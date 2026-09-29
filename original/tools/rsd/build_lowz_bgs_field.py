#!/usr/bin/env python3
"""Exploratory BGS BRIGHT 0.01<=z<0.40 field for the exact DR1 viewer.

Outputs are local candidate caches only. Scientific mask approval and PlayCanvas
deployment are separate gates.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
from extend_rsd_catalog import load_viewer_catalog  # noqa: E402
from rsd_pipeline import _deduplicate_targetids, _fits_rows, _query_radial_delta  # noqa: E402

SOURCE = ROOT / 'tools' / 'rsd' / 'lowz_sources'
OUT = SOURCE / 'field_candidate'
ZMIN, ZMAX = 0.01, 0.40


def source(region, kind):
    suffix = 'clustering.dat.fits' if kind == 'data' else '0_clustering.ran.fits'
    return SOURCE / f'BGS_BRIGHT_{region}_{suffix}'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def support_mask(random_path, candidates, percentile, batch):
    # The full catalogue is used for this proximity proxy, irrespective of the
    # deterministic stride used below to limit reconstruction memory.
    randoms = _deduplicate_targetids(_fits_rows(random_path, Planck18, 0, ZMIN, ZMAX))
    positions = np.asarray(randoms['position'], dtype='f4')
    train = positions[::2]
    holdout = positions[1::2][:250000]
    tree = cKDTree(train)
    holdout_distance = tree.query(holdout, k=1, workers=-1)[0]
    radius = float(np.percentile(holdout_distance, percentile))
    accepted = np.zeros(len(candidates), dtype=bool)
    for start in range(0, len(candidates), batch):
        end = min(start + batch, len(candidates))
        distance = tree.query(candidates[start:end], k=1,
                              distance_upper_bound=radius, workers=-1)[0]
        accepted[start:end] = np.isfinite(distance)
        print(f'support {random_path.name}: {end:,}/{len(candidates):,}', flush=True)
    report = {'fullSelectedRandomRows': len(positions), 'trainRows': len(train),
              'holdoutRows': len(holdout), 'holdoutP99RadiusMpc': radius,
              'acceptedBeforeRegionExclusivity': int(accepted.sum()),
              'meaning': '3D random proximity proxy, not an official mask'}
    del tree, positions, train, holdout, randoms
    return accepted, report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--random-stride', type=int, default=6)
    ap.add_argument('--batch-size', type=int, default=100000)
    ap.add_argument('--support-percentile', type=float, default=99.0)
    args = ap.parse_args()
    if args.random_stride < 1 or args.support_percentile != 99.0:
        raise ValueError('invalid reconstruction stride or support threshold')
    OUT.mkdir(exist_ok=True)
    manifest, positions, redshift, targetid = load_viewer_catalog(
        ROOT / 'assets', ROOT / 'ply' / 'chunks', ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    select = (redshift >= ZMIN) & (redshift < ZMAX)
    candidates = positions[select]
    ids = targetid[select]
    zs = redshift[select]
    del positions, redshift, targetid
    support, support_reports = {}, {}
    for region in ('NGC', 'SGC'):
        support[region], support_reports[region] = support_mask(
            source(region, 'random'), candidates, args.support_percentile, args.batch_size)
    both = support['NGC'] & support['SGC']
    by_region = {region: support[region] & ~both for region in support}
    run = {'status': 'running', 'catalog': 'DESI DR1 LSS iron LSScats v1.5 BGS_BRIGHT',
           'redshiftRange': [ZMIN, ZMAX], 'viewerSelectedRows': int(select.sum()),
           'ambiguousBothRegions': int(both.sum()), 'supportProxy': support_reports,
           'randomStride': args.random_stride, 'cellSizeMpcH': 16.0,
           'smoothingMpcH': 15.0, 'bias': 1.5, 'f': 0.682,
           'iterations': 3, 'boxpad': 1.2, 'field': 'rsd',
           'scientificMaskStatus': 'unvalidated', 'regions': []}
    from pyrecon import IterativeFFTReconstruction
    for region in ('NGC', 'SGC'):
        print(f'reconstructing BGS_BRIGHT {region}', flush=True)
        data_path, random_path = source(region, 'data'), source(region, 'random')
        data = _deduplicate_targetids(_fits_rows(data_path, Planck18, 0, ZMIN, ZMAX))
        randoms = _deduplicate_targetids(_fits_rows(random_path, Planck18, 0, ZMIN, ZMAX,
                                                    row_stride=args.random_stride))
        if not len(data['position']) or not len(randoms['position']):
            raise ValueError('empty source selection')
        setup = np.concatenate((data['position'], randoms['position']))
        recon = IterativeFFTReconstruction(f=0.682, bias=1.5, los=None,
            cellsize=16 / Planck18.h, boxpad=1.2, positions=setup)
        del setup
        recon.assign_data(data['position'], weights=data['weight'])
        recon.assign_randoms(randoms['position'], weights=randoms['weight'])
        recon.set_density_contrast(smoothing_radius=15 / Planck18.h)
        recon.run(niterations=3)
        mesh_center = np.asarray(recon.boxcenter, dtype='f8')
        mesh_halfsize = np.asarray(recon.boxsize, dtype='f8') / 2
        within_mesh = np.all(np.abs(candidates - mesh_center) <= mesh_halfsize, axis=1)
        chosen = by_region[region] & within_mesh
        values, _ = _query_radial_delta(recon, candidates[chosen], args.batch_size)
        if not np.isfinite(values).all() or np.max(np.abs(values), initial=0) > 150:
            raise ValueError('nonfinite or implausible low-z shift')
        output = OUT / f'BGS_BRIGHT_{region}_viewer_rsd_v1.npz'
        np.savez(output, targetid=ids[chosen], z=zs[chosen], delta_mpc=values)
        # Same-field exact-ID check: query source positions and compare the
        # resulting shift to viewer positions for shared source TARGETIDs.
        source_ids = np.asarray(data['targetid'], dtype='u8')
        sample_count = min(50000, len(source_ids))
        sample = np.linspace(0, len(source_ids)-1, sample_count, dtype='i8')
        source_delta, _ = _query_radial_delta(recon, data['position'][sample], args.batch_size)
        source_order = np.argsort(source_ids[sample])
        sample_ids = source_ids[sample][source_order]
        sample_delta = source_delta[source_order]
        at = np.searchsorted(sample_ids, ids[chosen])
        hit = (at < len(sample_ids)) & (sample_ids[np.minimum(at,len(sample_ids)-1)] == ids[chosen])
        difference = values[hit] - sample_delta[at[hit]]
        report = {'region': region, 'dataSource': str(data_path), 'dataSha256': sha256(data_path),
                  'randomSource': str(random_path), 'randomSha256': sha256(random_path),
                  'selectedDataRows': len(source_ids), 'selectedRandomRows': len(randoms['targetid']),
                  'dataStride': data['stride'], 'randomStride': randoms['stride'],
                  'duplicatesDiscardedData': data.get('duplicatesDiscarded',0),
                  'duplicatesDiscardedRandom': randoms.get('duplicatesDiscarded',0),
                  'meshN': np.asarray(recon.nmesh).tolist(),
                  'meshBoxsizeMpc': np.asarray(recon.boxsize).tolist(),
                  'meshCenterMpc': mesh_center.tolist(),
                  'candidateRowsAfterUniqueSupportAndMesh': int(chosen.sum()),
                  'exactIdValidationMatches': int(hit.sum()),
                  'exactIdValidationMedianAbsDifferenceMpc': float(np.median(np.abs(difference))) if len(difference) else None,
                  'deltaPercentilesMpc': np.percentile(values,[0,1,50,99,100]).tolist(),
                  'output': str(output), 'outputSha256': sha256(output)}
        run['regions'].append(report)
        (OUT / 'field-manifest.json').write_text(json.dumps(run,indent=2)+'\n')
        print(json.dumps({'region':region,'candidateRows':int(chosen.sum()),
                          'exactIdCheck':int(hit.sum())}),flush=True)
        del recon, data, randoms, values, chosen, within_mesh
    run['status'] = 'complete-for-technical-review'
    (OUT / 'field-manifest.json').write_text(json.dumps(run,indent=2)+'\n')


if __name__ == '__main__':
    main()
