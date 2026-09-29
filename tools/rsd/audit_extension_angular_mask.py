#!/usr/bin/env python3
"""Audit extension angular support against independent halves of DR1 LSS randoms.

This is a random-catalog diagnostic, not a certified DESI veto mask.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

from build_dr1_chunk_index import header, morton_order

ROOT = Path(__file__).resolve().parents[2]


def unit(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    c = np.cos(dec)
    return np.column_stack((c * np.cos(ra), c * np.sin(ra), np.sin(dec))).astype('f4')


def viewer_sample_positions(sample_ids):
    manifest = json.loads((ROOT / 'ply/chunks/manifest.json').read_text())
    units = float(manifest['unitsPerMpc'])
    positions = []
    for name in manifest['tracers']:
        path = ROOT / 'assets' / (name + '.desi')
        count, own_units, _, _ = header(path)
        raw = np.memmap(path, dtype='<i2', mode='r', offset=48, shape=(count, 3))
        positions.append(np.rint(raw.astype('f4') / own_units * units).astype('<i2'))
    pos = np.concatenate(positions)
    hashes = np.abs(pos[:, 0].astype('i8') * 17 + pos[:, 2].astype('i8') * 59 + pos[:, 1].astype('i8') * 101)
    wanted = np.sort(np.unique(np.concatenate(list(sample_ids.values()))))
    found = np.zeros(len(wanted), dtype=bool)
    xyz = np.zeros((len(wanted), 3), dtype='f4')
    for i, chunk in enumerate(manifest['chunks']):
        selected = np.flatnonzero(hashes % len(manifest['chunks']) == i)
        selected = selected[morton_order(pos[selected])]
        ids = np.load(ROOT / 'ply/reconstruction_dr1/private_index' / f'targetid_{i:03d}.npy', mmap_mode='r')
        if len(ids) != chunk['count'] or len(ids) != len(selected):
            raise ValueError(f'chunk {i}: geometry/index mismatch')
        slot = np.searchsorted(wanted, ids)
        hit = (slot < len(wanted)) & (wanted[np.minimum(slot, len(wanted) - 1)] == ids)
        xyz[slot[hit]] = pos[selected[hit]].astype('f4') / units
        found[slot[hit]] = True
    if not found.all():
        raise ValueError(f'{(~found).sum()} sampled extension IDs absent from viewer index')
    return wanted, xyz


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sample-per-field', type=int, default=10000)
    ap.add_argument('--random-sample', type=int, default=500000)
    ap.add_argument('--random-index', type=int, choices=(0, 1), default=0)
    ap.add_argument('--out', type=Path)
    a = ap.parse_args()
    if a.out is None:
        suffix = '' if a.random_index == 0 else '-index1'
        a.out = ROOT / f'tools/rsd/extensions/angular-mask-audit{suffix}.json'
    manifest = json.loads((ROOT / 'tools/rsd/extensions/extension-manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('extension manifest is incomplete')
    rng = np.random.default_rng(20260924)
    samples = {}
    for entry in manifest['fields']:
        key = (entry['tracer'], entry['region'])
        path = ROOT / 'tools/rsd/extensions' / f'{key[0]}_{key[1]}_uncached_extension_v1.npz'
        with np.load(path, allow_pickle=False) as z:
            ids = z['targetid']
            samples[key] = np.asarray(ids[rng.choice(len(ids), min(a.sample_per_field, len(ids)), replace=False)], dtype='u8')
    wanted, xyz = viewer_sample_positions(samples)
    norm = np.linalg.norm(xyz, axis=1)
    sky = (xyz / norm[:, None]).astype('f4')
    result = {'method': 'angular nearest-neighbor and eighth-neighbor support; disjoint random calibration/validation halves',
              'status': 'needs-independent-mask-validation', 'samplePerField': a.sample_per_field,
              'randomIndex': a.random_index, 'randomSampleCap': a.random_sample, 'fields': []}
    flagged = {}
    for entry in manifest['fields']:
        key = (entry['tracer'], entry['region'])
        source = entry['inputFiles']['randoms']['path']
        path = Path('D:/' + source[len('/mnt/d/'):]) if source.startswith('/mnt/d/') else Path(source)
        if a.random_index == 1:
            path = ROOT / 'tools/rsd/validation_sources' / path.name.replace('_0_clustering.ran.fits', '_1_clustering.ran.fits')
        lo, hi = entry['selection']['zMin'], entry['selection']['zMaxExclusive']
        with fits.open(path, memmap=True) as hdus:
            tab = hdus[1].data
            z = np.asarray(tab['Z'])
            weight = np.asarray(tab['WEIGHT']) * np.asarray(tab['WEIGHT_FKP'])
            valid = np.flatnonzero((z >= lo) & (z < hi) & np.isfinite(weight) & (weight > 0))
            if len(valid) > a.random_sample:
                valid = np.sort(rng.choice(valid, a.random_sample, replace=False))
            random_sky = unit(np.asarray(tab['RA'][valid]), np.asarray(tab['DEC'][valid]))
        train, holdout = random_sky[::2], random_sky[1::2]
        tree = cKDTree(train)
        hold_dist = tree.query(holdout, k=8, workers=-1)[0]
        thresholds = np.percentile(hold_dist, 99, axis=0)
        slot = np.searchsorted(wanted, samples[key])
        ext_dist = tree.query(sky[slot], k=8, workers=-1)[0]
        near = ext_dist[:, 0] <= thresholds[0]
        dense = ext_dist[:, 7] <= thresholds[7]
        flagged[f'{key[0]}_{key[1]}'] = samples[key][~(near & dense)]
        item = {'tracer': key[0], 'region': key[1], 'sampledExtensionRows': len(slot),
                'filteredRandomRowsSampled': len(random_sky),
                'calibrationHoldoutRows': len(holdout),
                'angularNearestP99Arcmin': float(np.rad2deg(2*np.arcsin(thresholds[0]/2))*60),
                'angularEighthP99Arcmin': float(np.rad2deg(2*np.arcsin(thresholds[7]/2))*60),
                'sampleNearFraction': float(near.mean()), 'sampleDenseFraction': float(dense.mean()),
                'sampleBothFraction': float((near & dense).mean()),
                'sampleRejectedNear': int((~near).sum()), 'sampleRejectedDensity': int((~dense).sum()),
                'sampleRejectedEither': int((~(near & dense)).sum())}
        result['fields'].append(item)
        print(json.dumps(item), flush=True)
    result['limitations'] = ['Random catalogue 0 samples the selection but is not an independent official angular veto-mask polygon.',
                             'Finite random sampling and p99 thresholds can reject valid edge points or accept holes.',
                             'This does not validate radial selection beyond the configured z interval or field physics.']
    result['sampledExtensionRowsTotal'] = sum(x['sampledExtensionRows'] for x in result['fields'])
    result['sampleRejectedEitherTotal'] = sum(x['sampleRejectedEither'] for x in result['fields'])
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    np.savez_compressed(a.out.with_name(a.out.stem + '-flagged-targetids.npz'), **flagged)


if __name__ == '__main__':
    main()
