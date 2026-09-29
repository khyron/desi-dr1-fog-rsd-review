#!/usr/bin/env python3
"""Generate local-only RSD extensions for uncached current DR1 viewer targets.

This tool reconstructs one tracer×hemisphere at a time from the complete local
DESI DR1 LSS v1.5 data/random selection. It writes only under tools/rsd and
never modifies viewer assets or the existing reference caches.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
sys.path.insert(0, str(ROOT))
from rsd_pipeline import (  # noqa: E402
    TRACERS, _deduplicate_targetids, _fits_rows, _physical_mpc,
    _query_radial_delta, _source_path, config, local_sources, read_desi_header,
)
from audit_rsd_catalog_coverage import FIELDS
from build_desi_chunks import morton_order


def rss_gib() -> float:
    # Linux ru_maxrss is KiB. This tool is expected to run under the verified WSL venv.
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def sorted_cache(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path) as obj:
        ids = np.asarray(obj['targetid'], dtype='u8')
        delta = np.asarray(obj['delta_mpc'], dtype='f4')
    order = np.argsort(ids, kind='stable')
    ids, delta = ids[order], delta[order]
    if np.any(ids[1:] == ids[:-1]) or not np.isfinite(delta).all():
        raise ValueError(f'{path}: duplicate TARGETID or nonfinite reference delta')
    return ids, delta


def lookup(ids: np.ndarray, ref_ids: np.ndarray, ref_delta: np.ndarray):
    at = np.searchsorted(ref_ids, ids)
    hit = (at < len(ref_ids)) & (ref_ids[np.minimum(at, len(ref_ids)-1)] == ids)
    delta = np.zeros(len(ids), dtype='f4')
    delta[hit] = ref_delta[at[hit]]
    return delta, hit


def compare(pred: np.ndarray, ref: np.ndarray) -> dict:
    if not len(pred):
        return {'count': 0}
    diff = pred.astype('f8') - ref.astype('f8')
    absdiff = np.abs(diff)
    nontrivial = np.abs(ref) >= 0.5
    sign_agree = None
    if nontrivial.any():
        sign_agree = float(np.mean(np.sign(pred[nontrivial]) == np.sign(ref[nontrivial])))
    corr = float(np.corrcoef(pred, ref)[0, 1]) if len(pred) > 1 else None
    return {
        'count': int(len(pred)), 'pearsonR': corr,
        'meanBiasMpc': float(diff.mean()), 'medianAbsDifferenceMpc': float(np.median(absdiff)),
        'p90AbsDifferenceMpc': float(np.percentile(absdiff, 90)),
        'p99AbsDifferenceMpc': float(np.percentile(absdiff, 99)),
        'signAgreementForAbsReferenceAtLeast0p5Mpc': sign_agree,
        'signComparedCount': int(nontrivial.sum()),
        'referenceDeltaPercentilesMpc': np.percentile(ref, [1, 50, 99]).tolist(),
        'extensionDeltaPercentilesMpc': np.percentile(pred, [1, 50, 99]).tolist(),
    }


def load_viewer_catalog(assets: Path, chunks_dir: Path, index_dir: Path):
    manifest = json.loads((chunks_dir / 'manifest.json').read_text(encoding='utf-8'))
    total = int(manifest['count'])
    index_report = json.loads((index_dir / 'index-report.json').read_text(encoding='utf-8'))
    if int(index_report.get('count', -1)) != total or not all(c.get('geometryExact') for c in index_report.get('chunks', [])):
        raise ValueError('private viewer TARGETID index is not certified exact against current geometry')
    units = float(manifest['unitsPerMpc'])
    raw_parts, z_parts = [], []
    for name in manifest['tracers']:
        path = assets / f'{name}.desi'
        h = read_desi_header(path); count = int(h['count'])
        if h['version'] != 2:
            raise ValueError(f'{path}: expected exact current DR1 geometry v2')
        raw = np.memmap(path, dtype='<i2', mode='r', offset=48, shape=(count, 3))
        raw_parts.append(np.rint(raw.astype('f4') / h['unitsPerMpc'] * units).astype('<i2'))
        zq = np.memmap(path, dtype='<u2', mode='r', offset=48 + count * 6, shape=(count,))
        z_parts.append(h['zMin'] + zq.astype('f4') * ((h['zMax'] - h['zMin']) / 65535.0))
    raw_xyz = np.concatenate(raw_parts)
    z_all = np.concatenate(z_parts).astype('f4')
    del raw_parts, z_parts
    hash_key = np.abs(raw_xyz[:, 0].astype('i8') * 17 + raw_xyz[:, 2].astype('i8') * 59 + raw_xyz[:, 1].astype('i8') * 101)
    pos_all = np.empty((total, 3), dtype='f4')
    z_view = np.empty(total, dtype='f4')
    ids_all = np.empty(total, dtype='u8')
    cursor = 0
    for chunk_i, chunk in enumerate(manifest['chunks']):
        chosen = np.flatnonzero(hash_key % len(manifest['chunks']) == chunk_i)
        chosen = chosen[morton_order(raw_xyz[chosen])]
        ids = np.load(index_dir / f'targetid_{chunk_i:03d}.npy', mmap_mode='r')
        if len(chosen) != len(ids) or len(ids) != int(chunk['count']):
            raise ValueError(f'chunk {chunk_i}: current geometry and exact TARGETID index disagree')
        end = cursor + len(chosen)
        pos_all[cursor:end] = raw_xyz[chosen].astype('f4') / units
        z_view[cursor:end] = z_all[chosen]
        ids_all[cursor:end] = ids
        cursor = end
    del raw_xyz, z_all, hash_key
    if cursor != total or len(np.unique(ids_all)) != total:
        raise ValueError('viewer catalog count or TARGETID uniqueness gate failed')
    return manifest, pos_all, z_view, ids_all


def mesh_matches(recon, report: dict, params: dict, cell_h: float, smoothing_mpc: float, tol: float):
    mesh = report.get('mesh', {})
    actual_nmesh = np.asarray(recon.nmesh, dtype='i8')
    actual_size = np.asarray(recon.boxsize, dtype='f8')
    actual_center = np.asarray(recon.boxcenter, dtype='f8')
    checks = {
        'nmeshExact': bool(np.array_equal(actual_nmesh, np.asarray(mesh.get('nmesh'), dtype='i8'))),
        'boxsizeWithinMpcTolerance': bool(np.allclose(actual_size, mesh.get('boxsizeMpc'), rtol=1e-8, atol=tol)),
        'boxcenterWithinMpcTolerance': bool(np.allclose(actual_center, mesh.get('boxcenterMpc'), rtol=1e-8, atol=tol)),
        'cellSizeMpcHMatches': abs(float(report.get('cellSizeMpcH', -999)) - cell_h) < 1e-12,
        'smoothingMpcMatches': abs(float(report.get('smoothingMpc', -999)) - smoothing_mpc) < tol,
        'biasMatches': abs(float(report.get('bias', -999)) - float(params['bias'])) < 1e-12,
        'growthRateMatches': abs(float(report.get('f', -999)) - float(params['f'])) < 1e-12,
        'fieldMatches': report.get('field') == 'rsd',
        'localLosMatches': report.get('los') == 'local',
    }
    return checks, actual_nmesh, actual_size, actual_center


def support_for_region(random_path: Path, params: dict, cosmo, candidates: np.ndarray,
                       batch: int, validation_n: int, percentile: float):
    print(f'[support {random_path.name}] reading full z/weight-selected randoms',flush=True)
    rows = _deduplicate_targetids(_fits_rows(random_path, cosmo, 0, params['zMin'], params['zMax']))
    rp = np.asarray(rows['position'], dtype='f4')
    print(f'[support {random_path.name}] selected={len(rp):,}; building random KD-tree',flush=True)
    if not len(rp):
        raise ValueError(f'{random_path}: no selected randoms')
    train = rp[::2]
    valid_all = rp[1::2]
    nvalid = min(validation_n, len(valid_all))
    valid = valid_all[:nvalid]
    tree = cKDTree(train, compact_nodes=True, balanced_tree=True)
    val_d, _ = tree.query(valid, k=1, workers=-1)
    radius = float(np.percentile(val_d, percentile))
    support = np.zeros(len(candidates), dtype=bool)
    nearest_percentiles = np.percentile(val_d, [50, 90, 95, 99, 100]).tolist()
    for start in range(0, len(candidates), batch):
        end = min(len(candidates), start + batch)
        dist, _ = tree.query(candidates[start:end], k=1, distance_upper_bound=radius, workers=-1)
        support[start:end] = np.isfinite(dist)
        del dist
        if end == len(candidates) or end % max(batch * 10, batch) == 0:
            print(f'[support {random_path.name}] queried {end:,}/{len(candidates):,} candidate positions',flush=True)
    result = {
        'randomRowsSourceAfterZWeightAndIDDedup': int(len(rp)),
        'randomRowsRead': int(rows['rowsRead']), 'randomStride': int(rows['stride']),
        'trainRandomRows': int(len(train)), 'independentValidationRows': int(nvalid),
        'validationNearestDistancePercentilesMpc': nearest_percentiles,
        'supportRadiusMpc': radius, 'supportPercentile': percentile,
        'validationAcceptedFraction': float(np.mean(val_d <= radius)),
        'interpretation': 'empirical random-cloud proximity proxy, not an official angular/radial mask polygon',
    }
    del tree, rp, train, valid, valid_all, val_d, rows
    return support, result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir', type=Path, default=Path('/mnt/d/desi/data_external/RSD'))
    ap.add_argument('--cache-dir', type=Path, default=Path('/mnt/d/desi/data_external/RSD/derived/rsd-v1-exploratory'))
    ap.add_argument('--config', type=Path, default=ROOT / 'tools/rsd/rsd_config.json')
    ap.add_argument('--assets', type=Path, default=ROOT / 'assets')
    ap.add_argument('--chunks', type=Path, default=ROOT / 'ply/chunks')
    ap.add_argument('--index', type=Path, default=ROOT / 'ply/reconstruction_dr1/private_index')
    ap.add_argument('--out', type=Path, default=ROOT / 'tools/rsd/extensions')
    ap.add_argument('--tracer', choices=TRACERS, action='append', help='default: all tracers sequentially')
    ap.add_argument('--region', choices=('NGC', 'SGC'), action='append', help='default: both hemispheres')
    ap.add_argument('--cell-size-mpc-h', type=float, default=None)
    ap.add_argument('--iterations', type=int, default=3)
    ap.add_argument('--validation-sample', type=int, default=50000)
    ap.add_argument('--mask-validation-sample', type=int, default=250000)
    ap.add_argument('--use-reference-random-stride', action='store_true',
                    help='reproduce each cache report random stride for reconstruction only; mask support still uses full randoms')
    ap.add_argument('--mask-percentile', type=float, default=99.0)
    ap.add_argument('--batch-size', type=int, default=50000)
    ap.add_argument('--max-rss-gib', type=float, default=5.5)
    ap.add_argument('--mesh-tolerance-mpc', type=float, default=1e-4)
    ap.add_argument('--min-sign-agreement', type=float, default=0.99)
    ap.add_argument('--min-correlation', type=float, default=0.995)
    ap.add_argument('--max-median-abs-difference-mpc', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=20260923)
    args = ap.parse_args()

    if not np.isfinite(args.max_rss_gib) or args.max_rss_gib <= 0:
        raise ValueError('--max-rss-gib must be a positive finite memory limit')
    if args.mask_percentile != 99.0:
        raise ValueError('the support threshold is fixed at p99 for this run; do not relax it')
    if args.min_sign_agreement < 0.99 or args.min_correlation < 0.995 or args.max_median_abs_difference_mpc > 0.1:
        raise ValueError('validation gates cannot be relaxed')

    from astropy.cosmology import Planck18
    from pyrecon import IterativeFFTReconstruction
    import pyrecon

    cfg = config(args.config); files = local_sources(args.data_dir)
    chosen = args.tracer or list(cfg['tracers'])
    regions = args.region or ['NGC', 'SGC']
    if len(set(chosen)) != len(chosen) or len(set(regions)) != len(regions):
        raise ValueError('do not repeat tracer or region selectors')
    args.out.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out / 'extension-manifest.json'
    previous_fields=[]
    if manifest_path.exists():
        try:
            old_manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
            previous_fields=[f for f in old_manifest.get('fields', [])
                if f.get('status')=='complete' and f.get('extension') and Path(f['extension']).exists()]
        except (OSError,json.JSONDecodeError):
            previous_fields=[]
    run = {
        'status': 'running', 'generatedAt': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'catalog': 'DESI DR1 LSS iron LSScats v1.5', 'viewer': 'current exact-indexed DR1 v2 catalog',
        'pyreconVersion': getattr(pyrecon, '__version__', 'unknown'),
        'config': str(args.config), 'cellSizeMpcH': float(args.cell_size_mpc_h or cfg['gridCellSizeMpcH']),
        'field': 'rsd', 'lineOfSight': 'local', 'radialProjection': 'delta_r = -dot(read_shifts(pos, field="rsd"), pos/|pos|)',
        'selectionPolicy': 'current viewer z in configured half-open interval; exact TARGETID cache union per same tracer excluded; each viewer ID is assigned once by disjoint z bins',
        'deduplicationPolicy': 'official data/random TARGETID duplicates retain highest WEIGHT*(WEIGHT_FKP if present), stable first row on tie',
        'reconstructionRandomPolicy': ('deterministic per-field stride copied from reference report; full randoms are used independently for support proxy' if args.use_reference_random_stride else 'full random catalog for reconstruction and support'),
        'supportPolicy': 'unique support by hemisphere: nearest distance <= p99 on disjoint random split and within reconstructed mesh box; proxy only, not official mask',
        'resourceGateGiB': float(args.max_rss_gib), 'fields': previous_fields,
        'requestedTracers':chosen,'requestedRegions':regions,
        'notProductionAsset': True,
    }
    def record_field(entry):
        for i, old in enumerate(run['fields']):
            if old.get('tracer')==entry.get('tracer') and old.get('region')==entry.get('region'):
                run['fields'][i]=entry
                break
        else:
            run['fields'].append(entry)
    manifest_path.write_text(json.dumps(run, indent=2) + '\n', encoding='utf-8')

    catalog_manifest, viewer_pos, viewer_z, viewer_ids = load_viewer_catalog(args.assets, args.chunks, args.index)
    if rss_gib() > args.max_rss_gib:
        run['status']='blocked_memory_limit_during_viewer_load'; run['peakRssGiB']=rss_gib()
        manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
        raise MemoryError(f'viewer catalog load peaked at {rss_gib():.3f} GiB, above the {args.max_rss_gib:.3f} GiB limit')
    run['viewerRows'] = int(len(viewer_ids)); run['viewerTargetIdsUnique'] = True
    run['viewerGeometryUnitsPerMpc'] = float(catalog_manifest['unitsPerMpc'])
    run['exactPrivateIndex'] = str(args.index)
    run['iterations'] = int(args.iterations)
    run['boxpad'] = 1.2
    manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
    ref_reports: dict[tuple[str, str], dict] = {}
    ref_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for tracer in chosen:
        ref_cache[tracer] = {}
        pieces=[]
        for region in ('NGC','SGC'):
            report_path=args.cache_dir/f'{tracer}_{region}_report.json'
            cache_path=args.cache_dir/f'{tracer}_{region}_rsd_query_v1.npz'
            report=json.loads(report_path.read_text(encoding='utf-8'))
            ids, delta=sorted_cache(cache_path)
            ref_reports[(tracer,region)]=report
            pieces.append((ids,delta))
        all_ids=np.concatenate([p[0] for p in pieces]); all_delta=np.concatenate([p[1] for p in pieces])
        order=np.argsort(all_ids,kind='stable'); all_ids,all_delta=all_ids[order],all_delta[order]
        if np.any(all_ids[1:]==all_ids[:-1]):
            raise ValueError(f'{tracer}: reference caches overlap by TARGETID across NGC/SGC')
        ref_cache[tracer]=(all_ids,all_delta)

    start_all=time.time()
    complete_keys={(f['tracer'],f['region']) for f in run['fields']}
    for tracer in chosen:
        params=cfg['tracers'][tracer]
        zmin,zmax=float(params['zMin']),float(params['zMax'])
        field_mask=(viewer_z>=zmin)&(viewer_z<zmax)
        field_indices=np.flatnonzero(field_mask)
        ref_ids,ref_delta=ref_cache[tracer]
        at=np.searchsorted(ref_ids,viewer_ids[field_indices])
        cached=(at<len(ref_ids))&(ref_ids[np.minimum(at,len(ref_ids)-1)]==viewer_ids[field_indices])
        missing_indices=field_indices[~cached]
        missing_pos=np.asarray(viewer_pos[missing_indices],dtype='f4')
        missing_ids=np.asarray(viewer_ids[missing_indices],dtype='u8')
        if len(np.unique(missing_ids)) != len(missing_ids):
            raise ValueError(f'{tracer}: uncached viewer candidate TARGETIDs are not unique')
        region_support={}; support_meta={}
        print(f'[{tracer}] eligible={len(field_indices):,}; exact cache={int(cached.sum()):,}; uncached candidates={len(missing_ids):,}',flush=True)
        for region in ('NGC','SGC'):
            run['supportProgress']={'tracer':tracer,'region':region,'candidateCount':int(len(missing_ids)),
                                    'stage':'read_randoms_build_tree_query_with_p99_upper_bound'}
            manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
            random_path=_source_path(files,tracer,region,'random')
            region_support[region],support_meta[region]=support_for_region(
                random_path,params,Planck18,missing_pos,args.batch_size,
                args.mask_validation_sample,args.mask_percentile)
            support_meta[region]['supportCandidatesBeforeMesh'] = int(region_support[region].sum())
            run['supportProgress'].update({'stage':'complete','supportRadiusMpc':support_meta[region]['supportRadiusMpc'],
                                           'supportedCandidates':int(region_support[region].sum()),'peakRssGiB':rss_gib()})
            manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
            print(f'[{tracer} {region}] support radius={support_meta[region]["supportRadiusMpc"]:.3f} Mpc; candidates={int(region_support[region].sum()):,}; RSS={rss_gib():.2f} GiB',flush=True)
            if rss_gib()>args.max_rss_gib:
                run['status']='blocked_memory_limit_during_support'
                run['peakRssGiB']=rss_gib()
                run['lastCompletedSupport']={'tracer':tracer,'region':region,'support':support_meta[region]}
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                raise MemoryError(f'RSS {rss_gib():.3f} GiB exceeded {args.max_rss_gib:.3f} GiB before reconstruction')
        both=region_support['NGC']&region_support['SGC']
        unique_region={r:region_support[r]&~region_support['SGC' if r=='NGC' else 'NGC'] for r in ('NGC','SGC')}
        del region_support

        for region in regions:
            if (tracer,region) in complete_keys:
                print(f'[{tracer} {region}] reusing completed extension from prior manifest',flush=True)
                continue
            field_t=time.time(); field_entry={'tracer':tracer,'region':region,'status':'running',
                'selection':{'zMin':zmin,'zMaxExclusive':zmax,'weight':'WEIGHT * WEIGHT_FKP if present'},
                'supportProxy':support_meta[region],
                'ambiguousBothRegionSupportCandidates':int(both.sum()),
                'inputFiles':{},'validation':{},'extension':None}
            # Ensure support for this region is unique, and later also lies within
            # the exact mesh that is rebuilt in this execution.
            region_candidate_mask=unique_region[region]
            region_entry=ref_reports[(tracer,region)]
            data_path=_source_path(files,tracer,region,'data')
            random_path=_source_path(files,tracer,region,'random')
            data=_deduplicate_targetids(_fits_rows(data_path,Planck18,0,zmin,zmax))
            ref_random_stride=int(region_entry['randoms'].get('stride',1)) if args.use_reference_random_stride else None
            randoms=_deduplicate_targetids(_fits_rows(random_path,Planck18,0,zmin,zmax,row_stride=ref_random_stride))
            # Source/runtime inputs must be identical to the full-data reference
            # reconstruction whose mesh and cache are being extended.
            source_gates={
                'dataRowsKeptMatchReport': int(data['rowsKept'])==int(region_entry['data']['rowsKept']),
                'randomRowsKeptMatchReport': int(randoms['rowsKept'])==int(region_entry['randoms']['rowsKept']),
                'randomRowsReadMatchReport': int(randoms['rowsRead'])==int(region_entry['randoms']['rowsRead']),
                'dataStrideIsFull': int(data['stride'])==1,
                'randomStrideMatchesReport': int(randoms['stride'])==int(region_entry['randoms']['stride']),
                'zRangeMatches': np.allclose([zmin,zmax],region_entry['zRange'],atol=1e-12,rtol=0),
                'cellSizeMatchesConfig': abs(float(region_entry['cellSizeMpcH'])-float(args.cell_size_mpc_h or cfg['gridCellSizeMpcH']))<1e-12,
                'smoothingMatchesConfig': abs(float(region_entry['smoothingMpcH'])-float(params['smoothingMpcH']))<1e-12,
                'dataDedupCountMatchesReport': int(data.get('duplicatesDiscarded',0))==int(region_entry['data'].get('duplicatesDiscarded',0)),
                'randomDedupCountMatchesReport': int(randoms.get('duplicatesDiscarded',0))==int(region_entry['randoms'].get('duplicatesDiscarded',0)),
            }
            field_entry['inputFiles']={
                'data':{'path':str(data_path),'rowsRead':int(data['rowsRead']),'rowsKept':int(data['rowsKept']),
                        'duplicatesDiscarded':int(data.get('duplicatesDiscarded',0)),'stride':int(data['stride'])},
                'randoms':{'path':str(random_path),'rowsRead':int(randoms['rowsRead']),'rowsKept':int(randoms['rowsKept']),
                           'duplicatesDiscarded':int(randoms.get('duplicatesDiscarded',0)),'stride':int(randoms['stride'])},
            }
            field_entry['runtimeParameters']={'bias':float(params['bias']),'f':float(params['f']),
                'cellSizeMpcH':float(args.cell_size_mpc_h or cfg['gridCellSizeMpcH']),
                'cellSizeMpc':_physical_mpc(float(args.cell_size_mpc_h or cfg['gridCellSizeMpcH']),Planck18),
                'smoothingMpcH':float(params['smoothingMpcH']),
                'smoothingMpc':_physical_mpc(float(params['smoothingMpcH']),Planck18),
                'iterations':int(args.iterations),'boxpad':1.2,'field':'rsd','los':'local'}
            field_entry['sourceRuntimeGates']=source_gates
            if not all(source_gates.values()):
                field_entry['status']='blocked_source_or_config_mismatch'
                record_field(field_entry); run['status']='blocked'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                print(f'[{tracer} {region}] blocked: source/config mismatch',flush=True)
                return 2

            cell_h=float(args.cell_size_mpc_h or cfg['gridCellSizeMpcH'])
            cell_mpc=_physical_mpc(cell_h,Planck18)
            smoothing_mpc=_physical_mpc(float(params['smoothingMpcH']),Planck18)
            positions=np.concatenate((data['position'],randoms['position']))
            recon=IterativeFFTReconstruction(f=float(params['f']),bias=float(params['bias']),los=None,
                cellsize=cell_mpc,boxpad=1.2,positions=positions)
            recon.assign_data(data['position'],weights=data['weight'])
            recon.assign_randoms(randoms['position'],weights=randoms['weight'])
            recon.set_density_contrast(smoothing_radius=smoothing_mpc)
            recon.run(niterations=args.iterations)
            del positions
            mesh_gates,actual_nmesh,actual_size,actual_center=mesh_matches(
                recon,region_entry,params,cell_h,smoothing_mpc,args.mesh_tolerance_mpc)
            field_entry['mesh']={'actualNmesh':actual_nmesh.tolist(),'referenceNmesh':region_entry['mesh']['nmesh'],
                'actualBoxsizeMpc':actual_size.tolist(),'referenceBoxsizeMpc':region_entry['mesh']['boxsizeMpc'],
                'actualBoxcenterMpc':actual_center.tolist(),'referenceBoxcenterMpc':region_entry['mesh']['boxcenterMpc'],
                'gates':mesh_gates}
            field_entry['peakRssGiBAfterReconstruction']=rss_gib()
            print(f'[{tracer} {region}] reconstruction complete; mesh gates={all(mesh_gates.values())}; RSS={rss_gib():.2f} GiB',flush=True)
            if rss_gib()>args.max_rss_gib:
                field_entry['status']='blocked_memory_limit'
                run['fields'].append(field_entry); run['status']='blocked_memory_limit'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                raise MemoryError(f'RSS {rss_gib():.3f} GiB exceeded {args.max_rss_gib:.3f} GiB')
            if not all(mesh_gates.values()):
                field_entry['status']='blocked_mesh_reference_mismatch'
                record_field(field_entry); run['status']='blocked'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                return 3

            half=np.asarray(recon.boxsize,dtype='f8')/2
            box_center=np.asarray(recon.boxcenter,dtype='f8')
            inside_box=np.all(np.abs(missing_pos-box_center)<=half,axis=1)
            region_candidates=region_candidate_mask&inside_box
            candidate_rows=int(region_candidates.sum())
            field_entry['supportProxy']['insideRebuiltMesh']=int((region_support_mask := (region_candidate_mask&inside_box)).sum())
            field_entry['supportProxy']['uniqueCandidatesBeforeMesh']=int(region_candidate_mask.sum())
            field_entry['supportProxy']['uniqueCandidatesInsideMesh']=candidate_rows
            # Sample data positions retained in the exact reference cache to
            # check direction, sign and amplitudes for the live reconstructed field.
            rng=np.random.default_rng(args.seed+sum(ord(c) for c in tracer+region))
            val_n=min(args.validation_sample,len(data['position']))
            val_idx=np.sort(rng.choice(len(data['position']),size=val_n,replace=False))
            val_pred,_=_query_radial_delta(recon,data['position'][val_idx],args.batch_size)
            val_ref,val_hit=lookup(data['targetid'][val_idx],ref_ids,ref_delta)
            val_cmp=compare(val_pred[val_hit],val_ref[val_hit])
            val_cmp.update({'sampledDataRows':int(val_n),'cacheMatches':int(val_hit.sum()),
                'cacheMatchFraction':float(val_hit.mean()) if len(val_hit) else 0.0,
                'meaning':'same-input full-data field comparison to existing cache; sensitivity check, not independent velocity truth'})
            sign_agreement=val_cmp.get('signAgreementForAbsReferenceAtLeast0p5Mpc')
            validation_gates={
                'hasCommonCacheTargets':int(val_hit.sum())>0,
                'correlationPass':val_cmp.get('pearsonR') is not None and val_cmp['pearsonR']>=args.min_correlation,
                'signAgreementPass':sign_agreement is not None and sign_agreement>=args.min_sign_agreement,
                'medianDifferencePass':val_cmp.get('medianAbsDifferenceMpc',float('inf'))<=args.max_median_abs_difference_mpc,
            }
            field_entry['validation']={'commonTargetsVsFullDataCache':val_cmp,'gates':validation_gates}
            field_entry['peakRssGiBAfterValidation']=rss_gib()
            if not all(validation_gates.values()):
                field_entry['status']='blocked_cache_validation_gate'
                record_field(field_entry); run['status']='blocked'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                print(f'[{tracer} {region}] cache QA failed: {json.dumps(validation_gates)}',flush=True)
                return 4

            qpos=missing_pos[region_candidates]
            qids=missing_ids[region_candidates]
            if len(qids) and np.any(np.diff(np.sort(qids))==0):
                raise ValueError(f'{tracer} {region}: duplicate extension TARGETID after all gates')
            delta,_=_query_radial_delta(recon,qpos,args.batch_size)
            if rss_gib()>args.max_rss_gib:
                field_entry['status']='blocked_memory_limit_during_query'
                record_field(field_entry); run['status']='blocked_memory_limit'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                raise MemoryError(f'RSS {rss_gib():.3f} GiB exceeded {args.max_rss_gib:.3f} GiB during query')
            ext_path=args.out/f'{tracer}_{region}_uncached_extension_v1.npz'
            # Write only after source, mesh, resource, identity and sign gates pass.
            np.savez(ext_path,targetid=qids,delta_mpc=delta)
            field_entry.update({
                'status':'complete','extension':str(ext_path),'extensionRows':int(len(qids)),
                'deltaPercentilesMpc':np.percentile(delta,[0,1,50,99,100]).tolist() if len(delta) else [],
                'deltaSignCounts':{'negative':int((delta<0).sum()),'zero':int((delta==0).sum()),'positive':int((delta>0).sum())},
                'finalPeakRssGiB':rss_gib(),'elapsedSeconds':time.time()-field_t,
                'referenceReport':str(args.cache_dir/f'{tracer}_{region}_report.json'),
                'noAssetOrReferenceCacheModified':True,
            })
            field_entry['outputSha256']=sha256_file(ext_path)
            record_field(field_entry)
            complete_keys.add((tracer,region))
            run['status']='running'
            manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
            print(f'[{tracer} {region}] wrote {len(qids):,} uncached proxy-supported deltas; peak RSS={rss_gib():.2f} GiB; validation={json.dumps(validation_gates)}',flush=True)
            # Free all large field-specific state before the next region.
            del recon,data,randoms,qpos,qids,delta,val_pred,val_ref,val_hit,val_idx,inside_box
            if rss_gib()>args.max_rss_gib:
                run['status']='blocked_memory_limit_after_field'
                manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
                raise MemoryError(f'RSS {rss_gib():.3f} GiB exceeded {args.max_rss_gib:.3f} GiB after field')
        del unique_region,both,missing_pos,missing_ids,missing_indices,field_indices
    run['status']='complete'
    run['elapsedSeconds']=time.time()-start_all
    run['peakRssGiB']=rss_gib()
    manifest_path.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':run['status'],'fields':len(run['fields']),'elapsedSeconds':run['elapsedSeconds'],'peakRssGiB':run['peakRssGiB']},indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
