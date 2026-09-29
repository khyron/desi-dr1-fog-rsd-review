#!/usr/bin/env python3
"""Audit local RSD extension files against their manifest and exact caches."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from audit_rsd_catalog_coverage import current_redshift_chunks

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / 'tools/rsd/extensions'
CACHE = Path('/mnt/d/desi/data_external/RSD/derived/rsd-v1-exploratory')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def main():
    manifest = json.loads((EXT / 'extension-manifest.json').read_text())
    report = {'status': 'complete', 'manifestStatus': manifest.get('status'), 'fields': [], 'totalRows': 0}
    by_tracer: dict[str, list[np.ndarray]] = {}
    extension_ids_by_tracer: dict[str, list[np.ndarray]] = {}
    extension_payloads: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    for entry in manifest['fields']:
        item = {'tracer': entry['tracer'], 'region': entry['region'], 'status': entry['status']}
        path = Path(entry.get('extension') or '')
        if not path.is_absolute():
            path = ROOT / path
        checks = {
            'fieldComplete': entry['status'] == 'complete',
            'sourceGatesPass': bool(entry.get('sourceRuntimeGates')) and all(entry['sourceRuntimeGates'].values()),
            'meshGatesPass': bool(entry.get('mesh', {}).get('gates')) and all(entry['mesh']['gates'].values()),
            'validationGatesPass': bool(entry.get('validation', {}).get('gates')) and all(entry['validation']['gates'].values()),
            'fileExists': path.is_file(),
        }
        if path.is_file():
            item['path'] = str(path)
            item['sha256Actual'] = sha256(path)
            item['sha256Expected'] = entry.get('outputSha256')
            checks['sha256MatchesManifest'] = item['sha256Actual'] == item['sha256Expected']
            with np.load(path, allow_pickle=False) as z:
                ids = np.asarray(z['targetid'], dtype='u8')
                delta = np.asarray(z['delta_mpc'], dtype='f4')
            extension_payloads[(entry['tracer'], entry['region'])] = (ids, delta)
            checks['rowCountMatchesManifest'] = len(ids) == int(entry.get('extensionRows', -1)) == len(delta)
            checks['idsUniqueWithinField'] = len(np.unique(ids)) == len(ids)
            checks['deltasFinite'] = bool(np.isfinite(delta).all())
            pieces = []
            for region in ('NGC', 'SGC'):
                with np.load(CACHE / f"{entry['tracer']}_{region}_rsd_query_v1.npz") as c:
                    pieces.append(np.asarray(c['targetid'], dtype='u8'))
            cache_ids = np.unique(np.concatenate(pieces))
            checks['outsideSameTracerCache'] = not bool(np.isin(ids, cache_ids, assume_unique=False).any())
            by_tracer.setdefault(entry['tracer'], []).append(ids)
            extension_ids_by_tracer.setdefault(entry['tracer'], []).append(ids)
            item.update({'rows': int(len(ids)), 'finiteDeltaRows': int(np.isfinite(delta).sum()),
                         'uniqueTargetIds': int(len(np.unique(ids))), 'checks': checks})
            report['totalRows'] += int(len(ids))
        else:
            item['checks'] = checks
        report['fields'].append(item)
    for tracer, groups in by_tracer.items():
        all_ids = np.concatenate(groups)
        report.setdefault('tracerChecks', {})[tracer] = {
            'extensionIdsUniqueAcrossRegions': bool(len(np.unique(all_ids)) == len(all_ids)),
            'rows': int(len(all_ids)),
        }
    # Global TARGETID union against all exact cache selections and the current viewer.
    all_extension_ids = np.concatenate([x for groups in extension_ids_by_tracer.values() for x in groups])
    unique_extension_ids, extension_counts = np.unique(all_extension_ids, return_counts=True)
    cross_tracer_ids = set()
    unique_by_tracer = {t: np.unique(np.concatenate(v)) for t, v in extension_ids_by_tracer.items()}
    tracer_names = list(unique_by_tracer)
    for i, a in enumerate(tracer_names):
        for b in tracer_names[i + 1:]:
            cross_tracer_ids.update(np.intersect1d(unique_by_tracer[a], unique_by_tracer[b], assume_unique=True).tolist())
    cache_parts = []
    cache_payloads = {}
    for tracer in extension_ids_by_tracer:
        for region in ('NGC', 'SGC'):
            with np.load(CACHE / f'{tracer}_{region}_rsd_query_v1.npz') as z:
                cids = np.asarray(z['targetid'], dtype='u8')
                cdelta = np.asarray(z['delta_mpc'], dtype='f4')
            order = np.argsort(cids, kind='stable')
            cids, cdelta = cids[order], cdelta[order]
            cache_parts.append(cids)
            cache_payloads[(tracer, region)] = (cids, cdelta)
    exact_cache_ids = np.unique(np.concatenate(cache_parts))
    ids_in_viewer = []
    index_dir = ROOT / 'ply/reconstruction_dr1/private_index'
    index_report = json.loads((index_dir / 'index-report.json').read_text())
    for i in range(len(index_report['chunks'])):
        ids_in_viewer.append(np.asarray(np.load(index_dir / f'targetid_{i:03d}.npy', mmap_mode='r'), dtype='u8'))
    viewer_ids = np.concatenate(ids_in_viewer)
    viewer_unique = np.unique(viewer_ids)
    catalog_manifest = json.loads((ROOT / 'ply/chunks/manifest.json').read_text())
    cached_current_viewer = 0
    cached_current_envelope = 0
    for chunk_i, z in current_redshift_chunks(ROOT / 'assets', catalog_manifest, index_dir):
        ids = np.asarray(np.load(index_dir / f'targetid_{chunk_i:03d}.npy', mmap_mode='r'), dtype='u8')
        hit = np.isin(ids, exact_cache_ids, assume_unique=True)
        eligible = np.zeros(len(ids), dtype=bool)
        for low, high in ((0.1, 0.4), (0.4, 0.8), (0.8, 1.1), (1.1, 1.6), (1.6, 2.1)):
            eligible |= (z >= low) & (z < high)
        cached_current_viewer += int(hit.sum())
        cached_current_envelope += int((hit & eligible).sum())
    extension_in_any_cache = np.isin(unique_extension_ids, exact_cache_ids, assume_unique=True)
    new_extension_ids = unique_extension_ids[~extension_in_any_cache]
    cross_cache_matrix = []
    cross_cache_overlap_ids = set()
    for (etracer, eregion), (eids, edelta) in extension_payloads.items():
        esort = np.argsort(eids, kind='stable')
        eids, edelta = eids[esort], edelta[esort]
        for (ctracer, cregion), (cids, cdelta) in cache_payloads.items():
            at = np.searchsorted(cids, eids)
            hit = (at < len(cids)) & (cids[np.minimum(at, len(cids) - 1)] == eids)
            n = int(hit.sum())
            if not n:
                continue
            cross_cache_overlap_ids.update(eids[hit].tolist())
            diff = edelta[hit].astype('f8') - cdelta[at[hit]].astype('f8')
            cross_cache_matrix.append({
                'extensionTracer': etracer, 'extensionRegion': eregion,
                'cacheTracer': ctracer, 'cacheRegion': cregion, 'count': n,
                'medianAbsDeltaDifferenceMpc': float(np.median(np.abs(diff))),
                'p90AbsDeltaDifferenceMpc': float(np.percentile(np.abs(diff), 90)),
                'p99AbsDeltaDifferenceMpc': float(np.percentile(np.abs(diff), 99)),
                'bitwiseEqualDeltaFraction': float(np.mean(diff == 0)),
                'interpretation': 'different-tracer cache hit; not presumed compatible with assigned extension selection',
            })
    prior_audit = json.loads((ROOT / 'tools/rsd/rsd-coverage-audit.json').read_text())
    envelope_rows = int(prior_audit['totals']['insideConfiguredRedshiftEnvelope'])
    extensions_in_viewer = np.isin(unique_extension_ids, viewer_unique, assume_unique=True)
    foreign_comparisons = sum(int(x['count']) for x in cross_cache_matrix)
    all_foreign_values_differ = bool(cross_cache_matrix) and all(
        float(x['bitwiseEqualDeltaFraction']) == 0.0 for x in cross_cache_matrix)
    report['globalCoverage'] = {
        'extensionRowsAcrossFields': int(len(all_extension_ids)),
        'extensionUniqueTargetIds': int(len(unique_extension_ids)),
        'duplicateExtensionRowsAcrossAnyTracer': int(len(all_extension_ids) - len(unique_extension_ids)),
        'crossTracerDuplicateTargetIds': int(len(cross_tracer_ids)),
        'crossTracerDuplicateExamples': [int(x) for x in sorted(cross_tracer_ids)[:20]],
        'cachedTargetIdsInCurrentViewerGlobalUnion': cached_current_viewer,
        'cacheIdsAcrossAllFields': int(len(exact_cache_ids)),
        'extensionUniqueIdsAlreadyInAnyTracerCache': int(extension_in_any_cache.sum()),
        'newExtensionTargetIdsUniqueAndNotCached': int(len(new_extension_ids)),
        'crossTracerCacheOverlapTargetIds': int(len(cross_cache_overlap_ids)),
        'crossTracerCacheDeltaComparisons': foreign_comparisons,
        'allComparedForeignCacheDeltasDifferBitwise': all_foreign_values_differ,
        'crossTracerCacheOverlapMatrix': cross_cache_matrix,
        'extensionIdsAllBelongToCurrentViewer': bool(extensions_in_viewer.all()),
        'currentViewerRows': int(len(viewer_unique)),
        'zEnvelopeRows': envelope_rows,
        'coveredUniqueRowsTotalCatalog': int(cached_current_viewer + len(new_extension_ids)),
        'coverageFractionTotalCatalog': float((cached_current_viewer + len(new_extension_ids)) / len(viewer_unique)),
        'cachedTargetIdsInCurrentViewerZEnvelope': cached_current_envelope,
        'coveredUniqueRowsZEnvelope': int(cached_current_envelope + len(new_extension_ids)),
        'coverageFractionZEnvelope': float((cached_current_envelope + len(new_extension_ids)) / envelope_rows),
        'selectionAlignedExactBaselineRowsFromPriorAudit': int(prior_audit['totals']['exactCacheMatches']),
        'selectionAlignedBaselinePlusAllExtensionUniqueRows': int(prior_audit['totals']['exactCacheMatches'] + len(unique_extension_ids)),
        'selectionAlignedCoverageFractionTotalCatalog': float((prior_audit['totals']['exactCacheMatches'] + len(unique_extension_ids)) / len(viewer_unique)),
        'selectionAlignedCoverageFractionZEnvelope': float((prior_audit['totals']['exactCacheMatches'] + len(unique_extension_ids)) / envelope_rows),
        'selectionAlignedBaselinePlusGloballyNovelExtensions': int(prior_audit['totals']['exactCacheMatches'] + len(new_extension_ids)),
        'selectionAlignedCoverageFractionTotalCatalogWithCrossCacheOverlapsExcluded': float((prior_audit['totals']['exactCacheMatches'] + len(new_extension_ids)) / len(viewer_unique)),
        'selectionAlignedCoverageFractionZEnvelopeWithCrossCacheOverlapsExcluded': float((prior_audit['totals']['exactCacheMatches'] + len(new_extension_ids)) / envelope_rows),
        'definition': 'reports both selection-aligned coverage and any-cache object union; foreign-tracer cache hits are not assumed scientifically interchangeable',
    }
    all_checks = [v for f in report['fields'] for v in f['checks'].values()]
    all_tracer_checks = [v for tr in report.get('tracerChecks', {}).values()
                         for k, v in tr.items() if isinstance(v, bool)]
    global_checks = [report['globalCoverage']['extensionIdsAllBelongToCurrentViewer'],
                     report['globalCoverage']['duplicateExtensionRowsAcrossAnyTracer'] == report['globalCoverage']['crossTracerDuplicateTargetIds'],
                     report['globalCoverage']['crossTracerDuplicateTargetIds'] == 0]
    report['status'] = 'complete' if report['manifestStatus'] == 'complete' and all(all_checks + all_tracer_checks + global_checks) else 'failed'
    (EXT / 'final-extension-audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'fields': len(report['fields']),
                      'totalRows': report['totalRows'], 'tracerChecks': report.get('tracerChecks'),
                      'globalCoverage': report['globalCoverage'],
                      'output': str(EXT / 'final-extension-audit.json')}, indent=2))
    return 0 if report['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())
