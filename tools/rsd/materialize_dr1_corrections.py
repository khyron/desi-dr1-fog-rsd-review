#!/usr/bin/env python3
"""Map previously reconstructed DR1 Gfinder/RSD TARGETIDs to exact DR1 chunks.

Unmatched objects retain their observed position (zero delta). The output report
states coverage explicitly; it does not claim that an LSS reconstruction field
extends beyond its source selection or redshift interval.
"""
import argparse
import gzip
import hashlib
import json
import struct
from pathlib import Path

import numpy as np

from build_dr1_chunk_index import header, morton_order

ROOT = Path(__file__).resolve().parents[2]
FIELDS = [('BGS', 0.1, 0.4), ('LRG', 0.4, 0.8),
          ('LRG_ELG', 0.8, 1.1), ('ELG', 1.1, 1.6), ('QSO', 1.6, 2.1)]


def lookup(ids, values, query):
    pos = np.searchsorted(ids, query)
    found = (pos < len(ids)) & (ids[np.minimum(pos, len(ids) - 1)] == query)
    result = np.zeros(len(query), dtype='<f4')
    result[found] = values[pos[found]]
    return result, found


def load_rsd(cache_dir):
    output = {}
    for field, _, _ in FIELDS:
        files = [cache_dir / f'{field}_{region}_rsd_query_v1.npz' for region in ('NGC', 'SGC')]
        files = [path for path in files if path.is_file()]
        if len(files) != 2:
            raise ValueError(f'{field}: expected NGC and SGC caches, found {len(files)}')
        sources = []
        for path in files:
            with np.load(path) as data:
                sources.append((np.asarray(data['targetid'], dtype='<u8'),
                                np.asarray(data['delta_mpc'], dtype='<f4')))
        ids = np.concatenate([x[0] for x in sources])
        delta = np.concatenate([x[1] for x in sources])
        order = np.argsort(ids, kind='stable')
        ids, delta = ids[order], delta[order]
        if np.any(ids[1:] == ids[:-1]) or not np.isfinite(delta).all():
            raise ValueError(f'{field}: duplicate TARGETID or nonfinite correction')
        output[field] = (ids, delta)
    return output


def add_rsd_extensions(rsd, extension_dir):
    """Merge only within the assigned tracer; foreign cache hits never take precedence."""
    manifest = json.loads((extension_dir / 'extension-manifest.json').read_text())
    if manifest.get('status') != 'complete' or len(manifest.get('fields', [])) != 10:
        raise ValueError('complete ten-field extension manifest required')
    grouped = {field: [] for field, _, _ in FIELDS}
    seen = set()
    for entry in manifest['fields']:
        field, region = entry['tracer'], entry['region']
        if field not in grouped or region not in ('NGC', 'SGC') or (field, region) in seen:
            raise ValueError(f'invalid or duplicate extension field: {field}/{region}')
        seen.add((field, region))
        if entry.get('status') != 'complete':
            raise ValueError(f'incomplete extension field: {field}/{region}')
        path = extension_dir / f'{field}_{region}_uncached_extension_v1.npz'
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != entry.get('outputSha256'):
            raise ValueError(f'extension hash mismatch: {path}')
        with np.load(path, allow_pickle=False) as data:
            ids = np.asarray(data['targetid'], dtype='<u8')
            values = np.asarray(data['delta_mpc'], dtype='<f4')
        if len(ids) != len(values) or len(ids) != entry.get('extensionRows') or not np.isfinite(values).all():
            raise ValueError(f'invalid extension rows: {path}')
        grouped[field].append((ids, values))
    if len(seen) != 10:
        raise ValueError('missing extension fields')
    for field, parts in grouped.items():
        base_ids, base_values = rsd[field]
        ids = np.concatenate([base_ids] + [p[0] for p in parts])
        values = np.concatenate([base_values] + [p[1] for p in parts])
        order = np.argsort(ids, kind='stable')
        ids, values = ids[order], values[order]
        if np.any(ids[1:] == ids[:-1]):
            raise ValueError(f'{field}: extension overlaps same-tracer cache or region')
        rsd[field] = ids, values
    return rsd


def redshift_by_chunk(assets, manifest, index):
    units = float(manifest['unitsPerMpc'])
    positions, redshifts = [], []
    for name in manifest['tracers']:
        path = assets / (name + '.desi')
        count, own_units, zmin, zmax = header(path)
        pos = np.memmap(path, dtype='<i2', mode='r', offset=48, shape=(count, 3))
        zq = np.memmap(path, dtype='<u2', mode='r', offset=48 + count * 6, shape=(count,))
        positions.append(np.rint(pos.astype('f4') / own_units * units).astype('<i2'))
        redshifts.append(zmin + zq.astype('f4') * ((zmax - zmin) / 65535.0))
    pos = np.concatenate(positions)
    z = np.concatenate(redshifts)
    hashes = np.abs(pos[:, 0].astype('i8') * 17 + pos[:, 2].astype('i8') * 59 + pos[:, 1].astype('i8') * 101)
    for i, spec in enumerate(manifest['chunks']):
        chosen = np.flatnonzero(hashes % len(manifest['chunks']) == i)
        chosen = chosen[morton_order(pos[chosen])]
        if len(chosen) != spec['count'] or len(np.load(index / f'targetid_{i:03d}.npy', mmap_mode='r')) != len(chosen):
            raise ValueError(f'chunk {i}: index/geometry count mismatch')
        yield i, z[chosen]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chunks', type=Path, default=ROOT / 'ply' / 'chunks')
    parser.add_argument('--assets', type=Path, default=ROOT / 'assets')
    parser.add_argument('--index', type=Path, default=ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    parser.add_argument('--fog', type=Path, default=Path(r'D:\desi\data_external\fog\analysis\targetid-delta.npy'))
    parser.add_argument('--rsd-cache', type=Path, default=Path(r'D:\desi\data_external\RSD\derived\rsd-v1-exploratory'))
    parser.add_argument('--rsd-extensions', type=Path, help='audited local extensions; requires a distinct explicit output directory')
    parser.add_argument('--bgs-lowz-candidate', type=Path, help='review-only BGS_BRIGHT 0.01<=z<0.40 field directory')
    parser.add_argument('--mask-approval-report', type=Path, help='scientific mask approval JSON for extension export')
    parser.add_argument('--stage-unvalidated-rsd', action='store_true', help='build local review assets before independent mask validation')
    parser.add_argument('--out', type=Path, default=ROOT / 'ply' / 'reconstruction_dr1' / 'payloads')
    args = parser.parse_args()
    if args.stage_unvalidated_rsd and not args.rsd_extensions:
        raise ValueError('--stage-unvalidated-rsd requires --rsd-extensions')
    if args.bgs_lowz_candidate and not args.stage_unvalidated_rsd:
        raise ValueError('low-z BGS candidate requires --stage-unvalidated-rsd')
    if args.rsd_extensions:
        if not args.stage_unvalidated_rsd:
            if not args.mask_approval_report:
                raise ValueError('extensions require --mask-approval-report or --stage-unvalidated-rsd')
            approval = json.loads(args.mask_approval_report.read_text())
            extension_manifest = args.rsd_extensions / 'extension-manifest.json'
            expected_hash = hashlib.sha256(extension_manifest.read_bytes()).hexdigest()
            if approval.get('status') != 'approved' or approval.get('extensionManifestSha256') != expected_hash:
                raise ValueError('mask approval must bind to the exact extension manifest')
    if args.rsd_extensions and (args.out == parser.get_default('out') or args.out.resolve() == Path(parser.get_default('out')).resolve()):
        raise ValueError('extensions require --out distinct from published baseline payloads')
    if args.rsd_extensions and args.out.exists() and any(args.out.iterdir()):
        raise ValueError('extension output directory must be empty')
    manifest = json.loads((args.chunks / 'manifest.json').read_text())
    exact = json.loads((args.index / 'index-report.json').read_text())
    if exact['count'] != manifest['count'] or not all(c['geometryExact'] for c in exact['chunks']):
        raise ValueError('exact geometry index validation is required')
    fog_source = np.load(args.fog, mmap_mode='r')
    order = np.argsort(fog_source['targetid'], kind='stable')
    fog_ids = np.asarray(fog_source['targetid'][order], dtype='<u8')
    fog_values = np.asarray(fog_source['delta_mpc'][order], dtype='<f4')
    if np.any(fog_ids[1:] == fog_ids[:-1]) or not np.isfinite(fog_values).all():
        raise ValueError('invalid FoG source')
    rsd = load_rsd(args.rsd_cache)
    if args.rsd_extensions:
        rsd = add_rsd_extensions(rsd, args.rsd_extensions)
    lowz_bgs = None
    if args.bgs_lowz_candidate:
        candidate_manifest = json.loads((args.bgs_lowz_candidate / 'field-manifest.json').read_text())
        if (candidate_manifest.get('status') != 'complete-for-technical-review' or
                candidate_manifest.get('redshiftRange') != [0.01, 0.4] or
                len(candidate_manifest.get('regions', [])) != 2):
            raise ValueError('complete two-region low-z BGS candidate required')
        pieces = []
        for entry in candidate_manifest['regions']:
            region = entry.get('region')
            if region not in ('NGC', 'SGC'):
                raise ValueError('invalid low-z BGS region')
            path = args.bgs_lowz_candidate / f'BGS_BRIGHT_{region}_viewer_rsd_v1.npz'
            if hashlib.sha256(path.read_bytes()).hexdigest() != entry.get('outputSha256'):
                raise ValueError(f'low-z BGS hash mismatch: {path}')
            with np.load(path, allow_pickle=False) as data:
                pieces.append((np.asarray(data['targetid'], dtype='<u8'),
                               np.asarray(data['delta_mpc'], dtype='<f4')))
        if {entry['region'] for entry in candidate_manifest['regions']} != {'NGC', 'SGC'}:
            raise ValueError('NGC and SGC low-z BGS candidates required')
        ids = np.concatenate([part[0] for part in pieces])
        values = np.concatenate([part[1] for part in pieces])
        order = np.argsort(ids, kind='stable')
        ids, values = ids[order], values[order]
        if np.any(ids[1:] == ids[:-1]) or not np.isfinite(values).all():
            raise ValueError('invalid low-z BGS candidate ID/value arrays')
        lowz_bgs = ids, values
    args.out.mkdir(parents=True, exist_ok=True)
    totals = {'rows': 0, 'fogMatched': 0, 'fogNonzero': 0, 'rsdEligible': 0,
              'rsdMatched': 0, 'rsdNonzero': 0}
    rows = []
    for i, z in redshift_by_chunk(args.assets, manifest, args.index):
        ids = np.load(args.index / f'targetid_{i:03d}.npy', mmap_mode='r')
        fog, fog_hit = lookup(fog_ids, fog_values, ids)
        correction = np.zeros(len(ids), dtype='<f4')
        matched = np.zeros(len(ids), dtype=bool)
        eligible = np.zeros(len(ids), dtype=bool)
        for field, low, high in FIELDS:
            mask = (z >= low) & (z < high)
            eligible |= mask
            if not np.any(mask):
                continue
            values, hit = lookup(*rsd[field], ids[mask])
            correction[mask] = values
            matched[mask] = hit
        if lowz_bgs is not None:
            # One BGS_BRIGHT field spans both sides of z=0.10. It replaces the
            # earlier BGS_BRIGHT-21.5 deltas over the entire BGS range, avoiding
            # a discontinuous switch of fields at that boundary.
            mask = (z >= 0.01) & (z < 0.4)
            eligible |= mask
            if np.any(mask):
                values, hit = lookup(*lowz_bgs, ids[mask])
                correction[mask] = values
                matched[mask] = hit
        if not np.isfinite(fog).all() or not np.isfinite(correction).all():
            raise ValueError(f'chunk {i}: nonfinite correction')
        if np.max(np.abs(fog)) > 300.001 or np.max(np.abs(correction)) > 150:
            raise ValueError(f'chunk {i}: implausible radial displacement')
        for label, values in [('fog', fog), ('rsd', correction)]:
            (args.out / f'{label}_{i:03d}.bin').write_bytes(gzip.compress(values.tobytes(), compresslevel=9))
        item = {'chunk': i, 'rows': len(ids), 'fogMatched': int(fog_hit.sum()),
                'fogNonzero': int(np.count_nonzero(fog)), 'rsdEligible': int(eligible.sum()),
                'rsdMatched': int(matched.sum()), 'rsdNonzero': int(np.count_nonzero(correction))}
        for key in totals:
            totals[key] += item[key]
        rows.append(item)
        print(json.dumps(item), flush=True)
    if totals['rows'] != manifest['count'] or totals['rsdMatched'] < manifest['count'] * 0.1:
        raise ValueError('count or RSD coverage validation failed')
    report = {'sourceGeometry': 'exact-validated-published-DR1-chunks',
              'fogSource': 'DESI DR1 Gfinder v1.0 group analysis',
              'rsdSource': ('DESI DR1 LSS v1.5 BGS_BRIGHT 0.01<=z<0.40 candidate plus other exploratory fields'
                            if lowz_bgs is not None else 'DESI DR1 LSS iron LSScats v1.5 exploratory 16 h^-1 Mpc field'),
              'rsdExtensions': str(args.rsd_extensions) if args.rsd_extensions else None,
              'bgsLowzCandidate': str(args.bgs_lowz_candidate) if args.bgs_lowz_candidate else None,
              'scientificMaskStatus': ('unvalidated-review-only' if args.stage_unvalidated_rsd else
                                       'approved' if args.rsd_extensions else 'baseline'),
              'rsdSelectionPolicy': 'half-open redshift interval assigns tracer; same-tracer extension fills uncached IDs; foreign-tracer caches ignored',
              'unmatchedPolicy': 'zero radial delta; retain observed position',
              'totals': totals, 'rsdMatchedFraction': totals['rsdMatched'] / totals['rows'],
              'rsdEligibleMatchedFraction': totals['rsdMatched'] / totals['rsdEligible'],
              'chunks': rows}
    (args.out / 'coverage-report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['totals']), flush=True)


if __name__ == '__main__':
    main()
