#!/usr/bin/env python3
"""Inventory local RSD extension assets for later scientific review and upload."""
import hashlib
import gzip
import json
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / 'ply/reconstruction_dr1/staged_extensions_v1'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, default=STAGE)
    parser.add_argument('--expected-rsd-matched', type=int, default=13_462_756)
    args = parser.parse_args()
    stage = args.stage.resolve()
    coverage = json.loads((stage / 'payloads/coverage-report.json').read_text())
    lod = json.loads((stage / 'lod_pairs/manifest.json').read_text())
    extension = ROOT / 'tools/rsd/extensions/extension-manifest.json'
    if coverage.get('scientificMaskStatus') != 'unvalidated-review-only':
        raise ValueError('expected review-only staged payloads')
    if coverage['totals']['rsdMatched'] != args.expected_rsd_matched:
        raise ValueError('unexpected RSD matched count')
    if sum(x['rsdNonzero'] for x in lod['catalogs'].values()) != coverage['totals']['rsdNonzero']:
        raise ValueError('LOD and point RSD totals differ')
    files = []
    for i in range(12):
        path = stage / 'payloads' / f'rsd_{i:03d}.bin'
        files.append({'role': 'point-rsd', 'chunk': i, 'path': str(path.relative_to(ROOT)),
                      'bytes': path.stat().st_size, 'sha256': sha256(path)})
        baseline_fog = ROOT / 'ply/reconstruction_dr1/payloads' / f'fog_{i:03d}.bin'
        staged_fog = stage / 'payloads' / f'fog_{i:03d}.bin'
        if gzip.decompress(staged_fog.read_bytes()) != gzip.decompress(baseline_fog.read_bytes()):
            raise ValueError(f'chunk {i}: staged FoG differs from baseline')
    for name in ('DR1_GALAXY', 'DR1_QSO'):
        path = stage / 'lod_pairs' / f'correction_{name}.bin'
        files.append({'role': 'lod-pair', 'catalog': name, 'path': str(path.relative_to(ROOT)),
                      'bytes': path.stat().st_size, 'sha256': sha256(path)})
    report = {'status': 'staged-for-review', 'readyForUpload': False,
              'reason': 'independent angular/radial mask and scientific field validation pending',
              'extensionManifestSha256': sha256(extension),
              'bgsLowzFieldManifestSha256': (sha256(ROOT / coverage['bgsLowzCandidate'] / 'field-manifest.json')
                                                if coverage.get('bgsLowzCandidate') else None),
              'rsdMatched': coverage['totals']['rsdMatched'],
              'fogUnchangedFromBaseline': True,
              'uploadFiles': files,
              'uploadSequence': 'after approval: upload 12 RSD chunk files and 2 LOD pairs; create new overlay with returned asset IDs; parse and test Launch before publication'}
    (stage / 'asset-staging-manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'files': len(files),
                      'totalBytes': sum(x['bytes'] for x in files),
                      'rsdMatched': report['rsdMatched']}))


if __name__ == '__main__':
    main()
