#!/usr/bin/env python3
"""Build a PlayCanvas review overlay from uploaded extension asset IDs."""
import json
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / 'ply/reconstruction_dr1/staged_extensions_v1'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, default=STAGE)
    parser.add_argument('--source-label', default='uncached extension v1, visual review only')
    parser.add_argument('--output-name', default='reconstruction-dr1-overlay-ext-review.json')
    args = parser.parse_args()
    stage = args.stage.resolve()
    ids = json.loads((stage / 'playcanvas-upload-ids.json').read_text())
    staged = json.loads((stage / 'payloads/coverage-report.json').read_text())
    base = json.loads((ROOT / 'ply/reconstruction_dr1/reconstruction-dr1-overlay.json').read_text())
    if ids['projectId'] != 1605307 or len(ids['rsdChunkAssetIds']) != len(base['chunks']) or len(base['chunks']) != 12:
        raise ValueError('project or chunk count mismatch')
    if staged['scientificMaskStatus'] != 'unvalidated-review-only':
        raise ValueError('unexpected staged extension status or count')
    base['coverage'] = staged['totals']
    base['rsd']['nonzero'] = staged['totals']['rsdNonzero']
    base['rsd']['source'] += '; ' + args.source_label
    base['scientificMaskStatus'] = 'unvalidated-review-only'
    for i, chunk in enumerate(base['chunks']):
        path = stage / 'payloads' / f'rsd_{i:03d}.bin'
        chunk['rsd'] = f"asset:{ids['rsdChunkAssetIds'][i]}"
        chunk['rsdBytes'] = path.stat().st_size
    out = stage / args.output_name
    out.write_text(json.dumps(base, indent=2) + '\n')
    print(json.dumps({'file': str(out), 'chunks': len(base['chunks']),
                      'rsdMatched': base['coverage']['rsdMatched']}))


if __name__ == '__main__':
    main()
