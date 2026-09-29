#!/usr/bin/env python3
"""Write a local final manifest with actual compressed catalogue byte sizes."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('desiV3/data/manifest.json'))
    parser.add_argument('--chunks', type=Path, default=Path('desiV3/data/full-photometry'))
    parser.add_argument('--out', type=Path, default=Path('desiV3/data/full-photometry/manifest.json'))
    args = parser.parse_args()
    manifest = json.loads(args.base.read_text())
    validation = json.loads((args.chunks / 'validation.json').read_text())
    if len(manifest['chunks']) != 12 or len(validation['chunks']) != 12:
        raise ValueError('expected twelve chunks')
    for number, spec in enumerate(manifest['chunks']):
        name = f'catalog_{number:03d}.bin.gz'
        row = validation['chunks'][number]
        actual = (args.chunks / name).stat().st_size
        if spec['file'] != name or row['chunk'] != number or row['count'] != spec['count'] or row['bytes'] != actual:
            raise ValueError(f'chunk {number}: manifest, validation or file mismatch')
        spec['bytes'] = actual
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Wrote {args.out}: {manifest["count"]:,} rows, '
          f'{sum(row["bytes"] for row in manifest["chunks"]):,} compressed bytes')


if __name__ == '__main__':
    main()
