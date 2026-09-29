#!/usr/bin/env python3
"""Check the four DR1 v2/targeting inputs against recorded production hashes."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets', type=Path, default=Path('assets'))
    parser.add_argument('--reference', type=Path, default=Path('desiV3/reference-inputs.json'))
    args = parser.parse_args()
    expected = json.loads(args.reference.read_text())
    for name, wanted in expected.items():
        path = args.assets / name
        sha = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                sha.update(block)
        if sha.hexdigest() != wanted:
            raise ValueError(f'{path}: SHA256 differs from historical input')
        print(f'{name}: exact SHA256 match')
    print('PASS: four exact historical inputs')


if __name__ == '__main__':
    main()
