#!/usr/bin/env python3
"""Verify all twelve deployed-format catalogue chunks without loading them together."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import struct
from pathlib import Path


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path('desiV3/data/full-photometry'))
    parser.add_argument('--manifest', type=Path, help='final manifest; defaults to DIRECTORY/manifest.json')
    parser.add_argument('--reference', type=Path, default=Path('desiV3/reference-current.sha256'))
    parser.add_argument('--require-reference-match', action='store_true',
                        help='also require byte-exact SHA256 matches to the historically deployed chunks')
    args = parser.parse_args()
    manifest = json.loads((args.manifest or args.directory / 'manifest.json').read_text())
    validation = json.loads((args.directory / 'validation.json').read_text())
    expected = {Path(line.split(None, 1)[1]).name: line.split()[0]
                for line in args.reference.read_text().splitlines() if line.strip()}
    if len(manifest['chunks']) != len(validation['chunks']) != 12:
        raise ValueError('expected twelve manifest and validation rows')
    if manifest['count'] != 15786217 or sum(row['count'] for row in manifest['chunks']) != manifest['count']:
        raise ValueError('catalogue row count differs from the recorded DR1 selection')
    total_bytes = 0
    for number, spec in enumerate(manifest['chunks']):
        filename = f'catalog_{number:03d}.bin.gz'
        if spec['file'] != filename or validation['chunks'][number]['chunk'] != number:
            raise ValueError(f'chunk {number}: unexpected order')
        path = args.directory / filename
        row = validation['chunks'][number]
        if path.stat().st_size != row['bytes'] or row['count'] != spec['count'] or spec['bytes'] != row['bytes']:
            raise ValueError(f'{filename}: size or row count differs from validation')
        actual = digest(path)
        if actual != row['sha256']:
            raise ValueError(f'{filename}: SHA256 differs from local validation')
        if args.require_reference_match and actual != expected[filename]:
            raise ValueError(f'{filename}: SHA256 differs from the historical release')
        with gzip.open(path, 'rb') as stream:
            header = stream.read(36)
        if header[:4] != b'DSC3':
            raise ValueError(f'{filename}: invalid DSC3 signature')
        version, count, *offsets, length = struct.unpack_from('<8I', header, 4)
        if version != 1 or count != spec['count'] or len(offsets) != 5:
            raise ValueError(f'{filename}: invalid header')
        sections = [count, count * 4, count * 6, count * 2, count * 4]
        if any(offset % 4 for offset in offsets) or offsets[0] != 64:
            raise ValueError(f'{filename}: invalid section alignment')
        if any(a + size > b for a, size, b in zip(offsets, sections, offsets[1:] + [length])):
            raise ValueError(f'{filename}: overlapping sections')
        total_bytes += path.stat().st_size
        print(f'{filename}: {count:,} rows, {path.stat().st_size:,} bytes, SHA256 OK')
    print(f'PASS: {manifest["count"]:,} rows in 12 chunks; {total_bytes:,} compressed bytes')


if __name__ == '__main__':
    main()
