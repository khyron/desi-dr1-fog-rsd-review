#!/usr/bin/env python3
"""Copy validated final chunks to the twelve current asset filenames."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path):
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('desiV3/data/full-photometry'))
    parser.add_argument('--out', type=Path, default=Path('desiV3/data/upload'))
    args = parser.parse_args()
    manifest = json.loads((args.source / 'manifest.json').read_text())
    validation = json.loads((args.source / 'validation.json').read_text())
    if len(manifest['chunks']) != 12 or len(validation['chunks']) != 12:
        raise ValueError('expected twelve final chunks')
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = []
    for number, spec in enumerate(manifest['chunks']):
        source = args.source / f'catalog_{number:03d}.bin.gz'
        name = f'desiV3_catalog_{number:03d}.bin.gz'
        destination = args.out / name
        expected = validation['chunks'][number]
        if spec['file'] != source.name or expected['chunk'] != number:
            raise ValueError(f'chunk {number}: unexpected ordering')
        if sha256(source) != expected['sha256'] or source.stat().st_size != expected['bytes']:
            raise ValueError(f'{source}: failed final-source validation')
        if not destination.exists() or sha256(destination) != expected['sha256']:
            shutil.copyfile(source, destination)
        if sha256(destination) != expected['sha256']:
            raise ValueError(f'{destination}: copy verification failed')
        spec['file'] = name
        receipt.append({'file': name, 'bytes': destination.stat().st_size, 'sha256': expected['sha256']})
    (args.out / 'desiV3-consolidated-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    (args.out / 'staging-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'Staged twelve named assets: {sum(row["bytes"] for row in receipt):,} bytes')


if __name__ == '__main__':
    main()
