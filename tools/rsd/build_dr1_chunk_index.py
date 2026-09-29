#!/usr/bin/env python3
"""Recover exact TARGETID order for the published DR1 point chunks.

The current DR1 .desi v2 files omit TARGETID. Recover identifiers in their
original FITS order, verify the stored redshift row by row, then replay the
published point partition and compare every encoded geometry chunk exactly.
The private .npy index is never a browser asset.
"""
import argparse
import gzip
import json
import struct
from pathlib import Path

import numpy as np
from astropy.io import fits

import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import build_dr1_full_catalog as catalog
from build_desi_chunks import encode_chunk, morton_order


def header(path):
    with path.open('rb') as stream:
        raw = stream.read(24)
    if raw[:4] != b'DESI':
        raise ValueError(f'{path}: invalid magic')
    version, count, units, zmin, zmax = struct.unpack_from('<IIfff', raw, 4)
    if version != 2:
        raise ValueError(f'{path}: expected published v2 catalog, found v{version}')
    return count, units, zmin, zmax


def recover_ids(source, assets, out, counts):
    catalog.MIN_RELIABLE_Z = 0.0  # Published v2 contains the original z>0 selection.
    ids = {}
    redshifts = {}
    positions = {}
    z_ranges = {}
    for kind in catalog.CLASSES:
        name = 'DR1_' + kind.decode()
        path = assets / (name + '.desi')
        count, units, zmin, zmax = header(path)
        if count != counts[name]:
            raise ValueError(f'{name}: {count} catalog rows vs {counts[name]} manifest rows')
        ids[kind] = np.lib.format.open_memmap(out / (name + '_targetid.npy'), mode='w+', dtype='<u8', shape=(count,))
        redshifts[kind] = np.memmap(path, dtype='<u2', mode='r', offset=48 + count * 6, shape=(count,))
        positions[kind] = (np.memmap(path, dtype='<i2', mode='r', offset=48, shape=(count, 3)), units)
        z_ranges[kind] = (zmin, zmax)
    offsets = {kind: 0 for kind in catalog.CLASSES}
    with fits.open(source, memmap=True) as hdul:
        data = hdul[1].data
        for start in range(0, len(data), catalog.ROW_BATCH):
            for kind in catalog.CLASSES:
                rows, z, _, _ = catalog.selected(data, start, min(start + catalog.ROW_BATCH, len(data)), kind)
                if not len(rows):
                    continue
                at = offsets[kind]
                end = at + len(rows)
                if end > len(ids[kind]):
                    raise ValueError(f'{kind.decode()}: more FITS rows than published')
                count = len(ids[kind])
                name = 'DR1_' + kind.decode()
                _, _, zmin, zmax = header(assets / (name + '.desi'))
                expected = np.rint((z - zmin) / (zmax - zmin) * 65535).astype('<u2')
                if np.any(np.abs(redshifts[kind][at:end].astype('i4') - expected.astype('i4')) > 1):
                    raise ValueError(f'{name}: redshift order differs at FITS rows {start}:{start + catalog.ROW_BATCH}')
                ids[kind][at:end] = np.asarray(rows['TARGETID'], dtype='<u8')
                offsets[kind] = end
            if start % (catalog.ROW_BATCH * 20) == 0:
                print(f'FITS {start:,}/{len(data):,}', flush=True)
    for kind in catalog.CLASSES:
        if offsets[kind] != len(ids[kind]):
            raise ValueError(f'{kind.decode()}: recovered {offsets[kind]} vs {len(ids[kind])}')
        ids[kind].flush()
    return ids, positions, redshifts, z_ranges


def make_index(ids, positions, redshifts, z_ranges, manifest, chunks, out):
    units = float(manifest['unitsPerMpc'])
    pos = np.concatenate([np.rint(src.astype('f4') / own_units * units).astype('<i2')
                          for src, own_units in (positions[kind] for kind in catalog.CLASSES)])
    target = np.concatenate([np.asarray(ids[kind]) for kind in catalog.CLASSES])
    z = np.concatenate([z_ranges[kind][0] + np.asarray(redshifts[kind], dtype='f4') *
                        ((z_ranges[kind][1] - z_ranges[kind][0]) / 65535.0)
                        for kind in catalog.CLASSES])
    if len(pos) != manifest['count']:
        raise ValueError('position count differs from manifest')
    hashes = np.abs(pos[:, 0].astype('i8') * 17 + pos[:, 2].astype('i8') * 59 + pos[:, 1].astype('i8') * 101)
    report = []
    for i, spec in enumerate(manifest['chunks']):
        chosen = np.flatnonzero(hashes % len(manifest['chunks']) == i)
        chosen = chosen[morton_order(pos[chosen])]
        if len(chosen) != spec['count']:
            raise ValueError(f'chunk {i}: row count mismatch')
        published = gzip.decompress((chunks / spec['file']).read_bytes())
        candidate = gzip.decompress(encode_chunk(pos[chosen]))
        if not np.array_equal(np.frombuffer(published, 'u1'), np.frombuffer(candidate, 'u1')):
            raise ValueError(f'chunk {i}: geometry/order differs from published payload')
        path = out / f'targetid_{i:03d}.npy'
        np.save(path, target[chosen])
        np.save(out / f'redshift_{i:03d}.npy', z[chosen])
        report.append({'chunk': i, 'rows': int(len(chosen)), 'geometryExact': True, 'index': path.name})
        print(f'chunk {i}: {len(chosen):,} exact rows', flush=True)
    (out / 'index-report.json').write_text(json.dumps({'count': len(pos), 'chunks': report}, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'data' / 'zall-pix-iron.fits')
    parser.add_argument('--assets', type=Path, default=ROOT / 'assets')
    parser.add_argument('--chunks', type=Path, default=ROOT / 'ply' / 'chunks')
    parser.add_argument('--out', type=Path, default=ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.chunks / 'manifest.json').read_text())
    ids, positions, redshifts, z_ranges = recover_ids(args.source, args.assets, args.out, manifest['tracerCounts'])
    make_index(ids, positions, redshifts, z_ranges, manifest, args.chunks, args.out)


if __name__ == '__main__':
    main()
