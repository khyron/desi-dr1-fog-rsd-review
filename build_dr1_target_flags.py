#!/usr/bin/env python3
"""Recover DESI target-selection labels without changing DR1 catalog row IDs.

Each byte is aligned with one row of DR1_GALAXY.desi or DR1_QSO.desi:
bit 0 BGS, bit 1 BGS_BRIGHT, bit 2 LRG, bit 3 ELG,
bit 4 QSO target, bit 5 other (none of the preceding target classes).
The spectral class remains the catalog filename; targeting labels may overlap.
"""

import json
import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits

import build_dr1_full_catalog as catalog
from build_dr1_full_catalog import ROW_BATCH, selected

DESI_COLUMNS = ('DESI_TARGET', 'SV1_DESI_TARGET', 'SV2_DESI_TARGET', 'SV3_DESI_TARGET')
BGS_COLUMNS = ('BGS_TARGET', 'SV1_BGS_TARGET', 'SV2_BGS_TARGET', 'SV3_BGS_TARGET')
LABELS = ('BGS', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO_TARGET', 'OTHER')


def flags_for_rows(rows):
    desi = np.zeros(len(rows), dtype=np.uint64)
    bgs = np.zeros(len(rows), dtype=np.uint64)
    for column in DESI_COLUMNS:
        desi |= np.asarray(rows[column], dtype=np.uint64)
    for column in BGS_COLUMNS:
        bgs |= np.asarray(rows[column], dtype=np.uint64)
    flags = np.zeros(len(rows), dtype=np.uint8)
    flags |= ((bgs & 15) != 0).astype(np.uint8)
    flags |= (((bgs & 2) != 0).astype(np.uint8) << 1)
    flags |= (((desi & 1) != 0).astype(np.uint8) << 2)
    flags |= (((desi & 2) != 0).astype(np.uint8) << 3)
    flags |= (((desi & 4) != 0).astype(np.uint8) << 4)
    flags |= ((flags == 0).astype(np.uint8) << 5)
    return flags


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=catalog.SOURCE)
    parser.add_argument('--dest', type=Path, default=catalog.DEST)
    parser.add_argument('--min-reliable-z', type=float, default=catalog.MIN_RELIABLE_Z)
    args = parser.parse_args()
    catalog.MIN_RELIABLE_Z = args.min_reliable_z
    args.dest.mkdir(parents=True, exist_ok=True)
    summary = {}
    with fits.open(args.source, memmap=True) as hdul:
        data = hdul[1].data
        required = set(DESI_COLUMNS + BGS_COLUMNS)
        missing = required.difference(data.names)
        if missing:
            raise ValueError(f'missing targeting columns: {sorted(missing)}')
        outputs = {}
        for kind in (b'GALAXY', b'QSO'):
            name = f'DR1_{kind.decode()}'
            catalogue_path = args.dest / f'{name}.desi'
            with catalogue_path.open('rb') as stream:
                header = stream.read(24)
            version = int(np.frombuffer(header, dtype='<u4', count=1, offset=4)[0])
            count = int(np.frombuffer(header, dtype='<u4', count=1, offset=8)[0])
            if version not in (2, 3):
                raise ValueError(f'{catalogue_path}: expected v2 or v3 catalog')
            z_min, z_max = np.frombuffer(header, dtype='<f4', count=2, offset=16)
            outputs[kind] = {
                'name': name, 'count': count, 'written': 0,
                'ids': np.memmap(catalogue_path, dtype='<u8', mode='r', offset=48 + count * 24,
                                 shape=(count,)) if version == 3 else None,
                'redshift': np.memmap(catalogue_path, dtype='<u2', mode='r', offset=48 + count * 6,
                                      shape=(count,)) if version == 2 else None,
                'z_min': float(z_min), 'z_max': float(z_max),
                'flags': np.memmap(args.dest / f'{name}.target-flags.bin', dtype='u1',
                                   mode='w+', shape=(count,)),
                'labels': np.zeros(len(LABELS), dtype=np.int64),
            }
        for start in range(0, len(data), ROW_BATCH):
            end = min(start + ROW_BATCH, len(data))
            for kind, output in outputs.items():
                rows, z, _, _ = selected(data, start, end, kind)
                n = len(rows)
                if not n:
                    continue
                sl = slice(output['written'], output['written'] + n)
                if output['ids'] is not None:
                    if not np.array_equal(output['ids'][sl], np.asarray(rows['TARGETID'], dtype='<u8')):
                        raise ValueError(f"{output['name']}: row order differs from the .desi catalog")
                else:
                    span = output['z_max'] - output['z_min']
                    expected_z = (np.rint((z - output['z_min']) / span * 65535)
                                  .astype('<u2') if span else np.zeros(n, dtype='<u2'))
                    if np.any(np.abs(output['redshift'][sl].astype(np.int32) -
                                     expected_z.astype(np.int32)) > 1):
                        raise ValueError(f"{output['name']}: redshift row order differs")
                flags = flags_for_rows(rows)
                output['flags'][sl] = flags
                output['labels'] += np.array([np.count_nonzero(flags & (1 << bit))
                                               for bit in range(len(LABELS))])
                output['written'] += n
            print(f'\rTarget flags {end:,}/{len(data):,}', end='', flush=True)
    print()
    for output in outputs.values():
        if output['written'] != output['count']:
            raise ValueError(f"{output['name']}: expected {output['count']}, got {output['written']}")
        output['flags'].flush()
        summary[output['name']] = {
            'count': output['count'],
            'targetCounts': dict(zip(LABELS, map(int, output['labels']))),
        }
    (args.dest / 'dr1_target_flags_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
