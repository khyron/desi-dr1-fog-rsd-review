#!/usr/bin/env python3
"""Build 3D-ready DR1 galaxy and quasar catalogs from the official zall-pix v1 FITS.

Only primary science targets with ZWARN=0 and positive finite redshift are used.
Stars are excluded because their spectroscopic redshift is not a distance.
The outputs use the existing .desi v3 layout; unavailable imaging fields are
marked unknown rather than invented.
"""

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from astropy.io import fits

SOURCE = Path(__file__).parent / 'data' / 'zall-pix-iron.fits'
DEST = Path(__file__).parent / 'assets'
CLASSES = (b'GALAXY', b'QSO')
ROW_BATCH = 250_000
MIN_RELIABLE_Z = 0.0045  # ≈20 Mpc in Planck18; nearer objects need independent distances.


def selected(data, start, end, kind):
    rows = data[start:end]
    z = np.asarray(rows['Z'], dtype=np.float64)
    ra = np.asarray(rows['TARGET_RA'], dtype=np.float64)
    dec = np.asarray(rows['TARGET_DEC'], dtype=np.float64)
    mask = (np.asarray(rows['ZCAT_PRIMARY'], dtype=bool)
            & (np.char.strip(rows['OBJTYPE'].astype('U')) == 'TGT')
            & (np.asarray(rows['ZWARN']) == 0)
            & (np.char.strip(rows['SPECTYPE'].astype('U')) == kind.decode())
            & np.isfinite(z) & (z >= MIN_RELIABLE_Z)
            & np.isfinite(ra) & (ra >= 0) & (ra < 360)
            & np.isfinite(dec) & (dec >= -90) & (dec <= 90)
            & (np.asarray(rows['TARGETID']) > 0))
    return rows[mask], z[mask], ra[mask], dec[mask]


def positions(z, ra, dec, z_grid, distance_grid):
    distance = np.interp(z, z_grid, distance_grid)
    ra_r = np.deg2rad(ra)
    dec_r = np.deg2rad(dec)
    flat = distance * np.cos(dec_r)
    return np.column_stack((flat * np.cos(ra_r),
                            flat * np.sin(ra_r), distance * np.sin(dec_r)))


def inspect(data):
    result = {name: {'count': 0, 'z_min': np.inf, 'z_max': 0.0}
              for name in CLASSES}
    for start in range(0, len(data), ROW_BATCH):
        end = min(start + ROW_BATCH, len(data))
        for name in CLASSES:
            rows, z, _, _ = selected(data, start, end, name)
            if len(rows):
                item = result[name]
                item['count'] += len(rows)
                item['z_min'] = min(item['z_min'], float(z.min()))
                item['z_max'] = max(item['z_max'], float(z.max()))
        print(f'\rCounting {end:,}/{len(data):,}', end='', flush=True)
    print()
    return result


def create_output(path, count, units, z_min, z_max, bounds):
    size = 48 + count * 32
    with path.open('wb') as stream:
        stream.truncate(size)
        stream.seek(0)
        stream.write(b'DESI')
        stream.write(struct.pack('<IIfff', 3, count, units, z_min, z_max))
        stream.write(struct.pack('<6f', *bounds))
    base = 48 + count * 8
    return {
        'pos': np.memmap(path, '<i2', 'r+', 48, (count, 3)),
        'redshift': np.memmap(path, '<u2', 'r+', 48 + count * 6, (count,)),
        'morph': np.memmap(path, 'u1', 'r+', base, (count,)),
        'shape': np.memmap(path, 'u1', 'r+', base + count, (count,)),
        'gr': np.memmap(path, 'u1', 'r+', base + count * 2, (count,)),
        'rz': np.memmap(path, 'u1', 'r+', base + count * 3, (count,)),
        'ra': np.memmap(path, '<f4', 'r+', base + count * 4, (count,)),
        'dec': np.memmap(path, '<f4', 'r+', base + count * 8, (count,)),
        'zwarn': np.memmap(path, '<u4', 'r+', base + count * 12, (count,)),
        'target': np.memmap(path, '<u8', 'r+', base + count * 16, (count,)),
    }


def main():
    global MIN_RELIABLE_Z
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--out', type=Path, default=DEST)
    parser.add_argument('--min-reliable-z', type=float, default=MIN_RELIABLE_Z,
                        help='minimum redshift included; use 0 for the historical full DR1 point selection')
    args = parser.parse_args()
    MIN_RELIABLE_Z = args.min_reliable_z
    if not args.source.is_file():
        parser.error(f'missing official DR1 FITS: {args.source}')
    args.out.mkdir(parents=True, exist_ok=True)

    with fits.open(args.source, memmap=True) as hdul:
        data = hdul[1].data
        required = {'TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'ZWARN',
                    'SPECTYPE', 'OBJTYPE', 'ZCAT_PRIMARY'}
        missing = required.difference(data.names)
        if missing:
            raise ValueError(f'missing FITS columns: {sorted(missing)}')
        stats = inspect(data)
        if not any(item['count'] for item in stats.values()):
            raise ValueError('no usable GALAXY or QSO rows in source FITS')
        largest_z = max(item['z_max'] for item in stats.values())
        z_grid = np.linspace(0, largest_z * 1.001, 16384)
        distance_grid = Planck18.comoving_distance(z_grid).value

        # Determine each output's quantization scale and exact bounding box.
        for name in CLASSES:
            stats[name]['bounds_min'] = np.full(3, np.inf)
            stats[name]['bounds_max'] = np.full(3, -np.inf)
        for start in range(0, len(data), ROW_BATCH):
            end = min(start + ROW_BATCH, len(data))
            for name in CLASSES:
                _, z, ra, dec = selected(data, start, end, name)
                if not len(z):
                    continue
                xyz = positions(z, ra, dec, z_grid, distance_grid)
                item = stats[name]
                item['bounds_min'] = np.minimum(item['bounds_min'], xyz.min(axis=0))
                item['bounds_max'] = np.maximum(item['bounds_max'], xyz.max(axis=0))
            print(f'\rMeasuring {end:,}/{len(data):,}', end='', flush=True)
        print()

        output = {}
        for name in CLASSES:
            item = stats[name]
            item['written'] = 0
            if not item['count']:
                continue
            extent = max(np.abs(item['bounds_min']).max(),
                         np.abs(item['bounds_max']).max())
            item['units'] = 32767 / extent
            item['written'] = 0
            path = args.out / f'DR1_{name.decode()}.desi'
            bounds = np.r_[item['bounds_min'], item['bounds_max']]
            output[name] = create_output(path, item['count'], item['units'],
                                         item['z_min'], item['z_max'], bounds)

        for start in range(0, len(data), ROW_BATCH):
            end = min(start + ROW_BATCH, len(data))
            for name, arrays in output.items():
                rows, z, ra, dec = selected(data, start, end, name)
                n = len(z)
                if not n:
                    continue
                item = stats[name]
                sl = slice(item['written'], item['written'] + n)
                xyz = positions(z, ra, dec, z_grid, distance_grid)
                arrays['pos'][sl] = np.rint(xyz * item['units']).astype('<i2')
                if item['z_max'] == item['z_min']:
                    arrays['redshift'][sl] = 0
                else:
                    arrays['redshift'][sl] = np.rint(
                        (z - item['z_min']) / (item['z_max'] - item['z_min']) * 65535
                    ).astype('<u2')
                arrays['morph'][sl] = 255
                arrays['shape'][sl] = 0
                arrays['gr'][sl] = 128
                arrays['rz'][sl] = 128
                arrays['ra'][sl] = ra.astype('<f4')
                arrays['dec'][sl] = dec.astype('<f4')
                arrays['zwarn'][sl] = 0
                arrays['target'][sl] = np.asarray(rows['TARGETID'], dtype='<u8')
                item['written'] += n
            print(f'\rWriting {end:,}/{len(data):,}', end='', flush=True)
        print()

    for arrays in output.values():
        for array in arrays.values():
            array.flush()

    galaxy_ids = None
    for name, arrays in output.items():
        ids = np.sort(np.asarray(arrays['target']))
        if np.any(ids[1:] == ids[:-1]):
            raise ValueError(f'duplicate TARGETID in {name.decode()} output')
        if galaxy_ids is not None:
            positions_in_galaxies = np.searchsorted(galaxy_ids, ids)
            overlaps = (positions_in_galaxies < len(galaxy_ids))
            if np.any(galaxy_ids[positions_in_galaxies[overlaps]] == ids[overlaps]):
                raise ValueError('TARGETID appears in both GALAXY and QSO outputs')
        else:
            galaxy_ids = ids

    summary = {}
    for name, item in stats.items():
        if item['written'] != item['count']:
            raise RuntimeError(f'{name.decode()}: count changed between passes')
        summary[name.decode()] = {key: item[key] for key in ('count', 'z_min', 'z_max')}
    (args.out / 'dr1_full_summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
