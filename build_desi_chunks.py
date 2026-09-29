#!/usr/bin/env python3
"""Package DESI point positions into chunks for progressive loading.

The objective is faster first display. Sorting by Morton order and splitting
into contiguous blocks compresses well, but the first block covers only a
corner of the volume. Spatial hash partitioning distributes each chunk across
the survey; Morton ordering within each partition retains compression.

manifest.json describes the chunks. Each chunk_NNN.bin is gzip-compressed
int16[count] dx | int16[count] dy | int16[count] dz, with per-axis deltas
relative to the previous point. The browser decompresses it directly.
Photometric color is excluded; the display shader derives color from radius.

Usage:
  python build_desi_chunks.py
  python build_desi_chunks.py --chunks 16
"""

import argparse
import gzip
import json
import struct
import sys
from pathlib import Path

import numpy as np

ASSETS = Path(__file__).parent / 'assets'
OUT_DIR = Path(__file__).parent / 'ply' / 'chunks'
LEGACY_TRACERS = ['BGS_BRIGHT', 'LRG', 'ELG', 'QSO']
DR1_TRACERS = ['DR1_GALAXY', 'DR1_QSO']


def read_desi(path: Path):
    raw = path.read_bytes()
    if raw[:4] != b'DESI':
        raise ValueError(f'{path.name}: missing DESI signature')
    version, count, units_per_mpc, z_min, z_max = struct.unpack('<IIfff', raw[4:24])
    pos = np.frombuffer(raw, np.int16, count * 3, 48).reshape(count, 3)
    base = 48 + count * 8
    gr = np.frombuffer(raw, np.uint8, count, base + count * 2)
    rz = np.frombuffer(raw, np.uint8, count, base + count * 3)
    if version < 3:
        raise ValueError(f'{path.name}: .desi v3 with TARGETID is required for the FoG index')
    target_id = np.frombuffer(raw, np.uint64, count, base + count * 16)
    return pos, float(units_per_mpc), float(z_max), gr, rz, target_id


def morton_order(pos: np.ndarray) -> np.ndarray:
    """Return Morton order on a 10-bit grid per axis.

    Ten bits preserve enough locality for delta coding. Positions remain
    stored at full int16 precision.
    """
    q = ((pos.astype(np.int32) + 32768).astype(np.uint32) >> 6)

    def spread(x):
        x = x.astype(np.uint64) & 0x3ff
        x = (x | (x << 16)) & 0x30000ff
        x = (x | (x << 8)) & 0x300f00f
        x = (x | (x << 4)) & 0x30c30c3
        x = (x | (x << 2)) & 0x9249249
        return x

    key = spread(q[:, 0]) | (spread(q[:, 1]) << 1) | (spread(q[:, 2]) << 2)
    return np.argsort(key, kind='stable')


def encode_chunk(pos: np.ndarray) -> bytes:
    """Gzip-compress axis-separated int16 deltas.

    Prepend zero rather than pos[0] so the first delta is an absolute position;
    otherwise the decoder shifts the whole chunk by its first point.
    """
    d = np.diff(pos.astype(np.int32), axis=0, prepend=0).astype(np.int16)
    return gzip.compress(np.ascontiguousarray(d.T).tobytes(), 9)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--chunks', type=int, default=12, help='number of chunks (default: 12)')
    ap.add_argument('--out', default=str(OUT_DIR))
    ap.add_argument('--tag', default='',
                    help='version suffix for cache invalidation, e.g. 20260819')
    ap.add_argument('--fog-index-dir', default='',
                    help='offline-only TARGETID index directory; never publish it')
    args = ap.parse_args()
    tag = f'_{args.tag}' if args.tag else ''

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    dr1_present = [(ASSETS / f'{name}.desi').exists() for name in DR1_TRACERS]
    if any(dr1_present) and not all(dr1_present):
        raise FileNotFoundError('DR1 requires both DR1_GALAXY.desi and DR1_QSO.desi')
    tracers = DR1_TRACERS if all(dr1_present) else LEGACY_TRACERS
    tracer_codes = {name: code for code, name in enumerate(tracers)}

    # Each .desi file has its own unitsPerMpc. Convert to Mpc before combining
    # tracers and quantizing them with one shared scale.
    parts, types, target_parts, flag_parts = [], [], [], []
    max_redshift = 0.0
    for tracer in tracers:
        src = ASSETS / f'{tracer}.desi'
        if not src.exists():
            print(f'  [skip]  missing {src.name}')
            continue
        pos_q, u, z_max, gr, rz, target_id = read_desi(src)
        max_redshift = max(max_redshift, z_max)
        mpc = pos_q.astype(np.float32) / u
        parts.append(mpc)
        types.append(np.full(len(mpc), tracer_codes[tracer], dtype=np.uint8))
        target_parts.append(target_id)
        if tracers == DR1_TRACERS:
            flags_path = ASSETS / f'{tracer}.target-flags.bin'
            if not flags_path.exists() or flags_path.stat().st_size != len(mpc):
                raise ValueError(f'{flags_path.name}: missing or wrong row count')
            flag_parts.append(np.memmap(flags_path, dtype='u1', mode='r'))
        print(f'  [read]  {tracer}: {len(mpc):,}  (source unitsPerMpc {u:.3f}, '
              f'max radius {np.linalg.norm(mpc, axis=1).max():.0f} Mpc)')

    if not parts:
        print('no data')
        return 1

    mpc = np.concatenate(parts)
    type_all = np.concatenate(types)
    target_flags = np.concatenate(flag_parts) if flag_parts else None
    target_all = np.concatenate(target_parts)
    n = len(mpc)

    # The maximum absolute coordinate sets the shared scale, with a 2% margin.
    extent = float(np.abs(mpc).max()) * 1.02
    units = 32767.0 / extent
    pos = np.rint(mpc * units).astype(np.int16)
    err = float(np.abs(pos.astype(np.float32) / units - mpc).max())
    print(f'\nTotal: {n:,} objects')
    print(f'Global scale: {units:.4f} units/Mpc  (extent {extent:.0f} Mpc, '
          f'max error {err * 1000:.1f} kpc)')

    order = morton_order(pos)
    pos_m = pos[order]

    # Reference size for contiguous blocks, which compress better but load unevenly.
    contiguous = sum(len(encode_chunk(c)) for c in np.array_split(pos_m, args.chunks))

    manifest = {
        'version': 3,
        'count': int(n),
        'zMax': max_redshift,
        'unitsPerMpc': units,
        'partition': 'position-hash-v1',
        'tracers': tracers,
        'tracerCounts': {name: int((type_all == code).sum())
                         for name, code in tracer_codes.items()},
        'chunks': []
    }
    if target_flags is not None:
        manifest['typeEncoding'] = 'dr1-target-bits-v1'
        manifest['targetCounts'] = {name: int(np.count_nonzero(target_flags & (1 << bit)))
                                    for bit, name in enumerate(
                                        ('BGS', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO_TARGET', 'OTHER'))}
        # Bits 0-5 retain overlapping targeting labels; bit 6 is spectral QSO.
        type_all = target_flags | (type_all << 6)
    fog_index_dir = Path(args.fog_index_dir) if args.fog_index_dir else None
    if fog_index_dir:
        fog_index_dir.mkdir(parents=True, exist_ok=True)
    total = 0
    ttotal = 0
    position_hash = np.abs(pos[:, 0].astype(np.int64) * 17 +
                           pos[:, 2].astype(np.int64) * 59 +
                           pos[:, 1].astype(np.int64) * 101)
    for i in range(args.chunks):
        # Deterministic spatial partition. Clients can infer a point's chunk
        # from its position alone. Morton ordering within each bucket aids gzip.
        chosen = np.flatnonzero(position_hash % args.chunks == i)
        local_order = morton_order(pos[chosen])
        chosen = chosen[local_order]
        chunk = pos[chosen]
        blob = encode_chunk(chunk)
        name = f'chunk_{i:03d}{tag}.bin'
        (out_dir / name).write_bytes(blob)

        # One byte per object preserves targeting/class information. It can be
        # fetched after positions to avoid delaying the first frame.
        tname = f'type_{i:03d}{tag}.bin'
        tblob = gzip.compress(type_all[chosen].tobytes(), 9)
        (out_dir / tname).write_bytes(tblob)

        manifest['chunks'].append({'file': name, 'count': int(len(chunk)),
                                   'bytes': len(blob),
                                   'type': tname, 'typeBytes': len(tblob)})
        if fog_index_dir:
            # Keep this index offline to align FoG with exact published geometry
            # without sending TARGETIDs or a large hash table to the browser.
            iname = f'targetid_{i:03d}{tag}.npy'
            np.save(fog_index_dir / iname, target_all[chosen])
            manifest['chunks'][-1]['fogTargetIndex'] = iname
        total += len(blob)
        ttotal += len(tblob)
        print(f'  [write] {name}: {len(chunk):,} points, {len(blob) / 1e6:.2f} MB'
              f'  + type {len(tblob) / 1e6:.2f} MB')

    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=1))

    first = manifest['chunks'][0]
    print(f'\nTotal      : {total / 1e6:.1f} MB in {args.chunks} chunks')
    print(f'Types      : {ttotal / 1e6:.1f} MB separately (asynchronous)')
    print(f'Contiguous : {contiguous / 1e6:.1f} MB '
          f'(+{(total - contiguous) / 1e6:.1f} MB partition overhead)')
    print(f'First      : {first["bytes"] / 1e6:.2f} MB, {first["count"]:,} points')
    print(f'Output     : {out_dir}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
