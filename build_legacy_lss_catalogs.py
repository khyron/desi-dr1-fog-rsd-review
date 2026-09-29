"""Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""

import argparse
import struct
import sys
from pathlib import Path
from urllib.request import urlopen, Request

import numpy as np

BASE_URL = 'https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5'
TRACERS = {
    'BGS_BRIGHT': {
        'parts': ['BGS_BRIGHT_NGC', 'BGS_BRIGHT_SGC'],
        'full': 'BGS_BRIGHT_full_HPmapcut.dat.fits',
    },
    'LRG': {
        'parts': ['LRG_NGC', 'LRG_SGC'],
        'full': 'LRG_full_HPmapcut.dat.fits',
    },
    'ELG': {
        'parts': ['ELG_LOPnotqso_NGC', 'ELG_LOPnotqso_SGC'],
        'full': 'ELG_LOPnotqso_full_HPmapcut.dat.fits',
    },
    'QSO': {
        'parts': ['QSO_NGC', 'QSO_SGC'],
        'full': 'QSO_full_HPmapcut.dat.fits',
    },
}

DATA_DIR = Path(__file__).parent / 'data'
OUT_DIR = Path(__file__).parent / 'assets'


def remote_size(url: str) -> int:
    req = Request(url, method='HEAD', headers={'User-Agent': 'desi-splats/1.0'})
    with urlopen(req) as r:
        return int(r.headers.get('Content-Length', 0))


def download(fname: str, attempts: int = 12) -> Path:
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    dest = DATA_DIR / fname
    url = f'{BASE_URL}/{fname}'
    total = remote_size(url)

    if dest.exists():
        have = dest.stat().st_size
        if have == total:
            print(f'  [cache] {fname} ({have / 1e6:.0f} MB)')
            return dest
        print(f'  [redo]  {fname} incomplete ({have/1e6:.0f}/{total/1e6:.0f} MB)')
        dest.unlink()

    print(f'  [get]   {url}  ({total / 1e6:.0f} MB)')
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix('.part')

    for attempt in range(1, attempts + 1):
        done = tmp.stat().st_size if tmp.exists() else 0
        if done >= total:
            break
        headers = {'User-Agent': 'desi-splats/1.0'}
        if done:
            headers['Range'] = f'bytes={done}-'
        try:
            with urlopen(Request(url, headers=headers)) as r:
                if done and r.status != 206:
                    done = 0
                    tmp.unlink(missing_ok=True)
                with open(tmp, 'ab' if done else 'wb') as f:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        f.write(chunk)
                        done += len(chunk)
                        print(f'\r          {100*done/total:5.1f}%  '
                              f'{done/1e6:6.0f}/{total/1e6:.0f} MB', end='', flush=True)
        except OSError as e:
            print(f'\n  [retry] attempt {attempt}: {e}')
            continue
        print()

    done = tmp.stat().st_size if tmp.exists() else 0
    if done != total:
        raise RuntimeError(f'{fname}: incomplete download ({done}/{total} bytes) '
                           f'after {attempts} attempts; rerun to resume')
    tmp.replace(dest)
    return dest


def comoving_distance_table(z_max: float):
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    from astropy.cosmology import Planck18

    grid = np.linspace(0.0, z_max * 1.01 + 1e-3, 4096)
    dist = Planck18.comoving_distance(grid).value  # Mpc
    return grid, dist


def read_catalog(path: Path):
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        tid = np.asarray(d['TARGETID'], dtype=np.int64)
        ra = np.asarray(d['RA'], dtype=np.float64)
        dec = np.asarray(d['DEC'], dtype=np.float64)
        z = np.asarray(d['Z'], dtype=np.float64)
    return tid, ra, dec, z


MORPH_CODES = {'PSF': 0, 'REX': 1, 'EXP': 2, 'DEV': 3, 'SER': 4}


def read_morphology(path: Path):
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        d = hdul[1].data
        tid = np.asarray(d['TARGETID'], dtype=np.int64)
        morph_raw = np.asarray(d['MORPHTYPE'])
        shape_r = np.asarray(d['SHAPE_R'], dtype=np.float32)
        flux_g = np.asarray(d['FLUX_G'], dtype=np.float32)
        flux_r = np.asarray(d['FLUX_R'], dtype=np.float32)
        flux_z = np.asarray(d['FLUX_Z'], dtype=np.float32)
        zwarn = np.asarray(d['ZWARN'], dtype=np.uint32)

    morph = np.full(len(tid), 255, dtype=np.uint8)
    for name, code in MORPH_CODES.items():
        morph[np.char.strip(morph_raw.astype('U8')) == name] = code

    order = np.argsort(tid, kind='stable')
    return (tid[order], morph[order], shape_r[order], flux_g[order], flux_r[order],
            flux_z[order], zwarn[order])


def nanomaggie_to_mag(flux):
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(flux > 0, 22.5 - 2.5 * np.log10(np.maximum(flux, 1e-12)), np.nan)


def attach_morphology(tracer: str, spec: dict, tid: np.ndarray, skip_download: bool):
    """Build historical v3 LSS catalogues with exact TARGETIDs and imaging metadata."""
    fname = spec['full']
    path = DATA_DIR / fname

    if skip_download:
        if not path.exists():
            print(f'  [skip]  missing {fname}: writing without morphology')
            return None, None, None, None, None
    else:
        path = download(fname)

    print(f'  [morph] leyendo {fname}')
    f_tid, f_morph, f_shape, f_g, f_r, f_z, f_zwarn = read_morphology(path)
    slot = np.searchsorted(f_tid, tid)
    slot = np.clip(slot, 0, len(f_tid) - 1)
    hit = f_tid[slot] == tid

    n = len(tid)
    morph = np.full(n, 255, dtype=np.uint8)
    shape_r = np.zeros(n, dtype=np.float32)
    mag_g = np.full(n, np.nan, dtype=np.float32)
    mag_r = np.full(n, np.nan, dtype=np.float32)
    mag_z = np.full(n, np.nan, dtype=np.float32)
    zwarn = np.full(n, np.uint32(0xffffffff), dtype=np.uint32)

    morph[hit] = f_morph[slot[hit]]
    shape_r[hit] = f_shape[slot[hit]]
    mag_g[hit] = nanomaggie_to_mag(f_g[slot[hit]])
    mag_r[hit] = nanomaggie_to_mag(f_r[slot[hit]])
    mag_z[hit] = nanomaggie_to_mag(f_z[slot[hit]])
    zwarn[hit] = f_zwarn[slot[hit]]

    print(f'  [morph] {hit.sum():,} de {n:,} emparejados ({100*hit.mean():.1f}%)')
    counts = {name: int((morph == code).sum()) for name, code in MORPH_CODES.items()}
    print('  [morph] ' + '  '.join(f'{k}={v:,}' for k, v in counts.items()))
    shape_r_q = np.clip(np.rint(shape_r * 10), 0, 255).astype(np.uint8)
    gr = mag_g - mag_r
    rz = mag_r - mag_z
    gr_q = np.where(np.isfinite(gr), np.clip(np.rint((gr + 0.5) * 85), 0, 255), 128).astype(np.uint8)
    rz_q = np.where(np.isfinite(rz), np.clip(np.rint((rz + 0.5) * 102), 0, 255), 128).astype(np.uint8)

    return morph, shape_r_q, gr_q, rz_q, zwarn


def build_tracer(tracer: str, spec: dict, skip_download: bool, with_morphology: bool) -> None:
    print(f'\n{tracer}')

    tids, ras, decs, zs = [], [], [], []
    for part in spec['parts']:
        fname = f'{part}_clustering.dat.fits'
        path = DATA_DIR / fname
        if skip_download:
            if not path.exists():
                print(f'  [skip]  missing {fname}')
                continue
        else:
            path = download(fname)
        tid, ra, dec, z = read_catalog(path)
        print(f'  [read]  {path.name}: {len(ra):,} objects')
        tids.append(tid)
        ras.append(ra)
        decs.append(dec)
        zs.append(z)

    if not ras:
        print(f'  [warn]  no data for {tracer}; skipping')
        return

    tid = np.concatenate(tids)
    ra = np.concatenate(ras)
    dec = np.concatenate(decs)
    z = np.concatenate(zs)
    ok = np.isfinite(z) & (z > 0.0) & (z < 6.0) & np.isfinite(ra) & np.isfinite(dec)
    tid, ra, dec, z = tid[ok], ra[ok], dec[ok], z[ok]
    print(f'  [clean] {len(z):,} valid objects (z {z.min():.3f} - {z.max():.3f})')

    morph = shape_r_q = gr_q = rz_q = zwarn = None
    if with_morphology:
        morph, shape_r_q, gr_q, rz_q, zwarn = attach_morphology(
            tracer, spec, tid, skip_download)

    grid, dist_grid = comoving_distance_table(float(z.max()))
    d = np.interp(z, grid, dist_grid)  # comoving Mpc

    ra_r = np.radians(ra)
    dec_r = np.radians(dec)
    cos_dec = np.cos(dec_r)
    x = d * cos_dec * np.cos(ra_r)
    y = d * cos_dec * np.sin(ra_r)
    zc = d * np.sin(dec_r)

    extent = float(max(np.abs(x).max(), np.abs(y).max(), np.abs(zc).max()))
    units_per_mpc = 32767.0 / extent  # aprovecha todo el rango int16

    pos = np.empty((len(z), 3), dtype=np.int16)
    pos[:, 0] = np.rint(x * units_per_mpc)
    pos[:, 1] = np.rint(y * units_per_mpc)
    pos[:, 2] = np.rint(zc * units_per_mpc)

    z_max = float(z.max())
    z_min = float(z.min())
    z_q = np.rint((z - z_min) / (z_max - z_min) * 65535.0).astype(np.uint16)

    version = 3 if morph is not None else 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f'{tracer}.desi'
    with open(out, 'wb') as f:
        f.write(b'DESI')
        f.write(struct.pack('<IIfff', version, len(z), units_per_mpc, z_min, z_max))
        f.write(struct.pack('<6f',
                            float(x.min()), float(y.min()), float(zc.min()),
                            float(x.max()), float(y.max()), float(zc.max())))
        f.write(pos.tobytes())
        f.write(z_q.tobytes())
        if morph is not None:
            f.write(morph.tobytes())
            f.write(shape_r_q.tobytes())
            f.write(gr_q.tobytes())
            f.write(rz_q.tobytes())
            f.write(ra.astype('<f4').tobytes())
            f.write(dec.astype('<f4').tobytes())
            f.write(zwarn.astype('<u4', copy=False).tobytes())
            f.write(tid.astype('<u8', copy=False).tobytes())

    mb = out.stat().st_size / 1e6
    print(f'  [write] {out.name}: v{version}, {len(z):,} objects, {mb:.1f} MB, '
          f'{extent:.0f} Mpc extent, {1/units_per_mpc*1000:.1f} kpc/unit')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tracer', choices=sorted(TRACERS), help='process only this tracer')
    ap.add_argument('--skip-download', action='store_true',
                    help='use only FITS files already present in data/')
    ap.add_argument('--no-morphology', action='store_true',
                    help='write v1 without downloading the _full catalogues (about 20 GB)')
    args = ap.parse_args()

    todo = {args.tracer: TRACERS[args.tracer]} if args.tracer else TRACERS
    for tracer, spec in todo.items():
        build_tracer(tracer, spec, args.skip_download, not args.no_morphology)

    print(f'\nComplete. Binary files in {OUT_DIR}')


if __name__ == '__main__':
    sys.exit(main())
