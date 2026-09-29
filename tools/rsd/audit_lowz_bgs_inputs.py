#!/usr/bin/env python3
"""Read-only input audit for a proposed DR1 BGS BRIGHT low-z reconstruction."""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
from extend_rsd_catalog import load_viewer_catalog  # noqa: E402

SOURCE = ROOT / 'tools' / 'rsd' / 'lowz_sources'
OUT = SOURCE / 'input-audit.json'


def stats(path):
    with fits.open(path, memmap=True) as hdus:
        table = hdus[1].data
        z = np.asarray(table['Z'])
        weight = np.asarray(table['WEIGHT'])
        if 'WEIGHT_FKP' in table.names:
            weight = weight * np.asarray(table['WEIGHT_FKP'])
        valid = np.isfinite(z) & np.isfinite(weight) & (weight > 0)
        result = {'file': path.name, 'rows': len(table), 'desidr': hdus[1].header.get('DESIDR'),
                  'columns': list(table.names), 'finitePositiveWeightRows': int(valid.sum()),
                  'zMin': float(np.nanmin(z)), 'zMax': float(np.nanmax(z)),
                  'z0p01to0p10': int(np.count_nonzero(valid & (z >= 0.01) & (z < 0.1))),
                  'z0p10to0p40': int(np.count_nonzero(valid & (z >= 0.1) & (z < 0.4)))}
        if path.name.endswith('.dat.fits'):
            ids = np.asarray(table['TARGETID'][valid & (z >= 0.01) & (z < 0.1)], dtype='u8')
            result['uniqueLowzTargetIds'] = int(len(np.unique(ids)))
        return result


def main():
    names = [f'BGS_BRIGHT_{region}_{kind}_clustering.ran.fits' if kind == '0'
             else f'BGS_BRIGHT_{region}_clustering.dat.fits'
             for region in ('NGC', 'SGC') for kind in ('data', '0')]
    files = [stats(SOURCE / name) for name in names]
    manifest, _, viewer_z, viewer_ids = load_viewer_catalog(
        ROOT / 'assets', ROOT / 'ply' / 'chunks', ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    viewer_low = (viewer_z >= 0.01) & (viewer_z < 0.1)
    low_ids = np.sort(viewer_ids[viewer_low])
    matched = np.zeros(len(low_ids), dtype=bool)
    by_region = {}
    for region in ('NGC', 'SGC'):
        path = SOURCE / f'BGS_BRIGHT_{region}_clustering.dat.fits'
        with fits.open(path, memmap=True) as hdus:
            d = hdus[1].data
            z = np.asarray(d['Z'])
            ids = np.sort(np.asarray(d['TARGETID'][(z >= 0.01) & (z < 0.1)], dtype='u8'))
            at = np.searchsorted(ids, low_ids)
            hit = (at < len(ids)) & (ids[np.minimum(at, len(ids)-1)] == low_ids)
            by_region[region] = int(hit.sum())
            matched |= hit
    report = {'catalog': 'DESI DR1 LSS iron LSScats v1.5 BGS_BRIGHT',
              'sourceFiles': files, 'viewerCount': manifest['count'],
              'viewerZ0p01to0p10': int(viewer_low.sum()),
              'viewerExactDataMatchesByRegion': by_region,
              'viewerExactDataMatchesUnion': int(matched.sum()),
              'viewerLowzWithoutExactDataMatch': int((~matched).sum()),
              'interpretation': 'Exact data identity only; unmatched viewer rows may still be queryable within survey support after field reconstruction. No correction generated.'}
    OUT.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('viewerZ0p01to0p10','viewerExactDataMatchesByRegion','viewerExactDataMatchesUnion','viewerLowzWithoutExactDataMatch')}, indent=2))


if __name__ == '__main__':
    main()
