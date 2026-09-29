#!/usr/bin/env python3
"""Compare two BGS BRIGHT random realizations at low-z candidate positions."""
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
from extend_rsd_catalog import load_viewer_catalog  # noqa: E402
from rsd_pipeline import _fits_rows  # noqa: E402

SOURCE = ROOT / 'tools' / 'rsd' / 'lowz_sources'
FIELD = SOURCE / 'field_candidate'


def flags(region, index, points):
    path = SOURCE / f'BGS_BRIGHT_{region}_{index}_clustering.ran.fits'
    rows = _fits_rows(path, Planck18, 0, 0.01, 0.4)
    positions = np.asarray(rows['position'], dtype='f4')
    train = positions[::2]
    holdout = positions[1::2][:250000]
    tree = cKDTree(train)
    radius = float(np.percentile(tree.query(holdout,k=1,workers=-1)[0],99))
    distance = tree.query(points,k=1,workers=-1)[0]
    return distance > radius, {'file':path.name,'selectedRandomRows':len(positions),
                               'trainingRows':len(train),'p99RadiusMpc':radius}


def main():
    _, positions, redshift, viewer_ids = load_viewer_catalog(
        ROOT / 'assets', ROOT / 'ply' / 'chunks', ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    order = np.argsort(viewer_ids)
    sorted_ids = viewer_ids[order]
    report = {'status':'two-random-proximity-diagnostic', 'seed':20260925,'regions':[],
              'limitation':'Two realizations share the same DESI selection recipe; neither validates the official mask independently.'}
    rng = np.random.default_rng(20260925)
    for region in ('NGC','SGC'):
        with np.load(FIELD / f'BGS_BRIGHT_{region}_viewer_rsd_v1.npz') as entry:
            z = np.asarray(entry['z'])
            ids = np.asarray(entry['targetid'])
        sample_indices=[]
        for low,high in ((0.01,0.1),(0.1,0.12)):
            available=np.flatnonzero((z>=low)&(z<high))
            sample_indices.extend(rng.choice(available,size=min(10000,len(available)),replace=False))
        sample_ids=ids[np.asarray(sample_indices)]
        at=np.searchsorted(sorted_ids,sample_ids)
        if np.any(at>=len(sorted_ids)) or not np.array_equal(sorted_ids[at],sample_ids):
            raise ValueError('candidate/viewer TARGETID mismatch')
        query=positions[order[at]]
        bad0,meta0=flags(region,0,query)
        bad1,meta1=flags(region,1,query)
        item={'region':region,'sampleRows':len(query),'index0':meta0,'index1':meta1,
              'flaggedIndex0':int(bad0.sum()),'flaggedIndex1':int(bad1.sum()),
              'flaggedBoth':int((bad0&bad1).sum()),'flaggedEither':int((bad0|bad1).sum())}
        report['regions'].append(item)
        (SOURCE / 'random-realization-audit.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(item),flush=True)


if __name__ == '__main__':
    main()
