#!/usr/bin/env python3
"""Collect complete Gfinder members only for FoG-eligible groups, streamed FITS."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[2]
FOG=ROOT/'data_external'/'fog'
GF=ROOT/'data_external'/'desi'/'dr1'/'gfinder'/'v1.0'
MEM=np.dtype([('igal','<i8'),('igrp','<i8'),('rank','<i2')])
GAL=np.dtype([('igal','<i8'),('ra','<f8'),('dec','<f8'),('z','<f4'),('zsrc','<i2')])

def contains(ids, values):
    slot=np.searchsorted(ids,values)
    return (slot<len(ids)) & (ids[np.minimum(slot,len(ids)-1)]==values)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--groups-dir',default=str(FOG/'groups'));p.add_argument('--out',default=str(FOG/'members'))
    p.add_argument('--relation',default=str(GF/'iDESIDR9.y1.v1_1.fits'));p.add_argument('--galaxy',default=str(GF/'DESIDR9.y1.v1_galaxy.fits'))
    p.add_argument('--minimum-richness',type=int,default=5);p.add_argument('--batch',type=int,default=1_000_000)
    a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    groups=np.concatenate([np.load(f,mmap_mode='r') for f in sorted(Path(a.groups_dir).glob('group_*.npy'))])
    eligible=np.unique(groups['igrp'][groups['rich']>=a.minimum_richness]);np.save(out/'eligible-igrps.npy',eligible)
    member_paths=[];count=0
    with fits.open(a.relation,memmap=True,lazy_load_hdus=True) as h:
        d=h[1].data
        for start in range(0,len(d),a.batch):
            stop=min(len(d),start+a.batch);b=d[start:stop];keep=contains(eligible,np.asarray(b['IGRP']))
            if not keep.any():continue
            rows=np.empty(int(keep.sum()),dtype=MEM);rows['igal']=np.asarray(b['IGAL'])[keep];rows['igrp']=np.asarray(b['IGRP'])[keep];rows['rank']=np.asarray(b['RANK'])[keep]
            f=out/f'member_{start:09d}_{stop:09d}.npy';np.save(f,rows);member_paths.append(f);count+=len(rows)
            print(f'relation {stop:,}/{len(d):,}: members {count:,}',flush=True)
    wanted=np.unique(np.concatenate([np.load(f,mmap_mode='r')['igal'] for f in member_paths]));np.save(out/'member-igals.npy',wanted)
    gal_paths=[];gal_count=0
    with fits.open(a.galaxy,memmap=True,lazy_load_hdus=True) as h:
        d=h[1].data
        for start in range(0,len(d),a.batch):
            stop=min(len(d),start+a.batch);b=d[start:stop];keep=contains(wanted,np.asarray(b['IGAL']))
            if not keep.any():continue
            rows=np.empty(int(keep.sum()),dtype=GAL)
            for key in rows.dtype.names:rows[key]=np.asarray(b[key.upper()])[keep]
            f=out/f'galaxy_{start:09d}_{stop:09d}.npy';np.save(f,rows);gal_paths.append(f);gal_count+=len(rows)
            print(f'galaxy {stop:,}/{len(d):,}: members {gal_count:,}',flush=True)
    summary={'minimumRichness':a.minimum_richness,'eligibleGroups':int(len(eligible)),'members':int(count),'memberGalaxiesResolved':int(gal_count)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
