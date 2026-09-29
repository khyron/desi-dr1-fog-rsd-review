#!/usr/bin/env python3
"""Conservative geometry FoG analysis. Outputs sorted TARGETID -> delta Mpc."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from astropy.cosmology import Planck18

ROOT=Path(__file__).resolve().parents[2]; FOG=ROOT/'data_external'/'fog'
def join(parts): return np.concatenate([np.load(x,mmap_mode='r') for x in sorted(parts)])
def robust(x):
    m=np.median(x); return 1.4826*np.median(np.abs(x-m))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--members',default=str(FOG/'members'));p.add_argument('--groups',default=str(FOG/'groups'));p.add_argument('--matches',default=str(FOG/'matches'));p.add_argument('--out',default=str(FOG/'analysis'));p.add_argument('--min-spec',type=int,default=5);p.add_argument('--elongation',type=float,default=1.5);p.add_argument('--alpha-min',type=float,default=.08);p.add_argument('--max-delta-mpc',type=float,default=300.);a=p.parse_args()
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 mem=join(Path(a.members).glob('member_*.npy')); gal=join(Path(a.members).glob('galaxy_*.npy')); grp=join(Path(a.groups).glob('group_*.npy'))
 # Exact IGAL join after stable sort; no positional fallback.
 mo=np.argsort(mem['igal']);go=np.argsort(gal['igal']);mem=mem[mo];gal=gal[go]
 if len(mem)!=len(gal) or not np.array_equal(mem['igal'],gal['igal']):raise RuntimeError('IGAL membership/galaxy mismatch')
 gp=np.argsort(grp['igrp']);grp=grp[gp]; slot=np.searchsorted(grp['igrp'],mem['igrp']); group_rows=grp[slot]
 zmax=float(max(gal['z'].max(),group_rows['z'].max())); zg=np.linspace(0,zmax*1.001+1e-4,16384); dg=Planck18.comoving_distance(zg).value
 chi=np.interp(gal['z'],zg,dg); gchi=np.interp(group_rows['z'],zg,dg); delta=np.zeros(len(mem),np.float32); stats=[]
 order=np.argsort(mem['igrp'],kind='stable'); starts=np.r_[0,np.flatnonzero(np.diff(mem['igrp'][order]))+1,len(order)]
 for lo,hi in zip(starts[:-1],starts[1:]):
  ix=order[lo:hi];g=group_rows[ix[0]];spec=ix[gal['zsrc'][ix]>0]
  if len(spec)<a.min_spec:continue
  radial=chi[spec]-gchi[spec];sp=robust(radial)
  ra=np.deg2rad(gal['ra'][spec]);dec=np.deg2rad(gal['dec'][spec]);ra0=np.deg2rad(g['ra']);dec0=np.deg2rad(g['dec'])
  cos=np.sin(dec)*np.sin(dec0)+np.cos(dec)*np.cos(dec0)*np.cos(ra-ra0);theta=np.arccos(np.clip(cos,-1,1));st=robust(gchi[spec]*theta)
  if not(np.isfinite(sp) and np.isfinite(st) and sp>0 and st>0):continue
  elong=sp/st
  if elong<a.elongation:continue
  alpha=max(a.alpha_min,min(1.,st/sp));local=(alpha*(chi[ix]-gchi[ix])-(chi[ix]-gchi[ix]));local=np.clip(local,-a.max_delta_mpc,a.max_delta_mpc);delta[ix]=local.astype(np.float32)
  stats.append((int(g['igrp']),int(g['rich']),len(spec),float(sp),float(st),float(elong),float(alpha)))
 matches=join(Path(a.matches).glob('match_*.npy')); io=np.argsort(mem['igal']);pos=np.searchsorted(mem['igal'][io],matches['igal']);hit=(pos<len(io))&(mem['igal'][io[np.minimum(pos,len(io)-1)]]==matches['igal']); outrows=np.empty(int(hit.sum()),dtype=[('targetid','<u8'),('delta_mpc','<f4')]);outrows['targetid']=matches['targetid'][hit];outrows['delta_mpc']=delta[io[pos[hit]]];outrows=outrows[np.argsort(outrows['targetid'],kind='stable')]
 np.save(out/'targetid-delta.npy',outrows);np.save(out/'group-stats.npy',np.array(stats,dtype=[('igrp','<i8'),('rich','<i4'),('spec','<i4'),('sigma_parallel','<f4'),('sigma_perp','<f4'),('elongation','<f4'),('alpha','<f4')]))
 summary={'eligibleGroups':int(len(grp)),'correctedGroups':len(stats),'matchedTargets':int(len(outrows)),'nonzeroDeltas':int(np.count_nonzero(outrows['delta_mpc'])),'minSpec':a.min_spec,'elongationThreshold':a.elongation,'alphaMin':a.alpha_min,'maxDeltaMpc':a.max_delta_mpc,'cosmology':'Planck18'};(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
