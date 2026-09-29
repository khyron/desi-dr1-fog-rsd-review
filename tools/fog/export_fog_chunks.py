#!/usr/bin/env python3
"""Create gzip Float32 delta chunks aligned with regenerated point chunks."""
from __future__ import annotations
import argparse,gzip,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser();p.add_argument('--chunks',required=True);p.add_argument('--target-index-dir',required=True);p.add_argument('--deltas',default=str(ROOT/'data_external/fog/analysis/targetid-delta.npy'));p.add_argument('--out',required=True);a=p.parse_args()
 chunks=Path(a.chunks);out=Path(a.out);out.mkdir(parents=True,exist_ok=True);man=json.loads((chunks/'manifest.json').read_text());d=np.load(a.deltas,mmap_mode='r'); ids=d['targetid']; vals=d['delta_mpc']; total=0;nonzero=0
 for i,c in enumerate(man['chunks']):
  targets=np.load(Path(a.target_index_dir)/c['fogTargetIndex'],mmap_mode='r');slot=np.searchsorted(ids,targets);hit=(slot<len(ids))&(ids[np.minimum(slot,len(ids)-1)]==targets);delta=np.zeros(len(targets),np.float32);delta[hit]=vals[slot[hit]]
  name=f'fog_{i:03d}.bin';blob=gzip.compress(delta.tobytes(),compresslevel=9);(out/name).write_bytes(blob);c['fog']=name;c['fogBytes']=len(blob);c.pop('fogTargetIndex',None);total+=len(blob);nonzero+=int(np.count_nonzero(delta));
  if len(delta)!=c['count']:raise RuntimeError('count mismatch '+name)
 man['fog']={'format':'gzip-float32le-delta-mpc-v1','algorithm':'gfinder-geometry-fog-v1','bytes':total,'nonzero':nonzero};(out/'manifest.json').write_text(json.dumps(man,indent=1));print(json.dumps({'chunks':len(man['chunks']),'bytes':total,'nonzero':nonzero},indent=2))
if __name__=='__main__':main()
