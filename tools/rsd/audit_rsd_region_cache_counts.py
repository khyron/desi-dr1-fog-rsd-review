#!/usr/bin/env python3
"""Append exact NGC/SGC cache identity counts to the support audit JSON."""
import argparse, json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'rsd'))
from audit_rsd_catalog_coverage import FIELDS, current_redshift_chunks
from rsd_pipeline import _fits_rows, local_sources

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--chunks',type=Path,default=ROOT/'ply/chunks')
    ap.add_argument('--assets',type=Path,default=ROOT/'assets')
    ap.add_argument('--index',type=Path,default=ROOT/'ply/reconstruction_dr1/private_index')
    ap.add_argument('--data-dir',type=Path,default=Path('/mnt/d/desi/data_external/RSD'))
    ap.add_argument('--cache-dir',type=Path,default=Path('/mnt/d/desi/data_external/RSD/derived/rsd-v1-exploratory'))
    ap.add_argument('--report',type=Path,default=ROOT/'tools/rsd/rsd-field-support-audit.json')
    a=ap.parse_args()
    field_cache={}; output={f:{r:{'cacheRows':0,'exactCatalogMatches':0} for r in ('NGC','SGC')} for f,_,_ in FIELDS}
    for f,_,_ in FIELDS:
        for r in ('NGC','SGC'):
            path=a.cache_dir/f'{f}_{r}_rsd_query_v1.npz'
            with np.load(path) as npz: ids=np.asarray(npz['targetid'],dtype='u8')
            ids.sort(kind='stable')
            field_cache[(f,r)]=ids; output[f][r]['cacheRows']=int(len(ids))
    overlaps={}
    for f,_,_ in FIELDS:
        a_ids=field_cache[(f,'NGC')]; b_ids=field_cache[(f,'SGC')]
        overlaps[f]=int(np.intersect1d(a_ids,b_ids,assume_unique=True).size)
    manifest=json.loads((a.chunks/'manifest.json').read_text())
    for k,z in current_redshift_chunks(a.assets,manifest,a.index):
        ids=np.load(a.index/f'targetid_{k:03d}.npy',mmap_mode='r')
        for f,lo,hi in FIELDS:
            ix=np.flatnonzero((z>=lo)&(z<hi))
            if not len(ix): continue
            q=ids[ix]
            for r in ('NGC','SGC'):
                cached=field_cache[(f,r)]; at=np.searchsorted(cached,q)
                hit=(at<len(cached))&(cached[np.minimum(at,len(cached)-1)]==q)
                output[f][r]['exactCatalogMatches']+=int(hit.sum())
    report=json.loads(a.report.read_text())
    report['exactCacheMatchesByRegion']=output
    report['withinFieldNgsSgcCacheTargetIdIntersections']=overlaps
    report['exactCacheMatchesByRegionInterpretation']='Exact TARGETID matches by current viewer redshift selection and hemisphere cache; NGC/SGC overlap is checked directly.'
    # Confirm repeated identities inside the combined LRG+ELG input are
    # duplicate sky rows, rather than distinct/conflicting positions.
    from astropy.cosmology import Planck18
    sources=local_sources(a.data_dir)
    consistency={}
    for f,lo,hi in FIELDS:
        for r in ('NGC','SGC'):
            source=next(Path(s['localPath']) for s in sources if s['tracer']==f and s['region']==r and s['kind']=='data')
            rows=_fits_rows(source,Planck18,0,lo,hi)
            order=np.argsort(rows['targetid'],kind='stable'); tid=rows['targetid'][order]; xyz=rows['position'][order]
            starts=np.r_[0,np.flatnonzero(tid[1:]!=tid[:-1])+1]; ends=np.r_[starts[1:],len(tid)]
            duplicate_groups=[(int(s),int(e)) for s,e in zip(starts,ends) if e-s>1]
            max_sep=0.0; max_dz=0.0; n_position_conflicts=0; n_z_spread=0
            from astropy.io import fits
            with fits.open(source,memmap=True,lazy_load_hdus=True) as hdul:
                table=next(h.data for h in hdul if getattr(h.data,'names',None))
                step=1; row=table[::step]
                w=np.asarray(row['WEIGHT'],dtype='f8')
                if 'WEIGHT_FKP' in row.names: w*=np.asarray(row['WEIGHT_FKP'],dtype='f8')
                rz=np.asarray(row['Z'],dtype='f8'); rra=np.asarray(row['RA'],dtype='f8'); rdec=np.asarray(row['DEC'],dtype='f8')
                keep=(np.isfinite(w)&(w>0)&np.isfinite(rz)&np.isfinite(rra)&np.isfinite(rdec)&(rz>=lo)&(rz<hi))
                raw_tid=np.asarray(row['TARGETID'][keep],dtype='u8'); raw_z=rz[keep]; raw_ra=rra[keep]; raw_dec=rdec[keep]
            raw_order=np.argsort(raw_tid,kind='stable'); raw_tid=raw_tid[raw_order]; raw_z=raw_z[raw_order]; raw_ra=raw_ra[raw_order]; raw_dec=raw_dec[raw_order]
            for s,e in duplicate_groups:
                group_sep=float(np.max(np.linalg.norm(xyz[s:e]-xyz[s],axis=1)))
                max_sep=max(max_sep,group_sep)
                if group_sep>1.e-5: n_position_conflicts+=1
                zspan=float(np.ptp(raw_z[s:e])); max_dz=max(max_dz,zspan)
                if zspan>1.e-8: n_z_spread+=1
            consistency[f'{f}_{r}']={'duplicateTargetIdGroups':len(duplicate_groups),
                'duplicateRowsBeyondFirst':int(sum(e-s-1 for s,e in duplicate_groups)),
                'maxWithinGroup3dSeparationMpc':max_sep,
                'groupsWith3dPositionConflictsAbove1e-5Mpc':n_position_conflicts,
                'maxWithinGroupRedshiftSpan':max_dz,'groupsWithRedshiftSpanAbove1e-8':n_z_spread,
                'meaning':'Duplicate TARGETIDs are deduplicated by highest total weight; position/redshift conflicts are quantified separately.'}
    report['inputClusteringDataIdentityAudit']['duplicateSkyRowConsistency']=consistency
    a.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'exactCacheMatchesByRegion':output,'withinFieldCacheIntersections':overlaps},indent=2))
if __name__=='__main__': main()
