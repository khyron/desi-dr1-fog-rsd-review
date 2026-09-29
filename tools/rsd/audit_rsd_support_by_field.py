#!/usr/bin/env python3
"""Full-catalogue, random-cloud support audit for all DR1 RSD tracer/regions.

This estimates where an in-memory field query could be attempted, not where a
correction is scientifically validated. No shifts are reconstructed or written.
"""
from __future__ import annotations
import argparse, json, resource, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
sys.path.insert(0, str(ROOT))
from rsd_pipeline import _fits_rows, _deduplicate_targetids, _source_path, local_sources, config, read_desi_header
from audit_rsd_catalog_coverage import FIELDS, read_rsd
from build_desi_chunks import morton_order

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir',type=Path,default=Path('/mnt/d/desi/data_external/RSD'))
    ap.add_argument('--assets',type=Path,default=ROOT/'assets')
    ap.add_argument('--chunks',type=Path,default=ROOT/'ply/chunks')
    ap.add_argument('--index',type=Path,default=ROOT/'ply/reconstruction_dr1/private_index')
    ap.add_argument('--cache-dir',type=Path,default=Path('/mnt/d/desi/data_external/RSD/derived/rsd-v1-exploratory'))
    ap.add_argument('--config',type=Path,default=ROOT/'tools/rsd/rsd_config.json')
    ap.add_argument('--out',type=Path,default=ROOT/'tools/rsd/rsd-field-support-audit.json')
    ap.add_argument('--validation-sample',type=int,default=250000)
    ap.add_argument('--support-percentile',type=float,default=99.0)
    ap.add_argument('--random-limit',type=int,default=1500000,
                    help='deterministic stride sample from each official random FITS; 0 reads all rows')
    ap.add_argument('--viewer-sample-per-field',type=int,default=100000,
                    help='uniform sample of eligible DR1 points per field, expanded with Wilson intervals')
    ap.add_argument('--field',choices=[x[0] for x in FIELDS],action='append',help='repeat to restrict fields')
    ap.add_argument('--seed',type=int,default=20260923)
    a=ap.parse_args(); t0=time.time()
    from astropy.cosmology import Planck18
    sources=local_sources(a.data_dir); cfg=config(a.config)
    manifest=json.loads((a.chunks/'manifest.json').read_text())
    n=int(manifest['count']); units=float(manifest['unitsPerMpc'])
    # Recreate exact current viewer partition order, so each point is paired to
    # the corresponding TARGETID index entry, while retaining xyz for support.
    all_names=list(manifest['tracers']); xyz_parts=[]; z_parts=[]; id_parts=[]
    for name in all_names:
        p=a.assets/f'{name}.desi'; h=read_desi_header(p); c=h['count']
        if h['version']!=2: raise ValueError(f'{p}: expected v2')
        raw=np.memmap(p,dtype='<i2',mode='r',offset=48,shape=(c,3))
        xyz_parts.append(np.rint(raw.astype('f4')/h['unitsPerMpc']*units).astype('<i2'))
        q=np.memmap(p,dtype='<u2',mode='r',offset=48+c*6,shape=(c,))
        z_parts.append(h['zMin']+q.astype('f4')*((h['zMax']-h['zMin'])/65535.0))
    xyz=np.concatenate(xyz_parts).astype('f4'); z=np.concatenate(z_parts).astype('f4')
    del xyz_parts,z_parts
    hashes=np.abs(xyz[:,0].astype('i8')*17+xyz[:,2].astype('i8')*59+xyz[:,1].astype('i8')*101)
    ids=np.empty(n,dtype='u8'); pos=np.empty((n,3),dtype='f4'); zz=np.empty(n,dtype='f4')
    cursor=0
    for k,ch in enumerate(manifest['chunks']):
        ix=np.flatnonzero(hashes%len(manifest['chunks'])==k); ix=ix[morton_order(xyz[ix])]
        ti=np.load(a.index/f'targetid_{k:03d}.npy',mmap_mode='r')
        if len(ix)!=len(ti) or len(ti)!=ch['count']: raise ValueError(f'chunk {k}: index mismatch')
        lo=cursor; hi=lo+len(ix); cursor=hi
        pos[lo:hi]=xyz[ix].astype('f4')/units; zz[lo:hi]=z[ix]; ids[lo:hi]=ti
    del xyz,z,hashes
    if len(ids)!=n: raise ValueError('catalog count mismatch')
    cache=read_rsd(a.cache_dir)
    config_fields={name:(float(p['zMin']),float(p['zMax'])) for name,p in cfg['tracers'].items()}
    # Audit duplicate identities in input data selections, including the
    # combined LRG+ELG catalogue and NGC/SGC overlap.
    data_ids={}; data_diag={}
    for f,lo,hi in FIELDS:
        for region in ('NGC','SGC'):
            path=_source_path(sources,f,region,'data')
            rows=_fits_rows(path,Planck18,0,lo,hi); dedup=_deduplicate_targetids(rows)
            key=f'{f}_{region}'; data_ids[key]=dedup['targetid']
            data_diag[key]={'rowsAfterZAndWeight':int(len(rows['targetid'])),
                'duplicatesDiscardedByTargetID':int(dedup.get('duplicatesDiscarded',0)),
                'rowsAfterIdentityDedup':int(len(dedup['targetid']))}
    data_overlaps={}
    for x,y in [('LRG_ELG_NGC','LRG_NGC'),('LRG_ELG_NGC','ELG_NGC'),
                ('LRG_ELG_SGC','LRG_SGC'),('LRG_ELG_SGC','ELG_SGC'),
                ('BGS_NGC','BGS_SGC'),('LRG_NGC','LRG_SGC'),('LRG_ELG_NGC','LRG_ELG_SGC'),
                ('ELG_NGC','ELG_SGC'),('QSO_NGC','QSO_SGC')]:
        data_overlaps[f'{x}__{y}']=int(np.intersect1d(data_ids[x],data_ids[y],assume_unique=True).size)
    by={}; rng=np.random.default_rng(a.seed)
    for field,lo,hi in FIELDS:
        if a.field and field not in a.field: continue
        field_t=time.time()
        eligible=(zz>=lo)&(zz<hi); eligible_ix=np.flatnonzero(eligible); nelig=len(eligible_ix)
        fids, _=cache[field]
        at=np.searchsorted(fids,ids[eligible_ix]); hit=(at<len(fids))&(fids[np.minimum(at,len(fids)-1)]==ids[eligible_ix])
        cached_exact=int(hit.sum())
        ns=min(a.viewer_sample_per_field,nelig)
        ix=np.sort(rng.choice(eligible_ix,size=ns,replace=False)) if ns<nelig else eligible_ix
        at_s=np.searchsorted(fids,ids[ix]); cache_hit=(at_s<len(fids))&(fids[np.minimum(at_s,len(fids)-1)]==ids[ix])
        # For each hemisphere, calibrate 3-D distance threshold on a disjoint
        # half split of the official z/weight-filtered random catalogue, then
        # query every eligible viewer point against the other half.
        regions={}; support_by_region=[]; box_by_region=[]
        for region in ('NGC','SGC'):
            reg_t=time.time()
            data_path=_source_path(sources,field,region,'data'); rand_path=_source_path(sources,field,region,'random')
            print(f'[{field} {region}] reading and filtering randoms: {rand_path.name}',flush=True)
            rows=_deduplicate_targetids(_fits_rows(rand_path,Planck18,a.random_limit,lo,hi))
            rp=np.asarray(rows['position'],dtype='f4'); nrand=len(rp)
            print(f'[{field} {region}] filtered randoms={nrand:,}; building support tree',flush=True)
            if nrand<4: raise ValueError(f'{field} {region}: insufficient filtered randoms')
            # Deterministic independent subsets without loading two copies.
            train=rp[::2]; valid=rp[1::2]
            nv=min(a.validation_sample,len(valid)); valid=valid[:nv]
            tree=cKDTree(train,compact_nodes=True,balanced_tree=True)
            vd,_=tree.query(valid,k=1,workers=-1)
            radius=float(np.percentile(vd,a.support_percentile))
            d,_=tree.query(pos[ix],k=1,workers=-1)
            report_path=a.cache_dir/f'{field}_{region}_report.json'; report=json.loads(report_path.read_text())
            mesh=report['mesh']; center=np.asarray(mesh['boxcenterMpc']); half=np.asarray(mesh['boxsizeMpc'])/2
            in_box=np.all(np.abs(pos[ix]-center)<=half,axis=1)
            supported=(d<=radius)&in_box
            support_by_region.append(supported); box_by_region.append(in_box)
            regions[region]={'officialRandomRowsAfterZWeightAndIDDedup':nrand,
                'randomStride':int(rows['stride']),
                'calibrationTrainRows':int(len(train)),
                'independentValidationRows':int(nv),'supportRadiusMpc':radius,
                'validationAcceptedFraction':float(np.mean(vd<=radius)),
                'viewerInsideReportedMesh':int(in_box.sum()),'viewerRandomSupportedAndInMesh':int(supported.sum()),
                'meshSource':'existing reconstruction report bounds; actual field mesh is not serialized here'}
            del tree,rp,rows,train,valid,vd,d
            print(f'[{field} {region}] support radius={radius:.3f} Mpc; completed in {time.time()-reg_t:.1f}s',flush=True)
        both=support_by_region[0]&support_by_region[1]
        any_support=support_by_region[0]|support_by_region[1]
        # Require a unique region to call a spatial query candidate; ambiguous
        # overlap is reported separately and is not counted as new candidate.
        unique=any_support&~both
        new=unique&~cache_hit
        no_support=~any_support&~cache_hit
        ambiguous=both&~cache_hit
        cache_n=cached_exact
        uncached_pop=nelig-cache_n
        def expand(mask):
            denom=int((~cache_hit).sum()); successes=int((mask&~cache_hit).sum())
            if denom==0: return {'sampleCount':0,'estimatedCatalogCount':0,'ci95':[0,0]}
            p=successes/denom; zc=1.959963984540054; den=1+zc*zc/denom
            center=(p+zc*zc/(2*denom))/den
            half=zc*np.sqrt(p*(1-p)/denom+zc*zc/(4*denom*denom))/den
            return {'sampleCount':successes,'estimatedCatalogCount':round(p*uncached_pop),
                'ci95':[round(max(0,center-half)*uncached_pop),round(min(1,center+half)*uncached_pop)],'sampleFraction':p}
        # Region attribution is limited to uniquely supported sample points;
        # each region expansion has its own sample count and Wilson interval.
        region_estimates={}
        for ri,region in enumerate(('NGC','SGC')):
            reg_unique=support_by_region[ri]&~support_by_region[1-ri]
            region_estimates[region]=expand(reg_unique)
        by[field]={'zRange':[lo,hi],'eligible':nelig,'cachedExactTargets':cache_n,
            'sampledEligibleRows':int(ns),'uncachedPopulationForExpansion':uncached_pop,
            'estimatedNewTargetsWithUniqueRandomSupportProxy':expand(unique),
            'estimatedUnmatchedWithoutSupportProxy':expand(~any_support),
            'estimatedUnmatchedRegionAmbiguous':expand(both),
            'estimatedNewTargetsByUniqueRegion':region_estimates,
            'sampledSupportedUniqueUncached':int(new.sum()),'sampledUnsupportedUncached':int(no_support.sum()),
            'sampledSupportedButRegionAmbiguousUncached':int(ambiguous.sum()),
            'sampledSupportedIncludingCache':int(any_support.sum()),
            'regions':regions,
            'cacheOnlyFieldRows':int(len(fids)),
            'interpretation':'candidate support diagnostic only; not validated corrections'}
        print(f'{field}: eligible={nelig:,} cacheExact={cache_n:,} sample={ns:,} newProxySample={new.sum():,} noSupportSample={no_support.sum():,} ambiguousSample={ambiguous.sum():,}; elapsed={time.time()-field_t:.1f}s',flush=True)
    # Detect duplicate identities between independent cache fields; LRG_ELG
    # combined selection gets explicit comparison against neighboring tracers.
    cache_ids={}
    for f in cache:
        cache_ids[f]=cache[f][0]
    overlaps={}
    names=list(cache_ids)
    for i,x in enumerate(names):
        for y in names[i+1:]:
            overlaps[f'{x}__{y}']=int(np.intersect1d(cache_ids[x],cache_ids[y],assume_unique=True).size)
    report={'catalog':'DESI DR1 current viewer catalog; exact-validated identity partition',
      'objects':n,'insideEnvelope':sum(x['eligible'] for x in by.values()),
      'outsideEnvelope':n-sum(x['eligible'] for x in by.values()),'byField':by,
      'cacheTargetIdIntersectionsBetweenFields':overlaps,
      'inputClusteringDataIdentityAudit':{'perCatalog':data_diag,'intersections':data_overlaps},
      'method':{'randomCatalog':'official local v1.5 random index 0, exact z envelope, finite positive WEIGHT*(WEIGHT_FKP if available), TARGETID deduplicated',
        'support':'nearest 3-D random distance <= p99 of disjoint random validation distances and within existing report mesh bounds',
        'regionPolicy':'NGC and SGC must give a unique supported region; both-supported points are reported ambiguous and excluded from new candidates',
        'cache':'exact TARGETID matches across existing NGC/SGC cache files',
        'notPerformed':'No reconstruction, no new field query, no deltas written; reported mesh bounds do not prove mesh mask validity'},
      'limitations':['Random-cloud proximity is a proxy, not DESI angular-mask polygon or formal volume mask.',
        'A supported new target would still require live-field query, holdout validation, and scientific review.',
        'Existing cached exact matches are counted even when region cache intersections occur; overlaps reported separately.',
        'LRG_ELG is its own configured selection; cache identity intersections are explicitly audited.'],
      'sampling':{'uniformWithoutReplacementPerField':True,'viewerSamplePerField':a.viewer_sample_per_field,
        'officialRandomStrideLimit':a.random_limit,'wilsonIntervalZ':1.959963984540054,
        'interpretation':('candidate counts are design-based viewer-sample expansions with 95% Wilson intervals; exact cache counts are from all catalog IDs; all z/weight-selected official random rows were used for support calibration' if a.random_limit==0 else 'candidate counts are design-based viewer-sample expansions with 95% Wilson intervals; exact cache counts are from all catalog IDs; random support uses a deterministic stride sample of official random rows')},
      'elapsedSeconds':time.time()-t0,'peakResidentMemoryGiB':float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2)}
    a.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8'); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
