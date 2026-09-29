"""Build checked, single-download chunks from the unchanged published catalogue."""
import gzip, hashlib, json, struct, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from build_desi_chunks import morton_order
OUT = ROOT / 'desiV3/data'
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source_manifest = ROOT/'ply/chunks/manifest.json'
    old = json.loads(source_manifest.read_text())
    catalogs, positions, redshifts, attributes, flags, ids = [], [], [], [], [], []
    inputs = []
    for code, name in enumerate(old['tracers']):
        path = ROOT/'assets'/f'{name}.desi'
        raw = path.read_bytes()
        version, n, units, zmin, zmax = struct.unpack_from('<IIfff', raw, 4)
        assert raw[:4] == b'DESI' and version == 2 and len(raw) == 48+12*n
        catalogs.append(dict(name=name+'.desi',count=n,unitsPerMpc=units,zMin=zmin,zMax=zmax,headerHex=raw[:48].hex()))
        positions.append(np.frombuffer(raw,'<i2',3*n,48).reshape(n,3))
        redshifts.append(np.frombuffer(raw,'<u2',n,48+6*n))
        attributes.append(np.frombuffer(raw,'u1',4*n,48+8*n).reshape(4,n).T)
        fp=ROOT/'assets'/f'{name}.target-flags.bin'
        f=np.frombuffer(fp.read_bytes(),'u1'); assert len(f)==n
        flags.append(f | (code << 6))
        ids.append(np.load(ROOT/'ply/reconstruction_dr1/private_index'/f'{name}_targetid.npy',mmap_mode='r'))
        inputs.extend([dict(file=str(path.relative_to(ROOT)),sha256=sha(path)),dict(file=str(fp.relative_to(ROOT)),sha256=sha(fp))])
    pos=np.concatenate(positions); zq=np.concatenate(redshifts); attrs=np.concatenate(attributes)
    types=np.concatenate(flags); target=np.concatenate(ids)
    rows=np.concatenate([np.arange(c['count'],dtype='<u4') for c in catalogs])
    globalpos=np.concatenate([np.rint(p.astype('f4')/c['unitsPerMpc']*old['unitsPerMpc']).astype('<i2') for p,c in zip(positions,catalogs)])
    hashes=np.abs(globalpos[:,0].astype('i8')*17+globalpos[:,2].astype('i8')*59+globalpos[:,1].astype('i8')*101)
    manifest={**old,'encoding':'desi-consolidated-1','catalogs':catalogs,'chunks':[]}
    if manifest.get('localPrecision'):
        manifest['localPrecision'] = dict(manifest['localPrecision'])
        if not manifest['localPrecision']['file'].startswith(('https://', 'asset:')):
            manifest['localPrecision']['file'] = 'https://www.tinyplanet.cl/desi/desi-chunks/' + manifest['localPrecision']['file']
    report=[]; seen=np.zeros(old['count'],dtype='u1')
    for i,spec in enumerate(old['chunks']):
        chosen=np.flatnonzero(hashes%len(old['chunks'])==i)
        chosen=chosen[morton_order(globalpos[chosen])]; n=len(chosen)
        assert n==spec['count'] and not np.any(seen[chosen]); seen[chosen]=1
        reference=np.load(ROOT/'ply/reconstruction_dr1/private_index'/f'targetid_{i:03d}.npy')
        assert np.array_equal(target[chosen],reference), 'Object row mapping changed'
        deltas=np.diff(globalpos[chosen].astype('i4'),axis=0,prepend=0).astype('<i2')
        assert deltas.T.copy().tobytes()==gzip.decompress((ROOT/'ply/chunks'/spec['file']).read_bytes())
        assert types[chosen].tobytes()==gzip.decompress((ROOT/'ply/chunks'/spec['type']).read_bytes())
        native=np.diff(pos[chosen].astype('i4'),axis=0,prepend=0).astype('<i2').T.copy()
        sections=[types[chosen].tobytes(),rows[chosen].tobytes(),native.tobytes(),zq[chosen].tobytes(),attrs[chosen].T.copy().tobytes()]
        raw=bytearray(64); offsets=[]
        for section in sections:
            raw.extend(b'\0'*((-len(raw))%4)); offsets.append(len(raw)); raw.extend(section)
        struct.pack_into('<4s8I',raw,0,b'DSC3',1,n,*offsets,len(raw))
        # Independently decode every section and compare to the native source rows.
        p=np.frombuffer(raw,'<i2',3*n,offsets[2]).reshape(3,n)
        decoded=np.cumsum(p.astype('i4'),axis=1).astype('<i2').T
        assert np.array_equal(decoded,pos[chosen])
        assert np.array_equal(np.frombuffer(raw,'<u4',n,offsets[1]),rows[chosen])
        assert np.array_equal(np.frombuffer(raw,'<u2',n,offsets[3]),zq[chosen])
        assert np.array_equal(np.frombuffer(raw,'u1',4*n,offsets[4]).reshape(4,n).T,attrs[chosen])
        blob=gzip.compress(raw,compresslevel=6,mtime=0)
        name=f'catalog_{i:03d}.bin.gz'; (OUT/name).write_bytes(blob)
        manifest['chunks'].append(dict(file=name,count=n,bytes=len(blob)))
        report.append(dict(chunk=i,count=n,sha256=hashlib.sha256(blob).hexdigest(),geometryExact=True,typesExact=True,nativeMetadataExact=True,targetOrderExact=True))
        print(f'{i}: {n:,} rows, {len(blob):,} bytes, exact',flush=True)
    assert np.all(seen) and sum(c['count'] for c in catalogs)==old['count']
    baseline=sum(c['bytes']+c['typeBytes'] for c in old['chunks'])+sum((ROOT/'assets'/c['name']).stat().st_size+c['count'] for c in catalogs)
    total=sum(c['bytes'] for c in manifest['chunks'])
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'validation.json').write_text(json.dumps(dict(status='passed',count=old['count'],inputs=inputs,chunks=report,previousTransferBytes=baseline,newTransferBytes=total,savedBytes=baseline-total),indent=2)+'\n')
if __name__=='__main__': main()
