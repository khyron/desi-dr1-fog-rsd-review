"""Recover imaging attributes by exact TARGETID from locally retained catalogues."""
from pathlib import Path
import struct, gzip, json
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'desiV3/data/morphology'
OUT.mkdir(exist_ok=True)
ids, attrs = [], []
for name in ['BGS_BRIGHT','LRG','ELG','QSO']:
    raw = (ROOT/'assets'/f'{name}.desi').read_bytes()
    version,n = struct.unpack_from('<II',raw,4)
    assert version == 3 and len(raw) == 48+32*n
    ids.append(np.frombuffer(raw,'<u8',n,48+24*n).copy())
    attrs.append(np.frombuffer(raw,'u1',4*n,48+8*n).reshape(4,n).T.copy())
target = np.concatenate(ids); photo = np.concatenate(attrs)
order = np.argsort(target,kind='stable'); target=target[order]; photo=photo[order]
unique, first, counts = np.unique(target,return_index=True,return_counts=True)
# A conflicting imaging fit is not a reliable classification: leave it unknown.
conflicts = np.zeros(len(unique),bool)
for axis in range(4):
    lo=np.minimum.reduceat(photo[:,axis],first); hi=np.maximum.reduceat(photo[:,axis],first)
    conflicts |= lo != hi
photo=photo[first]; unique=unique[~conflicts]; photo=photo[~conflicts]
print('Reliable IDs',len(unique),'conflicting IDs excluded',int(conflicts.sum()),flush=True)
rowids=[np.load(ROOT/'ply/reconstruction_dr1/private_index'/f'{name}_targetid.npy',mmap_mode='r') for name in ['DR1_GALAXY','DR1_QSO']]
report=[]
for i in range(12):
    path=ROOT/'desiV3/data'/f'catalog_{i:03d}.bin.gz'
    raw=bytearray(gzip.decompress(path.read_bytes()))
    version,n,o_type,o_row,o_pos,o_z,o_attrs,total=struct.unpack_from('<8I',raw,4)
    rows=np.frombuffer(raw,'<u4',n,o_row); types=np.frombuffer(raw,'u1',n,o_type)
    names=(types&64)!=0; keys=np.empty(n,dtype='<u8')
    for code in [0,1]:
        mask=names==code; keys[mask]=rowids[code][rows[mask]]
    idx=np.searchsorted(unique,keys); safe=np.minimum(idx,len(unique)-1)
    match=(idx<len(unique))&(unique[safe]==keys)
    before=bytes(raw)
    view=np.frombuffer(raw,'u1',4*n,o_attrs).reshape(4,n)
    view[:,match]=photo[safe[match]].T
    assert raw[:o_attrs]==before[:o_attrs] and raw[o_attrs+4*n:]==before[o_attrs+4*n:]
    blob=gzip.compress(raw,compresslevel=6,mtime=0)
    (OUT/path.name).write_bytes(blob)
    counts=dict(zip(*[a.tolist() for a in np.unique(view[0],return_counts=True)]))
    entry=dict(chunk=i,count=n,matched=int(match.sum()),galaxyMatched=int((match&~names).sum()),qsoMatched=int((match&names).sum()),morphologyCounts=counts,bytes=len(blob),geometryAndIdentityUnchanged=True)
    report.append(entry); print(entry,flush=True)
(OUT/'validation.json').write_text(json.dumps(dict(method='Exact TARGETID; ambiguous imaging tuples excluded',excludedConflicts=int(conflicts.sum()),chunks=report),indent=2))
