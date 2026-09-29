"""Fetch missing exact TARGETIDs from the public DESI DR1 SQL service.

Interval queries avoid huge ADQL IN lists. Retain only exact requested IDs.
"""
from pathlib import Path
import csv, gzip, hashlib, io, json, struct, urllib.parse, urllib.request, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'desiV3/data/missing-photometry-sql'
OUT.mkdir(exist_ok=True)
ENDPOINT = 'https://datalab.noirlab.edu/query/query'
ids = [np.load(ROOT / 'ply/reconstruction_dr1/private_index' / (name + '_targetid.npy')) for name in ['DR1_GALAXY', 'DR1_QSO']]
missing = []
for i in range(12):
    raw = gzip.decompress((ROOT / 'desiV3/data/morphology' / f'catalog_{i:03d}.bin.gz').read_bytes())
    version, n, ot, ori, op, oz, oa, total = struct.unpack_from('<8I', raw, 4)
    qso = (np.frombuffer(raw, 'u1', n, ot) & 64) != 0
    rows = np.frombuffer(raw, '<u4', n, ori)
    mask = np.frombuffer(raw, 'u1', n, oa) == 255
    for code in [0, 1]: missing.append(ids[code][rows[mask & (qso == code)]])
keys = np.unique(np.concatenate(missing))
batch = 50000
blocks = [keys[i:i+batch] for i in range(0, len(keys), batch)]

def process(item):
    number, targets = item
    name = f'tap-{number:04d}'
    compact, receipt = OUT / (name + '.npz'), OUT / (name + '.json')
    condition = (f'p.targetid BETWEEN {int(targets[0])} AND {int(targets[-1])} '
        "AND EXISTS (SELECT 1 FROM desi_dr1.zpix z WHERE z.targetid=p.targetid "
        "AND z.zcat_primary AND z.zwarn=0 AND z.spectype IN ('GALAXY','QSO'))")
    query = ('SELECT p.targetid,p.morphtype,p.shape_r,p.flux_g,p.flux_r,p.flux_z FROM desi_dr1.photometry p WHERE ' + condition)
    if receipt.exists() and compact.exists():
        r = json.loads(receipt.read_text())
        if r['sourceQuery'] == query and hashlib.sha256(compact.read_bytes()).hexdigest() == r['compactSHA256']:
            return r
    def request(sql):
        url = ENDPOINT + '?' + urllib.parse.urlencode(dict(sql=sql, ofmt='csv', out='', **{'async':'False'}))
        headers = {'X-DL-AuthToken':'anonymous.0.0.anon_access', 'X-DL-TimeoutRequest':'180'}
        return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=240)
    print('Query TAP', number+1, '/', len(blocks), flush=True)
    for attempt in range(4):
        try:
            with request(query) as response:
                raw = response.read()
            reader = csv.DictReader(io.StringIO(raw.decode()))
            assert reader.fieldnames == ['targetid','morphtype','shape_r','flux_g','flux_r','flux_z'], raw[:300]
            rows = list(reader)
            tid = np.array([int(r['targetid']) for r in rows], dtype='u8')
            assert len(np.unique(tid)) == len(tid)
            assert ((tid >= targets[0]) & (tid <= targets[-1])).all()
            count_query = 'SELECT COUNT(*) AS n FROM desi_dr1.photometry p WHERE ' + condition
            with request(count_query) as response:
                count_raw = response.read()
            expected = int(list(csv.DictReader(io.StringIO(count_raw.decode())))[0]['n'])
            assert len(rows) == expected, 'Incomplete TAP result'
            break
        except Exception as exc:
            print('Retry TAP', number+1, type(exc).__name__, str(exc)[:160], flush=True)
            if attempt == 3: raise
            time.sleep(2 * (attempt + 1))
    total_returned = len(rows)
    hit = np.isin(tid, targets)
    rows = [r for r, keep in zip(rows, hit) if keep]
    arrays = dict(TARGETID=tid[hit], TYPE=np.array([r['morphtype'] for r in rows], dtype='U4'))
    for col in ['SHAPE_R','FLUX_G','FLUX_R','FLUX_Z']:
        arrays[col] = np.array([float(r[col.lower()]) if r[col.lower()] else np.nan for r in rows], dtype='f4')
    part = OUT / (name + '.part')
    with part.open('wb') as stream: np.savez_compressed(stream, **arrays)
    part.replace(compact)
    (OUT / (name + '.csv.gz')).write_bytes(gzip.compress(raw, mtime=0))
    r = dict(sourceURL=ENDPOINT, sourceQuery=query, sourceFile=name, sourcePixel=-1,
        sourceSHA256=hashlib.sha256(raw).hexdigest(), sourceBytes=len(raw), rows=len(rows),
        requestedRows=len(targets), independentCount=expected, returnedRows=total_returned,
        compactSHA256=hashlib.sha256(compact.read_bytes()).hexdigest())
    receipt.write_text(json.dumps(r, indent=2))
    print('Verified TAP', number+1, '/', len(blocks), 'found', len(rows), 'of', len(targets), flush=True)
    return r

records = []
with ThreadPoolExecutor(max_workers=4) as pool:
    futures = [pool.submit(process, item) for item in enumerate(blocks)]
    for future in as_completed(futures):
        records.append(future.result())
        (OUT / 'download-progress.json').write_text(json.dumps(dict(complete=len(records),total=len(blocks),requestedUniqueIDs=len(keys),files=records),indent=2))
print('COMPLETE', len(records), 'batches', flush=True)
