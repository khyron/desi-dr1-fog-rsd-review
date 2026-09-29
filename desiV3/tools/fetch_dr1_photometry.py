"""Download verified official observed-target Tractor files, retaining compact columns.

Run with the workspace vendor directory on PYTHONPATH. Resumable per sky pixel.
Only the downloader's own temporary FITS file is removed after checksum verification
and atomic compact extraction; the source remains publicly reproducible.
"""
from pathlib import Path
import hashlib, json, re, urllib.request, shutil, time, gc, argparse
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'desiV3/data/official-photometry'
OUT.mkdir(parents=True, exist_ok=True)
BASE = 'https://data.desi.lbl.gov/public/dr1/vac/dr1/lsdr9-photometry/iron/v1.1/observed-targets/tractorphot/'

def fetch(url):
    for attempt in range(4):
        try:
            return urllib.request.urlopen(url, timeout=120)
        except Exception:
            if attempt == 3: raise
            time.sleep(2 * (attempt + 1))

index = fetch(BASE).read().decode()
checksum_name = re.search(r'href="([^"]+sha256sum)"', index).group(1)
checksums = fetch(BASE + checksum_name).read().decode()
(OUT / 'source.sha256sum').write_text(checksums)
expected = {line.split()[-1].split('/')[-1]: line.split()[0] for line in checksums.splitlines() if line.strip()}
names = sorted(set(re.findall(r'href="(tractorphot[^"/]+\.fits)"', index)))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pixel', type=int, help='Download only one audit pixel; omit for all official files')
options = parser.parse_args()
if options.pixel is not None:
    selected = f'tractorphot-nside4-hp{options.pixel:03d}-iron.fits'
    assert selected in names, 'Pixel not present in official manifest'
    names = [selected]
def process(item):
    number, name = item
    compact = OUT / (name + '.npz')
    receipt = OUT / (name + '.json')
    if compact.exists() and receipt.exists():
        r = json.loads(receipt.read_text())
        if hashlib.sha256(compact.read_bytes()).hexdigest() == r['compactSHA256']:
            (OUT / (name + '.download')).unlink(missing_ok=True)
            print('Cached', number + 1, name, flush=True)
            return r
    temp = OUT / (name + '.download')
    print('Downloading', number + 1, '/', len(names), name, flush=True)
    size = int(re.search(r'href="' + re.escape(name) + r'"[^\n]+?\s(\d+)\s*$', index, re.M).group(1))
    if shutil.disk_usage(OUT).free < size * 8 + 1024**3:
        raise RuntimeError('Insufficient free space for next temporary FITS and extraction')
    # The public server may close a large transfer early. Resume verified byte ranges.
    failures = 0
    while not temp.exists() or temp.stat().st_size < size:
        offset = temp.stat().st_size if temp.exists() else 0
        request = urllib.request.Request(BASE + name, headers={'Range': f'bytes={offset}-'})
        try:
            with fetch(request) as response:
                if offset:
                    assert response.status == 206 and response.headers['Content-Range'].startswith(f'bytes {offset}-'), 'Server ignored resume range'
                with temp.open('ab') as stream:
                    remaining = size - offset
                    while remaining:
                        block = response.read(min(4 * 1024**2, remaining))
                        if not block: break
                        stream.write(block)
                        remaining -= len(block)
            assert temp.stat().st_size > offset, 'Transfer made no progress'
            failures = 0
        except Exception:
            failures += 1
            if failures >= 4: raise
            time.sleep(2 * failures)
    assert temp.stat().st_size == size
    digest = hashlib.sha256()
    with temp.open('rb') as stream:
        while block := stream.read(4 * 1024**2): digest.update(block)
    assert digest.hexdigest() == expected[name], 'Official SHA256 mismatch: ' + name
    with fits.open(temp, memmap=True) as hdus:
        table = hdus[1].data
        columns = ['TARGETID', 'TYPE', 'SHAPE_R', 'FLUX_G', 'FLUX_R', 'FLUX_Z']
        arrays = {col: np.array(table[col]) for col in columns}
        # Preserve measured floats, including nonpositive fluxes and invalid values.
        # No proxy classification or colour defaults are fabricated here.
        part = OUT / (name + '.npz.part')
        with part.open('wb') as stream:
            np.savez_compressed(stream, **arrays)
        part.replace(compact)
        rows = len(table)
        del table
    del hdus
    gc.collect()
    r = dict(sourceURL=BASE + name, sourceSHA256=digest.hexdigest(),
             sourceBytes=temp.stat().st_size, rows=rows, columns=columns,
             compactSHA256=hashlib.sha256(compact.read_bytes()).hexdigest())
    receipt.write_text(json.dumps(r, indent=2))
    temp.unlink()
    print('Verified and extracted', rows, 'rows;', compact.stat().st_size, 'bytes', flush=True)
    return r

records = []
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = [pool.submit(process, item) for item in enumerate(names)]
    for future in as_completed(futures):
        records.append(future.result())
        (OUT / 'download-progress.json').write_text(json.dumps(dict(
            complete=len(records), total=len(names), files=records), indent=2))
(OUT / 'download-progress.json').write_text(json.dumps(dict(complete=len(records), total=len(names), files=records), indent=2))
print('COMPLETE', len(records), 'files', flush=True)
