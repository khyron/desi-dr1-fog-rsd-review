#!/usr/bin/env python3
"""Download selected official DR1 v1.5 BGS files with verified resumable ranges."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import time
from pathlib import Path

import requests

BASE = 'https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/'
DEST = Path(__file__).resolve().parent / 'lowz_sources'
ALLOWED = {
    'BGS_BRIGHT_NGC_clustering.dat.fits',
    'BGS_BRIGHT_SGC_clustering.dat.fits',
    'BGS_BRIGHT_NGC_0_clustering.ran.fits',
    'BGS_BRIGHT_SGC_0_clustering.ran.fits',
    'BGS_BRIGHT_NGC_1_clustering.ran.fits',
    'BGS_BRIGHT_SGC_1_clustering.ran.fits',
    'BGS_BRIGHT_full.dat.fits',
    'BGS_BRIGHT_0_full.ran.fits',
}
PART_BYTES = 32 * 1024 * 1024


def get_part(url, total, part_dir, index, start, end):
    path = part_dir / f'{index:04d}.part'
    expected = end - start + 1
    if path.exists() and path.stat().st_size == expected:
        return
    for attempt in range(12):
        try:
            with requests.get(url, headers={'Range': f'bytes={start}-{end}'}, stream=True,
                              timeout=(30, 120)) as response:
                response.raise_for_status()
                if response.status_code != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{total}':
                    raise ValueError(f'bad range response for {index}')
                with path.open('wb') as stream:
                    for block in response.iter_content(1024 * 1024):
                        if block:
                            stream.write(block)
            if path.stat().st_size != expected:
                raise ValueError(f'part {index}: wrong byte count')
            return
        except (OSError, ValueError, requests.RequestException):
            if attempt == 11:
                raise
            time.sleep(min(15, 2 * (attempt + 1)))


def download(name, workers):
    if name not in ALLOWED:
        raise ValueError(f'file not on allowlist: {name}')
    DEST.mkdir(parents=True, exist_ok=True)
    url = BASE + name
    response = requests.head(url, allow_redirects=True, timeout=30)
    response.raise_for_status()
    size = int(response.headers['Content-Length'])
    output = DEST / name
    if not output.exists() or output.stat().st_size != size:
        part_dir = DEST / (name + '.parts')
        part_dir.mkdir(exist_ok=True)
        ranges = [(i, i * PART_BYTES, min(size - 1, (i + 1) * PART_BYTES - 1))
                  for i in range((size + PART_BYTES - 1) // PART_BYTES)]
        print(f'{name}: {size} bytes, {len(ranges)} ranges', flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(get_part, url, size, part_dir, *item) for item in ranges]
            for job in concurrent.futures.as_completed(jobs):
                job.result()
        assembled = DEST / (name + '.assembled')
        with assembled.open('wb') as target:
            for i, _, _ in ranges:
                with (part_dir / f'{i:04d}.part').open('rb') as source:
                    shutil.copyfileobj(source, target, 8 * 1024 * 1024)
        if assembled.stat().st_size != size:
            raise ValueError('assembled size mismatch')
        os.replace(assembled, output)
        for i, _, _ in ranges:
            (part_dir / f'{i:04d}.part').unlink()
        part_dir.rmdir()
    digest = hashlib.sha256()
    with output.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            digest.update(block)
    record = {'url': url, 'file': str(output), 'bytes': size, 'sha256': digest.hexdigest()}
    (DEST / (name + '.source.json')).write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', choices=sorted(ALLOWED))
    parser.add_argument('--workers', type=int, default=6)
    args = parser.parse_args()
    for filename in args.files:
        download(filename, args.workers)
