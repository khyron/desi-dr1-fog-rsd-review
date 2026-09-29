#!/usr/bin/env python3
"""Download official DESI random index 1 files with checked HTTP byte ranges."""
import concurrent.futures
import json
import os
import shutil
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'tools/rsd/validation_sources'
BASE = 'https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/'
PART_BYTES = 64 * 1024 * 1024
WORKERS = 6


def download_part(url, directory, number, start, end):
    path = directory / f'{number:04d}.part'
    expected = end - start + 1
    if path.is_file() and path.stat().st_size == expected:
        return number
    for attempt in range(10):
        try:
            with requests.get(url, headers={'Range': f'bytes={start}-{end}'},
                              stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                if response.status_code != 206 or response.headers.get('Content-Range') != f'bytes {start}-{end}/{TOTAL_SIZE[url]}':
                    raise ValueError(f'invalid range response for part {number}: {response.status_code} {response.headers.get("Content-Range")}')
                with path.open('wb') as f:
                    for block in response.iter_content(1024 * 1024):
                        if block:
                            f.write(block)
            if path.stat().st_size != expected:
                raise ValueError(f'part {number} size mismatch')
            return number
        except (requests.RequestException, OSError, ValueError):
            if attempt == 9:
                raise
            time.sleep(min(10, 2 * (attempt + 1)))


TOTAL_SIZE = {}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    source = json.loads(Path(r'D:\desi\data_external\RSD\source-manifest.json').read_text())
    for entry in source['files']:
        if entry['kind'] != 'random':
            continue
        name = entry['file'].replace('_0_clustering.ran.fits', '_1_clustering.ran.fits')
        url = BASE + name
        response = requests.head(url, allow_redirects=True, timeout=30)
        response.raise_for_status()
        size = int(response.headers['Content-Length'])
        TOTAL_SIZE[url] = size
        out = DEST / name
        if out.is_file() and out.stat().st_size == size:
            print(f'verified size {name} {size}', flush=True)
            continue
        directory = DEST / (name + '.parts')
        directory.mkdir(exist_ok=True)
        ranges = [(i, i * PART_BYTES, min(size - 1, (i + 1) * PART_BYTES - 1))
                  for i in range((size + PART_BYTES - 1) // PART_BYTES)]
        print(f'downloading {name}: {size} bytes in {len(ranges)} parts', flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futures = [pool.submit(download_part, url, directory, i, start, end)
                       for i, start, end in ranges]
            for future in concurrent.futures.as_completed(futures):
                future.result()
        assembled = DEST / (name + '.assembled')
        with assembled.open('wb') as target:
            for i, _, _ in ranges:
                with (directory / f'{i:04d}.part').open('rb') as part:
                    shutil.copyfileobj(part, target, 8 * 1024 * 1024)
        if assembled.stat().st_size != size:
            raise ValueError(f'{name}: assembled size mismatch')
        os.replace(assembled, out)
        for i, _, _ in ranges:
            (directory / f'{i:04d}.part').unlink()
        directory.rmdir()
        print(f'verified size {name} {size}', flush=True)


if __name__ == '__main__':
    main()
