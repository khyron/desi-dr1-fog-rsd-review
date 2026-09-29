#!/usr/bin/env python3
"""Verify and inventory downloaded DESI DR1 random realization 1 FITS files."""
import hashlib
import json
from pathlib import Path

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'tools/rsd/validation_sources'
BASE = 'https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    source = json.loads(Path(r'D:\desi\data_external\RSD\source-manifest.json').read_text())
    records = []
    for item in source['files']:
        if item['kind'] != 'random':
            continue
        name = item['file'].replace('_0_clustering.ran.fits', '_1_clustering.ran.fits')
        path = DEST / name
        if not path.is_file():
            raise FileNotFoundError(path)
        with fits.open(path, memmap=True) as hdus:
            hdus.verify('exception')
            columns = set(hdus[1].columns.names)
            if not {'TARGETID', 'RA', 'DEC', 'Z', 'WEIGHT', 'WEIGHT_FKP'} <= columns:
                raise ValueError(f'{name}: required columns absent')
            rows = int(hdus[1].header['NAXIS2'])
            release = hdus[1].header.get('DESIDR') or hdus[0].header.get('DESIDR')
        record = {'tracer': item['tracer'], 'region': item['region'], 'randomIndex': 1,
                  'url': BASE + name, 'file': name, 'bytes': path.stat().st_size,
                  'sha256': sha256(path), 'rows': rows, 'DESIDR': release}
        records.append(record)
        print(json.dumps(record), flush=True)
    if len(records) != 10:
        raise ValueError(f'expected 10 random sources, got {len(records)}')
    report = {'catalog': 'DESI DR1 LSS iron LSScats v1.5', 'purpose': 'independent random realization for later support comparison',
              'status': 'downloaded-and-file-verified', 'files': records}
    (DEST / 'source-manifest-index1.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
