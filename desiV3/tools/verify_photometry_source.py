"""Cross-check SQL extraction against an independently SHA256-verified FITS pixel."""
from pathlib import Path
import json, hashlib
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
source = ROOT / 'desiV3/data/official-photometry/tractorphot-nside4-hp000-iron.fits.npz'
receipt = json.loads(source.with_suffix('.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest() == receipt['compactSHA256']
official_hash = (source.parent / 'source.sha256sum').read_text().splitlines()[0].split()[0]
assert receipt['sourceSHA256'] == official_hash
with np.load(source) as f:
    order = np.argsort(f['TARGETID'].astype('u8'))
    reference = {k: f[k][order] for k in f.files}
ids = reference['TARGETID'].astype('u8')
counts = dict(matches=0, profileDifferences=0, measuredFloatDifferences=0)
for path in sorted((ROOT / 'desiV3/data/missing-photometry-sql').glob('*.npz')):
    with np.load(path) as f:
        keys = f['TARGETID']
        slot = np.searchsorted(ids, keys)
        safe = np.minimum(slot, len(ids)-1)
        hit = (slot < len(ids)) & (ids[safe] == keys)
        counts['matches'] += int(hit.sum())
        counts['profileDifferences'] += int((np.char.strip(reference['TYPE'][safe[hit]]) != np.char.strip(f['TYPE'][hit])).sum())
        for col in ['SHAPE_R','FLUX_G','FLUX_R','FLUX_Z']:
            a, b = reference[col][safe[hit]], f[col][hit]
            counts['measuredFloatDifferences'] += int((~((a == b) | (np.isnan(a) & np.isnan(b)))).sum())
assert counts['matches'] > 0
assert counts['profileDifferences'] == counts['measuredFloatDifferences'] == 0, counts
(ROOT / 'desiV3/data/missing-photometry-sql/source-crosscheck.json').write_text(json.dumps(dict(referencePixel=0, floatComparison='exact float32, including paired NaNs', **counts), indent=2))
print(counts)
