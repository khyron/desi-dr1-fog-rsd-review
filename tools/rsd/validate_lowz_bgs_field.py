#!/usr/bin/env python3
"""Technical and seam audit of local BGS BRIGHT candidate against active RSD."""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools' / 'rsd'))
from extend_rsd_catalog import load_viewer_catalog  # noqa: E402

FIELD = ROOT / 'tools' / 'rsd' / 'lowz_sources' / 'field_candidate'
PAYLOAD = ROOT / 'ply' / 'reconstruction_dr1' / 'staged_extensions_v1' / 'payloads'


def main():
    report = json.loads((FIELD / 'field-manifest.json').read_text())
    if report['status'] != 'complete-for-technical-review' or len(report['regions']) != 2:
        raise ValueError('two completed fields required')
    ids_parts, values_parts = [], []
    for region in ('NGC', 'SGC'):
        with np.load(FIELD / f'BGS_BRIGHT_{region}_viewer_rsd_v1.npz') as data:
            ids_parts.append(np.asarray(data['targetid'], dtype='u8'))
            values_parts.append(np.asarray(data['delta_mpc'], dtype='f4'))
    ids = np.concatenate(ids_parts)
    values = np.concatenate(values_parts)
    order = np.argsort(ids)
    ids, values = ids[order], values[order]
    if np.any(ids[1:] == ids[:-1]) or not np.isfinite(values).all():
        raise ValueError('duplicate IDs or nonfinite candidate values')
    manifest, pos, z, viewer_ids = load_viewer_catalog(
        ROOT / 'assets', ROOT / 'ply' / 'chunks', ROOT / 'ply' / 'reconstruction_dr1' / 'private_index')
    old = np.concatenate([
        np.frombuffer(gzip.decompress((PAYLOAD / f'rsd_{i:03d}.bin').read_bytes()), dtype='<f4')
        for i in range(len(manifest['chunks']))])
    if len(old) != len(z):
        raise ValueError('active RSD payload order mismatch')
    at = np.searchsorted(ids, viewer_ids)
    hit = (at < len(ids)) & (ids[np.minimum(at, len(ids)-1)] == viewer_ids)
    new = np.zeros(len(z), dtype='f4')
    new[hit] = values[at[hit]]
    if np.any(hit & ((z < 0.01) | (z >= 0.4))):
        raise ValueError('candidate escaped BGS redshift range')
    updated = old.copy()
    bgs = (z >= 0.01) & (z < 0.4)
    updated[bgs] = new[bgs]
    bins = [(0.01,0.05),(0.05,0.08),(0.08,0.09),(0.09,0.095),
            (0.095,0.1),(0.1,0.105),(0.105,0.11),(0.11,0.12),(0.12,0.4)]
    shells = []
    for low, high in bins:
        m = (z >= low) & (z < high)
        shells.append({'z': [low,high], 'viewerRows': int(m.sum()),
                       'candidateMatched': int((hit&m).sum()),
                       'activeRsdNonzero': int(np.count_nonzero(old[m])),
                       'candidateRsdNonzero': int(np.count_nonzero(new[m])),
                       'activeMedianAbsMpc': float(np.median(np.abs(old[m]))),
                       'candidateMedianAbsMpc': float(np.median(np.abs(new[m])))})
    overlap = hit & (z >= 0.1) & (z < 0.4) & (old != 0)
    delta = new[overlap] - old[overlap]
    # Count the radial positions in 5 Mpc shells across the previous z=0.10 cut.
    radius = np.linalg.norm(pos, axis=1)
    r0 = float(Planck18.comoving_distance(0.1).value)
    radial_edges = np.arange(r0 - 30, r0 + 35, 5)
    in_window = (radius >= radial_edges[0] - 150) & (radius <= radial_edges[-1] + 150)
    histogram = {}
    for label, radial in [('observed',radius[in_window]),
                          ('active',radius[in_window]+old[in_window]),
                          ('candidate',radius[in_window]+updated[in_window])]:
        histogram[label] = np.histogram(radial,bins=radial_edges)[0].tolist()
    result = {'status':'technical-audit', 'scientificMaskStatus':'unvalidated',
              'viewerRows':len(z), 'candidateUniqueIds':len(ids),
              'candidateCoverageByObservedZ':shells,
              'oldVsCandidateOverlapRows':int(overlap.sum()),
              'oldVsCandidatePearsonR':float(np.corrcoef(old[overlap],new[overlap])[0,1]),
              'oldVsCandidateMedianAbsDifferenceMpc':float(np.median(np.abs(delta))),
              'oldVsCandidateP99AbsDifferenceMpc':float(np.percentile(np.abs(delta),99)),
              'previousCutComovingRadiusMpc':r0,
              'radialHistogramBinEdgesMpc':radial_edges.tolist(),
              'radialHistogramCounts':histogram,
              'note':'A continuous source field removes the algorithmic z=0.10 cutoff, but this count audit and random proximity are not an official mask or peculiar-velocity validation.'}
    (FIELD / 'validation-report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('candidateUniqueIds','oldVsCandidateOverlapRows','oldVsCandidatePearsonR','oldVsCandidateMedianAbsDifferenceMpc')},indent=2))


if __name__ == '__main__':
    main()
