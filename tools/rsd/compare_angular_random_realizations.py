#!/usr/bin/env python3
"""Compare angular support flags from independent random 0 and 1 realizations."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2] / 'tools/rsd/extensions'


def main():
    a = json.loads((ROOT / 'angular-mask-audit.json').read_text())
    b = json.loads((ROOT / 'angular-mask-audit-index1.json').read_text())
    if a['randomIndex'] != 0 or b['randomIndex'] != 1 or a['samplePerField'] != b['samplePerField']:
        raise ValueError('incompatible random realization audits')
    rows = []
    with np.load(ROOT / 'angular-mask-audit-flagged-targetids.npz') as flags0, \
         np.load(ROOT / 'angular-mask-audit-index1-flagged-targetids.npz') as flags1:
        for field0, field1 in zip(a['fields'], b['fields']):
            if (field0['tracer'], field0['region']) != (field1['tracer'], field1['region']):
                raise ValueError('field order differs')
            key = f"{field0['tracer']}_{field0['region']}"
            ids0 = np.sort(flags0[key])
            ids1 = np.sort(flags1[key])
            both = np.intersect1d(ids0, ids1, assume_unique=True)
            rows.append({'tracer': field0['tracer'], 'region': field0['region'],
                         'sampleRows': field0['sampledExtensionRows'],
                         'flaggedIndex0': len(ids0), 'flaggedIndex1': len(ids1),
                         'flaggedBoth': len(both),
                         'flaggedEither': len(ids0) + len(ids1) - len(both),
                         'flagAgreementFraction': (field0['sampledExtensionRows'] - len(ids0) - len(ids1) + 2*len(both)) / field0['sampledExtensionRows']})
    report = {'status': 'random-realization-proxy-comparison', 'fields': rows,
              'sampleRows': sum(x['sampleRows'] for x in rows),
              'flaggedBoth': sum(x['flaggedBoth'] for x in rows),
              'flaggedEither': sum(x['flaggedEither'] for x in rows),
              'limitations': 'Agreement between random realizations is not an independent official veto-mask or radial-selection validation.'}
    (ROOT / 'angular-random-realization-comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'sampleRows': report['sampleRows'],
                      'flaggedBoth': report['flaggedBoth'], 'flaggedEither': report['flaggedEither']}))


if __name__ == '__main__':
    main()
