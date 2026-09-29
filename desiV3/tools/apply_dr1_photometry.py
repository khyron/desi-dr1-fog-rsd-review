"""Apply verified official Tractor imaging columns by exact uint64 TARGETID.

Inputs and prior enriched chunks are preserved. Output is a separate candidate.
No network calls; fails unless every official source file has been processed.
"""
from pathlib import Path
import gzip, hashlib, json, struct
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'desiV3/data/missing-photometry-sql'
OUT = ROOT / 'desiV3/data/full-photometry'

def main():
    progress = json.loads((SOURCE / 'download-progress.json').read_text())
    assert progress['complete'] == progress['total'], 'Official download is incomplete'
    ids = [np.load(ROOT / 'ply/reconstruction_dr1/private_index' / (name + '_targetid.npy'))
           for name in ['DR1_GALAXY', 'DR1_QSO']]
    keys = np.concatenate(ids)
    order = np.argsort(keys, kind='stable')
    sorted_keys = keys[order]
    values = np.zeros((len(keys), 4), dtype='u1')
    found = np.zeros(len(keys), dtype=bool)
    ambiguous = np.zeros(len(keys), dtype=bool)
    raw_measurements = np.full((len(keys), 4), np.nan, dtype='f4')
    profiles = np.full(len(keys), b'', dtype='S4')
    source_pixel = np.full(len(keys), -1, dtype='i2')
    source_records = []
    for record in progress['files']:
        name = record.get('sourceFile', record['sourceURL'].split('/')[-1])
        path = SOURCE / (name + '.npz')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['compactSHA256']
        pixel = record.get('sourcePixel', int(name.split('-hp')[1][:3]) if '-hp' in name else -1)
        with np.load(path) as data:
            tid = data['TARGETID'].astype('u8')
            assert len(np.unique(tid)) == len(tid), 'Duplicate IDs within official pixel'
            slot = np.searchsorted(sorted_keys, tid)
            safe = np.minimum(slot, len(keys) - 1)
            hit = (slot < len(keys)) & (sorted_keys[safe] == tid)
            indices = order[safe[hit]]
            types = np.char.strip(data['TYPE'][hit].astype('S4'))
            code = np.full(len(indices), 255, dtype='u1')
            for label, number in [('PSF', 0), ('REX', 1), ('EXP', 2), ('DEV', 3), ('SER', 4)]:
                code[types == label.encode()] = number
            shape = data['SHAPE_R'][hit].astype('f4')
            flux = np.stack([data['FLUX_' + band][hit] for band in ['G', 'R', 'Z']], axis=1).astype('f4')
            mag = np.full(flux.shape, np.nan, dtype='f4')
            good = np.isfinite(flux) & (flux > 0)
            mag[good] = 22.5 - 2.5 * np.log10(flux[good])
            gr, rz = mag[:, 0] - mag[:, 1], mag[:, 1] - mag[:, 2]
            encoded = np.stack([code,
                np.clip(np.rint(np.nan_to_num(shape, nan=0, posinf=0, neginf=0) * 10), 0, 255).astype('u1'),
                np.clip(np.rint(np.where(np.isfinite(gr), (gr + .5) * 85, 128)), 0, 255).astype('u1'),
                np.clip(np.rint(np.where(np.isfinite(rz), (rz + .5) * 102, 128)), 0, 255).astype('u1')], axis=1)
            measures = np.column_stack([shape, flux])
            conflict = found[indices] & ((profiles[indices] != types) |
                ~np.all((raw_measurements[indices] == measures) |
                        (np.isnan(raw_measurements[indices]) & np.isnan(measures)), axis=1))
            ambiguous[indices[conflict]] = True
            values[indices] = encoded
            profiles[indices] = types
            raw_measurements[indices] = measures
            source_pixel[indices] = pixel
            found[indices] = True
        source_records.append(dict(**record, matchedCatalogueRows=int(hit.sum())))
        print(name, 'matched', int(hit.sum()), flush=True)
    valid = found & ~ambiguous & (values[:, 0] != 255)
    OUT.mkdir(exist_ok=True)
    reports = []
    for number in range(12):
        name = f'catalog_{number:03d}.bin.gz'
        source = ROOT / 'desiV3/data/morphology' / name
        raw = bytearray(gzip.decompress(source.read_bytes()))
        before = bytes(raw)
        version, n, o_type, o_row, o_pos, o_z, o_attrs, total = struct.unpack_from('<8I', raw, 4)
        assert version == 1 and len(raw) == total
        qso = (np.frombuffer(raw, 'u1', n, o_type) & 64) != 0
        row = np.frombuffer(raw, '<u4', n, o_row)
        global_rows = row.astype('i8') + qso.astype('i8') * len(ids[0])
        view = np.frombuffer(raw, 'u1', 4 * n, o_attrs).reshape(4, n)
        hit = valid[global_rows]
        previously_missing = view[0] == 255
        view[:, hit] = values[global_rows[hit]].T
        assert raw[:o_attrs] == before[:o_attrs] and raw[o_attrs + 4*n:] == before[o_attrs + 4*n:]
        blob = gzip.compress(raw, compresslevel=6, mtime=0)
        (OUT / name).write_bytes(blob)
        reports.append(dict(chunk=number, count=n, bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(),
            officialMatched=int(hit.sum()), newlyRecovered=int((hit & previously_missing).sum()),
            galaxyMissing=int(((view[0] == 255) & ~qso).sum()),
            qsoMissing=int(((view[0] == 255) & qso).sum()), geometryAndIdentityUnchanged=True))
    # Full precision measurements and provenance are retained for future API import.
    # These private row-aligned arrays are not additional browser downloads.
    np.savez_compressed(OUT / 'api-imaging.npz', TARGETID=keys[found],
        PROFILE=profiles[found], SHAPE_R=raw_measurements[found, 0],
        FLUX_G=raw_measurements[found, 1], FLUX_R=raw_measurements[found, 2],
        FLUX_Z=raw_measurements[found, 3], SOURCE_PIXEL=source_pixel[found], AMBIGUOUS=ambiguous[found])
    report = dict(method='Exact uint64 TARGETID join to public desi_dr1.photometry; only previously missing rows requested',
        sources=source_records, chunks=reports, officialFound=int(found.sum()),
        ambiguous=int(ambiguous.sum()), officialRecognizedProfile=int(valid.sum()),
        unsupportedProfiles=dict(zip(*[a.tolist() for a in np.unique(profiles[found & (values[:, 0] == 255)].astype('U4'), return_counts=True)])),
        colorQuantization='g-r: (mag+.5)*85; r-z: (mag+.5)*102; round and clamp to uint8',
        sizeQuantization='round(arcsec*10), clamp 0..255',
        invalidColorFallbackByte=128, rawValuesPreserved=True)
    (OUT / 'validation.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k not in ['sources', 'chunks']}), flush=True)

if __name__ == '__main__': main()
