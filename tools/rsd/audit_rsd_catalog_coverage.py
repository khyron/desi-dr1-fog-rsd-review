#!/usr/bin/env python3
"""Audit exact RSD-cache coverage against the current 15.8M DR1 point catalog.

The point catalog omits TARGETID in its published v2 bytes; this audit uses the
private exact-validated TARGETID index and redshifts recovered from those same
catalog bytes. It never fills values by nearest-neighbour interpolation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from build_desi_chunks import morton_order

FIELDS = [('BGS', 0.1, 0.4), ('LRG', 0.4, 0.8),
          ('LRG_ELG', 0.8, 1.1), ('ELG', 1.1, 1.6), ('QSO', 1.6, 2.1)]


def read_rsd(cache_dir: Path):
    output = {}
    for field, _, _ in FIELDS:
        arrays = []
        for region in ('NGC', 'SGC'):
            path = cache_dir / f'{field}_{region}_rsd_query_v1.npz'
            with np.load(path) as cache:
                arrays.append((np.asarray(cache['targetid'], dtype='<u8'),
                               np.asarray(cache['delta_mpc'], dtype='<f4')))
        ids = np.concatenate([a[0] for a in arrays]); delta = np.concatenate([a[1] for a in arrays])
        order = np.argsort(ids, kind='stable'); ids, delta = ids[order], delta[order]
        if np.any(ids[1:] == ids[:-1]) or not np.isfinite(delta).all():
            raise ValueError(f'{field}: duplicated IDs or nonfinite values in cache')
        output[field] = (ids, delta)
    return output


def current_redshift_chunks(assets: Path, manifest: dict, index: Path):
    units = float(manifest['unitsPerMpc'])
    positions, redshifts = [], []
    for name in manifest['tracers']:
        path = assets / f'{name}.desi'
        with path.open('rb') as stream:
            header = stream.read(24)
        version, count, own_units, zmin, zmax = __import__('struct').unpack('<IIfff', header[4:24])
        if header[:4] != b'DESI' or version != 2:
            raise ValueError(f'{path}: expected current DR1 v2 geometry')
        pos = np.memmap(path, dtype='<i2', mode='r', offset=48, shape=(count, 3))
        zq = np.memmap(path, dtype='<u2', mode='r', offset=48 + count * 6, shape=(count,))
        positions.append(np.rint(pos.astype('f4') / own_units * units).astype('<i2'))
        redshifts.append(zmin + zq.astype('f4') * ((zmax - zmin) / 65535.0))
    pos = np.concatenate(positions); z = np.concatenate(redshifts)
    hashes = np.abs(pos[:, 0].astype('i8') * 17 + pos[:, 2].astype('i8') * 59 + pos[:, 1].astype('i8') * 101)
    for chunk in range(len(manifest['chunks'])):
        chosen = np.flatnonzero(hashes % len(manifest['chunks']) == chunk)
        chosen = chosen[morton_order(pos[chosen])]
        ids = np.load(index / f'targetid_{chunk:03d}.npy', mmap_mode='r')
        if len(chosen) != len(ids) or len(ids) != manifest['chunks'][chunk]['count']:
            raise ValueError(f'chunk {chunk}: redshift/identity order mismatch')
        yield chunk, z[chosen]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunks", type=Path, default=Path("ply/chunks"))
    ap.add_argument("--assets", type=Path, default=Path("assets"))
    ap.add_argument("--index", type=Path, default=Path("ply/reconstruction_dr1/private_index"))
    ap.add_argument("--cache-dir", type=Path, default=Path(r"D:\desi\data_external\RSD\derived\rsd-v1-exploratory"))
    ap.add_argument("--out", type=Path, default=Path("tools/rsd/rsd-coverage-audit.json"))
    args = ap.parse_args()

    manifest = json.loads((args.chunks / "manifest.json").read_text(encoding="utf-8"))
    index_report = json.loads((args.index / "index-report.json").read_text(encoding="utf-8"))
    if index_report["count"] != manifest["count"] or not all(c["geometryExact"] for c in index_report["chunks"]):
        raise ValueError("current private TARGETID index is not validated against the published geometry")
    fields = read_rsd(args.cache_dir)
    totals = {"objects": 0, "insideConfiguredRedshiftEnvelope": 0,
              "exactCacheMatches": 0, "insideEnvelopeWithoutExactSample": 0,
              "outsideConfiguredRedshiftEnvelope": 0}
    by_field = {name: {"eligible": 0, "exactMatches": 0} for name, _, _ in FIELDS}
    matched_deltas = []
    for i, z in current_redshift_chunks(args.assets, manifest, args.index):
        ids = np.load(args.index / f"targetid_{i:03d}.npy", mmap_mode="r")
        if len(ids) != len(z):
            raise ValueError(f"chunk {i}: TARGETID/redshift count mismatch")
        eligible = np.zeros(len(ids), dtype=bool)
        matched = np.zeros(len(ids), dtype=bool)
        for field, low, high in FIELDS:
            mask = (z >= low) & (z < high)
            eligible |= mask
            by_field[field]["eligible"] += int(mask.sum())
            if mask.any():
                source_ids, deltas = fields[field]
                pos = np.searchsorted(source_ids, ids[mask])
                hit = (pos < len(source_ids)) & (source_ids[np.minimum(pos, len(source_ids) - 1)] == ids[mask])
                matched[np.flatnonzero(mask)[hit]] = True
                by_field[field]["exactMatches"] += int(hit.sum())
                if hit.any():
                    matched_deltas.append(deltas[pos[hit]])
        totals["objects"] += len(ids)
        totals["insideConfiguredRedshiftEnvelope"] += int(eligible.sum())
        totals["exactCacheMatches"] += int(matched.sum())
        totals["insideEnvelopeWithoutExactSample"] += int((eligible & ~matched).sum())
        totals["outsideConfiguredRedshiftEnvelope"] += int((~eligible).sum())

    n = totals["objects"]
    env = totals["insideConfiguredRedshiftEnvelope"]
    totals["exactCacheCoverageOfCatalog"] = totals["exactCacheMatches"] / n if n else 0
    totals["exactCacheCoverageWithinEnvelope"] = totals["exactCacheMatches"] / env if env else 0
    if matched_deltas:
        delta = np.concatenate(matched_deltas)
        delta_summary = {"medianAbsMpc": float(np.median(np.abs(delta))),
                         "p90AbsMpc": float(np.percentile(np.abs(delta), 90)),
                         "p99AbsMpc": float(np.percentile(np.abs(delta), 99)),
                         "maxAbsMpc": float(np.max(np.abs(delta)))}
    else:
        delta_summary = {}
    result = {
        "catalog": "DESI DR1 current 15,786,217-object point catalog",
        "identitySource": "exact-validated private targetid_*.npy index",
        "selection": [{"field": f, "zMin": low, "zMaxExclusive": high} for f, low, high in FIELDS],
        "maskAndRegion": "Only exact TARGETIDs in NGC/SGC query caches are counted as covered. Unmatched objects have no mask-support evidence in the saved files.",
        "spatialFieldAvailability": "Query caches contain TARGETID/delta samples, not serialized spatial fields. Spatial interpolation is not executable from current derived files.",
        "totals": totals,
        "byField": by_field,
        "matchedDeltaMagnitude": delta_summary,
        "interpretation": {
            "exactCacheMatches": "supported now by exact identity lookup",
            "insideEnvelopeWithoutExactSample": "upper bound on possible extra query targets only; actual valid footprint/mask support is unknown and a saved reconstruction mesh is absent",
            "outsideConfiguredRedshiftEnvelope": "requires a newly defined tracer/redshift reconstruction with corresponding LSS data and random catalogues",
            "unmatchedPolicy": "retain observed position (zero RSD delta); do not interpolate outside demonstrated valid mask/volume",
        },
        "validationBasis": "Redshifts and TARGETIDs share the reconstructed exact point partition; index report marks every chunk geometryExact.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
