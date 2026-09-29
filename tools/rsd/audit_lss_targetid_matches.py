#!/usr/bin/env python3
"""Exact TARGETID cross-match of current viewer objects with DR1 LSS source rows.

This is an identity audit only. It does not infer mask support or create deltas.
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np

SELECTIONS = {
    "BGS": (0.1, 0.4),
    "LRG": (0.4, 0.8),
    "LRG_ELG": (0.8, 1.1),
    "ELG": (1.1, 1.6),
    "QSO": (1.6, 2.1),
}


def load_eligible_ids(path: Path, tracer: str) -> tuple[np.ndarray, int]:
    from astropy.io import fits

    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        table = next(h.data for h in hdul if getattr(h.data, "names", None))
        required = ("TARGETID", "RA", "DEC", "Z", "WEIGHT")
        missing = [name for name in required if name not in table.names]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        zmin, zmax = SELECTIONS[tracer]
        weight = np.asarray(table["WEIGHT"], dtype="f8")
        if "WEIGHT_FKP" in table.names:
            weight *= np.asarray(table["WEIGHT_FKP"], dtype="f8")
        z = np.asarray(table["Z"], dtype="f8")
        ra = np.asarray(table["RA"], dtype="f8")
        dec = np.asarray(table["DEC"], dtype="f8")
        keep = (np.isfinite(weight) & (weight > 0) & np.isfinite(z) & np.isfinite(ra) &
                np.isfinite(dec) & (z >= zmin) & (z < zmax))
        ids = np.asarray(table["TARGETID"][keep], dtype="u8")
    return np.unique(ids), int(keep.sum())


def load_cache_ids(path: Path) -> np.ndarray:
    with np.load(path) as data:
        ids = np.asarray(data["targetid"], dtype="u8")
    return np.unique(ids)


def load_viewer_catalog_ids_z(index: Path, assets: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    path = assets / f"{name}.desi"
    with path.open("rb") as stream:
        header = stream.read(48)
    if header[:4] != b"DESI":
        raise ValueError(f"{path}: invalid DESI geometry")
    version, count, _, zmin, zmax = struct.unpack_from("<IIfff", header, 4)
    if version != 2:
        raise ValueError(f"{path}: expected v2 geometry")
    zq = np.memmap(path, dtype="<u2", mode="r", offset=48 + count * 6, shape=(count,))
    ids = np.load(index / f"{name}_targetid.npy", mmap_mode="r")
    if len(ids) != count:
        raise ValueError(f"{name}: TARGETID count differs from geometry")
    z = zmin + np.asarray(zq, dtype="f4") * ((zmax - zmin) / 65535.0)
    return np.asarray(ids, dtype="u8"), z


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=Path("ply/reconstruction_dr1/private_index"))
    parser.add_argument("--assets", type=Path, default=Path("assets"))
    parser.add_argument("--data-dir", type=Path, default=Path(r"D:\desi\data_external\RSD"))
    parser.add_argument("--cache-dir", type=Path, default=Path(r"D:\desi\data_external\RSD\derived\rsd-v1-exploratory"))
    parser.add_argument("--out", type=Path, default=Path("tools/rsd/lss-targetid-match-audit.json"))
    args = parser.parse_args()

    galaxy_ids, galaxy_z = load_viewer_catalog_ids_z(args.index, args.assets, "DR1_GALAXY")
    qso_ids, qso_z = load_viewer_catalog_ids_z(args.index, args.assets, "DR1_QSO")
    viewer_ids = np.concatenate((galaxy_ids, qso_ids))
    viewer_z = np.concatenate((galaxy_z, qso_z))

    source_manifest = json.loads((args.data_dir / "source-manifest.json").read_text(encoding="utf-8"))
    sources = [row for row in source_manifest["files"] if row["kind"] == "data"]
    field_caches: dict[tuple[str, str], np.ndarray] = {}
    for row in sources:
        field = row["tracer"]
        if field == "LRG_ELG":
            cache_name = "LRG_ELG"
        else:
            cache_name = field
        path = args.cache_dir / f"{cache_name}_{row['region']}_rsd_query_v1.npz"
        field_caches[(field, row["region"])] = load_cache_ids(path) if path.exists() else np.empty(0, dtype="u8")

    report_rows = []
    source_viewer_union: list[np.ndarray] = []
    source_cached_union: list[np.ndarray] = []
    for row in sorted(sources, key=lambda item: (item["tracer"], item["region"])):
        path = args.data_dir / row["file"]
        ids, selected_source_rows = load_eligible_ids(path, row["tracer"])
        zmin, zmax = SELECTIONS[row["tracer"]]
        viewer_in_field = np.unique(viewer_ids[(viewer_z >= zmin) & (viewer_z < zmax)])
        on_viewer = np.intersect1d(ids, viewer_in_field, assume_unique=True)
        cached = np.intersect1d(on_viewer, field_caches[(row["tracer"], row["region"])], assume_unique=True)
        report_rows.append({
            "tracer": row["tracer"], "region": row["region"], "file": row["file"],
            "selectedRowsAfterRSDRedshiftAndWeightCuts": selected_source_rows,
            "uniqueEligibleSourceTARGETIDs": int(ids.size),
            "eligibleViewerMatches": int(on_viewer.size),
            "eligibleViewerMatchesAlreadyInExactRSDCache": int(cached.size),
            "eligibleViewerMatchesMissingExactRSDCache": int(on_viewer.size - cached.size),
        })
        source_viewer_union.append(on_viewer)
        source_cached_union.append(cached)

    matched_union = np.unique(np.concatenate(source_viewer_union)) if source_viewer_union else np.empty(0, dtype="u8")
    cached_union = np.unique(np.concatenate(source_cached_union)) if source_cached_union else np.empty(0, dtype="u8")
    envelope_mask = (viewer_z >= 0.1) & (viewer_z < 2.1)
    envelope_ids = np.unique(viewer_ids[envelope_mask])
    outside_ids = np.unique(viewer_ids[~envelope_mask])
    result = {
        "catalog": "DESI DR1 LSS iron LSScats v1.5 exact TARGETID membership audit with RSD redshift/weight cuts",
        "identityOnly": True,
        "viewerObjects": {"galaxy": int(galaxy_ids.size), "qso": int(qso_ids.size), "total": int(viewer_ids.size)},
        "lssDataFiles": report_rows,
        "totals": {
            "viewerObjectsMatchingAnyLSSDataCatalog": int(matched_union.size),
            "thoseAlreadyInExactRSDCache": int(cached_union.size),
            "matchingLSSObjectsMissingExactRSDCache": int(matched_union.size - cached_union.size),
            "viewerObjectsWithinConfiguredRedshiftEnvelope": int(envelope_ids.size),
            "withinEnvelopeNotMatchingAnyEligibleLSSDataCatalog": int(envelope_ids.size - matched_union.size),
            "outsideConfiguredRedshiftEnvelope": int(outside_ids.size),
        },
        "interpretation": (
            "Exact membership identifies targets from an input LSS sample. It does not prove spatial-mask support for unmatched targets. "
            "Objects outside the LSS catalog remain candidates for field queries only after their tracer selection, mask and radial support are validated."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
