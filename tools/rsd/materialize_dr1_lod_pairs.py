#!/usr/bin/env python3
"""Build exact TARGETID-aligned FoG/RSD pairs for DR1 GalaxyLOD catalog assets.

The outputs match the row order in DR1_GALAXY.desi and DR1_QSO.desi, verified
against the exact per-catalog TARGETID indexes and their embedded v3 TARGETIDs.
Unmatched source IDs retain zero corrections. No coordinate/proximity matching.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
NAMES = ("DR1_GALAXY", "DR1_QSO")


def validate_catalog(path: Path) -> tuple[int, np.ndarray | None]:
    with path.open("rb") as f:
        header = f.read(48)
    if len(header) != 48 or header[:4] != b"DESI":
        raise ValueError(f"{path}: invalid DESI header")
    version = int.from_bytes(header[4:8], "little")
    count = int.from_bytes(header[8:12], "little")
    if version not in (2, 3):
        raise ValueError(f"{path}: unsupported .desi v{version}")
    expected = 48 + count * (32 if version >= 3 else 12)
    if path.stat().st_size != expected:
        raise ValueError(f"{path}: expected {expected} bytes, got {path.stat().st_size}")
    if version >= 3:
        return count, np.memmap(path, dtype="<u8", mode="r", offset=48 + count * 16, shape=(count,))
    return count, None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalogs", type=Path, default=ROOT / "assets")
    ap.add_argument("--index", type=Path, default=ROOT / "ply" / "reconstruction_dr1" / "private_index")
    ap.add_argument("--payloads", type=Path, default=ROOT / "ply" / "reconstruction_dr1" / "payloads")
    ap.add_argument("--out", type=Path, default=ROOT / "ply" / "reconstruction_dr1" / "lod_pairs")
    args = ap.parse_args()

    index_report = json.loads((args.index / "index-report.json").read_text(encoding="utf-8"))
    if int(index_report.get("count", 0)) <= 0 or not all(c.get("geometryExact") for c in index_report["chunks"]):
        raise ValueError("private chunk TARGETID index is not validated against geometry")
    source_ids, source_fog, source_rsd = [], [], []
    for spec in index_report["chunks"]:
        i = int(spec["chunk"])
        ids = np.load(args.index / spec["index"], mmap_mode="r")
        n = int(spec["rows"])
        if ids.dtype != np.dtype("uint64") or len(ids) != n:
            raise ValueError(f"chunk {i}: invalid TARGETID index dtype/count")
        layers = []
        for prefix in ("fog", "rsd"):
            path = args.payloads / f"{prefix}_{i:03d}.bin"
            with gzip.open(path, "rb") as f:
                raw = f.read()
            if len(raw) != n * 4:
                raise ValueError(f"{path}: expected {n * 4} bytes, got {len(raw)}")
            values = np.frombuffer(raw, dtype="<f4")
            if not np.isfinite(values).all():
                raise ValueError(f"{path}: non-finite delta")
            layers.append(values)
        source_ids.append(np.asarray(ids, dtype="<u8"))
        source_fog.append(layers[0])
        source_rsd.append(layers[1])

    ids = np.concatenate(source_ids)
    fog = np.concatenate(source_fog)
    rsd = np.concatenate(source_rsd)
    if len(ids) != int(index_report["count"]):
        raise ValueError("chunk index aggregate count mismatch")
    order = np.argsort(ids, kind="stable")
    ids, fog, rsd = ids[order], fog[order], rsd[order]
    if np.any(ids[1:] == ids[:-1]):
        raise ValueError("duplicate TARGETID in exact point index; ambiguous correction join")

    args.out.mkdir(parents=True, exist_ok=True)
    report = {"version": 1, "format": "gzip-float32le-fog-rsd-pairs-v1",
              "join": "exact TARGETID; unmatched rows zero", "catalogs": {}}
    for name in NAMES:
        target_path = args.index / f"{name}_targetid.npy"
        targetids = np.load(target_path, mmap_mode="r")
        cat_path = args.catalogs / f"{name}.desi"
        catalog_count, embedded = validate_catalog(cat_path)
        if len(targetids) != catalog_count:
            raise ValueError(f"{name}: private TARGETID index count differs from v2 catalog")
        if embedded is not None and not np.array_equal(targetids, embedded):
            raise ValueError(f"{name}: private TARGETID index differs from v3 catalog row order")
        if np.any(targetids[1:] == targetids[:-1]):
            raise ValueError(f"{name}: duplicate TARGETID in catalog index")
        slot = np.searchsorted(ids, targetids)
        found = (slot < len(ids)) & (ids[np.minimum(slot, len(ids) - 1)] == targetids)
        pairs = np.zeros((len(targetids), 2), dtype="<f4")
        pairs[found, 0] = fog[slot[found]]
        pairs[found, 1] = rsd[slot[found]]
        if not np.isfinite(pairs).all():
            raise ValueError(f"{name}: generated pair contains non-finite values")
        blob = gzip.compress(pairs.tobytes(order="C"), compresslevel=9)
        out_path = args.out / f"correction_{name}.bin"
        out_path.write_bytes(blob)

        # Validate the actual serialized asset, including row count, before reporting it.
        decoded = gzip.decompress(out_path.read_bytes())
        if len(decoded) != len(targetids) * 8:
            raise ValueError(f"{out_path}: serialized count/size mismatch")
        roundtrip = np.frombuffer(decoded, dtype="<f4").reshape((-1, 2))
        if not np.array_equal(roundtrip, pairs):
            raise ValueError(f"{out_path}: serialized values differ from exact join")
        report["catalogs"][name] = {
            "file": out_path.name, "rows": len(targetids), "bytes": len(blob),
            "targetIdsFoundInPointIndex": int(found.sum()),
            "fogNonzero": int(np.count_nonzero(pairs[:, 0])),
            "rsdNonzero": int(np.count_nonzero(pairs[:, 1])),
            "validated": True,
        }
        print(json.dumps(report["catalogs"][name]), flush=True)

    (args.out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
