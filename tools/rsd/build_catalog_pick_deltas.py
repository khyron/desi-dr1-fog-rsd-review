#!/usr/bin/env python3
"""Materialize FoG/RSD deltas in original .desi row order for browser picking.

The public point chunks are Morton/hash ordered, while desiGalaxyLod picks from the
four original .desi catalogues.  This creates one gzipped Float32 pair per catalogue:
`[fog_delta_mpc, rsd_delta_mpc]` for each row.  It contains no TARGETIDs.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TRACERS = ("BGS_BRIGHT", "LRG", "ELG", "QSO")


def target_ids(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    if raw[:4] != b"DESI":
        raise ValueError(f"{path}: not a DESI catalogue")
    version = int.from_bytes(raw[4:8], "little")
    count = int.from_bytes(raw[8:12], "little")
    if version < 3:
        raise ValueError(f"{path}: v{version} has no TARGETID table")
    header = 48
    table = header + count * 8  # positions (6 B) plus quantized redshift (2 B)
    return np.frombuffer(raw, dtype="<u8", count=count, offset=table + count * 16).copy()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalogs", default=str(ROOT / "assets"))
    ap.add_argument("--target-index", default=str(ROOT / "data_external" / "fog" / "chunk-targetids"))
    ap.add_argument("--fog-chunks", default=str(ROOT / "ply" / "chunks_fog_v1"))
    ap.add_argument("--rsd-chunks", default=str(ROOT / "ply" / "chunks_rsd_v1"))
    ap.add_argument("--out", default=str(ROOT / "ply" / "rsd_catalog_pick_v1"))
    args = ap.parse_args()

    fog_dir, rsd_dir, index_dir, out = map(Path, (args.fog_chunks, args.rsd_chunks, args.target_index, args.out))
    out.mkdir(parents=True, exist_ok=True)
    fog_manifest = json.loads((fog_dir / "manifest.json").read_text())
    ids_all, fog_all, rsd_all = [], [], []
    for n, chunk in enumerate(fog_manifest["chunks"]):
        index_name = chunk.get("fogTargetIndex")
        if index_name:
            index_path = index_dir / index_name
        else:
            matches = sorted(index_dir.glob(f"targetid_{n:03d}_*.npy"))
            if len(matches) != 1:
                raise FileNotFoundError(f"chunk {n}: private TARGETID index not found")
            index_path = matches[0]
        ids = np.load(index_path, mmap_mode="r")
        fog = np.frombuffer(gzip.decompress((fog_dir / chunk["fog"]).read_bytes()), dtype="<f4")
        rsd = np.frombuffer(gzip.decompress((rsd_dir / f"rsd_{n:03d}.bin").read_bytes()), dtype="<f4")
        if not (len(ids) == len(fog) == len(rsd) == chunk["count"]):
            raise ValueError(f"chunk {n}: count mismatch")
        ids_all.append(np.asarray(ids)); fog_all.append(fog); rsd_all.append(rsd)
    ids = np.concatenate(ids_all); fog = np.concatenate(fog_all); rsd = np.concatenate(rsd_all)
    order = np.argsort(ids, kind="stable")
    ids, fog, rsd = ids[order], fog[order], rsd[order]
    # Viewer geometry can contain repeated TARGETIDs.  Corrections are TARGETID-based, so one
    # stable representative serves every repeated catalogue row.
    unique = np.r_[True, ids[1:] != ids[:-1]]
    duplicates = int(len(ids) - np.count_nonzero(unique))
    if duplicates:
        print(f"deduplicated {duplicates:,} repeated TARGETIDs")
        ids, fog, rsd = ids[unique], fog[unique], rsd[unique]

    manifest = {"version": 1, "format": "gzip-float32le-fog-rsd-pairs-v1", "catalogs": {}}
    for tracer in TRACERS:
        cat_ids = target_ids(Path(args.catalogs) / f"{tracer}.desi")
        slot = np.searchsorted(ids, cat_ids)
        hit = (slot < len(ids)) & (ids[np.minimum(slot, len(ids) - 1)] == cat_ids)
        pair = np.zeros((len(cat_ids), 2), dtype="<f4")
        pair[hit, 0] = fog[slot[hit]]
        pair[hit, 1] = rsd[slot[hit]]
        name = f"correction_{tracer}.bin"
        blob = gzip.compress(pair.tobytes(), compresslevel=9)
        (out / name).write_bytes(blob)
        manifest["catalogs"][tracer] = {"file": name, "count": len(cat_ids), "bytes": len(blob),
            "fogNonzero": int(np.count_nonzero(pair[:, 0])), "rsdNonzero": int(np.count_nonzero(pair[:, 1]))}
        print(f"{tracer}: {len(cat_ids):,} rows, {len(blob) / 1e6:.2f} MB")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
