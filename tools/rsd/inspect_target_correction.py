#!/usr/bin/env python3
"""Print catalogue row and stored FoG/RSD deltas for one TARGETID."""
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for_name = ("BGS_BRIGHT", "LRG", "ELG", "QSO")

def main() -> int:
    ap = argparse.ArgumentParser(); ap.add_argument("targetid", type=int); args = ap.parse_args()
    target = np.uint64(args.targetid)
    for name in for_name:
        raw = (ROOT / "assets" / f"{name}.desi").read_bytes(); count = int.from_bytes(raw[8:12], "little")
        base = 48 + count * 8
        ids = np.frombuffer(raw, dtype="<u8", count=count, offset=base + count * 16)
        hits = np.flatnonzero(ids == target)
        if not len(hits): continue
        pair = np.frombuffer(gzip.decompress((ROOT / "ply" / "rsd_catalog_pick_v1" / f"correction_{name}.bin").read_bytes()), dtype="<f4").reshape(-1, 2)
        for row in hits:
            print({"tracer": name, "row": int(row), "fogMpc": float(pair[row, 0]), "rsdMpc": float(pair[row, 1]), "combinedMpc": float(pair[row].sum())})
        return 0
    print("TARGETID not found"); return 1

if __name__ == "__main__": raise SystemExit(main())
