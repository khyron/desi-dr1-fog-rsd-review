#!/usr/bin/env python3
"""Build offline TARGETID -> source-row index for current DESI .desi catalogues.

Output is never shipped to browser. Source row is concatenated BGS_BRIGHT, LRG, ELG,
QSO order, matching build_desi_chunks.py before its Morton/partition reorder.
"""
from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TRACERS = ("BGS_BRIGHT", "LRG", "ELG", "QSO")


def target_ids(path: Path) -> np.ndarray:
    with path.open("rb") as source:
        header = source.read(48)
    if header[:4] != b"DESI":
        raise ValueError(f"{path}: invalid magic")
    version, count = struct.unpack_from("<II", header, 4)
    if version < 3:
        raise ValueError(f"{path}: FoG requires v3 TARGETID data, got v{version}")
    base = 48 + count * 8
    return np.memmap(path, mode="r", dtype="<u8", offset=base + count * 16, shape=(count,))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets", default=str(ROOT / "assets"))
    parser.add_argument("--out", default=str(ROOT / "data_external" / "fog" / "viewer-target-index.npz"))
    args = parser.parse_args()
    assets, out = Path(args.assets), Path(args.out)
    parts, offsets, offset = [], {}, 0
    for tracer in TRACERS:
        ids = target_ids(assets / f"{tracer}.desi")
        parts.append(np.asarray(ids, dtype=np.uint64))
        offsets[tracer] = {"offset": offset, "count": int(ids.size)}
        offset += ids.size
    all_ids = np.concatenate(parts)
    order = np.argsort(all_ids, kind="stable")
    sorted_ids = all_ids[order]
    duplicate_count = int(np.count_nonzero(sorted_ids[1:] == sorted_ids[:-1]))
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, targetid=sorted_ids, source_index=order.astype(np.uint32),
             metadata=json.dumps({"tracers": list(TRACERS), "rows": int(all_ids.size),
                                  "duplicates": duplicate_count, "offsets": offsets}))
    print(json.dumps({"output": str(out), "rows": int(all_ids.size),
                      "duplicates": duplicate_count, "offsets": offsets}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
