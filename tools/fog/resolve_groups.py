#!/usr/bin/env python3
"""Resolve Gfinder memberships/groups for viewer-matched galaxies, streaming FITS only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
FOG = ROOT / "data_external" / "fog"
GF = ROOT / "data_external" / "desi" / "dr1" / "gfinder" / "v1.0"
REL_DTYPE = np.dtype([("igal", "<i8"), ("igrp", "<i8"), ("rank", "<i2")])
GROUP_DTYPE = np.dtype([("igrp", "<i8"), ("rich", "<i4"), ("ra", "<f8"), ("dec", "<f8"), ("z", "<f4"), ("logm", "<f4"), ("logl", "<f4")])


def contains(sorted_ids: np.ndarray, values: np.ndarray) -> np.ndarray:
    slot = np.searchsorted(sorted_ids, values)
    return (slot < len(sorted_ids)) & (sorted_ids[np.minimum(slot, len(sorted_ids) - 1)] == values)


def matched_igals(match_dir: Path) -> np.ndarray:
    paths = sorted(match_dir.glob("match_*.npy"))
    if not paths:
        raise RuntimeError("no match batches")
    return np.unique(np.concatenate([np.load(p, mmap_mode="r")["igal"] for p in paths]))


def run(args: argparse.Namespace) -> int:
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    matched = matched_igals(Path(args.matches))
    matched_path = out / "matched-igals.npy"; np.save(matched_path, matched)
    relevant_parts, relevant_count = [], 0
    with fits.open(args.relation, memmap=True, lazy_load_hdus=True) as hdul:
        data = hdul[1].data
        for start in range(0, len(data), args.batch):
            stop = min(len(data), start + args.batch); block = data[start:stop]
            keep = contains(matched, np.asarray(block["IGAL"]))
            if not keep.any(): continue
            rows = np.empty(int(keep.sum()), dtype=REL_DTYPE)
            rows["igal"] = np.asarray(block["IGAL"])[keep]
            rows["igrp"] = np.asarray(block["IGRP"])[keep]
            rows["rank"] = np.asarray(block["RANK"])[keep]
            path = out / f"relevant_{start:09d}_{stop:09d}.npy"; np.save(path, rows); relevant_parts.append(path)
            relevant_count += len(rows)
            print(f"relation {stop:,}/{len(data):,}: relevant {relevant_count:,}", flush=True)
    relevant = np.concatenate([np.load(p, mmap_mode="r") for p in relevant_parts])
    groups = np.unique(relevant["igrp"]); np.save(out / "relevant-igrps.npy", groups)
    group_parts = []
    with fits.open(args.groups, memmap=True, lazy_load_hdus=True) as hdul:
        data = hdul[1].data
        for start in range(0, len(data), args.batch):
            stop = min(len(data), start + args.batch); block = data[start:stop]
            keep = contains(groups, np.asarray(block["IGRP"]))
            if not keep.any(): continue
            rows = np.empty(int(keep.sum()), dtype=GROUP_DTYPE)
            rows["igrp"] = np.asarray(block["IGRP"])[keep]; rows["rich"] = np.asarray(block["RICH"])[keep]
            rows["ra"] = np.asarray(block["GRP_RA"])[keep]; rows["dec"] = np.asarray(block["GRP_DEC"])[keep]
            rows["z"] = np.asarray(block["GRP_Z"])[keep]; rows["logm"] = np.asarray(block["GRP_LOGM"])[keep]; rows["logl"] = np.asarray(block["GRP_LOGL"])[keep]
            path = out / f"group_{start:09d}_{stop:09d}.npy"; np.save(path, rows); group_parts.append(path)
    resolved = np.concatenate([np.load(p, mmap_mode="r") for p in group_parts])
    missing = int(len(groups) - len(resolved))
    (out / "summary.json").write_text(json.dumps({"matchedIgals": int(len(matched)), "relevantGroups": int(len(groups)), "resolvedGroups": int(len(resolved)), "missingGroups": missing}, indent=2), encoding="utf-8")
    print(json.dumps(json.loads((out / "summary.json").read_text()), indent=2))
    return 1 if missing else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matches", default=str(FOG / "matches")); parser.add_argument("--out", default=str(FOG / "groups"))
    parser.add_argument("--relation", default=str(GF / "iDESIDR9.y1.v1_1.fits")); parser.add_argument("--groups", default=str(GF / "DESIDR9.y1.v1_group.fits"))
    parser.add_argument("--batch", type=int, default=1_000_000)
    return run(parser.parse_args())


if __name__ == "__main__": raise SystemExit(main())
