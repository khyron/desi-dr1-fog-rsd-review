#!/usr/bin/env python3
"""Stream official Gfinder GALAXY rows and retain only TARGETID matches in viewer index.

No browser input. No FITS table copy. Output batches are resumable intermediate data.
"""
from __future__ import annotations

import argparse
import json
import sys
import types
from pathlib import Path

import numpy as np
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[2]
GFINDER = ROOT / "data_external" / "desi" / "dr1" / "gfinder" / "v1.0"
INDEX = ROOT / "data_external" / "fog" / "viewer-target-index.npz"
OUT = ROOT / "data_external" / "fog" / "matches"
DTYPE = np.dtype([
    ("targetid", "<u8"), ("igal", "<i8"), ("ra", "<f8"), ("dec", "<f8"),
    ("z", "<f4"), ("zsrc", "<i2"), ("source_hits", "<u2")
])


def official_encoder():
    # targets.py imports healpy for unrelated routines. encode_targetid itself does not use it.
    sys.modules.setdefault("healpy", types.SimpleNamespace())
    from desitarget.targets import encode_targetid
    return encode_targetid


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", default=str(INDEX))
    parser.add_argument("--gfinder", default=str(GFINDER / "DESIDR9.y1.v1_galaxy.fits"))
    parser.add_argument("--out", default=str(OUT))
    parser.add_argument("--batch", type=int, default=1_000_000)
    parser.add_argument("--limit", type=int, default=0, help="rows; 0 = full catalog")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    source = np.load(args.index, allow_pickle=False)
    viewer_ids = source["targetid"]
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    encode_targetid = official_encoder()
    matched_rows = matched_viewer_rows = scanned = 0
    with fits.open(args.gfinder, memmap=True, lazy_load_hdus=True) as hdul:
        data = hdul[1].data
        total = min(len(data), args.limit) if args.limit else len(data)
        for start in range(0, total, args.batch):
            stop = min(total, start + args.batch)
            path = out / f"match_{start:09d}_{stop:09d}.npy"
            if path.exists() and not args.force:
                rows = np.load(path, mmap_mode="r"); matched_rows += len(rows); scanned += stop - start
                matched_viewer_rows += int(rows["source_hits"].sum()); continue
            block = data[start:stop]
            targetid = np.asarray(encode_targetid(objid=np.asarray(block["OBJID"]),
                brickid=np.asarray(block["BRICKID"]), release=np.asarray(block["RELEASE"])), dtype=np.uint64)
            left = np.searchsorted(viewer_ids, targetid, side="left")
            right = np.searchsorted(viewer_ids, targetid, side="right")
            hit = right > left
            count = int(hit.sum())
            rows = np.empty(count, dtype=DTYPE)
            if count:
                indices = np.flatnonzero(hit)
                rows["targetid"] = targetid[indices]
                rows["igal"] = np.asarray(block["IGAL"])[indices]
                rows["ra"] = np.asarray(block["RA"])[indices]
                rows["dec"] = np.asarray(block["DEC"])[indices]
                rows["z"] = np.asarray(block["Z"])[indices]
                rows["zsrc"] = np.asarray(block["ZSRC"])[indices]
                rows["source_hits"] = (right[indices] - left[indices]).astype(np.uint16)
            np.save(path, rows)
            scanned += stop - start; matched_rows += count; matched_viewer_rows += int(rows["source_hits"].sum())
            print(f"{stop:,}/{total:,}: Gfinder matches {matched_rows:,}; viewer rows {matched_viewer_rows:,}", flush=True)
    summary = {"scannedGfinderRows": scanned, "matchedGfinderRows": matched_rows,
               "matchedViewerRows": matched_viewer_rows, "viewerRows": int(viewer_ids.size),
               "matchFraction": matched_viewer_rows / max(1, int(viewer_ids.size))}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
