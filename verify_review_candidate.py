#!/usr/bin/env python3
"""Check that the bundled V3 catalogue and review payloads are intact and aligned."""
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def main() -> None:
    candidate = ROOT / "review_candidate"
    records = json.loads((candidate / "manifest.json").read_text())["files"]
    for name, expected in records.items():
        path = candidate / name
        assert path.stat().st_size == expected["bytes"], name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected["sha256"], name
    catalogs = ROOT / "desiV3" / "catalogs"
    for line in (catalogs / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        path = ROOT / name if name.startswith("desiV3/") else catalogs / name.lstrip("*")
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
    manifest = json.loads((catalogs / "desiV3-consolidated-manifest.json").read_text())
    assert len(manifest["chunks"]) == 12 and sum(c["count"] for c in manifest["chunks"]) == manifest["count"]
    for cls, spec in zip(("DR1_GALAXY", "DR1_QSO"), manifest["catalogs"]):
        with gzip.open(candidate / "identity" / f"{cls}_targetid.npy.gz", "rb") as stream:
            ids = np.load(io.BytesIO(stream.read()), allow_pickle=False)
        assert ids.dtype == np.dtype("<u8") and len(ids) == spec["count"]
    for index, chunk in enumerate(manifest["chunks"]):
        for component in ("fog", "rsd"):
            path = candidate / "payloads" / f"{component}_{index:03d}.bin"
            values = np.frombuffer(gzip.decompress(path.read_bytes()), dtype="<f4")
            assert len(values) == chunk["count"] and np.isfinite(values).all(), path
    print(f"Verified 12 V3 chunks, 2 identity arrays, 24 correction payloads, {manifest['count']:,} rows.")


if __name__ == "__main__":
    main()
