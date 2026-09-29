#!/usr/bin/env python3
"""Convert DR1 .desi v3 files to the historical v2 layout without changing row order."""
from __future__ import annotations

import argparse
import struct
from pathlib import Path

BLOCK = 8 * 1024 * 1024


def compact(source: Path, output: Path) -> None:
    with source.open("rb") as input_stream:
        header = bytearray(input_stream.read(48))
        if len(header) != 48 or header[:4] != b"DESI":
            raise ValueError(f"invalid catalogue header: {source}")
        version, count = struct.unpack_from("<II", header, 4)
        if version != 3 or source.stat().st_size != 48 + count * 32:
            raise ValueError(f"expected complete v3 catalogue: {source}")
        struct.pack_into("<I", header, 4, 2)
        remaining = count * 12
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("wb") as result:
            result.write(header)
            while remaining:
                block = input_stream.read(min(BLOCK, remaining))
                if not block:
                    raise ValueError(f"short catalogue read: {source}")
                result.write(block)
                remaining -= len(block)
    if output.stat().st_size != 48 + count * 12:
        raise ValueError(f"short compact catalogue: {output}")
    print(f"{source.name}: {count:,} rows -> {output.stat().st_size:,} bytes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("assets"))
    parser.add_argument("--out-dir", type=Path, default=Path("assets/compact"))
    args = parser.parse_args()
    for name in ("DR1_GALAXY.desi", "DR1_QSO.desi"):
        compact(args.source_dir / name, args.out_dir / name)


if __name__ == "__main__":
    main()
