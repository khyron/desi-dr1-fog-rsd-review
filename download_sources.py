#!/usr/bin/env python3
"""Download selected public DESI inputs listed in sources.json.

This is an offline acquisition helper. It does not run a reconstruction.
The Gfinder downloader in tools/fog/fog.py additionally verifies official SHA256.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
BLOCK = 8 * 1024 * 1024


def remote_size(url: str) -> int:
    with urlopen(Request(url, method="HEAD"), timeout=60) as response:
        return int(response.headers["Content-Length"])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(BLOCK):
            digest.update(block)
    return digest.hexdigest()


def download(entry: dict, destination: Path, attempts: int = 8) -> dict:
    url = entry["url"]
    expected = remote_size(url)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    if destination.is_file() and destination.stat().st_size == expected:
        path = destination
    else:
        if partial.exists() and partial.stat().st_size > expected:
            partial.unlink()
        offset = partial.stat().st_size if partial.exists() else 0
        failures = 0
        while offset < expected:
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            try:
                with urlopen(Request(url, headers=headers), timeout=120) as response:
                    if offset and response.status != 206:
                        raise RuntimeError(f"server ignored Range for {url}")
                    if offset and not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                        raise RuntimeError(f"bad Content-Range for {url}")
                    with partial.open("ab" if offset else "wb") as output:
                        while block := response.read(BLOCK):
                            output.write(block)
                            offset += len(block)
                failures = 0
            except (HTTPError, URLError, TimeoutError, OSError) as error:
                failures += 1
                if failures >= attempts:
                    raise RuntimeError(f"download failed at byte {offset}: {url}") from error
                time.sleep(min(60, 2 ** failures))
        if partial.stat().st_size != expected:
            raise RuntimeError(f"size mismatch: {partial}")
        partial.replace(destination)
        path = destination
    result = {"file": entry["file"], "url": url, "bytes": path.stat().st_size, "sha256": sha256(path)}
    if result["bytes"] != expected:
        raise RuntimeError(f"size mismatch: {path}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", required=True, choices=(
        "dr1_zcatalog", "fog_gfinder", "rsd_base", "bgs_lowz", "rsd_validation_index1"))
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "sources.json").read_text(encoding="utf-8-sig"))
    entries = [entry for entry in manifest["files"] if entry["group"] == args.group]
    results = []
    for entry in entries:
        destination = args.dest / entry["file"]
        if args.dry_run:
            print(f"{entry['url']} -> {destination}")
        else:
            result = download(entry, destination)
            results.append(result)
            print(f"verified {result['file']}: {result['bytes']} bytes SHA256 {result['sha256']}", flush=True)
    if not args.dry_run:
        (args.dest / "download-receipt.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
