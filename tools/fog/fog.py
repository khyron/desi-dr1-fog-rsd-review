#!/usr/bin/env python3
"""Acquire and verify the official DESI DR1 Gfinder v1.0 inputs for offline FoG work.

Raw FITS files are deliberately kept outside version control under data_external/.
This command never loads a FITS table into memory; later pipeline stages open them with
astropy memmap.  Usage (PowerShell):

  $env:PYTHONPATH='D:\desi\.python-packages'
  & $python tools/fog/fog.py download
  & $python tools/fog/fog.py verify
  & $python tools/fog/fog.py inspect
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = Path(__file__).with_name("data_sources.json")
DEFAULT_DATA = ROOT / "data_external" / "desi" / "dr1" / "gfinder" / "v1.0"
BLOCK = 8 * 1024 * 1024
USER_AGENT = "desi-fog-offline/1.0"
MAX_RETRIES = 12


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def url_for(manifest: dict, name: str) -> str:
    return "https://" + manifest["host"] + manifest["basePath"] + "/" + name


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while block := source.read(BLOCK):
            digest.update(block)
    return digest.hexdigest()


def expected_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.strip().replace("*", "").split()
        if len(fields) >= 2 and len(fields[0]) == 64:
            checksums[Path(fields[-1]).name] = fields[0].lower()
    return checksums


def download_one(url: str, destination: Path, dry_run: bool) -> None:
    if dry_run:
        print("[dry-run] " + url + " -> " + str(destination))
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    try:
        with urlopen(Request(url, method="HEAD", headers={"User-Agent": USER_AGENT}), timeout=60) as head:
            total_bytes = int(head.headers["Content-Length"])
    except (HTTPError, URLError, KeyError, ValueError) as error:
        raise RuntimeError("could not determine remote size for " + url + ": " + str(error)) from error
    if partial.exists() and partial.stat().st_size > total_bytes:
        partial.unlink()
    offset = partial.stat().st_size if partial.exists() else 0
    attempts = 0
    retry_count = 0
    while offset < total_bytes:
        attempts += 1
        if attempts > 1024:
            raise RuntimeError("too many interrupted responses for " + url)
        headers = {"User-Agent": USER_AGENT}
        if offset:
            headers["Range"] = "bytes=" + str(offset) + "-"
        try:
            response = urlopen(Request(url, headers=headers), timeout=60)
        except HTTPError as error:
            # The public DESI host occasionally returns 503 while serving these
            # multi-gigabyte files.  Preserve the partial file and retry the same
            # byte range; never treat it as a completed download.
            if error.code not in (408, 429, 500, 502, 503, 504) or retry_count >= MAX_RETRIES:
                raise RuntimeError("download failed for " + url + ": " + str(error)) from error
            retry_count += 1
            delay = min(300, 5 * 2 ** (retry_count - 1))
            print("  HTTP {0}; retrying byte {1} in {2}s ({3}/{4})".format(
                error.code, offset, delay, retry_count, MAX_RETRIES), flush=True)
            time.sleep(delay)
            continue
        except URLError as error:
            if retry_count >= MAX_RETRIES:
                raise RuntimeError("download failed for " + url + ": " + str(error)) from error
            retry_count += 1
            delay = min(300, 5 * 2 ** (retry_count - 1))
            print("  network error; retrying byte {0} in {1}s ({2}/{3})".format(
                offset, delay, retry_count, MAX_RETRIES), flush=True)
            time.sleep(delay)
            continue
        status = getattr(response, "status", None)
        if offset and status != 206:
            response.close()
            raise RuntimeError("server ignored Range request at byte " + str(offset) + " for " + url)
        before = offset
        with response, partial.open("ab" if offset else "wb") as output:
            while block := response.read(BLOCK):
                output.write(block)
                offset += len(block)
                print("\r  {0:.1f}%  {1:.2f}/{2:.2f} GiB".format(
                    offset * 100 / total_bytes, offset / 2**30, total_bytes / 2**30), end="", flush=True)
        print()
        if offset <= before:
            raise RuntimeError("response made no progress for " + url)
        retry_count = 0
    partial.replace(destination)


def record_state(data_dir: Path, manifest: dict, checksums: dict[str, str]) -> None:
    rows = []
    for item in manifest["files"]:
        path = data_dir / item["name"]
        rows.append({
            "filename": item["name"],
            "url": url_for(manifest, item["name"]),
            "approximateBytes": item["approximateBytes"],
            "actualBytes": path.stat().st_size if path.exists() else None,
            "expectedSha256": checksums.get(item["name"]),
            "actualSha256": sha256(path) if path.exists() else None,
        })
    (data_dir / "download-state.json").write_text(json.dumps({
        "dataset": manifest["dataset"], "version": manifest["version"],
        "generatedAt": datetime.now(timezone.utc).isoformat(), "files": rows
    }, indent=2), encoding="utf-8")


def command_download(args: argparse.Namespace) -> int:
    manifest, data_dir = load_manifest(), Path(args.data_dir)
    checksum_path = data_dir / manifest["checksumFile"]
    if not checksum_path.exists() or args.force:
        print("[get] " + manifest["checksumFile"])
        download_one(url_for(manifest, manifest["checksumFile"]), checksum_path, args.dry_run)
    checksums = expected_checksums(checksum_path) if checksum_path.exists() else {}
    for item in manifest["files"]:
        path = data_dir / item["name"]
        if path.exists() and item["name"] in checksums and sha256(path) == checksums[item["name"]]:
            print("[valid cache] " + item["name"])
            continue
        print("[get] " + item["name"])
        download_one(url_for(manifest, item["name"]), path, args.dry_run)
        if not args.dry_run and item["name"] in checksums and sha256(path) != checksums[item["name"]]:
            raise RuntimeError("SHA256 mismatch: " + item["name"])
    if not args.dry_run:
        record_state(data_dir, manifest, checksums)
    return 0


def command_verify(args: argparse.Namespace) -> int:
    manifest, data_dir = load_manifest(), Path(args.data_dir)
    checksum_path = data_dir / manifest["checksumFile"]
    if not checksum_path.exists():
        raise RuntimeError("missing checksum file: " + str(checksum_path))
    checksums, failed = expected_checksums(checksum_path), False
    for item in manifest["files"]:
        path = data_dir / item["name"]
        actual = sha256(path) if path.exists() else None
        ok = actual is not None and actual == checksums.get(item["name"])
        print(("[ok]   " if ok else "[FAIL] ") + item["name"])
        failed |= not ok
    if not failed:
        record_state(data_dir, manifest, checksums)
    return 1 if failed else 0


def command_inspect(args: argparse.Namespace) -> int:
    from astropy.io import fits
    manifest, data_dir = load_manifest(), Path(args.data_dir)
    for item in manifest["files"]:
        path = data_dir / item["name"]
        if not path.exists():
            raise RuntimeError("missing FITS: " + str(path))
        with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
            table = next((h for h in hdul if getattr(h, "data", None) is not None and getattr(h.data, "names", None)), None)
            if table is None:
                raise RuntimeError("no table HDU in " + item["name"])
            columns = set(table.data.names)
            missing = set(item["requiredColumns"]) - columns
            print(item["name"] + ": rows=" + str(len(table.data)) + ", HDU=" + str(table.name))
            print("  columns: " + ", ".join(table.data.names))
            if missing:
                raise RuntimeError("missing required columns in " + item["name"] + ": " + ", ".join(sorted(missing)))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA))
    sub = parser.add_subparsers(dest="command", required=True)
    download = sub.add_parser("download")
    download.add_argument("--dry-run", action="store_true")
    download.add_argument("--force", action="store_true")
    sub.add_parser("verify")
    sub.add_parser("inspect")
    args = parser.parse_args()
    return {"download": command_download, "verify": command_verify, "inspect": command_inspect}[args.command](args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print("fog: " + str(error), file=sys.stderr)
        raise SystemExit(1)
