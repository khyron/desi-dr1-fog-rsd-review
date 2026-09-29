#!/usr/bin/env python3
"""Offline DESI DR1 large-scale RSD pipeline; never changes FoG files."""
from __future__ import annotations

import argparse
import json
import re
import struct
import gzip
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = ROOT / "data_external" / "RSD"
TRACERS = ("BGS", "LRG", "LRG_ELG", "ELG", "QSO")


def config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def local_sources(data_dir: Path) -> list[dict]:
    """Use actual downloaded v1.5 filenames; no guessed remote filenames."""
    out = []
    for path in sorted(data_dir.glob("*_clustering.*.fits")):
        m = re.match(r"(?P<name>.+)_(?P<region>NGC|SGC)(?:_(?P<random>\d+))?_clustering\.(?P<kind>dat|ran)\.fits$", path.name)
        if not m:
            continue
        name = m["name"]
        tracer = {"BGS_BRIGHT-21.5": "BGS", "LRG": "LRG", "LRG+ELG_LOPnotqso": "LRG_ELG", "ELG_LOPnotqso": "ELG", "QSO": "QSO"}.get(name)
        if not tracer:
            continue
        out.append({"tracer": tracer, "region": m["region"], "kind": "random" if m["kind"] == "ran" else "data",
                    "randomIndex": int(m["random"]) if m["random"] else None, "file": path.name,
                    "localPath": str(path), "bytes": path.stat().st_size,
                    "remotePath": "/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/" + path.name})
    return out


def read_desi_header(path: Path) -> dict:
    with path.open("rb") as f:
        head = f.read(48)
    if head[:4] != b"DESI":
        raise ValueError(f"{path}: invalid DESI magic")
    version, count, units, zmin, zmax = struct.unpack_from("<IIfff", head, 4)
    return {"version": version, "count": count, "unitsPerMpc": units, "zMin": zmin, "zMax": zmax}


def cmd_sources(args: argparse.Namespace) -> int:
    rows = local_sources(Path(args.data_dir))
    selected = [r for r in rows if not args.tracer or r["tracer"] == args.tracer]
    result = {"generatedAt": datetime.now(timezone.utc).isoformat(), "catalog": "DESI DR1 LSS iron LSScats v1.5",
              "files": selected, "totalBytes": sum(r["bytes"] for r in selected)}
    text = json.dumps(result, indent=2)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    return 0


def cmd_inspect(args: argparse.Namespace) -> int:
    from astropy.io import fits
    rows = []
    for source in local_sources(Path(args.data_dir)):
        if args.tracer and source["tracer"] != args.tracer:
            continue
        with fits.open(source["localPath"], memmap=True, lazy_load_hdus=True) as hdul:
            table = next(h for h in hdul if getattr(h.data, "names", None))
            data = table.data
            cols = list(data.names)
            missing = [c for c in ("RA", "DEC", "Z", "WEIGHT") if c not in cols]
            if missing:
                raise ValueError(f"{source['file']}: missing {missing}")
            sample = data[::max(1, len(data) // 100000)]
            rows.append({**source, "rows": len(data), "columns": cols,
                         "zRange": [float(np.nanmin(sample["Z"])), float(np.nanmax(sample["Z"]))],
                         "raRange": [float(np.nanmin(sample["RA"])), float(np.nanmax(sample["RA"]))],
                         "decRange": [float(np.nanmin(sample["DEC"])), float(np.nanmax(sample["DEC"]))],
                         "weightPercentiles": np.nanpercentile(sample["WEIGHT"], [0, 50, 100]).tolist()})
    print(json.dumps({"files": rows}, indent=2))
    return 0


def cmd_audit_viewer(args: argparse.Namespace) -> int:
    rows = {p.stem: read_desi_header(p) for p in sorted(Path(args.assets).glob("*.desi"))}
    print(json.dumps({"coordinate": "catalog equatorial x,y,z Mpc; scene x,z,y", "records": rows}, indent=2))
    return 0


def sky_to_cartesian(ra, dec, z, cosmo):
    """Viewer-compatible Planck18 RA/DEC/Z to Cartesian comoving Mpc."""
    r = np.asarray(cosmo.comoving_distance(z).value, dtype="f8")
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    cdec = np.cos(dec)
    return np.column_stack((r * cdec * np.cos(ra), r * cdec * np.sin(ra), r * np.sin(dec)))


def _physical_mpc(value_mpc_h: float, cosmo) -> float:
    """Convert published h^-1 Mpc parameters to the viewer's physical Mpc geometry."""
    return float(value_mpc_h) / float(cosmo.h)


def _fits_rows(path: Path, cosmo, limit: int, z_min: float, z_max: float, row_stride: int | None = None):
    """Read a deterministic FITS stride; no coordinate or weight convention is guessed."""
    from astropy.io import fits
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdul:
        table = next(h.data for h in hdul if getattr(h.data, "names", None))
        if row_stride is not None and row_stride < 1:
            raise ValueError("row_stride must be a positive integer")
        if row_stride is not None and limit:
            raise ValueError("specify either row_stride or limit, not both")
        step = int(row_stride) if row_stride is not None else (max(1, len(table) // limit) if limit else 1)
        row = table[::step]
        required = ("TARGETID", "RA", "DEC", "Z", "WEIGHT")
        missing = [name for name in required if name not in row.names]
        if missing:
            raise ValueError(f"{path.name}: missing {missing}")
        weight = np.asarray(row["WEIGHT"], dtype="f8")
        # WEIGHT is the catalogue's completeness/systematics weight. FKP is the supplied
        # reconstruction weighting and is applied only when the column is present.
        if "WEIGHT_FKP" in row.names:
            weight *= np.asarray(row["WEIGHT_FKP"], dtype="f8")
        z = np.asarray(row["Z"], dtype="f8")
        ra, dec = np.asarray(row["RA"], dtype="f8"), np.asarray(row["DEC"], dtype="f8")
        keep = (np.isfinite(weight) & (weight > 0) & np.isfinite(z) & np.isfinite(ra) &
                np.isfinite(dec) & (z >= z_min) & (z < z_max))
        return {
            "targetid": np.asarray(row["TARGETID"][keep], dtype="u8"),
            "position": sky_to_cartesian(ra[keep], dec[keep], z[keep], cosmo),
            "weight": weight[keep],
            "rowsRead": int(len(row)), "rowsKept": int(keep.sum()), "stride": int(step),
        }


def _fits_positions(path: Path, cosmo, limit: int, z_min: float = -np.inf, z_max: float = np.inf):
    rows = _fits_rows(path, cosmo, limit, z_min, z_max)
    return rows["position"], rows["weight"]


def _deduplicate_targetids(rows: dict) -> dict:
    """Keep one exact sky row per TARGETID, preferring the largest supplied total weight.

    The LRG+ELG combined clustering catalogue contains a small number of identical
    TARGETID/position rows selected by both tracer definitions. Counting both would
    double-weight one galaxy in the combined density field and would make the exact
    viewer-order export ambiguous. This is an identity deduplication, never a spatial
    merge. Ties preserve the earliest FITS row.
    """
    targetid, weight = rows["targetid"], rows["weight"]
    if len(targetid) < 2:
        rows["duplicatesDiscarded"] = 0
        return rows
    ranked = np.lexsort((-weight, targetid))
    sorted_ids = targetid[ranked]
    keep = np.empty(len(ranked), dtype=bool)
    keep[0] = True
    keep[1:] = sorted_ids[1:] != sorted_ids[:-1]
    winner = np.sort(ranked[keep], kind="stable")
    discarded = len(targetid) - len(winner)
    if not discarded:
        rows["duplicatesDiscarded"] = 0
        return rows
    result = {key: (value[winner] if isinstance(value, np.ndarray) and len(value) == len(targetid) else value)
              for key, value in rows.items()}
    result["rowsKept"] = int(len(winner))
    result["duplicatesDiscarded"] = int(discarded)
    return result


def cmd_prototype(args: argparse.Namespace) -> int:
    """Small, repeatable LRG-region field. Never writes production corrections."""
    from astropy.cosmology import Planck18
    from pyrecon import IterativeFFTReconstruction
    cfg = config(Path(args.config)); tracer = args.tracer or "LRG"
    if tracer not in cfg["tracers"]:
        raise ValueError("prototype requires configured tracer")
    files = local_sources(Path(args.data_dir))
    find = lambda kind: next(Path(x["localPath"]) for x in files if x["tracer"] == tracer and x["region"] == args.region and x["kind"] == kind)
    cosmo = Planck18
    p = cfg["tracers"][tracer]
    data, dw = _fits_positions(find("data"), cosmo, args.limit, p["zMin"], p["zMax"])
    randoms, rw = _fits_positions(find("random"), cosmo, args.limit * args.random_factor, p["zMin"], p["zMax"])
    cell_mpc = _physical_mpc(args.cell_size, cosmo)
    smoothing_mpc = _physical_mpc(p["smoothingMpcH"], cosmo)
    recon = IterativeFFTReconstruction(f=p["f"], bias=p["bias"], los=None,
        cellsize=cell_mpc, boxpad=1.2, positions=np.concatenate((data, randoms)))
    recon.assign_data(data, weights=dw)
    recon.assign_randoms(randoms, weights=rw)
    recon.set_density_contrast(smoothing_radius=smoothing_mpc)
    recon.run()
    shifts = recon.read_shifts(data[:min(len(data), 10000)], field="rsd")
    # pyrecon defines reconstructed position as x - read_shifts(x). The scalar stored by
    # the viewer is therefore delta_r = -dot(rsd_shift, rhat), not its opposite.
    radial = -np.sum(shifts * data[:len(shifts)] / np.linalg.norm(data[:len(shifts)], axis=1)[:, None], axis=1)
    report = {"tracer": tracer, "region": args.region, "data": len(data), "randoms": len(randoms),
              "cellSizeMpcH": args.cell_size, "cellSizeMpc": cell_mpc,
              "smoothingMpcH": p["smoothingMpcH"], "smoothingMpc": smoothing_mpc,
              "field": "rsd", "shiftAbsPercentilesMpc": np.percentile(np.linalg.norm(shifts, axis=1), [50,90,99]).tolist(),
              "deltaRadialPercentilesMpc": np.percentile(radial, [1,50,99]).tolist(), "finite": bool(np.isfinite(shifts).all())}
    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2)); return 0


def _source_path(files: list[dict], tracer: str, region: str, kind: str) -> Path:
    candidates = [Path(row["localPath"]) for row in files if row["tracer"] == tracer and
                  row["region"] == region and row["kind"] == kind]
    if not candidates:
        raise FileNotFoundError(f"missing {tracer} {region} {kind} clustering FITS")
    # The source manifest records all random indices. The reconstruction command currently
    # uses random-0 deliberately, rather than silently selecting a random subset.
    return candidates[0]


def _shift_statistics(shifts: np.ndarray, positions: np.ndarray) -> dict:
    radius = np.linalg.norm(positions, axis=1)
    radial_vec = np.sum(shifts * positions / radius[:, None], axis=1)
    transverse = np.sqrt(np.maximum(0, np.sum(shifts * shifts, axis=1) - radial_vec * radial_vec))
    magnitude = np.linalg.norm(shifts, axis=1)
    return {
        "shiftAbsMpc": {"mean": float(magnitude.mean()), "percentiles": np.percentile(magnitude, [50, 75, 90, 95, 99, 100]).tolist()},
        "radialAbsMpc": {"mean": float(np.abs(radial_vec).mean()), "percentiles": np.percentile(np.abs(radial_vec), [50, 75, 90, 95, 99, 100]).tolist()},
        "transverseAbsMpc": {"mean": float(transverse.mean()), "percentiles": np.percentile(transverse, [50, 90, 99, 100]).tolist()},
    }


def _query_radial_delta(recon, positions: np.ndarray, batch_size: int) -> tuple[np.ndarray, dict]:
    """Query only field='rsd' in bounded NumPy batches and return viewer-Mpc deltas."""
    out = np.empty(len(positions), dtype="f4")
    diagnostics = []
    for start in range(0, len(positions), batch_size):
        end = min(len(positions), start + batch_size)
        pos = positions[start:end]
        shifts = np.asarray(recon.read_shifts(pos, field="rsd"), dtype="f8")
        if not np.isfinite(shifts).all():
            raise ValueError(f"non-finite RSD shift in query rows {start}:{end}")
        radii = np.linalg.norm(pos, axis=1)
        if np.any(radii == 0):
            raise ValueError("zero-radius reconstruction query")
        out[start:end] = -np.sum(shifts * pos / radii[:, None], axis=1)
        if len(diagnostics) < 4:
            diagnostics.append(_shift_statistics(shifts, pos))
    if not np.isfinite(out).all():
        raise ValueError("non-finite projected RSD delta")
    return out, {"batches": int(np.ceil(len(positions) / batch_size)), "samples": diagnostics}


def cmd_reconstruct(args: argparse.Namespace) -> int:
    """Build per-tracer/per-region reusable TARGETID query caches from official LSS data."""
    from astropy.cosmology import Planck18
    from pyrecon import IterativeFFTReconstruction

    cfg, cosmo = config(Path(args.config)), Planck18
    files, out = local_sources(Path(args.data_dir)), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    chosen = [args.tracer] if args.tracer else list(cfg["tracers"])
    manifest_path = out / "reconstruction-manifest.json"
    previous_reports = []
    if manifest_path.exists():
        try:
            previous_reports = json.loads(manifest_path.read_text(encoding="utf-8")).get("reports", [])
        except (OSError, json.JSONDecodeError):
            previous_reports = []
    reports = []
    for tracer in chosen:
        p = cfg["tracers"][tracer]
        cell_mpc = _physical_mpc(args.cell_size or cfg["gridCellSizeMpcH"], cosmo)
        smoothing_mpc = _physical_mpc(p["smoothingMpcH"], cosmo)
        for region in ("NGC", "SGC"):
            data_path = _source_path(files, tracer, region, "data")
            random_path = _source_path(files, tracer, region, "random")
            data = _deduplicate_targetids(_fits_rows(data_path, cosmo, args.limit, p["zMin"], p["zMax"]))
            randoms = _deduplicate_targetids(_fits_rows(random_path, cosmo, args.random_limit, p["zMin"], p["zMax"]))
            if not len(data["position"]) or not len(randoms["position"]):
                raise ValueError(f"{tracer} {region}: empty valid data/random selection")
            recon = IterativeFFTReconstruction(f=p["f"], bias=p["bias"], los=None,
                cellsize=cell_mpc, boxpad=1.2,
                positions=np.concatenate((data["position"], randoms["position"])))
            recon.assign_data(data["position"], weights=data["weight"])
            recon.assign_randoms(randoms["position"], weights=randoms["weight"])
            recon.set_density_contrast(smoothing_radius=smoothing_mpc)
            recon.run(niterations=args.iterations)
            delta, query = _query_radial_delta(recon, data["position"], cfg["batchSize"])
            cache = out / f"{tracer}_{region}_rsd_query_v1.npz"
            np.savez(cache, targetid=data["targetid"], delta_mpc=delta)
            report = {
                "tracer": tracer, "region": region, "catalog": cfg["catalog"], "field": "rsd", "los": "local",
                "data": {k: data[k] for k in ("rowsRead", "rowsKept", "stride", "duplicatesDiscarded")},
                "randoms": {k: randoms[k] for k in ("rowsRead", "rowsKept", "stride", "duplicatesDiscarded")},
                "zRange": [p["zMin"], p["zMax"]], "bias": p["bias"], "f": p["f"],
                "cellSizeMpcH": args.cell_size or cfg["gridCellSizeMpcH"], "cellSizeMpc": cell_mpc,
                "smoothingMpcH": p["smoothingMpcH"], "smoothingMpc": smoothing_mpc,
                "mesh": {"nmesh": np.asarray(recon.nmesh).tolist(), "boxsizeMpc": np.asarray(recon.boxsize).tolist(), "boxcenterMpc": np.asarray(recon.boxcenter).tolist()},
                "query": query, "shift": query["samples"][0] if query["samples"] else {},
                "deltaPercentilesMpc": np.percentile(delta, [0, 1, 50, 99, 100]).tolist(),
                "cache": str(cache),
            }
            (out / f"{tracer}_{region}_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            reports.append(report)
            print(json.dumps({"tracer": tracer, "region": region, "data": len(data["position"]), "randoms": len(randoms["position"]), "cache": str(cache)}, indent=2))
    changed = {(report["tracer"], report["region"]) for report in reports}
    reports = [report for report in previous_reports if (report.get("tracer"), report.get("region")) not in changed] + reports
    manifest_path.write_text(json.dumps({"generatedAt": datetime.now(timezone.utc).isoformat(), "reports": reports}, indent=2) + "\n", encoding="utf-8")
    return 0


def cmd_preflight(args: argparse.Namespace) -> int:
    """Estimate per-region FFT size from deterministic samples before a full run."""
    from astropy.cosmology import Planck18
    cfg, cosmo = config(Path(args.config)), Planck18
    files, chosen, rows = local_sources(Path(args.data_dir)), ([args.tracer] if args.tracer else list(cfg["tracers"])), []
    for tracer in chosen:
        p = cfg["tracers"][tracer]
        cell_mpc = _physical_mpc(args.cell_size or cfg["gridCellSizeMpcH"], cosmo)
        for region in ("NGC", "SGC"):
            data = _fits_rows(_source_path(files, tracer, region, "data"), cosmo, args.sample, p["zMin"], p["zMax"])
            randoms = _fits_rows(_source_path(files, tracer, region, "random"), cosmo, args.sample, p["zMin"], p["zMax"])
            positions = np.concatenate((data["position"], randoms["position"]))
            lo, hi = positions.min(axis=0), positions.max(axis=0)
            box = (hi - lo) * 1.2
            nmesh = np.maximum(1, np.ceil(box / cell_mpc).astype("i8"))
            cells = int(np.prod(nmesh, dtype="i8"))
            # Lower bound: several real/complex FFT work arrays, before positions and MPI overhead.
            lower_gib = cells * 48 / 1024 ** 3
            rows.append({"tracer": tracer, "region": region, "sampleData": data["rowsKept"], "sampleRandoms": randoms["rowsKept"],
                         "cellSizeMpc": cell_mpc, "estimatedBoxMpc": box.tolist(), "estimatedNmesh": nmesh.tolist(),
                         "estimatedCells": cells, "minimumWorkingGiB": lower_gib})
    lower = max((row["minimumWorkingGiB"] for row in rows), default=0)
    result = {"note": "Sampled bounds; memory estimate is a lower bound, not authorization to run.", "fields": rows,
              "maximumMinimumWorkingGiB": lower, "recommendedAvailableGiB": lower * 3}
    print(json.dumps(result, indent=2)); return 0


def cmd_collect_reports(args: argparse.Namespace) -> int:
    """Rebuild the durable reconstruction manifest from completed per-region reports."""
    out = Path(args.out)
    reports = []
    for path in sorted(out.glob("*_report.json")):
        report = json.loads(path.read_text(encoding="utf-8"))
        cache = Path(report.get("cache", ""))
        if not cache.exists():
            raise FileNotFoundError(f"report without cache: {path}")
        reports.append(report)
    if len(reports) != args.expected_fields:
        raise ValueError(f"expected {args.expected_fields} region reports, found {len(reports)}")
    manifest = {"generatedAt": datetime.now(timezone.utc).isoformat(), "profile": args.profile,
                "reports": reports}
    (out / "reconstruction-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"reports": len(reports), "profile": args.profile}, indent=2)); return 0


def cmd_materialize(args: argparse.Namespace) -> int:
    """Project validated TARGETID query caches into the exact private point-chunk ordering."""
    cache_dir, target_dir, out = Path(args.cache_dir), Path(args.target_index_dir), Path(args.out)
    files = sorted(cache_dir.glob("*_rsd_query_v1.npz"))
    if not files:
        raise FileNotFoundError(f"no query caches in {cache_dir}")
    by_field: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
    for path in files:
        with np.load(path) as cache:
            field = path.name.split("_")[0] if not path.name.startswith("LRG_ELG_") else "LRG_ELG"
            by_field.setdefault(field, []).append((np.asarray(cache["targetid"], dtype="u8"),
                                                    np.asarray(cache["delta_mpc"], dtype="f4")))
    lookup = {}
    for field, chunks in by_field.items():
        ids = np.concatenate([chunk[0] for chunk in chunks])
        deltas = np.concatenate([chunk[1] for chunk in chunks])
        order = np.argsort(ids, kind="stable"); ids, deltas = ids[order], deltas[order]
        if np.any(ids[1:] == ids[:-1]):
            raise ValueError(f"{field}: duplicate TARGETID inside one reconstruction field")
        lookup[field] = (ids, deltas)

    def fill_from(field: str, target: np.ndarray, allowed: np.ndarray, result: np.ndarray, matched_mask: np.ndarray):
        if field not in lookup or not np.any(allowed):
            return
        ids, deltas = lookup[field]
        subset = target[allowed]
        pos = np.searchsorted(ids, subset)
        hit = (pos < len(ids)) & (ids[np.minimum(pos, len(ids) - 1)] == subset)
        indices = np.flatnonzero(allowed)
        result[indices[hit]] = deltas[pos[hit]]
        matched_mask[indices[hit]] = True

    source_manifest = json.loads((Path(args.chunks) / "manifest.json").read_text(encoding="utf-8"))
    if len(source_manifest.get("chunks", [])) != len(list(target_dir.glob("targetid_*.npy"))):
        raise ValueError("target index count differs from source chunk manifest")
    out.mkdir(parents=True, exist_ok=True)
    summary, total, matched = [], 0, 0
    for index, path in enumerate(sorted(target_dir.glob("targetid_*.npy"))):
        target = np.asarray(np.load(path, mmap_mode="r"), dtype="u8")
        type_file = Path(args.chunks) / source_manifest["chunks"][index]["type"]
        tracer_type = np.frombuffer(gzip.decompress(type_file.read_bytes()), dtype="u1")
        if len(tracer_type) != len(target):
            raise ValueError(f"chunk {index}: type/target count mismatch")
        delta, hit = np.zeros(len(target), dtype="f4"), np.zeros(len(target), dtype=bool)
        fill_from("BGS", target, tracer_type == 0, delta, hit)
        fill_from("LRG", target, (tracer_type == 1) & ~hit, delta, hit)
        fill_from("ELG", target, tracer_type == 2, delta, hit)
        # The combined field is the deliberate 0.8--1.1 fallback for either LRG or ELG.
        fill_from("LRG_ELG", target, ((tracer_type == 1) | (tracer_type == 2)) & ~hit, delta, hit)
        fill_from("QSO", target, tracer_type == 3, delta, hit)
        if not np.isfinite(delta).all():
            raise ValueError(f"{path}: non-finite materialized delta")
        np.save(out / f"rsd_{index:03d}.npy", delta)
        total += len(target); matched += int(hit.sum())
        summary.append({"chunk": index, "count": int(len(target)), "matched": int(hit.sum()), "nonzero": int(np.count_nonzero(delta))})
    values = np.concatenate([value[1] for value in lookup.values()])
    report = {"total": total, "matched": matched, "coverage": matched / total if total else 0, "chunks": summary,
              "deltaPercentilesMpc": np.percentile(values, [0, 1, 50, 99, 100]).tolist(),
              "identityPolicy": "exact TARGETID plus published viewer tracer byte; LRG_ELG only fallback for LRG/ELG"}
    (out / "materialize-report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if report["coverage"] < args.min_coverage:
        raise ValueError(f"RSD TARGETID coverage {report['coverage']:.4%} below required {args.min_coverage:.4%}")
    print(json.dumps(report, indent=2)); return 0


def cmd_validate_chunks(args: argparse.Namespace) -> int:
    """Static validation for production RSD assets before any PlayCanvas upload."""
    source = json.loads((Path(args.source_chunks) / "manifest.json").read_text())
    out = Path(args.rsd_chunks)
    manifest = json.loads((out / "manifest.json").read_text())
    if len(source["chunks"]) != len(manifest["chunks"]):
        raise ValueError("RSD manifest chunk count differs from source geometry")
    rows, all_values = [], []
    for index, (base, rsd) in enumerate(zip(source["chunks"], manifest["chunks"])):
        if base["count"] != rsd["count"] or not rsd.get("rsd"):
            raise ValueError(f"chunk {index}: missing RSD file or count mismatch")
        payload = gzip.decompress((out / rsd["rsd"]).read_bytes())
        if len(payload) != base["count"] * 4:
            raise ValueError(f"chunk {index}: payload bytes {len(payload)}")
        delta = np.frombuffer(payload, dtype="<f4")
        if not np.isfinite(delta).all():
            raise ValueError(f"chunk {index}: non-finite RSD values")
        rows.append({"chunk": index, "count": int(len(delta)), "nonzero": int(np.count_nonzero(delta)),
                     "positive": int((delta > 0).sum()), "negative": int((delta < 0).sum())})
        all_values.append(delta)
    values = np.concatenate(all_values)
    report = {"objects": int(len(values)), "valid": int(np.isfinite(values).sum()), "nonzero": int(np.count_nonzero(values)),
              "zero": int((values == 0).sum()), "coverage": float(np.count_nonzero(values) / len(values)),
              "meanMpc": float(values.mean()), "medianMpc": float(np.median(values)),
              "absPercentilesMpc": np.percentile(np.abs(values), [50, 75, 90, 95, 99, 100]).tolist(),
              "positive": int((values > 0).sum()), "negative": int((values < 0).sum()), "chunks": rows}
    if report["coverage"] < args.min_coverage:
        raise ValueError(f"RSD coverage {report['coverage']:.4%} below required {args.min_coverage:.4%}")
    (out / "validation-report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2)); return 0


def cmd_validate_coherence(args: argparse.Namespace) -> int:
    """Check that neighboring corrected points vary more smoothly than a shuffled control."""
    from scipy.spatial import cKDTree

    source = Path(args.source_chunks)
    manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    rsd = Path(args.rsd_chunks)
    positions, deltas = [], []
    for index, spec in enumerate(manifest["chunks"]):
        count = int(spec["count"])
        raw = gzip.decompress((source / spec["file"]).read_bytes())
        encoded = np.frombuffer(raw, dtype="<i2")
        if len(encoded) != count * 3:
            raise ValueError(f"chunk {index}: invalid source geometry payload")
        xyz = np.column_stack([np.cumsum(encoded[axis * count:(axis + 1) * count], dtype="i4")
                               for axis in range(3)]).astype("f4") / float(manifest["unitsPerMpc"])
        delta = np.frombuffer(gzip.decompress((rsd / f"rsd_{index:03d}.bin").read_bytes()), dtype="<f4")
        useful = np.flatnonzero(delta != 0)
        if len(useful) > args.samples_per_chunk:
            step = len(useful) / args.samples_per_chunk
            useful = useful[np.floor(np.arange(args.samples_per_chunk) * step).astype("i8")]
        positions.append(xyz[useful]); deltas.append(delta[useful])
    positions, deltas = np.concatenate(positions), np.concatenate(deltas)
    tree = cKDTree(positions)
    rng = np.random.default_rng(20260911)
    shuffled = rng.permutation(deltas)
    bands = []
    for radius in args.radii:
        pairs = tree.query_pairs(radius, output_type="ndarray")
        if not len(pairs):
            bands.append({"radiusMpc": radius, "pairs": 0})
            continue
        a, b = pairs[:, 0], pairs[:, 1]
        actual_diff = np.abs(deltas[a] - deltas[b])
        shuffled_diff = np.abs(shuffled[a] - shuffled[b])
        corr = float(np.corrcoef(deltas[a], deltas[b])[0, 1]) if len(pairs) > 1 else None
        bands.append({"radiusMpc": radius, "pairs": int(len(pairs)), "pearson": corr,
                      "medianAbsDifferenceMpc": float(np.median(actual_diff)),
                      "shuffledMedianAbsDifferenceMpc": float(np.median(shuffled_diff)),
                      "smoothnessRatio": float(np.median(actual_diff) / max(np.median(shuffled_diff), 1e-12))})
    report = {"samples": int(len(deltas)), "method": "same-position radial deltas against deterministic shuffle control",
              "bands": bands}
    (rsd / "coherence-report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2)); return 0


def cmd_build_chunks(args: argparse.Namespace) -> int:
    """Write only already-validated per-chunk Float32 Mpc deltas; never touches FoG."""
    manifest = json.loads(Path(args.chunks).joinpath("manifest.json").read_text())
    source, out = Path(args.delta_dir), Path(args.out); out.mkdir(parents=True, exist_ok=True)
    total = nonzero = 0
    for index, chunk in enumerate(manifest["chunks"]):
        delta_path = source / f"rsd_{index:03d}.npy"
        if not delta_path.exists(): raise FileNotFoundError(delta_path)
        delta = np.load(delta_path, mmap_mode="r")
        if delta.dtype != np.float32 or delta.ndim != 1 or len(delta) != chunk["count"]: raise ValueError(f"{delta_path}: count/dtype")
        if not np.isfinite(delta).all(): raise ValueError(f"{delta_path}: non-finite")
        name = f"rsd_{index:03d}.bin"; blob = gzip.compress(delta.tobytes(), compresslevel=9); (out / name).write_bytes(blob)
        chunk["rsd"], chunk["rsdBytes"] = name, len(blob); total += len(blob); nonzero += int(np.count_nonzero(delta))
    manifest["rsd"] = {
        "version": 1, "format": "gzip-float32le-delta-mpc-v1", "algorithm": "iterativeFFT-rsd-v1",
        "field": "rsd", "lineOfSight": "local", "bytes": total, "nonzero": nonzero,
        "source": "DESI DR1 LSS iron LSScats v1.5", "cosmology": "Planck18 physical comoving Mpc",
        "independentOf": "fog"
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    # This is deliberately an overlay, not a mutation of chunks_fog_v1/manifest.json.
    # Copy its fields into the deployed geometry manifest and keep rsdBaseUrl relative to it.
    overlay = {"rsdBaseUrl": args.runtime_rsd_base_url, "rsd": manifest["rsd"],
               "chunks": [{"rsd": c["rsd"], "rsdBytes": c["rsdBytes"]} for c in manifest["chunks"]]}
    (out / "manifest-rsd-overlay.json").write_text(json.dumps(overlay, indent=1) + "\n")
    print(json.dumps({"chunks":len(manifest["chunks"]),"bytes":total,"nonzero":nonzero}, indent=2)); return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=str(DEFAULT_DATA))
    parser.add_argument("--config", default=str(Path(__file__).with_name("rsd_config.json")))
    parser.add_argument("--tracer", choices=TRACERS)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("sources"); p.add_argument("--out"); p.set_defaults(func=cmd_sources)
    sub.add_parser("inspect").set_defaults(func=cmd_inspect)
    p = sub.add_parser("audit-viewer"); p.add_argument("--assets", default=str(ROOT / "assets")); p.set_defaults(func=cmd_audit_viewer)
    p = sub.add_parser("prototype"); p.add_argument("--region", choices=("NGC", "SGC"), default="NGC"); p.add_argument("--limit", type=int, default=100000); p.add_argument("--random-factor", type=int, default=5); p.add_argument("--cell-size", type=float, default=4.0); p.add_argument("--out", default=str(DEFAULT_DATA / "prototype-lrg.json")); p.set_defaults(func=cmd_prototype)
    p = sub.add_parser("reconstruct", help="build per-region reusable TARGETID RSD query caches")
    p.add_argument("--out", default=str(DEFAULT_DATA / "derived" / "rsd-v1"))
    p.add_argument("--limit", type=int, default=0, help="deterministic FITS data limit; 0 means full catalog")
    p.add_argument("--random-limit", type=int, default=0, help="deterministic FITS random limit; 0 means full catalog")
    p.add_argument("--cell-size", type=float, default=0, help="h^-1 Mpc; 0 uses rsd_config.json")
    p.add_argument("--iterations", type=int, default=3)
    p.set_defaults(func=cmd_reconstruct)
    p = sub.add_parser("preflight", help="estimate grid and lower-bound memory from sampled FITS bounds")
    p.add_argument("--sample", type=int, default=10000)
    p.add_argument("--cell-size", type=float, default=0, help="h^-1 Mpc; 0 uses rsd_config.json")
    p.set_defaults(func=cmd_preflight)
    p = sub.add_parser("collect-reports", help="rebuild reconstruction manifest from completed region reports")
    p.add_argument("--out", default=str(DEFAULT_DATA / "derived" / "rsd-v1"))
    p.add_argument("--expected-fields", type=int, default=10)
    p.add_argument("--profile", default="exploratory")
    p.set_defaults(func=cmd_collect_reports)
    p = sub.add_parser("materialize", help="match query caches to exact offline chunk TARGETID order")
    p.add_argument("--cache-dir", default=str(DEFAULT_DATA / "derived" / "rsd-v1"))
    p.add_argument("--target-index-dir", default=str(ROOT / "data_external" / "fog" / "chunk-targetids"))
    p.add_argument("--chunks", default=str(ROOT / "ply" / "chunks_fog_v1"))
    p.add_argument("--out", default=str(DEFAULT_DATA / "derived" / "rsd-v1" / "chunk-deltas"))
    p.add_argument("--min-coverage", type=float, default=0.10)
    p.set_defaults(func=cmd_materialize)
    p = sub.add_parser("build-chunks"); p.add_argument("--chunks", default=str(ROOT / "ply" / "chunks_fog_v1")); p.add_argument("--delta-dir", required=True); p.add_argument("--out", default=str(ROOT / "ply" / "chunks_rsd_v1")); p.add_argument("--runtime-rsd-base-url", default="../chunks_rsd_v1/"); p.set_defaults(func=cmd_build_chunks)
    p = sub.add_parser("validate-chunks", help="verify RSD payloads and coverage before deployment")
    p.add_argument("--source-chunks", default=str(ROOT / "ply" / "chunks_fog_v1"))
    p.add_argument("--rsd-chunks", default=str(ROOT / "ply" / "chunks_rsd_v1"))
    p.add_argument("--min-coverage", type=float, default=0.10)
    p.set_defaults(func=cmd_validate_chunks)
    p = sub.add_parser("validate-coherence", help="compare nearby RSD deltas with a shuffled control")
    p.add_argument("--source-chunks", default=str(ROOT / "ply" / "chunks_fog_v1"))
    p.add_argument("--rsd-chunks", default=str(ROOT / "ply" / "chunks_rsd_v1"))
    p.add_argument("--samples-per-chunk", type=int, default=4000)
    p.add_argument("--radii", type=float, nargs="+", default=[15.0, 30.0, 60.0])
    p.set_defaults(func=cmd_validate_coherence)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
