#!/usr/bin/env python3
"""Run an isolated BGS NGC RSD query/coverage prototype; never writes payloads.

Reconstructs with a deterministic 20% data holdout and all official randoms.
Random-cloud support is a conservative local 3-D proximity diagnostic calibrated
on a disjoint random subset. It is not a substitute for the DESI angular mask.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from rsd_pipeline import (  # noqa: E402
    _deduplicate_targetids,
    _fits_rows,
    _physical_mpc,
    _query_radial_delta,
    _source_path,
    config,
    local_sources,
    read_desi_header,
)


def cache_lookup(path: Path):
    with np.load(path) as cache:
        ids = np.asarray(cache["targetid"], dtype="u8")
        delta = np.asarray(cache["delta_mpc"], dtype="f4")
    order = np.argsort(ids, kind="stable")
    ids, delta = ids[order], delta[order]
    if np.any(ids[1:] == ids[:-1]) or not np.isfinite(delta).all():
        raise ValueError("reference BGS NGC cache has duplicate IDs or non-finite deltas")
    return ids, delta


def lookup(ids: np.ndarray, source_ids: np.ndarray, source_delta: np.ndarray):
    at = np.searchsorted(source_ids, ids)
    hit = (at < len(source_ids)) & (source_ids[np.minimum(at, len(source_ids) - 1)] == ids)
    out = np.zeros(len(ids), dtype="f4")
    out[hit] = source_delta[at[hit]]
    return out, hit


def agreement(pred: np.ndarray, ref: np.ndarray) -> dict:
    if not len(pred):
        return {"count": 0}
    diff = pred.astype("f8") - ref.astype("f8")
    corr = float(np.corrcoef(pred, ref)[0, 1]) if len(pred) > 1 else None
    return {
        "count": int(len(pred)),
        "pearsonR": corr,
        "meanDeltaDifferenceMpc": float(diff.mean()),
        "medianAbsDeltaDifferenceMpc": float(np.median(np.abs(diff))),
        "p90AbsDeltaDifferenceMpc": float(np.percentile(np.abs(diff), 90)),
        "p99AbsDeltaDifferenceMpc": float(np.percentile(np.abs(diff), 99)),
        "referenceMedianAbsMpc": float(np.median(np.abs(ref))),
        "prototypeMedianAbsMpc": float(np.median(np.abs(pred))),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, default=Path(r"D:\desi\data_external\RSD"))
    ap.add_argument("--config", type=Path, default=ROOT / "tools/rsd/rsd_config.json")
    ap.add_argument("--assets", type=Path, default=ROOT / "assets")
    ap.add_argument("--index", type=Path, default=ROOT / "ply/reconstruction_dr1/private_index")
    ap.add_argument("--cache", type=Path, default=Path(r"D:\desi\data_external\RSD\derived\rsd-v1-exploratory\BGS_NGC_rsd_query_v1.npz"))
    ap.add_argument("--reference-report", type=Path, default=Path(r"D:\desi\data_external\RSD\derived\rsd-v1-exploratory\BGS_NGC_report.json"))
    ap.add_argument("--out", type=Path, default=ROOT / "tools/rsd/bgs-ngc-coverage-prototype.json")
    ap.add_argument("--viewer-sample", type=int, default=50000)
    ap.add_argument("--random-limit", type=int, default=1500000, help="deterministic random sample used for this pilot field")
    ap.add_argument("--mask-validation-sample", type=int, default=250000)
    ap.add_argument("--mask-percentile", type=float, default=99.0)
    ap.add_argument("--holdout-fraction", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--query-batch-size", type=int, default=50000)
    args = ap.parse_args()

    from astropy.cosmology import Planck18
    from astropy.io import fits
    from pyrecon import IterativeFFTReconstruction

    p = config(args.config)["tracers"]["BGS"]
    sources = local_sources(args.data_dir)
    data_path = _source_path(sources, "BGS", "NGC", "data")
    random_path = _source_path(sources, "BGS", "NGC", "random")
    ref_ids, ref_delta = cache_lookup(args.cache)
    old_report = json.loads(args.reference_report.read_text(encoding="utf-8"))
    if old_report.get("tracer") != "BGS" or old_report.get("region") != "NGC":
        raise ValueError("reference report is not for BGS NGC")

    t0 = time.time()
    data = _deduplicate_targetids(_fits_rows(data_path, Planck18, 0, p["zMin"], p["zMax"]))
    randoms = _deduplicate_targetids(_fits_rows(random_path, Planck18, args.random_limit, p["zMin"], p["zMax"]))
    if not len(data["position"]) or not len(randoms["position"]):
        raise ValueError("empty BGS NGC data or random selection")
    print(f"Loaded BGS NGC data={len(data['position']):,} random sample={len(randoms['position']):,} stride={randoms['stride']}", flush=True)
    random_rows_used = len(randoms["position"])
    random_stride = int(randoms["stride"])

    rng = np.random.default_rng(args.seed)
    n_hold = max(1, int(round(len(data["position"]) * args.holdout_fraction)))
    shuffled = rng.permutation(len(data["position"]))
    hold_idx = np.sort(shuffled[:n_hold])
    train_idx = np.sort(shuffled[n_hold:])

    # Calibrate support from the full local random catalog even though the pilot
    # reconstruction uses a bounded random subsample. The alternating rows are
    # disjoint and the held-out subset is not used to build the support tree.
    mask_randoms = _deduplicate_targetids(_fits_rows(random_path, Planck18, 0, p["zMin"], p["zMax"]))
    mask_train = np.asarray(mask_randoms["position"][::2], dtype="f4")
    mask_valid_all = mask_randoms["position"][1::2]
    n_mask_valid = min(args.mask_validation_sample, len(mask_valid_all))
    mask_valid = np.asarray(mask_valid_all[:n_mask_valid], dtype="f4")
    mask_random_count = int(len(mask_randoms["position"]))
    del mask_valid_all, mask_randoms
    mask_tree = cKDTree(mask_train, compact_nodes=True, balanced_tree=True)
    mask_dist, _ = mask_tree.query(mask_valid, k=1, workers=-1)
    support_radius = float(np.percentile(mask_dist, args.mask_percentile))
    random_validation = {
        "sourceRandoms": mask_random_count,
        "trainRandoms": int(len(mask_train)),
        "independentValidationRandoms": int(len(mask_valid)),
        "nearestDistancePercentilesMpc": np.percentile(mask_dist, [50, 90, 95, 99, 99.5, 100]).tolist(),
        "supportRadiusMpc": support_radius,
        "supportRadiusPercentile": args.mask_percentile,
        "validationAcceptedFraction": float(np.mean(mask_dist <= support_radius)),
        "interpretation": "empirical local random-cloud support in 3-D; finite random sampling makes this a diagnostic, not a formal DESI angular-mask polygon",
    }
    print(f"Random support calibrated: radius={support_radius:.3f} Mpc; valid sample={len(mask_valid):,}", flush=True)
    del mask_train, mask_valid, mask_dist

    cell_mpc = _physical_mpc(16.0, Planck18)
    smoothing_mpc = _physical_mpc(p["smoothingMpcH"], Planck18)
    train_pos = data["position"][train_idx]
    recon = IterativeFFTReconstruction(
        f=p["f"], bias=p["bias"], los=None,
        cellsize=cell_mpc, boxpad=1.2,
        positions=np.concatenate((train_pos, randoms["position"])),
    )
    recon.assign_data(train_pos, weights=data["weight"][train_idx])
    recon.assign_randoms(randoms["position"], weights=randoms["weight"])
    recon.set_density_contrast(smoothing_radius=smoothing_mpc)
    recon.run(niterations=3)
    print(f"BGS NGC reconstruction complete: nmesh={np.asarray(recon.nmesh).tolist()}", flush=True)
    del train_pos, randoms

    # Held-out LSS galaxies must have random-supported positions and are compared
    # to the existing full-data TARGETID cache. This measures data-removal sensitivity,
    # not independent truth accuracy.
    hold_pos = data["position"][hold_idx]
    hold_ids = data["targetid"][hold_idx]
    hold_dist, _ = mask_tree.query(hold_pos, k=1, workers=-1)
    box_half = np.asarray(recon.boxsize, dtype="f8") / 2
    box_center = np.asarray(recon.boxcenter, dtype="f8")
    in_box = np.all(np.abs(hold_pos - box_center) <= box_half, axis=1)
    hold_supported = (hold_dist <= support_radius) & in_box
    hold_query, _ = _query_radial_delta(recon, hold_pos[hold_supported], args.query_batch_size)
    hold_ref, hold_hit = lookup(hold_ids[hold_supported], ref_ids, ref_delta)
    holdout_cmp = agreement(hold_query[hold_hit], hold_ref[hold_hit])
    holdout_cmp.update({
        "heldOutRows": int(len(hold_idx)),
        "supportedByRandomsAndMesh": int(hold_supported.sum()),
        "supportFraction": float(hold_supported.mean()),
        "referenceCacheMatches": int(hold_hit.sum()),
        "referenceCacheMatchFractionOfSupported": float(hold_hit.mean()) if len(hold_hit) else 0.0,
        "meaning": "cross-fit sensitivity against the existing full-data field cache; this is not an independent velocity truth test",
    })

    # Draw an unbiased point sample from the current DR1 galaxy geometry inside
    # the BGS z interval; query only positions supported by randoms and mesh box.
    geometry = args.assets / "DR1_GALAXY.desi"
    h = read_desi_header(geometry)
    if h["version"] != 2:
        raise ValueError("expected current DR1_GALAXY v2 geometry")
    count = h["count"]
    zq = np.memmap(geometry, dtype="<u2", mode="r", offset=48 + count * 6, shape=(count,))
    z = h["zMin"] + np.asarray(zq, dtype="f8") * ((h["zMax"] - h["zMin"]) / 65535.0)
    eligible_idx = np.flatnonzero((z >= p["zMin"]) & (z < p["zMax"]))
    eligible_count = int(len(eligible_idx))
    n_candidate = min(args.viewer_sample, len(eligible_idx))
    sampled_idx = np.sort(rng.choice(eligible_idx, size=n_candidate, replace=False))
    del z, zq, eligible_idx
    pos_data = np.memmap(geometry, dtype="<i2", mode="r", offset=48, shape=(count, 3))
    query_pos = np.asarray(pos_data[sampled_idx], dtype="f8") / h["unitsPerMpc"]
    ids_path = args.index / "DR1_GALAXY_targetid.npy"
    viewer_ids = np.load(ids_path, mmap_mode="r")
    if len(viewer_ids) != count:
        raise ValueError("private current DR1 galaxy TARGETID index does not align with geometry rows")
    query_ids = np.asarray(viewer_ids[sampled_idx], dtype="u8")
    query_dist, _ = mask_tree.query(query_pos, k=1, workers=-1)
    in_box = np.all(np.abs(query_pos - box_center) <= box_half, axis=1)
    viewer_supported = (query_dist <= support_radius) & in_box
    print(f"DR1 sample selected={n_candidate:,} support={int(viewer_supported.sum()):,}; querying holdout and viewer positions", flush=True)
    viewer_delta, _ = _query_radial_delta(recon, query_pos[viewer_supported], args.query_batch_size)
    viewer_ref, viewer_hit = lookup(query_ids[viewer_supported], ref_ids, ref_delta)
    viewer_cmp = agreement(viewer_delta[viewer_hit], viewer_ref[viewer_hit])
    viewer_cmp.update({
        "catalog": "DR1_GALAXY.desi",
        "redshiftRange": [p["zMin"], p["zMax"]],
        "catalogRows": int(count),
        "catalogEligibleRows": eligible_count,
        "sampledEligibleRows": int(n_candidate),
        "randomAndMeshSupportedSampleRows": int(viewer_supported.sum()),
        "supportedFractionOfSample": float(viewer_supported.mean()) if len(viewer_supported) else 0.0,
        "supportedSampleQueries": int(viewer_delta.size),
        "referenceCacheMatches": int(viewer_hit.sum()),
        "referenceCacheMatchFractionOfSupported": float(viewer_hit.mean()) if len(viewer_hit) else 0.0,
        "sampleIsCatalogCoverageEstimate": False,
    })

    result = {
        "prototype": "BGS NGC; local-only; no payload or cache writes",
        "selection": {"source": str(data_path), "randoms": str(random_path), "redshift": [p["zMin"], p["zMax"]],
                      "weight": "WEIGHT * WEIGHT_FKP if present", "randomIndex": 0, "region": "NGC"},
        "reconstruction": {"algorithm": "pyrecon IterativeFFTReconstruction", "field": "rsd", "los": None,
            "dataRows": int(len(data["position"])), "trainingRows": int(len(train_idx)), "heldOutRows": int(len(hold_idx)),
            "randomRowsSourceFull": int(old_report["randoms"]["rowsKept"]),
            "randomRowsUsed": random_rows_used,
            "randomStride": random_stride, "cellSizeMpcH": 16.0, "cellSizeMpc": cell_mpc,
            "smoothingMpcH": p["smoothingMpcH"], "smoothingMpc": smoothing_mpc,
            "nmesh": np.asarray(recon.nmesh).tolist(), "boxsizeMpc": np.asarray(recon.boxsize).tolist(),
            "boxcenterMpc": box_center.tolist()},
        "randomSupport": random_validation,
        "holdoutVsExistingFullDataTargetidCache": holdout_cmp,
        "currentDr1ViewerSample": viewer_cmp,
        "limitations": [
            "A random-cloud proximity threshold is a finite-sampling support proxy; it does not replace the official angular mask or prove unbiased velocities.",
            "The existing full-data cache is a sensitivity reference, not ground truth. The holdout comparison quantifies the effect of omitting 20% of data.",
            "Viewer sample results are diagnostic only and are not extrapolated to a catalog-wide coverage claim.",
            "Only positions inside both the random support radius and reconstructed mesh box were queried; no correction asset was written.",
        ],
        "elapsedSeconds": time.time() - t0,
        "peakResidentMemoryGiB": float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
