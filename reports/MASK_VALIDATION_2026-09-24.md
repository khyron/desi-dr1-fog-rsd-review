# DESI DR1 v1.5 RSD mask audit — 2026-09-24

## Status and inputs

The historical extension payloads were staged for visual inspection without scientific approval. `asset-staging-manifest.json` retained `readyForUpload: false` as a scientific gate. The current viewer release disables corrections. `materialize_dr1_corrections.py --stage-unvalidated-rsd` permits local review output; approved output requires a `--mask-approval-report` bound by hash to the exact extension manifest.

The reconstruction used DESI DR1 LSS `iron/LSScats/v1.5`: five tracers in NGC/SGC, matching index-0 clustering randoms, half-open redshift intervals from `rsd_config.json`, and finite positive `WEIGHT × WEIGHT_FKP`. Ten data and ten index-0 random FITS files were available. Ten index-1 random files were subsequently downloaded (**11,208,458,880 bytes**) for an independent random-realization comparison. `source-manifest-index1.json` recorded URLs, sizes, SHA256, row counts and `DESIDR=dr1`; FITS columns were checked. Index-1 randoms do not provide an independent official survey mask.

The [DESI DR1 documentation](https://data.desi.lbl.gov/doc/releases/dr1/) and [LSS data model](https://desidatamodel.readthedocs.io/en/latest/DESI_ROOT/survey/catalogs/RELEASE/LSS/SPECPROD/LSScats/VERSION/index.html) describe the clustering products and angular veto variants. `full`, `full_noveto`, `full_HPmapcut` and HEALPix property maps require a scientifically specified selection before use. Property maps alone do not define a tracer-, hemisphere- and redshift-specific acceptance mask.

## Audit performed

1. The extension manifest recorded ten complete fields and output hashes. `audit_rsd_extensions.py` returned `complete`, checking field reproducibility and integrity rather than scientific correctness.
2. `audit_rsd_support_by_field.py` used complete index-0 randoms, an even/odd split, p99 nearest-neighbor distances and FFT box limits. This three-dimensional support test is a heuristic, not the DESI footprint boundary.
3. `audit_extension_angular_mask.py` sampled 10,000 TARGETIDs per field with seed 20260924. Positions came from exact chunk geometry and index. Up to 500,000 valid randoms per field were split into training/holdout sets; first- and eighth-neighbor angular distances were compared with holdout p99 thresholds. Eighth-neighbor distance detects low angular density but cannot uniquely identify a veto hole.
4. Of 100,000 sampled positions, **5,798** exceeded at least one index-0 p99 threshold. Flags per 10,000: BGS NGC 445, SGC 505; LRG NGC 545, SGC 1,058; LRG+ELG NGC 334, SGC 584; ELG NGC 390, SGC 541; QSO NGC 564, SGC 832. These are sample rates, not exact totals for full extensions.
5. The same sample with index-1 randoms flagged **5,738/100,000**. Exact-TARGETID comparison found **3,447** flagged by both realizations and **8,089** by at least one. LRG SGC had 761 persistent flags per 10,000. Persistent flags deserve review but are not automatically outside the official mask.

The historical generated reports are `extensions/angular-mask-audit.json`, `extensions/angular-mask-audit-index1.json` and `extensions/angular-random-realization-comparison.json`. They are not included in this source-only repository.

## Reproduction

After downloading the exact LSS inputs and installing the Python dependencies, run from the repository root:

```bash
python tools/rsd/verify_validation_sources.py
python tools/rsd/audit_rsd_extensions.py
python tools/rsd/audit_extension_angular_mask.py
python tools/rsd/audit_extension_angular_mask.py --random-index 1
python tools/rsd/compare_angular_random_realizations.py
```

The historical Windows angular audit completed with NumPy 2.4.6, SciPy 1.17.1 and Astropy 8.0.1. WSL was unavailable in that session (`E_ACCESSDENIED`), so the global audit was not rerun there. The historical downloader checked HTTP length and ranges, reused complete files and assembled partial transfers; the verifier checked FITS columns, rows and SHA256. Private `.npz` flag files contain TARGETIDs and should stay offline.

## Scientific review still required

1. Specify the official angular and radial selection for all ten clustering fields, including vetoes, `HPmapcut`, completeness and tracer/hemisphere dependence. Record versions, code, parameters, hashes and provenance. Download additional `full` variants only after selecting the precise products.
2. Classify flagged cases and unflagged controls against the official boundary: outer edge, internal hole/veto, sparse random sampling, quantized position, or accepted region. Measure distance to the actual mask edge.
3. Apply the official mask to all **8,298,917** extended TARGETIDs, reporting acceptance by field, hemisphere, redshift and edge distance. Compare displacements with an independent reference in spatial holdouts and redshift bins; test random stride, weights and mesh sensitivity.
4. If approved, produce JSON with `status: "approved"`, exact `extensionManifestSha256`, reviewer, date, mask version and per-field results. Rematerialize into an empty output directory using `--rsd-extensions` and `--mask-approval-report`; compare decompressed hashes and counts. The assigned-tracer precedence retains extension deltas for **4,166** cross-hits against foreign-tracer caches.

Both random realizations share the same selection construction and vetoes. Sampling 500,000 randoms per field broadens neighbor distances compared with full randoms. Neither persistent flags nor unflagged cases certify validity. Quantized `.desi` coordinates and cache agreement do not establish true peculiar velocities or authorize scientific publication.
