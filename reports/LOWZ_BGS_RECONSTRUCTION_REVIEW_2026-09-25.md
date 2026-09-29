# Low-redshift BGS field: technical review (2026-09-25)

## Why the z=0.10 gap appeared

The active exploratory RSD assets use the DESI DR1 LSS v1.5 `BGS_BRIGHT-21.5` clustering sample for BGS. The local NGC and SGC data files each contain zero rows below z=0.10. The materializer therefore sets RSD displacement to zero below the configured BGS interval `[0.10, 0.40)`, while objects immediately above 0.10 can move when the slider is used. The resulting radial density deficit is a selection-boundary artifact.

## Replacement source and download

Official DESI DR1 `iron/LSScats/v1.5` **BGS_BRIGHT** clustering data and matching index-0 randoms were downloaded to `tools/rsd/lowz_sources`. Index-1 BGS_BRIGHT randoms were also downloaded for an independent-realization support diagnostic. Each download is recorded by source URL, expected byte count and local SHA-256 in its `.source.json` sidecar. `download_lowz_sources.py` implements checked HTTP ranges and is resumable. The NGC data file was first downloaded through a resumable direct request and then locally hashed; its sidecar uses the same format.

`audit_lowz_bgs_inputs.py` verifies FITS headers/columns, finite positive weights, redshift ranges and exact `TARGETID` overlap with the current geometry index. Index-0 source counts:

| Region | Data rows, 0.01 ≤ z < 0.10 | Random rows, 0.01 ≤ z < 0.10 |
| --- | ---: | ---: |
| NGC | 430,166 | 1,910,677 |
| SGC | 148,435 | 815,998 |

The current viewer has 882,086 objects in `[0.01, 0.10)`, of which 578,458 match these low-z data by exact `TARGETID`. The remaining viewer objects require a position query inside verified field support; an ID match is not invented.

## Candidate field and technical checks

`build_lowz_bgs_field.py` reconstructs one continuous BGS_BRIGHT field for `[0.01, 0.40)` independently in NGC and SGC using `IterativeFFTReconstruction`, local line of sight, `field="rsd"`, 16 h⁻¹ Mpc nominal cells, 15 h⁻¹ Mpc smoothing, boxpad 1.2, three iterations, b=1.5 and f=0.682. It uses all selected data rows and every sixth selected index-0 random row for mesh construction. The full index-0 randoms determine a p99 nearest-neighbor support proxy; candidate points must have unique NGC or SGC support and lie inside the corresponding FFT box. These values are provisional: a flux-limited BGS_BRIGHT sample need not share the optimal bias or growth approximation of BGS_BRIGHT-21.5. Random-stride convergence and an official survey mask check remain necessary.

The candidate produced 4,714,019 NGC and 1,612,246 SGC viewer deltas, 6,326,265 in total for 6,643,924 viewer positions in `[0.01, 0.40)`. At 49,473 NGC and 49,625 SGC sample IDs present in the source catalogues, the median absolute difference between querying the source coordinate and the quantized viewer coordinate was 0.00376 and 0.00369 Mpc respectively. This tests numerical/geometry consistency, not physical truth.

`validate_lowz_bgs_field.py` compared the candidate against the currently active RSD payloads. On 5,343,619 shared corrected objects above z=0.10, the displacement correlation is 0.8501 and the median absolute difference is 1.363 Mpc. Replacing the BGS field is therefore a material scientific change, even though it removes the artificial z=0.10 field switch. In a 5 Mpc radial bin immediately after the old boundary, the observed count is 19,024, the active reconstruction yields 13,912 and the continuous candidate yields 18,834. The adjacent bin immediately before the boundary changes from 22,598 active to 19,216 candidate, compared with 19,329 observed. These are descriptive point counts, not a cosmological goodness-of-fit test.

Local review-only payloads and exact-ID LOD pairs were built under `ply/reconstruction_dr1/staged_lowz_bgs_v1`. The 12 RSD point arrays and two LOD pairs pass Float32 finiteness, count, gzip round-trip, and exact point/LOD join checks; the FoG layer is unchanged. Their coverage report records 15,786,217 total viewer objects, 15,268,054 eligible within the configured fields including low-z BGS, and 14,345,346 matched/nonzero RSD deltas. `asset-staging-manifest.json` lists hashes and retains `readyForUpload: false` as a scientific gate. This flag prohibits interpreting the files as scientifically approved; the user authorized an exploratory PlayCanvas upload for visual review.

## PlayCanvas review deployment and duplicate cleanup

The active scene is project 1605307, scene 2602187, branch `main`. `DESI_Points.rsdOverlayUrl` points to overlay asset 308591082, which references twelve new RSD chunk assets 308591055, 308591056, 308591057, 308591058, 308591059, 308591060, 308591061, 308591062, 308591064, 308591066, 308591069 and 308591071. `GalaxyLOD.reconstructionDeltaAssets` points to the two reuploaded pairs 308622658 (`correction_DR1_GALAXY.bin`) and 308622659 (`correction_DR1_QSO.bin`). Their remote sizes and MD5 hashes match the local staged files: 57,825,350 bytes / `65511b5dc7a89fe80e5f60fc739bf0e9` and 5,129,512 bytes / `ab1d9c31765158b23d46108a6e8e5aa8`.

The earlier LOD IDs 308589367, 308589368, 308591072 and 308591073 were checked by asset ID and are absent from PlayCanvas after duplicate cleanup. The first two are retained in the local backup at `ply/reconstruction_dr1/backup_before_lowz_2026_09_25`; the latter two are reproduced by the staged low-z pair files. The old RSD chunks and overlay were also removed after their verified local backup. The twelve existing FoG assets remain in use.

A fresh WebGL2 Launch loaded all 15,786,217 points. With the slider at 100%, runtime logs reported 14,140,375 GALAXY paired deltas and 1,645,842 QSO paired deltas ready. The slider returned to 0% observed mode. This verifies asset binding and application loading, not the official mask or physical correctness of the reconstructed field. The numerical radial occupancy diagnostic above supports removal of the z=0.10 artifact; astronomer visual review and independent mask validation remain pending.

## Scientific review still needed

1. Reconstruct the tracer-specific official angular and radial mask for BGS_BRIGHT v1.5, including veto and completeness rules. Classify every candidate ID, especially those without an exact BGS_BRIGHT data-row match, against it. Random-cloud proximity is not equivalent to official mask membership.
2. Test sensitivity to the every-sixth-random reconstruction, `b`, `f`, smoothing and mesh size. Validate near z=0.01, z=0.10 and z=0.40 with mocks or an independent reference. The current exact-ID comparison only checks coordinate consistency.
3. Check the physical appropriateness of replacing the BGS_BRIGHT-21.5 field with BGS_BRIGHT across all `[0.01,0.40)` and whether the chosen FoG layer should be added independently.
4. If deployed for visual review, verify slider 0 → 100 → 0, both LOD pairs, point/LOD agreement and the disappearance of the reported z=0.10 void. Do not describe a visual-review deployment as scientific approval.

## Reproduction

From `D:\desiV2`, using the installed WSL DESI environment and official source files already in `tools/rsd/lowz_sources`:

```bash
cd /mnt/d/desiV2
/home/khyro/.venvs/desi-rsd/bin/python tools/rsd/audit_lowz_bgs_inputs.py
/home/khyro/.venvs/desi-rsd/bin/python tools/rsd/build_lowz_bgs_field.py --random-stride 6
/home/khyro/.venvs/desi-rsd/bin/python tools/rsd/validate_lowz_bgs_field.py
/home/khyro/.venvs/desi-rsd/bin/python tools/rsd/audit_lowz_random_realizations.py
```

Materialization uses `materialize_dr1_corrections.py --rsd-extensions tools/rsd/extensions --bgs-lowz-candidate tools/rsd/lowz_sources/field_candidate --stage-unvalidated-rsd --out <new-empty-dir>`. The candidate flag is restricted to review mode. `materialize_dr1_lod_pairs.py` then rebuilds the two LOD pair files against the same payload directory.

Official references: https://data.desi.lbl.gov/doc/releases/dr1/ and https://desidatamodel.readthedocs.io/en/latest/DESI_ROOT/survey/catalogs/RELEASE/LSS/SPECPROD/LSScats/VERSION/index.html .
