# Rebuild the twelve current DESI V3 catalogue assets

This guide rebuilds the twelve files currently named `desiV3_catalog_000.bin.gz` through `desiV3_catalog_011.bin.gz`. They contain **15,786,217 observed DR1 objects** (14,140,375 galaxies and 1,645,842 quasars). The currently deployed files total **190,266,622 compressed bytes**. Their geometry is the observed redshift-space catalogue; **FoG and RSD deltas are not baked into these twelve files**. The disputed reconstruction is a separate, currently disabled product.

All commands below run from this repository's root. This is a Python data pipeline and does not require the visualization engine. The final copy step gives the files their deployed names. The local final manifest describes the same rows, byte sizes and encoding; a host that uses asset IDs may replace its file references after upload.

## Inputs, software and storage

| Input | Acquisition | Why it is needed |
|---|---|---|
| `data/zall-pix-iron.fits` | `download_sources.py --group dr1_zcatalog`; [official DR1 zcatalog](https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/zcatalog/v1/) | Exact DR1 GALAXY/QSO selection, redshift, sky position, TARGETID and targeting bits |
| Eight NGC/SGC `clustering.dat.fits` files and four matching `full_HPmapcut.dat.fits` files | `build_legacy_lss_catalogs.py`; [official LSS v1.5 directory](https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/) | Previously available imaging metadata joined by exact TARGETID |
| Missing-target Tractor columns | `desiV3/tools/fetch_missing_photometry_tap.py`; [DESI database access documentation](https://data.desi.lbl.gov/doc/access/database/) | Supported TYPE, SHAPE_R and G/R/Z fluxes for IDs absent from the legacy imaging stage |
| Optional independent source audit | `desiV3/tools/fetch_dr1_photometry.py --pixel 0`; [official DR1 imaging VAC](https://data.desi.lbl.gov/doc/releases/dr1/vac/lsdr9-photometry/) | Compare selected SQL values against an official checksum-verified FITS pixel |

The SQL stage queries public `desi_dr1.photometry` and checks each response against an independent `COUNT(*)` request. It retains raw CSV, exact SQL, source SHA256 and compact NPZ receipts. It is a large, multi-request operation; network access to the public service is necessary. Keep receipts because public data services and software environments can change. Allow substantial temporary storage: the zcatalog alone is about 22.4 GB, and the LSS FITS and intermediate arrays add many GB. **At least 100 GB of free working space and 64 GB RAM are practical planning targets**, not measured hard minima. Run the large Python stages sequentially.

### macOS

Install Python 3.11 or newer. For a Homebrew installation:

```bash
brew install python@3.12
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy astropy
python -m pip freeze > environment-catalog.txt
```

### Linux

On Debian/Ubuntu, install a Python 3.11+ interpreter and its venv package; use equivalent packages on other distributions:

```bash
sudo apt-get update
sudo apt-get install python3 python3-venv
python3 --version  # require 3.11 or newer
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy astropy
python -m pip freeze > environment-catalog.txt
```

`pyrecon`, MPI and the FoG/RSD dependencies in the top-level `requirements.txt` are **not needed** to rebuild these catalogue files. NumPy and Astropy are enough for the data pipeline. Keep `environment-catalog.txt` with the output for a reproducibility record.

## Stage 1 — reproduce the exact observed DR1 selection and point partition

Download the official FITS, build the original full selection, recover targeting flags, create the twelve historical point/type chunks, compact the DR1 catalogues to their v2 layout, then validate the private TARGETID index. The historical selection used `--min-reliable-z 0`; the current builder default of 0.0045 would change the 15.8-million-row selection.

```bash
python download_sources.py --group dr1_zcatalog --dest data
python build_dr1_full_catalog.py --source data/zall-pix-iron.fits \
  --out assets/v3 --min-reliable-z 0
python build_dr1_target_flags.py --source data/zall-pix-iron.fits \
  --dest assets/v3 --min-reliable-z 0
python build_desi_chunks.py --assets-dir assets/v3 \
  --out ply/chunks --chunks 12 --tag dr1v2
python compact_dr1_catalogs.py --source-dir assets/v3 --out-dir assets
cp assets/v3/DR1_GALAXY.target-flags.bin assets/
cp assets/v3/DR1_QSO.target-flags.bin assets/
python tools/rsd/build_dr1_chunk_index.py --source data/zall-pix-iron.fits \
  --assets assets --chunks ply/chunks \
  --out ply/reconstruction_dr1/private_index
```

The DR1 builder keeps the v3 files in `assets/v3`; compact v2 files go to `assets/`. The point builder **must** read the v3 files because they contain TARGETID. The indexer **must** read the compact v2 files and exact point chunks. Do not put either private TARGETID index on a public server. `index-report.json` must show twelve `geometryExact: true` rows and 15,786,217 total objects. Compare `ply/chunks/manifest.json` with [`reference-point-manifest.json`](reference-point-manifest.json): 12 chunk counts, tracer counts, type/target counts, partition name and `unitsPerMpc` must agree. Compressed byte sizes can vary with zlib version even when decoded geometry agrees.

The two compact v2 catalogues and two targeting files can also be checked against [`reference-inputs.json`](reference-inputs.json):

```bash
python desiV3/tools/verify_reference_inputs.py --assets assets
```

Exact SHA256 matches are strong evidence of an identical historical input; a different hash requires investigating selection, FITS version, quantization or row order before proceeding.

## Stage 2 — consolidate geometry and retain known imaging

First build the earlier four LSS v3 catalogues with exact TARGETID and imaging fields. The script downloads eight hemisphere clustering files and four `full_HPmapcut` files from the versioned DESI directory. Do **not** use `--no-morphology`; `enrich_morphology.py` requires v3 files with the four imaging byte arrays. If sources are already present in `data/`, `--skip-download` avoids network transfer, but it does not verify official checksums; record local hashes separately.

```bash
python build_legacy_lss_catalogs.py
python desiV3/tools/build_consolidated.py
python desiV3/tools/enrich_morphology.py
```

`build_consolidated.py` validates the old chunk's decoded position and type bytes against the newly consolidated data, verifies the exact TARGETID order, then writes `desiV3/data/catalog_000.bin.gz` … `catalog_011.bin.gz`, `manifest.json` and `validation.json`. `enrich_morphology.py` joins retained LSS metadata by exact uint64 TARGETID, excludes conflicting imaging tuples and writes a separate candidate to `desiV3/data/morphology/`. It asserts that bytes outside the four imaging arrays remain unchanged. The earlier morphology pass recovered imaging for 8,403,421 galaxies and 1,225,483 quasars; this is an intermediate state, **not** the current release.

The `DSC3` decompressed chunk has a 64-byte header followed by five aligned sections: one targeting/class byte per row, original uint32 row index, axis-separated native int16 position deltas, uint16 redshift and four uint8 imaging arrays (profile, angular radius, g−r, r−z). The header stores version, row count and section offsets. Each gzip member is generated with a fixed timestamp (`mtime=0`), but cross-platform zlib versions may still change compressed hashes. The underlying row/geometry/attribute validation is the primary scientific-data check.

## Stage 3 — complete imaging from the official photometry table

The SQL helper enumerates exact TARGETIDs whose intermediate profile byte is `255`, queries them in batches of 50,000 with four workers and resumes verified batches. It filters returned rows back to the requested IDs, checks uniqueness, and checks every broad interval response against an independent `COUNT(*)`. The final application refuses an incomplete progress manifest or changed compact hashes. It encodes only supported Tractor profiles: PSF, REX, EXP, DEV and SER. `TYPE` is a fitted imaging profile, **not** a validated spiral/elliptical classification.

```bash
python desiV3/tools/fetch_missing_photometry_tap.py
python desiV3/tools/apply_dr1_photometry.py
python desiV3/tools/finalize_manifest.py
python desiV3/tools/verify_final_catalog.py
python desiV3/tools/stage_current_assets.py
```

The output is `desiV3/data/full-photometry/catalog_000.bin.gz` … `catalog_011.bin.gz`, `validation.json`, `manifest.json` and a private `api-imaging.npz` containing full-precision source measurements for possible later API work. The staging command copies and verifies the twelve files under their **current asset filenames** in `desiV3/data/upload/`, with `desiV3-consolidated-manifest.json` and `staging-receipt.json`. Those twelve `desiV3_catalog_***.bin.gz` files are the deliverable. Upload only them and the matching manifest when deploying to a new visualization host; the private TARGETID index and source FITS stay offline.

The checked [active-manifest reference](reference-active-manifest.json) also records a separate nearby-position precision overlay for 43,922 objects. That overlay is outside the twelve catalogue files and is not generated by this pipeline. If a visualization host uses it, retain or rebuild that separate asset and its manifest entry when deploying; removing the entry changes nearby-position behavior even if all twelve catalogue hashes match.

The historical run queried 6,157,313 previously missing IDs in 124 verified batches and recovered 6,120,687 supported profiles. The current files provide supported profiles for 14,104,500 galaxies and 1,645,091 quasars; 35,875 galaxies and 751 quasars remain unknown. Missing or unsupported profiles are never silently classified. Geometry, redshifts, exact identity/order, targeting and spectral class remain unchanged through imaging enrichment.

### Optional independent FITS cross-check

The SQL workflow is sufficient to reproduce the current files. For an additional source check, download and verify one original official FITS pixel, then compare overlapping SQL values exactly:

```bash
python desiV3/tools/fetch_dr1_photometry.py --pixel 0
python desiV3/tools/verify_photometry_source.py
```

The historical comparison found 43,965 exact matches with zero profile or float32 differences. Downloading every official FITS pixel would be much larger and is unnecessary for the selected-column pipeline.

## Validation and interpretation

The final verifier checks twelve SHA256 hashes against the generated `validation.json`, the manifest counts and sizes, and the `DSC3` section layout. To demand **byte-for-byte agreement** with the currently deployed local originals, add `--require-reference-match`; the recorded hashes are in [`reference-current.sha256`](reference-current.sha256):

```bash
python desiV3/tools/verify_final_catalog.py --require-reference-match
```

An output can be scientifically identical while gzip hashes differ across zlib versions; if strict hash comparison fails, inspect the stage reports, decoded sections and original source receipts rather than replacing files blindly. The checked local current files had twelve matching SHA256 hashes and totalled 190,266,622 bytes. Matching these assets does not validate the separate RSD/FoG reconstruction.

The official query uses public infrastructure. If an endpoint, table schema or product path changes, record the replacement source and rerun source consistency checks. The fixed selection, exact TARGETID order, imaging-code mapping and row counts are required to claim reproduction of the current catalogue. The original four Python stages under `desiV3/tools/` were copied from the local production workflow; package changes and reference hashes are described in [PROVENANCE.md](../PROVENANCE.md).
