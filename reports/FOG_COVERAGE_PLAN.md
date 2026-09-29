# DESI DR1 FoG coverage audit and extension plan

Date: 2026-09-23. Scope: FoG coverage of the exact 15,786,217-TARGETID DR1 catalogue using Gfinder DR1 v1.0. This historical audit did not change any active correction payload.

## Measured coverage

The DR1 materializer matched **595,821** TARGETIDs to the calculated FoG table (3.77% of the catalogue); **195,200** had nonzero radial deltas (1.24%). The source analysis had 595,823 IDs, two absent from the current twelve-chunk DR1 catalogue. These counts measure the existing delta file's coverage, not all possible Gfinder matches. By class, 589,608 `DR1_GALAXY` rows (194,006 nonzero) and 6,213 `DR1_QSO` rows (1,194 nonzero) were matched. QSO interpretation needs separate review because Gfinder is a galaxy group catalogue.

For the 195,200 nonzero absolute deltas, median was 3.06 Mpc; p90 9.27, p95 11.84, p99 17.94, p99.9 29.09 and maximum 175.81 Mpc. Of 28,080 corrected groups, 13,645 (48.6%) reached the `alpha=0.08` floor. This is a transverse/radial dispersion ratio, not a 0.08 Mpc displacement limit. `maxDeltaMpc=300` excluded no generated value: 117 objects exceeded 35 Mpc, 52 exceeded 50 Mpc, 11 exceeded 100 Mpc, two exceeded 150 Mpc and none exceeded 200 Mpc.

The maximum belongs to TARGETID `39628423219380830`, `ZSPEC=0.21545`, Gfinder `IGRP=3129`, RICH=71 and group redshift 0.1679. The analysis retained 18 spectroscopic members, `sigma_parallel=7.07 Mpc`, `sigma_perp=0.377 Mpc`, elongation 18.76 and `alpha=0.08`. Its 175.81 Mpc delta corresponds to roughly 10,900 km/s under `v ≈ H(z)|delta_chi|/(1+z)`. Investigate membership, spectral quality and redshift errors; this is not a confirmed peculiar velocity.

## Why coverage remains 595,821

1. `match_gfinder.py` used the earlier 9,751,955-row LSS viewer index. It scanned 134,731,580 Gfinder rows and obtained 5,507,091 unique Gfinder TARGETIDs, corresponding to 5,549,432 earlier viewer rows (56.91%, including duplicates). Almost all also match the current catalogue: 5,388,525 galaxies plus 117,430 QSOs = 5,505,955 rows. Reindexing those same matched IDs against current DR1 would add only 1,136 matches.
2. `collect_members.py` excluded RICH<5 groups. Among 5,137,633 relevant groups, 3,197,763 have RICH=1 and cannot supply internal dispersion. The 352,351 groups with RICH>=5 yielded 3,593,784 collected members and 595,823 viewer-matched rows in `targetid-delta.npy`. `analyze_fog.py` additionally required five `ZSRC>0` members and radial/transverse elongation >=1.5. Then 28,080 groups and 195,200 rows received nonzero deltas. Eligible groups failing geometry retained zero deltas.

Sensitivity tests on the same members did not create new matches. With `min-spec=5`, lowering elongation 1.5→1.0 raised corrected groups 28,080→28,217 and nonzero rows 195,200→195,769. With `min-spec=3`, elongation 1.5 gave 65,117 groups and 283,000 nonzero rows; elongation 1.0 gave 65,907 and 285,051. All remained within the same 595,821 DR1-matched IDs. These were exploratory tests, not adopted parameters.

The local group table also had 1,017,771 RICH=2, 390,517 RICH=3 and 179,231 RICH=4 groups. Their members were not materialized under RICH>=5. Matching the full 5,505,955 current IDs is feasible, but RICH 1–4 does not automatically provide reliable dispersion. RICH>=3 merits a validation study; RICH=1/2 should not automatically be corrected.

## Abell 2162 check

At approximate catalogue coordinates RA=243.125°, Dec=+29.533°, z≈0.0310 ([reference](https://www.icc.dur.ac.uk/~tt/Lectures/Galaxies/LocalGroup/Back/nearsc.html)), nearby Gfinder `IGRP=4641737` has RICH=2, z=0.0320 and angular separation 0.045°, so fails `minimum-richness=5`. `IGRP=52780607` has RICH=1, z=0.0295 and separation 0.016°. Nearby RICH=5 `IGRP=964226` has z=0.0549 and only one ZSPEC member among five collected members, so fails `min-spec=5`. This explains local exclusions but does not establish physical membership in A2162. Check individual Gfinder memberships, ZSPEC, the mask and an external cluster reference before assigning corrections.

## Physical displacement checks

A universal 35 Mpc cap is not directly supported by the FoG references. For a hypothetical maximum peculiar speed, the first-order comoving radial scale is `|delta_chi|max(z) ≈ (1+z) v_max / H_Planck18(z)`. Under Planck18, 2,000–2,500 km/s corresponds to 29.7–37.1 Mpc at z=0.008, 31.0–38.7 Mpc at z=0.1 and 33.2–41.4 Mpc at z=1. A 35 Mpc threshold roughly resembles 2,300 km/s but is not redshift invariant. Such a hard individual-velocity limit needs evidence and was not adopted. Review the largest outlier with group membership, `ZERR`, spectra and dispersion.

Tegmark et al. compress FoG using group finding and radial versus transverse dispersion, with threshold sensitivity tests. The local analysis instead uses Gfinder membership and robust MAD. It is an exploratory visual approximation, not an exact reproduction of their procedure, and the reference does not authorize universal 35 or 300 Mpc caps.

## Extension and validation plan

1. Build a unique sorted TARGETID index for all 15,786,217 DR1 rows, retaining class and redshift. Rescan all 134,731,580 Gfinder rows, encoding IDs from `RELEASE`, `BRICKID`, `OBJID`; avoid angular or rounded-coordinate matching. Record unique hits, duplicates and tracer/redshift/ZSRC breakdown. Reproduce the earlier 5,507,091 Gfinder-ID / 5,549,432 viewer-row result as a regression check.
2. Join official `IGAL → IGRP` relations from `iDESIDR9.y1.v1_1.fits` and group properties from `DESIDR9.y1.v1_group.fits`. Test RICH>=3 separately from RICH>=5. Keep RICH=1/2 and missing ZSPEC apart. Use only `ZSRC>0` to estimate radial dispersion; the official Gfinder documentation identifies `ZSRC=0` as photometric redshift.
3. With Planck18 distances, compare scaled MAD to robust/truncated standard deviation; inspect `ZERR`, spectral quality and interlopers. Publish a sensitivity matrix for RICH>=3/5, `min-spec`=3/5 and elongation=1.0/1.5/2.0. Test group stability under leave-one-spectrum-out and compare observed/residual radial dispersion. Keep insufficient or nonpositive-dispersion groups at zero delta.
4. Report displacement percentiles and counts above 35/50/100 Mpc. If testing a speed cap, convert it by redshift using Planck18 and compare 2,000 and 2,500 km/s with mocks and spectral distributions before adoption. The current 300 Mpc cap has no effect on existing values.
5. Require every nonzero delta to trace to verifiable `IGAL`, `IGRP` relation and group properties. Report separate denominators for current TARGETID, Gfinder match, group member, spectroscopic member, geometry-qualified group and nonzero correction, by class and redshift. Preserve exact TARGETID/chunk order, finite values and zero for ineligible rows. Recheck a known rich cluster, including failed cuts around Abell 2162. Accept a profile only after membership, dispersion, resampling and independent comparison pass.

## Sources and limits

- [Official Gfinder DR1 documentation](https://data.desi.lbl.gov/doc/releases/dr1/vac/gfinder/) defines files and `IGAL`, `IGRP`, `RICH`, `GRP_Z`, `ZSRC`, including `ZSRC=0` photometric and `ZSRC>0` spectroscopic. Gfinder's selection is limited to `MAG_Z<21`.
- [Yang et al., extended halo-based group method](https://arxiv.org/abs/2012.14998) studies mocks with z<=1 and MAG_Z<=21, including groups with at least three members. This motivates a sensitivity test, not direct DR1 calibration.
- [Tegmark et al. 2004, section 3.3 and figure 7](https://web.physics.rutgers.edu/grad/690/Tegmark-etal-2004.pdf) supports radial/transverse compression and threshold testing, not this implementation's exact cuts.

Historical outputs included `coverage-report.json`, Gfinder `matches/`, `groups/`, `members/`, `analysis/summary.json`, `targetid-delta.npy` and `group-stats.npy`. The three sensitivity runs used `analyze_fog.py` with `(min-spec, elongation)` = `(5, 1.0)`, `(3, 1.5)` and `(3, 1.0)` in separate temporary output directories; they did not replace the published deltas.
