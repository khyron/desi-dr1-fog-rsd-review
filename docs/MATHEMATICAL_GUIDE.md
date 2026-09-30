# Mathematical methodology of the FoG and RSD review candidate

This document specifies the **exploratory correction candidate stored in this repository** for scientific review. It distinguishes measured inputs, modelling choices, code-level operations and unresolved validation. The resulting coordinates are hypotheses about real-space positions, not measured peculiar velocities or a DESI-approved reconstructed catalogue. The corrections are disabled in the current visualization release. Reproduction commands and source versions are in the [README](../README.md); the [worked regions](../examples/README.md) show diagnostics from the bundled candidate.

## 1. Observed coordinates and conventions

DESI measures an object's right ascension $\alpha$, declination $\delta$ and spectroscopic redshift $z$. The observed redshift combines cosmological expansion and line-of-sight peculiar motion. A redshift-derived radial distance therefore places an object in **redshift space**. We denote its observed catalogue position by $\mathbf{s}$ and its unknown real-space position by $\mathbf{r}$. The procedure constructs a candidate radial displacement of $\mathbf{s}$; it does not observe $\mathbf{r}$.

To first order in peculiar speed, $(1+z_{\rm obs})\simeq(1+z_{\rm cos})(1+v_{\parallel}/c)$. This relation explains why a velocity component parallel to the line of sight changes an inferred radial distance while leaving the measured sky angles essentially unchanged. The pipeline does **not** invert this expression object by object; FoG and RSD estimate displacements from group geometry and a reconstructed large-scale field, respectively.

With the Planck18 cosmological model, the code converts redshift to a comoving radial distance in physical megaparsecs (Mpc):

$$
\chi(z)=c\int_0^z\frac{dz'}{H(z')}, \qquad
\mathbf{s}=\chi(z)
\begin{pmatrix}
\cos\delta\cos\alpha\\
\cos\delta\sin\alpha\\
\sin\delta
\end{pmatrix}.
$$

Here $c$ is the speed of light and $H(z)$ is the expansion rate assumed by Planck18. The code uses `Planck18.comoving_distance(z).value` in **physical comoving Mpc**. This is a coordinate conversion, not a correction. Peculiar motions are particularly relevant to redshift-derived distances nearby. The catalogue contains 15,786,217 selected DR1 objects, rather than a complete census of the underlying galaxy distribution.

For any nonzero position, define the direction from Earth to that object as $\hat{\mathbf{n}}=\mathbf{s}/|\mathbf{s}|$. Both correction files store a **scalar radial displacement** $\Delta$ in physical Mpc. Applying one moves the object along that direction:

$$
\mathbf{s}_{\rm corrected}=\mathbf{s}+\Delta\hat{\mathbf{n}}.
$$

A positive $\Delta$ moves it farther away; a negative $\Delta$ brings it closer. The displayed examples multiply Mpc by approximately $3.26156$ to label distances in million light-years (Mly). A value quoted in $h^{-1}\mathrm{Mpc}$ is converted before geometry calculations: $d_{\rm Mpc}=d_{h^{-1}\mathrm{Mpc}}/h$, where $h$ comes from Planck18.

## 2. FoG: group-based radial compression

**Fingers of God** are galaxy groups that appear extended along the line of sight because their members have different internal velocities. The implementation uses DESI DR1 Gfinder v1.0 galaxy, group and membership relations. It joins those relations to the selected DESI objects by exact identifiers, ultimately `TARGETID`; it never assigns a correction by nearest sky position. Only spectroscopic members (`zsrc > 0` in the imported Gfinder table) enter the group-spread estimates, though the resulting scale factor is evaluated on the group's matched members.

For galaxy $i$ in a group, the observed radial separation from the Gfinder group centre is

$$
d_{\parallel,i}=\chi(z_i)-\chi(z_g),
$$

where $z_g$ is the group's redshift. To estimate the group's transverse scale, the code calculates angular separation $\theta_i$ between each spectroscopic member and the group centre, then uses $\chi(z_g)\theta_i$. The angle is in radians. Each distribution gets a robust spread estimate

$$
R(x)=1.4826\cdot\mathrm{median}\left(\left|x-\mathrm{median}(x)\right|\right),
\quad
\sigma_{\parallel}=R(d_{\parallel}),
\quad
\sigma_{\perp}=R\left(\chi(z_g)\theta\right).
$$

The angular separation is calculated from

$$
\cos\theta_i=\sin\delta_i\sin\delta_g+
\cos\delta_i\cos\delta_g\cos(\alpha_i-\alpha_g).
$$

The factor $1.4826$ makes the median absolute deviation comparable to a standard deviation for a Gaussian distribution; no Gaussian group model is imposed. The code requires at least **five spectroscopic members**, finite positive spreads, and an elongation ratio $E=\sigma_{\parallel}/\sigma_{\perp}$ of at least **1.5**. If those tests fail, that group receives zero FoG displacement. The transverse scale $\chi(z_g)\theta_i$ is a projected separation around the Gfinder centre, rather than a three-dimensional halo fit.

For an eligible group, the radial scale factor is

$$
a=\max\left(0.08,\min\left(1,\frac{\sigma_{\perp}}{\sigma_{\parallel}}\right)\right),
\qquad
\Delta_{{\rm FoG},i}=\mathrm{clip}\left[(a-1)d_{\parallel,i},-300,+300\right]\mathrm{Mpc}.
$$

The factor $a$ shrinks each radial offset toward the group centre; the $0.08$ floor prevents complete collapse and the $300$ Mpc cap limits an individual move. Those numerical values are **exploratory model choices**, not DESI-approved thresholds. The calculation estimates group *shape*, not an individual galaxy's measured velocity.

**Small example:** if a member is $20$ Mpc behind its group centre, and the group's transverse/radial spread ratio gives $a=0.5$, then $\Delta_{\rm FoG}=(0.5-1)20=-10$ Mpc. Its radial position moves $10$ Mpc toward Earth. A member $20$ Mpc in front gets $+10$ Mpc, also moving toward the group centre. The cap does not matter in this example.

## 3. RSD: density-field reconstruction and radial projection

Large-scale redshift-space distortion is associated with coherent motion toward or away from overdense regions. The reconstruction uses the DESI DR1 LSS galaxy and quasar *data* catalogues plus matching *random* catalogues. Random points trace where the survey could have observed targets; comparing data with randoms helps separate clustering from the survey's angular and radial selection. The catalogue weights are `WEIGHT`, multiplied by `WEIGHT_FKP` when that column exists. Data are split by tracer, redshift interval and north/south Galactic-cap region.

A simplified picture of the field supplied to the reconstruction is the weighted overdensity

$$
\delta_g(\mathbf{x})\sim
\frac{n_{\rm data}(\mathbf{x})-A\cdot n_{\rm random}(\mathbf{x})}
{A\cdot n_{\rm random}(\mathbf{x})},
$$

where $A$ normalizes the total weighted random counts to the data counts. This is **an intuition for the data-minus-randoms step**, not the exact `pyrecon` internal estimator. The actual script calls `pyrecon.IterativeFFTReconstruction`, assigns weighted data and random positions, smooths the density field, and runs three iterations for the recorded candidate. It passes tracer bias $b$, growth rate $f$, a local line of sight (`los=None`), and a mesh cell size to `pyrecon`. The adopted $b$, $f$, smoothing and grid values are listed in [`rsd_config.json`](../tools/rsd/rsd_config.json).

The installed `pyrecon` API returns a shift vector $\mathbf{q}_{\rm RSD}(\mathbf{s})$ from `read_shifts(position, field="rsd")`. In the convention used by this script, reconstructed position is $\mathbf{s}-\mathbf{q}$. The exported correction is therefore **the negative radial projection**:

$$
\Delta_{\rm RSD}(\mathbf{s})
=-\mathbf{q}_{\rm RSD}(\mathbf{s})\cdot\hat{\mathbf{n}},
\qquad
\mathbf{s}_{\rm RSD}
=\mathbf{s}+\Delta_{\rm RSD}\hat{\mathbf{n}}.
$$

Only the radial component goes into the viewer payload; the transverse part of `pyrecon`'s vector is not applied to the observed sky direction. This sign convention is important: reversing it would move objects in the opposite direction. A positive stored value moves an object outward, even though it is obtained from a *negative* dot product with the returned shift.

The reference field profile and low-redshift candidate used the following parameters. Intervals are half-open, $[z_{\min},z_{\max})$; NGC and SGC are reconstructed separately.

| Field | Redshift range | Bias $b$ | Growth $f$ | Smoothing ($h^{-1}$ Mpc) |
|---|---:|---:|---:|---:|
| BGS_BRIGHT-21.5 base BGS | $[0.10,0.40)$ | 1.5 | 0.682 | 15 |
| LRG | $[0.40,0.80)$ | 2.0 | 0.834 | 15 |
| LRG+ELG | $[0.80,1.10)$ | 1.6 | 0.870 | 15 |
| ELG | $[1.10,1.60)$ | 1.2 | 0.900 | 15 |
| QSO | $[1.60,2.10)$ | 2.1 | 0.928 | 30 |
| BGS_BRIGHT low-redshift replacement | $[0.01,0.40)$ | 1.5 | 0.682 | 15 |

The recorded runs used **16 $h^{-1}$ Mpc mesh cells** and three iterations because available RAM limited the analysis. The low-redshift run used `boxpad=1.2` and every sixth reconstruction random; the *full* selected random set was still used for the support heuristic below. A finer grid or all reconstruction randoms creates a **different** candidate and requires separate validation. The low-redshift BGS_BRIGHT field replaces the base BGS field throughout $0.01\le z<0.40$, including both sides of $z=0.10$, to avoid a hard field switch at that boundary. Random realization 1 was reserved for diagnostics; it does not certify the selection mask.

### 3.1 Query support for objects outside the LSS input

The core LSS data provide exact-ID shift queries for their own targets. The extension stage also queries the reconstructed mesh at positions of selected DR1 objects absent from the original clustering input. For the low-redshift BGS extension, a viewer position is queried only if its redshift is in $[0.01,0.40)$, it lies inside the reconstruction mesh, and it has unique support in either NGC or SGC. The support test takes the full selected random catalogue, builds a 3D nearest-neighbour tree from one half, and sets a radius to the **99th percentile** nearest-neighbour distance of the held-out half. A candidate passes if its distance to the training randoms is within that radius. Objects supported by both NGC and SGC are excluded as ambiguous. This is a **proximity heuristic**, not the official DESI angular/radial footprint or veto mask. Its behaviour near survey edges and holes remains an open scientific question.

## 4. Exact-ID materialization and combined candidate

Each correction source is keyed by the exact unsigned 64-bit `TARGETID`. The materializer validates the mapping from point order to original DR1 catalogue rows, joins FoG deltas, selects one RSD tracer field by the half-open redshift intervals above, and fills IDs absent from the base query cache with same-tracer extension values. Foreign-tracer cache hits do not take precedence. When supplied, the low-redshift BGS candidate **replaces** the earlier BGS values throughout $[0.01,0.40)$. An unmatched source contributes zero; zero does not establish that the true correction is zero. The two stored scalars can be inspected separately or summed for a *visual comparison*:

$$
\Delta_{\rm both}=\Delta_{\rm FoG}+\Delta_{\rm RSD},
\qquad
\mathbf{s}_{\rm both}=\mathbf{s}+\Delta_{\rm both}\hat{\mathbf{n}}.
$$

This is an additive review convention. It does **not** prove the two effects are independent or that their sum is an unbiased real-space estimate. The materializer verifies finite values and rejects $|\Delta_{\rm FoG}|>300$ Mpc or $|\Delta_{\rm RSD}|>150$ Mpc. The twelve V3 catalogue files retain observed positions. Separate gzip-compressed float32 arrays hold one Mpc delta per V3 row, with exact chunk and row correspondence. The included identity arrays map each V3 class/original-row pair back to `TARGETID` without rebuilding the catalogue. `verify_review_candidate.py` checks hashes, lengths and finite contents.

The recorded materialization report for the bundled candidate gives the following coverage. "Matched" means an exact ID joined to a source value; "nonzero" means a stored displacement other than zero. Neither means that the displacement is accurate.

| Measure | Rows |
|---|---:|
| Catalogue objects | 15,786,217 |
| FoG exact-ID matches | 595,821 |
| Nonzero FoG deltas | 195,200 |
| RSD interval-eligible objects | 15,268,054 |
| RSD exact-ID matches | 14,345,346 |
| Nonzero RSD deltas | 14,345,346 |

The local preview chooses an observed-space box around a requested `TARGETID`. Width and height are perpendicular to the Earth-to-target line of sight; depth is parallel to it. It draws the **same selected IDs** before and after correction, so a point may leave the plotted box after moving. Colors identify inward RSD, outward RSD, FoG-only, and unchanged objects. Histograms show signed RSD and absolute displacement, report exact zeros separately, and retain numeric bin counts in JSON. This makes the candidate inspectable without the visualization engine:

```bash
python verify_review_candidate.py
python visualize_target_region.py 39627715854207566 \
  --width-mly 500 --height-mly 500 --depth-mly 500 \
  --component both --out local-review/midrange
```

Open `local-review/midrange/comparison.svg` and `histograms.svg`, then inspect `summary.json`, `histograms.json`, and `largest_displacements.csv`. Repeat with `--component rsd` and `--component fog` for the same box. See the [four included comparisons](../examples/README.md) for ready-made plots.

## 5. Interpretation, diagnostics and scientific acceptance

File hashes, finite numbers, exact IDs and plausible-looking plots establish **computational integrity**, not scientific validity. An independent review should examine at least:

1. **Survey selection:** Do data and randoms represent the same angular and radial mask in each tracer/redshift region, including boundaries and the low-redshift BGS extension?
2. **Parameter sensitivity:** How do results change with finer mesh cells, smoothing, bias, growth rate, random sampling and FoG thresholds? A stable answer should not depend strongly on a resource-saving setting.
3. **External truth tests:** On simulations or mocks where real-space positions are known, do corrected positions actually move closer to truth? Compare clustering statistics and group shapes before and after correction.
4. **Outliers and boundaries:** Are the largest FoG/RSD shifts tied to reliable group membership and supported field regions? Are there discontinuities near $z=0.01$, $0.10$, $0.40$ and other field boundaries?
5. **Sign and frame:** Does an independent implementation reproduce the sign of the `pyrecon` shift, the Mpc versus $h^{-1}$ Mpc conversion, and the Earth-centred line of sight?

Until those tests pass, the files should be described as an **unvalidated reconstruction candidate**. The plots are designed to help a reviewer find problems, including a result that simply *looks* wrong.

## 6. Implementation trace for independent audit

| Operation described here | Implementation and audit record |
|---|---|
| Planck18 Cartesian conversion, weighted LSS inputs, `pyrecon` configuration and negative radial projection | [`tools/rsd/rsd_pipeline.py`](../tools/rsd/rsd_pipeline.py), [`tools/rsd/rsd_config.json`](../tools/rsd/rsd_config.json) |
| Exact Gfinder membership joins, robust group spreads and clipped radial FoG delta | [`tools/fog/analyze_fog.py`](../tools/fog/analyze_fog.py) and preceding scripts in [`tools/fog/`](../tools/fog/) |
| Low-$z$ BGS data/random selection, support proxy, mesh query and cross-check against shared source IDs | [`tools/rsd/build_lowz_bgs_field.py`](../tools/rsd/build_lowz_bgs_field.py), [mask audit](../reports/MASK_VALIDATION_2026-09-24.md) |
| Redshift field assignment, low-redshift replacement, exact-ID joins, zero-unmatched policy and payload export | [`tools/rsd/materialize_dr1_corrections.py`](../tools/rsd/materialize_dr1_corrections.py) |
| Current candidate integrity and observed-versus-corrected comparisons | [`verify_review_candidate.py`](../verify_review_candidate.py), [`visualize_target_region.py`](../visualize_target_region.py), [examples](../examples/README.md) |

The repository's [`PROVENANCE.md`](../PROVENANCE.md) distinguishes preserved historical scripts from adaptations made for reproducibility. An independent recalculation should record the exact DESI source-file versions and hashes, software versions, parameter file, field reports, identity/index validation and materialization report. A matching plot alone is weaker evidence than matching the underlying numerical arrays and source receipts.
