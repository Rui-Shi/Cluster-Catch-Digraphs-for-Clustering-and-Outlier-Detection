# WP8 experiment 4 — false positives under violated CSR, no outliers (R3.4, R3.5)

**Script:** `revision_experiments/tr1/88_wp8_csr_violation.R`. **Grid:** launched
2026-09-05 in the background (RK methods U-MCCD/SU-MCCD at 40 reps, all other
methods at 100 reps, one process), complete as of the files analysed here
(2026-09-06). Analysis only — no script, protocol, or manuscript file was
touched; `results/tr1/wp8/86_*`, `87_*`, and `results/tr1/wp7/` were not read
or written.

## R3.4 and R3.5 (verbatim)

> **R3.4** The radius-selection procedure and the final within-cluster
> validation rely substantially on Complete Spatial Randomness or a
> homogeneous Poisson process. However, the target scenarios emphasized in
> the manuscript include Gaussian clusters, non-uniform densities, density
> gradients, and irregular cluster structures. This creates a potential
> mismatch between the underlying statistical assumption and the intended
> applications. Shape-adaptive multi-ball coverage may improve the geometric
> representation of clusters, but it does not fundamentally eliminate the
> local homogeneity assumption embedded in the spatial-randomness test. The
> authors should explain the statistical interpretation of the proposed
> tests under non-uniform intensity, investigate the false-positive behavior
> when the CSR assumption is violated, and consider alternatives such as
> local density normalization or an inhomogeneous spatial-process null model.
>
> **R3.5** The NND-based MC-SRT examines multiple nested candidate radii for
> each observation and repeatedly performs Monte Carlo tests at different
> radii. Holm's procedure is used to combine the tests based on the mean and
> median nearest-neighbor distances, but the manuscript does not sufficiently
> explain how multiplicity across candidate radii or across all observations
> is handled. In addition, the same data are used to select the covering
> radii, construct the cluster structure, and determine whether observations
> are outliers, which may raise concerns about the validity of inference
> after adaptive selection. The authors should clarify which error
> probability is controlled by the reported significance level and provide
> further discussion or theoretical results concerning false-positive
> control, consistency, or detection probability under repeated testing and
> data-adaptive selection.

## Design recap (from `WP8_PROTOCOL.md`, experiment 4, and its dated notes)

- **Settings:** generator in {uniform_control, beta_gradient, gauss_equal,
  gauss_unequal} x d in {3, 10}, n = 200, **n0 = 0 throughout** (`Y = rep(1,
  n)`). 8 settings. `gauss_equal` (isotropic, sigma = 1) was added on
  2026-09-05 so that `gauss_unequal - gauss_equal` isolates anisotropy
  specifically from a Gaussian's own radial density gradient.
- **Reps:** RK-based methods (U-MCCD, SU-MCCD) at **40 reps**; every other
  method (UN-MCCD, SUN-MCCD, LOF, DBSCAN, MST, ODIN, iForest, and the 4-member
  MST threshold sweep MST@1.05/1.40/1.60/2.00) at **100 reps** — a launch
  decision recorded in `WP8_REVERIFICATION.md` to control the RK bracket
  search's 8-140x cost blow-up on single-population data (`WP8_PROTOCOL.md`'s
  "88 cost finding"). 13 methods total (9 default + 4 MST-sweep members not
  duplicating the fixed `MST`=1.2 row).
- **Metric:** because `n0 = 0`, `evaluate()`/`count_scores2` cannot be called
  (it indexes `label_pred[(n-n0+1):n]`, which is `(n+1):n` at n0=0, reading
  past the vector). The reported quantity is the raw **flag rate** =
  `mean(score >= threshold)` over all n points, using the same
  `REAL_DATA_THRESHOLDS`/`WP0_THRESHOLDS` as every other WP8 script. With no
  true outliers, TNR = 1 − flag_rate exactly, so flag_rate is the whole
  story — it **is** the empirical false-positive rate of the whole detection
  procedure on this generator.
- **alpha resolvers** (`shared/harness.R`, unchanged, never overridden):
  RK (`U-MCCD`, `SU-MCCD`): alpha = 1% for d < 10, 0.1% for d >= 10.
  `UN-MCCD`: 15%/10%/5%/1%/0.1% at d <= 2/4/9/19/>=19 respectively (10% at
  d=3, **1%** at d=10). `SUN-MCCD`: identical to UN-MCCD below d=10, but
  forces 0.1% already at d=10 (one step earlier than UN-MCCD). **Note:** the
  task brief's "10% at d=3, 0.1% at d=10 for the NND methods" holds for
  SUN-MCCD but not for UN-MCCD, which is still at 1% (not 0.1%) at d=10 per
  the resolver above — both are reported below.

## Integrity results

`count.fields` pass on both output files — every row matches its header's
field count:

| File | Header fields | Data rows | Bad rows |
|---|---|---|---|
| `88_csr_violation.csv` | 10 | 9,440 | 0 |
| `88_csr_violation_clusters.csv` | 8 | 448,000 | 0 |

- **Expected metric rows:** 8 settings x (40 reps x 2 RK methods + 100 reps x
  11 other method/variant rows) = 8 x (80 + 1,100) = **9,440** — matches
  exactly, before *and* after de-duplication (there is nothing to
  de-duplicate).
- **Duplicates on `(setting_id, d, rep, method)`:** 0.
- **`status != "ok"` rows:** 0 — every one of the 9,440 cells succeeded on
  its first attempt; no error rows, nothing retried.
- **Completeness:** all 104 distinct `(setting_id, method)` cells (8 settings
  x 13 methods) have exactly their expected rep count (40 for U-MCCD/SU-MCCD,
  100 for the other 11). No duplicate seeds within any `(setting_id,
  method)`.
- **No `88_csr_violation_done.csv` file exists** — this matches the protocol:
  88 is single-row-per-cell and gated directly by `has_result()` on the main
  CSV; only 85 (which writes 20 sub-rows per cell) uses a companion done
  file.
- **Cluster file:** 448,000 rows = 8 settings x (40 reps x 2 RK + 100 reps x
  2 NND MCCD methods) x 200 points = 8 x 280 x 200 — matches exactly.

The grid completed cleanly. Nothing was excluded from the tables below on
integrity grounds.

## Script's own `--summarize`

`Rscript revision_experiments/tr1/88_wp8_csr_violation.R --summarize` wrote
`results/tr1/wp8/88_csr_violation_summary.csv` (104 rows: one per
`(generator, d, method)`), plus a console MST-threshold report. Its
`flag_rate_mean`/`flag_rate_se` columns were cross-checked against this
analysis's independent computation and match to full precision (spot-checked
on all 8 `uniform_control x {U/SU/UN/SUN}-MCCD` rows). Its `delta_vs_control`
column, however, is computed by **merging on `rep`** — i.e. as a *paired*
difference between a generator's rep-*k* draw and the control's rep-*k*
draw. Per `WP8_REVERIFICATION.md`'s analysis-time item, **this pairing is
not meaningful**: each generator/d combination gets its own `s_i`-indexed
seed offset (`BASE_SEED + 100000*s_i + rep`), so "rep 5 of beta_gradient"
and "rep 5 of uniform_control" are statistically independent draws that
merely share an integer label — there is no shared randomness to exploit by
pairing. The tables below instead report the delta as a **difference of two
independent means**, `SE = sqrt(SE1^2 + SE2^2)`, computed separately in
`wp8_88_summary.R` (scratchpad). The point estimates of the delta are
numerically close between the two methods (pairing two independent samples
by an arbitrary shared index does not bias the mean), but only the
independent-means SE is the statistically correct one to quote.

Its "best MST threshold" console report (`which.min(|delta_vs_control_mean|)`
among the swept thresholds) is also reproduced and corrected below — at d =
10, `MST@2.00`'s control flag rate is 0.0007 (effectively flags nothing), so
it wins "closest to zero delta" by flagging almost nothing under any
generator, not by tracking the CSR violation. `WP8_REVERIFICATION.md` flagged
exactly this failure mode.

## R3.5 answer: alpha (a per-ball test level) vs. the empirical procedure-level false-positive rate

**Control (uniform_control) flag rate — the empirical false-positive rate of
the whole procedure, with no outliers present:**

| method | d | n_reps | flag_rate (mean) | SE | per-ball alpha | flag_rate / alpha |
|---|---|---|---|---|---|---|
| SU-MCCD | 3 | 40 | 0.0090 | 0.00185 | 1% | 0.90 |
| SUN-MCCD | 3 | 100 | 0.00955 | 0.00112 | 10% | 0.096 |
| U-MCCD | 3 | 40 | 0.01963 | 0.00272 | 1% | 1.96 |
| UN-MCCD | 3 | 100 | 0.02330 | 0.00195 | 10% | 0.233 |
| SU-MCCD | 10 | 40 | 0.04662 | 0.00417 | 0.1% | **46.6** |
| SUN-MCCD | 10 | 100 | 0.02655 | 0.00179 | 0.1% | **26.6** |
| U-MCCD | 10 | 40 | 0.08325 | 0.00959 | 0.1% | **83.3** |
| UN-MCCD | 10 | 100 | 0.05560 | 0.00313 | 1% | **5.56** |

**Reading.** At d = 3 the per-ball alpha is a loose but not absurd guide to
the procedure-level rate: RK methods (U-MCCD ~2x alpha, SU-MCCD ~0.9x alpha)
track it within a factor of 2, while the NND methods run noticeably more
conservatively than their nominal 10% (SUN-MCCD at ~1/10 of alpha, UN-MCCD at
~1/4). At d = 10 the relationship collapses for every method: the observed
flag rate is **5.6 to 83 times** the declared per-ball alpha, worst for
U-MCCD (83x) and SU-MCCD (47x), and even the least-inflated method (UN-MCCD)
is still 5.6x its nominal level. This is a direct, numeric answer to R3.5:
alpha is the significance level of one covering-ball's spatial-randomness
test, not the false-positive rate of the *combined* clustering-plus-testing
procedure applied to n = 200 points and (for the RK/NND density calibration)
many candidate balls per cluster; the two diverge by an order of magnitude
or more in higher dimension. The reviewer's "does the same data used to
select radii and test for outliers invalidate the reported significance
level" concern is confirmed empirically, not just in principle.

## R3.4 answer: flag-rate rise above control under CSR violation

**Delta vs. uniform_control at the same d, difference of independent means**
(full per-generator/method/d table in `wp8_88_delta.csv`; MCCD + baselines
only shown here, MST sweep in its own section below):

| generator | d | method | control | gen. rate | delta | SE(delta) | delta/SE |
|---|---|---|---|---|---|---|---|
| beta_gradient | 3 | U-MCCD | 0.0196 | 0.1179 | **+0.098** | 0.0078 | 12.6 |
| beta_gradient | 3 | UN-MCCD | 0.0233 | 0.0863 | +0.063 | 0.0042 | 15.1 |
| beta_gradient | 3 | SU-MCCD | 0.0090 | 0.0461 | +0.037 | 0.0045 | 8.3 |
| beta_gradient | 3 | SUN-MCCD | 0.0096 | 0.0402 | +0.031 | 0.0030 | 10.3 |
| beta_gradient | 10 | U-MCCD | 0.0833 | 0.2313 | **+0.148** | 0.0166 | 8.9 |
| beta_gradient | 10 | UN-MCCD | 0.0556 | 0.1722 | +0.117 | 0.0079 | 14.7 |
| beta_gradient | 10 | SU-MCCD | 0.0466 | 0.1601 | +0.114 | 0.0114 | 9.9 |
| beta_gradient | 10 | SUN-MCCD | 0.0266 | 0.0756 | +0.049 | 0.0041 | 11.9 |
| gauss_equal | 3 | U-MCCD | 0.0196 | 0.1616 | **+0.142** | 0.0089 | 15.9 |
| gauss_equal | 3 | UN-MCCD | 0.0233 | 0.1249 | +0.102 | 0.0057 | 17.7 |
| gauss_equal | 10 | U-MCCD | 0.0833 | 0.3028 | **+0.220** | 0.0184 | 11.9 |
| gauss_equal | 10 | UN-MCCD | 0.0556 | 0.2147 | +0.159 | 0.0093 | 17.1 |
| gauss_unequal | 3 | U-MCCD | 0.0196 | 0.2284 | **+0.209** | 0.0108 | 19.4 |
| gauss_unequal | 3 | UN-MCCD | 0.0233 | 0.1714 | +0.148 | 0.0065 | 22.9 |
| gauss_unequal | 10 | U-MCCD | 0.0833 | 0.3103 | **+0.227** | 0.0194 | 11.7 |
| gauss_unequal | 10 | UN-MCCD | 0.0556 | 0.2390 | +0.183 | 0.0086 | 21.4 |
| — | — | DBSCAN | 0 | 0 | 0 | 0 | — (by construction, see below) |

Every non-DBSCAN core method's flag rate rises well above its own control
under every CSR-violating generator, at both d = 3 and d = 10, by 8-37 SEs —
none of these deltas are noise. iForest is the one exception with a
*negative* delta at every generator/d combination (−0.013 to −0.029): its
flag rate *drops* under the gradient/anisotropic regimes relative to the
uniform control (its threshold is a fixed score cut, and these regimes
concentrate score mass below it rather than above).

## (c) gauss_unequal − gauss_equal: anisotropy isolated from the Gaussian's own radial gradient

| method | d | gauss_equal | gauss_unequal | delta (anisotropy only) | SE | delta/SE |
|---|---|---|---|---|---|---|
| U-MCCD | 3 | 0.1616 | 0.2284 | **+0.067** | 0.0134 | 5.0 |
| SU-MCCD | 3 | 0.0725 | 0.1210 | +0.049 | 0.0081 | 6.0 |
| UN-MCCD | 3 | 0.1249 | 0.1714 | +0.046 | 0.0082 | 5.7 |
| SUN-MCCD | 3 | 0.0469 | 0.0935 | +0.047 | 0.0051 | 9.2 |
| ODIN | 3 | 0.1227 | 0.1016 | −0.021 | 0.0024 | −8.8 |
| U-MCCD | 10 | 0.3028 | 0.3103 | +0.008 | 0.0231 | 0.3 |
| SU-MCCD | 10 | 0.1974 | 0.2246 | +0.027 | 0.0190 | 1.4 |
| UN-MCCD | 10 | 0.2147 | 0.2390 | +0.024 | 0.0119 | 2.1 |
| SUN-MCCD | 10 | 0.0780 | 0.1146 | **+0.037** | 0.0053 | 6.9 |
| ODIN | 10 | 0.2849 | 0.2325 | −0.052 | 0.0033 | −15.7 |
| LOF | 10 | 0.0024 | 0.0081 | +0.006 | 0.0008 | 7.5 |

(full table, all 13 methods x 2 d, in `wp8_88_aniso.csv`.) At d = 3, all four
MCCD variants respond to anisotropy alone with a further +4.6 to +6.7 pp
flag-rate rise beyond the isotropic Gaussian's own gradient effect (5-9 SEs).
At d = 10 the RK-based U-MCCD/SU-MCCD's anisotropy-specific effect is no
longer resolvable at 40 reps (0.3 and 1.4 SEs respectively — this is exactly
the reduced-rep cost the launch decision accepted), but the 100-rep NND
methods still show a clear, resolvable anisotropy effect (UN-MCCD +2.1 SEs,
SUN-MCCD +6.9 SEs). ODIN moves in the *opposite* direction at both d (−9 to
−16 SEs): its flag rate is lower under anisotropic than isotropic Gaussian
data, i.e. anisotropy specifically makes ODIN flag less, not more.

## (d) Most / least sensitive methods (core 9, excluding DBSCAN)

Ranked by the largest observed |delta| against control across the three
violating generators and both d (full per-cell table in
`wp8_88_sensitivity.csv`):

| rank | method | max |delta| | at | mean delta (3 generators x 2 d pooled) |
|---|---|---|---|---|
| 1 (most sensitive) | U-MCCD | +0.227 | gauss_unequal, d=10 | +0.174 |
| 2 | UN-MCCD | +0.183 | gauss_unequal, d=10 | +0.129 |
| 3 | SU-MCCD | +0.178 | gauss_unequal, d=10 | +0.109 |
| 4 | ODIN | +0.132 | gauss_equal, d=10 | +0.082 |
| 5 | SUN-MCCD | +0.088 | gauss_unequal, d=10 | +0.057 |
| 6 | MST (1.2) | +0.080 | gauss_equal, d=10 | +0.032 |
| 7 | LOF | +0.055 | gauss_unequal, d=3 | +0.024 |
| 8 (least sensitive) | iForest | −0.029 | gauss_equal, d=3 | **−0.022** |

**Reading.** The two RK-based MCCD variants (U-MCCD, SU-MCCD) and the
1%-alpha UN-MCCD are the most sensitive to CSR violation by a wide margin —
consistent with their per-ball test running closer to its nominal alpha at
d=3 (Section R3.5 above) and therefore having more headroom to inflate.
**SUN-MCCD, the paper's headline recommendation, is the least sensitive of
the four MCCD variants** (mean delta +0.057 vs. +0.109 to +0.174 for the
other three) and sits below ODIN, though still above MST and well above
LOF. In absolute terms SUN-MCCD's flag rate never exceeds 0.115 across any
of the 6 non-control settings, versus 0.31 for U-MCCD and 0.285 for ODIN —
so while its *relative* rise above its own (already low) control is not the
smallest in the ranking, its *absolute* false-positive exposure under CSR
violation is the second-lowest of the 8 non-DBSCAN methods, behind only
iForest and LOF. iForest is the outlier of the ranking in kind, not just
degree: its delta is negative at all 6 settings — it is not "robust" to the
CSR violation so much as pointed the opposite way, presumably because its
isolation-path score threshold (0.55, fixed across all generators) happens
to sit above more of the score mass once density concentrates away from the
box's uniform fill.

## MST threshold sweep

Control (uniform_control) flag rate per swept threshold — this is what
determines eligibility for "best":

| method | d | control flag rate | SE |
|---|---|---|---|
| MST (1.2, fixed) | 3 | 0.388 | 0.0053 |
| MST@1.05 | 3 | 0.680 | 0.0054 |
| MST@1.40 | 3 | 0.185 | 0.0040 |
| MST@1.60 | 3 | 0.096 | 0.0029 |
| MST@2.00 | 3 | 0.035 | 0.0017 |
| MST (1.2, fixed) | 10 | 0.130 | 0.0026 |
| MST@1.05 | 10 | 0.492 | 0.0042 |
| MST@1.40 | 10 | 0.022 | 0.0011 |
| MST@1.60 | 10 | 0.005 | 0.0005 |
| **MST@2.00** | **10** | **0.0007** | **0.0002** |

**Excluded from "best" selection: MST@2.00 at d = 10** (control flag rate
0.0007, i.e. it flags essentially nothing on the control regardless of
generator — `WP8_REVERIFICATION.md`'s named failure mode). No other
sweep member falls below the 0.005 exclusion bar.

"Best" = swept threshold whose delta-vs-control is closest to zero, among
eligible members, vs. the always-reported fixed MST(1.2) row:

| generator | d | fixed MST(1.2) rate / delta | best-eligible member | rate / delta | its own control |
|---|---|---|---|---|---|
| beta_gradient | 3 | 0.382 / −0.005 | MST@2.00 | 0.042 / +0.007 | 0.035 |
| beta_gradient | 10 | 0.183 / +0.053 | MST@1.60 | 0.014 / +0.009 | 0.005 |
| gauss_equal | 3 | 0.380 / −0.008 | MST@2.00 | 0.045 / +0.010 | 0.035 |
| gauss_equal | 10 | 0.209 / +0.080 | MST@1.05 | 0.497 / +0.005 | 0.491 |
| gauss_unequal | 3 | 0.381 / −0.007 | MST@1.40 | 0.196 / +0.011 | 0.185 |
| gauss_unequal | 10 | 0.206 / +0.076 | MST@1.05 | 0.490 / −0.002 | 0.491 |

(see `wp8_88_summary.R` console output / `wp8_88_delta.csv` for the exact
figures; the table above matches the script's own console report once
MST@2.00 at d=10 is excluded — at d=10, once MST@2.00 is excluded, "best"
falls to MST@1.05 or MST@1.60 depending on generator, both with
non-trivial control flag rates of 0.005-0.49, i.e. a genuine, not
degenerate, comparison point.) **Without the exclusion**, the script's raw
`--summarize` output picks `MST@2.00` as "best" for every d=10 generator
(delta magnitudes 0.0002-0.001) purely because it flags ~0% of everything —
this would misleadingly read as "MST at threshold 2.0 is nearly immune to
CSR violation," when in fact it has no discriminating power left at that
threshold in d=10 to be sensitive to anything.

## DBSCAN (by construction)

DBSCAN's flag rate is **0.000 at every one of the 8 settings** (all 800
"ok" rows, `note = "oracle contamination = 0; flag rate 0 by construction"`).
This is fed by `Y = rep(1, n)` (no outliers), which sets DBSCAN's internal
oracle-contamination-based quantile to 0 — it is not evidence that DBSCAN
detected the CSR violation or resisted it; it never had the opportunity to
flag anything. DBSCAN is excluded from every sensitivity comparison and
delta table above and must never be read as "0 false positives" in the R3.4
sense.

## Caveats

- **40-rep RK reduction.** U-MCCD/SU-MCCD ran at 40 reps (vs. 100 for every
  other method), a launch-time cost decision (`WP8_PROTOCOL.md`'s "88 cost
  finding": the RK density-calibration bracket search is 8-140x slower on
  88's single-population generators than on the two-cluster settings the
  original per-cell cost estimate assumed). Observed control-flag-rate SEs
  for these two methods range 0.0018-0.0096 (~0.2-1 pp), and every reported
  R3.4 delta for them is 8-19 SEs from zero — the reduced rep count does not
  threaten the qualitative finding, but it does widen individual confidence
  intervals visibly relative to the 100-rep methods (e.g., the gauss_unequal
  d=10 anisotropy-only effect for U-MCCD is only 0.3 SEs and cannot be
  distinguished from zero at this rep count, unlike the same contrast for
  the 100-rep SUN-MCCD, which resolves at 6.9 SEs).
- **Independent-means delta, not paired.** The script's own `--summarize`
  merges generator and control rows on `rep`, which looks like a paired
  comparison but is not one — each generator/d setting has its own
  seed offset, so "rep k" carries no shared randomness across generators.
  All delta/SE figures reported here use the independent-means formula
  `SE = sqrt(SE1^2 + SE2^2)` instead; point estimates barely move, but the
  SE is the one that is actually valid to quote.
- **MST fixed threshold vs. per-data-set tuning.** The paper's real-data
  pipeline tunes MST's `thresh` per data set (1.05-1.6); this experiment's
  fixed `MST` row uses 1.2 throughout, which is well inside that tuned
  range but is not itself tuned to any of these four synthetic generators.
  The threshold sweep exists precisely to show how sensitive the "MST is
  fine/not fine under CSR violation" reading is to that choice — see the
  sweep table above, especially the exclusion of near-zero-control members.
- **RK-degeneracy provenance notes are not this experiment's concern.**
  `HANDOFF_FROM_TR2.md`'s two RK high-dimensional-degeneracy provenance
  notes concern the *real-data* quantile tables and are unrelated to this
  synthetic, no-outlier CSR experiment; not referenced further here.
- This file reports flag rates and their standard errors as descriptive
  statistics of the observed sample; no additional hypothesis test, p-value,
  or multiplicity correction beyond what is stated above has been computed
  or is implied.

## Source files

- Raw data (not committed, gitignored): `88_csr_violation.csv`,
  `88_csr_violation_clusters.csv`, script's own summary
  `88_csr_violation_summary.csv` — all in this directory.
- Analyst's independent computation: `wp8_88_summary.R`, and its CSV outputs
  `wp8_88_agg.csv`, `wp8_88_delta.csv`, `wp8_88_control_vs_alpha.csv`,
  `wp8_88_aniso.csv`, `wp8_88_sensitivity.csv` — temporary
  scratch directory, not retained and not part of this repository.
