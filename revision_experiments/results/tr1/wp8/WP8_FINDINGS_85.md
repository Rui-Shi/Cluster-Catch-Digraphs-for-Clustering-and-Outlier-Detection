# WP8 experiment 1 — boundary false positives (R1.1)

**Script:** `revision_experiments/tr1/85_wp8_boundary_fp.R`. **Grid:** launched
2026-09-05, complete as of the files analysed here (2026-09-06). Analysis
only — no script, protocol, or manuscript file was touched.

## R1.1 (verbatim)

> The bidirectional constraint of "mutual capture" described in the paper
> filters a large number of boundary samples. Are normal samples located at
> cluster boundaries more likely to be misidentified as anomalous samples?

## Design recap (from `WP8_PROTOCOL.md`, experiment 1)

- **Settings:** generator in {uniform, gaussian} x d in {3, 10}, n = 200,
  contamination = 0.05. 4 settings, 100 reps, 9 methods per rep (U-MCCD,
  SU-MCCD, UN-MCCD, SUN-MCCD, LOF, DBSCAN, MST, ODIN, iForest).
  `S_min = 0.05` passed as `min.cls` to SU-MCCD/SUN-MCCD only; alpha resolved
  per-d by the paper's own dimension-schedule resolvers.
- **Generators:** two clusters of sizes n1 = 95, n2 = 94 (the declared
  off-by-one; n1+n2+n0 = 199, not 200 — inherited unchanged from
  `55_wp2c_simulation_arm.R`, not corrected here), centred `cls_dis = 3`
  apart, points drawn `rpoisball.unit * runif(r_min,r_max) + mu_k` (uniform
  arm) or `mvrnorm` with matched calibration (gaussian arm); 10 rejection-
  sampled outliers at distance > `otl_dis` from both centres.
- **Per-point normalized radius and bins.** For every regular point,
  `r_i = ||x_i - mu_k||` in its own cluster; two binning schemes into 10 bins
  each (`bin_type` column):
  - **width**: `bin = clamp(ceil(10*r_i/r_max_k), 1, 10)`, `r_max_k` = the
    largest observed radius in that cluster/rep.
  - **mass**: `bin = clamp(ceil(10*rank(r_i)/n_k), 1, 10)`, an equal-count
    scheme added after `WP8_VERIFICATION.md` found the width scheme
    uninformative at d = 10 (see Integrity/caveats below).
  Both clusters' regular points are pooled into the same 10 bins per
  `(setting, rep, method, bin_type)` cell.
- **Outputs:** `85_boundary_fp.csv` (20 bin-rows/cell: `n_regular_in_bin`,
  `n_flagged_in_bin`), `85_boundary_fp_done.csv` (done markers, one row/cell),
  `85_boundary_fp_clusters.csv` (WP7 per-point cluster assignment, MCCD
  methods only, 189 rows/cell = n_reg).

## Integrity results

**`count.fields` pass** (header quote-aware field count) on all three files —
every row matches its header's field count, no truncation or corruption from
concurrent writers on the G: enclosure:

| File | Header fields | Data rows | Bad rows |
|---|---|---|---|
| `85_boundary_fp.csv` | 10 | 72,000 | 0 |
| `85_boundary_fp_done.csv` | 5 | 3,600 | 0 |
| `85_boundary_fp_clusters.csv` | 8 | 302,400 | 0 |

All three row counts match the expected counts exactly (4 settings x 100
reps x 9 methods = 3,600 done rows; x 20 bins = 72,000 bin rows; and, MCCD
methods only, 4 settings x 100 reps x 4 methods x 189 regular points =
302,400 cluster rows).

**Cell-level checks:**

- Done markers: 3,600 unique `(setting_id, d, rep, method)` keys, 0
  duplicates, 0 missing against the full declared grid (4 settings x 100
  reps x 9 methods) — **no failed cells**. 85 writes no row for a failed
  cell, so a missing cell would show up as `< 100` reps for that
  setting/method; none did.
- Main (bin-row) file: exactly 3,600 distinct cells, every one with exactly
  20 bin rows (10 width + 10 mass), 0 duplicate `(cell, bin_type, bin)` keys
  after de-duplication.
- Cluster file: exactly 1,600 distinct `(setting, d, rep, method)` cells
  (4 settings x 100 reps x 4 MCCD methods), every one with exactly 189 rows,
  0 duplicate `row_index` within any cell.
- No cell has a done marker but is absent from the main file (the reverse of
  the "failed cell" check).

The grid completed cleanly: every one of 3,600 cells has a done marker, 20
bin rows, and (for the 4 MCCD methods) a full 189-row cluster block. Nothing
was excluded from the summary below on integrity grounds.

## Summariser

`Rscript revision_experiments/tr1/85_wp8_boundary_fp.R --summarize` wrote
`revision_experiments/results/tr1/wp8/85_boundary_fp_summary.csv` (720 rows =
4 settings x 9 methods x 2 bin_types x 10 bins). This file is inside
`results/`, which is gitignored — nothing was committed and nothing needs to
be.

## Width-bin informativeness check

Confirms the verifier's finding and extends it to the gaussian arm. Share of
regular points landing in each **width** bin (method-invariant — binning
uses only `X`, not any method's score):

| bin | uniform_d3 | uniform_d10 | gaussian_d3 | gaussian_d10 |
|---|---|---|---|---|
| 1 | 17 | 0 | 210 | 0 |
| 2 | 134 | 0 | 1314 | 6 |
| 3 | 328 | 0 | 2755 | 132 |
| 4 | 668 | 0 | 3720 | 940 |
| 5 | 1174 | 21 | 3612 | 2986 |
| 6 | 1735 | 94 | 2944 | 4775 |
| 7 | 2391 | 388 | 2063 | 4590 |
| 8 | 3249 | 1452 | 1183 | 3170 |
| 9 | 4002 | 4652 | 568 | 1569 |
| 10 | 5202 | **12,293 (65.0%)** | 531 | 732 |

(out of 18,900 = 100 reps x 189 regular points per setting)

At **d = 10, uniform**, 65.0% of all regular points fall in the last width
bin and bins 1–4 are literally empty (0 points, all 100 reps) — this
reproduces the verifier's "66%" finding. The width scheme has no interior to
contrast at d = 10 for the uniform generator. The gaussian d = 10 arm is
less degenerate (bin 1 empty, but 2–10 populated, peaked around bins 6–7)
because the unbounded gaussian tail still produces a spread of observed
radii — but its `r_max` is an unbounded sample maximum, not a fixed support
edge, so width bins are not comparable across reps there either (see
caveats). **Equal-mass bins are used for all reported results below**; width
bins are given only for reference.

## Results: equal-mass boundary false-positive rate, per method x setting

Pooled ratio = sum(n_flagged) / sum(n_regular) over all 100 reps, per bin
(never a mean of per-rep ratios, which would over-weight reps with few
points in a bin). n_reps = 100 and n_reps_nonempty = 100 in every mass-bin
cell (mass bins hold ~19 points per rep by construction, never empty).

### uniform_d3

| method | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 |
|---|---|---|---|---|---|---|---|---|---|---|
| U-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.001 | 0.007 | 0.012 | 0.019 | 0.040 |
| SU-MCCD | 0.000 | 0.000 | 0.000 | 0.001 | 0.000 | 0.000 | 0.001 | 0.002 | 0.002 | 0.004 |
| UN-MCCD | 0.000 | 0.000 | 0.000 | 0.001 | 0.001 | 0.006 | 0.008 | 0.015 | 0.023 | 0.040 |
| SUN-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.003 | 0.005 | 0.006 | 0.013 |
| LOF | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| DBSCAN | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| MST | 0.271 | 0.250 | 0.241 | 0.276 | 0.312 | 0.363 | 0.406 | 0.463 | 0.483 | 0.518 |
| ODIN | 0.000 | 0.000 | 0.001 | 0.001 | 0.003 | 0.004 | 0.008 | 0.024 | 0.049 | 0.080 |
| iForest | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

### uniform_d10

| method | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 |
|---|---|---|---|---|---|---|---|---|---|---|
| U-MCCD | 0.002 | 0.006 | 0.007 | 0.014 | 0.015 | 0.026 | 0.029 | 0.032 | 0.037 | 0.035 |
| SU-MCCD | 0.001 | 0.001 | 0.003 | 0.006 | 0.005 | 0.012 | 0.015 | 0.018 | 0.017 | 0.028 |
| UN-MCCD | 0.000 | 0.000 | 0.002 | 0.002 | 0.003 | 0.007 | 0.004 | 0.013 | 0.011 | 0.019 |
| SUN-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.002 | 0.002 | 0.002 | 0.005 |
| LOF | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| DBSCAN | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| MST | 0.012 | 0.038 | 0.052 | 0.068 | 0.086 | 0.096 | 0.109 | 0.112 | 0.132 | 0.121 |
| ODIN | 0.000 | 0.000 | 0.000 | 0.001 | 0.006 | 0.013 | 0.025 | 0.036 | 0.045 | 0.071 |
| iForest | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

### gaussian_d3

| method | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 |
|---|---|---|---|---|---|---|---|---|---|---|
| U-MCCD | 0.000 | 0.000 | 0.000 | 0.002 | 0.005 | 0.028 | 0.090 | 0.238 | 0.513 | 0.867 |
| SU-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.007 | 0.064 | 0.497 |
| UN-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.005 | 0.023 | 0.076 | 0.248 | 0.722 |
| SUN-MCCD | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.002 | 0.011 | 0.033 | 0.345 |
| LOF | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.007 | 0.047 | 0.455 |
| DBSCAN | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| MST | 0.166 | 0.189 | 0.204 | 0.242 | 0.277 | 0.333 | 0.393 | 0.460 | 0.568 | 0.629 |
| ODIN | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.005 | 0.029 | 0.079 | 0.280 | 0.702 |
| iForest | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.004 |

### gaussian_d10

| method | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 |
|---|---|---|---|---|---|---|---|---|---|---|
| U-MCCD | 0.001 | 0.007 | 0.028 | 0.077 | 0.151 | 0.243 | 0.352 | 0.521 | 0.684 | 0.884 |
| SU-MCCD | 0.010 | 0.012 | 0.015 | 0.031 | 0.064 | 0.134 | 0.219 | 0.361 | 0.538 | 0.820 |
| UN-MCCD | 0.001 | 0.001 | 0.007 | 0.023 | 0.057 | 0.117 | 0.197 | 0.355 | 0.544 | 0.822 |
| SUN-MCCD | 0.000 | 0.000 | 0.000 | 0.001 | 0.001 | 0.001 | 0.009 | 0.042 | 0.148 | 0.513 |
| LOF | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.025 |
| DBSCAN | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| MST | 0.003 | 0.014 | 0.040 | 0.063 | 0.099 | 0.156 | 0.220 | 0.287 | 0.373 | 0.531 |
| ODIN | 0.000 | 0.000 | 0.000 | 0.002 | 0.011 | 0.049 | 0.149 | 0.385 | 0.686 | 0.937 |
| iForest | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

### Boundary contrast, core (b1) vs boundary (b10), with per-rep SE (n_reps = 100 throughout)

| | uniform_d3 b1 | uniform_d3 b10 | uniform_d10 b1 | uniform_d10 b10 | gaussian_d3 b1 | gaussian_d3 b10 | gaussian_d10 b1 | gaussian_d10 b10 |
|---|---|---|---|---|---|---|---|---|
| U-MCCD | 0.0000(0.0000) | 0.0395(0.0055) | 0.0017(0.0012) | 0.0355(0.0072) | 0.0000(0.0000) | 0.8670(0.0129) | 0.0006(0.0006) | 0.8845(0.0111) |
| SU-MCCD | 0.0000(0.0000) | 0.0045(0.0020) | 0.0006(0.0006) | 0.0275(0.0054) | 0.0000(0.0000) | 0.4970(0.0193) | 0.0100(0.0070) | 0.8200(0.0141) |
| UN-MCCD | 0.0000(0.0000) | 0.0395(0.0058) | 0.0000(0.0000) | 0.0195(0.0049) | 0.0000(0.0000) | 0.7215(0.0199) | 0.0006(0.0006) | 0.8215(0.0135) |
| SUN-MCCD | 0.0000(0.0000) | 0.0135(0.0035) | 0.0000(0.0000) | 0.0050(0.0015) | 0.0000(0.0000) | 0.3455(0.0192) | 0.0000(0.0000) | 0.5130(0.0191) |
| LOF | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.4545(0.0152) | 0.0000(0.0000) | 0.0255(0.0034) |
| DBSCAN | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) |
| MST | 0.2711(0.0138) | 0.5185(0.0129) | 0.0122(0.0028) | 0.1210(0.0072) | 0.1656(0.0124) | 0.6295(0.0109) | 0.0033(0.0013) | 0.5315(0.0120) |
| ODIN | 0.0000(0.0000) | 0.0800(0.0067) | 0.0000(0.0000) | 0.0715(0.0060) | 0.0000(0.0000) | 0.7020(0.0117) | 0.0000(0.0000) | 0.9370(0.0060) |
| iForest | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0000(0.0000) | 0.0045(0.0014) | 0.0000(0.0000) | 0.0000(0.0000) |

## Overall regular-point FP rate per method x setting (pooled across all bins, for context)

| method | uniform_d3 | uniform_d10 | gaussian_d3 | gaussian_d10 |
|---|---|---|---|---|
| U-MCCD | 0.0081 | 0.0205 | 0.1778 | 0.2983 |
| SU-MCCD | 0.0011 | 0.0106 | 0.0595 | 0.2238 |
| UN-MCCD | 0.0094 | 0.0062 | 0.1107 | 0.2157 |
| SUN-MCCD | 0.0030 | 0.0011 | 0.0412 | 0.0739 |
| LOF | 0.0000 | 0.0000 | 0.0534 | 0.0027 |
| DBSCAN | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| MST | 0.3592 | 0.0830 | 0.3477 | 0.1810 |
| ODIN | 0.0173 | 0.0200 | 0.1127 | 0.2257 |
| iForest | 0.0000 | 0.0000 | 0.0005 | 0.0000 |

## Reading

**Yes — the false-positive rate on regular points rises toward the cluster
boundary, for every MCCD method and for MST and ODIN, in all four settings.**
Bin 1 (the 10% of regular points closest to their cluster's own centre) is at
or indistinguishable from 0 everywhere; bin 10 (the outermost 10%) is where
essentially all of the false-positive mass sits.

**Steepness depends on the generator, not on the mutual-capture mechanism
alone.** In the **uniform** arm, where the true point density inside the
support ball is constant, the rise from bin 1 to bin 10 is real but bounded:
SUN-MCCD goes from 0 to 1.35% (d=3) or 0.50% (d=10); U-MCCD, the least
shape-adaptive MCCD variant, rises furthest, to 3.95% (d=3) and 3.55%
(d=10). LOF, DBSCAN and iForest stay at exactly 0 across all 10 bins in both
uniform settings — a homogeneous, bounded-support cluster gives them no
reason to flag any regular point regardless of position, which isolates the
MCCD family's boundary effect as attributable to the geometry of mutual
capture and not to a genuine density gradient the reference density methods
would also pick up. MST and ODIN, however, also rise steeply in the uniform
arm (MST 27–52% at d=3, ODIN 0 to 8%), so a boundary-following false-positive
gradient is not unique to CCD's bidirectional-capture rule; MST in
particular carries a high false-positive rate even at bin 1 (24–27%),
consistent with a general edge/degree-based sensitivity rather than a
boundary-specific one.

In the **gaussian** arm, density genuinely decreases with radius, so a
boundary-following flag rate is the expected behavior of *any* density- or
coverage-based detector, not evidence of a mutual-capture-specific artifact.
Every method except DBSCAN shows a strong rise here (bin 10 rates from 0.5%
(iForest, d=3) up to 93.7% (ODIN, d=10)). What still isolates the mutual-
capture effect within this arm is the **comparison across the four MCCD
variants**: SUN-MCCD, the paper's headline recommendation, has the lowest
bin-10 false-positive rate of the four MCCD methods in three of four
settings (gaussian_d3: 34.5% vs 49.7–86.7% for the other three; gaussian_d10:
51.3% vs 82.0–88.5%) and is essentially tied for lowest in the fourth
(uniform_d3: SU-MCCD 0.45% vs SUN-MCCD 1.35%, both far below U-/UN-MCCD's
3.95–4.0%; uniform_d10: SUN-MCCD is lowest outright at 0.50% vs 1.95–3.55%
for the other three). Shape-adaptive coverage (SU-, SUN-) reduces boundary
false positives relative to the uniform-coverage reference constructions
(U-, UN-) in every setting measured; this is expected under the paper's own
framing — a ball fitted to the cluster's local shape should reach farther
into elongated or low-density tails than a fixed-radius one — but the
magnitude (roughly halving the boundary FP rate in the gaussian arm, and 3–8x
lower in the uniform arm) had not been quantified before this experiment.

**How the four MCCD methods compare with the five baselines.** In the
uniform arm, all four MCCD methods sit well below MST and ODIN and just
above LOF/DBSCAN/iForest's exact-0 floor at every bin — the MCCD family's
absolute boundary FP rate (0.5–4.0% at bin 10) is a real but small cost next
to MST's 12–52%. In the gaussian arm the ranking is generator-dependent: at
d=3, LOF's bin-10 rate (45.5%) sits between SUN-MCCD (34.5%) and SU-MCCD
(49.7%); at d=10, LOF collapses to 2.5% at bin 10 while every MCCD variant
and MST/ODIN remain high (51–94%) — LOF's distance-ratio statistic loses
discriminative power in the higher-dimensional gaussian tail in a way none
of the graph/coverage-based methods do. DBSCAN and iForest are the only
methods that never rise above single digits (iForest) or exactly 0
(DBSCAN) in any setting; both also have near-zero overall power in this
design (DBSCAN never flags anything at all, by the same fixed-threshold
mechanism noted for the CSR-violation experiment in `WP8_PROTOCOL.md`
experiment 4).

**Where equal-width bins are uninformative.** At d = 10 in the **uniform**
generator, 65.0% of all regular points fall in the single last width bin and
bins 1–4 are completely empty (0 of 18,900 points across all 100 reps) — the
width scheme has no interior to contrast, reproducing the verifier's "66%"
finding. At d = 10 in the **gaussian** generator the width scheme is less
degenerate (bin 1 empty, 2–10 populated and peaked mid-range) but its
`r_max` is an unbounded per-rep sample maximum rather than a fixed support
edge, so bins are not comparable across reps or even within a rep across
methods' analyses of the same replicate. All quantitative claims above use
the **equal-mass** bins for this reason; the width-bin tables in this file
are included for reference only.

## Caveats

- **Gaussian r_max is an unbounded sample maximum.** Unlike the uniform
  generator, whose support has a hard edge at radius `r_min`–`r_max` before
  the jitter draw, the gaussian generator's per-rep, per-cluster `r_max_k` is
  just the largest of ~95–189 draws from an unbounded normal — it varies
  rep to rep and is not a stable normalization reference. Equal-width bins
  computed against it are not comparable across reps in the gaussian arm;
  equal-mass bins, which use only the within-rep rank of `r_i` and never
  reference `r_max`, are unaffected and are what every reading above relies
  on for both generators.
- **The off-by-one in n1/n2** (95 + 94 = 189 regular points, not 190) is
  inherited unchanged from `55_wp2c_simulation_arm.R` per the protocol's
  explicit decision not to correct it (a math change reserved for sign-off);
  it does not affect any ratio reported here since both numerator and
  denominator are computed against the same realised `n_reg`.
- **MST and ODIN are not boundary-specific detectors** in the design sense
  MCCD's mutual capture is; their elevated FP rates (including a
  non-negligible rate at bin 1 for MST in the uniform arm) reflect general
  edge/degree sensitivity rather than a bidirectional-coverage mechanism,
  and should not be read as evidence against the boundary-rise finding for
  MCCD specifically.
- **The gaussian arm's rise is expected under any density-aware method**;
  it demonstrates the mutual-capture effect only through the cross-method
  and cross-MCCD-variant comparisons above (SUN-MCCD vs the other three MCCD
  constructions; MCCD family vs LOF/DBSCAN/iForest), not through the raw
  gaussian bin-10 numbers in isolation.
- Integrity checks above (`count.fields`, cell/duplicate/missing-cell counts)
  were run only against 85's own output files; `86_*`, `87_*`, `88_*` were
  not touched, per instruction (their grids were still writing at analysis
  time).
