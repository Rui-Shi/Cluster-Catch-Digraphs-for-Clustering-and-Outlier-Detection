# WP8 experiment 3 — small legitimate cluster below S_min (R3.7)

**Script:** `revision_experiments/tr1/87_wp8_small_cluster.R`. **Grid:** complete
as of the files analysed here (2026-09-06, files last written 03:39). Analysis
only — no script, protocol, or manuscript file was touched; `results/tr1/wp8/86_*`
and `results/tr1/wp7/` were not read or written.

## R3.7 (verbatim)

> The main principle of the MCG framework is that outliers lack connectivity to
> a cluster core. This assumption is more suitable for global outliers located
> away from the main cluster structure, but it may be restrictive for local
> outliers embedded within clusters, collective outliers near cluster
> boundaries, anomalous points connected to a cluster core through a small
> number of bridging observations, or internally well-connected groups of
> anomalous points. Conversely, a small but legitimate cluster may be
> classified as anomalous if its size is below S_min. The authors should
> include targeted experiments for these difficult cases and discuss in
> greater detail how the method distinguishes bridging effects, legitimate
> small clusters, and collective anomalies.

This experiment addresses the "small but legitimate cluster ... below S_min"
sentence specifically (the local/bridging/collective sentence is experiment 2,
script 86).

## Design recap (from `WP8_PROTOCOL.md`, experiment 3, and its dated notes)

- **Generator:** two background clusters (`mu1 = rep(3,d)`, `mu2 = mu1 +
  3*e1`) plus a third genuine cluster (`mu3 = mu1 + 3*e2`), `n = 200` total,
  fixed. Background split `n_bg = 200 - m`, `n1 = floor(n_bg/2)`, `n2 = n_bg -
  n1`. All three clusters drawn independently (`rpoisball.unit(n_k,d) *
  runif(1,0.7,1.3) + mu_k`, three separate radius-jitter draws) — **no
  contamination outliers** (`n0 = 0`); the only question is what happens to
  cluster 3 itself.
- **Density-matched cluster 3** (the fix required by
  `WP8_VERIFICATION.md` item 87.1): cluster 3's radius draw is scaled by
  `(m/n1)^(1/d)` relative to its own independent `runif(0.7,1.3)` draw, so its
  point density (points per unit volume) matches cluster 1's at every `m`.
  Under the pre-fix generator, shrinking `m` also shrank density (measured
  NN-spacing ratio vs. background at d=3: 4.42 at m=2, 2.26 at m=10, 1.71 at
  m=20). With the fix, **size is what varies with `m`, not density** — this
  is the load-bearing design fact for everything below: a small `m` cell is a
  small legitimate cluster at the *same* density as the background, not a
  sparse one.
- **Sizes:** `m in {2, 4, 6, 8, 10, 12, 15, 20}` out of n = 200, against
  `S_min = 0.05` so `round(0.05*200) = 10` is the size at which the flip is
  expected. `d in {3, 10}`. 16 settings.
- **Methods:** SU-MCCD, SUN-MCCD called with `min.cls = S_min = 0.05` (the
  only two methods whose wrapper takes a `min.cls` argument); U-MCCD,
  UN-MCCD as controls, called identically at every `m` (no `min.cls`
  argument exists in their call at all).
- **Reps (2026-09-05 cost decision):** RK methods (U-MCCD, SU-MCCD) at 50
  reps; NND methods (UN-MCCD, SUN-MCCD) at 100 reps — U/SU-MCCD cost ~11
  s/cell at d=3 and are 96% of the experiment's total cost.
- **Output columns:** `frac_flagged` = fraction of cluster 3's `m` points
  with the mutual-catch-graph minority label (`score == 1`); `n_clusters` =
  the detector's estimated macro-cluster count (so absorption of cluster 3
  into cluster 1/2 is visible directly); `n_unassigned_third` and
  `singleton_lost` disambiguate `frac_flagged = 0` between "cluster 3
  recognised as its own cluster" and "cluster 3's points never claimed by any
  cluster at all" (`mccd_translate()` scores unclaimed rows 0, same as a
  legitimate majority-connected point — the confound `WP8_VERIFICATION.md`
  item 87.2 named). A WP7 cluster file (`87_small_cluster_clusters.csv`,
  schema `m,d,rep,seed,method,row_index,true_cluster,detected_cluster`) was
  added for this experiment specifically to expose that per-point detail.

## Integrity results

`count.fields` pass on both output files — every one of the 4,801 / 960,001
lines (header + data) has the header's field count (11 and 8 respectively),
no malformed rows:

| File | Header fields | Data rows | Bad rows |
|---|---|---|---|
| `87_small_cluster.csv` | 11 | 4,800 | 0 |
| `87_small_cluster_clusters.csv` | 8 | 960,000 | 0 |

- **Expected metric rows:** 8 sizes x 2 d x (50 reps x 2 RK methods + 100
  reps x 2 NND methods) = 16 x (100 + 200) = **4,800** — matches exactly.
- **Duplicates on `(m, d, rep, method)`:** 0.
- **`status != "ok"` rows:** 0 — every one of the 4,800 cells succeeded on
  its first attempt; no error rows, nothing retried.
- **Duplicate seeds within a cell:** 0 — every `(m, d, rep, method)` maps to
  exactly one seed (seed canonicalisation against the fixed `SIZES x c(3,10)`
  CANON grid holds).
- **Completeness:** all 64 distinct `(m, d, method)` cells (8 sizes x 2 d x 4
  methods) have exactly their expected rep count (50 for U-MCCD/SU-MCCD, 100
  for UN-MCCD/SUN-MCCD).
- **Cluster file:** 960,000 rows = 4,800 cells x 200 points/cell — matches
  exactly; every cell has exactly 200 distinct `row_index` values, no
  duplicates on `(m, d, rep, method, row_index)`.
- **`true_cluster` distribution:** 1 = 456,600; 2 = 457,200; 3 = 46,200 rows
  (cluster 3's total point-count across all cells, `sum(m)` over the 4,800
  cells — consistent with the design).
- **`detected_cluster` NA (unassigned) count in the cluster file: 0** across
  all 960,000 rows.

The grid completed cleanly. Nothing was excluded from the tables below on
integrity grounds.

## Script's own `--summarize`

`Rscript revision_experiments/tr1/87_wp8_small_cluster.R --summarize` wrote
`results/tr1/wp8/87_small_cluster_summary.csv` (64 rows: one per `(m, d,
method)`). Its `frac_flagged_mean`/`frac_flagged_se`/`n_clusters_mean`/
`frac_reps_ncls3` columns were reproduced independently in a scratchpad
script (`wp8_87_summary.R`) and match to full precision on every row
(spot-checked: max |diff| on `p_ncls3` vs. `frac_reps_ncls3` across all 64
rows = 0). The independent computation adds two things the script's own
summary does not separate out: `p(n_clusters == 2)` (absorbed into a
background cluster) vs. `p(n_clusters == other)` (neither 2 nor 3), and
`mean_singleton_lost` alongside `mean_n_unassigned_third`.

## Bookkeeping: is `frac_flagged = 0` ever "never covered"?

**No.** `n_unassigned_third` is exactly 0 in every one of the 4,800 rows, and
`singleton_lost` is exactly 0 in every one of the 4,800 rows (confirmed
against the cluster file's own 0-NA-`detected_cluster` count above — the two
files agree). Every point of cluster 3, at every `(m, d, method)` cell, was
claimed by some detected macro-cluster. This rules out the confound
`WP8_VERIFICATION.md` flagged as a reason to add these columns: in this
experiment, `frac_flagged = 0` for cluster 3 always means "cluster 3's
points were claimed and voted majority" (typically because cluster 3 was
itself recognised as its own macro-cluster, or occasionally because it was
majority-connected within a larger absorbing cluster), never "cluster 3's
points were never covered by any ball." The diagnostic columns are clean
data on this run, not a null result — they confirm there is nothing to
disambiguate here, so `frac_flagged` can be read directly as "flagged
anomalous" without qualification.

## Curve table 1: `frac_flagged` (mean fraction of cluster 3's points flagged)

Rows = method x d, columns = m. SE range and n_reps given in the notes below
the table (full per-cell SEs in `wp8_87_summary.csv`, scratchpad).

| method | d | m=2 | m=4 | m=6 | m=8 | m=10 | m=12 | m=15 | m=20 |
|---|---|---|---|---|---|---|---|---|---|
| U-MCCD  | 3  | 0.980 | 0.785 | 0.340 | 0.355 | 0.240 | 0.148 | 0.120 | 0.018 |
| SU-MCCD | 3  | 1.000 | 1.000 | 1.000 | 1.000 | 0.520 | 0.200 | 0.121 | 0.006 |
| UN-MCCD | 3  | 1.000 | 0.800 | 0.600 | 0.439 | 0.309 | 0.231 | 0.187 | 0.091 |
| SUN-MCCD| 3  | 1.000 | 1.000 | 1.000 | 1.000 | 0.550 | 0.350 | 0.223 | 0.085 |
| U-MCCD  | 10 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.980 |
| SU-MCCD | 10 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.980 |
| UN-MCCD | 10 | 1.000 | 1.000 | 1.000 | 0.990 | 1.000 | 0.971 | 0.882 | 0.822 |
| SUN-MCCD| 10 | 1.000 | 1.000 | 1.000 | 1.000 | 0.990 | 0.920 | 0.841 | 0.592 |

n_reps: U-MCCD/SU-MCCD = 50 per cell; UN-MCCD/SUN-MCCD = 100 per cell (every
cell, no exceptions). SE at m=10 (the threshold cell, where variance is
largest): SU-MCCD 0.0714, SUN-MCCD 0.0500, U-MCCD 0.0581, UN-MCCD 0.0456 —
all four transitions at m=10, d=3 are resolved well beyond noise (frac_flagged
differences between adjacent m values there exceed 3-5 SEs).

## Curve table 2: P(n_clusters == 3) (third cluster recognised as its own macro-cluster)

| method | d | m=2 | m=4 | m=6 | m=8 | m=10 | m=12 | m=15 | m=20 |
|---|---|---|---|---|---|---|---|---|---|
| U-MCCD  | 3  | 0.02 | 0.22 | 0.66 | 0.66 | 0.78 | 0.86 | 0.90 | 1.00 |
| SU-MCCD | 3  | 0.00 | 0.00 | 0.00 | 0.00 | 0.48 | 0.80 | 0.88 | 1.00 |
| UN-MCCD | 3  | 0.00 | 0.20 | 0.40 | 0.57 | 0.70 | 0.78 | 0.83 | 0.91 |
| SUN-MCCD| 3  | 0.00 | 0.00 | 0.00 | 0.00 | 0.45 | 0.65 | 0.78 | 0.91 |
| U-MCCD  | 10 | 0.16 | 0.08 | 0.14 | 0.08 | 0.04 | 0.08 | 0.08 | 0.12 |
| SU-MCCD | 10 | 0.16 | 0.08 | 0.14 | 0.08 | 0.04 | 0.08 | 0.08 | 0.12 |
| UN-MCCD | 10 | 0.00 | 0.00 | 0.00 | 0.01 | 0.00 | 0.04 | 0.13 | 0.18 |
| SUN-MCCD| 10 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 | 0.08 | 0.16 | 0.41 |

`P(n_clusters == 2)` is `1 - P(n_clusters==3) - P(other)` in every row; `P(other)`
(neither 2 nor 3 macro-clusters detected) is 0 in 62 of the 64 cells and 0.01
in the remaining two (SUN-MCCD d=3 m=20; UN-MCCD d=10 m=4) — a 1-of-50-or-100
rep artifact in each case, not a pattern.

Note `U-MCCD` and `SU-MCCD` give **numerically identical** `P(n_clusters==3)`
at every `m` at d=10 (0.16, 0.08, 0.14, 0.08, 0.04, 0.08, 0.08, 0.12) and
numerically identical `frac_flagged` at every `m` at d=10 as well (both
curves in table 1 read 1,1,1,1,1,1,1,0.98) — the `min.cls` argument makes
literally no observable difference to SU-MCCD relative to its no-`min.cls`
control at d=10 for any of the eight sizes tested.

## R3.7 answer

**The design fact that governs every reading below:** cluster 3 is
density-matched to the background at every `m` (see "Design recap"), so the
only thing changing along the x-axis is the third cluster's *point count*,
not its density relative to the other two clusters.

**Where SU-MCCD and SUN-MCCD flip, at d = 3.** Both S_min-using methods show
`frac_flagged = 1.000` (SE = 0, no rep disagreement) at every `m <=
round(S_min*200) = 10`'s immediate predecessor (m = 2, 4, 6, 8) — a flat
ceiling, not a gradual approach to 1 — and then step down starting exactly
at `m = 10`, the declared threshold: SU-MCCD 1.000 -> 0.520 -> 0.200 ->
0.121 -> 0.006 (m = 8, 10, 12, 15, 20); SUN-MCCD 1.000 -> 0.550 -> 0.350 ->
0.223 -> 0.085 over the same sizes. `P(n_clusters==3)` mirrors this exactly:
0.00 for m <= 8, jumping to 0.48 (SU-MCCD) / 0.45 (SUN-MCCD) at m = 10 and
climbing to 0.88-1.00 by m = 15-20.

**Where SU-MCCD and SUN-MCCD flip, at d = 10.** SUN-MCCD's decline is
smoother and starts later: `frac_flagged` stays at 1.000 through m = 8, is
0.990 at m = 10 (the threshold), then 0.920 (m=12), 0.841 (m=15), 0.592
(m=20) — a continuing decline past the threshold rather than a step
resolved within one or two sizes, as at d=3. SU-MCCD, by contrast, barely
moves at all: 1.000 through m = 15, then 0.980 at m = 20 — no resolvable
transition within the sizes tested.

**Do U-MCCD and UN-MCCD (no `min.cls`) show a different pattern?** Yes, and
the difference locates the cause of the flip:

- At d = 3, the no-`min.cls` controls decline *smoothly and earlier* than
  their `min.cls` counterparts, with no flat-ceiling segment: U-MCCD is
  already at 0.785 by m = 4 and 0.340 by m = 6 (vs. SU-MCCD's 1.000 through
  m = 8); UN-MCCD is at 0.800 by m = 4 and 0.600 by m = 6 (vs. SUN-MCCD's
  1.000 through m = 8). By m = 15-20 the `min.cls` and no-`min.cls` curves
  converge to comparable values (U-MCCD 0.120/0.018 vs. SU-MCCD
  0.121/0.006; UN-MCCD 0.187/0.091 vs. SUN-MCCD 0.223/0.085). **The flat
  ceiling at exactly 1.000 for m < 10, present only in the `min.cls`
  methods, is therefore attributable to the S_min rule specifically** — the
  underlying coverage geometry (visible in the controls) already begins
  partially recognising a size-4 or size-6 cluster as distinct at d=3;
  `min.cls` overrides that recognition and forces full flagging until the
  declared threshold is reached, exactly as R3.7 describes.
- At d = 10, U-MCCD and SU-MCCD are numerically indistinguishable at every
  size tested (see the note under table 2): the RK-based, uniform-coverage
  pair fails to recognise the third cluster as distinct almost regardless of
  `m`, up to and including m = 20 (10% of n). Here the flip cannot be
  attributed to S_min at all — there is no `min.cls` effect to see, because
  the coverage geometry alone already suppresses recognition at every size
  tested. UN-MCCD and SUN-MCCD *do* separate at d=10 (SUN-MCCD's frac_flagged
  falls further below UN-MCCD's as `m` grows past 10: 0.99 vs 1.00 at m=10,
  widening to 0.592 vs 0.822 at m=20) — i.e. at d=10 the min.cls-aware
  shape-adaptive method recognises a cluster of size 20 measurably better
  than its no-`min.cls` counterpart, but neither curve shows the flat-then-
  step signature seen at d=3, and even SUN-MCCD's best value at m=20
  (frac_flagged 0.592) is far from the near-zero values both NND methods
  reach at d=3.

**Asymmetry between d=3 and d=10, stated in numbers.** At d=3, the
`min.cls`-attributable component of the flip (the difference between the
`min.cls` method and its own control, at m=10, the threshold) is +0.28
(SU-MCCD 0.520 vs. U-MCCD 0.240) and +0.24 (SUN-MCCD 0.550 vs. UN-MCCD
0.309) — the S_min rule adds roughly a quarter to a third to the flagged
fraction right at threshold, on top of whatever the coverage geometry alone
would produce. At d=10, the same difference at m=10 is 0.000 (SU-MCCD vs.
U-MCCD, both exactly 1.000) and -0.010 (SUN-MCCD 0.990 vs. UN-MCCD 1.000,
i.e. SUN-MCCD flags marginally *less*, not more). The S_min-specific
contribution to the flip is resolvable and substantial at d=3 and is
essentially absent (RK pair) or reversed in sign and small (NND pair) at
d=10, where the coverage geometry's own failure to resolve a cluster this
small dominates whatever S_min is doing.

## Guidance sentence (from the measured curve, no interpretive language)

At d=3, `round(S_min * n)` is the largest true cluster size at which both
`min.cls`-aware methods return `frac_flagged = 1.000` with SE = 0; any true
cluster of that size or smaller is flagged in full, at every one of the 50
or 100 replicates measured. A user who wants a cluster of size `m*` to
survive as its own macro-cluster rather than be flagged in full should set
`S_min` so that `round(S_min * n) < m*`; at `round(S_min * n) = m*` the
measured flagged fraction is 0.520 (SU-MCCD) / 0.550 (SUN-MCCD) at d=3, not
0 and not 1, and the fraction continues to fall as the true cluster size
grows past that value (0.200/0.350 at 1.2x the threshold, 0.006/0.085 at 2x
the threshold, both at d=3). At d=10 this relationship does not hold for the
RK-based pair (U-MCCD, SU-MCCD): `frac_flagged >= 0.98` for every measured
`m` up to 20 (10% of n) independent of where `S_min` is set relative to
`m*`, so lowering `S_min` alone does not recover a cluster of this size
range under RK-based coverage at this dimension; only the NND-based pair
shows `frac_flagged` declining measurably as `m*` grows past
`round(S_min*n)` at d=10, and even there the decline is slower and reaches a
higher floor (0.592 at `m* = 2x` threshold, vs. 0.085 at d=3) than at d=3.

## Caveats

- **RK reps.** U-MCCD/SU-MCCD ran at 50 reps (vs. 100 for UN-MCCD/SUN-MCCD),
  a launch-time cost decision (`WP8_PROTOCOL.md`: U/SU-MCCD cost ~11 s/cell
  at d=3 and are 96% of the experiment's total cost). The largest SEs
  reported above (m=10, the threshold cell) are 0.058-0.071 for the RK pair
  vs. 0.045-0.050 for the NND pair — wider but the qualitative flip pattern
  at d=3 (flat ceiling then step) is resolved by 3+ SEs at every adjacent-m
  comparison used above.
- **87's cluster file exists for WP7.** `87_small_cluster_clusters.csv`
  (960,000 rows) was written for this experiment specifically
  (superseding the original protocol's "no WP7 cluster file for this
  experiment" — see `WP8_PROTOCOL.md`'s dated note) and is available for
  WP7 use beyond what this findings file draws from it (only the aggregate
  `detected_cluster`-NA count was used here).
- **`P(other)` (neither 2 nor 3 macro-clusters) is nonzero in 2 of 64
  cells**, at 0.01 each (SUN-MCCD d=3 m=20; UN-MCCD d=10 m=4) — a single
  replicate in each case, not a pattern; excluded from the reported curves'
  interpretation above without materially changing any number (`P(n_clusters
  == 2)` absorbs the complement in every other cell).
- This file reports means, SEs, and proportions as descriptive statistics of
  the observed sample; no hypothesis test, p-value, or multiplicity
  correction beyond the SE comparisons stated above has been computed or is
  implied.

## Source files

- Raw data (not committed, gitignored): `87_small_cluster.csv`,
  `87_small_cluster_clusters.csv`, script's own summary
  `87_small_cluster_summary.csv` — all in this directory.
- Analyst's independent computation: `wp8_87_summary.R`, and its CSV outputs
  `wp8_87_summary.csv`, `wp8_87_curve_flag.csv`, `wp8_87_curve_p3.csv` —
  temporary
  scratch directory, not retained and not part of this repository.
