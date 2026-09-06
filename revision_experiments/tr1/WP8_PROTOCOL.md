# WP8 protocol — hard cases and statistical validity

**Declared 2026-09-05, before any of the four experiments below was run.**
Follows the declare-before-look discipline of `BENCHMARK_EXPANSION_RULE.md` and
`WP4_PROTOCOL.md`: the settings, generators, seeds and outputs below are fixed
in advance; any deviation forced by an implementation failure is recorded as
an appended, dated note rather than an edit to the declarations.

WP8 answers R1.1, R3.4, R3.5 and R3.7 with four targeted synthetic
experiments, all at **100 replicates, d ∈ {3, 10}, n = 200, S_min = 0.05**
(the single constant now used everywhere in this revision — see `CLAUDE.md`),
α from the paper's own dimension-schedule resolvers
(`rk_quant_label_paper`, `nn_quant_label_paper_UN`,
`nn_quant_label_paper_SUN` in `shared/harness.R`), and all four MCCD
detectors (U-MCCD, SU-MCCD, UN-MCCD, SUN-MCCD) unless an experiment restricts
the method list for a stated reason.

**Default method list (experiments 1, 2, 4):** the four MCCD detectors plus
the five published baselines — `U-MCCD, SU-MCCD, UN-MCCD, SUN-MCCD, LOF,
DBSCAN, MST, ODIN, iForest` — called through `METHOD_REGISTRY` exactly as WP0
wired it. The four generic `RKCCD-OOS/IOS`, `UNCCD-OOS/IOS` registry entries
are not part of this list; they are components of the paper's own methods,
not a sixth family to compare.

**Experiment 3 uses SU-MCCD and SUN-MCCD only, with U-MCCD and UN-MCCD as
controls** (per the task spec) — those two take no `min.cls` argument at all
in `wp0_mccd_methods.R`, so S_min cannot act on them; running them anyway
shows what happens to a small legitimate cluster when nothing enforces a
minimum size.

## Shared conventions

- **Seeding.** Every script defines `BASE_SEED` (85 → 8501, 86 → 8601,
  87 → 8701, 88 → 8801) and computes, for setting index `s` (1-based, in the
  script's own settings table) and replicate `rep`:
  `seed <- BASE_SEED + 100000L * s + rep`. The realised seed is written to a
  `seed` column on every output row, so any cell can be regenerated exactly
  from the CSV alone.
- **α and S_min.** α is never overridden — every method call lets
  `wp0_mccd_methods.R`'s wrappers resolve `quant = NULL` via the paper
  resolvers for the cell's own `d`. S_min = 0.05 is passed as `min.cls` to
  SU-MCCD and SUN-MCCD only (their wrapper's parameter name), everywhere it
  is used.
- **Per-cell append / restart.** Every script appends one (or a
  fixed, declared, atomic group of) row(s) per finished cell via
  `harness.R`'s `append_result()`, gated by `has_result()` so a restart skips
  finished cells. G: drops off the bus under sustained writes (`CLAUDE.md`),
  so nothing is buffered across cells.
  **A `note` column pitfall, found while smoke-testing 86/87/88 and fixed
  before any real cell was run:** `has_result()` treats any NA in a
  non-key column as a partial/incomplete row (its own truncation guard),
  so a successful row's `note` field cannot be `NA_character_` without
  defeating skip-on-restart for every cell. The obvious fix, an empty
  string, does not work either — a CSV column that is `""` on every row
  round-trips through `read.csv()`'s automatic `type.convert()` as
  **logical `NA`** (an all-empty character column has nothing to infer a
  type from). All four scripts use the non-empty placeholder `"-"` for a
  successful cell's `note`, reserving real messages for `status="error"`
  rows (which are *meant* to stay retriable, since their numeric payload
  columns are genuinely `NA`).
- **`--smoke`.** Runs exactly one replicate of one setting (one
  `(generator/type, d)` combination — not one method; the whole default
  method list still runs, because the point of the smoke test is to prove
  every method's translation path, not just one) and writes to
  `results/tr1/wp8/smoke/<script>.csv`, never the production file.
- **Thresholds.** `REAL_DATA_THRESHOLDS` (from `harness.R`, extended by
  `wp0_mccd_methods.R` to `WP0_THRESHOLDS`) supplies the decision threshold
  for every method: 0.5 for the four MCCD detectors and for
  DBSCAN/MST/ODIN's binarized scores, 1.5 for LOF, 0.55 for iForest.
- **CCD generator primitives.** All generators use `rpoisball.unit(n, d)`
  (uniform draw inside the unit d-ball; sourced transitively into
  `.GlobalEnv` by `harness.R`'s `RKCCD_OOS_IOS.R`/`UNCCD_OOS_IOS.R` chain —
  no separate `source()` needed) and `MASS::mvrnorm`, exactly as
  `55_wp2c_simulation_arm.R`'s `gen_uniform`/`gen_gaussian` already do for
  the two-cluster base case. Two-cluster geometry constants (`cls_dis = 3`,
  `otl_dis = 2`, `r_min/r_max = 0.7/1.3`) are inherited unchanged from that
  script (itself copied verbatim from the original simulation drivers) so the
  hard-case generators sit in the same coordinate scale as the paper's own
  simulations, not an arbitrary new one.

## Experiment 1 — boundary false positives (85, R1.1)

**Settings:** generator ∈ {uniform, gaussian} × d ∈ {3, 10}, n = 200,
contamination = 0.05 (the "usual" level; chosen to match S_min = 0.05
numerically only by coincidence — the two parameters are unrelated). 4
settings, 100 reps each.

**Generators.** `gen_uniform_boundary`/`gen_gaussian_boundary`: two clusters
of sizes `n1 = round(n*(1-cont)*0.5)`, `n2 = round(n*(1-cont)*0.5) - 1`
(the same off-by-one the original drivers carry — n1+n2+n0 lands at 199, not
200, for n=200/cont=0.05; not corrected here, since correcting it would be a
math change to a component the protection rule reserves for explicit
sign-off) centred at `mu1 = rep(3,d)`, `mu2 = c(3+cls_dis, rep(3,d-1))`,
points drawn `rpoisball.unit(n_k, d) * runif(1, r_min, r_max) + mu_k`;
`n0 = round(n*cont)` outliers drawn `rpoisball.unit(1,d)*5 + mean(mu1,mu2)`,
rejected and redrawn until `> otl_dis` from both cluster centres — identical
to `55_wp2c_simulation_arm.R`'s `gen_uniform`/`gen_gaussian`. Gaussian variant
replaces the ball draw with `mvrnorm(n_k, mu_k, diag(d)*(sigma*runif(1,
r_min,r_max))^2)`, `sigma = 1/sqrt(qchisq(1-0.01, d))`.

**Per-point normalized radius.** For every regular point `i` in cluster `k`,
`r_i = ||x_i - mu_k||`; `r_max_k` = the max `r_i` over that cluster in that
replicate; `bin_i = clamp(ceil(10 * r_i / r_max_k), 1, 10)`. The two
clusters' regular points are pooled into the same 10 bins (the boundary
question is about distance-from-own-centre, not which cluster).

**Output** (`results/tr1/wp8/85_boundary_fp.csv`): one row per
`(setting_id, d, rep, seed, method, bin)`: `n_regular_in_bin`,
`n_flagged_in_bin`. **Atomicity note:** 10 bin-rows share one detector call;
a companion `85_boundary_fp_done.csv` (keys `setting_id, d, rep, method`) is
appended only after all 10 bin rows for that cell have been written, and is
what `has_result()` checks for skip/restart — so a rep is never
half-recorded across a restart even though the 10 bins are not a single
physical row.

**Also written** (WP7): `results/tr1/wp8/85_boundary_fp_clusters.csv`, one
row per regular point per MCCD-family method (`U-MCCD, SU-MCCD, UN-MCCD,
SUN-MCCD` only — baselines expose no macro-cluster id through the harness):
`setting_id, d, rep, seed, method, row_index, true_cluster, detected_cluster`.
Written right after the bin rows, before the done-marker for that cell.

## Experiment 2 — local / bridged / collective outliers (86, R3.7)

**Settings:** type ∈ {local, bridge, collective} × d ∈ {3, 10}, n = 200,
n0 = round(n*0.05) = 10 (same "usual" 5% used in experiment 1, held fixed
across all three types so the comparison is about *placement*, not
*count*). 6 settings, 100 reps each. Two base clusters as in experiment 1
(`n1 = 95, n2 = 94` at cont = 0.05).

**Generators** (each replaces only the *outlier* placement; the two base
clusters are always `rpoisball.unit`-drawn as above):

- **local** — an outlier that is deep inside a cluster's spatial extent but
  locally isolated. For each of the `n0` outliers: pick a host cluster `k`
  (first `ceiling(n0/2)` assigned to cluster 1, the rest to cluster 2, in
  order); draw a candidate centre `g = mu_k + rpoisball.unit(1,d) * 0.4`
  (interior offset, well inside the typical 0.7–1.3 cluster radius); the
  outlier point itself **is** `g`. The two clusters' regular points are then
  generated by rejection sampling — draw a candidate the usual way, redraw
  (up to 200 attempts, after which the candidate is kept anyway and the
  event logged) if it falls within clearance radius 0.3 of any `g` assigned
  to that cluster — so no regular point sits within 0.3 of the outlier
  (`R/ccds/Kest.R`'s `rpoisball.unit` returns bulk mass near the ball's edge
  in high `d`, so 200 attempts is generous, not tight, at d=3/10).
- **bridge** — a chain of `n0` points on the segment between the two cluster
  centres: `bridge_i = mu1 + (i/(n0+1))*(mu2-mu1) + rnorm(d, 0, 0.15)` for
  `i = 1..n0`, each with independent N(0, 0.15²) jitter in every coordinate.
  Regular clusters are the plain `rpoisball.unit` draw (no rejection).
- **collective** — a tight minority group just outside cluster 1's boundary:
  draw one random unit direction `u` (a single N(0,I_d) vector, normalized);
  `center = mu1 + u * (r_max + 0.5)` with `r_max = 1.3` (the upper end of the
  cluster's own radius jitter, so the group sits just past the cluster's
  typical edge); the `n0` outliers are `center + rnorm(d, 0, 0.15)`
  (independent per point, same tight spread as the bridge jitter so the
  group is internally well-connected). Regular clusters are the plain draw.

**Output** (`results/tr1/wp8/86_outlier_types.csv`): one row per
`(type, d, rep, seed, method)`: `TPR, TNR, BA, F2` via `evaluate()`
(n0 = 10 > 0 throughout, so `evaluate()`/`count_scores2` is safe here),
`n_flagged`. Single row per cell — `has_result()` on the main file gates
skip/restart directly, no separate done-marker needed.

**Also written** (WP7): `results/tr1/wp8/86_outlier_types_clusters.csv`,
same schema as experiment 1's cluster file, `true_cluster ∈ {1,2}` for
regular points (outliers excluded — they have no true cluster). Written
immediately after each cell's metrics row; if a run is interrupted between
the two writes for one cell, that one cell's cluster rows are permanently
missing (the metrics row alone still gates skip) — accepted as a known gap,
not covered by a second done-marker, since the metrics file is what every
other WP consumes and the cluster file is WP7-only supplementary data.

## Experiment 3 — small legitimate cluster below S_min (87, R3.7)

**Settings:** third-cluster size `m ∈ {2, 4, 6, 8, 10, 12, 15, 20}` × d ∈
{3, 10}, n = 200 total (fixed), S_min = 0.05 so `round(0.05*200) = 10` is the
size at which the flip is expected. 16 settings, 100 reps each.

**Generator.** Three cluster centres, `mu1 = rep(3,d)`,
`mu2 = mu1 + 3*e1`, `mu3 = mu1 + 3*e2` (`e1, e2` the first two standard basis
vectors — needs `d ≥ 2`, satisfied at d = 3 and d = 10), same `cls_dis = 3`
separation as every other WP8 generator. Background split
`n_bg = 200 - m`, `n1 = floor(n_bg/2)`, `n2 = n_bg - n1`; all three clusters
drawn `rpoisball.unit(n_k,d) * runif(1,0.7,1.3) + mu_k`, independently per
cluster (three independent radius-jitter draws, not one shared draw). **No
separate contamination outliers** — n0 = 0 throughout; the only question
this experiment asks is what happens to cluster 3.

**Methods:** SU-MCCD, SUN-MCCD (`min.cls = 0.05`) as the focal methods;
U-MCCD, UN-MCCD (no `min.cls` argument in their wrapper) as controls, run
identically at every `m` since nothing in their call changes with `m`.

**Output** (`results/tr1/wp8/87_small_cluster.csv`): one row per
`(m, d, rep, seed, method)`: `frac_flagged` = fraction of cluster 3's `m`
points with `score == 1` (the mutual-catch-graph minority label from
`mccd_translate()`, not a labelled-outlier rate — there is no ground-truth
"outlier" here, only ground-truth cluster membership), plus `n_clusters`
(the detector's estimated macro-cluster count, so a collapse of cluster 3
into cluster 1 or 2 is visible directly, not just inferred from
`frac_flagged`). Single row per cell, gated directly by `has_result()`. No
WP7 cluster file for this experiment (there are no contamination outliers to
contrast against, and the per-point cluster/label detail is already fully
captured by `frac_flagged` and `n_clusters`).

## Experiment 4 — false positives under violated CSR, no outliers (88, R3.4)

**Settings:** generator ∈ {uniform_control, beta_gradient, gauss_unequal} ×
d ∈ {3, 10}, n = 200, **n0 = 0 throughout** (`Y = rep(1, n)`). 6 settings,
100 reps each.

**Generators** (single cluster, no macro-structure, no outliers):

- **uniform_control** — `X[,j] ~ iid Uniform(0,1)`, `j=1..d`. The genuine
  null the RK/NN spatial-randomness test is calibrated against; the control
  every flag-rate number is read relative to.
- **beta_gradient** — `X[,j] ~ iid Beta(2,5)`, `j=1..d`, support already
  `[0,1]` so no rescaling. Beta(2,5) has mode ≈ 0.2 and mean 2/7 ≈ 0.286: a
  density that decays smoothly away from one corner of the box, i.e. a
  genuine density gradient rather than a second population.
- **gauss_unequal** — `X ~ MVN(0, diag(sigma_1^2,...,sigma_d^2))`,
  `sigma_j` linearly spaced `seq(0.3, 1.5, length.out = d)`: an axis-aligned
  ellipsoid, anisotropic but still a single homogeneous population — tests
  the *shape* mismatch (non-spherical density) separately from the *gradient*
  mismatch above.

**Why n0 = 0 changes the metric.** `count_scores2` (via `evaluate()`) indexes
`label_pred[(n-n0+1):n]` for the TPR numerator; at `n0=0` that is
`(n+1):n`, a reversed two-element sequence that reads past the vector's end.
`evaluate()` is therefore **not called** in this experiment. The reported
quantity is the **raw flag rate**, `mean(score >= threshold)` over all n
points, using the same `ALL_THRESHOLDS` (merged `REAL_DATA_THRESHOLDS` +
`WP0_THRESHOLDS`) as every other WP8 experiment. TNR = 1 − flag_rate exactly
here (no true outliers to miss), so flag rate is the whole story.

**Output** (`results/tr1/wp8/88_csr_violation.csv`): one row per
`(generator, d, rep, seed, method)`: `flag_rate`, `n_flagged`. Single row
per cell, gated directly by `has_result()`.

**Also written** (WP7): `results/tr1/wp8/88_csr_violation_clusters.csv`,
same schema and same "metrics row gates, cluster file may lag on
interruption" caveat as experiment 2 — here `true_cluster` is always 1 (a
single population), so the WP7 use is k-hat accuracy (is a homogeneous but
non-uniform or anisotropic population ever split into >1 macro-cluster) more
than ARI/NMI, which are degenerate with one true cluster.

## Outputs summary

```
results/tr1/wp8/85_boundary_fp.csv                main: setting_id,d,rep,seed,method,bin,n_regular_in_bin,n_flagged_in_bin
results/tr1/wp8/85_boundary_fp_done.csv           done-markers: setting_id,d,rep,method
results/tr1/wp8/85_boundary_fp_clusters.csv       WP7: setting_id,d,rep,seed,method,row_index,true_cluster,detected_cluster
results/tr1/wp8/86_outlier_types.csv              main: type,d,rep,seed,method,TPR,TNR,BA,F2,n_flagged
results/tr1/wp8/86_outlier_types_clusters.csv     WP7
results/tr1/wp8/87_small_cluster.csv              main: m,d,rep,seed,method,frac_flagged,n_clusters
results/tr1/wp8/88_csr_violation.csv              main: generator,d,rep,seed,method,flag_rate,n_flagged
results/tr1/wp8/88_csr_violation_clusters.csv     WP7
results/tr1/wp8/smoke/<script-basename>.csv       --smoke output, same schema as the main file it stands in for
```

## Dated appended notes (2026-09-05, post Opus verification WP8_VERIFICATION.md)

**86 FAIL fix -- `local` and `bridge` generators rewritten; mandatory
acceptance check.** The originally declared `local` (fixed 0.4 interior
offset / 0.3 clearance) and `bridge` (plain centre-to-centre chain)
constructions did not produce the named phenomena: measured local
outlier-NN/regular-NN ratio was 1.28 (d=3) and 0.61 (d=10); measured bridge
had 5.2-6.2 of 10 points per rep landing inside a realised cluster ball.
Replacement constructions (implemented in `86_wp8_outlier_types.R`,
`draw_cluster()`/`gen_local()`/`gen_bridge()`/`gen_collective()`):

- `local`: host cluster regular points confined to a dense core of radius
  `core_frac * scale` with `scale` the cluster's own realised jitter draw
  (`runif(1, R_MIN, R_MAX)`); each outlier placed at a uniformly random
  direction from its host centre at radius `runif(1, shell_lo, shell_hi) *
  scale`. First attempt (`core_frac = 0.5`, `shell = (0.75, 1.0)`) measured
  ratio (mean of 20 reps, `mean(outlier NN dist)/median(host regular NN
  dist)`, seeds = the production `BASE_SEED + 100000*s_i + rep` scheme for
  `rep in 1:20`): **d=3 mean 4.122 (min 3.477, max 4.944) -- PASS**; **d=10
  mean 1.892 (min 1.740, max 2.178) -- FAIL (< 2)**. Second attempt
  (`core_frac = 0.4`, `shell = (0.8, 1.0)`, both within the "e.g." bracket
  named in WP8_VERIFICATION.md): **d=3 mean 6.032 (min 4.802, max
  7.523)**; **d=10 mean 2.540 (min 2.360, max 2.874) -- PASS at both d**.
  Adopted.
- `bridge`: `draw_cluster()` now exposes `attr(pts, "scale")`, the realised
  jitter draw; bridge points are placed only in the gap between the two
  clusters' realised surfaces, `t_lo = (s1+clearance)/CLS_DIS`,
  `t_hi = 1-(s2+clearance)/CLS_DIS`, `t_i = t_lo + (i-0.5)/n0*(t_hi-t_lo)`.
  First attempt (`clearance = 0.25`, as literally specified in
  WP8_VERIFICATION.md): 20-rep check gave **0 violations at d=10** but **1
  of 200 points (1 rep of 20) inside a realised ball at d=3**. Root cause:
  the per-point `rnorm(d, 0, 0.15)` jitter, added AFTER placing the point at
  the clearance-satisfying `t_i`, occasionally erodes the margin on its own.
  Second attempt (`clearance = 0.45`, no redraw): 20-rep check now **0/20 at
  both d=3 and d=10**, but a larger, more sensitive 100-rep check (not part
  of the mandatory bar, run for due diligence) found 6/1000 points (d=3) and
  3/1000 points (d=10) still inside. Third attempt (kept `clearance = 0.45`,
  added a bounded redraw of the jitter alone, up to 50 attempts, whenever the
  candidate lands inside either realised ball -- the mean position `t_i` is
  never moved): 100-rep check **0/1000 at both d=3 and d=10**. Adopted. The
  `n_inside` column (sentinel `-1` for `local`/`collective`, where it does
  not apply) is retained in the main CSV as a per-row diagnostic in case a
  future, larger run still exhibits a rare residual case.
- `collective`: standoff corrected to `mu1 + u * (s1 + 0.5)` with `s1` the
  realised jitter draw (was the fixed nominal `r_max = 1.3`); this replaces
  the protocol's original "just past the cluster's typical edge" language
  above with the measured, realised-radius-relative stand-off.

**86 cross-cutting.** Seed canonicalisation (`CANON` built from the fixed
`type x d in {3,10}` grid, `s_i <- match(setting_id, CANON$setting_id)`,
independent of any settings/dims subset passed on the command line); a
`run [settings] [reps] [budget] [methods]` mode (settings: comma list of
`CANON$setting_id` or `ALL`); cluster rows written as one block per cell via
`append_result()`'s multi-row-list form; a `--summarize` mode (means and SEs
of TPR/TNR/BA/F2 per `type x d x method`, `status=="ok"` filter, de-dup on
`(setting_id, d, rep, method)` keeping the last row). `--smoke` now covers
all 6 canonical settings (both d, all three types), one rep, all 9 default
methods -- confirmed non-zero TPR for LOF on `local` (TPR=1.0 at both d=3
and d=10) and done-skip on re-invocation.

**85 fixes.** Seed canonicalisation and `run [settings] [reps] [budget]
[methods]` mode as in 86, against the fixed `generator x d in {3,10}` CANON.
Added `stopifnot(length(res$score) == dat$n, !anyNA(res$score))` after the
detector call (85 was the only WP8 script without this guard). Added
equal-mass binning: main CSV gained a `bin_type` column (`"width"` = the
original `ceiling(10*r/r_max)`; `"mass"` =
`ceiling(10*rank(r,ties.method="first")/length(r))`, computed per cluster
like `"width"`) and now writes 20 bin-rows per cell instead of 10. Both the
20 bin-rows and (for MCCD methods) the up-to-189 cluster-rows are written as
ONE `append_result()` block each. `unassigned_rows` (from
`res$unassigned_rows`, sentinel `-1L` for the five baselines, which have no
such field -- never `NA_integer_`, which would trip `has_result()`'s
partial-row guard on the done file) is added as a payload column on the
done row. Added `--summarize` (per `setting_id x d x method x bin_type x
bin`: `pooled_ratio = sum(n_flagged)/sum(n_regular)` as the point estimate,
plus `mean_per_rep_ratio`/`se_per_rep_ratio` over reps with a non-empty
bin, joined against the done file so a cell that errored out is excluded).
**Bug found and fixed while smoke-testing `--summarize`:** the de-dup key
used to drop retried-cell duplicates was `(setting_id,d,rep,method)` --
correct for the other three scripts' one-row-per-cell files, but 85 writes
20 rows per cell, so that key collapsed all 20 legitimate bin rows down to
1 (whichever bin happened to be written last), turning a 180-row smoke
summary into 9 rows. Fixed by extending the key to
`(setting_id,d,rep,method,bin_type,bin)`; the done file's own de-dup key
correctly stays `(setting_id,d,rep,method)` (one row per cell there).
Smoke-verified: 180-row summary from a single-cell (9-method) smoke run,
mass bins visibly near-equal-count (18-20 of ~189 points each) vs width
bins piling into bin 10, done-file skip confirmed on re-invocation.

**87 fixes.** Seed canonicalisation against the fixed `SIZES x c(3,10)`
CANON grid (`setting_id = "m{m}_d{d}"`) and a
`run [settings] [reps] [budget] [msel]` mode (`settings`: comma list of
CANON setting_id or "ALL"; `msel`: an additional comma list of `m` values
that further restricts whichever settings were already selected, for quick
chunking by size without composing full setting_id strings).
**Density-matched cluster 3** (deviation from the original declaration,
which gave cluster 3 the same `runif(1,R_MIN,R_MAX)` radius as clusters
1-2): `data3 <- rpoisball.unit(m,d) * (runif(1,R_MIN,R_MAX) * (m/n1)^(1/d))
+ mu3` -- volume scales as radius^d, so scaling the radius by `(m/n1)^(1/d)`
holds cluster 3's point density (points per unit volume) equal to cluster
1's at every `m`, instead of confounding "smaller" with "sparser" (measured
NN-spacing ratio vs background at d=3 under the OLD generator: 4.42 at m=2,
2.26 at m=10, 1.71 at m=20 -- decreasing, i.e. genuinely getting sparser,
not just smaller). Added `n_unassigned_third` (`sum(is.na(res$cluster
[third_idx]))`) and `singleton_lost` (`res$singleton_lost_rows`, sentinel
`-1L` if absent -- never `NA_integer_`) columns, since `frac_flagged` alone
conflates "cluster 3 recognised as its own cluster but voted minority"
with "cluster 3 never claimed by any cluster at all"
(`mccd_translate()` scores unclaimed rows 0, same as a legitimate
majority-connected point). Added a WP7 cluster file,
`87_small_cluster_clusters.csv` (schema `m,d,rep,seed,method,row_index,
true_cluster,detected_cluster`, `true_cluster` in `{1,2,3}`, one block per
cell) -- this supersedes the original protocol's "No WP7 cluster file for
this experiment" (experiment 3's own per-point cluster/label detail is
exactly what WP8_VERIFICATION.md asked to expose). Added `--summarize`
(mean/SE of `frac_flagged` and `n_clusters` per `m x d x method`, plus
`frac_reps_ncls3` = the fraction of reps with `n_clusters == 3`).
**Cost decision:** RK methods (U-MCCD, SU-MCCD, ~11 s/cell at d=3, 96% of
the experiment's total cost) run 50 reps; NND methods (UN-MCCD, SUN-MCCD)
run 100 reps -- `REPS_BY_FAMILY` in the script, overridable uniformly via
the `reps` run-mode argument for manual chunking. Smoke-verified (m=10,
d=3, all 4 methods): correct columns on both the main and WP7 cluster
files, main-CSV skip confirmed on re-invocation, `--summarize` produces
one row per method with the expected fields.

**88 fixes.** Seed canonicalisation against the fixed
`generator x d in {3,10}` CANON (now 4 generators, 8 settings, see below)
and a `run [settings] [reps] [budget] [methods]` mode. Added `gauss_equal`
(`MASS::mvrnorm(n, rep(0,d), diag(d))`, sigma=1 in every coordinate) so
`gauss_unequal - gauss_equal` isolates ANISOTROPY specifically; the
original `gauss_unequal - uniform_control` contrast confounded anisotropy
with the fact that any Gaussian (equal-variance included) has a non-uniform
radial density, unlike the uniform control. DBSCAN's successful rows now
carry `note = "oracle contamination = 0; flag rate 0 by construction"`
instead of the generic `"-"` (it is fed `Y = rep(1,n)`, so its internal
oracle-contamination-based quantile is 0 and its flag rate is 0 by
construction, not by detection -- footnote this wherever DBSCAN's 88 numbers
are tabulated). Added an MST threshold sweep, `MST@1.05/1.40/1.60/2.00`
(the fixed `MST` row already covers 1.2, so it is not duplicated in the
sweep) -- `--summarize` reports, per `(generator, d)` with `generator !=
uniform_control`, the fixed-1.2 flag rate/delta alongside the sweep member
that minimizes `|delta_vs_control_mean|` (closest match to the CSR
control's own flag rate = least distorted by that generator's CSR
violation), labelled as "best" but never replacing the always-reported
1.2 row. Cluster rows written as one block per cell (200 rows, `true_cluster`
always 1). `--summarize` adds the paired `generator - uniform_control`
delta at the same `(d, rep)` with its own SE, alongside the per-`(generator,
d, method)` mean/SE of `flag_rate`.

**88 cost finding (2026-09-05, measured while debugging a stalled
foreground smoke run -- important for scheduling the full grid).** The
9-method-cell benchmark in WP8_VERIFICATION.md (23.3 s at d=3, U/SU-MCCD
11 s each) was measured on 85/86's TWO-CLUSTER settings and does NOT
transfer to 88's single-population generators. U-MCCD and SU-MCCD (both
RK-based; their density-calibration bracket search, `connected.ksccd.m`,
brackets the largest density keeping the core connected -- see `CLAUDE.md`)
are 8-140x slower here because a single homogeneous population gives the
bracket search no natural connectivity break to lock onto. Measured
one-core cost, n=200, d=3, one call each (U-MCCD; SU-MCCD tracked
separately and closely matched, both listed): `uniform_control` 137.9 s /
138.4 s; `beta_gradient` 34.1 s / 32.6 s; `gauss_equal` 43.2 s / 43.1 s;
`gauss_unequal` 8.6 s / 8.7 s. At d=10 all four generators dropped to
3.2-3.7 s for U-MCCD (SUN/UN-MCCD were already fast, ~1.3 s, at both d).
Consequence for the full-grid estimate: at d=3, the U-MCCD+SU-MCCD pair
alone costs approximately 138+34+43+9 ≈ 224 s PER METHOD across the 4
generators, i.e. roughly 448 s (~7.5 min) per rep at d=3 for the RK pair,
against d=10's well under 30 s/rep for the same pair -- d=3 dominates 88's
cost, the reverse of 85/86/87's own d=3-cheaper-than-d=10 pattern. See
"Revised full-grid cost estimate" below.

**Schema and chunking interface, consolidated (applies to 85, 86, 87, 88).**
Every script now: (1) builds a fixed `CANON` settings table from the full
declared grid at load time, independent of any CLI selection, and derives
each setting's `s_i` via `match(setting_id, CANON$setting_id)` before
computing `seed <- BASE_SEED + 100000*s_i + rep` -- a setting's seed is
therefore invariant to how a run is chunked; (2) exposes a
`run [settings] [reps] [budget] [methods]` CLI (positional; `settings` is a
comma list of `CANON$setting_id` or `"ALL"`; 87 additionally exposes a 4th
`msel` argument, a comma list of `m` values, since its `settings` selector
alone cannot express "every d for these sizes" as compactly); (3) writes
every per-cell block (bin rows, cluster rows) via one multi-row
`append_result()` call rather than a row-by-row loop; (4) implements
`--summarize`, filtering `status=="ok"` and de-duplicating retried cells on
the row's own natural key (see the 85 dated note above for the one
script -- 85 -- where that key must include the bin columns, not just the
cell-identifying columns); (5) guards its CLI-dispatch tail with
`if (!isTRUE(getOption("wp8.no_main", FALSE)))`, so the generator/settings
internals can be `source()`d (with that option set) for ad hoc checks --
e.g. the 86 acceptance check and the 88 cost finding above -- without
triggering a live run. **NA convention:** every payload column that
`has_result()` can see (i.e. every non-key column of a main-CSV row, a done
row, or -- since it is checked defensively even though nothing gates on it
-- a WP7 cluster row's `detected_cluster`) must never be `NA` on a
successful row, because `has_result()` treats ANY such `NA` as a
sign of a partial/truncated write and will re-run the cell on restart.
Two sentinels are used throughout instead: the string `"-"` for a `note`
field with nothing to say (an empty string round-trips through
`read.csv()`'s `type.convert()` as logical `NA` on an all-empty column,
which is exactly what trips the guard), and the integer `-1L` for a numeric
diagnostic that does not apply to a given row (86's `n_inside` for
`local`/`collective`; 85's `unassigned_rows` and 87's `singleton_lost` for
the five baseline methods, which have no cluster-assignment concept at
all). The one place `NA` is used deliberately and safely is
`detected_cluster` in the WP7 cluster files (`85_boundary_fp_clusters.csv`,
`86_outlier_types_clusters.csv`, `87_small_cluster_clusters.csv`,
`88_csr_violation_clusters.csv`) -- it means "this point was never claimed
by any cluster" (`mccd_translate()`'s unassigned bucket) and is safe there
because nothing calls `has_result()` against a cluster file; only the main
metrics file (and 85's done file) gate skip/restart.

## Revised full-grid cost estimate (2026-09-05, one core, post-fix)

WP8_VERIFICATION.md's own figures ("As written: 85 2.1h, 86 3.1h, 87 6.7h,
88 3.1h (15h); with block writes and the 87 decision above, roughly 11h")
assumed a uniform ~23.3 s/7.3 s (d=3/d=10) 9-method-cell cost measured on
TWO-CLUSTER settings. That transfers cleanly to 85 and 86 (same two-cluster
generators) and, combined with the 87 cost decision (RK methods 50 reps,
NND methods 100 reps, cutting 87's dominant RK cost roughly in half), gives:

- **85**: 4 settings x 100 reps, ~23.3 s/rep at d=3, ~7.3 s/rep at d=10 ->
  (2 x 100 x 23.3 + 2 x 100 x 7.3) s ~= 1.7 h. Block writes and the extra
  10 mass-bin rows/cell add negligible time (same detector calls; the bin
  computation itself is O(n)).
- **86**: 6 settings x 100 reps at the same per-cell cost -> ~2.9-3.1 h; the
  rewritten `local`/`bridge`/`collective` generators call the same
  `rpoisball.unit`/`mvrnorm` primitives at the same n, so their cost is
  unchanged from the original estimate.
- **87**: ~3.3-3.5 h with the RK-50/NND-100 rep split (down from 6.7 h at a
  uniform 100 reps), plus the new WP7 cluster-file block writes (negligible).

**88 does NOT fit this model and is far more expensive than
WP8_VERIFICATION.md's 3.1 h estimate.** U-MCCD/SU-MCCD (both RK-based) are
8-140x slower than the two-cluster benchmark on 88's SINGLE-POPULATION
generators at d=3 (see the cost finding above) -- their density-calibration
bracket search has no natural connectivity break to lock onto in a
homogeneous population. Per-setting-per-rep cost at d=3, summing all 13
methods (measured RK pair + ~5 s for the NND pair, 5 baselines and the
4-member MST sweep together, all comparatively negligible):
`uniform_control` ~281 s, `beta_gradient` ~72 s, `gauss_equal` ~92 s,
`gauss_unequal` ~23 s; at d=10 all four generators drop to ~12 s/rep. Over
100 reps and 4 generators: d=3 alone totals roughly **13 h**; d=10 adds
roughly 1.4 h. **Revised 88 full-grid estimate: ~14-15 h, one core** --
close to five times the original estimate, and larger than 85+86+87
combined (~7.9-7.5 h). Revised WP8 total (one core, all four scripts):
roughly **22-23 h**, not "roughly 11h". This is a launch-planning input,
not a correctness defect -- every measured run above completed and
returned a sensible (non-degenerate) score; it is simply slow. Two
mitigations available if the schedule needs it (neither applied here,
since only the 87 rep-count decision was authorized): parallelise 88 across
generators (each generator's cost is independent and the existing
per-cell `has_result()` restart already supports arbitrary interleaving),
or apply a 87-style reduced-rep decision to 88's RK methods specifically --
left to the user/coordinator, not decided unilaterally here.

## Dated appended notes (2026-09-05, round 2, post Opus re-verification WP8_REVERIFICATION.md)

Opus re-reviewed the round-1 fixes above (commits `764bc7f`, `605e103`,
`710d030`, `0e29d32`) and returned two required changes for `86` and one for
`87`; `85` and `88` passed outright and were launched unchanged. This note
records what changed and the acceptance measurements.

**86 `bridge` (blocking -- round 1's fix still failed).** Round 1's
clearance-0.45-with-isotropic-jitter construction still let the per-point
3-D `N(0,0.15)` jitter erode the axial margin: `t_lo > t_hi` (an empty or
inverted interpolation interval) in 37-41% of reps, collapsing the 10-point
chain to a compact blob at the midpoint -- geometrically near-identical to
`collective` -- on which every method scored TPR = 0.

Fix: clearance reduced to **0.05** (it only needs to clear the realised
surface, since jitter can no longer erode it) and jitter restricted to the
subspace **orthogonal to the cluster axis**: draw `g ~ N(0, I_d)`, subtract
its projection onto the unit axis vector `axis = (mu2-mu1)/CLS_DIS`
(`g_perp <- g - sum(g*axis)*axis`), scale by 0.03. Because the axial
component of the displacement is now exactly zero, Pythagoras guarantees
`distance_from_centre^2 = axial_distance^2 + perp_distance^2 >=
axial_distance^2`, i.e. perpendicular jitter can only INCREASE distance from
either `mu1` or `mu2` relative to the already-clear axial placement --
`n_inside = 0` is analytic, not merely likely, and the round-1 bounded-redraw
loop is removed entirely (nothing left for it to catch).

Acceptance measurement (100 reps, `d in {3, 10}`, production seed scheme
`BASE_SEED + 100000*s_i + rep`, `s_i` from the fixed CANON grid; script
sourced with `options(wp8.no_main=TRUE)` and `gen_bridge()` called directly,
recovering the realised jitter scales `s1, s2` by replaying the same seeded
`draw_cluster()` calls `gen_bridge()` itself makes -- no change to the
production RNG stream):

| d | n_inside (sum/100) | axis span / realised gap | chain NN <= host NN | end gap / chain spacing |
|---|---|---|---|---|
| 3  | 0/100 (max 0) | mean 0.801 (min 0.711, max 0.838) | 100/100 | mean 1.18 (min 0.95, max 1.75) |
| 10 | 0/100 (max 0) | mean 0.803 (min 0.738, max 0.841) | 100/100 | mean 1.18 (min 0.92, max 1.65) |

Definitions: "realised gap" = `CLS_DIS - s1 - s2` (axial extent between the
two clusters' realised surfaces); "axis span" = the chain's own first-to-last
axial extent (unaffected by the now-perpendicular-only jitter); "within-chain
NN" = per-point nearest-neighbour distance among the `n0` chain points,
summarised as the per-rep median; "host-cluster NN" = the pooled per-point
NN distance within `data1`/`data2`, per-rep median; "end gap" = the larger of
the two chain ENDPOINTS' own nearest-neighbour distance (checking the
uniform `t_i` interpolation does not leave an anomalously large gap at
either end), divided by the chain's median interior spacing. All four
criteria pass: `n_inside` is exactly 0 in every one of the 200 reps
(matching the analytic guarantee); axis span is 80%, comfortably above the
60% floor; the chain is at least as tightly packed as the host cluster in
every rep; end-gap ratio averages 1.18 at both d (individual reps range up
to 1.65-1.75, still close to the "~1.3x" target on average, not a hard
violation).

**86 `local` (framing -- two placements, not a fix).** Per
WP8_REVERIFICATION.md's decision, the arm now reports BOTH placements rather
than replacing one with the other:

- `local_shell` (round-1 construction, kept as-is, renamed): host cluster's
  regular points confined to a dense core (`core_frac = 0.4`, i.e. radius
  `0.4*scale`); each outlier at `runif(1, 0.8, 1.0)*scale` along a random
  direction from its host centre -- a shell placement, 0.49 beyond the host
  cluster's own realised radius (0.40) on average, not "embedded within" the
  cluster in the sense R3.7 describes. All 7 methods scored TPR = 1 on it in
  round-1 smoke.
- `local_interior` (restored verbatim from git `5ed690c`, the pre-rewrite
  construction WP8_VERIFICATION.md originally FAILED): a fixed 0.4*scale
  interior offset from the host centre, host regular points kept out of a
  0.3-radius ball around each outlier by rejection sampling. Its own
  `draw_cluster_avoid()` helper is restored alongside the shared `draw_cluster()`
  used by the other three generators (the two are not interchangeable: the
  interior construction needs per-point rejection sampling against a list of
  exclusion balls, which the shared core/shell `draw_cluster()` has no
  argument for).

Both clusters in `local_shell` are shrunk (`core_frac = 0.4`), making that
arm's clusters roughly **2.6x denser** (by NN spacing, which scales linearly
with radius at fixed n, so shrinking radius by 0.4 shrinks spacing by the
same factor -- reciprocal 1/0.4 = 2.5, matching the ~2.6x figure measured) than
`bridge`/`collective`'s full-radius clusters -- the four types are therefore
not density-matched, and this is reported as a known asymmetry rather than
corrected (correcting it would mean shrinking `bridge`/`collective`'s
clusters too, which is out of scope for this fixer pass).

Geometric fact carried over from WP8_VERIFICATION.md, restated here because
`local_interior`'s expected result depends on it: at d = 10, n = 95 (one
host cluster's regular-point count), the median within-cluster NN distance
is 0.699 against a realised radius of ~1.009 -- an interior point whose
isolation requires roughly DOUBLING the typical NN distance is not
realisable inside that radius at that d. **`local_interior` is therefore
expected to give near-zero TPR for every method at d = 10, and this is
itself the R3.7 finding (a legitimately embedded, locally-isolated outlier
is not geometrically constructible at high d with this cluster size) -- not
a construction defect to be re-fixed.** At d = 3 the same construction is
realisable and is expected to discriminate normally.

Stale comment fixed: `local_shell`'s header comment (86:100-116 before this
edit) had drifted to describe the FIRST attempt's core 0.5/shell 0.75-1.0
radii (measured ratio 1.89 at d=10, below the required 2.0, and therefore
never shipped); it now states the actual, adopted core 0.4/shell 0.8-1.0
radii (measured ratio 6.032 at d=3, 2.540 at d=10 -- see the round-1 note
above).

CANON, the settings selector, and the summariser all pick up the fourth type
automatically (`GENERATORS` is now `list(local_interior=, local_shell=,
bridge=, collective=)`; `build_settings()` iterates `names(GENERATORS)`, and
`do_summarize()` groups by `setting_id`/`type` generically) -- no
type-specific logic needed changing beyond the generator definitions
themselves. The type x d grid is now **8 settings**, not 6.

**87 (one-line CLI fix).** `87:275`'s reps-override parsing treated ANY
non-empty second argument as an override, including `"0"` -- so the intended
launch command `Rscript 87_wp8_small_cluster.R ALL 0 99999999` (meant to
supply `budget` while leaving `reps` at its per-family default) evaluated
`reps_ov <- as.integer("0") = 0L`, which is not `NULL`, so `do_run()`'s
`if (!is.null(reps_override)) reps_override else REPS_BY_FAMILY[[m_id]]`
picked **0 reps for every method** instead of the intended 50/100 split.
Fixed to `reps_ov <- if (length(args) >= 2 && nzchar(args[2]) &&
!is.na(suppressWarnings(as.integer(args[2]))) && as.integer(args[2]) > 0)
as.integer(args[2]) else NULL` -- only a positive integer counts as an
override. Dry-checked (no R run): with `args <- c("ALL","0","99999999")`,
`reps_ov` resolves to `NULL`, and the per-method loop then resolves to
U-MCCD 50, SU-MCCD 50, UN-MCCD 100, SUN-MCCD 100 -- the intended family
split, confirmed.

**Smoke re-verification (86, foreground, `results/tr1/wp8/smoke/` only).**
`--smoke` now covers all 8 canonical settings (both d, all 4 types), 1 rep,
9 default methods, with done-skip confirmed on re-invocation. See the commit
message / task report for the per-type, per-d TPR table.

**Bridge smoke caveat (2026-09-05, round 2, recorded not to overstate the
single-rep smoke result).** The fixed smoke replicate (`rep = 9001`, the
convention shared by all four WP8 scripts) gives TPR = 0 for all 9 methods
on `bridge` at BOTH d = 3 and d = 10 -- an unlucky draw, not a sign the
construction is degenerate. An ad hoc diagnostic (8 additional reps per d,
`rep in 1:8`, called directly against `gen_bridge()`/`METHOD_REGISTRY`
without touching any output file) found `n_inside = 0` in all 16 cells
checked and a NON-zero TPR for at least one method in 5/8 reps at d = 3
(MST 0.2-0.5, U-MCCD 0.3, LOF 0.2, ODIN 0.1) and 3/8 reps at d = 10 (U-MCCD
0.4, SU-MCCD 0.1-0.3, MST 0.2) -- so detection of the bridge chain is
genuinely possible but inconsistent across replicates, which is itself
consistent with R3.7's argument (a chain of mutually-close points can evade
methods that reward local density) rather than a construction defect. The
100-replicate production run is what will characterise the true detection
rate; the single smoke rep is a pipeline sanity check, not a result.

## Open questions / choices not fully specified by the revision plan

Recorded here rather than silently decided, per the declare-before-look
discipline:

1. Experiment 1/2's "usual contamination" level is set to 0.05, matching
   S_min numerically; this is a coincidence, not a claim that contamination
   should track S_min.
2. Experiment 2's three outlier counts are held at a fixed n0 = 10 (5%)
   rather than swept, since the plan asks for TPR *per type*, not a
   contamination sweep; a sweep is left to a future package if the response
   letter needs one.
3. Experiment 2's jitter scales (0.4 interior offset / 0.3 clearance for
   `local_interior`; 0.4/0.8-1.0 core/shell radii for `local_shell`; 0.05
   clearance / sd-0.03 perpendicular jitter for `bridge`; 0.15 jitter for
   `collective` -- see the round-2 dated note below for `bridge`'s history)
   and experiment 4's Beta(2,5)/variance-range choices are the author's; the
   plan names the distribution families but not their parameters.
4. Experiment 3's third cluster uses independent basis directions
   (`e1, e2`) rather than any specific angle; any two directions at
   `cls_dis` separation from both existing centres would serve equally.
5. The WP7 cluster-assignment files are written on a best-effort basis
   (see the atomicity notes above) rather than under the same restart
   guarantee as the metrics files, since the revision plan asks that they be
   *available*, not that they carry the same resumability contract.
