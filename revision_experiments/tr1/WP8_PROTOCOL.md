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
   "local", 0.15 for "bridge" and "collective") and experiment 4's
   Beta(2,5)/variance-range choices are the author's; the plan names the
   distribution families but not their parameters.
4. Experiment 3's third cluster uses independent basis directions
   (`e1, e2`) rather than any specific angle; any two directions at
   `cls_dis` separation from both existing centres would serve equally.
5. The WP7 cluster-assignment files are written on a best-effort basis
   (see the atomicity notes above) rather than under the same restart
   guarantee as the metrics files, since the revision plan asks that they be
   *available*, not that they carry the same resumability contract.
