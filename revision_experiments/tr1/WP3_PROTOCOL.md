# WP3 protocol — parameter-selection protocol statement, variability, paired tests

**Declared 2026-09-05, before `tr1/92_wp3_stats.R` is written.**

WP3 answers R3.9, R5.3 and R5.4: the manuscript compares seventeen methods
(the five original baselines, the four proposed MCG detectors, and the eight
WP4 competitors) without a single stated rule for how each method's
parameters were chosen, and reports no variability or significance test
behind any of the real-data numbers. This file fixes, in advance of running
anything, (a) the design of the paired significance tests on the sixteen
real data sets, (b) what the manuscript's Monte Carlo simulation tables would
need for a standard-deviation / Monte Carlo standard-error column and
whether the data to compute it exist, and (c) the shape of the parameter
provenance table that the manuscript's protocol subsection (Section 5.1)
will be built from. Nothing here is re-derived from scratch: the per-cell
real-data table WP3 tests is the exact one WP4 already assembled and the
manuscript already prints (`tab:Real_Data_Aggregate`, `tab:Real_Data_Result_Summary`),
via `tr1/82_wp4_metrics.R`'s `fx` (nine incumbent methods, alpha- and
S_min-corrected) and `combined("T1")` (fx plus the eight WP4 competitors at
the contamination-0.1 threshold, T1 = the main-table setting). WP3 sources
that script and reuses its objects; it does not re-read raw scores or
re-derive labels.

## (a) Paired-test design on the sixteen real data sets

**Population.** The sixteen data sets of `tab:Real_Data`, one paired
observation per data set, for each of the two metrics reported in
`tab:Real_Data_Aggregate`: $F_2$ and BA. The values compared are the exact
rounded (3 dp) numbers the manuscript prints — `fx` for the nine incumbents
(after the vertebral/SU-MCCD S_min patch and the alpha-repair overrides) and
`combined("T1")` for the eight WP4 competitors, i.e., DIF and LUNAR at their
5-seed mean, mutual-kNN and SNN at their pre-declared oracle-$k$, and the
remaining four (ECOD, COPOD, GLOSH, OPTICS) at their single deterministic
fit. Using the published rounded values, not intermediate unrounded ones,
means the test operates on the numbers a reader can check against the
tables.

**Primary comparison.** SUN-MCCD vs. each of the other sixteen methods (the
eight incumbents — U-MCCD, SU-MCCD, UN-MCCD, LOF, DBSCAN, MST, ODIN,
iForest — and the eight WP4 competitors), separately for $F_2$ and for BA.
That is 16 comparisons × 2 metrics = 32 tests.

**Secondary comparison.** The same design run from the perspective of each
of U-MCCD, SU-MCCD, and UN-MCCD in turn: each vs. the other sixteen methods
(which now includes SUN-MCCD), on both metrics. 3 methods × 16 comparisons ×
2 metrics = 96 tests. These are secondary because SUN-MCCD is the
manuscript's headline recommendation (see `CLAUDE.md`); the other three are
didactic reference constructions, and the point of running the same design
for them is to show whether the shape-adaptive / NND-based gains they
individually contribute are themselves distinguishable from noise.

**Test.** The exact two-sided paired Wilcoxon signed-rank test,
`wilcox.test(x, y, paired = TRUE, exact = NULL)` from base R (`stats`).
`coin::wilcoxsign_test` was considered — it offers an explicit
`zero.method`/`ties.method` argument — but the `coin` package is **not
installed** in this R installation (`C:/Program Files/R/R-4.6.1`, checked
with `requireNamespace("coin", quietly = TRUE)` on 2026-09-05, returns
`FALSE`). WP3 therefore uses base R's `wilcox.test`, which is already exact
for this sample size whenever the data allow it.

**Tie and zero handling — R's own default, not overridden.**
`wilcox.test(..., exact = NULL)` (the default) uses the exact null
distribution automatically when the number of *non-zero* paired differences
is below 50 and none of the non-zero differences share the same absolute
value; both hold trivially for the non-zero part of a 16-observation sample,
so the test is exact unless the data themselves produce a tied absolute
difference. Concretely:

- A **zero difference** (the primary and opponent method print the identical
  rounded value on a data set) is dropped before ranking — this is
  `wilcox.test`'s own behavior, not a WP3-specific choice — and reduces the
  effective sample size for that comparison. The dropped count is reported
  as `n_tied` in the output and the comparison is still reported (not
  skipped) whenever at least one non-zero difference remains.
- A **tied absolute value among the non-zero differences** makes the exact
  null distribution inapplicable; `wilcox.test` then falls back to the
  normal approximation with a continuity correction and would ordinarily
  emit a warning under `exact = TRUE`. WP3 does not force `exact = TRUE`; it
  leaves `exact = NULL` so R silences that warning and makes the fallback
  itself, and records **whether the exact or the approximate branch fired**
  in an `exact_used` column, so this is auditable per comparison rather than
  asserted globally.
- We do **not** force `exact = FALSE` (all-approximate) or hand-implement a
  mid-p correction; whichever branch base R's own default takes is reported
  as-is, because the point of the test here is a standard, reproducible
  significance check, not a bespoke correction for the small sample.

**Multiplicity.** Holm's step-down adjustment (`p.adjust(method = "holm")`)
is applied separately within each block of 16 raw p-values — one block per
(primary method, metric) pair. So SUN-MCCD's 16 $F_2$ comparisons are
Holm-adjusted among themselves, its 16 BA comparisons among themselves, and
likewise for each of the three secondary primary methods. Blocks are not
pooled across metrics or across primary methods, because each block answers
a self-contained question ("does SUN-MCCD beat the field on $F_2$") and
pooling would make the correction depend on how many secondary analyses were
run, which is an analysis-design choice rather than a property of the
comparison itself.

**Ahead/behind/tied count.** For every comparison, alongside the test
statistic, WP3 reports the number of data sets (of 16) on which the primary
method's rounded value is strictly greater (`n_ahead`), strictly less
(`n_behind`), or equal (`n_tied`) to the opponent's, for that metric. This is
the sign-test-level summary that motivates and cross-checks the Wilcoxon
result; a large `n_ahead`/`n_behind` imbalance with a non-significant
Wilcoxon test (e.g., because the wins are small in magnitude) is exactly the
distinction the manuscript's protocol subsection needs to state honestly.

**What this does not claim.** Sixteen data sets is a small, heterogeneous,
non-random sample (ODDS/ELKI benchmarks of convenience), so a significant
paired Wilcoxon test here supports "the ordering observed on these sixteen
benchmarks is not readily explained by chance ordering of wins and losses of
this magnitude" — it is not a claim about performance on some superpopulation
of tasks. The manuscript's protocol subsection states this scope limitation
alongside the numbers.

## (b) Simulation-variability plan

**Rule.** For every simulation cell the manuscript reports as a mean (a
number in prose, a figure line/point, or a table entry) over its 1000 (main
synthetic grid) or fewer (robustness/Neyman–Scott) replicates, WP3 computes
the standard deviation over replicates and the 95% Monte Carlo standard
error, $\mathrm{SE} = \mathrm{SD}/\sqrt{R}$ with a normal-approximation 95%
half-width $1.96\,\mathrm{SE}$ — **conditional on the per-replicate raw
output existing**. This sub-package does not re-run any simulation; WP1–WP2
already established the write-per-cell discipline this repository requires
(`CLAUDE.md`: "every long run must append its results per cell"), and WP3
only reads what is already on disk.

**What exists, checked 2026-09-05.** Every driver behind the manuscript's
Section 5 synthetic experiments (`sec:synthetic-experiments`: uniform and
Gaussian clusters over $d\in\{2,3,5,10,20,50,100\}$,
$n\in\{50,100,200,500,1000\}$, 1000 replicates; `sec:MC-focus-settings`: the
six robustness checks; `sec:rand-clust-proc`: the Matérn/Thomas/mixed
Neyman–Scott experiments) lives under
`Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection/simulations/outlier_detection/{RU-MCCDs,SU-MCCDs,UN-MCCDs,SUN-MCCDs,Algo_Compare_OutlierDetection}/Simulation*/`.
Every one of those directories holds only two file types: the cluster-job
driver script (`*.R`) and its SLURM stdout log (`slurm-*.out`). A full
recursive search (`find simulations/outlier_detection -type f`) over 2,223 R
files and 1,973 `.out` logs returns **zero** `.csv`, `.rds`, or `.RData`
files anywhere in the tree. Each driver script ends with
`save.image(...)` to an `.RData` path in the same directory — that is where
the raw per-replicate TPR/TNR/BA/$F_2$ vectors would live — but `*.RData` is
listed in this repo's `.gitignore` and the files themselves are not present
in this checkout; they were evidently written on the compute cluster and
never copied back, or were cleaned up afterward. The SLURM logs print only a
single aggregate line per job (e.g. `"The mean success rate is 1 , and, the
mean True positive rate is 0.9595..."`), not a per-replicate vector, so no
SD is recoverable from the logs either.

**Conclusion for the manuscript's own headline simulation numbers: no
per-replicate data exist in this repository, for any of the three simulation
families (main synthetic grid, six robustness checks, Neyman–Scott).**
`92_wp3_stats.R` therefore cannot compute a single SD or Monte Carlo SE for
any number currently printed as a mean in Sections 5.2–5.4 or their
supporting figures. This is recorded as `results/tr1/wp3/sim_variability_INVENTORY.md`,
not asserted from memory — the inventory file lists every directory searched
and the exact zero-file result.

**What does exist, and why it does not substitute.** `55_wp2c_simulation_arm.R`
(a WP2(c) script, not owned by WP3) wrote genuine per-replicate rows to
`results/tr1/wp2c_simulation_arm.csv` (10,609 rows) while comparing the
S_min rule candidates. Its cells are real per-replicate data, but they answer
a different, narrower question — S_min rule sensitivity — with a design that
does not match the manuscript's reported cells: only $d=10$; only
$n\in\{200,500\}$; only SU-MCCD and SUN-MCCD (not the other proposed methods,
and no baselines); 100 replicates per cell, not 1000; and contamination and
generator combinations chosen for the S_min study, not for reproducing a
specific figure point. A direct check (gaussian, $d=10$, contamination 5%,
$n=200$, SUN-MCCD) gives mean $F_2=0.787$ against the $F_2=0.830$ the main
text reports at $d=10$ for the Gaussian setting — the two numbers do not
match, confirming this file is not the source of any manuscript-printed
cell. `92_wp3_stats.R` still computes the per-cell SD and MC-SE from this
file, in `results/tr1/wp3/sim_variability.csv`, labeled throughout as **not
corresponding to any manuscript-reported value** — it is the only real
evidence in this repository of the magnitude of Monte Carlo noise the
proposed methods' simulation cells carry (SD of $F_2$ around 0.01–0.07 per
cell at $n=100$–$200$ replicates in that file), offered as an indicative
lower bound on what the true (unavailable) per-cell SDs would look like at
the manuscript's own, larger, 1000-replicate cells.

**Response-letter framing.** R5.4's request for "standard deviations,
confidence intervals, or statistical significance tests" is answered for the
real-data comparison in full (part (a)); for the synthetic simulation grid,
it is answered by disclosure that the underlying per-replicate output no
longer exists in this repository and by the indicative WP2(c) noise-level
evidence, not by a re-run — a full-grid re-run (7 dimensions × 5 sample
sizes × 1000 replicates × 2 generators × 4 proposed methods, plus the six
robustness checks and three Neyman–Scott settings) is outside WP3's
same-day budget and is not requested by name in R3.9/R5.3/R5.4 the way the
real-data paired test is. This is flagged as a follow-up decision for the
user: re-run the grid to get simulation SDs (a multi-hour to multi-day
compute job), or state in the response letter that per-replicate variance
was not retained from the original compute run and report the indicative
WP2(c) figures instead.

## (c) Parameter provenance table (draft)

One row per (method, parameter): all seventeen methods used in the real-data
study, with DBSCAN's `Eps` broken out as its own row because it is the one
parameter R3.9 and R5.3 name directly ("If ... the parameters of baseline
methods are selected according to the actual anomaly proportion"). The
`how set` column uses a fixed, closed vocabulary: `library default`,
`fixed constant`, `fixed rule (function of n or d)`, `label-free rule`,
`uses true contamination`, `oracle-k` (or `oracle` generally, where the true
labels are consulted). The full table is generated by `92_wp3_stats.R` from
the values already fixed and cited in `CCD_OutlierDetection_Neurocomputing.tex`
Section 6 (real-data baseline settings, referring back to Section 5.3 for
LOF/ODIN/iForest) and `tr1/WP4_PROTOCOL.md` §2, and written to
`results/tr1/wp3/parameter_provenance.csv`; a rendering is embedded below for
reference at declaration time.

| Method | Parameter | Value | How set | Source |
|---|---|---|---|---|
| LOF | $k$ | $\{11,\ldots,30\}$, largest LOF retained per point | fixed constant (literature-recommended envelope) | main text, `sec:rand-clust-proc`, reused unchanged for real data |
| LOF | outlier threshold | $1.5$ | fixed constant (matches Breunig et al.'s reported optimum) | main text, `sec:rand-clust-proc` |
| DBSCAN | MinPts | $4$ | fixed constant | main text, `sec:rand-clust-proc` and `sec:Real-Data-Examples` |
| **DBSCAN** | **Eps** | **4-distance heuristic evaluated at each data set's own true contamination rate** | **uses true contamination** | main text, `sec:Real-Data-Examples`: "sets Eps from the 4-distance heuristic at the known contamination rate of each data set" |
| MST | inconsistent-edge threshold | $1.2$ (real data); $1.7\to1.1$ decreasing with $d$ (Neyman–Scott simulation) | fixed constant, held across all 16 real data sets (differs by design from the simulation study, which varies it with $d$ but not with the data) | main text, `sec:Real-Data-Examples` and `sec:rand-clust-proc` |
| MST | minimum cluster size to avoid an outlier flag | $2\%$ of $n$ | fixed constant | main text, `sec:Real-Data-Examples` |
| ODIN | $k$ | $\lceil\sqrt{n}\rceil$ | fixed rule (function of $n$, label-free) | main text, `sec:rand-clust-proc` |
| ODIN | in-degree threshold $T$ | $\lceil n^{1/3}\rceil$ | fixed rule (function of $n$, label-free) | main text, `sec:rand-clust-proc` |
| iForest | number of trees | $1000$ | fixed constant (library-recommended) | main text, `sec:rand-clust-proc` |
| iForest | sub-sample size | $256$ | fixed constant (library-recommended) | main text, `sec:rand-clust-proc` |
| iForest | outlier-score threshold | $0.55$ | fixed constant | main text, `sec:rand-clust-proc` |
| U-MCCD, SU-MCCD (RK-based) | MC-SRT significance level $\alpha$ | step function of $d$: $1\%$ ($d<10$), $0.1\%$ ($d\geq10$) | fixed rule (function of $d$, pre-declared, label-free) | main text, `sec:Uniform_clusters`; realized values per real data set in `tab:alpha_real` |
| U-MCCD, SU-MCCD | KS-CCD density-search tolerance $\tau$ | $0.01$ | fixed constant (sensitivity checked over $10^{-3}$–$1$; WP2(a)) | Supplementary Material Algorithm boxes; sensitivity result in supplement around line 507 |
| U-MCCD, SU-MCCD | Monte Carlo replications behind the RK quantile tables | $5000$ (99% level tables) / $10000$ (99.9% level tables) | fixed constant, computed once and stored (no randomness at detection time) | `R/RK-test_quantile/*d_99%.R`, `*d_999%.R` |
| SU-MCCD | $S_{\min}$ (minimum cluster fraction) | $0.05$ | fixed constant, label-free (same value for every real data set and for the simulations) | main text, `sec:Real-Data-Examples`: "we set $S_{\min}=0.05$ ... so that no parameter depends on the contamination level" |
| UN-MCCD, SUN-MCCD (NND-based) | MC-SRT significance level $\alpha$ | step function of $d$, own schedule (UN-MCCD: $15\%,10\%,5\%,1\%,0.1\%$ at $d=2,3,5,10,\{20,50,100\}$; SUN-MCCD: $0.1\%$ already at $d=10$) | fixed rule (function of $d$, pre-declared, label-free) | main text, `sec:Uniform_clusters`; realized values per real data set in `tab:alpha_real` |
| UN-MCCD, SUN-MCCD | Monte Carlo replications behind the NND quantile tables | $5000$ / $10000$, by level, as above | fixed constant | `R/NN-test_quantile/*.R` |
| SUN-MCCD | $S_{\min}$ | $0.05$ | fixed constant, label-free | main text, `sec:Real-Data-Examples` |
| ECOD | none | — (parameter-free by construction) | library default | `tr1/WP4_PROTOCOL.md` §2 |
| COPOD | none | — (parameter-free by construction) | library default | `tr1/WP4_PROTOCOL.md` §2 |
| DIF | architecture/training defaults (`hidden_neurons=[500,100]`, `n_ensemble=50`, `n_estimators=6`, `max_samples=256`, ...) | pyod 3.6.1 `get_params()` defaults | library default | `tr1/WP4_PROTOCOL.md` §5 |
| DIF | random seed | $\{1,\ldots,5\}$, mean over seeds reported | fixed constant (5 seeds, pre-declared) | `tr1/WP4_PROTOCOL.md` §2 |
| LUNAR | architecture/training defaults (`n_neighbours=5`, `n_epochs=200`, `model_type='WEIGHT'`, ...) | pyod 3.6.1 `get_params()` defaults | library default | `tr1/WP4_PROTOCOL.md` §5 |
| LUNAR | random seed | $\{1,\ldots,5\}$, mean over seeds reported | fixed constant (5 seeds, pre-declared) | `tr1/WP4_PROTOCOL.md` §2 |
| GLOSH (HDBSCAN) | `min_cluster_size` | $5$ | library default | `tr1/WP4_PROTOCOL.md` §2, §5 |
| GLOSH (HDBSCAN) | `min_samples` | `None` (falls back to `min_cluster_size`) | library default | `tr1/WP4_PROTOCOL.md` §5 |
| OPTICS | `min_samples` | $5$ | library default | `tr1/WP4_PROTOCOL.md` §5 |
| OPTICS | `xi` | $0.05$ | library default | `tr1/WP4_PROTOCOL.md` §5 |
| mutual-kNN | $k$ | best of $\{5,10,15,20,30\}$ by $F_2$ under the true labels, per data set | **oracle-$k$ (uses true labels)** | `tr1/WP4_PROTOCOL.md` §2, "the oracle concession on $k$" |
| SNN | $k$ | best of $\{5,10,15,20,30\}$ by $F_2$ under the true labels, per data set | **oracle-$k$ (uses true labels)** | `tr1/WP4_PROTOCOL.md` §2 |
| all 8 WP4 competitors (score-based label rule, T1) | contamination fraction used to threshold the score | $0.1$, identical for every data set | library default, label-free | `tr1/WP4_PROTOCOL.md` §3, main text `sec:Real-Data-Examples` |

The table intentionally separates a method's *fitting* parameters (rows
above the last one) from the *thresholding* contamination used only to turn
a continuous score into a label (the last row) — the two are conflated in
some of the reviewer language and the manuscript's protocol subsection needs
to keep them apart, since a `library default` fitting parameter can still sit
behind a labelling step that is `label-free` (all 8 WP4 competitors at T1)
or that is not (DBSCAN's Eps; mutual-kNN and SNN's oracle-$k$).

## Deliverables and outputs

```
results/tr1/wp3/wilcoxon_real.csv              32 primary + 96 secondary paired-test rows
results/tr1/wp3/sim_variability.csv            per-cell SD/MC-SE from wp2c_simulation_arm.csv (not manuscript cells; see (b))
results/tr1/wp3/sim_variability_INVENTORY.md   what per-replicate simulation output exists and does not, and why
results/tr1/wp3/parameter_provenance.csv       the table in (c), all rows
```

Gate: before any test is run, `92_wp3_stats.R` re-asserts (does not merely
trust) that the nine-method aggregate it sources from `82_wp4_metrics.R`
reproduces the same 36 values printed in `tab:Real_Data_Aggregate`, by
literally re-running the PRINTED-vs-computed comparison 82 performs, on the
`fx` object 82 already built. If that gate fails, nothing downstream is
computed.
