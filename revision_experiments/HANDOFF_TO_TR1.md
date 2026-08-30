# Handoff to the TR1 / Neurocomputing session

Written 2026-07-26 by the TR2 session (Pattern Recognition PR-D-26-05767, the
outlyingness-score paper). Everything below lives in the shared code repo
`Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection`, which both
projects check out.

Read section 1 first. It contains the only items that can change numbers TR1
has already published.

---

## 1. Bugs found in shared code that may touch TR1's published tables

These were found while running the TR2 revision experiments. All fixes are
committed. Each one sits in code that TR1's pipeline also uses, so TR1 should
check whether its own tables inherit the problem.

**1.1 Glass label-column bug — highest priority for an outlier-detection paper.**
`data/outlier_detection/RealData_Collection.R:34` sorts Glass by column 9
(the feature Fe) instead of column 10 (the label). Glass's 9 outliers therefore
land at rows 182-190, which breaks the outliers-last positional convention that
`count_scores2` assumes. Every metric computed for Glass through that path is
suspect. No other dataset is affected; all 14 exports were checked.

Measured consequence on the TR2 side: the published Glass LOF row
(TPR .778 / TNR .618 / BA .697 / F2 .289) keeps its TPR but its TNR does not
survive the corrected pipeline; the guarded reference is
.778 / .794 / .786 / .412. Part of that shift also comes from 1.2 below.
Separately, the published run used the deduplicated 9-outlier Glass while the
dataset-description table claims n0 = 10.

If TR1 reports Glass, regenerate the row and disclose the change in the response
letter. The harness fix is in `evaluate()` (commit 39126d8), which reorders
(Y, score) jointly rather than trusting position.

**1.2 LOF upper bound changed from 50 to 30** (commit 38013d9). Affects every
LOF baseline row, not only Glass.

**1.3 Vicinity-density overflow at high d** (commit history around 2026-07-16).
`methods/outlyingness_scores/Outlyingness_Score.R:12` computed `(size/R^d)^(1/d)`,
which evaluates `R^d` first and overflows to Inf once R > exp(709/d) — about 13.3
at d = 274, about 1212 at d = 100. Fixed to the algebraically identical
`size^(1/d)/R`. Regression test: `revision_experiments/tr2/test_outlyingness_density.R`.
Low-d agreement with the original formula is 2.1e-16, so nothing published at
d <= 100 moves. This only matters to TR1 if it runs the scores at high d.

**1.4 `Kest.R:449` gamma() overflow.** The constant term is finite through
d = 341, zero at 342, NaN at 343 and above. A log-space override is selected
automatically at d >= 342 in `revision_experiments/tr2/01_gen_quantile_table.R`;
the original file is untouched. Only relevant above d = 341.

**1.5 std_MADN fallback** (commit e138b85): MADN = 0 now falls back to SD, then
to 0, instead of zeroing every score in the cluster. Checked exhaustively in
`FINDINGS.md` section "High-d std_MADN recheck" — the fix moves AUC by at most
0.03 and never in the direction that flatters the scores.

---

## 2. Structural result TR1 may want to cite

**RK-CCD's spatial-randomness test loses its signal as d grows, and this is
measurable, not conjectural.** The share of zero entries in the Ripley-K CSR
null envelope rises from 0.1% at d = 5 to 7.5% at d = 10, 72% at d = 50, 91.4%
at d = 100, and 100% at d = 166, 274 and 1555. A zero envelope short-circuits
the radius search at the first candidate, so RK-CCD gets *faster* in higher
dimensions while doing less statistical work. Tables at d >= 166 can be
generated but carry no information.

Two consequences TR1 can use directly: any RK-CCD-based claim at high d needs
this caveat, and the degeneracy is a concrete, quantified answer to a reviewer
asking under what conditions the method fails.

Evidence: `FINDINGS.md` section 3; probe script in the TR2 scratchpad, and the
envelope percentages are reproducible from the tables in `R/RK-test_quantile/`.

---

## 3. Quantile-table inventory and its coverage limit

`R/RK-test_quantile/` (51 files) and `R/NN-test_quantile/` (81 files) hold the
production tables. The recovered tables were verified against the published WBC
row: 6 of 9 methods reproduce to the third decimal, max drift 0.011 (LOF).

Conventions, confirmed by that gate:
- RK quantile 99% for d <= 5, 999% for d >= 10.
- NN per-d anchors: 2d→85%, 3d→90%, 5d→95%, 10d→99%, 20d→999%.
- Score threshold 2; labels 1 = regular, 0 = outlier, outliers last.

**Coverage limit, easy to violate silently.** Low-d tables cover n up to 5000;
high-d tables (d >= 50) cover n up to 1000. Past coverage, `harness.R:153`
clamps to the last entry, i.e. silently uses the wrong critical values. TR1's
datasets are low-d (largest is Wilt at n = 4839), so they are inside the 5000
coverage — but check before adding anything larger.

New tables built for TR2, available for reuse: NN at d = 166, 274, 400 and 500.
Generation cost is real: d = 166 took about 15 h, d = 274 about 32 h, d = 400
about 43 h on 20 cores. `USER_RUN_TABLES.md` is the command sheet.

---

## 4. What can be reused as-is

**`revision_experiments/shared/harness.R`** — sourceable. Gives a quantile-table loader
`get_simul(variant, d, quant)`, a `METHOD_REGISTRY` of 9 scoring functions
(4 CCD-based outlyingness scores plus LOF, DBSCAN, MST, ODIN, iForest), an
`evaluate(Y, score, threshold)` returning TPR/TNR/BA/F2, and CSV checkpointing
so an interrupted run resumes. It only reads `R/`, `methods/`, `simulations/`
and `data/`; where a defensive fix was needed it redefines the function in the
global environment after sourcing, leaving the originals untouched.

Note for TR1: the registry does **not** include the MCCD detectors
(U-MCCD, SU-MCCD, UN-MCCD, SUN-MCCD) that are TR1's own methods. Those live in
`methods/outlier_detection/`. Adding them to the registry is the cheapest way to
get TR1 runtime numbers, AUCs, or new-baseline comparisons on the same footing.

**`revision_experiments/results/datasets_csv/`** — 16 datasets in one uniform
format with a `manifest.csv`: the ten original real datasets plus Arrhythmia,
Musk, Speech, InternetAds and n = 1000 subsamples of Musk and Speech.

**`revision_experiments/results/scores_cache/`** — 43 cached score vectors, so
AUCs and threshold sweeps can be recomputed without rerunning anything.

**PyOD environment** at `revision_experiments/.venv`: pyod 3.6.1, torch
2.13.0+cpu. ECOD, LUNAR and AutoEncoder are already run on the ten original
datasets and the high-dimensional ones — 154 of 154 fits, no failures. If TR1
also needs recent baselines, the results exist: `results/tr2/wp6_pyod_metrics.csv`,
with both contamination = 0.1 and true-rate labelings.

---

## 5. Experiments run for TR2, with result files

All under `revision_experiments/results/tr2/`. `FINDINGS.md` (272 lines, now
at `revision_experiments/tr2/FINDINGS.md`) is the running log and is the
authority where a number is disputed.

**WP3, cutoff sensitivity.** 18-point sweep from 0.5x to 2x the calibrated
cutoff. Three synthetic settings at 200 Monte Carlo replications each, plus WBC,
Thyroid and Arrhythmia. Balanced accuracy moves by at most 0.04 across the whole
grid in all twelve synthetic curves. Arrhythmia (d = 274) is the exception: no
plateau, BA falls from 0.66 to 0.53 as the cutoff grows.
Files: `wp3_sensitivity_synthetic.csv`, `wp3_sensitivity_real.csv`,
`wp3_synthetic_raw.csv`; figures in `results/tr2/figures/`.
Scripts: `09_wp3_synthetic.R`, `10_wp3_real.R`, `11_wp3_lines_plots.R`, `11b_wp3_synthetic_f2_plot.R` (all under `revision_experiments/tr2/`).

**WP4, runtime and scalability.** Wall-clock seconds for every method on two
grids: n in {100, 250, 500, 1000, 2000} at d = 10, and d in {5, 10, 50, 100, 500}
at n = 500. Ten repetitions per cell, three at n = 2000, on an idle machine
(24-core Intel Core Ultra 7 270K Plus, 64 GB DDR5, Windows 11, R 4.6.1).
Headline: UN-CCD construction dominates — 0.58 s at n = 100 rising to 9427 s at
n = 2000, a local slope near n^4, steeper than the O(n^3 log n) bound the papers
quote. Cost is flat in d, about 36 s from d = 10 to d = 500. Only the UN-CCD
radius search is parallelized (22 workers); everything else is single-threaded.
Files: `wp4_runtime2_{n,d,pyod}.csv` and the `_raw` per-repetition versions.
The older `wp4_runtime_*.csv` grids were measured under CPU contention and are
superseded — do not use them.
Scripts: `04_wp4_runtime.R`, `05_wp4_runtime_pyod.py`.

**WP5, high-dimensional real data.** Arrhythmia (452 x 274), Musk (3062 x 166),
Speech (3686 x 400), full data and n = 1000 subsamples. Reports rank AUC beside
the fixed-cutoff rates, because the calibrated cutoff does not transfer at these
dimensions. Two distinct failure modes are recorded: on Musk the IOS ranking
inverts (AUC 0.09 against 0.74 for OOS) because the outlier cluster is denser
than the genuine data; on Speech the scores go degenerate. InternetAds
(d = 1555) was dropped — no CCD variant can be built there.
Files: `wp5_highdim_metrics.csv`, `wp5_fulldata_raw.csv`, `wp5_subsample_raw.csv`.
Scripts: `06_wp5_highdim.R`, `07_wp5_subsample_ccd.R`, `07b_wp5_fulldata_ccd.R`.

**WP6, recent baselines.** ECOD, LUNAR, AutoEncoder at PyOD defaults, five seeds
for the two stochastic ones. Files: `wp6_pyod_metrics.csv`, `wp6_fit_log.csv`.
Scripts: `08_wp6_pyod_baselines.py`, `08b_wp6_metrics.R` (under `revision_experiments/tr2/`).

**High-d std_MADN recheck** (2026-07-24). Ten cells, 100 replications each,
testing whether "IOS beats OOS at high d" was an artifact of the standardization
bug. It was not — and the by-product is that OOS is the stronger score at high d
in every non-masking Gaussian and uniform cell. IOS's advantage is specific to
masking, not to dimension. Files: `highd_madn_recheck.csv` plus per-cell
per-replication dumps. Script: `12_highd_madn_recheck.R`.

---

## 6. Operating rules learned the hard way

- **Never edit a script while a detached Rscript is executing it.** R's
  incremental top-level parser reads by byte offset; editing mid-run produced a
  parse error hours into a job.
- **Never gate persistence of expensive compute on a verification assertion.**
  A 5.65 h Monte Carlo run was discarded by a diagnostic that was itself wrong.
  Save first, verify after.
- **Check `$LASTEXITCODE` in every chain step.** A PowerShell chain logged
  "done" for a step that had failed and marched on.
- **PowerShell `*>` writes UTF-16 logs.** bash `grep` fails silently on them.
- Invoke R as `C:\Program Files\R\R-4.6.1\bin\Rscript.exe` from PowerShell; the
  Bash tool segfaults on it.
- R was upgraded 4.4.1 to 4.6.1 mid-pipeline (2026-07-11, after a power outage).
  Results before and after that date were spot-checked and agree, but it is worth
  knowing if a number looks off by a hair.

---

## 7. What is not done

- The MCCD detectors have no runtime measurements, no AUCs, and no comparison
  against the three recent PyOD baselines. Everything needed to produce them is
  in place; only the registry entries are missing.
- InternetAds has no CCD-based rows at all.
- `wp4_runtime2_pyod.csv` aggregates all ten repetitions, while the TR2
  manuscript tables use repetitions 2-10 to discard a library-initialization
  warmup. Both statistics are reproducible from `wp4_runtime2_pyod_raw.csv`.
