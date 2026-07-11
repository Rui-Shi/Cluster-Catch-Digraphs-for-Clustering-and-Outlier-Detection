# Running findings log — revision experiments (orchestrator-maintained)

Feeds the final EXPERIMENT_RESULTS.md (T8). Each entry: what was found, evidence, consequence.

## 1. Recovered quantile tables are the production tables (T2, 2026-07-10)
WBC reproduction gate: 6 of 9 methods reproduce the manuscript's real-data WBC row to the 3rd decimal; max drift 0.011 (LOF). Conventions confirmed: RK quantile 99% for d<=5 / 999% for d>=10; NN per-d anchors (2d→85%, 3d→90%, 5d→95%, 10d→99%, 20d→999%); OS threshold = 2; label polarity 1=regular/0=outlier, outliers last.

## 2. Glass loader bug + suspect published Glass LOF row (T2 follow-up, 2026-07-10)
`data/outlier_detection/RealData_Collection.R:34` sorts glass by column 9 (feature Fe), not column 10 (label) → glass's 9 outliers sit at rows 182-190, violating the outliers-last convention `count_scores2` assumes positionally. Harness `evaluate()` now jointly reorders (Y, score) regulars-first (commit 39126d8); no other dataset affected (all 14 exports checked).
Consequence for the manuscript: published Glass LOF row (TPR .778/TNR .618/BA .697/F2 .289) has matching TPR but TNR inconsistent with the corrected pipeline (guarded reference: .778/.794/.786/.412; also affected by the LOF UB=50→30 fix in commit 38013d9). The Glass row (at least LOF; possibly other methods) should be regenerated for the revised manuscript, and the change disclosed in the response letter if numbers move.
Also: published run used the deduped 9-outlier glass while the paper's dataset table claims n0=10 — check the dataset-description table when revising.

## 3. RK-CCD quantile machinery breaks down at high d (T3a, complete 2026-07-11, commit 0160701)
Numerical: `Kest.R:449` cons term uses gamma() → finite through d=341, 0 at d=342, NaN at d>=343. Fixed via log-space stable override in 01_gen_quantile_table.R (auto-selected at d>=342; validated to 2.8e-15 relative diff vs original at d=20, identical at d=100; originals untouched).
Statistical (the deeper problem): the RK CSR envelope is already 91.4% zero quantiles in the production d=100 table and 100% zero at d=166/274/1555 probes — the Monte-Carlo spatial-randomness test carries no signal at these dimensions. RK tables at d>=166 are generatable but scientifically vacuous.
Consequence: WP5 CCD rows at d>=166 should be UN-CCD-based only (UNCCD-OOS/IOS); this is also manuscript-relevant — it is a concrete "condition under which the approach degrades" answering R2's failure-conditions comment, worth a sentence in the paper.
Production parameters reconciled: two tiers in recovered tables — low/mid d: RK m=5000, rn=10, niter=2000; NN n=5000. High d (50, 100): m=n=1000, rn=10, niter=10000 (matches committed 50100d driver). New dims adopt the high-d tier.
NN table costs (measured probes, cores=20): d=166 ≈ 15.4 h, 274 ≈ 31.5 h; d>=400 needs two validated generator optimizations (mvrnorm(0, diag(d)) → matrix(rnorm) — distributionally identical, kills an O(d^3) eigen per call; streaming to cut RAM) → 400 ≈ 43 h, 1555 ≈ 257 h at niter=10000 or ≈ 51 h at niter=2000 (matches the low-d production tier; needs sign-off).

## 3b. Unattributed files in R/NN-test_quantile/ (2026-07-11)
`35-99d_999%.R` (new, loops NN table generation d=35:99 at niter=10000), whitespace-only edit to `20d_999%.R`, and a fresh `.Rhistory` — user's own exploratory R session (confirmed 2026-07-11). Left uncommitted; not part of the pipeline; not running.

## 3c. User scope decisions (2026-07-11)
- WP5 CCD rows at d>=166: **UN-CCD only** (UNCCD-OOS/IOS). RKCCD rows dropped there on degeneracy evidence; the paper gains a failure-conditions sentence (answers R2).
- **InternetAds (d=1555) CCD experiments skipped for now** — baselines + PyOD cover it; NN-1555 cost analysis (11 d full / 2 d thinned) retained as documentation for the response letter.
- WP4 runtime d-grid: **{10, 50, 100, 166, 274, 400}** (500/1000 dropped; 166-400 reuse WP5 tables; RKCCD cells auto-skip where no table).
- Phase B therefore: NN tables at d ∈ {166, 274, 400} only (~15 + 32 + 43 h), after implementing + validating T3a's two generator optimizations (rnorm swap — distributionally identical; streaming for RAM).

## 4. WP6 PyOD baselines complete (T6, 2026-07-10, commit 42caaed)
154/154 fits, no failures. ECOD strongest of the three (Musk raw-score ROC-AUC 0.956, matches ADBench ~0.95). Speech near-chance for all (BA≈0.51) — consistent with literature. pyod 3.6.1, torch 2.13.0+cpu, defaults, LUNAR/AE 5 seeds. Metrics: results/wp6_pyod_metrics.csv (both contamination=0.1 and true-rate labelings).

## 4b. WP5 baseline rows complete (T5, 2026-07-11, commit 2b56120)
20 rows (5 baselines × 4 new datasets) in results/wp5_highdim_metrics.csv; UNCCD rows pending user-generated NN tables (6 skip markers). Paper-relevant observations, all verified non-bugs:
- iForest near-ceiling on Musk (BA 0.943, rank-AUC 1.000 — matches ADBench) but TPR≈0 on Arrhythmia/InternetAds/Speech purely because the fixed 0.55 threshold sits above its whole score range there (ranking still good: AUC 0.804 on Arrhythmia). Threshold-calibration artifact worth a sentence.
- MST degenerates on Musk (TPR=1, TNR=0): pruned-MST components all fall below the min-cluster size, so everything is flagged — distance concentration at d=166 flattens edge ratios. Genuine high-d failure mode of a first-cycle baseline.
- Speech (3686×400): every method near chance (best BA 0.527; AUCs 0.475-0.529), consistent with PyOD rows and the literature.
- Coherence with WP6 PyOD numbers confirmed dataset-by-dataset; no contradictions.

## 4c. WP3 sensitivity: real sweeps + probe complete (T7a, 2026-07-11)
wp3_sensitivity_real.csv: WBC + Thyroid × 4 OS methods × 18 cutoffs (144 rows); WBC reproduction gate at cutoff=2 PASSED (all 4 methods match published values). All 5 figures rendered (2 real + 3 synthetic-probe). Synthetic probe: 30 reps in 34.7 min on 20 cores (gaussian 107 s/rep, uniform 89, matern 11); production sweep at 200 reps/setting ≈ 11.5 h — queued after the WP4 timing grids. Thyroid timings: RKCCD-OOS 587 s construct / RKCCD-IOS 514 s; UNCCD completed overnight (caches present).

## 4d. NN fast generator validated; user command sheet ready (T3b, 2026-07-11)
01e_nn_fast.R: rnorm swap (rotation-invariance argument, kills O(d^3) eigen per draw) + streaming SimuOnce. Validation (nn_fast_validation*.log): statistical agreement original-vs-fast within the different-seed MC yardstick (d=50: cor 0.995, mean rel diff 0.4%); per-worker RAM at d=400 2.29 GB → 0.27 GB; measured speedups 1.08x (d=166) to 1.29x (d=400) plus the cores=20 unlock. A worker-export bug (rpoisball.unit.fast not exported to PSOCK workers) was caught in validation step C and fixed before it could reach a production run. USER_RUN_TABLES.md: three commands (NN 166/274/400, niter=10000, cores=20), est. ~14 h / ~28 h / ~37-42 h, no checkpointing (a killed run restarts from zero), full-R-path invocations.

## 4e. R upgraded 4.4.1 → 4.6.1 mid-pipeline (discovered 2026-07-11, after power outage)
R 4.4.1 was uninstalled (stale PATH entry remains); only 4.6.1 exists. 9 of 17 required packages were missing from 4.6.1's library and scales/ggplot2 were corrupted ("unknown input format", plausibly the power outage) — all reinstalled; 00_env_check.R passes under 4.6.1; the exact user table command completed a toy run. Consequences: all metrics so far (computed under 4.4.1) remain valid; ALL WP4 timing runs execute uniformly under 4.6.1 (none had run yet, so no mixed-version timings); report R 4.6.1 in the manuscript's runtime section provenance.

## 5. Environment quirks (T0/T2/T6)
- Bash tool segfaults invoking Rscript → use PowerShell.
- pip into the venv needs `subst X:` long-path workaround.
- `Rscript --%` (PowerShell stop-parse) breaks argument handling — invoke without it.
- Windows 260-char path limit: repo needs `core.longpaths true` (set locally) for 3 deep SVDD files.
- UNCCD registry entries pay construction twice when t_construct is timed separately — skip separate timing for long cells.
