# Running findings log — revision experiments (orchestrator-maintained)

Feeds the final EXPERIMENT_RESULTS.md (T8). Each entry: what was found, evidence, consequence.

## 1. Recovered quantile tables are the production tables (T2, 2026-07-10)
WBC reproduction gate: 6 of 9 methods reproduce the manuscript's real-data WBC row to the 3rd decimal; max drift 0.011 (LOF). Conventions confirmed: RK quantile 99% for d<=5 / 999% for d>=10; NN per-d anchors (2d→85%, 3d→90%, 5d→95%, 10d→99%, 20d→999%); OS threshold = 2; label polarity 1=regular/0=outlier, outliers last.

## 2. Glass loader bug + suspect published Glass LOF row (T2 follow-up, 2026-07-10)
`data/outlier_detection/RealData_Collection.R:34` sorts glass by column 9 (feature Fe), not column 10 (label) → glass's 9 outliers sit at rows 182-190, violating the outliers-last convention `count_scores2` assumes positionally. Harness `evaluate()` now jointly reorders (Y, score) regulars-first (commit 39126d8); no other dataset affected (all 14 exports checked).
Consequence for the manuscript: published Glass LOF row (TPR .778/TNR .618/BA .697/F2 .289) has matching TPR but TNR inconsistent with the corrected pipeline (guarded reference: .778/.794/.786/.412; also affected by the LOF UB=50→30 fix in commit 38013d9). The Glass row (at least LOF; possibly other methods) should be regenerated for the revised manuscript, and the change disclosed in the response letter if numbers move.
Also: published run used the deduped 9-outlier glass while the paper's dataset table claims n0=10 — check the dataset-description table when revising.

## 3. RK-CCD quantile generation produces NaN at d=1555 (T3a, in progress)
Probe failed with NaN — original RK math breaks at very high d (suspected overflow/underflow in d-dependent terms). Diagnosis in progress; fallback if unfixable: UN-CCD-based scores only on InternetAds, disclosed in response letter.

## 4. WP6 PyOD baselines complete (T6, 2026-07-10, commit 42caaed)
154/154 fits, no failures. ECOD strongest of the three (Musk raw-score ROC-AUC 0.956, matches ADBench ~0.95). Speech near-chance for all (BA≈0.51) — consistent with literature. pyod 3.6.1, torch 2.13.0+cpu, defaults, LUNAR/AE 5 seeds. Metrics: results/wp6_pyod_metrics.csv (both contamination=0.1 and true-rate labelings).

## 5. Environment quirks (T0/T2/T6)
- Bash tool segfaults invoking Rscript → use PowerShell.
- pip into the venv needs `subst X:` long-path workaround.
- `Rscript --%` (PowerShell stop-parse) breaks argument handling — invoke without it.
- Windows 260-char path limit: repo needs `core.longpaths true` (set locally) for 3 deep SVDD files.
- UNCCD registry entries pay construction twice when t_construct is timed separately — skip separate timing for long cells.
