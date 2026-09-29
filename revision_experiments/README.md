# revision_experiments

This directory is shared by two paper revisions that both build on the same
CCD/MCCD codebase (`data/`, `methods/`, `R/`, `simulations/` at the repo
root):

- **tr1/** -- Neurocomputing, "Shape-Adaptive Outlier Detection Using Cluster
  and Mutual Catch Digraphs" (NEUCOM-D-26-15191).
- **tr2/** -- Pattern Recognition, "Outlyingness Scores with Cluster Catch
  Digraphs..." (PR-D-26-05767).
- **shared/** -- infrastructure both projects source.

Reorganized 2026-08-30 into a full split: each project's scripts and docs now
live in their own subdirectory (`tr1/`, `tr2/`), not just their `results/`.
This supersedes the 2026-08-09 partial reorg, which split only `results/`
into `results/tr1/` / `results/tr2/` and left every script at the top level
of `revision_experiments/` -- see git history for that rationale if needed;
it no longer describes the current layout.

2026-09-27 tr1 consolidation: every tr1 script and document now sits in
`tr1/`. Moved from `revision_experiments/` with names unchanged: the WP2
scripts `40_*` .. `78_*`, `01i_nn_multiquant_table.R` (not the same file as
`tr2/01i_nn_d500_quant.R`), `DELETED_DUPLICATE_TABLES.md` and
`INSTALLED_REGEN999_TABLES.md`. Renamed:

| Old name | New name |
|---|---|
| `revision_experiments/79_rerun_sun_smin_fixed_tables.R` | `tr1/78b_rerun_sun_smin_fixed_tables.R` |
| `tr1/92_validate_incremental_radi.R` | `tr1/92b_validate_incremental_radi.R` (output `results/tr1/wp6_incr/92b_validate.csv`) |

So `tr1/79_harness_guard_tests.R` and `tr1/92_wp3_stats.R` are now the only
`79_*` and `92_*` files. New: `tr1/95_wp6_incr_summary.R`. The WP2 console
logs (`revision_experiments/logs/`) left git and now sit in
`results/tr1/_logs/`, which is ignored. Paths inside the scripts were updated
to match; no computation changed.

## Running

Scripts are run from the repo clone root (paths use `here::here()`). R
package prerequisites are checked by `tr2/00_env_check.R`. On this project's
Windows setup, `Rscript` must be launched from PowerShell, not Git-Bash/MSYS
(known segfault). The Python steps (`tr2/02b_convert_highdim.py`,
`tr2/05_wp4_runtime_pyod.py`, `tr2/08_wp6_pyod_baselines.py`) use the bundled
`revision_experiments/.venv` (interpreter at `.venv/python.exe`).

The tr1 scripts follow the same rules: they run from the clone root, their R
prerequisites are checked by `tr2/00_env_check.R`, and the four tr1 Python
steps (`tr1/81_wp4_baselines.py`, `tr1/84a_wp5_fetch_convert.py`,
`tr1/90b_wp7_hdbscan.py`, `tr1/91b_wp6_runtime_py.py`) use the same `.venv`.

Reruns. The drivers write one row per finished cell (`append_result()`) and
skip every cell already in their output (`has_result()`, or a keys-only
`*_done.csv` file where one cell spans several rows). Rerunning a driver
against a committed result CSV therefore skips the finished cells and
computes nothing new. To recompute, move the file aside first.

Order of the tr1 steps that depend on each other:

- Run `tr1/80_export_datasets_wp4.R` before `81`, `82` and `92`. It writes
  the data CSVs they read, which are not committed.
- WP5: `84a`, then `84`, then `81` with
  `--data-dir revision_experiments/results/tr1/wp5/data` and
  `--out-dir revision_experiments/results/tr1/wp5` (paths relative to the
  working directory), then `84b`.
- `90_wp7_clustering_quality.R` reads the WP8 cluster files uncompressed.
  Run `gunzip -k revision_experiments/results/tr1/wp8/*_clusters.csv.gz`
  first.
- `18_table_consistency.R` and `39_verify_manuscript_tables.R` read the
  manuscript sources, which are not in this repository. Set
  `TR1_MANUSCRIPT_DIR` to the folder that holds them (default: a sibling
  folder `TR1_Neurocomputing_resubmit`).

Environment variables read by tr1 scripts:

| Variable | Read by | Effect |
|---|---|---|
| `WP0_GATE_OUT_CSV` | `13_wp0_gate.R`; set by both `regen_wilt_*_launcher.ps1` | Output CSV (default `results/tr1/wp0_gate_v2.csv`). |
| `WP0_GATE_TIMEOUT_SEC` | `13_wp0_gate.R`; set by both launchers | Per-cell timeout in seconds (default 600; the launchers set 5400). |
| `REGEN_FINAL_OUT_CSV` | `28_regen_final.R` | Output CSV (default `results/tr1/regen_final.csv`). |
| `FINAL_BASE_OUT_CSV` | `33_final_baselines.R` | Output CSV (default `results/tr1/final_baselines.csv`). |
| `WP2A_SAT_SUFFIX` | `49_wp2a_rk_saturation.R` | Suffix added to its two saturation CSVs. |
| `WP5_HIGHD_TIMEOUT_SEC` | `84_wp5_highd.R` | Per-cell timeout in seconds (default 540). |
| `TR1_MANUSCRIPT_DIR` | `18_table_consistency.R`, `39_verify_manuscript_tables.R` | Folder with the manuscript sources. |
| `RSCRIPT` | both launchers | Rscript executable (default `Rscript.exe`). |
| `CCD_QUANTILE_*` | every loader of a quantile table | See "Quantile tables" below. |

Several committed result files were written with an output variable set:
`wp0_gate_v3.csv` and `regen_proposed_*.csv` by `13`, and
`regen_final_<set>.csv` and `wp2c_rerun_alpha001.csv` by `28`.

## shared/

- **`harness.R`** -- the 9(+)-method registry, `evaluate()`,
  `load_real_dataset()`, `get_simul()`, checkpointing
  (`has_result()`/`append_result()`), and the manuscript's alpha schedule
  (`rk_quant_label_paper()`, `nn_quant_label_paper_UN()`,
  `nn_quant_label_paper_SUN()`, the only definitions in the repository).
  tr1's `wp0_mccd_methods.R` *appends* the four MCCD detectors to its
  registries. It reads `R/`, `methods/`, `simulations/` and `data/` from the
  repo root and writes nothing there except generated quantile tables (see
  "Quantile tables"), so it works unchanged from either `tr1/` or `tr2/`
  callers.
- **`get_simul(variant, d, quant, n)`** needs an explicit level (since
  2026-09-05): pass `rk_quant_label_paper(d)` for U-MCCD/SU-MCCD and the RK
  wrappers, `nn_quant_label_paper_UN(d)` for UN-MCCD and the UNCCD wrappers,
  or `nn_quant_label_paper_SUN(d)` for SUN-MCCD. The tr2 scripts pass
  `rk_quant_label_paper(d)` / `nn_quant_label_paper_UN(d)` (one passes the
  literal `"999"`). A table that is missing, or shorter than `n`, is
  generated on the spot through `R/ccds/quantile_table.R`. With
  `CCD_QUANTILE_GENERATE=FALSE` (or `options(ccd.quantile.generate = FALSE)`)
  it stops with an error instead, as before.

### Harness contract (unchanged by this reorg)

- Every method in `METHOD_REGISTRY` has the signature
  `function(X, d, Y = NULL, ...)` and returns
  `list(score = <numeric vector>, t_construct = <seconds>, t_total = <seconds>, ...)`.
- Labels: `Y == 1` is regular, `Y == 0` is an outlier.
- Always score via `evaluate(Y, score, threshold)` -- never via the
  positional `count_scores2()`/`count_DBSCAN()`/`count_MST2()`/`count_ODIN()`
  helpers used by the original (pre-revision) drivers. `evaluate()` reorders
  `(Y, score)` jointly before counting; the positional helpers slice the
  first `n-n0` / last `n0` rows without checking `Y`, which is silently wrong
  whenever a loader's final `order(...)` call sorts by the wrong column (see
  `tr1/26_loader_sort_audit.R`: glass and ecoli are MIS-SORTED this way).

## Quantile tables

The RK-based and NND-based spatial-randomness tests read their critical
values from Monte Carlo tables, `R/RK-test_quantile/RK-test-simul_<d>d_<level>%.RData`
and `R/NN-test_quantile/NN-test-simul_<d>d_<level>%.RData`. No table is
published for tr1: the RK tables are 100-700 MB each. The only tables in the
repository are five NN tables committed for tr2, all at the 0.1% level:
`NN-test-simul_166d_999%.RData`, `NN-test-simul_274d_999%.RData`,
`NN-test-simul_400d_999%.RData`, and the spliced
`NN-test-simul_166d_999%_n3062_spliced.RData` and
`NN-test-simul_400d_999%_n3686_spliced.RData`.

`R/ccds/quantile_table.R` (read its header) supplies a missing table:

- A table file that exists is loaded as before, and a shipped table is
  used exactly as before, so published results do not change when the
  tables are present.
- `get_simul()` in `shared/harness.R` generates a missing table at the data
  size `n`, and completes a shipped table shorter than `n` with generated
  entries past its end.
- The 1197 original drivers under `simulations/` load tables through
  `load_quantile_table()`, a drop-in for `load()`. When the file is missing,
  `simul` becomes a placeholder that records the variant, d and level, and
  `nnccd.radi()` (`R/ccds/UN_CCD.R`) or `ccd.Kest.edge.quantile()`
  (`R/ccds/RK_CCD_New.R`) fills it at the data's own size on first use.
- `ccd.Kest.edge.quantile()` grows an RK table only when its search reaches
  a row past the end of the table (where it used to stop with "subscript out
  of bounds"): first to `ccd.quantile.n` points, then to all `n`.
- The estimator is the original one: the per-iteration bodies are
  transcriptions of the Monte Carlo loops in `NNDestP.simpois.lower.quant()`
  (`R/ccds/NN_Dist_Est.R`) and `KestP.simpois.edge.quantile()`
  (`R/ccds/Kest.R`), with the same `quantile()` reductions. For RK at
  d >= 342 it uses the log-space body of `tr2/01_gen_quantile_table.R`.
- Generated tables are cached as
  `R/<RK|NN>-test_quantile/generated/<RK|NN>-test-simul_<d>d_<level>%_n<n>.RData`
  (ignored by git) and reused by any later call that needs no more than `n`
  points. One process generates a given table at a time: parallel workers
  and concurrent jobs wait for it (a `.lock_*` folder in the cache) and then
  reuse its file.

Settings (R option, else environment variable, else default):

| R option | Environment variable | Default | Meaning |
|---|---|---|---|
| `ccd.quantile.niter` | `CCD_QUANTILE_NITER` | as the shipped tables: NN 10000; RK 2000 for d < 50, 10000 for d >= 50 | Monte Carlo iterations. |
| `ccd.quantile.n` | `CCD_QUANTILE_N` | 1000 | Points covered when the data size is not known; first extent of a grown RK table. |
| `ccd.quantile.cores` | `CCD_QUANTILE_CORES` | all cores - 1 (1 inside a parallel worker) | Worker processes. |
| `ccd.quantile.seed` | `CCD_QUANTILE_SEED` | 20260927, plus d | Base seed. |
| `ccd.quantile.cache` | `CCD_QUANTILE_CACHE` | `R/<RK\|NN>-test_quantile/generated/` | Cache folder. |
| `ccd.quantile.generate` | `CCD_QUANTILE_GENERATE` | TRUE | FALSE stops with an error instead of generating. |

Cost. Generation is the expensive part of a first run. Per iteration, an NN
table of n points costs roughly n^3 d operations and an RK table roughly n^2
numerical integrals, so a table for a few hundred points takes minutes on 16
cores and one for several thousand points can take many hours at the default
iterations. For a quick run, lower `CCD_QUANTILE_NITER`; the 0.1% quantiles
are the most affected (see `tr1/73_niter_sensitivity.R`). Generated tables
are cached, so the cost is paid once per table.

Monte Carlo error. A generated table reproduces a published number up to
Monte Carlo error, not bit for bit, because its draws differ from those of
the tables the papers used. A generated RK table also leaves out the raw
draw matrix `Kest.m` (the detectors read only `$quan` and `$r`), so
`tr1/42`, `tr1/46` and `tr1/49`, which re-reduce `Kest.m` at other levels,
need the original RK tables. To use original tables, put them in
`R/RK-test_quantile/` or `R/NN-test_quantile/` under their original names;
a file on disk is always loaded as it is.

Several tr1 WP2 scripts read tables from their own folders
(`R/NN-test_quantile_n200/`, `R/NN-test_quantile_d21_regen/`,
`R/NN-test_quantile_provenance/`, `R/NN-test_quantile_regen999/`), written
by `48`, `56`, `58`, `64`, `71` and `74`. These folders are not published.

## .venv/

The pinned PyOD 3.6.1 / torch 2.13.0+cpu Python environment, used by tr2's
`05_wp4_runtime_pyod.py` and `08_wp6_pyod_baselines.py`, and (since
2026-09-05) by tr1's `81_wp4_baselines.py`, `84a_wp5_fetch_convert.py`,
`90b_wp7_hdbscan.py` and `91b_wp6_runtime_py.py`. Unmoved by this reorg
(still `revision_experiments/.venv`).

## results/

- `results/datasets_csv/` -- 16 exported real-data CSVs + `manifest.csv`
  (two of the 16, `Musk_sub1000.csv` and `Speech_sub1000.csv`, come from the
  WP5 subsample path, `07_wp5_subsample_ccd.R`, not the main loader). Written
  by tr2; read by tr2 scripts that need a dataset as a flat CSV rather than
  through `load_real_dataset()`. tr1's regeneration pipeline does not read
  this directory -- it loads data only via `load_real_dataset()`, per
  `tr1/REGEN_SPEC.md`; the one exception is tr1's WP5 step
  `84a_wp5_fetch_convert.py`, which copies musk/arrhythmia from here into
  `results/tr1/wp5/data/`. Not moved or touched by this reorg.
- `results/scores_cache/` -- cached `.rds` score vectors keyed by
  `<dataset>_<method>.rds`. Shared in principle, written by tr2's
  `06_wp5_highdim.R` / `07_wp5_subsample_ccd.R` / `07b_wp5_fulldata_ccd.R`
  and read by `10_wp3_real.R`. Not moved or touched by this reorg.
- `results/tr1/` -- tr1's result files. Layout:
  - Flat files: the WP0 gate and real-data regeneration (`wp0_gate_v2.csv`, `wp0_gate_v3.csv`,
    `regen_*.csv`, `final_baselines.csv`, `final_comparison.csv`,
    `lof_original8.csv`, `waveform_*.csv`, `diag_waveform_ecoli.csv`,
    `dataset_inventory.csv`, `row_order_invariance.csv`), and the WP2 files
    (`wp2a_*`, `wp2b_*`, `wp2c_*`: CSVs and `*_findings.md`).
    `row_order_invariance.csv` has no writer among the present scripts.
  - `wp3/` .. `wp9/`: one folder per package (see the tr1 tables below).
    `wp6_incr/` holds the incremental radius search checks and the paired
    UN-/SUN-MCCD timings (`91` with `--resdir=wp6_incr`, `92b`-`96`).
  - `ns_sun_d10/` (`97`) and `su_driver_fix/` (`98`): the simulation-driver
    checks.
  - `cost_probe/` (`42 --cost`) and `_logs/` (console logs of the WP2
    background jobs): not committed.

  `results/` is gitignored, as for tr2; the files that are committed were
  force-added (`git ls-files revision_experiments/results/tr1` lists them).
  Committed: result tables and summaries, findings files (`*_findings.md`,
  `WP*_FINDINGS*.md`, `wp3/sim_variability_INVENTORY.md`), verification
  records (`*_verify*.csv`, `wp6_incr/92b_validate.csv`,
  `wp6_incr/93_validate_extended.csv`), raw per-cell CSVs, the WP4 and WP5
  score vectors (`wp4/scores/`, `wp5/scores/`) with their `fit_log.csv` and
  `versions.txt`, the converted WP5 data
  (`wp5/data/*.csv`), the WP6 timing inputs (`wp6/data/*.csv`, the rep-1
  data sets `91b` timed), the WP8 cluster files gzipped
  (`wp8/*_clusters.csv.gz`), and one console log
  (`wp6_incr/94_ab_run.log`). Not committed: console logs (`*.log`, `*.err`,
  `*.DONE`, `_logs/`, `logs/`), checkpoint files (`*_done.csv`), smoke-test
  output (`smoke/`), scratch files (`wp7/tmp/`, `wp6/GO`), data a script
  rebuilds (`wp4/data/<set>.csv` from `80`, `wp5/data/raw/` from `84a`,
  `wp6_incr/data/` from `91`), the uncompressed WP8 cluster files, and
  superseded or abandoned runs (`wp0_gate.csv`, `final_tables.csv`, `wp2a_rk_saturation*_v2.csv`,
  `wp2a_rk_saturation_v1_incomplete_abandoned.csv`,
  `wp5/letter_with_duplicates/`).

  Never run `git clean -x` or `git clean -X` in a clone that holds
  unpublished results. Both delete ignored files, which here means every
  uncommitted result, the `.venv`, and every quantile table (`*.RData` is
  ignored too).
- `results/tr2/` -- tr2's result files (csv, `figures/`, `probes/`,
  `wp4_data/`, `wp4_data2/`, `wp6_scores/`). CSVs and `figures/` stay directly
  under `results/tr2/`.
  - **`results/tr2/_logs/`** (new in this reorg) -- every `*.log`, `*.flag`
    and `*DONE*` file that had accumulated directly under `results/tr2/`
    (run logs, shutdown-watcher logs, completion flags) was moved here to
    de-clutter the results directory. Nothing was deleted. CSVs, `.RData`,
    `.err`, `.lock`, `.bak*` files and the subdirectories above were left in
    place. A handful of `*_log.txt` text logs (e.g. `wp4_run_log.txt`,
    `wp3_real_log.txt`) were not moved and remain directly under
    `results/tr2/`.

## tr1/

NEUCOM-D-26-15191 revision experiments. Every script runs from the clone
root, e.g. `Rscript revision_experiments/tr1/28_regen_final.R`. Purpose and
main output are condensed from each script's own header comment. Paths in
the Main output column are under `results/tr1/` unless they start with `R/`
or `tr1/`. "(not committed)" marks files that are not in git; "console only"
means the script prints and writes no file.

### WP0: real-data regeneration pipeline

| Script | Purpose | Main output |
|---|---|---|
| `wp0_mccd_methods.R` | Wires the four MCCD detectors (U-MCCD, SU-MCCD, UN-MCCD, SUN-MCCD) into the harness by appending to `METHOD_REGISTRY` and `REAL_DATA_THRESHOLDS`. Source it after `shared/harness.R`. | n/a (sourced) |
| `published_realdata_truth.csv` | Input: published TPR/TNR/BA/F2 of the nine methods on the eight original real data sets, transcribed from the submitted supplement tables (288 rows). | n/a (input) |
| `published_datasets_truth.csv` | Input: published n, d and outlier count of the eight original real data sets. No present script reads it. | n/a (input) |
| `13_wp0_gate.R` | Runs the four MCCD detectors on real data sets and compares TPR/TNR/BA/F2 with `published_realdata_truth.csv` to 3 decimals; sweeps the S_min readings for SU-MCCD and SUN-MCCD. | `wp0_gate_v2.csv` (default, not committed); through `WP0_GATE_OUT_CSV`: `wp0_gate_v3.csv`, `regen_proposed_small.csv`, `regen_proposed_mid.csv` |
| `regen_wilt_nn_launcher.ps1`, `regen_wilt_rk_launcher.ps1` | Run `13` on wilt with a 5400 s cell timeout, one for the NND-based pair (UN-MCCD, SUN-MCCD) and one for the RK-based pair (U-MCCD, SU-MCCD). | `regen_proposed_wilt_nn.csv`, `regen_proposed_wilt_rk.csv`; `regen_wilt_*.log` (not committed) |
| `14_wp0_probe.R` | Two probes on small data sets: SUN-MCCD at the 1% NN table instead of 0.1%, and a constant S_min against the contamination-based readings. | console only |
| `15_wp0_constant_smin.R` | Extends the constant-S_min probe to hepatitis, ecoli and vowels over a small grid. | console only |
| `16_wp0_su_highsmin.R` | Tests whether SU-MCCD on glass and stamps reproduces the published rows at S_min above 0.0625. | console only |
| `17_dataset_audit.R` | Prints n, d, n0 and contamination of the eight original real data sets as `load_real_dataset()` returns them. | console only |
| `18_table_consistency.R` | Checks BA = (TPR + TNR)/2 over the two supplement result tables (72 method and data set pairs), reading the LaTeX source from `TR1_MANUSCRIPT_DIR`. | console only |
| `19_baseline_recheck.R` | Reruns three internally inconsistent baseline rows (DBSCAN on wilt, LOF on ecoli and on vertebral). | console only |
| `20_scale_audit.R` | Reports the feature scale of each data set after the loader's preprocessing. | console only |
| `21_regen_smin_grid.R` | S_min sweep over {0.01, 0.02, 0.03, 0.05, 0.0625, 0.10, 0.15, 0.20} for SU-MCCD and SUN-MCCD on hepatitis, glass, vertebral, ecoli and stamps. | `regen_smin_grid.csv` |
| `23_regen_baselines.R` | Label-free DBSCAN, MST, ODIN and iForest rows for the eight original data sets, plus reference rows that reproduce the published drivers (DBSCAN given the label column, oracle eps, per-set MST threshold). One data set per call. | `regen_baselines.csv` |
| `24_diag_waveform_ecoli.R` | Diagnoses why waveform and ecoli did not reproduce the published rows; one task per call. | `diag_waveform_ecoli.csv` |
| `25_smin_match.R` | Finds which S_min values of `21`'s grid reproduce the published rows. | console only |
| `26_loader_sort_audit.R` | Measures the loader's mis-sort of glass and ecoli (sorted by a feature column instead of the label). | console only |
| `27_positional_vs_correct.R` | Scores the same detector output on glass and ecoli both correctly (`evaluate()`) and positionally (`count_scores2()`), and compares both with the published rows. | console only |
| `28_regen_final.R` | Runs the four MCCD detectors at one fixed S_min (third argument, default 0.0625). | `regen_final.csv` (default, not committed); through `REGEN_FINAL_OUT_CSV`: `regen_final_{small,vowels,waveform,wilt,new_small,pendigits,thyroid,pageblocks}.csv`, and `wp2c_rerun_alpha001.csv` (the six cells rerun with the repaired 0.1% NN tables, S_min 0.05) |
| `29_waveform_alpha.R` | Runs UN-MCCD and SUN-MCCD on waveform at the 1% NN table, against the published row and the 0.1% run. | `waveform_alpha.csv` |
| `30_waveform_scaling.R` | Runs U-MCCD and SUN-MCCD on waveform under raw, z-score and median/MADN scaling. | `waveform_scaling.csv` |
| `31_headline_impact.R` | Compares the regenerated and published rows of the original eight data sets: change per method, and the winner per data set. | console only |
| `32_dataset_inventory.R` | Structural inventory (n, d, contamination, sort correctness) of every data set the loader builds; runs no detector. | `dataset_inventory.csv` |
| `33_final_baselines.R` | The five baselines at the configuration fixed in `BENCHMARK_EXPANSION_RULE.md`, on the eight added data sets. | `final_baselines.csv` |
| `34_final_summary.R` | Assembles the real-data comparison of 16 data sets by 9 methods. | `final_comparison.csv` |
| `35_print_tables.R` | Prints `final_comparison.csv` in the row and column order of the manuscript tables. | console only |
| `37_arith_audit.R` | Checks BA = (TPR + TNR)/2 and rebuilds F2 from TPR, TNR, n and n0 for every row of the 16 by 9 table. | console only |
| `38_lof_original8.R` | LOF on the original eight data sets with the published settings, scored with `evaluate()`. | `lof_original8.csv` |
| `39_verify_manuscript_tables.R` | Parses the real-data tables back out of the LaTeX source (`TR1_MANUSCRIPT_DIR`) and compares them with `final_comparison.csv`. | console only |

### NN multi-quantile generator

| Script | Purpose | Main output |
|---|---|---|
| `01i_nn_multiquant_table.R` | NN quantile-table generator that runs the Monte Carlo once and reduces the same draws at several levels. Writes one drop-in table per level plus the raw draws. Its default engine reproduces the original generator's draws bit for bit when run serially under the same seed. Sourced by `48`, `58`, `64`, `65`, `71`, `73` and `74`. | tables and `NN-draws_<d>d_n<n>.RData` in the folder given as its fifth argument (not committed) |

### WP2a: alpha, tolerance, search direction and S_min sensitivity

| Script | Purpose | Main output |
|---|---|---|
| `40_wp2a_tolerance_sweep.R` | Sweeps the tolerance of the KS-CCD connectivity search in `connected.ksccd.m()` (fixed at 0.01 in the code) through a copy of the function, and records the recovered delta and the flagged set per cell; `--summary`. | `wp2a_tolerance_sweep.csv`, `wp2a_tolerance_summary.csv` |
| `41_verify_tolerance_sweep.R` | Checks `40`: recomputes stability from the raw rows, compares the tol = 0.01 rows with `final_comparison.csv`, and reruns with the unmodified function. | console only |
| `42_wp2a_alpha_sweep.R` | Alpha sweep of the four MCCD detectors on the 16 real data sets. RK levels are derived from the raw draws `Kest.m` stored in the original RK tables; NN levels are limited to the tables on disk. Modes `--inventory`, `--summary`, `--cost`. | `wp2a_alpha_sweep.csv`, `wp2a_alpha_summary.csv`, `wp2a_alpha_findings.md`; `cost_probe/` (not committed) |
| `43_wp2a_direction_sweep.R` | Ascending against descending radius search for UN-MCCD and SUN-MCCD on the real data sets. | `wp2a_direction_sweep.csv`, `wp2a_direction_findings.md` |
| `44_wp2a_smin_sweep16.R` | S_min sweep for SU-MCCD and SUN-MCCD on all 16 real data sets over a 10-point grid, in chunks with a time budget; `--fill-not-run` records the cells not attempted. | `wp2a_smin_sweep16.csv` |
| `46_verify_wp2a_alpha.R` | Checks `42`: which level the duplicated NN tables follow, that RK `$quan` equals the quantile of `Kest.m`, and the effect sizes from the raw rows. | console only |
| `47_alpha_operating_range.R` | Splits `42`'s alpha effect into the full swept range and the operating range (alpha 1% and 0.1%). | `wp2a_alpha_operating_range.csv` |
| `48_wp2a_nn_saturation.R` | Tests whether the NND test stops responding to alpha at high d, with NN tables built by `01i` at n = 200 at the same four levels for d in {2, 3, 5, 10, 20, 50}; modes `--validate`, `--crosscheck`, `--summary`. | `wp2a_nn_saturation.csv`, `wp2a_nn_saturation_summary.csv`, `wp2a_nn_crosscheck.csv`, `wp2a_nn_crosscheck_tables.csv`, `wp2a_nn_cost.csv`, `wp2a_nn_zero_quantiles.csv`, `wp2a_saturation_findings.md`; tables in `R/NN-test_quantile_n200/` (not committed) |
| `49_wp2a_rk_saturation.R` | RK counterpart: zero-quantile fraction of every RK table from its raw draws, and U-MCCD and SU-MCCD at n = 200 for d from 2 to 100 at four alpha levels on uniform and Gaussian clusters. | `wp2a_rk_table_degeneracy.csv`, `wp2a_rk_saturation.csv`, `wp2a_rk_saturation_summary.csv`, `wp2a_rk_saturation_findings.md` |
| `50_wp2a_nn_alpha_bigsets.R` | Alpha sensitivity of UN-MCCD and SUN-MCCD on wilt (d = 5, levels 10% and 5%) and pageblocks (d = 10, 1% and 0.1%), the two real data sets whose two NN tables differ; refuses to run on identical tables. | `wp2a_nn_alpha_bigsets.csv` |
| `51_wp2a_smin_bigfour.R` | Completes `44`: pageblocks, thyroid, waveform and wilt on a reduced 4-point grid, and the missing PenDigits points. | `wp2a_smin_bigfour.csv` |
| `52_direction_summary.R` | Summarizes `43` from its per-cell CSV, after checking that the direction argument reached the radius search. | `wp2a_direction_summary.csv` |
| `53_verify_nn_saturation.R` | Checks `48`: table distinctness by content, design completeness, effect sizes, and the replicate noise floor. | console only |
| `56_verify_d21_duplication.R` | Tests at d = 21 whether the duplicate NN tables came from the generator or from a file copy: one set of draws reduced at two levels, and two independent runs. | `wp2a_d21_duplication.csv`; tables in `R/NN-test_quantile_d21_regen/` (not committed) |
| `57_waveform_d21_rerun.R` | Reruns UN-MCCD and SUN-MCCD on waveform (d = 21) with a new same-draws 1%/0.1% NN pair and with the shipped table. | `wp2a_waveform_d21_rerun.csv` |
| `58_gen_d21_chunked.R` | Generates the d = 21, n = 5000 same-draws NN pair in resumable chunks of raw draws, reduced once after pooling. | `R/NN-test_quantile_d21_regen/samedraws/` (not committed) |

### WP2b: stability

| Script | Purpose | Main output |
|---|---|---|
| `45_wp2b_seed_stability.R` | Stability of the four MCCD detectors under the R random-number stream (`--determinism`) and under small data perturbations (jitter of 0.01 SD, leave 1% out; 50 reps), at S_min = 0.05; `--summary`. | `wp2b_determinism_smin005.csv`, `wp2b_seed_stability_smin005.csv`, `wp2b_stability_summary_smin005.csv`; the same names without `_smin005` hold the earlier run at S_min = 0.0625 |

### WP2c: label-free S_min

| Script | Purpose | Main output |
|---|---|---|
| `54_wp2c_smin_rule.R` | Real-data side. The detectors see S_min only through k = round(S_min * n), so plateaus are reported in k and only new k values are run; evaluates the candidate label-free rules. | `wp2c_plateau_table.csv`, `wp2c_rule_candidates.csv`, `wp2c_extra_cells.csv`, `wp2c_adaptive_real.csv` |
| `55_wp2c_simulation_arm.R` | Simulation side: SU-MCCD and SUN-MCCD under the contamination-based S_min rule and each label-free candidate, paired on the same replicates (uniform and Gaussian clusters, d = 5 and 10, four contamination levels, n = 500, 100 reps). | `wp2c_simulation_arm.csv`, `wp2c_simulation_summary.csv` |
| `59_verify_wp2c.R` | Checks `55` from the raw rows: pairing, paired differences, single-cluster collapse by rule, worst cell. | console only |
| `60_wp2c_smin_zero_real.R` | S_min = 0 (the detectors' own default) on the 11 real data sets with n < 1500. | `wp2c_smin_zero_real.csv` |
| `61_aggregate_impact_smin005.R` | Lists the printed aggregate values that change when the real-data S_min moves from 0.0625 to 0.05 (31 of 32 cells are identical). | console only |
| `78b_rerun_sun_smin_fixed_tables.R` | Reruns the SUN-MCCD S_min cells with the repaired 0.1% NN tables: the 11 small sets at S_min in {0, 0.03, 0.05, 0.0625}, and waveform and PenDigits at {0.01, 0.05, 0.0625, 0.15, 0.30}. Supersedes the SUN-MCCD column of `wp2c_smin_zero_real.csv` and the SUN-MCCD waveform and PenDigits cells of `44` and `51` at those values. | `wp2c_sun_smin_rerun_fixed.csv` |

### Alpha provenance and NN table repair

| Script | Purpose | Main output |
|---|---|---|
| `62_audit_nn_duplicates.R` | Audits NN table duplication per dimension by content (`$average`, `$median`), reporting every dimension. | console only |
| `63_nn_provenance_recon.R` | For each affected dimension: the level token each method loads, the shipped table's extent, and whether the two shipped files differ. | `wp2a_provenance_targets.csv` |
| `64_nn_provenance_gen.R` | Generates replicate clouds of fresh NN tables at the candidate levels for each duplicated dimension (one Monte Carlo run reduced at every level). | `wp2a_provenance_gen_manifest.csv`; tables in `R/NN-test_quantile_provenance/` (not committed) |
| `65_nn_provenance_test.R` | Which level each shipped duplicated table holds: pointwise sign test against the replicate clouds, positive controls at d = 5, 10 and 20, and a behavioural check. | `wp2a_provenance_signtest.csv`, `wp2a_provenance_control.csv`, `wp2a_provenance_behavioural.csv`, `wp2a_provenance_behavioural_summary.csv` (`wp2a_provenance_findings.md` is written by hand from these) |
| `66_alpha_claim_audit.R` | Checks every alpha-level statement of the manuscript against the table each detector loads and whether that file differs from its sibling; classes each statement as verified, unverifiable or false. | `wp2a_alpha_claim_audit.csv` |
| `67_verify_provenance.R` | Re-derives `65`'s verdicts with a direct level estimate (the share of fresh draws at or below each table entry, median over subsample sizes), with six control files of known level. | `wp2a_provenance_verify.csv` |
| `68_verify_claim_audit.R` | Checks `66`: RK table distinctness by content at the real-data dimensions, and the simulation grid, calling the live resolvers. | `wp2a_verify_rk_distinct.csv`, `wp2a_verify_sim_grid.csv` |
| `69_verify_table_inventory.R` | Checks that every table the resolvers can ask for exists, and that restored files keep the known duplicate pattern. | console only |
| `70_remove_duplicate_tables.R` | Deletes the redundant member of each duplicated NN pair where no resolver loads it (dry run by default, `--apply` to delete). | `tr1/DELETED_DUPLICATE_TABLES.md` |
| `71_regen_nnd_alpha001.R` | Generates genuine 0.1% NN tables at d = 12, 16, 18, 19 and 21, each to its data set's own n, as chunks of raw draws. | `R/NN-test_quantile_regen999/` (not committed) |
| `72_verify_regen999.R` | Acceptance test of the staged 0.1% tables: the new 1% table against the shipped one, the ordering of the levels, and a direct level read against independent draws from `64`. | `wp2c_regen999_verify.csv` |
| `73_niter_sensitivity.R` | How far the 0.1% table at d = 12 moves when the number of Monte Carlo iterations is cut (random subsets of 10000). | `wp2c_niter_sensitivity.csv` |
| `74_pool_partial.R` | Builds a table from the chunks of a `71` run that exist so far; the pooled iteration count goes into the metadata and the file name. | `R/NN-test_quantile_regen999/` (not committed) |
| `75_install_regen999.R` | Installs verified 0.1% tables into `R/NN-test_quantile/`, keeping the shipped content under the name of the level it holds (dry run by default, `--apply`). | `tr1/INSTALLED_REGEN999_TABLES.md` |
| `76_aggregate_after_alpha_fix.R` | Recomputes the aggregate table after the six rerun cells (SUN-MCCD on vowels, PenDigits, lymphography, hepatitis and waveform; UN-MCCD on waveform). | `wp2c_aggregate_after_fix.csv` |
| `77_purge_wrong_tables.R` | Deletes the remaining duplicate or mislabelled NN tables: the 0.1% member of identical pairs at d = 11, 13, 14, 15 and 17, and the 0.1%-named single files holding 1% at d = 22-28 (dry run by default, `--apply`). | `tr1/DELETED_DUPLICATE_TABLES.md` |
| `78_manuscript_tables_after_fix.R` | Rebuilds the two main-text real-data tables (per data set, and aggregate) from the corrected per-cell table; ranks and winners use the rounded values. | `wp2c_manuscript_tables.csv` |

The two manifests `DELETED_DUPLICATE_TABLES.md` and
`INSTALLED_REGEN999_TABLES.md` are the only record in git of which table
files were deleted, renamed or installed (`*.RData` is ignored).

### Harness guard tests

| Script | Purpose | Main output |
|---|---|---|
| `79_harness_guard_tests.R` | Guard tests for five harness defects found on 2026-09-05 (two competing alpha schedules, tables shorter than the data, NA or short score vectors, appends that ignore the header, truncated rows read as complete), and for on-the-spot table generation: the generator matches the original Monte Carlo under the same seed, the caller's random-number stream is left untouched, missing tables are generated, cached and reused, a short shipped table keeps its entries and is completed, a missing file becomes a placeholder that `nnccd.radi()` fills at its own level, and `ccd.Kest.edge.quantile()` grows a short RK table only when its search needs more rows. | console only (temporary files under `tempdir()`) |

### WP4: modern baselines on the real data

| Script | Purpose | Main output |
|---|---|---|
| `80_export_datasets_wp4.R` | Exports the 16 real data sets, exactly as the R loader prepares them, to CSV for the Python competitors; stops if (n, d, outliers) differ from the manuscript's data table. | `wp4/data/<set>.csv` (not committed), `wp4/data/manifest.csv` |
| `81_wp4_baselines.py` | Eight competitors on the 16 sets: ECOD, COPOD, HDBSCAN (GLOSH), OPTICS, mutual-kNN, SNN, DIF and LUNAR (DIF and LUNAR with seeds 1-5). Mutual-kNN and SNN run at k in {5, 10, 15, 20, 30} and at k = round(sqrt(n)). Writes raw scores and native labels only, one cell at a time under a time budget (`--budget`); `--data-dir` and `--out-dir` reuse it for WP5. | `wp4/scores/<set>_<method>.csv` (and `_labels.csv`), `wp4/fit_log.csv`, `wp4/versions.txt`; `wp4/fit_errors.log` (not committed) |
| `82_wp4_metrics.R` | Turns `81`'s scores into TPR/TNR/BA/F2 with the harness `evaluate()` under PyOD's strict threshold rule, at contamination 0.1 and at the true rate, and rebuilds the two main-text real-data tables with all 17 methods. Mutual-kNN and SNN use k = round(sqrt(n)); the k chosen on the labels gives upper-bound rows. | `wp4/wp4_metrics_long.csv`, `wp4/wp4_metrics_main.csv`, `wp4/WP4_FINDINGS.md` |
| `83_wp4_latex_rows.R` | Prints the LaTeX rows of the WP4 tables from `82`'s output, rebuilding the integer counts so that third decimals match exactly. | console only |

### WP5: real data above d = 21

| Script | Purpose | Main output |
|---|---|---|
| `84a_wp5_fetch_convert.py` | Downloads letter and mnist from the ADBench mirror, checks (n, d, outliers), flips the label convention, drops the two duplicate row pairs of letter, scales letter and mnist by median/MADN, draws mnist's n = 1000 subsample (seed 20260905), and copies musk and arrhythmia from `results/datasets_csv/`. | `wp5/data/{letter,mnist,musk,arrhythmia}.csv`, `wp5/data/manifest.csv`; `wp5/data/raw/*.npz` (not committed) |
| `84_wp5_highd.R` | The four MCCD detectors and five baselines on letter (d = 32), mnist (d = 100), musk (d = 166) and arrhythmia (d = 274). U-MCCD and SU-MCCD are recorded n/a at musk and arrhythmia, where no RK table exists (checked by file existence); the RK zero-quantile fraction is recorded at each d. | `wp5/wp5_highd_results.csv`, `wp5/rk_degeneracy.csv`; `wp5/wp5_highd_done.csv` (not committed) |
| `84b_wp5_metrics.R` | Scores the eight competitors (`81` run on the WP5 folders) with `82`'s thresholding, and merges them with `84`'s rows into one table of 17 methods. | `wp5/wp5_metrics_long.csv`, `wp5/wp5_metrics_main.csv`, `wp5/WP5_FINDINGS.md`; `81` writes `wp5/scores/`, `wp5/fit_log.csv`, `wp5/versions.txt` |

### WP8: hard cases

All four take `--smoke` (output in `wp8/smoke/`, not committed) and
`--summarize`. The nine methods are the four MCCD detectors and LOF, DBSCAN,
MST, ODIN and iForest.

| Script | Purpose | Main output |
|---|---|---|
| `85_wp8_boundary_fp.R` | Boundary false positives: flag rate of regular points by normalized distance from their cluster centre (10 equal-width and 10 equal-mass bins); uniform and Gaussian two-cluster data, d = 3 and 10, n = 200, 5% contamination, 100 reps, nine methods. | `wp8/85_boundary_fp.csv`, `wp8/85_boundary_fp_summary.csv`, `wp8/85_boundary_fp_clusters.csv.gz`; `wp8/85_boundary_fp_clusters.csv`, `wp8/85_boundary_fp_done.csv` (not committed) |
| `86_wp8_outlier_types.R` | Local (interior and shell), bridged and collective outliers placed in the same two-cluster base; d = 3 and 10, n = 200, 10 outliers, 100 reps, nine methods. | `wp8/86_outlier_types.csv`, `wp8/86_outlier_types_summary.csv`, `wp8/86_outlier_types_clusters.csv.gz` |
| `87_wp8_small_cluster.R` | A third genuine cluster of m in {2, 4, 6, 8, 10, 12, 15, 20} points out of n = 200 against S_min = 0.05; SU-MCCD and SUN-MCCD, with U-MCCD and UN-MCCD as controls; 50 reps for the RK-based methods, 100 for the NND-based ones. | `wp8/87_small_cluster.csv`, `wp8/87_small_cluster_summary.csv`, `wp8/87_small_cluster_clusters.csv.gz` |
| `88_wp8_csr_violation.R` | Flag rate on data without outliers, under complete spatial randomness and three departures from it (Beta(2,5) density gradient, isotropic and anisotropic Gaussian); d = 3 and 10, n = 200, 100 reps, nine methods plus an MST threshold sweep. | `wp8/88_csr_violation.csv`, `wp8/88_csr_violation_summary.csv`, `wp8/88_csr_violation_clusters.csv.gz` |

### WP9: SUN-MCCD ablation

| Script | Purpose | Main output |
|---|---|---|
| `wp9_sun_variants.R` | Ablation machinery: a copy of `nnccd.radi()` with three switches (NND statistic, centre-point removal, search direction), substituted through private closure environments so no global binding changes; the synthetic generators. Sourced by `89`. | n/a (sourced) |
| `89_wp9_ablation.R` | Ablation driver: 4 settings (uniform and Gaussian, d = 3 and 10) by 5 variants by 100 reps; `--smoke` also checks that the stock variant reproduces SUN-MCCD exactly; `--summarize`. | `wp9/wp9_ablation.csv`, `wp9/wp9_ablation_summary.csv`; `wp9/smoke/` (not committed). `wp9/wp9_paired_diff.csv` and `wp9/WP9_FINDINGS.md` are derived from `wp9_ablation.csv`; no script writes them |

### WP7: clustering quality

| Script | Purpose | Main output |
|---|---|---|
| `90_wp7_clustering_quality.R` | Scores the MCCD partitions in WP8's four `*_clusters.csv` files against the true clusters (ARI, NMI, AMI, estimated number of clusters, unassigned fraction), and runs DBSCAN and HDBSCAN on each cell's regenerated data for comparison. Runs no MCCD detector. Start it only after `85`-`88` have finished, and unzip their cluster files first. | `wp7/wp7_metrics_long.csv`, `wp7/wp7_summary_metrics.csv`, `wp7/wp7_khat_table.csv`; `wp7/tmp/`, `wp7/smoke/`, run logs (not committed) |
| `90b_wp7_hdbscan.py` | Fits HDBSCAN (WP4 defaults, `min_cluster_size = 5`) on one regenerated cell and writes its native labels; called by `90` through a CSV round trip. | `wp7/tmp/` (not committed) |

### WP6: runtime and memory

| Script | Purpose | Main output |
|---|---|---|
| `91_wp6_runtime.R` | Runtime and memory of the four MCCD detectors and five baselines on uniform two-cluster data (5% contamination): n in {100, 250, 500, 1000, 2000} at d = 10, and d in {5, 10, 50, 100} at n = 500; 10 reps, single thread. Refuses to start unless the machine is idle (`--idle-check`; `--force` overrides). Exports the rep-1 data sets for `91b`. `--summarize` writes per-cell medians and log-log slopes, and summarizes `91b`. Other options: `--grid=`, `--reps=`, `--methods=`, `--resdir=`, `--nn-direction=`, `--smoke`. | `wp6/91_wp6_runtime_raw.csv`, `wp6/91_wp6_runtime_n.csv`, `wp6/91_wp6_runtime_d.csv`, `wp6/91_wp6_slope.csv`, `wp6/91b_wp6_runtime_py_summary.csv`, `wp6/data/*.csv`; `wp6/91_wp6_runtime_done.csv` (not committed). With `--resdir=wp6_incr`: `wp6_incr/91_wp6_runtime_{raw,n,d}.csv` and `wp6_incr/91_wp6_slope.csv`, the UN-MCCD and SUN-MCCD re-timing of 2026-09-23, superseded by `94` |
| `91b_wp6_runtime_py.py` | Runtime and memory of the eight WP4 competitors on `91`'s exported data sets, single-threaded; times `fit()` only, memory as the `tracemalloc` peak. | `wp6/91b_wp6_runtime_py_raw.csv`; `wp6/91b_wp6_runtime_py_done.csv` (not committed) |

Follow-ups on the incremental nearest-neighbour radius search in
`R/ccds/UN_CCD.R` (2026-09-23); the earlier per-step recomputation is kept
as `R/ccds/UN_CCD_radi_recompute_reference.R`:

| Script | Purpose | Main output |
|---|---|---|
| `92b_validate_incremental_radi.R` | Checks that the incremental search returns exactly the radii of the per-step recomputation, in both directions, on synthetic and real data, and reports the speed-up. | `wp6_incr/92b_validate.csv` |
| `93_validate_incremental_extended.R` | Extended check: radii over 5 seeds, d = 2-100 and six generators, both directions; end-to-end UN-MCCD and SUN-MCCD output on the 16 real data sets and the four above d = 21. | `wp6_incr/93_validate_extended.csv` |
| `94_wp6_ab_old_vs_new.R` | Paired timing: each replicate times the per-step and the incremental search inside UN-MCCD and SUN-MCCD on the same data (cells and seeds of `91`, alternating order), and times the radius search alone. Supersedes the `wp6_incr/91_*` re-timing, which ran on another day. | `wp6_incr/94_ab_raw.csv`, `wp6_incr/94_ab_run.log` |
| `95_wp6_incr_summary.R` | UN-MCCD and SUN-MCCD runtime rows from `94`'s paired run (arm `impl == "new"`, `part == "detector"`) with `91`'s definitions: median and CV of time and memory per cell, and log-log slopes over the d = 10 sweep. Also summarizes `96`. | `wp6_incr/95_ab_runtime_n.csv`, `wp6_incr/95_ab_runtime_d.csv`, `wp6_incr/95_ab_slope.csv`, `wp6_incr/95_direction_summary.csv` |
| `96_wp6_direction_timing.R` | Ascending against descending radius search of UN-MCCD and SUN-MCCD at n = 500, d = 10 (cell and seeds of `91`), 10 reps, alternating order. | `wp6_incr/96_direction_raw.csv` |

Outcome: 42 + 420 radius cases and 66 end-to-end runs identical; paired timing (94) shows the radius search 44-51% faster at n = 2000 and the whole detector 2-5% faster at d = 10, n >= 500.

### WP3: paired tests, simulation variability, parameter provenance

| Script | Purpose | Main output |
|---|---|---|
| `92_wp3_stats.R` | Exact paired Wilcoxon signed-rank tests of SUN-MCCD against each other method on F2 and BA over the 16 real data sets, Holm-adjusted per block, and the same with U-MCCD, SU-MCCD and UN-MCCD as the primary method; an inventory of per-replicate simulation output; the parameter-provenance table. Sources `82`. | `wp3/wilcoxon_real.csv`, `wp3/sim_variability.csv`, `wp3/sim_variability_INVENTORY.md`, `wp3/parameter_provenance.csv` |

### Simulation-driver checks

| Script | Purpose | Main output |
|---|---|---|
| `97_sun_ns_d10_rerun.R` | Reruns SUN-MCCD in the three Neyman-Scott settings (Matern, Thomas, mixed) at d = 10 with a chosen NN level: `--level=99` reproduces the published drivers, which load the 1% table; `--level=999` uses the 0.1% table of the stated schedule. Everything else is copied from the drivers. | `ns_sun_d10/97_<process>_d<d>_<level>_<variant>.csv` |
| `98_su_mccd_driver_check.R` | Reruns one SU-MCCD two-cluster simulation cell with an explicit d and RK table, to check two groups of drivers that do not run what their file names say: the d = 3, n = 1000 drivers set d = 2, and the d = 50, n = 50 drivers load the 2d table, which has no 0.1% level. | `su_driver_fix/98_<setting>_d<d>_n<n>_<table>_<impl>.csv` |

### tr1 documents

- `REGEN_SPEC.md` -- rules and output schema for the real-data regeneration
  (`13`-`39`).
- `BENCHMARK_EXPANSION_RULE.md` -- inclusion rule for the 16 real data sets
  and the baseline configuration, fixed before any result was computed.
- `REGENERATION_REPORT.md` -- report of the 16 by 9 regeneration as of
  2026-08-09, with a status note on the parts superseded since.
- `WP3_PROTOCOL.md` .. `WP9_PROTOCOL.md` -- the design of each package,
  written before its run; later changes are appended as dated notes.
- `WP6_VERIFICATION.md`, `WP7_VERIFICATION.md`, `WP8_VERIFICATION.md`,
  `WP8_REVERIFICATION.md` -- independent code reviews before launch, with
  verdicts and required changes.
- `WP5_INVENTORY.md` -- reconnaissance of candidate real data sets above
  d = 21, made before WP5 ran.
- `DELETED_DUPLICATE_TABLES.md` -- every NN table deleted by `70` and `77`:
  paths, sizes, MD5 sums, and the content comparison made at deletion.
- `INSTALLED_REGEN999_TABLES.md` -- the 0.1% NN tables installed by `75`,
  and what happened to the files they replaced.
- `FINDINGS.md` -- the running log of tr1 results, and the authority where a
  number is disputed.

## tr2/

PR-D-26-05767 revision-experiments pipeline. Scripts were renumbered to close
the gap left by moving tr1's `13`-`35` range into `tr1/`: the two duplicate
prefixes in the old flat layout (`01h_*` x2, `07_*` x2) were resolved, and the
old `36`-`42` scripts (added after tr1's numbering was already in use) were
renumbered down to `14`-`20`. Old name -> new name for every renamed file is
in the reorg's commit/session history; every internal `source()`,
`here::here()` and usage-comment path was updated to match.

Run order (02b precedes 02: it produces the raw high-dimensional CSVs 02
loads); purpose and main output condensed from each script's own header
comment:

| Script | Purpose | Main output |
|---|---|---|
| `00_env_check.R` | Load required packages; sanity-check the RK/NN quantile lookup tables; print `sessionInfo()`. | console only |
| `01_gen_quantile_table.R` | Production RK/NN quantile-table generator (the numerically-stable RK override auto-selects at d >= 342). | `.RData` tables in `R/RK-test_quantile/`, `R/NN-test_quantile/` |
| `01b_nn_component_probe.R` | Times the per-iteration cost components of the NN quantile generator at a given d, to size probes/extrapolate cost. | console only |
| `01c_validate_rk_stable.R` | Validates the log-space-stable RK weight computation against the original at dimensions where both work. | console only |
| `01d_nn_serial_iteration_probe.R` | Real-code-path timing check for NN at very high d (validates the component model from `01b`). | console only |
| `01e_nn_fast.R` | Optimized NN generator override (sourced by `01_gen_quantile_table.R`, not run standalone). | n/a (sourced) |
| `01f_validate_nn_fast.R` | Validates `01e_nn_fast.R`'s two optimizations against the original generator (correctness + timing + RAM). | console only |
| `01g_nn_sizes_quant.R` | Size-targeted NN null-quantile generation (knot sizes + interpolation) for full-data Musk/Speech. | `.RData` tables |
| `01h_nn_mc_yardstick.R` | Monte-Carlo-noise yardstick used to judge whether `01f`'s original-vs-fast differences are within seed noise. | console only |
| `01i_nn_d500_quant.R` | NN(UN-CCD) quantile table for the WP4 grid-2 cell d=500 (n=500). | `.RData` table |
| `02b_convert_highdim.py` | Converts raw ODDS/ADBench Musk/Speech/InternetAds (.npz) and Arrhythmia (.mat) into plain CSVs for `02_load_data.R`. | `data/outlier_detection/*_raw.csv` |
| `02_load_data.R` | Loads the 10 existing real-data benchmarks plus 4 new high-dimensional ones; exports all 14 to CSV. | `results/datasets_csv/*.csv`, `manifest.csv` |
| `03_smoke_test.R` | Smoke test + reproduction gate for the 9-method registry (WBC + synthetic Gaussian; Part C compares against the published manuscript table). | `results/tr2/smoke_wbc.csv`, `smoke_synthetic.csv` |
| `04_wp4_runtime.R` | WP4 wall-clock runtime/scalability grid (all 9 R methods, two grids over n and d). | `results/tr2/wp4_runtime2_{n,d}.csv` (+ `_raw`) |
| `04b_validate_parallel_radii.R` | Bit-exactness check for the parallel `nnccd.radi` override used by `04_wp4_runtime.R`. | console only |
| `05_wp4_runtime_pyod.py` | WP4 runtime grid for the PyOD baselines (ECOD/LUNAR/AutoEncoder fits, single-threaded). | `results/tr2/wp4_runtime2_pyod{,_raw}.csv` |
| `06_wp5_highdim.R` | WP5 high-dimensional real data: baselines on all four datasets (Arrhythmia, InternetAds, Musk, Speech); UN-CCD OOS/IOS on Arrhythmia, Musk, Speech only (InternetAds deliberately excluded). | `results/tr2/wp5_highdim_{metrics,raw}.csv` |
| `07_wp5_subsample_ccd.R` | WP5 follow-up: UNCCD-OOS/IOS on n=1000 subsamples of Musk and Speech. | `results/tr2/wp5_subsample_raw.csv` |
| `07b_wp5_fulldata_ccd.R` | Full-data Musk/Speech UNCCD via parallel radius search + spliced exact-n quantile tables. Supersedes the `07_*` n=1000 subsample workaround. | `results/tr2/wp5_fulldata_raw.csv`, cached `.rds` scores |
| `08_wp6_pyod_baselines.py` | WP6: fits ECOD, LUNAR, AutoEncoder on all 14 real datasets (5 seeds for the two stochastic methods). | `results/tr2/wp6_scores/*.csv`, `wp6_fit_log.csv` |
| `08b_wp6_metrics.R` | Computes TPR/TNR/BA/F2 from `08`'s raw PyOD scores, at two label thresholds. | `results/tr2/wp6_pyod_metrics.csv` |
| `09_wp3_synthetic.R` | WP3 cutoff-sensitivity study, synthetic part (3 settings x 4 CCD-based OS methods, 18-point cutoff-multiplier grid). | `results/tr2/wp3_synthetic_raw.csv`, `wp3_sensitivity_synthetic.csv` |
| `10_wp3_real.R` | WP3 cutoff-sensitivity study, real-data part (absolute cutoff grid; reproduces the WBC and Vowels manuscript rows as gates). | `results/tr2/wp3_sensitivity_real.csv` |
| `11_wp3_lines_plots.R` | Renders the WP3 BA/F2-vs-cutoff line plots (one PNG per synthetic setting and per real dataset). | `results/tr2/figures/wp3_{synthetic,real}_*.png` |
| `11b_wp3_synthetic_f2_plot.R` | Manuscript Appendix C figure: three-panel F2-only curve for the synthetic settings. | `results/tr2/figures/wp3_synthetic_F2.png` |
| `11c_wp3_real_f2_plot.R` | Manuscript Appendix C figure: three-panel F2-only curve for the real datasets (WBC, Thyroid, vowels). | `results/tr2/figures/wp3_real_F2.png` |
| `12_highd_madn_recheck.R` | Re-measures high-d synthetic cells under both the buggy and fixed `std_MADN`, to check whether the "IOS > OOS at high d" claim was a standardization artifact. | `results/tr2/highd_madn_recheck*.csv` |
| `13_glass_os_regen.R` | Regenerates the four OS-method rows for glass only, using `evaluate()` (fixes the positional-scoring bug on glass's mis-sorted outliers). | `results/tr2/glass_os_regen.csv` |
| `14_patch_manuscript_tables.ps1` | Replaces TR2's four MCCD rows in its two real-data tables with values regenerated for TR1 (one-time, author-authorized; see the script's own header for scope/authorization notes). | edits outside this repo |
| `15_os_repro_audit.R` | Audits whether the four OS rows of the manuscript's real-data tables reproduce under the current harness, dataset by dataset. | `results/tr2/os_repro_audit.csv` |
| `16_speech_ios_clusters.R` | Recovers per-cluster labels to explain why standardized IOS collapses to zero on 57% of Speech. | `results/tr2/speech_ios_clusters.csv` |
| `17_tiebreak_impact.R` | Bounds how much the tie-breaking rule (legacy vs. manuscript Eq. 10, whole-dataset vs. per-cluster) moves the published real-data OS numbers. | `results/tr2/tiebreak_impact*.csv` |
| `18_wp2_assumption_checks.R` | WP2 checklist: empty-outbound-neighbourhood handling, and similarity-equivariance of the radius search under a random rigid+scale transform. | `results/tr2/wp2_assumption_checks.csv` |
| `19_musk_geometry_audit.R` | Provenance audit (Ceyhan checklist item 8), Musk half: confirms Appendix D numbers come from the same construction as the reported AUC. | `results/tr2/musk_geometry_audit.csv` |
| `20_wp_collinearity.R` | Ad hoc collinearity-sensitivity sweep (equicorrelated Gaussian clusters, rho in {0, .3, .5, .7, .9}) backing the manuscript's collinearity-stability claim. | `results/tr2/wp_collinearity_{agg,raw}.csv` |
| `test_outlyingness_density.R` | Regression test for the vicinity-density overflow fix (`(size/R^d)^(1/d)` -> `size^(1/d)/R`). | console only |
| `test_std_madn.R` | Regression test for the `std_MADN` fallback (MADN=0 -> SD -> 0) behavior. | console only |
| `verify_endpoint_ties.R` | Verifies the tie-breaking implementation retains ties that reach either endpoint of the ranked sample. | console only |

`FINDINGS.md` is the running log and the authority where a number is
disputed. `USER_RUN_TABLES.md` is the command sheet for the NN quantile-table
production runs.
