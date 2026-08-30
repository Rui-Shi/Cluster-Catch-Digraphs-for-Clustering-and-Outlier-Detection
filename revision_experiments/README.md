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

## Running

Scripts are run from the repo clone root (paths use `here::here()`). R
package prerequisites are checked by `tr2/00_env_check.R`. On this project's
Windows setup, `Rscript` must be launched from PowerShell, not Git-Bash/MSYS
(known segfault). The Python steps (`tr2/02b_convert_highdim.py`,
`tr2/05_wp4_runtime_pyod.py`, `tr2/08_wp6_pyod_baselines.py`) use the bundled
`revision_experiments/.venv` (interpreter at `.venv/python.exe`).

## shared/

- **`harness.R`** -- the 9(+)-method registry, `evaluate()`,
  `load_real_dataset()`, `get_simul()`, checkpointing
  (`has_result()`/`append_result()`). Read-only by convention: neither
  project edits it directly; tr1's `wp0_mccd_methods.R` only *appends* to its
  registries. It only reads `R/`, `methods/`, `simulations/` and `data/` from
  the repo root, so it works unchanged from either `tr1/` or `tr2/` callers.

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

## .venv/

The pinned PyOD 3.6.1 / torch 2.13.0+cpu Python environment, used by tr2's
`05_wp4_runtime_pyod.py` and `08_wp6_pyod_baselines.py`. Nothing tr1-specific
uses Python. Unmoved by this reorg (still `revision_experiments/.venv`).

## results/

- `results/datasets_csv/` -- 16 exported real-data CSVs + `manifest.csv`
  (two of the 16, `Musk_sub1000.csv` and `Speech_sub1000.csv`, come from the
  WP5 subsample path, `07_wp5_subsample_ccd.R`, not the main loader). Written
  by tr2; read by tr2 scripts that need a dataset as a flat CSV rather than
  through `load_real_dataset()`. tr1 does not read this directory -- it loads
  data only via `load_real_dataset()`, per `tr1/REGEN_SPEC.md`. Not moved or
  touched by this reorg.
- `results/scores_cache/` -- cached `.rds` score vectors keyed by
  `<dataset>_<method>.rds`. Shared in principle, written by tr2's
  `06_wp5_highdim.R` / `07_wp5_subsample_ccd.R` / `07b_wp5_fulldata_ccd.R`
  and read by `10_wp3_real.R`. Not moved or touched by this reorg.
- `results/tr1/` -- tr1's flat result files (csv/log). Contents untouched by
  this reorg.
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

NEUCOM-D-26-15191 regeneration pipeline. Script names and numbering
(`13_wp0_gate.R` .. `39_verify_manuscript_tables.R`, plus `wp0_mccd_methods.R`
and the `regen_wilt_*_launcher.ps1` pair) are unchanged from before this reorg
-- only their location moved, from `revision_experiments/` to
`revision_experiments/tr1/`.
See `tr1/REGEN_SPEC.md`, `tr1/BENCHMARK_EXPANSION_RULE.md` and
`tr1/REGENERATION_REPORT.md` for tr1's own documentation; this file does not
duplicate it.

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
