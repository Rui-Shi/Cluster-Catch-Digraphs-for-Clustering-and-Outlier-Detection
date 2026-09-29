# tr1 findings index

Shi, R., Ceyhan, E., Billor, N. *Shape-Adaptive Outlier Detection Using
Cluster and Mutual Catch Digraphs.* Neurocomputing (in revision, 2026).

This file maps every table and headline number of the paper and its
supplementary material that comes from a tr1 experiment to the script that
produced it and the committed file that holds it. Script names are in
`revision_experiments/tr1/`; result paths are relative to
`revision_experiments/results/tr1/`. `revision_experiments/README.md` lists
every script with its purpose and outputs.

## 0. How to read this file

- Authority: a committed result file beats a findings note (`*_findings.md`,
  `WP*_FINDINGS*.md`), which beats this index.
- Most drivers checkpoint per cell, so rerunning one against a committed
  file computes nothing new; move the file aside first.
- The Monte Carlo quantile tables are not in the repository. A missing table
  is generated on the spot by `R/ccds/quantile_table.R`, so a rerun
  reproduces a published number up to Monte Carlo error, not bit for bit.
- Tables and results not listed here (the legacy simulation grid of the
  paper's Section 5, apart from the SUN-MCCD Neyman-Scott d = 10 cells) come
  from the original drivers under `simulations/outlier_detection/`, whose
  SLURM logs (`slurm-*.out`) hold the 1000-replicate means.

## 1. Main text

| Item | Reports | Script(s) | Committed file(s) | Notes |
|---|---|---|---|---|
| `tab:alpha_sim` | MC-SRT level by method family and simulation dimension | `shared/harness.R` resolvers; `68` (grid check) | `wp2a_verify_sim_grid.csv` | SUN-MCCD Neyman-Scott d = 10 drivers now load the 0.1% table (see `97`). |
| `tab:alpha_real` | MC-SRT level per real data set | `63`-`68`, `71`-`78` | `wp2a_provenance_verify.csv`, `wp2a_verify_rk_distinct.csv`, `wp2c_manuscript_tables.csv` | `tr1/INSTALLED_REGEN999_TABLES.md` records the regenerated 0.1% NN tables. |
| `tab:Real_Data` | n, d, contamination of the 16 real data sets | `32`, `80` | `dataset_inventory.csv`, `wp4/data/manifest.csv` | Inclusion rule: `tr1/BENCHMARK_EXPANSION_RULE.md`. |
| `tab:Real_Data_Result_Summary` | F2 ranks over 17 methods | `34`, `78`, `82` | `final_comparison.csv`, `wp2c_manuscript_tables.csv`, `wp4/wp4_metrics_main.csv` | Mutual-kNN and SNN at k = round(sqrt(n)). |
| `tab:Real_Data_Aggregate` | Mean/median F2 and BA, two method groups | `82`, `78` | `wp4/wp4_metrics_main.csv`, `wp2c_manuscript_tables.csv` | Label-chosen k rows: `MutualKNN-oracle`, `SNN-oracle`. |
| `tab:wp5_highd` | Results above d = 21 (letter, mnist, musk, arrhythmia) | `84`, `84b` | `wp5/wp5_metrics_main.csv`, `wp5/wp5_highd_results.csv` | U-/SU-MCCD n/a at musk and arrhythmia. |
| `tab:wp6_runtime` | Median time and memory | `91`, `91b`, `95` | `wp6/91_wp6_runtime_n.csv`, `wp6/91_wp6_runtime_d.csv`, `wp6/91b_wp6_runtime_py_summary.csv`, `wp6_incr/95_ab_runtime_n.csv`, `wp6_incr/95_ab_runtime_d.csv` | UN-/SUN-MCCD rows from `94`'s paired run through `95`. |
| `tab:Complex_Cluster_Compact` | Neyman-Scott F2 ranks | legacy drivers; `97` | `ns_sun_d10/97_{Matern,Thomas,Mix}_d10_999_sun.csv` | SUN-MCCD at d = 10 only; other cells from the legacy logs. |
| Wilcoxon, 5 of 32 | SUN-MCCD against 16 methods on F2 and BA, Holm per block | `92` | `wp3/wilcoxon_real.csv` | Significant: UN-MCCD on F2; DBSCAN, DIF, ODIN, MST on BA. |
| Mutual-kNN, 11 of 16 | SUN-MCCD F2 above mutual-kNN (both k rules) | `82` | `wp4/wp4_metrics_main.csv`, `wp4/WP4_FINDINGS.md` | UN-MCCD above it on 7 of 16. |
| Runtime exponents | log-log slopes of time and memory in n | `91`, `95` | `wp6/91_wp6_slope.csv` (U-/SU-MCCD 1.927/1.899), `wp6_incr/95_ab_slope.csv` (UN-/SUN-MCCD 1.988/1.932) | |
| Search direction | n = 500, d = 10: descending 6.55/7.43 s, ascending 0.94/2.07 s (UN-/SUN-MCCD) | `96`, `95` | `wp6_incr/96_direction_raw.csv`, `wp6_incr/95_direction_summary.csv` | Medians of 10 reps. |
| Neyman-Scott d = 10 SUN-MCCD | F2 0.942 (Matern), 0.891 (Thomas), 0.903 (mixed) at 0.1% | `97` | `ns_sun_d10/97_*_d10_999_sun.csv` | 1000 reps each; the 1% pass (`*_d10_99_sun.csv`) reproduces the drivers' logs. |
| RK zero-quantile share | 51.3%, 91.4%, 100%, 100% at d = 32, 100, 166, 274 | `84` | `wp5/rk_degeneracy.csv` | |

The four remaining main-text tables (`tab:delta-lineage`, `tab:acronym`,
`tab:mccd-design-family`, `tab:space_time`) hold definitions and proved
bounds, not results.

## 2. Supplementary material

| Label | Reports | Script(s) | Committed file(s) |
|---|---|---|---|
| `SM-tab:Real_Data_Result2.1`, `2.2` | TPR/TNR and BA/F2 of the nine original methods on 16 sets | `28`, `33`, `34`, `78` | `regen_final_*.csv`, `final_baselines.csv`, `regen_baselines.csv`, `final_comparison.csv`, `wp2c_manuscript_tables.csv`, `lof_original8.csv` |
| `SM-tab:WP4_Aggregate17`, `WP4_T2`, `WP4_Result1`, `WP4_Result2`, `WP4_mkNN`, `WP4_native`, `WP4_seeds` | The 17-method comparison, at contamination 0.1 and at the true rate; native labels; DIF/LUNAR seed spread | `81`, `82`, `83` | `wp4/wp4_metrics_main.csv`, `wp4/wp4_metrics_long.csv`, `wp4/scores/`, `wp4/WP4_FINDINGS.md` |
| `SM-tab:stability` | Label changes and BA spread under perturbation | `45` | `wp2b_stability_summary_smin005.csv`, `wp2b_determinism_smin005.csv`, `wp2b_seed_stability_smin005.csv` |
| `SM-tab:alpha_rk_real` | RK-based detectors over the alpha range | `42`, `46`, `47` | `wp2a_alpha_sweep.csv`, `wp2a_alpha_summary.csv`, `wp2a_alpha_operating_range.csv` |
| `SM-tab:alpha_nn_real` | NND-based detectors, wilt and pageblocks | `50` | `wp2a_nn_alpha_bigsets.csv` |
| `SM-tab:rk_degeneracy` | Share of zero RK critical values, d = 5-50 | `49` | `wp2a_rk_table_degeneracy.csv` |
| `SM-tab:rk_behaviour` | RK-based detectors on a synthetic dimension grid | `49` | `wp2a_rk_saturation_summary.csv`, `wp2a_rk_saturation.csv` |
| `SM-tab:smin_plateau` | S_min sensitivity, five largest sets | `51`, `54`, `78b` | `wp2a_smin_bigfour.csv`, `wp2c_plateau_table.csv`, `wp2c_sun_smin_rerun_fixed.csv` |
| `SM-tab:smin_rules` | Six S_min rules on simulated data | `55`, `59` | `wp2c_simulation_summary.csv`, `wp2c_simulation_arm.csv` |
| `SM-tab:smin_zero_real` | Real-data means at four S_min values | `60`, `78b` | `wp2c_smin_zero_real.csv`, `wp2c_sun_smin_rerun_fixed.csv` |
| `SM-tab:param_provenance` | Parameter source per method | `92` | `wp3/parameter_provenance.csv` |
| `SM-tab:wp9-ablation`, `wp9-paired` | SUN-MCCD component ablation | `89` | `wp9/wp9_ablation_summary.csv`, `wp9/wp9_paired_diff.csv`, `wp9/wp9_ablation.csv` |
| `SM-tab:wilcoxon`, `wilcoxon_U`, `wilcoxon_SU`, `wilcoxon_UN` | Paired Wilcoxon tests, Holm per block | `92` | `wp3/wilcoxon_real.csv` |
| `SM-tab:wp5_datasets` | Data sets above d = 21 | `84a` | `wp5/data/manifest.csv` |
| `SM-tab:wp5_rk` | Zero RK critical values above d = 21 | `84` | `wp5/rk_degeneracy.csv` |
| `SM-tab:wp5_full` | All methods above d = 21 | `84`, `84b` | `wp5/wp5_metrics_long.csv`, `wp5/wp5_metrics_main.csv` |
| `SM-tab:Matern_*`, `Thomas_*`, `Matern-Thomas_*`, `Complex_Cluster_Ranking` | Neyman-Scott results and ranks | legacy drivers; `97` | `ns_sun_d10/` (SUN-MCCD, d = 10) |
| `SM-tab:wp8-boundary-gradient`, `wp8-boundary-overall` | Boundary false positives | `85` | `wp8/85_boundary_fp_summary.csv`, `wp8/85_boundary_fp.csv` |
| `SM-tab:wp8-outliertypes` | Outlier types | `86` | `wp8/86_outlier_types_summary.csv` |
| `SM-tab:wp8-smallcluster-flag`, `wp8-smallcluster-p3` | A small third cluster | `87` | `wp8/87_small_cluster_summary.csv` |
| `SM-tab:wp8-csr-alpha`, `wp8-csr-flagmatrix` | Departures from spatial randomness | `88` | `wp8/88_csr_violation_summary.csv` |
| `SM-tab:wp7-ari`, `wp7-khat` | Clustering quality | `90` | `wp7/wp7_summary_metrics.csv`, `wp7/wp7_khat_table.csv`, `wp7/wp7_metrics_long.csv` |
| `SM-tab:wp6-n-r`, `wp6-d-r` | Runtime of the R methods | `91`, `95` | `wp6/91_wp6_runtime_n.csv`, `wp6/91_wp6_runtime_d.csv`, `wp6_incr/95_ab_runtime_n.csv`, `wp6_incr/95_ab_runtime_d.csv` |
| `SM-tab:wp6-n-py`, `wp6-d-py` | Runtime of the Python competitors | `91b` | `wp6/91b_wp6_runtime_py_summary.csv` |
| `SM-tab:wp6-slope` | Log-log exponents | `91`, `95` | `wp6/91_wp6_slope.csv`, `wp6_incr/95_ab_slope.csv` |

`SM-tab:focused-mc-settings` summarizes the legacy simulation settings.

## 3. Supersessions and overrides

- `wp2c_sun_smin_rerun_fixed.csv` (`78b`, repaired 0.1% NN tables)
  supersedes the SUN-MCCD cells of `wp2c_smin_zero_real.csv` and the SUN-MCCD
  waveform and PenDigits cells of `wp2a_smin_sweep16.csv` and
  `wp2a_smin_bigfour.csv` at the values it covers.
- The published stability results are the S_min = 0.05 vintage
  (`wp2b_*_smin005.csv`); the files without the suffix are the earlier
  S_min = 0.0625 run.
- `wp2c_rerun_alpha001.csv` (six cells rerun with the regenerated 0.1% NN
  tables) and `wp2c_manuscript_tables.csv` (`78`) supersede the matching
  cells of `final_comparison.csv`.
- `wp2a_alpha_claim_audit.csv` is the record before the table repair;
  `wp2a_provenance_verify.csv` is the check after it.
- `wp6_incr/94_ab_raw.csv` (paired timing) supersedes the `wp6_incr/91_*`
  re-timing of 2026-09-23, which ran on another day; `95` summarizes `94`.
- The tr1 table-provenance scripts (`62`-`78`) record which duplicated NN
  tables were deleted and which 0.1% tables were installed in
  `tr1/DELETED_DUPLICATE_TABLES.md` and `tr1/INSTALLED_REGEN999_TABLES.md`.

## 4. Not in the repository, and how to get it

| What | How |
|---|---|
| Quantile tables (`R/RK-test_quantile/*.RData`, `R/NN-test_quantile/*.RData`) | Generated on the spot by `R/ccds/quantile_table.R` (settings and costs in `revision_experiments/README.md`). The original RK tables also hold the raw draws `Kest.m` that `42`, `46` and `49` re-reduce; those three need the original files. |
| Uncompressed WP8 cluster files | `gunzip -k wp8/*_clusters.csv.gz` (needed by `90`). |
| WP4 data CSVs | `80` rebuilds them from `data/outlier_detection/`. |
| Console logs, checkpoint files, smoke runs | Rerun the script. |
| Per-replicate output of the legacy Section 5 grid | Never kept; the SLURM logs hold the 1000-replicate means only. `wp3/sim_variability_INVENTORY.md` records the search. |

## 5. Log

- **2026-08-09.** The nine original methods were regenerated on the 16 real
  data sets through the harness (`13`-`39`). Two loaders sort by the wrong
  column (glass, ecoli; `26`), which breaks positional scoring; all scoring
  now uses `evaluate()`.
- **2026-08-10.** Several NN tables were byte-identical at two levels (`62`,
  `70`); their true levels were estimated (`63`-`69`). S_min changed from
  0.0625 to the label-free constant 0.05 (`54`, `55`, `60`, `61`).
- **2026-08-11 to 12.** Genuine 0.1% NN tables were generated for d = 12,
  16, 18, 19 and 21 (`71`-`75`) and the six affected cells rerun
  (`wp2c_rerun_alpha001.csv`).
- **2026-09-05.** Harness audit: the NND test is a union of two marginal
  level-alpha tests (no Holm step); `get_simul()` requires an explicit
  level; tables shorter than the data are refused (`79`). WP3-WP9 ran under
  their protocols; the stability runs were repeated at S_min = 0.05.
- **2026-09-22 to 23.** SUN-MCCD S_min cells rerun on the repaired tables
  (`78b`). The NND radius search became incremental; radii and detector
  outputs are identical to the per-step search (`92b`, `93`) and the search
  is faster (`94`).
- **2026-09-24.** Ascending against descending radius search timed (`96`).
- **2026-09-26.** SUN-MCCD in the three Neyman-Scott settings at d = 10
  rerun at the stated 0.1% level (`97`; the drivers had loaded the 1%
  table). Mutual-kNN and SNN take k = round(sqrt(n)) as the primary rule.
  Two groups of SU-MCCD simulation drivers were fixed and checked with `98`
  (`su_driver_fix/`):
  - the d = 3, n = 1000 drivers had set d = 2. Rerun at d = 3 (22
    replicates): uniform TPR 1.000, TNR 1.000; Gaussian TPR 1.000,
    TNR 0.903;
  - the d = 50, n = 50 drivers had loaded the d = 2, 1% RK table. With that
    table the rerun reproduces the logged values (uniform TNR 0.216). With
    the d = 50, 0.1% table: uniform TPR 1.000, TNR 0.591, F2 0.379; Gaussian
    TPR 1.000, TNR 0.588, F2 0.370 (1000 replicates each).
- **2026-09-27.** tr1 scripts consolidated in `tr1/`; results committed.
- **2026-09-28.** Missing quantile tables are generated on the spot
  (`R/ccds/quantile_table.R`).

## 6. Protocols

Each package's design was committed before, or with, its driver:

| Protocol | Committed | Driver | Driver committed |
|---|---|---|---|
| `BENCHMARK_EXPANSION_RULE.md` | 294f8a5, 2026-08-09 00:29 | `33_final_baselines.R` (baselines on the eight added data sets) | 03fc090, 00:30 |
| `WP4_PROTOCOL.md` | 7653a17, 2026-09-05 17:18 | `81_wp4_baselines.py` | 8380363, 17:27 |
| `WP9_PROTOCOL.md` | c3afa38, 2026-09-05 18:45 | `89_wp9_ablation.R` | c3afa38 (same commit) |
| `WP3_PROTOCOL.md` | 29d1bdf, 2026-09-05 18:55 | `92_wp3_stats.R` | c9c23ea, 18:59 |
| `WP5_PROTOCOL.md` | d0594ab, 2026-09-05 18:55 | `84_wp5_highd.R` | c64985d, 19:00 |
| `WP8_PROTOCOL.md` | 5ed690c, 2026-09-05 18:59 | `85_wp8_boundary_fp.R` | 5ed690c (same commit) |
| `WP6_PROTOCOL.md` | caecfe8, 2026-09-05 23:06 | `91_wp6_runtime.R` | bd0681d, 23:12 |
| `WP7_PROTOCOL.md` | a0ee266, 2026-09-05 23:12 | `90_wp7_clustering_quality.R` | 7b1d350, 23:18 |

Changes made after a run are appended to each protocol as dated notes.
