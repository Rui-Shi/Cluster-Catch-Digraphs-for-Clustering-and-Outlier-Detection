# tr3: UN-CCD clustering paper

Revision experiments for the UN-CCD clustering paper (Shi, Ceyhan and Billor), now titled "Cluster Catch
Digraphs: Nearest-Neighbor Covering Radii and Label-Free Parameter Selection". It was submitted to KBS as
"Clustering with Uniformity- and Neighbor-Based Random Geometric Graphs" and rejected (KNOSYS-D-26-08010).
The manuscript is at `G:\Submissions\TR3\TR3_10182026\`. The map from each paper
result to its original script is `CLUSTERING_REPRODUCTION_GUIDE.md` at the clone root.

These drivers replace the 120 copy-paste simulation scripts and the three
real-data scripts for new work. The originals print averages to SLURM logs and
save the workspace. The drivers run the same code through one interface, with
the same settings and data, and write every replicate to a CSV. Nothing in
`methods/`, `R/` or `data/` is edited, except the descending radius search in `R/ccds/UN_CCD.R` (index fix
2026-10-23, sorted neighbour lists 2026-11-02; see their sections). In `simulations/`, the following settings were
changed on 2026-10-18 and 2026-10-23 by the authors' decision (see the audit section):
- the table line of the four UN-CCD Gaussian d = 10 scripts, from 0.1% to 1%;
- `iteN` of the two RK-CCD d = 3, n = 500 scripts, from 50 to 500;
- the three UN-CCD real-data table lines, changed 2026-10-23 (asymmetric `2d_80%` to `2d_85%`,
  user_knowledge `6d_95%` to `5d_95%`, wholesale `8d_95%` to `6d_95%`);
- the two RK-CCD real-data table lines, changed 2026-10-23 (user_knowledge `6d_99%` to `5d_99%`,
  wholesale `8d_99%` to `6d_99%`; `quant` stays 0.99);
- the two UN-CCD real-data search directions, changed 2026-10-23 (seeds and breast_cancer, `method` from
  `"descend"` to `"ascend"`; tables `7d_95%` and `9d_95%` stay);
- the five baseline real-data scripts (DBSCAN, MST, K-means++, Spectral, Louvain), changed 2026-10-23: a
  `set.seed(1)` per block, K-means++ on wholesale `k = 3`, Spectral on the data with a Gaussian kernel,
  DBSCAN with an elbow rule for Eps (see "Baseline fixes (2026-10-23)" below).

## Status

On `main` since 2026-11-03: the drivers and scripts behind the published numbers, their inputs and
results, the shared-code check `38_descend_pointer_check.R`, the settings check `00_check_settings.R`
and every protocol. The checks and audits that no published number rests on were removed from `main`
on 2026-11-03 and stay on branch `tr3-revision`: the reproduction gate `20_gate.R` (87 of 87 passed,
last run 2026-11-03) and its reference script `12_real_reference.R`, `21_iteN500_check.R`,
`25_greedy_prefix.R`, and their outputs (`gate.csv`, `real_reference.csv`, `iteN500_check.csv`,
`greedy_prefix/`), the audits `audit_cells.csv`, `audit_logs.csv`, `audit_repro.csv`,
`real_verify_ccd.csv`, `real_verify_baselines.csv`, and the superseded `real_baselines_eps.csv`.
The sections below still describe them as they were when written. Each new experiment got a
`*_PROTOCOL.md` here before it ran.

## Files

| File | What it does |
|---|---|
| `settings_simulation.csv` | The settings of the 120 simulation drivers, as the code ran them for the published numbers (not as the paper states them): seed, iteN, cls_dis, sig, radius range, size of the first of three clusters, quantile table, RK level, UN direction, KS grid. |
| `settings_realdata.csv` | Table, level, direction and grid per method and real data set (27 rows). |
| `00_check_settings.R` | Re-reads all 120 drivers and the three real-data scripts and stops if any value differs from the two tables, or if the drivers differ in anything else. |
| `tr3_lib.R` | Shared code: `ccd_cluster()`, data generation in the scripts' random-number order, the patched loader, CSV checkpointing. |
| `10_sim.R` | One simulation cell group (method, distribution, d, n). |
| `11_real.R` | CCD rows of the real-data table. |
| `ITEN500_PROTOCOL.md` | Variant `iteN500`: RK-CCDs at d = 3, n = 500 with 500 repetitions (its check script, `21_iteN500_check.R`, is on branch `tr3-revision`). |
| `13_real_baselines.R` | Baseline rows of the real-data table (DBSCAN, MST, K-means++, Spectral, Louvain): runs each data-set block of the original script verbatim, once or over seeds. |
| `23_noise.R` | The noise experiment (d = 3, 10; l = 3% to 20%), written 2026-10-24 because the original driver is lost. Its rows are in the paper; it does not reproduce the old table, see "Noise experiment (2026-10-24)". Variant `allpoints` is experiment 8 (k chosen on all points, 2026-10-31). |
| `24_alpha_supplement.R` | Supplement Section S2.1 (2026-10-25): writes `results/tr3/supp_alpha.tex`, the rows of the α-sweep table (config A, seeds from config S, the paper's current settings) and of the shipped-against-regenerated table (configs A and G, seeds left out), from `real_alpha.csv`; checks that the rule's cells equal `real.csv`. Nothing is rerun. |
| `22_runtime.R`, `22_runtime_run.ps1`, `22_runtime_monitor.ps1`, `22_runtime_summary.R` | Re-timing of `tab:runtime_comparison` (d = 5, uniform, n = 50 to 500, three methods) on this machine, with a load guard; see "Runtime comparison (2026-10-24)". |
| `ALPHA_PROTOCOL.md`, `14_alpha_tables.R`, `15_real_alpha.R`, `16_alpha_summary.R` | One α for UN-CCDs on all real data sets (20% to 1%): generates the missing NN tables, runs every cell, checks the reproduction and scores the levels. |
| `NEW_EXPERIMENTS_PROTOCOL.md` | Design, decisions and outcomes of the new experiments for the KBS rebuttal (2026-10-25): 1 real data with d ≥ 13, 2 Type I error, 3 search direction, 5 paired tests. See "New experiments (2026-10-25)". |
| `30_hd_load.R`, `settings_realdata_hd.csv` | Experiment 1: `tr3_load_hd()` for wine, hepatitis, lymphography, wdbc, vehicle and segmentation (raw files in `revision_experiments/tr3_data/` and `data/clustering/`), and the CCD settings by the rule for d ≥ 10 (`30_hd_load.R settings` rewrites the CSV). |
| `31_hd_real.R` | Experiment 1: UN-, RK- or KS-CCDs on those six sets (the counterpart of `11_real.R`). |
| `32_hd_baselines.R` | Experiment 1: DBSCAN, MST, K-means++, Spectral or Louvain (or `all`) on those six sets: each original script's `seeds` block with the data object (and k) substituted. `check` reproduces `real_baselines_fixed.csv`; `show` prints the evaluated block. |
| `33_type1.R` | Experiment 2: size of one MC-SRT under CSR (two nulls) and UN-CCDs on uniform data without clusters. |
| `34_ablation_direction.R` | Experiment 3: ascending against descending radius search at d = 10 and 20, n = 200, 500 replicates. |
| `35_paired_tests.R` | Experiment 5: paired Wilcoxon tests of UN- against RK- and KS-CCDs at n = 50 and 200, Holm within each n (reads existing results; see "Paired tests at n = 200 (2026-11-01/02)"). |
| `36_new_tables.R` | The LaTeX tabular bodies of experiments 1, 2, 3, 5, 7 and 8 (`results/tr3/tex/new_*.tex`) and the CSV summaries of experiments 7 and 8, with checks; nothing is rerun. |
| `37_sim_baselines.R` | Experiment 7: DBSCAN, MST, K-means++, Spectral or Louvain (or `all`) on the simulated data (`sim`) and the noise data (`noise`) of the CCDs; `check` reproduces the stored real-data baseline rows first. |
| `38_descend_pointer_check.R` | Checks that the descending radius search of 2026-11-02 (sorted neighbour lists) returns exactly the radii of the 2026-10-23 version, on simulated and real data with d ≥ 10 and edge cases; see "Descending search with sorted neighbour lists (2026-11-02)". |

## Running

From the clone root, in PowerShell (Rscript segfaults under Git Bash):

```powershell
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
& $R revision_experiments/tr3/00_check_settings.R
& $R revision_experiments/tr3/10_sim.R UN-CCDs uniform 2 50          # [clusters] [variant]
& $R revision_experiments/tr3/11_real.R RK-CCDs iris,seeds,R15
```

`TR3_CORES` sets the number of workers (default: all). Results do not depend on it.
A rerun skips every finished cell; move a result file aside to recompute.

## How the drivers stay faithful to the scripts

- **One method per process.** Each process loads only its method, in the script's
  order: wrapper, `Uni-Gau_cls.R`, then `quantile_table.R`. The real-data runs
  load the loader's packages too. `Uni-Gau_cls.R` attaches spatstat and
  hypercube/plotly, which mask igraph, cluster and MASS functions.
- **Same parallel call.** Workers run
  `foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar%`, called
  at top level. foreach exports only what the calling environment holds, so
  the call cannot sit inside a function.
- **Data.** Each script group calls `set.seed(2024)` once and then generates the
  2-, 3- and 5-cluster lists, in that order. Reruns regenerate all three even
  when they skip finished cells.
- **No table generation.** Quantile-table generation is off
  (`CCD_QUANTILE_GENERATE=FALSE`), so a missing table stops the run. RK tables
  are sent to the workers without their raw draws (`Kest.m`), which the RK
  code never reads.
- **The one code fix.** `UNMCCD_outlier()` (undefined) becomes
  `UNCCD_clustering()` with the same arguments.
- **The loader fix.** The real-data loader's `setwd()` is fixed in memory, and
  `readxl` is needed only for `user_knowledge`.

## Results (`revision_experiments/results/tr3/`)

`results/` is gitignored. As for tr1, the result files were force-added on 2026-11-03
(90 files; 78 on `main` since the removals listed under "Status", `git ls-files revision_experiments/results/tr3`): the files of this table, the
`*_done.csv` files (the readers keep only the cells listed there), `_before_alpha_adopt_20261023/`
(read by `16_alpha_summary.R`) and the files the sections below cite as results. Not committed: the
35 quantile tables in `alpha_tables/` (`*.RData`, written by `14_alpha_tables.R`). Removed from the
clone on 2026-11-03 (408 files, kept outside the repository): the console logs (`_logs/`,
`descend_ptr_20261102/logs/`), the backup folders `_before_*` other than `_before_alpha_adopt_20261023/`,
the part files (`real_alpha_parts/`, `real_baselines_fixed_<method>_run<n>[_params].csv`,
`real_baselines_fixed_v2_*`, `real_baselines_fixed_mst2_after_revert_*`), the MST Sil-max runs
(`real_baselines_fixed_mstrule_*`) and the driver copies in `descend_ptr_20261102/driver_copies/`.
The sections below this table still name some of these files as they were when written.

| File | One row per |
|---|---|
| `sim_replicates.csv` | replicate: variant, method, dist, d, n, n_clusters, iteN, rep, ari, sil, k_hat, seconds, param (KS: chosen m), data_fp |
| `sim_done.csv` | finished cell (what reruns skip) |
| `sim_ks_search.csv` | KS grid value per cell: delta, mean_sil, mean_ari, sr, seconds |
| `real.csv`, `real_done.csv` | method and data set: n, d, ari, sil, k_hat, seconds, param, table |
| `real_ks_search.csv` | KS grid value per data set: si |
| `real_baselines.csv`, `real_baselines_seeds.csv` | output of `13_real_baselines.R`: verbatim runs, and seeds 1 to 30 (both from before the baseline fixes of 2026-10-23) |
| `real_baselines_fixed.csv`, `real_baselines_fixed_run2.csv`, `real_baselines_fixed_params.csv` | the baselines after the fixes of 2026-10-23: first run, second run (identical), and the Eps / sigma each block chose (the per-method part files `real_baselines_fixed_<method>_run<1,2>[_params].csv` were removed on 2026-11-03) |
| `real_alpha.csv`, `alpha_summary.csv` | α sweep: config, level, data set (ARI, Sil, k_hat, table, source); the criteria per config and level. Generated tables in `alpha_tables/` |
| `real_hd.csv`, `real_hd_done.csv`, `real_hd_ks_search.csv` | experiment 1, CCDs: variant (`hd`), method, data set, n, d, ari, sil, k_hat, seconds, param (KS: chosen m), table; KS grid value per data set: si |
| `real_hd_baselines.csv`, `real_hd_baselines_params.csv`, `real_hd_baselines_done.csv` | experiment 1, baselines: method, data set, n, d, ari, sil, k_hat, seed, script; the Eps / sigma / bandwidth each block chose |
| `type1_test.csv`, `type1_test_summary.csv`, `type1_ccd_replicates.csv`, `type1_ccd_summary.csv` | experiment 2: rejection counts per (d, m) and null; UN-CCDs per replicate on uniform data; summaries |
| `ablation_direction_replicates.csv`, `ablation_direction_summary.csv`, `ablation_direction_repro.csv` | experiment 3: ascending replicates, the comparison per cell, the reproduction check of the descending rows |
| `paired_tests.csv`, `paired_tests_check.csv` | experiment 5: one row per n (50, 200), cell and comparison; cell means against the supplement |
| `simbase_replicates.csv`, `simbase_done.csv`; `noisebase_replicates.csv`, `noisebase_done.csv` | experiment 7: one row per replicate (method, setting or noise cell, mode, rep, ari, sil, k_hat, k_hat_block, n_label0, param, seconds, error, data_fp, fp_check); finished cells |
| `simbase_summary.csv`, `simbase_by_d.csv`, `simbase_comparisons.csv`, `noisebase_summary.csv` | experiments 7 and 8, from `36_new_tables.R`: per setting (eight methods; CCD rows read from the supplement), per distribution and d, the comparisons by the 0.02 rule, and the noise cells (since experiment 8, "UN-CCDs" there means k chosen on all points; "UN-CCDs (cluster points)" is the published choice) |
| `noise_allpoints_replicates.csv`, `noise_allpoints_done.csv` | experiment 8: per replicate the published columns of `noise_replicates.csv` (gate: equal) plus ari_all, sil200_all and nlab_all for k chosen on all points |
| `tex/new_*.tex` | output of `36_new_tables.R` |
| `descend_pointer_check.csv` | output of `38_descend_pointer_check.R`: case, n, d, low.num, scores, pass, seconds of both versions |
| `descend_ptr_20261102/` | the driver runs that checked the 2026-11-02 descending search against the published files (see its section): their outputs (`real_hd_UN.csv`, `sim_replicates_UN_d10_d20.csv` and the done files) and `compare_summary.json`, written by `compare_verify.py` |

Numbers are written with 17 significant digits, so means recomputed from the
CSV equal the in-memory ones exactly. Each cell is appended as soon as it
finishes, because G: drops off under sustained writes. Readers keep only the
cells listed in the matching `*_done.csv`.

Definitions follow the scripts:

- **Simulations:** Sil and k_hat are the method's own `si` and `si.ind`.
- **Real data:** Sil is `mean(silhouette(label, dist(X))[, 3])` and k_hat is
  `length(unique(label))`.

`data_fp` fingerprints the data set (`sum(x * seq_along(x))`), so replicates
can be paired across methods.

All 120 cells now use iteN = 500, so every cell can be paired across methods.
The RK-CCD d = 3, n = 500 scripts had iteN = 50 until 2026-10-18.
- Those published 50-repetition means (slurm-1442239, -1443038) are not
  paired with the other methods' data in the 3- and 5-cluster cells.
- Tables S3 and S7 now carry the 500-repetition values from variant
  `iteN500` (`ITEN500_PROTOCOL.md`). Those rows are the published values for
  this cell: running the cell as `published` would compute the same numbers
  again, in about 7 hours.

## Settings where the code and the paper differ

These are listed in the guide, section 7. The drivers use the code's values,
and changing any of them is the authors' decision:

- the Gaussian cluster centres, 6 apart, and the radius²·I covariance (resolved
  2026-10-31: the authors chose to describe the code's values in the paper, no rerun);
- the per-data-set real-data levels;
- the "or" combination of the mean and median NND tests, where the paper
  describes Holm.

## Descending search index (2026-10-24)

Commit ff343d6 (cherry-picked from 3289f34): `nnccd.radi.descend.one` compares a ball of n-j points with
`env.*.rev[j + 1]`, the quantile entry for n-j points (it read `[j + 2]`, the entry for one point fewer). Only the
descending search (d >= 10) changes; ascending searches, all real-data UN-CCD blocks and the 2d to 5d cells are unaffected.
`UN_CCD_radi_recompute_reference.R` keeps the old index as the historical reference.
- Measured first, paired (same data and tables, only the index changed; the function overridden in the session): 36 cells
  (n = 50, 100, 200) with 500 replicates, and the Gaussian n = 500 cells with the first 50 replicates. The as-is arm equalled
  the published logs to 15 digits in all 36 cells. Scripts and raw outputs are outside the repo (the session's scratchpad `offbyone\`).
- Rerun of the 48 UN-CCD cells with d = 10, 20 (Gaussian and Uniform, n = 50, 100, 200, 500, 2/3/5 clusters, 500 replicates) with
  `10_sim.R UN-CCDs <dist> <d> <n>` (variant `published`; `TR3_CORES` 10 to 24, the results do not depend on it). Logs:
  `_logs/descend_fix_*_20261023.log`. The old rows and files are in `_before_descend_fix_20261023/` (the n = 50 `published` rows
  of these four folders were removed from `sim_replicates.csv` and `sim_done.csv` before the rerun; the `scripts-0.1%` rows,
  also descending with the old index, were left). For the 42 cells the paired run covered, the rerun equals the paired
  "fixed" arm replicate by replicate (ARI, Sil, k_hat exactly).
- Against the published values at 3 decimals, 44 of the 48 cells change in at least one of ARI, Sil, SR. Median absolute change
  0.003 (ARI), 0.001 (Sil), 0.004 (SR); the largest is uniform d = 20, n = 200, 5 clusters: ARI 0.750 to 0.856, Sil 0.460 to 0.490,
  SR 0.284 to 0.516. At 3 decimals ARI rises in 24 cells, falls in 12 and is unchanged in 12.
- `20_gate.R` holds no d >= 10 simulation cell and no descending real-data cell, so its expected values do not change; after the
  fix `00_check_settings.R` passes and the gate passes, 87 of 87 (`_logs/*_after_descend_fix_20261024.log`).
- The noise experiment (d = 3 and 10) was rerun with a new driver whose rows replace the published table; see
  "Noise experiment (2026-10-24)".
- Not rerun: the alpha sweep's config A rows for seeds and breast_cancer (descending; superseded by config S).

## Descending search with sorted neighbour lists (2026-11-02)

Author decision 2026-11-01: adopt the pointer search, so that the paper can state O(n^3) per data set for d >= 10.
In `R/ccds/UN_CCD.R`, when a point leaves the ball, a member whose nearest neighbour it was took a minimum over all other
members, O(n) each, so a data set cost O(n^4) in the worst case. Now `nnccd.sorted.nbrs()` sorts the other points of each
point by distance once per data set (O(n^2 log n) time, n(n-1) integers), and `nnccd.radi.descend.one(..., S)` keeps per
member a position in its list that only moves forward while one centre is searched: O(n^2) per centre, O(n^3) per data set.
The comparison line with the table (`env.*.rev[j + 1]`) is unchanged, so `23_noise.R oldindex` still finds it. The
ascending search, RK-CCDs and KS-CCDs are untouched. Checks, all on this machine, R 4.6.1:
- `38_descend_pointer_check.R`: radii identical to the 2026-10-23 version (copied verbatim into the script) in 115 of 115
  cases: 96 simulated data sets (d = 10, 20; n = 50 to 500), the six real data sets of `31_hd_real.R`, 13 edge cases
  (duplicates, lattice data, n = 4 to 8, low.num 2 and 4, scores). Radius time, all real sets: 511.4 s before, 145.7 s after
  (segmentation 471.0 s to 124.9 s).
- Through the drivers, with outputs in `results/tr3/descend_ptr_20261102/` (the published files' MD5s did not change):
  `31_hd_real.R UN-CCDs` with `TR3_OUT_DIR`: the six data sets equal `real_hd.csv` in n, d, ari, sil and k_hat;
  a copy of `10_sim.R` with its three output paths redirected: the 48 UN-CCD cells with d = 10, 20 (24,000 replicates)
  equal `sim_replicates.csv` in ari, sil, k_hat and data_fp; a copy of `23_noise.R 10 allpoints` with its done file
  redirected and `TR3_SMOKE=1000`: its gate finds 0 differences in the 4,000 published d = 10 noise replicates.
- `00_check_settings.R` passes and `20_gate.R` passes, 87 of 87 (neither has a descending cell).
- Times. Segmentation, one R process: 616.7 s before (`real_hd.csv`), 255.7 s after. Simulations, 24 workers, median
  seconds per data set per (distribution, d) over 1,500 data sets: n = 200, 1.18-1.25 before and 1.23-1.28 after;
  n = 500, 10.97-12.88 before and 9.33-9.78 after; exponent log(t500 / t200) / log(2.5) 2.38-2.58 before, 2.21-2.23 after.
  The "before" times come from other runs (10 to 24 workers), so they are not a controlled comparison.

## Paired tests at n = 200 (2026-11-01/02)

Author request 2026-11-01: paired tests at n = 200 beside n = 50. Per-replicate rows at n = 200 existed only for UN-CCDs
with d = 10, 20 (rerun 2026-10-23), so `10_sim.R <method> <dist> <d> 200` (variant `published`) was run for the other 26
folders, one at a time on this machine: UN-CCDs d = 2, 3, 5 (2026-11-01 18:46-18:47), RK-CCDs (18:47-20:35) and KS-CCDs
(2026-11-02 10:07-18:41; their code does not load `UN_CCD.R`). All exit codes 0. Logs `_logs/n200_<method>_<dist>_<d>d_20261101.log`,
queue `_logs/n200_queue_20261101.log`; the CSVs as they were before are in `_before_n200_20261101/`.
- Check: in each of the 26 logs, the mean ARI, Sil and SR (for KS-CCDs also the chosen delta) for 2, 3 and 5 clusters equal
  the published SLURM log of the folder (`audit_logs.csv`; KS-CCDs batch 145) in all printed digits.
- `35_paired_tests.R` now tests n = 50 and 200 (`TR3_PAIRED_N`), with Holm over the 60 tests of each n
  (`_logs/paired_tests_20261102.log`). The 180 cell means equal the supplement tables at 3 decimals; the n = 50 rows equal
  the earlier run (`_before_n200_20261101/paired_tests_n50.csv`) in all columns.
- n = 200, Holm 0.05: UN vs RK 13 higher, 9 lower, 8 not significant (d >= 10: 12 of 12 higher, mean ARI difference
  +0.024 to +0.326; d <= 5: 9 lower, -0.037 to -0.006, 1 higher, 8 not significant). UN vs KS 2 higher (Gaussian, d = 20,
  3 and 5 clusters), 18 lower (14 with d <= 5, -0.063 to -0.001), 10 not significant. n = 50, unchanged: UN vs RK
  9 / 4 / 17 (higher / lower / not significant), UN vs KS 3 / 19 / 8.
- `36_new_tables.R` writes the paired table with both n (`_logs/new_tables_20261102.log`).

## Gate

2026-10-18: `20_gate.R` passes, 87 of 87 comparisons.
- Simulations: UN-, RK- and KS-CCDs, uniform clusters, d = 2, n = 50. Mean ARI, Sil
  and SR, and the KS-CCD m, match slurm-1442992, -1432233 and -1450256 in all 15
  printed digits. Run times on 24 cores: 2.5 s, 47 s and 2.8 min.
- Real data: iris, seeds and R15 for all three methods.
  - Against the manuscript table: ARI, Sil and k_hat match to 3 decimals.
  - UN- and RK-CCDs against the original script blocks: match to 15 digits.
  - KS-CCDs against Real_Data_KSCCD.html: match to 15 digits, including m.
- The data fingerprints of the three methods agree for all 1500 replicates.

## Audit of the published simulation cells (2026-10-18)

`results/tr3/audit_cells.csv` (one row per method, distribution, d, n and number
of clusters: 360 rows) and `audit_logs.csv` (164 logs) compare Tables S1 to S8
with the SLURM logs. `audit_repro.csv` holds the n = 50 rerun of all 30 folders.
- Every published value has its log. 340 of 360 rows matched to the printed
  3 decimals. In the other 20, one value per row was mistyped or truncated.
  All 20 were corrected in red on 2026-10-18, along with three prose numbers,
  the runtime table and Table S7's label (see the manuscript's ChangeLog.txt).
  A re-audit of the edited supplement shows every row that comes from a log
  equal to it.
- At n = 50, the published log of every folder reproduces to 15 digits,
  including the KS-CCD m. All KS-CCD cells come from the 145xxxx batch (the
  current grid). The 144xxxx KS-CCD logs are from an older grid and are not
  published.
- UN-CCDs, Gaussian, d = 10: the paper's text and Table S8 use the 1% run
  (slurm-1443096 to -99).
  - The scripts loaded the 0.1% table, which gives the unpublished run
    slurm-1451720 to -23. I checked this by rerunning n = 50 with each table.
  - The authors chose 1%, and on 2026-10-18 the four scripts' table line was
    changed to `NN-test-simul_10d_99%.RData`.
  - The 0.1% rerun is kept in the results as variant `scripts-0.1%`.
- Rerun needed: RK-CCDs d = 3, n = 500, uniform and Gaussian (iteN = 50).
  Protocol: `ITEN500_PROTOCOL.md`. Checks: `21_iteN500_check.R`.
  Done 2026-10-18: all checks pass. The new means differ from the old table values by
  at most 0.003, except the Gaussian 5-cluster SR, 0.992 against 1.000 (see the
  protocol's Outcome section). The authors adopted them the same day: Tables S3/S7
  are updated in red, and both scripts now set iteN = 500.

## Real-data rerun (2026-10-19)

All nine data sets with all eight methods. The earlier CCD results are in
`results/tr3/_before_realrerun_20261019/`. Afterwards `20_gate.R` still passes, 87 of 87.
- CCDs: in all 27 cells the driver equals the original script blocks (UN-, RK-CCDs)
  or `Real_Data_KSCCD.html` (KS-CCDs, including m) to 15 digits. Every paper cell
  matches at 3 decimals except UN-CCDs on D31: 0.749 / 0.444 / 31 against the paper's
  0.799 / 0.474 / 31 (main table, Supplement table and its D31 paragraph). Levels 80
  to 95%, z-scoring D31 and the radius search from before 2026-09 do not give 0.799.
  The cause is open; the paper's value predates the git history.
- RK-CCDs run on D31 (3100 points): no covering ball reaches the table's 1000 rows.
- 2026-10-23: the authors adopted a d-based α for UN-CCDs (see "One α" below), and `real.csv`
  now holds the new UN-CCD rows for user_knowledge, wholesale and asymmetric. The old rows are in
  `results/tr3/_before_alpha_adopt_20261023/`. The other bullets of this section describe the 2026-10-19 run.
- 2026-10-23: the RK-CCD rows for user_knowledge and wholesale were rerun with the 5d and 6d tables (the data's own d);
  the old rows are in `results/tr3/_before_rk_dim_20261023/`. user_knowledge 0.058 / 0.123 / 2 (was 0.052 / 0.119 / 2),
  wholesale -0.001 / 0.031 / 44 (was -0.002 / 0.344 / 2). Driver equals script block to 15 digits; the gate passes 87 of 87.
- `user_knowledge` loads as 258 x 5; the paper says d = 6.
- Baselines: DBSCAN, MST and K-means++ (seed 1 by default) are deterministic and
  equal their .html to 15 digits, except DBSCAN on asymmetric, whose .html used
  unscaled data and eps = 50. Spectral and Louvain set no seed.
- Paper cells the code does not give: Louvain asymmetric (the .html's aggregation
  row), MST asymmetric, MST knowledge and wholesale ARI (truncated, not rounded),
  K-means++ wholesale (the script sets k = 6; the paper's values match a Spectral
  run), DBSCAN asymmetric (needs MinPts = 4), and Spectral knowledge, wholesale,
  asymmetric Sil, R15, D31 and Louvain knowledge (not reached by 30 seeds).
- Code against Supplement S2: DBSCAN uses MinPts = 5, not 4; Spectral clusters the
  distance matrix, not a Gaussian similarity. Changing either is the authors' decision.

## One α for UN-CCDs on the real data (2026-10-20)

`ALPHA_PROTOCOL.md` (written before the run; its Outcome section has all tables).
UN-CCDs on all nine data sets at α = 20, 15, 10, 5 and 1%, with tables of the data's
own d. The missing tables were generated once per d (niter 10000, 22 workers, seed
20260927 + d). The published UN-CCD cells reproduce to 15 digits, and no cell failed.
- Mean ARI over the seven paper data sets: 20% 0.523, 15% 0.499, 10% 0.494, 5% 0.492,
  1% 0.483. The script's per-data-set levels give 0.523, with 4 of 7 k̂ correct against
  at most 3 for any single α.
- Mean Sil is highest at 1%. With the script's table d for knowledge and wholesale,
  15% comes first. New table draws move single cells by up to 0.05 ARI (D31) and
  0.14 (breast_cancer at 5%), more than the gaps between levels.
- α was chosen on the data it is scored on, so these numbers describe these nine sets.

2026-10-23: the authors adopted a rule that sets α by the data's d, as in Supplement S1: 15% for
d = 2, 10% for d = 3 and 4, 5% for 5 <= d <= 9 (`ALPHA_PROTOCOL.md`, "Decision"). `real.csv` now holds the
new UN-CCD rows for user_knowledge, wholesale and asymmetric; the old rows are in
`results/tr3/_before_alpha_adopt_20261023/`.

## Search direction by d (2026-10-23)

The authors' rule for UN-CCDs: the radius search is ascending for d < 10 and descending for d >= 10. The
simulation scripts already follow it (2d, 3d and 5d ascend; 10d and 20d descend). Two real-data blocks still
descended: seeds (d = 7) and breast_cancer (d = 9, not in the paper's tables). Both are now `method="ascend"`
in `Real_Data_UNCCD.R`; the tables stay `7d_95%` and `9d_95%`. `settings_realdata.csv` follows, and
`00_check_settings.R` passes. No other direction changes: all nine real-data UN-CCD blocks now ascend.
- Rerun with `11_real.R UN-CCDs seeds,breast_cancer` and `12_real_reference.R UN-CCDs seeds,breast_cancer`
  (logs `_logs/real_unccd_direction_{driver,reference}_20261023.log`). Driver and script block agree to 15
  digits in ARI, Sil and k_hat. The old rows are in `results/tr3/_before_direction_20261023/`; the two stale
  UN-CCD rows were removed from `real.csv`, `real_done.csv` and `real_reference.csv` before the rerun.
- Cross-check: `real_alpha.csv`, config S, level 5 (ascending, same shipped tables) is identical to 17 digits.

| data set | old (descend) ARI / Sil / k_hat | new (ascend) ARI / Sil / k_hat |
|---|---|---|
| seeds | 0.509 / 0.467 / 2 | 0.239 / 0.286 / 2 |
| breast_cancer | 0.075 / 0.208 / 2 | 0.049 / 0.034 / 34 |

New values at 17 digits: seeds 0.23942707729765705 / 0.28604275979153082 / 2; breast_cancer
0.048594143008411979 / 0.033911508115121738 / 34.
- Gate: the first run gave 85 of 87. The two failures were UN-CCDs seeds ARI and Sil against the old
  manuscript values hard-coded in `20_gate.R` (0.509, 0.467). That table was changed to the new seeds values
  (0.239 / 0.286 / 2), which the manuscript now carries, and the gate passes 87 of 87
  (`_logs/gate_after_direction_20261023*.log`).
- The 2026-10-19 bullet "every paper cell matches ... except UN-CCDs on D31" describes the run before this
  change; UN-CCDs on seeds no longer matches the paper's old value either (the manuscript was updated).
- The mean ARI of the adopted d-based α rule changes (see `ALPHA_PROTOCOL.md`, note under "Decision").

## Baseline fixes (2026-10-23)

The authors' decisions, fixed before any new number was seen. Code only; no `.tex` file was touched.
Scripts: `simulations/clustering/Algo_Compare_Clustering/<Method>/Simulation/Real_Data/*.R`
(DBSCAN, K-means++, Louvain, MST fix-cutoff, Spectral). Driver: `13_real_baselines.R`.
The old results are in `results/tr3/_before_baseline_fix_20261023/` (`real_baselines.csv`,
`real_baselines_seeds.csv`, `real_baselines_eps.csv`, `real_verify_baselines.csv`). The old
`real_verify_baselines.csv` was not regenerated; it describes the 2026-10-19 run.

What changed, and why
1. Seeds. Every dataset block of the five scripts starts with `set.seed(1)`, so no cell depends on the
   block order (before, Spectral and Louvain drew from whatever RNG state the previous block left).
   K-means++ also passes `seed = 1` to `KMeans_rcpp()` (its default).
2. K-means++ on wholesale: `k = 3`, the true number of clusters, as the paper's table states (was 6).
3. Spectral. `specc(as.matrix(dist(X)), ...)` is now `specc(as.matrix(X), centers = k)`; the dead
   `similarity_matrix` lines and the header `sigma=1` were removed. The old call: kernlab 0.9.33's
   `specc` method for a plain `matrix` (`getMethod("specc", "matrix")`) starts with `x <- as.matrix(x)`
   and `m <- nrow(x)`, and the documentation says `x` is "the matrix of data to be clustered, ... or a
   kernel Matrix of class 'kernelMatrix'". Only a `kernelMatrix` object is read as a similarity. So the
   rows of the distance matrix were clustered as n-dimensional points with an RBF kernel. With
   `kpar = "automatic"` the source draws `sam <- sample(1:m, floor(mod.sample * m))` (mod.sample = 0.75),
   builds the pairwise distances of that subsample, tries a grid of widths, runs `kmeans(yi, centers,
   iterations)` for each (random starts), and keeps the width with the smallest total within-cluster
   sum of squares: `kernel <- rbfdot((tmpsig[ms]^(-2))/2)`. The width is therefore random (subsample and
   k-means starts) and fixed here by `set.seed(1)`. The block stores specc's inverse width in
   `sigma_used` (`kernlab::kpar(kernlab::kernelf(clusters))$sigma`): k(x, y) = exp(-sigma ||x - y||^2),
   i.e. a Gaussian with bandwidth s = 1/sqrt(2 sigma).
4. DBSCAN. minPts = 5 (`k = 4`) is kept. The typed-in Eps is replaced by `knee_eps(X, k)`, defined in
   the script's header: sort the 4-NN distances d_(1) <= ... <= d_(n) (`dbscan::kNNdist`, the point
   itself excluded; checked against the 5th smallest entry of each row of `dist(X)`), x_i = (i-1)/(n-1),
   y_i = (d_(i) - d_(1))/(d_(n) - d_(1)), i* = argmax (x_i - y_i) (`which.max`, first index on ties),
   Eps = d_(i*) * (1 + 1e-9). The factor removes a knife-edge. d_(i*) is the 4-NN distance of one
   point, and `dbscan` counts neighbours at distance <= Eps, so with Eps = d_(i*) exactly, that
   neighbour sits on the boundary and rounding in the distance decides whether it counts.
   On seeds, point 121 is not a core point at d_(i*) but is at d_(i*)(1 + 1e-9). With
   d_(i*)(1 - 1e-9) the labels change on user_knowledge, wholesale, aggregation and D31 (the number of
   clusters does not); with d_(i*)(1 + 1e-9) no labels change on any of the nine data sets, against
   d_(i*) exactly. A first run with Eps = d_(i*) (the files before the factor) gave the same ARI, Sil
   and k_hat in all nine cells as the runs with the factor, to 17 digits; those files are in
   `results/tr3/_before_eps_margin_20261023/`.
5. Louvain (sigma = 2) and MST (thresh = 2): only the seed. (MST: a label-free threshold rule was tried the same day and not adopted; threshold 2 stays, see "MST threshold rule (2026-10-23)" at the end of this file.) All five methods get the loader's data.
   The loader z-scores (`scale`) iris, seeds, user_knowledge, wholesale, breast_cancer and asymmetric;
   it leaves aggregation, R15 and D31 raw (the `apply(X, 2, scale)` lines are commented out). Louvain
   therefore uses sigma = 2 on z-scores for the six scaled sets and on raw coordinates for the three
   synthetic ones (column sd about 3 to 10). Column sd after the loader: 1.000 for seeds,
   user_knowledge, wholesale, asymmetric; 0.99 to 1.01 for iris (duplicates removed after scaling);
   1.02 to 1.19 for breast_cancer; 8 to 10 (aggregation), 3.3 (R15), 6.8 to 7.3 (D31).

Driver (`13_real_baselines.R`). Header lines are now evaluated as complete expressions (a line that does
not parse alone is joined with the next ones), because the DBSCAN header defines a multi-line function;
before, it evaluated line by line. The sweep mode now replaces the block's own `set.seed(1)` by
`set.seed(s)`, and the kmseed mode replaces `seed = 1)` by `seed = s)`; both stop if the text is not
there exactly once. The sweep mode varies only Spectral and Louvain (K-means++ passes its own seed to
`KMeans_rcpp()`, DBSCAN and MST use no random numbers); only kmseed varies K-means++. The driver's own
variables live in `local()`, so that a loop variable (`i`, `ds`, `s`, `e`, `text`, `src`, `hdr`, `ends`,
`mp`) cannot collide with a variable of a script; the header and the blocks are still evaluated in the
global environment. A header line is joined with the next one only while the parse error is
"unexpected end of input"; any other parse error stops the run. `TR3_BASELINES_OUT` is used as the path
when it starts with a drive letter and a slash or with a slash, otherwise it names a file in
`results/tr3/`. New side file `<output stem>_params.csv` (method, dataset, param, value, seed) with
DBSCAN's `eps` and Spectral's `sigma` and `bandwidth`. The header comments describe all of this.

Runs. Each method on all nine data sets in verbatim mode, one process per method (5 at a time), twice:
`TR3_BASELINES_OUT=real_baselines_fixed_<method>_run<1|2>.csv`; logs in
`results/tr3/_logs/real_baselines_fixed_<method>_run<1|2>.log` (stderr in `.err.log`). No exit codes
were recorded for this batch. Each log has all nine result lines, and the `.err.log` files hold only
package start-up messages (no line with "Error", "Warning" or "Execution halted"). Run-to-run check: the 45 rows of run 2 equal run 1 in ARI, Sil and k_hat as 17-digit
text (so also to 15 digits; max absolute difference 0), and the DBSCAN and Spectral parameter files
are identical too. `real_baselines_fixed.csv` is run 1, `real_baselines_fixed_run2.csv` run 2.

Eps (DBSCAN, 4-NN elbow) and Spectral sigma (specc's chosen inverse width; bandwidth s in brackets)

| data set | Eps old | Eps new | sigma (s) |
|---|---|---|---|
| iris | 0.8 | 0.7890 | 7.280 (0.262) |
| seeds | 0.9 | 1.0983 | 5.651 (0.297) |
| user_knowledge | 1.5 | 1.4568 | 3.281 (0.390) |
| wholesale | 1 | 1.7884 | 0.1076 (2.156) |
| breast_cancer | 2 | 2.2182 | 0.8299 (0.776) |
| aggregation | 1.2 | 1.1236 | 0.1563 (1.789) |
| asymmetric | 0.15 | 0.1262 | 98.10 (0.0714) |
| R15 | 0.35 | 0.3271 | 0.006985 (8.461) |
| D31 | 0.6 | 0.6376 | 0.000388 (35.88) |

Results (ARI / Sil / k_hat, 3 decimals): `results/tr3/real_baselines_fixed.csv`. Against the paper's Table
S11 the cells that match at 3 decimals are: DBSCAN iris; MST iris, seeds, R15, D31; K-means++ iris, seeds,
knowledge, asymmetric, R15, D31 (K-means++ wholesale is now 0.014 / 0.543 / 3 against the paper's
-0.007 / 0.446 / 3); Spectral iris; Louvain iris, seeds, wholesale. The other cells differ. The 2026-10-19 bullets above (Baselines, paper cells the code does not give,
code against Supplement S2) describe the code before these fixes. Left as it was: DBSCAN uses minPts = 5
while the supplement said MinPts = 4 (the supplement was changed to MinPts 5 the same day, manuscript
commit 1baffda); MST and Louvain are not changed.

`20_gate.R` has no baseline check: it compares the CCD cells (simulation cells, iris, seeds and R15) with
the logs, the knit output and the manuscript. That it passes, 87 of 87, after these changes
(`_logs/gate_after_baseline_fix_20261023.log`) shows only that the CCD cells were not disturbed.

Follow-up of the same day (DBSCAN margin, driver changes)
- `knee_eps()` returns d_(i*) * (1 + 1e-9) (see item 4 above); `13_real_baselines.R` got the changes listed
  under "Driver". The DBSCAN 4-decimal Eps values in the table above did not change; the 17-digit values
  are in `real_baselines_fixed_params.csv`.
- Rerun with the changed script and driver, at most 2 processes at a time: DBSCAN twice, MST, K-means++,
  Louvain and Spectral once (`TR3_BASELINES_OUT=real_baselines_fixed_v2_<method>_run<n>.csv`, logs
  `_logs/real_baselines_fixed_v2_*`). Exit codes are in `_logs/real_baselines_fixed_v2_exitcodes.csv`: all
  six are 0; the `.err.log` files have no "Error", "Warning" or "Halted" line.
- ARI, Sil and k_hat of all 45 cells equal the earlier `real_baselines_fixed.csv` as 17-digit text. DBSCAN
  run 1 equals run 2 (rows and parameter file). The new Eps equals the old Eps times (1 + 1e-9) in all nine
  data sets (relative difference from 1e-9 at most 1e-16). The Spectral parameter file is unchanged.
- `real_baselines_fixed.csv` and `real_baselines_fixed_params.csv` now come from the v2 runs. In
  `real_baselines_fixed_run2.csv` the DBSCAN rows are the v2 run 2 and the other four methods' rows are the
  original second run. The earlier versions of these files are in `_before_eps_margin_20261023/`.
- `00_check_settings.R` passes (120 simulation drivers, 27 real-data blocks), and `20_gate.R` passes,
  87 of 87 (`_logs/check_settings_after_baseline_followups_20261023.log`,
  `_logs/gate_after_baseline_followups_20261023.log`).

## MST threshold rule (2026-10-23)

Tried and not adopted. The paper keeps MST with the fixed threshold 2 on every data set (the setting of
the original script), citing Zahn's ratio criterion. How 2 was first chosen is not recorded: the
supplement says it was "chosen after a trial-and-error process ... best overall performance", and nobody
remembers the criterion, so class labels may have been used.

What was done
- The authors specified a label-free rule before any run: on each data set, the threshold is the value of
  `seq(1.00, 3.00, by = 0.01)` (201 values) that maximizes the average silhouette index of the MST partition,
  `mean(silhouette(label, dist(X))[,3])`; ties go to the first maximum (`which.max`); a partition with a
  single cluster gets Sil = -Inf for the selection (so does k = n, where `cluster::silhouette()` returns
  NULL; neither was ever chosen). No class labels.
- The rule was implemented in `Real_Data_MST_silmax_cutoff.R` (a copy of `Real_Data_MST_fix_cutoff.R` with
  `thresh_grid`, `mst_prep()`, `mst_cut()` and `mst_sil()` in the header; the pruning step split from
  `MST_Clustering()` so that the adjacent-edge averages are computed once per data set). Its labels were
  `identical()` to `MST_Clustering()`'s at the chosen threshold in every block, and on iris (all 201 grid
  values), seeds (51), R15 (5), asymmetric (4) and D31 (3). Two runs on all nine data sets (two processes,
  separate part files, exit codes 0 and 0) gave byte-identical files: ARI, Sil, k_hat and the threshold agree
  to 17 digits. The script was removed on 2026-11-03 (author decision); it is in the git history (added in 6027db4).
- Result: it over-segments. Threshold, number of clusters and ARI / Sil, against threshold 2:

| data set | thresh (Sil-max) | k_hat, thresh 2 | k_hat, Sil-max | ARI / Sil, thresh 2 | ARI / Sil, Sil-max |
|---|---|---|---|---|---|
| iris | 1.03 | 6 | 74 | 0.530 / -0.145 | 0.045 / 0.184 |
| seeds | 1.03 | 4 | 101 | 0.000 / -0.364 | 0.036 / 0.179 |
| user_knowledge | 1.04 | 3 | 130 | 0.004 / -0.089 | 0.021 / 0.168 |
| wholesale | 1.00 | 13 | 227 | 0.005 / -0.491 | -0.000 / 0.103 |
| breast_cancer | 1.00 | 6 | 226 | 0.003 / -0.276 | 0.017 / 0.020 |
| aggregation | 1.00 | 12 | 404 | 0.786 / -0.081 | 0.014 / 0.212 |
| asymmetric | 1.04 | 87 | 465 | 0.344 / -0.277 | 0.014 / 0.341 |
| R15 | 1.00 | 65 | 282 | 0.536 / -0.071 | 0.077 / 0.346 |
| D31 | 1.00 | 276 | 1495 | 0.511 / -0.199 | 0.032 / 0.332 |

  The chosen threshold is the lower grid edge (1.00) on 5 of the 9 data sets (wholesale, breast_cancer,
  aggregation, R15, D31); none is at 3.00. The grid was not widened. The undefined-silhouette grid values were
  22 (iris), 41 (seeds), 57 (user_knowledge) and 0 on the other six.
- The result files of this experiment (`real_baselines_fixed_mstrule_MST_run1.csv`, `_run2.csv` and their
  `_params.csv`, param `thresh`) and its logs were removed from the clone on 2026-11-03 (author decision) and
  are not in the repository; the table above is the record. They were never part of `real_baselines_fixed*.csv`.

What the repository has now
- `Real_Data_MST_fix_cutoff.R` is the committed version (`thresh = 2`), unchanged by this experiment, and
  `13_real_baselines.R` is unchanged too (the driver changes made for the experiment, a `thresh` params row
  and a log field, were reverted: the committed script has no `thresh_op`).
- The MST rows of `real_baselines_fixed.csv`, `real_baselines_fixed_run2.csv` and the params file are the
  threshold-2 ones as before (restored from `_before_mst_rule_20261023/`, byte-identical to the files before
  the experiment). A new MST run with the committed script, `real_baselines_fixed_mst2_after_revert_MST_run1.csv`
  (exit code 0; removed with the other part files on 2026-11-03), equals the MST rows of `real_baselines_fixed.csv` in ARI, Sil, k_hat as 17-digit text on all
  nine data sets. The other methods' rows were never changed.
- `00_check_settings.R` and `20_gate.R` pass after the revert (`_logs/check_settings_after_mst_revert_20261023.log`,
  `_logs/gate_after_mst_revert_20261023.log`).
- The paper's Section S2 still says the threshold was chosen by trial and error; the text is the authors' to
  change (Zahn's ratio criterion). The `.tex` files were not touched.

## Noise experiment (2026-10-24)

Result: **the published noise table is not reproduced. On 2026-10-24 the authors decided to adopt this driver's rows for the
paper** (d = 3 from `23_noise.R 3 published`; d = 10 from `23_noise.R 10 published`, the committed code). Both manuscript tables
(`tab:Noise_Level_Clustering`, `tab:Noise_Level_Clustering_Apdx`), the settings in S1.3 and the text quoting the rows were updated, in red.

New rows (ARI / Sil / SR at l = 3, 5, 10, 20%; 1000 replicates, 3 decimals):
- d = 3: 0.973 / 0.617 / 0.998; 0.974 / 0.618 / 0.998; 0.975 / 0.619 / 0.998; 0.965 / 0.616 / 0.985.
- d = 10: 0.953 / 0.384 / 0.971; 0.947 / 0.384 / 0.955; 0.948 / 0.383 / 0.958; 0.947 / 0.385 / 0.955.

Why these settings (all are choices of this driver, not of the original run, which is lost; they are now stated in the paper):
r_noise = 7 and the noise ball centred on the mean of the cluster centres (the demo figures suggest about 6; the d = 3, 20% row
needs 7 or more to stay near the old values, d = 10 does not depend on it); cluster sizes 66/67/67 (the non-noise 3-cluster
rule); ARI, Sil and the choice of k on the 200 cluster points only (the published Sil is flat in l, which the all-points
silhouette does not give). The generator's noise branch (`uniform_cls`, `Gaussian_cls`) adds `colMeans(centers)` to an
n_noise x d matrix by recycling, which misplaces the noise cloud unless all coordinates of the mean are equal (d = 3 and 10 here:
they are not); `23_noise.R` draws the same random numbers but centres the noise with `sweep()` (`center = recycle` gives the
generator's behaviour). `R/general_functions/` was not edited. The details of the search and the failed validation follow.

- Driver: none exists. Searched: `simulations/clustering/`, `revision_experiments/`, all of `git log --all` (names with
  noise/cont/contamination), SLURM logs. The only noise drivers (`simulations/outlier_detection/*/Simulation/Noise/`) belong to the
  outlier-detection papers (contamination 5%, success and false-positive rates). `Gaussian_cls(noise = TRUE, r_noise, n_noise)`
  in `R/general_functions/Uni-Gau_cls.R` is the generator, and no committed script calls it with `noise = TRUE`.
  `23_noise.R` was written from Supplement S1.3. `CLUSTERING_REPRODUCTION_GUIDE.md` section 8.4 already noted the gap.
- Settings. From S1.3: d = 3, 10; l = 3, 5, 10, 20% (n_noise = round(200 l)); 200 cluster points; 3 equal clusters; centres
  (3,..,3), (9,3,..,3), (3,9,3,..,3); covariance Delta I, 1000 repetitions. From the non-noise driver
  (`Gaussian/<d>d/<d>d_235clusters_n200.R`, via `settings_simulation.csv`): seed 2024, sigma 1, radius factor U[0.8, 1.2],
  sizes 66/67/67, table `NN-test-simul_3d_90%` ascending at d = 3, `NN-test-simul_10d_99%` descending at d = 10, low.num 3,
  niter 1000. Not stated anywhere, chosen: r_noise = 7, the noise ball centred on the mean of the cluster centres (`sweep`; the
  generator as written recycles that vector down an n_noise x d matrix and misplaces the noise, available as `center = recycle`),
  and the scoring (below). The 2-d demo figures suggest r_noise of about 6.
- Scoring. The published Sil is flat in l (0.619 to 0.617 at d = 3) and ARI barely moves. The clustering runs on all points, but only
  scoring on the 200 cluster points reproduces that: `nnccd.silhouette(ind = 1:200)` for the choice of k and for Sil, and ARI on the
  200 points. With the silhouette over all points Sil falls from 0.60 to 0.51 (d = 3) as l goes 3% to 20%; with the noise as a fourth
  class ARI falls to 0.72. Both are in the CSV (`sil_all`, `k_all`); no committed script uses `ind`, so this is an inference.
- Validation (1000 replicates per cell, `results/tr3/noise_replicates.csv`; seeds of the published run unknown, so the criterion is
  Monte Carlo error: z = difference / (SE sqrt 2), SE from the replicates, the published run taken to have the same size).
  Same data for old and fixed index.

  | cell | ARI mine (pub, z) | Sil mine (pub, z) | SR mine (pub, z) |
  |---|---|---|---|
  | d = 3, 3% | 0.973 (0.977, -2.0) | 0.617 (0.619, -1.1) | 0.998 (1.000, -1.0) |
  | d = 3, 5% | 0.974 (0.977, -1.8) | 0.618 (0.619, -0.7) | 0.998 (1.000, -1.0) |
  | d = 3, 10% | 0.975 (0.975, +0.1) | 0.619 (0.618, +0.9) | 0.998 (0.998, 0.0) |
  | d = 3, 20% | 0.965 (0.971, -1.9) | 0.616 (0.617, -0.6) | 0.985 (0.982, +0.6) |
  | d = 10 old index, 3% | 0.950 (0.959, -2.2) | 0.383 (0.386, -1.8) | 0.965 (0.973, -1.0) |
  | d = 10 old index, 5% | 0.949 (0.961, -3.0) | 0.384 (0.384, +0.1) | 0.961 (0.977, -1.8) |
  | d = 10 old index, 10% | 0.946 (0.954, -1.8) | 0.383 (0.384, -0.7) | 0.956 (0.960, -0.4) |
  | d = 10 old index, 20% | 0.946 (0.947, -0.3) | 0.385 (0.385, -0.1) | 0.955 (0.942, +1.4) |

  Verdict: not reproduced. ARI is below the published value in 7 of 8 cells (z -1.8 to -3.0 in 6), and at d = 10 the published SR
  falls from 0.977 to 0.942 with l while this run stays near 0.96. Differences are small in size (ARI up to 0.012, SR up to 0.016).
  Also the published d = 3, 3% row is above the published no-noise 3-cluster cell (0.972 / 0.618 / 0.996, `slurm-1451714.out`). The
  d = 10 results hardly depend on r_noise (3% to 20% rows with r_noise 4 to 10, ARI 0.94 to 0.95 in all, old index); at d = 3 the 20%
  row needs r_noise of 7 or more to match (r_noise 3: ARI 0.91). Exploration grid (r_noise 3 to 10, both centrings at d = 3, r_noise
  4 to 10 asis at d = 10) was run with a scratch script outside the repo and is not saved here.
- For information only, not for the manuscript: the same data with the committed (fixed) index, d = 10, l = 3/5/10/20%:
  ARI 0.953 / 0.947 / 0.948 / 0.947, Sil 0.384 / 0.384 / 0.383 / 0.385, SR 0.971 / 0.955 / 0.958 / 0.955. Paired against the old
  index, k_hat differs in 16 to 20 of 1000 replicates per level and the mean ARI changes by +0.002, -0.002, +0.001, +0.001 (SE 0.0017):
  the index fix does not move the d = 10 noise row beyond Monte Carlo error.
- Runs: `23_noise.R 3 published`, `23_noise.R 10 oldindex`, `23_noise.R 10 published` (defaults r_noise 7, sweep, 1000);
  49 s, 5.1 min, 5.1 min on 24 cores. Logs `_logs/noise_d{3,10}_*_20261024.log`. The old index is an override of
  `nnccd.radi.descend.one` in the global environment (`[j + 2]`), checked on the workers by a column in the output; `R/ccds/` is untouched.
- An exact reproduction of the old table would need what the original run did (r_noise and centre, the scoring, the cluster
  sizes and the seed); that is unknown, so the authors adopted this driver's rows instead (see the top of this section).

## Runtime comparison (2026-10-24)

Re-timing of `tab:runtime_comparison` of the manuscript (total run time in minutes of the d = 5 uniform simulation, 500 repetitions,
2/3/5 clusters, n = 50, 100, 200, 500; UN-, RK- and KS-CCDs).

**Where the published times came from.** There was no timing script. The published values are the `Time difference` lines of the
d = 5 uniform SLURM logs (`simulations/clustering/<M>_CCDs/Simulations/Uniform/5d/slurm-*.out`: UN 1443072-75, RK 1442240-42 and
1443104, KS 1450264-67), which are in the repository since `82b22db` (2025-09-12), so they were produced with the code before the
incremental radius search of `059c798` (2026-09-23); the UN-CCD row used the recompute search, RK- and KS-CCD code is unchanged
since. The author confirmed on 2026-10-24 that these runs were on cluster nodes with 2x Intel Xeon Gold 6248R (24 cores per CPU, 48
per node), not on the i9-13900K workstation that the manuscript named. The number of cores allocated to each job is unknown: the
repository has no sbatch files, and the drivers start `detectCores()` workers whatever the allocation was.
What the logs time: from `t1 = Sys.time()` (after the method, the packages and the quantile table are loaded) to the end of the driver:
data generation, three `makeCluster(detectCores())` start-ups, the clustering of 3 x 500 data sets with ARI and Sil, and for KS-CCDs
the 200-value grid search for each number of clusters. Table loading is excluded. The drivers use 2, 3 and 5 clusters (the text says
"2 to 5").

**Script.** `22_runtime.R <method> <n> <rep> [iteN]` runs that interval once in a fresh R process, with the settings of the driver's row
in `settings_simulation.csv` and the methods called as the drivers call them (`UNCCD_clustering`, `RKCCD_clustering`,
`KSCCD_clustering`). The RK table is loaded with its raw draws (`Kest.m`) and exported to the workers inside the timed interval, as in
the drivers (`TR3_RUNTIME_DROP_KEST=1` drops them). It appends one row per run to `results/tr3/runtime_reps.csv` (seconds, split into data
generation and the three clustering blocks, the means of ARI, Sil and SR, KS-CCDs' best delta, cores, R version, CPU, RAM, OS, power
plan, `valid`). `22_runtime_monitor.ps1` runs beside it and samples every 3 s the R processes that are not the driver's (the driver,
its launcher chain and the workers are excluded); the script stops at the start if one is running, and a run is `valid = FALSE`
if one appears during it (none did). `22_runtime_run.ps1 -Reps 3` runs pass r over all 12 cells (cheapest first) before pass r + 1;
KS-CCDs at n = 500 is run in pass 1 only (`-Once`, the author's decision). `22_runtime_summary.R` writes `runtime_summary.csv` (median, min,
max, published value, ratio, and whether the means of ARI, Sil and SR, and KS's best delta, equal the log's).

**Settings and machine.** 24 workers (`detectCores()`, as the drivers), one timing process at a time,
nothing else running (checked by the monitor). Intel Core Ultra 7 270K Plus (24 cores, no hyper-threading), 64 GB RAM (63.5 GB usable), Windows 11 Pro
10.0.26300, R 4.6.1 (2026-06-24 ucrt), power plan Balanced. Quantile tables from disk (`NN-test-simul_5d_95%`, `RK-test-simul_5d_99%`).
Loop started 2026-10-24 02:27; pass 1 ran to 15:51 (KS-CCDs n = 500 alone took 11.5 h), passes 2 and 3 ended at 19:30. Log
`_logs/runtime_run_20261024.log`; the first two start attempts were stopped by the guard because the monitor counted the Rscript.exe
launcher (Rscript.exe, cmd.exe, x64 Rscript.exe) as foreign, which the monitor now handles.

**Results (minutes; median of 3 runs, KS-CCDs n = 500 one run; the three runs in brackets; published = SLURM log).**

| Method | n | new median | runs | published | new / published |
|---|---|---|---|---|---|
| UN-CCDs | 50 | 0.055 | 0.055, 0.055, 0.054 | 0.153 | 0.36 |
| UN-CCDs | 100 | 0.088 | 0.088, 0.088, 0.087 | 0.294 | 0.30 |
| UN-CCDs | 200 | 0.239 | 0.242, 0.239, 0.236 | 0.649 | 0.37 |
| UN-CCDs | 500 | 1.921 | 1.926, 1.917, 1.921 | 4.379 | 0.44 |
| RK-CCDs | 50 | 0.849 | 0.854, 0.849, 0.835 | 4.603 | 0.18 |
| RK-CCDs | 100 | 1.277 | 1.330, 1.277, 1.277 | 6.457 | 0.20 |
| RK-CCDs | 200 | 3.844 | 4.019, 3.825, 3.844 | 16.990 | 0.23 |
| RK-CCDs | 500 | 39.283 | 39.613, 39.283, 39.207 | 78.848 | 0.50 |
| KS-CCDs | 50 | 2.686 | 2.738, 2.686, 2.682 | 38.552 | 0.07 |
| KS-CCDs | 100 | 9.285 | 9.465, 9.230, 9.285 | 41.297 | 0.22 |
| KS-CCDs | 200 | 49.259 | 49.376, 49.133, 49.259 | 86.163 | 0.57 |
| KS-CCDs | 500 | 692.699 | 692.699 | 987.513 | 0.70 |

- Every valid run reproduces its log's mean ARI, Sil and SR (and KS-CCDs' best delta) to 12 digits: the speed changes, the results
  do not (`runtime_summary.csv`, column `results_equal_log`, TRUE in all 12 cells).
- The runs differ by at most 2% from the median except RK-CCDs n = 100 and n = 200 (first run 4% to 5% slower).
- This machine is 1.4 to 14 times faster than the cluster runs. RK- and KS-CCD code is the same as then, so those ratios (0.07 to
  0.70) measure the machine and the allocation, which is unknown; the UN-CCD ratios (0.30 to 0.44) also contain the incremental search
  of 059c798 and cannot be split from them. The ratio rises with n (0.07 to 0.70 for KS-CCDs, 0.18 to 0.50 for RK-CCDs): the small
  cells are dominated by the cluster start-ups and, on the cluster nodes, by the unknown allocation.
- The ordering UN < RK < KS holds in all four columns. At n = 500: RK/UN = 20, KS/UN = 361, KS/RK = 18 (published 18, 225 and 12.5).
- The per-dataset d = 10, 20 timings quoted in Supplement S3 (`sim_replicates.csv`, `seconds`, one worker of a parallel run) are
  a different measurement and are not part of this table.

## Greedy selection as the code runs it (2026-10-25)

The paper described the CCD pipeline as: find an approximate minimum dominating set with the outdegree
greedy algorithm, build the intersection graph of those balls, then select centers with the score-based
greedy algorithm. The code does not run the first step as described.
- All three wrappers (`methods/clustering/{UN,RK,KS}-CCDs.R`) pass `dom.method = "greedy2"`.
  `dominate.mat.greedy2` (`R/ccds/ccdfunctions.R`) marks only the chosen vertex as covered, so it returns
  every vertex ranked by outdegree. On a 12-vertex random digraph it returned all 12; the covering prefix
  of the same ranking had 4.
- `nnccd_clustering_quantile` (`UN_CCD.R`), `RKCCD_correct_quant` (`RK_CCD_New.R`) and `ksccd.clustering`
  (`ccd_ks_NEW.R`) then build the intersection graph on all n balls and select centers with
  `dominate.mat.ks(M.dom, score[D])`: the uncovered ball with the largest fixed score (points caught) is
  taken and every ball sharing a point with it is marked covered.
- `25_greedy_prefix.R` reruns the seven paper data sets with the shipped step and with the covering prefix
  (the step as the paper described it), with `dominate.mat.greedy2` replaced in the global environment.
  - UN-CCDs: identical labels on all seven; the shipped run equals `real.csv` within 5e-16. Centers passed
    to the silhouette step shrink (iris 21 to 13, asymmetric 49 to 28, D31 162 to 152) without changing a
    partition. The prefix only removes balls that `dominate.mat.ks` processes after all the others.
  - RK-CCDs (prefix only): ARI and Sil equal the published rows of `real.csv` within 4.4e-16, same k-hat
    on all seven (labels not compared).
  - KS-CCDs: not run (the 200-value grid search takes hours); the code path is the same.
- The manuscript now describes the step as the code runs it (manuscript commits c82e6de, 555666a).
  Nothing in `R/` or `methods/` was edited.

## New experiments (2026-10-25)

Experiments 1, 2, 3 and 5 of the revision plan, asked for by the authors on 2026-10-25 for the
KBS rebuttal. Design, decisions and outcomes are in `NEW_EXPERIMENTS_PROTOCOL.md` (written before
the runs). Each driver was written by one agent and checked by a separate verifier agent before
its full run, and the experiments ran one at a time. The manuscript reports them in red (blue
until 2026-10-31): Table 6 (experiment 1) and Supplement Section S4 (experiments 2, 3 and 5).
Experiment 7 (the five baselines on the simulated and noise data) was added on 2026-10-31 under
the same rules.

From the clone root, in PowerShell:

```powershell
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
& $R revision_experiments/tr3/31_hd_real.R UN-CCDs             # or RK-CCDs, KS-CCDs; [datasets]
& $R revision_experiments/tr3/32_hd_baselines.R all            # or one method; [datasets] [run|check|show]
& $R revision_experiments/tr3/33_type1.R test                  # then ccd, then summary
& $R revision_experiments/tr3/34_ablation_direction.R run      # [reps] [cells]; then summary
& $R revision_experiments/tr3/35_paired_tests.R
& $R revision_experiments/tr3/37_sim_baselines.R all check     # the gate; then: all sim, all noise
& $R revision_experiments/tr3/36_new_tables.R                  # tex/new_*.tex and a summary
```

- `TR3_OUT_DIR` sends the output of 31 to 36 to another folder (36 then reads its inputs there
  too). `TR3_SMOKE_GRID` (31 only) runs KS-CCDs on a coarse grid and labels the rows `hd_smoke`;
  31 stops if it is set while the output folder is `results/tr3`.
- Experiment 1 data: `30_hd_load.R` reads wine, vehicle, segmentation and wdbc from
  `revision_experiments/tr3_data/` (UCI downloads of 2026-10-25, md5 in `MANIFEST.csv`; the
  zip files and the files extracted from them, which the loader reads, are committed since
  2026-11-03) and hepatitis and lymphography from `data/clustering/`.
- Experiment 1 runs: the CCDs from 2026-10-25 23:26 to 2026-10-26 00:37 (UN-CCDs 11.1 min,
  RK-CCDs 1.7, KS-CCDs 58.0), then the five baselines (1.8 min); all exit codes 0.
  `32_hd_baselines.R ... check` reproduces `real_baselines_fixed.csv` on the nine older data sets.
- DBSCAN's k_hat: the original block counts clusters as distinct labels minus one (the noise
  label), so it would be one too low on a data set without noise points; all six new sets have
  noise points, and so do the nine older ones (`37_sim_baselines.R all check`, 2026-10-31).
  Experiment 7 counts the labels other than 0; DBSCAN found no noise point in 8827 of its 60000
  simulated data sets.
- Experiment 7 runs (2026-10-31): `all sim` 19:28 to 20:05 (36.5 min, 24 workers), `all noise`
  to 20:07 (2.5 min); exit codes 0 (`_logs/simbase_exitcodes_20261031.csv`). The verifier's own
  code, the 510-value real-data check and the fingerprints are described in the protocol outcome.
- Experiment 8 (2026-10-31): `23_noise.R 3 allpoints` then `23_noise.R 10 allpoints`, 21:53 to
  22:00, exit codes 0 (`_logs/noise_allpoints_20261031.log`); `TR3_SMOKE=3` runs the first three
  replicates of each level on the full data and writes nothing. Main Table 4 now shows UN-CCDs with
  k chosen on all points, as the baselines; Supplement Table S9 shows the same rule, one row per d
  (the cluster-points rows left the supplement on 2026-11-01, manuscript commit e4cc060; they stay in
  `noisebase_summary.csv`).
