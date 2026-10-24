# tr3: UN-CCD clustering paper

Revision experiments for "Clustering with Uniformity- and Neighbor-Based Random
Geometric Graphs" (Shi, Ceyhan and Billor; rejected by KBS as KNOSYS-D-26-08010).
The manuscript is at `G:\Submissions\TR3\TR3_10182026\`. The map from each paper
result to its original script is `CLUSTERING_REPRODUCTION_GUIDE.md` at the clone root.

These drivers replace the 120 copy-paste simulation scripts and the three
real-data scripts for new work. The originals print averages to SLURM logs and
save the workspace. The drivers run the same code through one interface, with
the same settings and data, and write every replicate to a CSV. Nothing in
`methods/`, `R/` or `data/` is edited. In `simulations/`, the following settings were
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

Reproduction gate (`20_gate.R`): see the end of this file. No new experiment runs
until the gate passes, and each one gets a `*_PROTOCOL.md` here before it runs.

## Files

| File | What it does |
|---|---|
| `settings_simulation.csv` | The settings of the 120 simulation drivers, as the code ran them for the published numbers (not as the paper states them): seed, iteN, cls_dis, sig, radius range, size of the first of three clusters, quantile table, RK level, UN direction, KS grid. |
| `settings_realdata.csv` | Table, level, direction and grid per method and real data set (27 rows). |
| `00_check_settings.R` | Re-reads all 120 drivers and the three real-data scripts and stops if any value differs from the two tables, or if the drivers differ in anything else. |
| `tr3_lib.R` | Shared code: `ccd_cluster()`, data generation in the scripts' random-number order, the patched loader, CSV checkpointing. |
| `10_sim.R` | One simulation cell group (method, distribution, d, n). |
| `11_real.R` | CCD rows of the real-data table. |
| `12_real_reference.R` | Runs the original real-data script blocks verbatim, to get full-precision reference values. |
| `20_gate.R` | Compares the drivers with the original logs, the knit output and the manuscript. |
| `ITEN500_PROTOCOL.md`, `21_iteN500_check.R` | Variant `iteN500`: RK-CCDs at d = 3, n = 500 with 500 repetitions, and its checks. |
| `13_real_baselines.R` | Baseline rows of the real-data table (DBSCAN, MST, K-means++, Spectral, Louvain): runs each data-set block of the original script verbatim, once or over seeds. |
| `ALPHA_PROTOCOL.md`, `14_alpha_tables.R`, `15_real_alpha.R`, `16_alpha_summary.R` | One α for UN-CCDs on all real data sets (20% to 1%): generates the missing NN tables, runs every cell, checks the reproduction and scores the levels. |

## Running

From the clone root, in PowerShell (Rscript segfaults under Git Bash):

```powershell
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
& $R revision_experiments/tr3/00_check_settings.R
& $R revision_experiments/tr3/10_sim.R UN-CCDs uniform 2 50          # [clusters] [variant]
& $R revision_experiments/tr3/11_real.R RK-CCDs iris,seeds,R15
& $R revision_experiments/tr3/12_real_reference.R UN-CCDs iris,seeds,R15
& $R revision_experiments/tr3/20_gate.R
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

## Results (`revision_experiments/results/tr3/`, not tracked by git)

| File | One row per |
|---|---|
| `sim_replicates.csv` | replicate: variant, method, dist, d, n, n_clusters, iteN, rep, ari, sil, k_hat, seconds, param (KS: chosen m), data_fp |
| `sim_done.csv` | finished cell (what reruns skip) |
| `sim_ks_search.csv` | KS grid value per cell: delta, mean_sil, mean_ari, sr, seconds |
| `real.csv`, `real_done.csv` | method and data set: n, d, ari, sil, k_hat, seconds, param, table |
| `real_ks_search.csv` | KS grid value per data set: si |
| `real_reference.csv` | output of `12_real_reference.R` |
| `real_baselines.csv`, `real_baselines_seeds.csv` | output of `13_real_baselines.R`: verbatim runs, and seeds 1 to 30 (both from before the baseline fixes of 2026-10-23; copies in `_before_baseline_fix_20261023/`) |
| `real_baselines_fixed.csv`, `real_baselines_fixed_run2.csv`, `real_baselines_fixed_params.csv` | the baselines after the fixes of 2026-10-23: first run, second run (identical), and the Eps / sigma each block chose. Part files `real_baselines_fixed_<method>_run<1,2>[_params].csv` |
| `real_verify_ccd.csv`, `real_verify_baselines.csv` | real-data rerun of 2026-10-19 against the paper and the knit output, per cell and metric |
| `real_alpha.csv`, `alpha_summary.csv` | α sweep: config, level, data set (ARI, Sil, k_hat, table, source); the criteria per config and level. Generated tables in `alpha_tables/` |
| `gate.csv` | gate comparison |

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

- the Gaussian cluster centres, 6 apart, and the radius²·I covariance;
- the per-data-set real-data levels;
- the "or" combination of the mean and median NND tests, where the paper
  describes Holm.

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
5. Louvain (sigma = 2) and MST (thresh = 2): only the seed. All five methods get the loader's data.
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
while the supplement says MinPts = 4; MST and Louvain are not changed.

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
