# Reproducing the UN-CCD clustering paper with this codebase

This guide is for reproducing the experiments of *Clustering with Uniformity- and Neighbor-Based Random Geometric Graphs* (Shi and Ceyhan, arXiv:2501.06268; the current manuscript is `G:\Submissions\TR3\TR3_10182026\UNCCD4Clustering.tex` with its `Supplementary_Material.tex`).
It was written on 2026-10-18 against `main` at `bdcf27e`, the commit this clone is on.
The file is not tracked by git.

Section numbers such as "Table S4" refer to that manuscript and supplement.

---

## 1. Summary

- The clustering code, data and drivers are in `methods/clustering/`, `R/ccds/`, `R/general_functions/`, `data/clustering/` and `simulations/clustering/`.
  Nothing else in the repository is needed; `revision_experiments/` belongs to two other papers.
- Every quantile table that any clustering driver loads is in this clone (section 5).
- The printed results of every simulation run are kept next to its driver as `slurm-*.out` files.
  These logs are the source of the numbers in the paper: the d = 2, n = 50 log of the UN-CCD uniform runs is exactly the n = 50, d = 2 block of Table S4.
- Three things stop the committed drivers from running (section 6.1):
  1. The 40 UN-CCD simulation drivers call `UNMCCD_outlier()`, which no file they source defines.
     Use `UNCCD_clustering()` instead.
  2. `data/clustering/Real_Data_Collection.R` changes into a folder that no longer exists.
     Point it at `data/clustering/`.
  3. Four R packages were missing on this machine (installed 2026-10-19): `readxl` (the real-data loader), and `factoextra`, `ClusterR` and `kernlab` (the baselines).
- The noise experiment (main-text noise table, Supplement S1.3) has no driver in the repository, and no figure script is in the repository (section 2).
- Several settings in the drivers differ from what the paper says (section 7).
  Reproduce the code as it is first; decide afterwards which side to change.

- Checked on this machine (section 9): with the current code and tables, the UN-CCD driver for uniform clusters at d = 2, n = 50 (after the one-word fix) and the RK-CCD driver for the same cell (unchanged) reproduce their logs to every printed digit, and the UN-CCD and RK-CCD results on iris, seeds and R15 match the paper's real-data table.

---

## 2. Paper item → code

All paths are relative to the clone root.
`<d>` is one of 2, 3, 5, 10, 20 and `<n>` one of 50, 100, 200, 500.

| Paper item | Driver(s) | Where the published numbers are |
|---|---|---|
| KS-CCDs, uniform clusters: Table S2, best density Table S1, line-plot figure | `simulations/clustering/KS_CCDs/Simulations/Uniform/<d>d/<d>d_235clusters_n<n>.R` | `slurm-*.out` in the same folder (prints ARI, Sil, SR and "best delta" for 2, 3 and 5 clusters) |
| RK-CCDs, uniform: Table S3, figure | `simulations/clustering/RK_CCDs/Simulations/Uniform/<d>d/...` | same |
| UN-CCDs, uniform: Table S4, figure | `simulations/clustering/UN_CCDs/Simulations/Uniform/<d>d/...` | same |
| KS-CCDs, Gaussian: Table S6, best density Table S5, figure | `simulations/clustering/KS_CCDs/Simulations/Gaussian/<d>d/...` | same |
| RK-CCDs, Gaussian: Table S7, figure | `simulations/clustering/RK_CCDs/Simulations/Gaussian/<d>d/...` | same |
| UN-CCDs, Gaussian: Table S8, figure | `simulations/clustering/UN_CCDs/Simulations/Gaussian/<d>d/...` | same |
| Runtime table (`tab:runtime_comparison`, d = 5) | none; the values are the `Time difference` lines of the d = 5 uniform logs (KS-CCDs: the 145xxxx runs). Since 2026-10-18 all twelve match to two decimals. Before that, RK-CCDs n = 500 and KS-CCDs n = 200 and 500 were hours rounded to two decimals, then times 60 | `simulations/clustering/*_CCDs/Simulations/Uniform/5d/slurm-*.out` |
| Noise experiment (`tab:Noise_Level_Clustering`, S1.3) | **no driver in the repository** | not in the repository |
| Real data (`Real_Data_Result_Clustering`), CCD rows | `simulations/clustering/{KS,RK,UN}_CCDs/Simulations/Real_Data/Real_Data_{KSCCD,RKCCD,UNCCD}.R` | `Real_Data_KSCCD.html` (knit output) for KS-CCDs; none kept for RK and UN |
| Real data, baseline rows | `simulations/clustering/Algo_Compare_Clustering/{DBSCAN,K-means++,Louvain,MST,Spectral_Clustering}/Simulation/Real_Data/*.R` | the `.html` knit output next to each script (from before the baseline fixes of 2026-10-23; the new numbers are in `revision_experiments/results/tr3/real_baselines_fixed.csv`, see `revision_experiments/tr3/README.md`) |
| Line-plot figures (`Lineplots_Clustering_*.png`) | not in the repository; `G:\Submissions\TR3\TR3_10182026\figure_unused\Line_Plots_Clustering.ipynb` draws the same six figures (it saves them as `..._LargerFont.png`) | |
| Demo figures (`Demo_*Cls.png`, `Demo_Noise_Clustering_*.png`, `UNCCD1-4.png`, `asymmetric/R15/D31.png`) | not in the repository | |

Each simulation folder holds four drivers (one per n) and four or eight SLURM logs.
The log names do not say which n they belong to; match a log to its n by content (the paper's table value) or by run time (larger n takes longer).
Folders with eight logs were run twice (job numbers 144xxxx and 145xxxx); which run the paper reports has to be read off by comparing with the table.

---

## 3. Code layout

### 3.1 Method wrappers (`methods/clustering/`)

| File | Function | What it calls |
|---|---|---|
| `UN-CCDs.R` | `UNCCD_clustering(datax, simul, method = "ascend", niter = 1000, ...)` | `nnccd_clustering_quantile()` with `low.num = 3`, `dom.method = "greedy2"`, then `nnccd.silhouette()` (both in `R/ccds/UN_CCD.R`) |
| `RK-CCDs.R` | `RKCCD_clustering(datax, simul = NULL, quant = 0.99, ...)` | `RKCCD_correct_quant()` with `r.seq = 10`, `dom.method = "greedy2"` (`R/ccds/RK_CCD_New.R`) |
| `KS-CCDs.R` | `KSCCD_clustering(datax, m = 10, ...)` | `ksccd.clustering(m = m, sequential = FALSE)` and `ksccd.silhouette()` (`R/ccds/ccd_ks_NEW.R`, `R/ccds/ccdfunctions.R`) |

All three return a list with `label` (cluster of every point), `si` (the average silhouette of the chosen partition) and `si.ind` (number of clusters chosen).
The drivers compute ARI with `mclust::adjustedRandIndex(x$label, label_T)`, Sil as `mean(x$si)` and SR as `mean(x$si.ind == N)`.

`methods/outlier_detection/UN-MCCD.R` defines `UNMCCD_outlier()`, the outlier detector of the Neurocomputing paper; it is not part of this paper.

### 3.2 The algorithm, step by step (`R/ccds/UN_CCD.R`)

1. `nnccd.radi()` finds a covering-ball radius for every point (section 4).
2. Every point catches the points within its radius; `dominate.mat.greedy2()` (`R/ccds/ccdfunctions.R`) finds an approximate minimum dominating set of the resulting digraph.
3. An intersection graph joins two dominating balls that catch a common point; `dominate.mat.ks()` orders its dominating set by the number of points each ball catches.
4. `nnccd.silhouette()` takes the first 2, 3, ... of those balls as cluster centres, assigns every point to the ball with the smallest distance-to-radius ratio (`nnccd.clustering.nonvalid()`), and keeps the partition with the largest average silhouette.
   Every point gets a cluster; there is no noise label.

RK-CCDs (`R/ccds/RK_CCD_New.R`) and KS-CCDs (`R/ccds/ccd_ks_NEW.R`) share steps 2 to 4 and differ only in step 1.

- RK-CCDs: `ccd.Kest.edge.quantile()` grows the ball one neighbour at a time and stops when the edge-corrected K-function estimate, on a grid of `r.seq = 10` radii, exceeds the table's upper quantile at any grid point.
- KS-CCDs: `ccd.nonseq_scale()` chooses, for each point, the radius r that maximises (number of points within r) − (m·r)^d.
  So `m` is the paper's scaled density √[d]δ, and δ = m^d.
  The drivers search m over `c(seq(0.01, 1, 0.01), seq(1.1, 10, 0.1), seq(11, 20, 1))`, 200 values, and keep the one with the largest mean silhouette over the 500 data sets.

### 3.3 Data generation (`R/general_functions/Uni-Gau_cls.R`)

- `uniform_cls(n_cls, d, centers, rad_range, noise = FALSE, r_noise, n_noise)`: cluster i is a uniform sample of `n_cls[i]` points in a d-ball of radius U(`rad_range`) around `centers[i, ]`.
- `Gaussian_cls(n_cls, d, centers, sigma, rad_range, noise = FALSE, ...)`: cluster i is N(`centers[i, ]`, sigma · radius² · I) with radius ~ U(`rad_range`).
- With `noise = TRUE` both add `n_noise` uniform points in a ball of radius `r_noise` around the mean of the centres.
  No committed driver uses this option; it is the natural starting point for the missing noise driver.
- The file also loads `spatstat` and `hypercube` and defines `Uni.Gau_cls()` (Matérn/Thomas processes), which the clustering paper does not use.

### 3.4 Real-data loader (`data/clustering/Real_Data_Collection.R`)

It reads iris, seeds, wholesale, breast cancer, aggregation, asymmetric, R15, D31 and user knowledge, z-scores the columns (not for aggregation, R15 and D31; asymmetric *is* scaled although the paper says the synthetic sets are not), drops rows with missing values and duplicated rows, and leaves each data set as a matrix whose last column is the class.
The data files in `data/clustering/` are byte-identical to the ones in the original 2025 commit (`82b22db`).

What the loader produces (with its folder fixed, section 6.1):

| Object | n | d | Classes | Paper |
|---|---|---|---|---|
| `iris` | 147 | 4 | 3 | n = 150; three duplicate rows are dropped |
| `seeds` | 210 | 7 | 3 | same |
| `user_knowledge` | 258 | 5 (the code keeps columns 1 to 5) | 4 levels of `UNS` | d = 6 |
| `wholesale` | 440 | 6 (Channel and Region dropped; Region is the class) | 3 | d = 8 |
| `breast_cancer` | 449 | 9 | 2 | not reported |
| `aggregation` | 788 | 2 | 7 | not reported |
| `asymmetric` | 1000 | 2 | 5 | same |
| `R15` | 600 | 2 | 15 | same |
| `D31` | 3100 | 2 | 31 | same |

---

## 4. The NND-based Monte Carlo test as implemented

`nnccd.radi()` (`R/ccds/UN_CCD.R`) takes the data, the table `simul`, a direction `method` and `low.num = 3`.
For a centre x and a candidate radius r (the distance to its j-th nearest neighbour), the ball holds the j − 1 points nearest to x, *without x itself*.
The statistics are the mean and the median of those points' nearest-neighbour distances within the ball, divided by r.

- `method = "ascend"`: try j = 3, 4, ... in turn; at the first j where the mean is below `simul$average[j - 1]` **or** the median is below `simul$median[j - 1]`, stop and return the previous radius (0 if that happens at j = 3).
- `method = "descend"`: start from the whole data set and shrink; return the first radius at which the mean is above its quantile **and** the median is above its quantile.
- The authors' rule (2026-10-23): ascend for d < 10, descend for d >= 10. The simulation scripts and, since that date, all nine real-data blocks follow it.

`simul$average[k]` and `simul$median[k]` are the lower α-quantiles of the mean and median NND of k uniform points in the unit d-ball.
So the test rejects complete spatial randomness when either statistic falls below its own critical value: two marginal level-α tests joined by "or".
The code forms no p-values and applies no Holm correction, although `alg:NNtest` in the paper describes a Holm step-down procedure.
The size of the joined test is at most 2α.
(The same discrepancy was found and corrected in the Neurocomputing paper, which uses the same function.)

The incremental nearest-neighbour update in `nnccd.radi()` dates from 2026-09-23; it returns exactly the radii of the original per-step computation, which is kept as `R/ccds/UN_CCD_radi_recompute_reference.R` (checked by `revision_experiments/tr1/92b_validate_incremental_radi.R` and `93_validate_incremental_extended.R`).

---

## 5. Quantile tables

### 5.1 Names and contents

- NN tables: `R/NN-test_quantile/NN-test-simul_<d>d_<L>%.RData`; RK tables: `R/RK-test_quantile/RK-test-simul_<d>d_<L>%.RData`.
- `<L>` encodes the level: 80 → α = 20%, 85 → 15%, 90 → 10%, 95 → 5%, 99 → 1%, 999 → 0.1%.
- Each file holds one object, `simul`.
  An NN table has `simul$average` and `simul$median`, one entry per ball size k.
  An RK table has `simul$r` (the radius grid), `simul$quan` (a list keyed by the level, e.g. `"0.99"`, each a matrix with one row per ball size) and `simul$Kest.m` (the raw Monte Carlo draws, which the detectors do not read).
- The scripts that made them are `R/NN-test_quantile/*.R`, `R/RK-test_quantile/*.R` and `R/ccds/{NN,RK}-test_quantile/*.R`.

### 5.2 The tables the clustering drivers load

| Driver group | d | Table | Level | Search direction (UN) |
|---|---|---|---|---|
| UN-CCDs, uniform and Gaussian | 2 | `NN-test-simul_2d_85%` | 15% | ascend |
| | 3 | `NN-test-simul_3d_90%` | 10% | ascend |
| | 5 | `NN-test-simul_5d_95%` | 5% | ascend |
| | 10 | `NN-test-simul_10d_99%` (Gaussian: `_999%` until 2026-10-18, see section 7) | 1% | descend |
| | 20 | `NN-test-simul_20d_999%` | 0.1% | descend |
| RK-CCDs, uniform and Gaussian | 2, 3, 5 | `RK-test-simul_<d>d_99%`, `quant = 0.99` | 1% | |
| | 10, 20 | `RK-test-simul_<d>d_999%`, `quant = 0.999` | 0.1% | |
| UN-CCDs, real data | | `2d_85%` (aggregation, asymmetric, R15, D31), `4d_90%` (iris), `5d_95%` (knowledge), `6d_95%` (wholesale), `7d_95%` (seeds), `9d_95%` (breast cancer); the table has the data's own d, and the level follows the authors' rule of 2026-10-23: 15% for d = 2, 10% for d = 3 and 4, 5% for 5 <= d <= 9. Before that, knowledge, wholesale and asymmetric loaded `6d_95%`, `8d_95%` and `2d_80%` | | ascend for all nine data sets since 2026-10-23, the rule "ascend for d < 10, descend for d >= 10" (seeds and breast cancer descended before) |
| RK-CCDs, real data | | `2d_99%` (aggregation, asymmetric, R15, D31), `4d_99%` (iris), `5d_99%` (knowledge), `6d_99%` (wholesale), `7d_99%` (seeds) (`quant = 0.99`), `9d_999%` (breast cancer, `quant = 0.999`); the table has the data's own d since 2026-10-23. Before that, knowledge and wholesale loaded `6d_99%` and `8d_99%` | | |
| KS-CCDs | | none | | |

All of these files are present.
The NN tables at d = 2, 3 and 4 cover ball sizes up to 1000 and the others up to 5000.
D31 (n = 3100) uses the 1000-entry 2d table: the search reads past the end of the table only if some ball grows beyond 1000 points, and it would then stop with an error rather than give a wrong radius.

How far the level in a file name has been checked: the Neurocomputing revision measured, against independent draws, the level of every table behind that paper's real-data results (`revision_experiments/tr1/63`–`68`).
Of the tables used here, that covers the NN tables `5d_95%`, `6d_95%`, `7d_95%`, `9d_95%`, `10d_99%` and `10d_999%`, and the RK tables `5d_99%` to `7d_99%` and `10d_999%` (`8d_99%` was measured too, but no driver loads it since 2026-10-23).
The NN tables `2d_85%`, `3d_90%`, `4d_90%` and `20d_999%`, and the RK tables `2d_99%`, `3d_99%`, `4d_99%`, `9d_999%` and `20d_999%`, have not been measured; their levels come from the file names only.

### 5.3 What is in this clone

`*.RData` is ignored by git, so none of these files comes with a clone of the GitHub repository (the only tracked tables are five NN tables at d = 166, 274 and 400).
This clone holds every table file that exists in the TR1 (Neurocomputing) and TR2 (Pattern Recognition) working copies, apart from the four mislabelled ones in section 5.4:

- `R/NN-test_quantile/`: 56 tables, 3 MB.
  d = 2 (80, 85, 90, 95%); 3 and 4 (80, 85, 90%); 5 (90, 95%); 6 to 9 (95%); 10, 12, 16, 18, 19, 20, 21 (99 and 99.9%); 11, 13, 14, 15, 17 and 22 to 28 (99%); 29 to 34, 50, 100, 166, 274, 400, 500 (99.9%); and `_spliced` versions of the d = 166 and 400 tables that cover more points.
- `R/RK-test_quantile/`: 43 tables, 17.4 GB (108 to 711 MB each).
  d = 2, 3, 4 (99%); 5 to 10 and 20 (99 and 99.9%); 11 to 19, 21 to 35, 50 and 100 (99.9%).
- `revision_experiments/results/tr2/probes/` and `results/tr1/cost_probe/`: small diagnostic tables from the other two papers (NN d = 12 and 166; RK d = 166, 274 and 1555; a cost probe at d = 5).
  The RK tables at d = 166 and 274 come from a 20-iteration run; nothing in the clustering code reads any of these.

Checked by MD5 on 2026-10-18: each of the 673 table files in the TR1 working copy (the `R/` table folders and `revision_experiments/results/`) has a byte-identical copy here.
TR2 had eleven files that TR1 lacks: the seven d = 22 to 28 tables, which were copied here, and the four mislabelled d = 6 to 9 files, which were not.
A search of all of `G:\` and `D:\` found one other set, `G:\code_working_folder\general functions\{NN,RK}-test_quantile\`, the original working folder from before the Neurocomputing repairs: its 43 RK tables and 45 of its 66 NN tables are the same files as here, and the other 21 are the mislabelled originals of section 5.4 (d = 6 to 9 `_99%`, and the 1% runs saved as `_999%` at d = 11 to 19 and 21 to 28).
Do not take tables from that folder.
The Box copy of the project was not mounted when this was checked; the only table known to be there is a copy of `RK-test-simul_2d_99%.RData`.

The RK probe tables at d = 166, 274 and 1555 are entirely zero quantiles (the K-function test carries no signal at those dimensions), which is why the d = 274 and d = 1555 files are byte-identical.
The NN probe `results/tr2/probes/NN-test-simul_12d_999%.RData` is a timing probe: its values lie above even the 1% table, so it does not hold a 0.1% quantile despite its name.
Never use it as a production table.

### 5.4 Tables with a wrong level in their name

These were found during the Neurocomputing revision (records: `revision_experiments/tr1/DELETED_DUPLICATE_TABLES.md` and `INSTALLED_REGEN999_TABLES.md`).
None of them is used by a clustering driver.

- `NN-test-simul_{6,7,8,9}d_99%.RData` were byte-identical copies of the 5% tables `..._95%.RData` (the generator scripts `R/ccds/NN-test_quantile/<d>d_9599%.R` saved one Monte Carlo run under both names).
  They are not in this clone and should not be copied back; the `_95%` files are correct.
- At d = 11 to 19 and 21 the shipped `_99%` and `_999%` files were one 1% run saved twice.
  The `_99%` files are correct; the wrong `_999%` copies were removed, and genuine 0.1% tables were built for d = 12, 16, 18, 19 and 21.
- At d = 22 to 28 the Neurocomputing cleanup deleted both names on purpose (2026-08-13), because the 1% level of these runs is inferred, from the generator scripts and the five measured pairs, and was never measured at these dimensions.
  The files `NN-test-simul_<d>d_99%.RData` were copied into this clone from the Pattern Recognition copy, and restored to the Neurocomputing repository on 2026-10-18, so that a copy survives; their level is still unmeasured.
  No clustering driver uses them.

### 5.5 Missing tables are generated on the spot

`load_quantile_table()` (`R/ccds/quantile_table.R`) loads a file that exists exactly as `load()` would.
If the file is missing it leaves a placeholder, and the first call that needs the table generates it with the original Monte Carlo estimator at the data's own size and caches it under `R/<RK|NN>-test_quantile/generated/`.
Generation of an RK table takes hours, so set `CCD_QUANTILE_GENERATE=FALSE` to make a missing table stop the run instead.
A generated table reproduces published numbers only up to Monte Carlo error.
The settings (`CCD_QUANTILE_NITER`, `CCD_QUANTILE_CORES`, `CCD_QUANTILE_SEED`, ...) are listed under "Quantile tables" in `revision_experiments/README.md`.

The other table folders in this clone (`R/NN-test_quantile_n200/`, `_d21_regen/`, `_provenance/`, `_regen999/`, `_seedstab/`) belong to the Neurocomputing revision's own checks and are not needed here.

---

## 6. Problems in the committed code, with fixes

### 6.1 Problems that stop a run

1. **UN-CCD simulation drivers call an undefined function.**
   All 40 files `simulations/clustering/UN_CCDs/Simulations/{Uniform,Gaussian}/*/` call `UNMCCD_outlier(datax = x, simul = simul, method = method)` inside `foreach`, but they source `methods/clustering/UN-CCDs.R`, which defines only `UNCCD_clustering()`.
   The call has been there since the first commit.
   Replace `UNMCCD_outlier(` by `UNCCD_clustering(`; the arguments are the same.
   With this one change, `simulations/clustering/UN_CCDs/Simulations/Uniform/2d/2d_235clusters_n50.R` reproduces `slurm-1442992.out`, and so the n = 50, d = 2 block of Table S4, to every printed digit.
2. **The real-data loader changes into a folder that is gone.**
   `data/clustering/Real_Data_Collection.R` line 5 runs `setwd(here::here("simulations/clustering/Algo_Compare_Clustering/Real_Datasets"))`; the 2026-05-10 reorganisation moved those files to `data/clustering/` and did not update this line.
   Change it to `setwd(here::here("data/clustering"))`.
   The loader leaves that folder as the working directory afterwards.
3. **Missing packages.**
   `readxl` (needed by the loader for the knowledge data), `factoextra`, `ClusterR` and `kernlab` (baselines) are not installed in R 4.6.1 on this machine.

### 6.2 Things to know before running

- Every simulation driver runs `makeCluster(detectCores())` three times and ends with `save.image("<d>d_235clusters_n<n>.RData")` in its own folder (`*.RData` is ignored by git).
- The data sets are generated on the main process with `set.seed(2024)` before the parallel step, and the clustering is deterministic, so the results do not depend on the number of cores.
- `Real_Data_KSCCD.R` sets `cores = 20` by hand.
- On this Windows machine `Rscript` must be started from PowerShell, not from Git Bash (known segfault), and it is not on `PATH`: use `"C:\Program Files\R\R-4.6.1\bin\Rscript.exe"`.
- The paper links to `.../tree/main/CCDs_Outlier_Detection`, which does not exist; the repository root is the right link.

---

## 7. Settings in the code that differ from the paper

Reproduce with the code as it is; these are for the authors to decide.

| Item | Paper / supplement | Code |
|---|---|---|
| MC-SRT combination | Holm step-down on p-values (`alg:NNtest`) | either statistic below its own α-quantile; no p-values (section 4) |
| UN-CCD level at d = 10 | 1% (S1: 15/10/5/1/0.1% for d = 2/3/5/10/20) | resolved 2026-10-18. The Gaussian scripts loaded the 0.1% table (logs slurm-1451720 to -23), but Table S8 comes from the earlier 1% run (slurm-1443096 to -99). The authors kept 1%, and the four scripts now load `NN-test-simul_10d_99%` |
| Gaussian cluster centres | "all other settings remain the same as in the uniform settings" (centres 3 apart) | `cls_dis = 6`: centres 6 apart |
| Gaussian covariance | Δ·I with Δ ~ U(0.8, 1.2) | radius²·I with radius ~ U(0.8, 1.2), so Δ ranges over (0.64, 1.44) and is not uniform |
| Repetitions | 500 everywhere | resolved 2026-10-18. `RK_CCDs/Simulations/{Uniform,Gaussian}/3d/3d_235clusters_n500.R` had `iteN = 50`, so the published cells (slurm-1442239, -1443038) were 50-repetition means. They were rerun with 500 (`revision_experiments/tr3/ITEN500_PROTOCOL.md`), Tables S3 and S7 were updated, and the scripts now set `iteN = 500` |
| Real-data level | "α is set to 1% for all datasets" (S2) | UN-CCDs: 15% (d = 2), 10% (d = 4) or 5% (d >= 5) depending on the data set, the authors' rule since 2026-10-23 (asymmetric was 20% before); RK-CCDs: 1%, and 0.1% for breast cancer |
| Real-data dimensions | knowledge d = 6, wholesale d = 8, iris n = 150 | knowledge has d = 5 and wholesale d = 6; the RK- and UN-CCD scripts ran them with the 6d and 8d tables until 2026-10-23 (they now load the 5d and 6d tables: `5d_99%` and `6d_99%` for RK-CCDs, `5d_95%` and `6d_95%` for UN-CCDs); iris has n = 147 (section 3.4) |
| Real-data scaling | synthetic sets (asymmetric, R15, D31) not normalised | asymmetric is z-scored |
| Real-data sets | seven sets | the scripts also run breast cancer and aggregation, which the paper does not report |
| DBSCAN | S2: "MinPts to 4, and find the value at the first ``elbow'' of the sorted 4-dist ... plot", Eps set to it | `minPts = k + 1` with k = 4, so minPts = 5 in R's `dbscan`, which counts the point itself (a core point has 4 other points within Eps), against the supplement's MinPts = 4. Eps: since 2026-10-23 not "the first elbow" but the point of the sorted 4-NN distances farthest below the chord, found by `knee_eps()` in the script's header (both axes scaled to [0, 1], first index on ties) and multiplied by (1 + 1e-9) so that the point defining it lies inside the radius. The supplement text is to be updated to match. Before, it was typed in per data set after looking at `kNNdistplot()` (iris 0.8, seeds 0.9, knowledge 1.5, wholesale 1, breast cancer 2, aggregation 1.2, asymmetric 0.15, R15 0.35, D31 0.6) |
| Spectral | Gaussian similarity, K-means on the embedding | resolved 2026-10-23. `kernlab::specc(as.matrix(X), centers = k)`: the Gaussian (`rbfdot`) kernel on the data, width chosen by `specc` (`kpar = "automatic"`, which subsamples 75% of the rows and runs `kmeans`, so it depends on the seed). Before, the script passed `as.matrix(dist(X))`, which `specc` clusters as a data matrix (RBF kernel on the rows of the distance matrix) |
| K-means++ on wholesale | k = 3 (the true number of clusters) | resolved 2026-10-23: `k = 3`; the script had `k = 6` |
| Baseline seeds | not stated | Since 2026-10-23 every dataset block of the five baseline scripts starts with `set.seed(1)`, and K-means++ passes `seed = 1` to `KMeans_rcpp()`, so no cell depends on the block order |
| Runtime hardware | Linux workstation, i9-13900K | the runtime table equals the run times in the SLURM job logs, and the original drivers read their files from `/mmfs1/home/...` |

---

## 8. Running

### 8.1 Setup

```powershell
cd G:\Submissions\TR3\Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages(c('readxl','factoextra','ClusterR','kernlab'))"
```

Packages used: `here`, `mclust`, `parallel`, `doParallel`, `foreach`, `MASS`, `cluster`, `igraph`, `spatstat`, `hypercube` (simulations); `R.matlab`, `foreign`, `dplyr`, `readxl` (data); `dbscan`, `ggplot2`, `factoextra`, `ClusterR`, `kernlab` (baselines).
`here::here()` finds the clone root from the `.git` folder, so a driver can be started from any working directory inside the clone.

### 8.2 One simulation cell

```powershell
$env:CCD_QUANTILE_GENERATE = "FALSE"
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "simulations/clustering/UN_CCDs/Simulations/Uniform/2d/2d_235clusters_n50.R"
```

(after fixing `UNMCCD_outlier`, section 6.1).
It prints three lines (2, 3 and 5 clusters) with mean ARI, mean Sil and the SR in percent, then the run time.
Compare them with the `slurm-*.out` files in the same folder and with the table in the supplement.

Run times in the original logs (500 repetitions, all cores), for example: UN-CCDs 15 seconds at d = 2, n = 50 and 13 minutes at d = 2, n = 500; RK-CCDs 5 minutes and 1.3 days for the same two cells; KS-CCDs up to 16.5 hours (d = 5, n = 500), because of the 200-value density search.
This machine is faster: the d = 2, n = 50 cells took 4 seconds (UN-CCDs, 8 cores) and 1.4 minutes (RK-CCDs, 12 cores).
Start with n = 50 and n = 100.

### 8.3 Real data

After fixing the loader (section 6.1), run each script in `simulations/clustering/*/Simulations/Real_Data/` and `simulations/clustering/Algo_Compare_Clustering/*/Simulation/Real_Data/`.
Each prints one line per data set: ARI, average silhouette and number of clusters.

### 8.4 The noise experiment

It has to be written.
The supplement (S1.3) gives the settings: three Gaussian clusters with centres (3, ..., 3), (9, 3, ..., 3), (3, 9, 3, ..., 3), covariance Δ·I with Δ ~ U(0.8, 1.2), 200 cluster points plus 3%, 5%, 10% or 20% noise "randomly distributed within the scope of the entire dataset", d = 3 and 10, 1000 repetitions, UN-CCDs only.
`Gaussian_cls(..., noise = TRUE, r_noise, n_noise)` generates such data, but the paper does not give `r_noise`, the level α, or the search direction used; ask the authors before treating any rerun as a reproduction.
Note that the noise branch of both generators adds the mean centre with `+ center_noise` to an `n_noise` × d matrix; R recycles that vector down the columns, so unless all coordinates of the mean centre are equal the noise cloud is not centred where intended.
`sweep(noise_data, 2, center_noise, "+")` does what was meant.

---

## 9. What was checked for this guide

All checks were run on 2026-10-18 with R 4.6.1 on Windows, from this clone, without writing into it (drivers were patched in memory, with `setwd()` and `save.image()` removed, and `CCD_QUANTILE_GENERATE=FALSE`).

| Check | Result |
|---|---|
| UN-CCDs, uniform, d = 2, n = 50, `UNMCCD_outlier` → `UNCCD_clustering`, 8 cores | identical to `slurm-1442992.out` in every printed digit (ARI 0.98000252747998 / 0.926033622988254 / 0.610065166255842; SR 99.6 / 90.2 / 21.4%); 4 s |
| RK-CCDs, uniform, d = 2, n = 50, driver unchanged, 12 cores | identical to `slurm-1432233.out` (ARI 0.981634855531438 / 0.906270958459637 / 0.517320463966977; SR 99.8 / 84.4 / 8.2%); 1.4 min |
| Real data, loader with the folder fixed | iris UN-CCDs 0.599 / 0.445 / 3, RK-CCDs 0.557 / 0.393 / 3; seeds UN-CCDs 0.509 / 0.467 / 2 (descending; since 2026-10-23 the seeds block ascends and gives 0.239 / 0.286 / 2, and the paper was updated); R15 UN-CCDs 0.975 / 0.742 / 15 (ARI / Sil / clusters): all equal to the paper |
| `methods/clustering/UN-CCDs.R` | defines `UNCCD_clustering` only; `UNMCCD_outlier` does not exist after sourcing it |
| The 22 tables in section 5.2 | all load; each holds one `simul` object with the fields described in 5.1 |
| Table copies | MD5-identical to the TR1 and TR2 working copies (section 5.3) |
| RK-CCDs d = 3, n = 500 | `iteN = 50` in the original runs, confirmed by the script and the run times. Rerun with 500 on 2026-10-18: replicates 1-50 of the 2-cluster data reproduce both logs to 15 digits |
| Real data, all nine data sets and all eight methods (2026-10-19) | CCDs: every paper cell except UN-CCDs on D31 (0.749 / 0.444 / 31 against 0.799 / 0.474 / 31, cause open). Baselines: several paper cells differ from the code. Details in `revision_experiments/tr3/README.md`, "Real-data rerun". The three UN-CCD cells whose table changed on 2026-10-23 (knowledge 0.147 / 0.165 / 2, wholesale -0.006 / 0.112 / 2, asymmetric 0.883 / 0.604 / 5) differ from the paper's table and were rerun; they equal the original script blocks and the α sweep to 15 digits. The seeds and breast cancer UN-CCD cells were rerun ascending on 2026-10-23 (seeds 0.239 / 0.286 / 2, breast cancer 0.049 / 0.034 / 34; before: 0.509 / 0.467 / 2 and 0.075 / 0.208 / 2); they equal the script blocks and the sweep's config S, level 5, to 15 digits |

Not checked: the noise experiment (no driver). The 2026-10-18 audit reran n = 50 in all 30 simulation folders; the larger cells were compared with their logs, not rerun.
The exact matches show that the code changes made since the original runs (the `here::here()` paths, the on-the-spot table loader, and the incremental radius search) did not change these results.
