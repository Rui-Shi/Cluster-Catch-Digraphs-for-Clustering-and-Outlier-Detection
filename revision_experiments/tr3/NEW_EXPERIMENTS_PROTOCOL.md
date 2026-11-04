# Protocol: new experiments for the KBS rebuttal (2026-10-25)

Written 2026-10-25, before any full run of these experiments. The authors asked for
experiments 1, 2, 3 and 5 of the revision plan on 2026-10-25; experiments 4 (failure
cases, extended noise) and 6 (sensitivity to N) are not run. Each script was written by
one agent and checked by a separate verifier agent before its full run.

Note of 2026-11-03: the console logs this file cites (`_logs/...`) were removed from the clone
in the clean-up before the merge to `main`; see `README.md`, "Results".

## Decisions taken before the runs (authors, 2026-10-25)

- UN-CCDs level for real data with d >= 10: alpha = 1% for 10 <= d < 20 (tables
  `NN-test-simul_<d>d_99%.RData`), 0.1% for d >= 20 (`..._999%`), as in the d = 10 and
  d = 20 simulations; search direction descend (d >= 10).
- RK-CCDs level for real data with d >= 10: 0.1% (`RK-test-simul_<d>d_999%.RData`, quant
  0.999), the rule of the RK-CCD simulations.
- KS-CCDs: the delta grid of the other real data sets.
- Ablation: only the search direction, at d >= 10, on a medium setting.
- All six candidate data sets are run (wine, vehicle, image segmentation, breast cancer
  Wisconsin diagnostic, hepatitis, lymphography).

## Experiment 1: real data with d >= 10 (`30_hd_load.R`, `31_hd_real.R`, `32_hd_baselines.R`)

Data in `revision_experiments/tr3_data/` (UCI downloads of 2026-10-25, md5 in
`MANIFEST.csv`) and `data/clustering/` (hepatitis, lymphography). Preprocessing as
`Real_Data_Collection.R`, plus one generic rule: zero-variance columns are dropped before
`scale()` (image segmentation's REGION-PIXEL-COUNT is constant over all 2310 rows).
Facts found before the runs: the vehicle files hold 846 rows (the UCI page says 946);
lymphography has 4 classes (81/61/4/2), not the 8 of `data/clustering/Datasets_Summary.txt`;
hepatitis has 80 complete rows of 155; hepatitis (12 of 19 features binary) and
lymphography (all 18 features coded with 2 to 8 values) give lattice distances, which the
NND test, built for continuous positions, does not model.
All tables needed cover n <= 5000 (checked 2026-10-25).
Methods: UN-, RK-, KS-CCDs (settings in `settings_realdata_hd.csv`, filled by the rules
above from the final d) and DBSCAN, MST, K-means++, spectral clustering and Louvain with
the setting rules of Supplement Section S2. Metrics as `11_real.R`.

## Experiment 2: Type I error (`33_type1.R`)

- Part test: size of one MC-SRT under CSR, d in {2, 3, 5, 10, 20} (simulation tables),
  m in {5, 10, 20, 50, 100, 200, 500} members, B = 10000. Two nulls: the table generator's
  (m points iid in the unit ball) and the exact null of the clustering's ball (m - 1 points
  in the ball and the farthest member on its boundary). Rates for mean only, median only and
  the rule (either).
- Part ccd: UN-CCDs (published settings) on n uniform points in [0,1]^d without clusters,
  d in {2, 3, 5, 10, 20}, n in {50, 100, 200, 500}, 500 replicates; proportion with
  k_hat = 1.

## Experiment 3: search-direction ablation (`34_ablation_direction.R`)

UN-CCDs at d = 10 and 20, n = 200 ("a medium simulated data set"), uniform and Gaussian,
2/3/5 clusters, 500 replicates on the paper's datasets: the ascending search against the
published descending runs in `sim_replicates.csv` (paired by replicate and data
fingerprint); Wilcoxon signed-rank on ARI, Holm over the 12 cells.

## Experiment 5: paired tests (`35_paired_tests.R`)

UN-CCDs against RK- and KS-CCDs at n = 50 (the only n with per-replicate rows for all three
methods), 30 cells, Wilcoxon signed-rank on per-replicate ARI, Holm over 60 tests.
This one reads existing results only; it was verified and run before this protocol was written
(commit 8696b5c).

## Outcome

Experiment 5 (run 2026-10-25, log `_logs/paired_tests_20261025.log`): 60 tests; Holm 0.05:
UN vs RK 9 higher, 4 lower, 17 not significant (for d >= 10: 8 / 1 / 3); UN vs KS 3 higher
(Gaussian, d = 20), 19 lower (15 with d <= 5), 8 not significant. The 90 cell means equal the
supplement's n = 50 values. Supplement Table S14.

Experiment 5 at n = 200 (author request 2026-11-01; README, "Paired tests at n = 200 (2026-11-01/02)"): after the
rerun of the 26 n = 200 folders without per-replicate rows, the same tests at n = 200, Holm within each n (60 tests
per n; log `_logs/paired_tests_20261102.log`). n = 200: UN vs RK 13 higher, 9 lower (all d <= 5), 8 not significant
(d >= 10: 12 of 12 higher); UN vs KS 2 higher (Gaussian, d = 20), 18 lower (14 with d <= 5), 10 not significant.
The n = 50 results above are unchanged; the 180 cell means equal the supplement. Table `tab:paired_Apdx`.

Experiment 3 (run 2026-10-25, 1.7 min, log `_logs/ablation_direction_20261025.log`): 12 cells x
500 replicates; the descending reruns of replicates 1-3 equal the stored rows in all 12 cells.
Mean ARI difference (ascending minus descending) from -0.012 to 0.0045; descending higher in
10 of 12; raw p < 0.05 in 5 cells, Holm-adjusted p >= 0.192 in all. The seconds columns are
not comparable (the descending rows were timed in another run). Supplement Table S15.

Experiment 2 (run 2026-10-25, 15.2 min, log `_logs/type1_20261025.log`; part check: the table
index replays the clustering's radii for 120 of 120 centres at d = 2 and d = 10). Size of one
test divided by alpha (B = 10000): table null 1.30-1.56 for d <= 5, up to 3.50 at d = 20;
radius-search null 0.88-1.48 for d <= 5, up to 2.70 at d = 20 (at most 35 rejections at
d = 20). UN-CCDs on uniform data without clusters: P(k_hat = 1) 0.006-0.518 at n = 50, at most
0.090 for n >= 200, where k_hat = 2 in 93.8% of replicates; mean silhouette at most 0.219
(0.065 for d >= 10). The silhouette refinement keeps any k >= 2 with positive silhouette, so
UN-CCDs do not test for the presence of clusters. Supplement Tables S16 and S17.

Experiment 1 (run 2026-10-25 23:26 to 2026-10-26 00:39, one method at a time; logs
`_logs/real_hd_<method>_20261025.log` and `_logs/real_hd_baselines_<method>_20261026.log`;
exit codes 0 in `real_hd_exitcodes_20261025.csv` and `real_hd_baselines_exitcodes_20261026.csv`):
UN-CCDs 11.1 min, RK-CCDs 1.7, KS-CCDs 58.0 (segmentation 54.5), the five baselines 1.8; 18 CCD
and 30 baseline cells, all in the done files. Final n/d/k: hepatitis 80/19/2, lymphography
148/18/4, wine 178/13/3, wdbc 569/30/2, vehicle 846/18/4, segmentation 2086/18/7 (224 duplicate
rows removed). ARI of UN-CCDs: 0.112, -0.004, 0.342, 0.168, 0.069, 0.094 (same order); higher than
RK-CCDs on wdbc, vehicle and segmentation, lower on the other three; higher than KS-CCDs on three.
K-means++ has a higher ARI than UN-CCDs on all six, spectral clustering on five (not on wdbc); the
best ARI is K-means++ on wine (0.915), lymphography (0.169) and wdbc (0.677), spectral clustering
on hepatitis (0.430), vehicle (0.101) and segmentation (0.511). k_hat = k: KS-CCDs 3 (hepatitis,
wine, wdbc), UN-CCDs 1 (wdbc), MST 1, Louvain 1, RK-CCDs 0, DBSCAN 0. RK-CCDs give 154 clusters on
wdbc; DBSCAN one cluster plus noise on five sets and 11 on segmentation; MST one cluster on wine
(Sil undefined); spectral clustering and KS-CCDs have equal ARI and Sil on wdbc to 17 digits
(labels are not stored, so the partitions were not compared). UN-CCDs took 616.72 s on
segmentation (one R process). Table 6 of the manuscript (`tex/new_realhd.tex`; manuscript
commits 21ea0e9, d420d64).
Correction to the facts written before the runs: hepatitis has 13 binary features of 19
(attributes 3 to 14 and 20, counted by script), not 12. The lymphography count (18 features with 2
to 8 values) was checked the same way.

## Experiment 7: baselines on the simulated data (`37_sim_baselines.R`)

Written 2026-10-31, before the driver. The authors approved it on 2026-10-31 (manuscript
`TODO.md`). Until now the five baselines ran on real data only; this run puts them on the
simulated data and the noise data of the CCDs.

Data
- Simulations: the 40 cells of `settings_simulation.csv` (uniform and Gaussian; d = 2, 3, 5,
  10, 20; n = 50, 100, 200, 500), each with 2, 3 and 5 clusters and 500 replicates: 120
  settings, 60000 data sets. Generated by `tr3_sim_data()` from the UN-CCDs row; the seed and
  data columns of the three method rows are equal (checked 2026-10-31). Raw coordinates, as the
  CCDs used.
- Fingerprints (`tr3_fp()`) are compared with the published rows of `sim_replicates.csv` where
  they exist: n = 50 at every d, and d = 10, 20 at every n (66 settings). A mismatch stops the
  cell. The other 54 settings (d <= 5, n >= 100) have no stored fingerprint; they come from the
  same generator and seed.
- Noise: d = 3, 10; l = 3%, 5%, 10%, 20%; 1000 replicates per level, generated as `23_noise.R`
  (published variant: r_noise 7, sweep, the seed and cluster settings of the Gaussian n = 200
  row, all four levels in order from one seed). Every fingerprint is compared with
  `noise_replicates.csv`.
- All data are generated on the main process. The workers load only the baselines' packages
  (`Uni-Gau_cls.R` masks functions of igraph, cluster and MASS).

Methods: the rules of Supplement Section S2, the same on every data set, no tuning per setting
- DBSCAN: minPts 5, Eps = `knee_eps(X, 4)`. MST: `MST_Clustering(X, 2)`. K-means++:
  `KMeans_rcpp(X, k, initializer = "kmeans++", seed = 1)`. Spectral: `specc(X, centers = k)`,
  automatic kernel width. Louvain: Gaussian-kernel graph with sigma = 2, `cluster_louvain()` at
  resolution 1. `set.seed(1)` before each method on each data set, as each real-data block does.
  K-means++ and spectral clustering get the true k (2, 3 or 5; noise data: 3).
- The driver calls the same functions as the real-data blocks, through their package
  namespaces. Check before any full run: the driver's code reproduces the stored real-data cells
  of `real_baselines_fixed.csv` (9 sets) and `real_hd_baselines.csv` (6 sets): ari, sil, k_hat,
  DBSCAN eps and spectral sigma equal as fmt17 text, with their `_params` files. The full run
  starts only if all are equal.
- DBSCAN k_hat: the number of distinct labels other than 0 (noise). The real-data block computes
  `length(unique(label)) - 1`, which is one too low when no point is noise. Both counts are
  written. The check also gives the number of noise points on each of the 15 real data sets; a
  set with none means its published k_hat is one too low (to be reported to the authors).
- KS-CCDs' delta is chosen per setting (largest mean silhouette); the baselines use fixed rules.
  Say so wherever they are compared.

Metrics per replicate
- ARI against the true labels (DBSCAN noise counts as one group, as in the real-data tables);
  Sil = `mean(silhouette(label, dist(X))[, 3])`, NA for 1 or n clusters; k_hat; seconds.
- Noise data (default chosen here; UN-CCDs pick k by the silhouette on the cluster points, the
  baselines have no such step): each method runs on all 200 + n_noise points; ARI and Sil are
  computed on the 200 cluster points, as for UN-CCDs; k_hat counts the clusters over all points.
- An error of a method on one replicate is written as NA with its message and counted; it is not
  rerun with other settings.

Summaries and comparisons (fixed now)
- Per setting: mean ARI, mean Sil (over defined values, with the NA count) and SR, the share of
  replicates with k_hat = k. No SR for K-means++ and spectral clustering (k is given).
- Main table: for each distribution and d, the mean over its 12 settings (4 n x 3 k) of mean ARI
  and SR, eight methods. CCD values are read from Supplement Tables S2-S4 and S6-S8 in
  `Supplementary_Material.tex`, the only complete source (`36_new_tables.R` already reads it).
- Wording: "comparable" when two such means differ by at most 0.02; otherwise "higher" or
  "lower". With it, the number of settings (of 12) where UN-CCDs have the higher or lower mean ARI.
- K-means++ and spectral clustering are the reference given k; DBSCAN, MST and Louvain are the
  direct competitors (they estimate k).
- Noise: mean ARI, Sil and SR for each d and l, the five baselines next to UN-CCDs (RK- and
  KS-CCDs were not run on the noise data).
- Not planned (authors' decision): paired tests of UN-CCDs against the baselines. The
  per-replicate rows allow them in the 66 settings with UN-CCD rows and in all noise cells.

Run
- The driver is written by one agent and checked by an independent verifier before any full run.
- Smoke test: 3 replicates per cell, every method; it gives the time per method. The full run
  starts if the projected total is under 12 hours; otherwise the authors decide.
- One Rscript per method, one after the other (never two runs on G: at once), from PowerShell;
  PSOCK workers over the replicates; each cell (method x setting, 500 or 1000 rows) is appended
  when it finishes, then marked in the done file.

Output (`results/tr3/`): `simbase_replicates.csv` and `noisebase_replicates.csv` (one row per
replicate: keys, rep, ari, sil, k_hat, DBSCAN's block count, seconds, error, data_fp), their
`_done.csv` files, and `simbase_check.csv` (the real-data check). Tables by `36_new_tables.R`.

Outcome of experiment 7 (driver 75110bf; run 2026-10-31, `all sim` 19:28 to 20:05 with 24 workers,
`all noise` to 20:07; logs `_logs/simbase_sim_20261031.log`, `_logs/simbase_noise_20261031.log`,
exit codes 0 in `_logs/simbase_exitcodes_20261031.csv`; tables `_logs/new_tables_20261031.log`):
- Gate: the real-data check reproduced 510 of 510 stored values (writer and verifier separately).
  Every one of the 15 real sets has DBSCAN noise points, so no published DBSCAN k_hat is too low.
- 600 sim cells x 500 and 40 noise cells x 1000, all in the done files. Fingerprints equal in the
  66 settings with stored ones (330 cells) and in all 40 noise cells; one fingerprint per setting
  and replicate across the five baselines. One error in 340000 replicates: spectral clustering,
  uniform, d = 2, n = 50, k = 2, replicate 4 ("NA/NaN/Inf in foreign function call (arg 1)", as in
  the original block). Sil undefined (one cluster) in 11209 MST and 747 DBSCAN replicates of 60000.
  DBSCAN found no noise point in 8827 of 60000 data sets (there the block's count is one too low).
- Mean ARI over the 12 settings of each distribution and d, by the 0.02 rule (10 groups):
  UN-CCDs against DBSCAN higher in 9 (lower at uniform d = 20, by 0.023), against MST higher in
  10, against Louvain higher in 6 (uniform, all d; Gaussian d = 20) and lower in 4 (Gaussian
  d <= 10, by 0.039 to 0.080); against K-means++ (given k) higher in 4 (uniform d >= 3), comparable
  in 2 (Gaussian d = 3, 5), lower in 4 (uniform d = 2; Gaussian d = 2, 10, 20); against spectral
  clustering (given k) higher in 4 (uniform d <= 10), comparable in 2 (Gaussian d = 3, 5), lower in
  4 (uniform d = 20; Gaussian d = 2, 10, 20). KS-CCDs (delta per setting): higher than DBSCAN and
  MST in all 120 settings; higher than K-means++ and spectral clustering in 8 groups, lower at
  Gaussian d = 10, 20; against Louvain higher in 6, comparable in 2, lower in 2.
- Noise (ARI and Sil on the 200 cluster points): UN-CCDs SR 0.955 to 0.998; DBSCAN SR at most
  0.667 (d = 3) and 0.046 (d = 10); MST SR at most 0.284; Louvain ARI at least 0.988 but SR 0.092
  to 0.563 at d = 3 and 0 at d = 10 (the noise points form extra clusters); K-means++ and spectral
  clustering (given k = 3) ARI 0.881 to 0.991. The UN-CCD values equal the supplement's noise table.
- Independent recomputation in Python from the replicate files (scratch `check_simbase.py`): the
  80 rows of `simbase_by_d.csv`, the 140 of `simbase_comparisons.csv` and the 48 of
  `noisebase_summary.csv` agree.

## Experiment 8: UN-CCDs on the noise data with k chosen on all points (`23_noise.R allpoints`)

Written 2026-10-31, before the code change. The authors chose this on 2026-10-31, after a
cross-check of the manuscript found that, on the noise data, UN-CCDs choose k by the silhouette
on the 200 cluster points (`nnccd.silhouette(ind = 1:200)`, which needs to know which points are
noise), while the five baselines of experiment 7 count k over all points. Neither document said
the rules differ. `noise_replicates.csv` already has k chosen on all points (`k_all`), but not the
ARI or Sil of that choice.

Data: as the published variant of `23_noise.R` (d = 3, 10; l = 3%, 5%, 10%, 20%; 1000 replicates;
r_noise 7, sweep; same seed and settings). Nothing in `R/ccds/` changes.

Per replicate, one UN-CCD graph and two selections of the number of clusters from it
- published: as before (silhouette on the 200 cluster points): ari, sil, k_hat;
- all points: the silhouette over all 200 + n_noise points, as `UNCCD_clustering()` does
  (`ind = NULL`; already stored as `sil_all`, `k_all`). New columns, scored as the baselines:
  `ari_all` = ARI of the labels of the 200 cluster points against the truth; `sil200_all` = mean
  silhouette of those labels on those points, NA unless they form 2 to 199 groups; `nlab_all` = the
  number of distinct labels over all points (the baselines' k_hat). SR = share with `nlab_all` = 3.
  Any replicate with `nlab_all` different from `k_all` is counted and reported.

Gate (before any number is used): in every replicate, data_fp, ari, sil, k_hat, sil_all and k_all
equal the stored published row as fmt17 text. Any difference stops the cell.

Run: smoke test (3 replicates per cell), then d = 3 and d = 10 one after the other (never two R
runs on G: at once), from PowerShell. Output, separate from the published files:
`noise_allpoints_replicates.csv` and `noise_allpoints_done.csv`.

Use in the paper (fixed now, whatever the numbers): main Table 4 shows UN-CCDs under the all-points
rule, the rule of the five baselines; Supplement Table S9 shows UN-CCDs under both rules; both
documents state both rules. The abstract, highlights and conclusion take their noise claims from
the all-points numbers.

Outcome of experiment 8 (driver fde65e5; run 2026-10-31 21:53 to 22:00 with 24 workers, exit codes 0;
log `_logs/noise_allpoints_20261031.log`; tables `_logs/new_tables_20261031b.log`):
- Gate: in all 8000 replicates, data_fp, ari, sil, k_hat, sil_all and k_all equal the published
  rows; `nlab_all` equals `k_all` in every replicate. Independent recomputation in Python (scratch
  `check_exp8.py`) gives the same summaries.
- UN-CCDs, k chosen on all points (ARI, Sil on the cluster points; SR), l = 3%, 5%, 10%, 20%:
  d = 3: ARI 0.973, 0.974, 0.975, 0.963; SR 0.998, 0.997, 0.996, 0.967.
  d = 10: ARI 0.951, 0.946, 0.939, 0.928; SR 0.967, 0.951, 0.935, 0.904.
  Published choice (cluster points): d = 3 SR 0.998, 0.998, 0.998, 0.985; d = 10 SR 0.971, 0.955,
  0.958, 0.955. The misses under the all-points choice are mostly too few clusters (d = 10, l = 20%:
  94 below 3, 2 above).
- Against the baselines that estimate k (same rule): DBSCAN SR at most 0.667, MST at most 0.284,
  Louvain 0.092 to 0.563 at d = 3 and 0 at d = 10 (ARI at least 0.988). UN-CCDs keep ARI above 0.92
  and SR above 0.90 at every level; no baseline that estimates k reaches SR 0.67.
