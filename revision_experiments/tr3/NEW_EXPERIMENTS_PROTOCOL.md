# Protocol: new experiments for the KBS rebuttal (2026-10-25)

Written 2026-10-25, before any full run of these experiments. The authors asked for
experiments 1, 2, 3 and 5 of the revision plan on 2026-10-25; experiments 4 (failure
cases, extended noise) and 6 (sensitivity to N) are not run. Each script was written by
one agent and checked by a separate verifier agent before its full run.

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

Filled in after the runs.
