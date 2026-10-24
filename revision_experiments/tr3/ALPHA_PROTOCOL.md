# Protocol: one significance level for UN-CCDs on all real data sets

Written 2026-10-20, before any run. The authors asked for it on 2026-10-20:
"rerun all the UN-CCDs with 1%, 5%, 10%, 15% and 20%, check which universal
alpha produce the best results".

## Question

Supplement S2 says "α is set to 1% for all datasets for simplicity". The
UN-CCD real-data script instead loads a table per data set, at levels from
5% to 20%:

| data set | n x d (loaded) | true k | script table | α | direction |
|---|---|---|---|---|---|
| iris | 147 x 4 | 3 | `4d_90%` | 10% | ascend |
| seeds | 210 x 7 | 3 | `7d_95%` | 5% | descend |
| user_knowledge | 258 x 5 | 4 | `6d_95%` | 5% | ascend |
| wholesale | 440 x 6 | 3 | `8d_95%` | 5% | ascend |
| breast_cancer | 449 x 9 | 2 | `9d_95%` | 5% | descend |
| aggregation | 788 x 2 | 7 | `2d_85%` | 15% | ascend |
| asymmetric | 1000 x 2 | 5 | `2d_80%` | 20% | ascend |
| R15 | 600 x 2 | 15 | `2d_85%` | 15% | ascend |
| D31 | 3100 x 2 | 31 | `2d_85%` | 15% | ascend |

For user_knowledge and wholesale the script's table is for a different d
than the data has (6 and 8 against 5 and 6). This run asks which single α,
used for every data set, gives the best results.

## Levels and configurations

Levels: α = 20%, 15%, 10%, 5%, 1% (tables `80%`, `85%`, `90%`, `95%`, `99%`).
The search direction stays the script's per data set, so the answer is a
universal α, not a universal setting. Everything else is as in `11_real.R`:
the loader, `UNCCD_clustering(X, simul, method = direction)`,
ARI = `adjustedRandIndex(label, Y)`, Sil = `mean(silhouette(label, dist(X))[, 3])`,
k̂ = `length(unique(label))`.

| config | table d | table source | data sets |
|---|---|---|---|
| A (primary) | d of the loaded data | the shipped table if it exists, else generated | all 9 |
| P | the script's d (user_knowledge 6, wholesale 8) | as A | user_knowledge, wholesale (the other 7 are A) |
| G | as A | generated at every level | all 9 |
| S | as A, with ascend for seeds and breast_cancer | as A | seeds, breast_cancer |

Shipped tables at the needed d: d = 2 at 80-95%; d = 4 at 80-90%; d = 5 at
90 and 95%; d = 6 to 9 at 95%. The others are generated.

## Generated tables

One Monte Carlo run per (d, extent) gives all five levels:
(2, 1000), (4, 147), (5, 258), (6, 440), (7, 210), (8, 440), (9, 449).
The extent is the largest n among the data sets that use that d; d = 2 has
1000 entries, as the shipped d = 2 tables do.
- Per-iteration body: `.ccd_qt_nn_once(n, d)` from `R/ccds/quantile_table.R`.
  It is the transcription of `NNDestP.simpois.lower.quant()`. The draws come
  from `.ccd_qt_draws()` with niter = 10000, as the shipped NN tables have.
- Seed 20260927 + d, with **22 workers fixed**. `.ccd_qt_draws()` uses
  `clusterSetRNGStream()` and `parLapply()`, so the draws depend on the
  number of workers; another worker count gives other tables.
- Reduction as in `ccd_generate_quantile_table()`: lower quantile at
  1 − level of the mean and the median NND, per number of points.
- Saved in `results/tr3/alpha_tables/` (not under `R/`), with the settings in
  an attribute. In memory the run strips the `ccd_table` and `generated`
  attributes. The tables are then used exactly as shipped tables are, and are
  never grown.
- Check: wherever a shipped table exists at the same d and level, report the
  largest relative difference between generated and shipped entries
  (`alpha_table_check.csv`).

## Reproduction check (before any analysis)

At the script's own level, config P must equal the `published` rows of
`results/tr3/real.csv` to 15 significant digits (ARI, Sil, k̂) in all 9 data
sets. If any cell differs, stop and report; do not analyse.

## Failures

A cell that stops with an error, for example a covering ball past the end of a
table, is recorded with its error message. It is not dropped silently.

If a cell fails, the primary ranking uses only the data sets where all five
levels ran. The failed cells are listed next to it.

## Criteria (fixed before the run)

- Primary: mean ARI over the seven data sets of the paper's table (iris,
  seeds, user_knowledge, wholesale, asymmetric, R15, D31), config A. The
  universal α is the level with the highest mean; levels within 0.001 of it
  are reported as tied.
- Secondary, per level:
  - mean Sil over the same seven data sets; k̂ = 1 counts as Sil 0 and is
    flagged;
  - the number of data sets with k̂ equal to the true k;
  - the mean rank of the level across the seven data sets, ranking by ARI
    with rank 1 for the best and ties averaged;
  - mean ARI over all nine data sets.
- Reference: the script's per-data-set levels (config P at the script's
  levels, that is, the published rerun) scored by the same criteria.
- Sensitivity: whether the primary choice changes under P, G or S.

## Caveats

- α is chosen with the labels of the same data sets it is scored on. The
  result describes these data sets; it is not an out-of-sample tuning rule.
- The generated tables carry Monte Carlo error. Config G against A, at the
  levels where both exist, shows how much that error moves the results.

## Files

| file | contents |
|---|---|
| `14_alpha_tables.R` | generates the seven tables |
| `15_real_alpha.R` | runs every (config, level, data set) cell |
| `16_alpha_summary.R` | reproduction check, then the criteria |
| `results/tr3/alpha_tables/` | generated tables |
| `results/tr3/alpha_table_check.csv` | generated against shipped |
| `results/tr3/real_alpha.csv` | one row per cell: config, level, dataset, table, source, n, d, ari, sil, k_hat, seconds, error |
| `results/tr3/alpha_summary.csv` | the criteria per config and level |

## Outcome

Run 2026-10-20 with `14_alpha_tables.R`, `15_real_alpha.R` and `16_alpha_summary.R`,
started once by `run_alpha.ps1`. All three steps exited 0. The criteria above were
not changed. The tables below are from `alpha_summary.csv` and `real_alpha.csv`.

### Readings of points the protocol left open

- Level column: α in percent (20, 15, 10, 5, 1) for tables 80, 85, 90, 95, 99.
- P and S are scored over all data sets: their own rows where they have them, config
  A's rows for the rest. P differs from A in user_knowledge and wholesale only, S in
  seeds and breast_cancer only.
- k̂ correct: k̂ equal to the true k of the protocol table, over the seven paper data sets.
- Mean rank: per data set, ranks of the five levels by ARI (1 = best, ties averaged),
  then the mean over the seven data sets. The row "script" has no rank among the five.
  `alpha_summary.csv` also has `mean_rank_with_ref`, which ranks the five levels and
  the script row together (6 candidates).
- No cell had k̂ = 1, so no Sil was set to 0.
- Reference ("script"): iris 10%; seeds, user_knowledge, wholesale and breast_cancer 5%;
  aggregation, R15 and D31 15%; asymmetric 20%; with the script's tables (config P for
  user_knowledge and wholesale, A for the rest). It is the same in all four configs.

### Reproduction check

Passed: 111 of 111 comparisons equal at 15 significant digits
(`alpha_repro_check.csv`).
- 27 published-vs-rerun comparisons: ARI, Sil and k̂ of all 9 data sets at the script's
  own level (config P for user_knowledge and wholesale, A for the other seven),
  against the `published` UN-CCD rows of `real.csv`.
- 84 G-equals-A comparisons: in the 21 cells where config A used a generated table, config G
  used the same file and gave the identical ARI, Sil, k̂ and error field.

### Tables (`results/tr3/alpha_tables/`, `alpha_table_check.csv`)

Seven tables, five levels each: 35 files. niter = 10000, 22 workers, seed 20260927 + d.
Draw times: d = 4 35 s, 7 82 s, 5 133 s, 6 615 s, 8 632 s, 9 685 s, 2 7483 s (2 h 05 min).
Step 14 took 161.2 min in all.

Largest relative difference between generated and shipped entries 2..min(extent, shipped
length), per (d, level):

| d | level | entries | average | median |
|---|---|---|---|---|
| 4 | 80% | 146 | 0.0055 | 0.0070 |
| 4 | 85% | 146 | 0.0084 | 0.0084 |
| 4 | 90% | 146 | 0.0095 | 0.0104 |
| 5 | 90% | 257 | 0.0146 | 0.0158 |
| 5 | 95% | 257 | 0.0584 | 0.0584 |
| 6 | 95% | 439 | 0.0214 | 0.0301 |
| 7 | 95% | 209 | 0.0118 | 0.0183 |
| 8 | 95% | 439 | 0.0177 | 0.0197 |
| 9 | 95% | 448 | 0.0236 | 0.0232 |
| 2 | 80% | 999 | 0.0111 | 0.0186 |
| 2 | 85% | 999 | 0.0133 | 0.0186 |
| 2 | 90% | 999 | 0.0201 | 0.0201 |
| 2 | 95% | 999 | 0.0337 | 0.0195 |

The largest is 0.058 (d = 5, 95%, average and median); every other one is at most 0.034.

### Cells and failures

110 cells (A 45, P 10, G 45, S 10). No cell failed: 0 error messages, and no covering
ball went past the end of a table. A used the shipped table in 24 cells and a generated
one in 21. Config G against A in the 24 cells where A used a shipped table: 18 are
identical in ARI and k̂ and 6 differ. The differences come from the Monte Carlo error of
the generated tables:

| cell | ARI, A (shipped) | ARI, G | k̂ A / G |
|---|---|---|---|
| D31 20% | 0.761 | 0.773 | 28 / 29 |
| D31 15% | 0.749 | 0.759 | 31 / 31 |
| D31 10% | 0.634 | 0.683 | 23 / 26 |
| D31 5% | 0.520 | 0.550 | 16 / 18 |
| wholesale 5% | -0.006 | -0.006 | 2 / 2 |
| breast_cancer 5% | 0.075 | 0.215 | 2 / 2 |

### Criteria: seven paper data sets (mean ARI, mean Sil, number of k̂ = true k out of 7, mean rank), and mean ARI over all nine

Config A (primary):

| α | mean ARI | mean Sil | k̂ correct | mean rank | mean ARI, 9 sets |
|---|---|---|---|---|---|
| 20% | 0.5226 | 0.3959 | 2 | 2.43 | 0.4803 |
| 15% | 0.4994 | 0.3793 | 3 | 2.50 | 0.4440 |
| 10% | 0.4943 | 0.3985 | 2 | 3.36 | 0.4465 |
| 5% | 0.4915 | 0.4234 | 2 | 3.14 | 0.4170 |
| 1% | 0.4833 | 0.4437 | 2 | 3.57 | 0.4486 |
| script | 0.5232 | 0.4251 | 4 | | 0.4606 |

Config P:

| α | mean ARI | mean Sil | k̂ correct | mean rank | mean ARI, 9 sets |
|---|---|---|---|---|---|
| 20% | 0.4995 | 0.4374 | 2 | 3.43 | 0.4623 |
| 15% | 0.5095 | 0.4234 | 3 | 2.43 | 0.4519 |
| 10% | 0.5056 | 0.3972 | 2 | 3.29 | 0.4553 |
| 5% | 0.4848 | 0.4227 | 3 | 3.00 | 0.4118 |
| 1% | 0.4901 | 0.4041 | 2 | 2.86 | 0.4539 |
| script | 0.5232 | 0.4251 | 4 | | 0.4606 |

Config G:

| α | mean ARI | mean Sil | k̂ correct | mean rank | mean ARI, 9 sets |
|---|---|---|---|---|---|
| 20% | 0.5244 | 0.3959 | 2 | 2.43 | 0.4816 |
| 15% | 0.5007 | 0.3806 | 3 | 2.50 | 0.4451 |
| 10% | 0.5013 | 0.4005 | 2 | 3.36 | 0.4519 |
| 5% | 0.4957 | 0.4240 | 2 | 3.14 | 0.4358 |
| 1% | 0.4833 | 0.4437 | 2 | 3.57 | 0.4486 |
| script | 0.5232 | 0.4251 | 4 | | 0.4606 |

Config S:

| α | mean ARI | mean Sil | k̂ correct | mean rank | mean ARI, 9 sets |
|---|---|---|---|---|---|
| 20% | 0.5226 | 0.3959 | 2 | 2.43 | 0.4937 |
| 15% | 0.4970 | 0.3807 | 3 | 2.43 | 0.4740 |
| 10% | 0.4943 | 0.3985 | 2 | 3.14 | 0.4465 |
| 5% | 0.4530 | 0.3976 | 2 | 3.71 | 0.3842 |
| 1% | 0.4879 | 0.4470 | 2 | 3.29 | 0.4505 |
| script | 0.5232 | 0.4251 | 4 | | 0.4606 |

### Config A per data set (ARI / Sil / k̂; true k in brackets)

| data set | 20% | 15% | 10% | 5% | 1% |
|---|---|---|---|---|---|
| iris [3] | 0.593 / 0.418 / 3 | 0.313 / 0.280 / 10 | 0.599 / 0.445 / 3 | 0.572 / 0.453 / 3 | 0.626 / 0.449 / 3 |
| seeds [3] | 0.414 / 0.346 / 2 | 0.470 / 0.363 / 2 | 0.454 / 0.373 / 2 | 0.509 / 0.467 / 2 | 0.402 / 0.403 / 2 |
| user_knowledge [4] | 0.201 / 0.160 / 2 | 0.104 / 0.155 / 2 | 0.104 / 0.155 / 2 | 0.147 / 0.165 / 2 | 0.083 / 0.156 / 2 |
| wholesale [3] | 0.003 / 0.035 / 31 | 0.002 / 0.066 / 23 | -0.000 / 0.055 / 28 | -0.006 / 0.112 / 2 | -0.003 / 0.376 / 2 |
| breast_cancer [2] | 0.258 / 0.108 / 2 | 0.094 / 0.009 / 16 | 0.168 / 0.247 / 2 | 0.075 / 0.208 / 2 | 0.342 / 0.278 / 2 |
| aggregation [7] | 0.406 / 0.414 / 3 | 0.407 / 0.414 / 3 | 0.391 / 0.409 / 3 | 0.238 / 0.392 / 2 | 0.313 / 0.444 / 2 |
| asymmetric [5] | 0.736 / 0.606 / 4 | 0.883 / 0.604 / 5 | 0.725 / 0.612 / 4 | 0.740 / 0.603 / 4 | 0.926 / 0.623 / 5 |
| R15 [15] | 0.951 / 0.729 / 15 | 0.975 / 0.742 / 15 | 0.943 / 0.724 / 15 | 0.959 / 0.732 / 15 | 0.795 / 0.643 / 13 |
| D31 [31] | 0.761 / 0.476 / 28 | 0.749 / 0.444 / 31 | 0.634 / 0.426 / 23 | 0.520 / 0.432 / 16 | 0.554 / 0.456 / 17 |

Config A equals the published rerun in the script's own cell for iris (10%), seeds,
breast_cancer (5%), aggregation, R15 and D31 (15%) and asymmetric (20%), where the
script's table is the shipped table of the same d. The A cells for user_knowledge and
wholesale use d = 5 and d = 6 tables, not the script's d = 6 and 8 (that is config P).

### Result

- Chosen universal α (primary, config A): 20%, mean ARI 0.5226 over the seven data
  sets. The next level, 15%, has 0.4994, a gap of 0.023. No level is within 0.001 of the
  best, so there is no tie.
- 20% also has the best mean rank (2.43; 15% 2.50) and the best mean ARI over all nine
  data sets (0.4803).
- Mean Sil is highest at 1% (0.4437) and k̂ correct is highest at 15% (3 of 7). The
  primary and secondary criteria do not agree on one level.
- Reference: the script's per-data-set levels give mean ARI 0.5232, Sil 0.4251 and
  4 of 7 k̂ correct. That is 0.0006 above the best single α in config A, and its k̂
  correct count is above every single α (at most 3).
- Sensitivity: G (0.5244) and S (0.5226) also pick 20%, with no tie. P picks 15%
  (0.5095), with 10% second (0.5056) and 20% fourth (0.4995). The choice changes
  under P only.
- Under G the six cells above that differ from A move D31 by up to +0.049 and
  breast_cancer 5% by +0.140; the 20% choice does not change.

### Run times

Step 14: 161.2 min (09:06 to 11:48). Step 15: 4.2 min for 110 tasks on 20 workers
(sum of task times 2768 s, longest task 250 s). Step 16: about 1 s. Runner: 09:06:41 to
11:52:05 (2 h 45 min).

## Post hoc: levels by dimension (not pre-registered)

The authors asked on 2026-10-20 for rules that set α by d. Each is scored from the
cells above; nothing was rerun. They were proposed after the results were seen.
- 5% for every d < 10: the 5% row above (mean ARI 0.492).
- 10% for d < 5, 5% for 5 <= d < 10: config A 0.507, P 0.501, G 0.514; 2 of 7 k̂.
- 15% for d = 2, 10% for d = 3 and 4, 5% for 5 <= d < 10: these are the simulation
  levels of Supplement S1 at d = 2, 3 and 5. Config A 0.551, P 0.544, G 0.552; k̂ correct
  4 (A, G) and 5 (P) of 7. Against the script's levels, only asymmetric changes
  (20% to 15%: 0.736 / 0.606 / 4 to 0.883 / 0.604 / 5). In config A, user_knowledge and
  wholesale also change, because the tables have their own d (0.147 / 0.165 / 2 and
  -0.006 / 0.112 / 2).

## Decision (2026-10-23)

The authors adopted the d-based rule above (15% for d = 2, 10% for d = 3 and 4, 5% for
5 <= d <= 9, with tables of the data's own d) for the UN-CCD real-data script on 2026-10-23. The three script
tables were changed: asymmetric `2d_80%` to `2d_85%`, user_knowledge `6d_95%` to `5d_95%`, wholesale
`8d_95%` to `6d_95%`; the other six data sets already followed the rule. `11_real.R` and
`12_real_reference.R` were rerun for these three data sets. They equal each other and the
sweep's config-A cells (user_knowledge 5%, wholesale 5%, asymmetric 15%) to 15 significant
digits in ARI, Sil and k̂. The old rows are in `results/tr3/_before_alpha_adopt_20261023/`.

Note (2026-10-23, later the same day): the authors also set the UN-CCD search direction by d (ascending for
d < 10, descending for d >= 10; README, "Search direction by d"). The sweep above ran the real data with the
script's directions, so seeds and breast_cancer descended in config A; config S is their ascending run. With
seeds ascending, the adopted rule's primary score changes. It is the mean ARI over the seven paper data sets,
config A at the adopted levels (iris 10%, seeds 5%, user_knowledge 5%, wholesale 5%, asymmetric, R15 and D31
15%), with seeds taken from config S at level 5 (`real_alpha.csv`):
- seeds at 5%: 0.5087 / 0.4670 / 2 (descend, config A) to 0.2394 / 0.2860 / 2 (ascend, config S);
- mean ARI 0.5509 (0.550878) to 0.5124 (0.512413); the gap is the seeds change, 0.2693 / 7 = 0.0385;
- mean Sil 0.4255 to 0.3997; k̂ correct stays 4 of 7 (seeds k̂ = 2, true k = 3);
- the other six cells are unchanged, and the new mean equals the mean of the seven UN-CCD rows now in
  `real.csv` (0.512413).
Against the single levels the ranking changes. Config S, with seeds ascending, has 20% 0.5226, 15% 0.4970,
10% 0.4943, 5% 0.4530, 1% 0.4879. The rule's 0.5124 is 0.0102 below the best single α (20%); with seeds
descending it was 0.0283 above it (0.5509 against 0.5226). Nothing else in this file was recomputed; the
"Result" and "Post hoc" sections describe the sweep as run.
