# Protocol: stability of the NND-based detectors to the seed of the Monte Carlo quantile table

Written 2026-10-10, before any run of `99_table_seed_stability.R`.

## Question

Given a stored Monte Carlo quantile table the four detectors are deterministic
(`45_wp2b_seed_stability.R --determinism`: same data, different `set.seed()`,
bit-identical scores, radii, clusters and deltas; the random-number stream is
never consumed). The only randomness in the NND-based detectors is in building
the NND quantile tables (mean and median nearest-neighbour distance of a
uniform sample in the unit ball, 10,000 draws per column). Reviewer 3 asked how
the radii, clusters and outlier labels vary across random seeds. This
experiment rebuilds the NND tables under independent seeds, at the full
production depth, reruns two NND-based detectors on real data sets and compares
each result with the production table.

## Data sets and methods

| data set | d | n | UN-MCCD table (alpha = 1%, label `99`) | SUN-MCCD table (alpha = 0.1%, label `999`) |
|---|---|---|---|---|
| hepatitis | 19 | 74 | yes | yes |
| lymphography | 18 | 148 | yes | yes |
| vowels | 12 | 1452 | yes | yes |

Loaded with `load_real_dataset()`, called through `METHOD_REGISTRY` exactly as
`run_cell()` in `45_wp2b_seed_stability.R`: UN-MCCD with no `min.cls`;
SUN-MCCD with `min.cls = 0.05` and the default ascending radius search. The
level of each method comes from the harness resolvers
(`nn_quant_label_paper_UN`, `nn_quant_label_paper_SUN`), not from this script.
RK-based methods (U-MCCD, SU-MCCD) are out of scope.

## Tables

For every cell (data set x method) the comparison base is the **production**
table, the file that produced the published numbers.

| id | what | where |
|---|---|---|
| `prod` | production table | `R/NN-test_quantile/NN-test-simul_<d>d_<label>%.RData` |
| `regen` | independent fresh 1% table (UN-MCCD only; the production 0.1% tables *are* the installed copies of the `regen999` 0.1% tables, byte-identical, so there is nothing independent to compare for SUN-MCCD there) | `R/NN-test_quantile_regen999/<d>d/NN-test-simul_<d>d_99%.RData` |
| `seed<k>` | new table, seed index k = 1..5 at d = 19 and d = 18; k = 1 at d = 12 | `R/NN-test_quantile_seedstab/<d>d/seed<k>/NN-test-simul_<d>d_{99,999}%.RData` |

New tables: 10,000 Monte Carlo draws each (the production depth, matching
`71_regen_nnd_alpha001.R`), engine `fast_stream`, 22 worker processes, chunks
under `.../seed<k>/chunks/`. One run yields both the 1% and the 0.1% table from
the same draws (same estimator as `71`: lower-tail quantile at 1 - quant of the
pooled draws, once, after pooling). Entry j of a table depends only on j
points, so a table built at the data set's own n is, in distribution, the
first n entries of a longer one (the argument of `71`). No file in
`R/NN-test_quantile/` or `R/NN-test_quantile_regen999/` is written.

The d = 12 table (n = 1452, about 7.7 h on 22 cores) is generated later by the
main session, not by the run that writes this protocol; the vowels cells that
need it are run afterwards.

## Seeds

Chunk i of seed index k at dimension d uses

    seed(k, d, i) = 20261010 + 100000 * k + 1000 * d + i

All such seeds are at least 20373010, which is above every seed used by
`71_regen_nnd_alpha001.R` (`20260811 + 1000 d + i`, at most 20281851),
`64_nn_provenance_gen.R` (`20260810 + 1000 r + d`), the default
`20260927 + d` of `R/ccds/quantile_table.R` and `20260809 + d` of
`01i_nn_multiquant_table.R`. The script prints the set-difference check
(`--generate`, including `--dry-run`) and aborts on any overlap. Each table's
metadata and `tables.csv` record n, d, niter, cores, engine, the per-chunk seed
vector and the number of chunks; a `PLAN.txt` in each `seed<k>` directory fixes
the chunk plan, and a resume with a different number of cores is refused.

## Metrics

For each (data set, method, table), with the production table as base:

* share of the n points whose outlier label differs from production;
* Jaccard index of the flagged sets (1 if both are empty);
* adjusted Rand index of the `cluster` partitions, and the number of clusters of each;
* alternative minus production for BA and F2; TPR, TNR, BA, F2 of each;
* median and maximum relative change of the per-point radii (radius tap on `nnccd.radi`, original row order),
  taken over the points whose production radius is positive; the number of points
  whose radius is zero in one table and positive in the other is reported separately.

Summary per (data set, method): number of tables compared, range of F2 and BA
over all tables including production, mean and minimum Jaccard, maximum share
of labels changed, minimum ARI.

## Observability of the table used

* On-the-spot generation is switched off (`ccd.quantile.generate = FALSE` and
  `CCD_QUANTILE_GENERATE=FALSE`, set before any `source()`, and asserted after).
* `get_simul` is wrapped after all `source()` calls. For each cell the wrapper
  records the file path it resolved, asserts it equals the intended file (paths
  normalised), asserts exactly one call per cell, and asserts the loaded object
  carries no `ccd_table` or `generated` attribute (so `.ccd_qt_grow` cannot fire).
  A violation stops the script; it is not turned into an error row.
* Each result row holds the path, the file's md5 and the sum of the table's
  `average` and `median` vectors. A comparison counts as evidence only if the
  md5 of the alternative differs from the production md5 (hard stop otherwise).

## Reporting

Results are reported whatever they show. No published number changes: the
production-table cells are first reproduced to three decimals (`--anchor`)
against the published values, and generation starts only if all six match.
The experiment is a sensitivity statement about the table seed, not a change
to any table or figure of the manuscript.

## Outputs

`revision_experiments/results/tr1/table_seed/`: `table_seed_cells.csv` (one row
per cell and table), `table_seed_comparisons.csv`, `table_seed_summary.csv`,
`tables.csv` (manifest with md5 of every table used).
