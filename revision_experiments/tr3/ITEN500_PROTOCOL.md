# Protocol: RK-CCDs at d = 3, n = 500 with 500 repetitions

Written 2026-10-18, before any run of the `iteN500` variant. The authors
asked for this rerun on 2026-10-18.

## Question

The supplement says every simulation cell is a mean over 500 repetitions.
Two RK-CCD cells are means over only 50: the scripts
`simulations/clustering/RK_CCDs/Simulations/{Uniform,Gaussian}/3d/3d_235clusters_n500.R`
set `iteN = 50`. Their published logs are `slurm-1442239.out` (uniform, 3.05 h)
and `slurm-1443038.out` (Gaussian, 24 min). These are the d = 3, n = 500
entries of Tables S3 and S7. All the other 118 cells, and the UN-CCD and
KS-CCD cells at d = 3, n = 500, use 500.

This run gives the two cells 500 repetitions, as the paper states.

## Cells and settings

| method | distribution | d | n | clusters | repetitions |
|---|---|---|---|---|---|
| RK-CCDs | uniform | 3 | 500 | 2, 3, 5 | 500 |
| RK-CCDs | Gaussian | 3 | 500 | 2, 3, 5 | 500 |

Every other setting is the published one, from these two cells' rows of
`settings_simulation.csv` (`00_check_settings.R` checks them against the
scripts):
- seed 2024;
- table `RK-test-simul_3d_99%.RData`, `quant = 0.99`;
- `cls_dis` 3 (uniform) or 6 (Gaussian), `sig = 1` (Gaussian);
- radii in (0.8, 1.2);
- the first of three clusters has `round(n/3) - 1` points.

`10_sim.R` applies the variant: `VARIANTS$iteN500` sets `iteN = "500"` and
refuses any other cell.

## Data

The script order is kept: `set.seed(2024)`, then 500 data sets with 2
clusters, then 500 with 3, then 500 with 5. As a result:

- the 1500 data sets are the ones the UN-CCD and KS-CCD cells at d = 3,
  n = 500 used, so the three methods can be compared replicate by replicate;
- the 2-cluster data sets 1-50 are the 50 published ones. The published 3-
  and 5-cluster data sets do not occur in this run, because with iteN = 50 they
  were drawn from a different point of the random stream.

## Runs

From the clone root, in PowerShell:

```powershell
$R = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
& $R revision_experiments/tr3/00_check_settings.R
& $R revision_experiments/tr3/10_sim.R RK-CCDs gaussian 3 500 2,3,5 iteN500
& $R revision_experiments/tr3/10_sim.R RK-CCDs uniform 3 500 2,3,5 iteN500
& $R revision_experiments/tr3/21_iteN500_check.R
```

Expected time on this machine (24 cores) is about 40 min for the Gaussian
cell and about 4 h for the uniform cell. These are ten times the published
50-repetition times, divided by the 6.5 to 8.4 times speed-up measured at
n = 50.

Each cluster count is written when it finishes. If a run stops (the G: drive
can drop out), the same command resumes it at the first unfinished cluster
count.

## Checks (all must pass before any number is used)

`21_iteN500_check.R` runs them and stops at the first failure:

1. **Complete:** 500 replicates for each cluster count, with no NA.
2. **Reproduction:** for 2 clusters, the means over replicates 1-50 (ARI,
   Sil, SR) equal the first result line of the published log to all 15
   printed digits.
3. **Pairing:** every data fingerprint equals that of the same replicate of
   the UN-CCD and KS-CCD data for the cell.

If a check fails, stop and report it. Nothing is tuned or rerun to force a
match.

## Output

- `results/tr3/sim_replicates.csv` and `sim_done.csv`: rows with
  `variant = iteN500`.
- `results/tr3/iteN500_check.csv`: the new means (17 digits and 3 dp) beside
  the 50-repetition values now in Tables S3 and S7.

## What this does not decide

The authors decide two things:
- whether the new means replace the d = 3, n = 500 entries of Tables S3 and
  S7. Any such edit is in red (AGENT.md), and the line plot redrawn from the
  notebook would change at that point;
- whether the two scripts' `iteN = 50` becomes 500.

## Outcome (added 2026-10-18, after the runs)

Both runs finished without errors, retries or resumes on 24 cores. All three
checks pass: `21_iteN500_check.R` was run by the run's agent and again
independently.

| Cell | Run time | Estimate above |
|---|---|---|
| Gaussian | 1.33 h | 40 min |
| Uniform | 5.28 h | 4 h |

At n = 500 this machine is about 3 times faster than the cluster, not 6.5 to
8.4 times as at n = 50.

New means (500 repetitions) against the values now in the tables. The full
17-digit values are in `results/tr3/iteN500_check.csv`.

| Table | Clusters | ARI new / table | Sil new / table | SR new / table |
|---|---|---|---|---|
| S3 uniform | 2 | 1.000 / 1.000 | 0.666 / 0.669 | 1.000 / 1.000 |
| S3 uniform | 3 | 1.000 / 1.000 | 0.656 / 0.654 | 1.000 / 1.000 |
| S3 uniform | 5 | 0.998 / 0.999 | 0.598 / 0.600 | 1.000 / 1.000 |
| S7 Gaussian | 2 | 0.989 / 0.990 | 0.636 / 0.633 | 1.000 / 1.000 |
| S7 Gaussian | 3 | 0.987 / 0.987 | 0.623 / 0.624 | 1.000 / 1.000 |
| S7 Gaussian | 5 | 0.974 / 0.976 | 0.582 / 0.581 | 0.992 / 1.000 |

Every change is at most 0.003, except the Gaussian 5-cluster SR (0.992
against 1.000). The 3-cluster ARIs, which are the points the line plots use,
are unchanged at 3 dp.

Decision (2026-10-18): the authors adopted the new values. Tables S3 and S7 were
updated in red, and both scripts now set `iteN = 500`.
