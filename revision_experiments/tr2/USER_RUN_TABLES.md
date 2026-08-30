# NN quantile-table production runs — command sheet

T3 Phase B. Three NN tables, d ∈ {166, 274, 400}, q = 0.999, niter = 10000,
n = 1000, cores = 20, canonical outdir `R/NN-test_quantile/`. You run these;
nothing here is scheduled automatically.

The generator (`revision_experiments/01_gen_quantile_table.R`) now routes every
NN run through an optimized override (`01e_nn_fast.R`): `matrix(rnorm(n*d))` in
place of `mvrnorm(n, 0, diag(d))` (distributionally identical — the mvrnorm
rotation is orthogonal and the isotropic Gaussian is rotation-invariant; kills
one O(d^3) eigendecomposition per draw call), and a streaming SimuOnce that
never materializes the per-iteration dataset list (measured per-worker peak at
d=400: ~2.2 GB original → ~0.2 GB fast, so cores=20 is safe at every d).
Validation evidence: `revision_experiments/results/nn_fast_validation.log`
(equivalence checks A/B), `nn_fast_validation_cd.log` (timing at d=166, RAM +
timing at d=400), `nn_mc_yardstick.log` (Monte-Carlo-noise tolerance baseline).
Measured statistical agreement original-vs-fast is within (at d=10, better
than) the original-vs-original different-seed baseline at the same scale.

## The three commands (PowerShell, from the repo root)

**R note (2026-07-11):** R on this machine was upgraded to 4.6.1 (4.4.1 is gone)
and plain `Rscript` is not on PowerShell's PATH. The commands below therefore
call the full path. The 4.6.1 library has been verified end-to-end (all 17
packages reinstalled where needed, `00_env_check.R` prints ENV CHECK OK, and
the exact NN-166 command path below completed a toy run at niter=2).

Run from
`G:\Submissions\TR2\Pattern Recognition (Elsevier) - Resubmit\Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection`.
One at a time — a single job at cores=20 finishes the set faster than two
jobs splitting the same 24 cores, and keeps each log's timing clean.

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "revision_experiments/01_gen_quantile_table.R" NN 166 0.999 10000 20 | Tee-Object -FilePath "revision_experiments/results/table_gen_166.log"
```

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "revision_experiments/01_gen_quantile_table.R" NN 274 0.999 10000 20 | Tee-Object -FilePath "revision_experiments/results/table_gen_274.log"
```

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "revision_experiments/01_gen_quantile_table.R" NN 400 0.999 10000 20 | Tee-Object -FilePath "revision_experiments/results/table_gen_400.log"
```

No sixth argument = the generator's production default outdir, which is the
canonical `R/NN-test_quantile/`. Do not pass `--%` anywhere (breaks Rscript
argument handling on this machine).

## Duration estimates (fresh, from post-optimization measurements)

| d   | est. wall clock at cores=20 | basis |
|-----|------------------------------|-------|
| 166 | ~14 h | Phase A measured marginal rate (5.55 s/iter -> 15.4 h) / measured head-to-head speedup 1.082x (2026-07-11) |
| 274 | ~28 h | Phase A component model (31.5 h) / model speedup ratio 1.128x; the model's ratio predictions were verified against measurement at d=166 (pred 1.065 vs meas 1.082) and d=400 (pred 1.241 vs meas 1.287) |
| 400 | ~37-42 h | measured fast serial iteration 155.3 s (single worker, 2026-07-11) x 1.69 self-contention x 500 rounds = 36.5 h; idle-model-scaled upper figure 41.7 h. cores=20 now viable (was RAM-capped at 16) |

Total ~80-84 h (~3.5 days) sequential. These assume the machine is otherwise
lightly loaded; heavy concurrent compute inflates them roughly proportionally
(the Phase A probes measured ~1.9x under a competing 20-worker job). The 1.69x
self-contention factor is the Phase A measurement at d=166/20 workers; actual
contention at d=274/400 may differ somewhat in either direction.

## No checkpointing — plan accordingly

The generator does NOT checkpoint. All niter=10000 iterations live inside one
`foreach %dopar%` call; the only save is a single `save(simul, ...)` at the
very end. If the run dies at hour N (power, reboot, OOM, closed console), you
get nothing and rerun from scratch. Corollaries:

- Run each table in a console that survives you logging off (or just leave the
  session up). A `Start-Process`-detached shell or a plugged-in laptop policy
  that blocks sleep both work; the important thing is the R process must live
  the full duration.
- Disable automatic Windows-update restarts for the window, or at least check
  the pending-restart state before launching the d=400 run.
- The three runs are independent — losing one doesn't affect the others.

## What "done" looks like

Each run ends with lines of the form

```
Result shapes: average length 1000 | median length 1000
Saved: .../R/NN-test_quantile/NN-test-simul_<d>d_999%.RData
TOTAL_ELAPSED_SECONDS=...
=== DONE ===
```

and produces exactly one file:

| d   | file |
|-----|------|
| 166 | `R/NN-test_quantile/NN-test-simul_166d_999%.RData` |
| 274 | `R/NN-test_quantile/NN-test-simul_274d_999%.RData` |
| 400 | `R/NN-test_quantile/NN-test-simul_400d_999%.RData` |

Object contract: one object `simul`, a list with `$average` and `$median`,
both numeric length 1000, NA-free, non-negative, element 1 == 0 (the n=1
placeholder), entries decreasing-ish in j (lower-tail 0.001 quantile of mean/
median NN distance among j points in the unit ball). One-line verification per
table (swap the d):

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "load('R/NN-test_quantile/NN-test-simul_166d_999%.RData'); stopifnot(is.list(simul), lengths(simul[c('average','median')])==1000, !anyNA(simul$average), !anyNA(simul$median), simul$average[1]==0); cat('OK: 166d table sane. avg[1000]=', simul$average[1000], ' med[1000]=', simul$median[1000], '\n')"
```

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "load('R/NN-test_quantile/NN-test-simul_274d_999%.RData'); stopifnot(is.list(simul), lengths(simul[c('average','median')])==1000, !anyNA(simul$average), !anyNA(simul$median), simul$average[1]==0); cat('OK: 274d table sane. avg[1000]=', simul$average[1000], ' med[1000]=', simul$median[1000], '\n')"
```

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "load('R/NN-test_quantile/NN-test-simul_400d_999%.RData'); stopifnot(is.list(simul), lengths(simul[c('average','median')])==1000, !anyNA(simul$average), !anyNA(simul$median), simul$average[1]==0); cat('OK: 400d table sane. avg[1000]=', simul$average[1000], ' med[1000]=', simul$median[1000], '\n')"
```

## Scheduling notes

- **Not while WP4 timing grids run.** `04_wp4_runtime.R` (and the PyOD twin)
  are wall-clock measurements; a 20-worker table generation on the same box
  poisons them. Either order is fine, just not simultaneous.
- Ordinary desktop use is fine alongside these runs: cores=20 of 24 leaves 4
  logical processors free, and per-worker RAM is no longer a constraint
  (measured peak ~0.2 GB/worker at d=400 with the streaming path, ~4 GB total
  for 20 workers on the 64 GB machine; the original path needed ~2.2 GB/worker
  = ~44 GB, which is what capped d=400 at 16 cores before).
- Overnight sequential is the intended pattern: launch d=166 one evening,
  d=274 the next, d=400 over a weekend (or chain them in one console with
  `&&` between the three commands and a single long-lived session).

## After generation

Tell the assistant/orchestrator the tables are in. The pipeline picks them up
without further edits:

- `06_wp5_highdim.R` — rerun fills the UNCCD-OOS/UNCCD-IOS rows for
  Musk (166), Arrhythmia (274), Speech (400); the existing SKIPPED_NO_TABLE
  marker rows are replaced automatically.
- `04_wp4_runtime.R` with mode `r_only` — fills the CCD cells at
  d = 166/274/400 in the runtime d-grid.

## Appendix (currently skipped): InternetAds, d = 1555

Skipped per scope decision 2026-07-11 (FINDINGS.md 3c) — run only if we
revisit InternetAds CCD experiments. At niter=2000 (the low-d production
tier's iteration count; needs explicit sign-off before use since the high-d
tier is niter=10000):

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "revision_experiments/01_gen_quantile_table.R" NN 1555 0.999 2000 20 | Tee-Object -FilePath "revision_experiments/results/table_gen_1555.log"
```

Rough estimate ~45 h (~2 days) wall clock at cores=20 with the fast path
(idle iteration model 3598.9 s x 0.31 after removing the eigen share, x ~0.85
measured model-to-reality correction, x 1.69 self-contention x 100 rounds).
The optimization helps most here — at d=1555 the eigendecomposition was ~69%
of the original per-iteration cost, and streaming lifts the old 5-worker RAM
cap to the full 20. Produces
`R/NN-test_quantile/NN-test-simul_1555d_999%.RData`.
