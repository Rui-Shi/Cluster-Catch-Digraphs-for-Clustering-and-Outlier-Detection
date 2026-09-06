# WP8 re-verification (Opus, 2026-09-05) — after the fixer commits 764bc7f, 605e103, 710d030, 0e29d32

| Script | Verdict |
|---|---|
| 85 boundary FP | PASS — launched 2026-09-05 (background, `ALL 100 99999999`) |
| 86 outlier types | PASS WITH REQUIRED CHANGES (bridge construct collapsed; local arm framing) |
| 87 small cluster | PASS WITH REQUIRED CHANGES (one-line CLI fix at 87:275) |
| 88 CSR violation | PASS — launched 2026-09-05 (background, RK at 40 reps then the rest at 100, one process) |

Seed canonicalisation PASS in all four (same (setting, rep) -> same seed under any selection).
Block writes PASS; summarisers PASS (85 de-dup key bug confirmed fixed). No per-cell timeout
exists anywhere, so a slow cell cannot be misrecorded as an error. No shared-file hazard between
the four scripts. Peak R memory ~47 MB per process.

## 86 required changes

**bridge (blocking).** With clearance 0.45 and 3-D jitter, `t_lo > t_hi` in 37-41% of reps and the
chain spans ~10% of the realised gap: the 10 points are a compact blob at the midpoint, geometrically
near-identical to `collective`, and every method scores TPR = 0 on it. Fix (tested by the verifier,
0/600 violations, 80% axis span): clearance **0.05**, jitter sd **0.03 projected onto the
orthogonal complement of the cluster axis**, no redraw (perpendicular moves strictly increase the
distance from both centres, so `n_inside = 0` is analytic). Record in the protocol the acceptance
criteria and measured values: n_inside = 0 over 100 reps at both d; axis span >= 60% of the realised
gap; within-chain NN <= host cluster median NN; end gap <= ~1.3x that spacing; smoke at both d with
non-zero TPR for at least one method.

**local (framing).** Acceptance ratio holds (mean 6.20 at d=3, 2.58 at d=10) but the point sits 0.49
beyond the host cluster's realised radius (0.40): it is a sparse-outer-shell outlier, not one
"embedded within" the cluster, and all 7 methods score TPR = 1. Also both clusters are shrunk
(core_frac 0.4), so this arm's clusters are 2.6x denser than bridge/collective's. Geometric fact: at
d = 10, n = 95, median within-cluster NN 0.699 vs radius 1.009, so an interior point with doubled NN
distance is not realisable at this n; at d = 3 it is. Decision: run the arm at TWO placements — the
original interior construction (whose all-methods-zero TPR at d = 3 is itself the answer to R3.7)
and the current shell construction — and report both, stating the stand-off and the density
difference. Fix the stale comment at 86:100-107 (says core 0.5 / shell 0.75-1.0; code is 0.4 /
0.8-1.0).

## 87 required change

`87:275` — parse the reps override so that an empty/zero second argument leaves the 50/100 family
split intact: `reps_ov <- if (length(args) >= 2 && nzchar(args[2]) && !is.na(suppressWarnings(as.integer(args[2]))) && as.integer(args[2]) > 0) as.integer(args[2]) else NULL`.
Launch after the fix: `Rscript revision_experiments/tr1/87_wp8_small_cluster.R ALL 0 99999999`.

## 88 analysis-time items (do not gate the launch)

- `88:236` "best MST threshold": at d = 10, MST@1.60 and MST@2.00 flag 0.000 on the control, so
  `which.min(|delta|)` picks a threshold that flags nothing. Exclude sweep members with ~0 control
  flag rate, or report the whole sweep.
- The paired delta merges on `rep`, but generators have different seeds; the SE is a difference of
  independent means, not a paired one. Correct the protocol wording.
- RK methods at 40 reps: worst-case SE of a flag rate ~1 pp; the R3.4 effect (UN-MCCD 1.6% on CSR vs
  9.3% under the gradient) is ~8 SEs. State the reduction.

## Analysis-time items for all four

- WP7 cluster files: de-duplicate on the natural key keeping the last, and assert each cell has
  exactly n_reg (85/86) or n (87/88) distinct row_index; a crash between the metrics row and the
  cluster block can leave a done cell with no cluster rows.
- 85 writes no row for a failed cell, so a failing method silently lowers n_reps; compare expected vs
  summarised cell counts.
- 86 summariser: `n_inside_mean` is NaN for local/collective; a method failing on every rep vanishes
  from the summary instead of showing n_error > 0.
- After the run: `count.fields` integrity pass on every output CSV (four concurrent writers on the
  USB enclosure).

## Wall clock (four processes in parallel)

85 ~1.9 h; 86 ~3.7 h with the second local arm; 87 ~4.3 h; 88 ~5.0 h at 40 RK reps. Critical path
~5 h. Do not co-schedule WP6.
