# WP9 protocol — SUN-MCCD component ablation (R3.10a)

Written before `wp9_sun_variants.R` / `89_wp9_ablation.R`, per the cut-list
scope (REVISION_PLAN §WP9, cut list item 1): **three toggles**, on SUN-MCCD
only, run across the full plan's grid of two dimensions (d = 3, 10) and two
generators (uniform, Gaussian) -- four settings, not a single dimension pair.
Four toggles were originally planned; Holm correction is dropped (below).

## Dropped toggle: Holm correction

CLAUDE.md ("The MC-SRT is not Holm-combined — corrected 2026-09-05"): the
live code (`R/ccds/UN_CCD.R:252-258`, critical values from
`R/ccds/NN_Dist_Est.R:41-51,101-111`) never formed p-values or applied Holm's
step-down procedure. It compares the observed mean and median NND against
their own per-statistic Monte Carlo lower α-quantiles and rejects CSR if
**either** falls below its critical value — the marginal per-statistic level
is exactly what the quantile-table filenames (`..._{95,99,999}%.RData`)
encode. There is nothing to toggle: "uncorrected" is what the code already
does, and for two hypotheses (mean, median) Holm's step-down procedure
reduces to a single comparison round at α/2 -- reject CSR if the smaller of
the two p-values is below α/2, which for the code's either/or rule is
exactly the same decision as running the existing rule at the α/2 quantile
table. That α/2 table exists at d = 3 (the α/2 = 5% table, filename suffix
`95`) but not at d = 10 (α/2 = 0.05%, no `9995`-suffixed table was ever
generated -- only `{95,99,999}` exist). More to the point, the effect of
halving α is already measured directly: the WP2(a) α sweep
(`42_wp2a_alpha_sweep.R`) varies α across the same marginal-quantile tables
this method already uses. So a dedicated Holm toggle here would be
redundant with that sweep, not impossible to build -- it is dropped for
that reason, and the response letter states it this way rather than as an
unavailable-table claim.

## The three toggles

All three act on `nnccd.radi()` (`R/ccds/UN_CCD.R:233-323`), the per-point
radius search inside `nnccd_clustering_quantile()`, which is what
`SUNMCCD_outlier()` (`methods/outlier_detection/SUN-MCCD.R`) calls.

| # | Toggle | Variants | Baseline (= stock SUN-MCCD) |
|---|---|---|---|
| 1 | NND statistic | `both` (either mean or median below its quantile rejects CSR — current) / `mean` only / `median` only | `both` |
| 2 | Centre-point removal | `on` (current: the ball's own centre point is excluded from the observed NN-distance statistic) / `off` (centre included) | `on` |
| 3 | Radius search direction | `ascend` (current: grow the ball from `low.num` neighbours outward) / `descend` (shrink from all `n-1` neighbours inward) | `ascend` |

**Design: one-factor-at-a-time, not a full 3×2×2 factorial.** Each toggle is
varied alone against the `(both, on, ascend)` baseline, matching how
`43_wp2a_direction_sweep.R` treated the direction toggle in isolation and how
the cut-list table itself is laid out (one row per toggle, not a cross
table). This gives **5 distinct variants** per setting, not 12:

| variant_id | stat | remove_centre | method |
|---|---|---|---|
| `stock` | both | TRUE | ascend | (= baseline, reproduces published SUN-MCCD bit-for-bit) |
| `stat_mean` | mean | TRUE | ascend |
| `stat_median` | median | TRUE | ascend |
| `centre_off` | both | FALSE | ascend |
| `direction_descend` | both | TRUE | descend |

`nnccd.radi.ablate()` (below) is written to accept the full `stat ×
remove_centre × method` cross so nothing prevents a later factorial run, but
the 5-cell one-at-a-time table above is what `89_wp9_ablation.R` schedules
and what the manuscript/response-letter reports.

## The centre-off indexing decision

This is the one place the cut-list spec leaves open, and it must be decided
before the code is written, not discovered by trial.

**Background.** The Monte Carlo envelope (`NN.envelop$average`,
`NN.envelop$median`, built by `NNDest.simpois.lower.quant()`,
`R/ccds/NN_Dist_Est.R:17-49`) is indexed so that **entry k is the simulated
NN-distance statistic (mean or median) of a k-point i.i.d. CSR cloud** —
`NN.dist.temp.ave = c(0, NN.dist.temp[1,])` pads entry 1 with a dummy 0 (a
1-point cloud has no NN distance) and entry k≥2 is the statistic over a
genuine k-point simulated cloud.

**Ascend branch** (`UN_CCD.R:245-264`). At step `j` (ball holding the
`j`-nearest neighbours of point `i`), centre-removed, the observed statistic
is computed over `ddx[o.d[2:j], o.d[2:j]]` — this excludes `o.d[1]` (point
`i` itself) and is a cloud of **`j-1` points**. The code compares it against
`NN.envelop$average[j-1]` / `[j-1]`. **The envelope index already equals the
observed cloud size** (`j-1` points → entry `j-1`) — the on-variant's index
choice is not arbitrary, it is size-matched.

Turning centre removal off means the observed statistic is computed over
`ddx[o.d[1:j], o.d[1:j]]` instead — the same ball, but now **including**
point `i`, a cloud of **`j` points** (one more than the on-variant, because
the centre itself is now counted). Preserving the same size-matching
invariant means comparing against entry **`j`**, not entry `j-1`.

**Decision: the off-variant reads envelope entry `[j]`, one entry higher
than the on-variant's `[j-1]`, because the observed cloud gained exactly one
point (the centre) and the envelope's convention is "entry k models a
k-point cloud."** Reusing `[j-1]` for a `j`-point observed statistic would
silently compare a `j`-point empirical statistic against a simulated
`(j-1)`-point null — a size mismatch that would bias the test toward
appearing *more* significant than it should (smaller simulated clouds have
larger NN-distance quantiles), which is exactly the kind of silent unit
error CLAUDE.md's δ₀/Δ and S_min notes warn about. It is not left as a free
parameter.

**Caveat.** The Monte Carlo envelope is built by simulating i.i.d. CSR
clouds (`NNDest.simpois.lower.quant()`). With centre removal on, the
observed cloud genuinely is a set of points free to fall anywhere, matching
that null. With centre removal off, one of the j (or size) observed points
is not free -- it sits, by construction, exactly at the ball's own centre,
since it is the point the ball was grown around. Reading the size-matched
envelope entry corrects the *count* but not this placement constraint, so
the j-point i.i.d. CSR envelope is an approximate null for the centre-off
variant, not an exact one. `centre_off`'s results should be read with that
in mind.

**Descend branch** (`UN_CCD.R:265-279`), for completeness even though the
5-cell grid above never combines `centre_off` with `descend`.
`nnccd.radi.ablate()` implements the general rule rather than special-casing
ascend: whatever entry the on-variant's formula would read at this step, the
off-variant reads **one entry higher** (in envelope-index terms, not in
`rev()`-subscript terms — `rev()`'s subscript *decreases* by 1 to move the
underlying index *up* by 1). Concretely, the stock descend line reads
`rev(NN.envelop$average)[j+2]`; the off-variant reads `rev(...)[j+1]`. This
was verified algebraically before writing code: `rev(v)[j+2] = v[n-j-1]`
against an observed centre-removed cloud of size `n-j`, i.e. index =
`size - 1` (the descend branch's own, different, pre-existing offset
convention — not something this ablation changes); the off-variant's cloud
grows to `n-j+1` points, so its index should be `(n-j+1)-1 = n-j =
rev(v)[j+1]`. The `+1` shift is identical in both directions; only the base
formula differs, because ascend and descend already used different
size-to-index offsets before this toggle existed.

## Settings

d ∈ {3, 10}, n = 200 (nominal), uniform and Gaussian cluster generators, 100
replicates, S_min = 0.05 (CLAUDE.md's current real-data/simulation constant,
not the old half-contamination oracle rule), α resolved per d by
`nn_quant_label_paper_SUN(d)` (harness.R) — "90" (10%) at d=3, "999" (0.1%)
at d=10, low.num = 3 (SUN-MCCD's own default).

**Contamination fixed at 5%** ("the study's usual contamination"): this is
the level baked into the original simulation drivers' own filenames
(`.../10d_2cls_n500_cont5%.R`, cited verbatim in
`55_wp2c_simulation_arm.R`'s header) and is the drivers' un-parameterised
default before WP2c promoted contamination to a swept argument. WP9 is an
ablation of algorithmic components, not a contamination sweep (that is
WP2a/WP2c's job), so one representative level is used throughout.

Generators (`gen_uniform`, `gen_gaussian` in `wp9_sun_variants.R`) are a
verbatim copy of `55_wp2c_simulation_arm.R`'s copies (2 clusters, inter-centre
distance `cls_dis=3`, outlier stand-off `otl_dis=2`, radius jitter
`r_min/r_max=0.7/1.3`, Gaussian `noise_level=0.01`, outlier radius 5), with
only `n`, `d`, `cont` promoted to arguments — same lineage note as that file.

Four settings: `{uniform, gaussian} × {d=3, d=10}`.

## Cells

4 settings × 5 variants × 100 reps = 2000 detector calls for the full run.
`--smoke` runs 1 setting (`uniform_d3`) × 5 variants × 1 rep = 5 calls, plus
one bit-identity assertion on the `stock` cell.

Seeds: `seed = WP9_BASE_SEED + rep_id`, `WP9_BASE_SEED = 9100` (distinct from
`55_wp2c_simulation_arm.R`'s `123` so the two experiments' replicate streams
never collide if ever compared side by side). Recorded per output row.

## Substitution mechanism (no global override)

`SUNMCCD_outlier()` calls `nnccd_clustering_quantile()` unqualified, which
itself calls `nnccd.radi()` unqualified; both are ordinary lexically-scoped
lookups resolved in the calling function's closure environment (`.GlobalEnv`,
since both were defined by a top-level `source()`). Rather than
`assign("nnccd.radi", ..., envir = .GlobalEnv)` — which is exactly the
already-diagnosed failure mode ("sourcing the detector files re-sources the
originals and silently reverts global overrides") — `run_sunmccd_ablation()`
builds a **child environment chain** and rebinds a *copy* of each dispatcher
function's `environment()` to point into it:

```
env_radi  <- new.env(parent = environment(nnccd_clustering_quantile))  # parent = .GlobalEnv
env_radi$nnccd.radi <- <closure over nnccd.radi.ablate with stat/remove_centre baked in>

nccq2 <- nnccd_clustering_quantile;  environment(nccq2) <- env_radi

env_sun <- new.env(parent = environment(SUNMCCD_outlier))              # parent = .GlobalEnv
env_sun$nnccd_clustering_quantile <- nccq2

sun2 <- SUNMCCD_outlier;  environment(sun2) <- env_sun
sun2(datax, simul = simul, min.cls = min.cls, method = method, low.num = low.num)
```

`SUNMCCD_outlier` and `nnccd_clustering_quantile` in `.GlobalEnv` are never
touched — `sun2`/`nccq2` are new function objects (R functions are values;
`environment(f) <- e` on a copy does not mutate the original binding).
Every other symbol either dispatcher references (`nnccd.silhouette`,
`connected.ksccd.m`, `ksccd.connected`, `dominate.mat.greedy2`, `dist`, …)
resolves by walking the parent chain up to `.GlobalEnv` exactly as it always
did, so nothing besides the one targeted call is affected. Re-sourcing
`SUN-MCCD.R` / `UN_CCD.R` (as `wp0_mccd_methods.R`'s idempotent
`if (!exists(...))` guards do) cannot revert this, because there is nothing
in `.GlobalEnv` to revert.

The `stock` variant (`stat="both", remove_centre=TRUE, method="ascend"`)
must reproduce plain `SUNMCCD_outlier()`'s output **bit-for-bit** on the
smoke data — same score vector, same radii — because `nnccd.radi.ablate()`
at those settings is definitionally identical to `nnccd.radi()`. The smoke
run asserts `identical()` on both the translated score vector and the raw
radius vector between `sunmccd_ablate_method(..., stat="both",
remove_centre=TRUE, method="ascend")` and stock `sunmccd_method(...)`
(`wp0_mccd_methods.R`), at matching `min.cls`, `low.num`, `simul`, on the
same data. A mismatch would mean the substitution mechanism itself is broken
and every other variant's result would be meaningless.

## Outputs

- `results/tr1/wp9/smoke/wp9_ablation_smoke.csv` — 5 rows from `--smoke`.
- `results/tr1/wp9/wp9_ablation.csv` — full run, one row per (setting,
  variant, rep): `setting_id, generator, d, n, n0, rep, seed, variant_id,
  stat, remove_centre, method, TPR, TNR, BA, F2, n_flagged, n_clusters,
  mean_radius, elapsed_sec, status, timestamp`. Appended per cell
  (`append_result`/`has_result` from `harness.R`), so a restart after a G:
  drive drop skips completed cells rather than recomputing or losing them.
- `results/tr1/wp9/wp9_ablation_summary.csv` (built by `--summarize`, not run
  in this session) — one row per (setting, variant): mean/sd over reps of
  TPR, TNR, BA, F2, plus mean number of clusters found and mean radius, so
  the mechanism (not only the score) is visible per the task spec.

All three files live under the gitignored `revision_experiments/results/`
tree; only the two `.R` scripts and this protocol are committed.

## What this WP9 pass does NOT do

Per REVISION_PLAN §WP9, the remaining three R3.10-named components (S_min,
silhouette-based k selection, approximate MDS) are handled elsewhere (WP2a;
an oracle-k comparison; noted as inherited from `manukyan2019parameter`) and
are out of scope for this ablation.
