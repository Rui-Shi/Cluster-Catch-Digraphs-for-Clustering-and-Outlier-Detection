# WP2(b) — Stability under randomness and under data perturbation

## RERUN UNDER THE FINAL CONFIGURATION — 2026-09-05

The 2026-08-09 run below was made under two settings that the manuscript no
longer uses:

1. `S_MIN <- 0.0625` in `45_wp2b_seed_stability.R:152`, whereas the revised
   paper uses one label-free constant `S_min = 0.05` for every data set and for
   the simulations (WP2c, user-approved 2026-08-10).
2. The NN quantile tables at d = 12, 16, 18, 19, 21 held **1%** values in files
   labelled 0.1%; they were regenerated at a genuine 0.1% level on 2026-08-13.
   *hepatitis* is d = 19, so SUN-MCCD on *hepatitis* ran against the wrong
   table.

The whole study was therefore re-run on 2026-09-05 at `S_MIN <- 0.05` against
the repaired tables. **The new outputs carry a `_smin005` suffix; the
2026-08-09 files are kept unchanged for comparison.**

* `wp2b_determinism_smin005.csv` — 20 cells, 5 data sets × 4 methods
* `wp2b_seed_stability_smin005.csv` — 1,616 rows, 30 columns, 0 non-`ok`,
  `radii_aligned` TRUE everywhere, `s_min` column holds only `0.05` / `NA`
* `wp2b_stability_summary_smin005.csv` — 32 summary rows

Wall clock, foreground, per chunk: determinism 1m47s; hepatitis (4 methods)
7m03s; glass U-MCCD 5m15s; glass SU/UN/SUN 5m28s; stamps U-MCCD 7m25s; stamps
SU-MCCD 5m13s; stamps UN/SUN 1m08s; pima U-MCCD 7m37s; pima SU-MCCD 7m22s; pima
UN/SUN 3m18s; summary 2s. Total ≈ 52 min.

### Determinism — unchanged, 20/20 on all three checks

Outcome identical across seeds 20/20; `.Random.seed` untouched 20/20; `runif(1)`
after the call identical 20/20. Nothing about the determinism verdict depends on
`S_min` or on which quantile table is loaded.

### Anchor — 16/16 PASS, now against the manuscript

`results/tr1/final_comparison.csv` is a 2026-08-09 artefact and is **stale for
two cells**: *hepatitis* × SUN-MCCD (still the pre-repair 0.446 / 0.686 / 0.714 /
0.657) and *vertebral* × SU-MCCD (pre-WP2c). `build_summary()` now overrides
those two rows from `SupplementaryMaterial.tex`
(`SM-tab:Real_Data_Result2.1` / `2.2`), which is the authority; the other 14 of
the 16 WP2(b) cells were checked one by one against the supplement and agree
with the CSV. With that correction the unperturbed cells reproduce the published
values to 3 dp on all four of BA, F2, TPR, TNR, **16/16**. In particular
*hepatitis* × SUN-MCCD now returns BA 0.705 / F2 0.469 (was 0.686 / 0.446) and
*glass* × SUN-MCCD returns 0.770 / 0.324, both matching the supplement exactly.

### What moved

Averaged over the four data sets:

| Perturbation | Method | Jaccard old → new | flip rate old → new | BA sd old → new |
|---|---|---|---|---|
| jitter | U-MCCD | 0.428 → **0.428** | 0.2377 → **0.2377** | 0.0863 → **0.0863** |
| | SU-MCCD | 0.752 → 0.754 | 0.1538 → 0.1529 | 0.0709 → 0.0714 |
| | UN-MCCD | 0.740 → **0.740** | 0.0860 → **0.0860** | 0.0603 → **0.0603** |
| | SUN-MCCD | 0.914 → 0.921 | 0.0401 → 0.0371 | 0.0403 → 0.0384 |
| dropout | U-MCCD | 0.822 → **0.822** | 0.0567 → **0.0567** | 0.0395 → **0.0395** |
| | SU-MCCD | 0.931 → 0.930 | 0.0413 → 0.0391 | 0.0457 → 0.0494 |
| | UN-MCCD | 0.906 → **0.906** | 0.0362 → **0.0362** | 0.0543 → **0.0543** |
| | SUN-MCCD | 0.968 → 0.965 | 0.0144 → 0.0169 | 0.0253 → 0.0240 |

Bold = bit-identical to the 2026-08-09 run. **U-MCCD and UN-MCCD reproduce
exactly**, which is the expected control: neither takes `min.cls`, and UN-MCCD
at d = 19 reads the 99% table, which was not among the regenerated files. Every
change is confined to five cells — *hepatitis* × {SU-MCCD, SUN-MCCD} × {jitter,
dropout} and *pima* × SUN-MCCD × jitter — and no aggregate moves by more than
0.007.

**The ranking survives intact.** SUN-MCCD is still the most stable of the four on
Jaccard, flip rate, BA sd and F2 sd under both perturbations, and U-MCCD is still
the least stable on every one of them. The gap widens marginally: SUN-MCCD's
jitter Jaccard rises 0.914 → 0.921 and its flip rate falls 0.040 → 0.037, while
U-MCCD is unchanged at 0.428 / 0.238.

The one structural change is *hepatitis* × SUN-MCCD, which now returns **2**
macro-clusters instead of 1 (both `S_min` and the repaired d = 19 table push in
that direction), so its ARI is no longer degenerate for that cell. That single
cell is what moves the SU-MCCD jitter ARI mean from 1.000 to 0.845 and the
SUN-MCCD jitter ARI mean from 0.950 to 0.984.

### The two standing qualifications still hold

* **ARI remains vacuous for the single-cluster methods.** 7 of the 8 SU-/SUN-MCCD
  cells still have `n_clusters_base = 1` (all but *hepatitis* × SUN-MCCD), and
  ARI between two one-cluster partitions is 1 by definition. **Quote Jaccard and
  the flip rate for SU-/SUN-MCCD, never ARI.**
* **The NND radius CV is still inflated by exact zeros.** Unperturbed zero-radius
  fractions: RK family 0/74, 0/213, 0/340, 0/555 (0.0% everywhere); NND family
  *hepatitis* UN-MCCD 26/74 (35.1%) and SUN-MCCD 15/74 (20.3%), *glass* 66/213
  (31.0%), *stamps* 103/340 (30.3%), *pima* 83/555 (15.0%). `radius_max_pointcv`
  saturates at √49 = 7.071 in nearly every NND cell, the signature of a point
  toggling across the zero boundary in exactly one of 50 replicates. RK
  per-point CV under jitter is ≤ 0.073 throughout. **Radii are stable; the
  connectivity built on them is not.**

  New detail: at d = 19 the UN and SUN radius statistics **no longer coincide**
  (0.173 vs 0.120 mean CV under jitter, 26 vs 15 zero radii). They coincided in
  the 2026-08-09 run only because the 0.1% and 1% d = 19 files were byte-identical
  duplicates. The RK pair still coincides exactly, as it must.

### Where this went in the manuscript

* `SupplementaryMaterial.tex` §S3.3 *Stability of the Outlier Labels Under Data
  Perturbation* — Table `SM-tab:stability` plus design and reading.
* `CCD_OutlierDetection_Neurocomputing.tex` §6 — four sentences at the end of the
  real-data discussion.
* `Response_to_Reviewers.tex` — `\point{R3.1}` with the verbatim referee text.

---

# The 2026-08-09 run (S_min = 0.0625, pre-repair d = 19 table)


Reviewer point **R3.1**: *"report how the radii, clustering results, and outlier
labels vary across parameter settings, random seeds, and numbers of Monte Carlo
simulations"*, with the observation that **graph connectivity is a discontinuous
function of the data**, so a small perturbation can flip a component and hence a
label.

Script: `revision_experiments/45_wp2b_seed_stability.R`
Outputs: `wp2b_determinism.csv`, `wp2b_seed_stability.csv` (1,616 per-run rows),
`wp2b_stability_summary.csv` (32 summary rows).

---

## Step 1 — Where does randomness enter? Verdict: **nowhere**

**The four MCCD detectors are exactly deterministic given fixed data and a fixed
quantile table, and they never touch the R random-number stream.**

### Static trace

`grep -nE 'sample\(|runif|rnorm|set\.seed|kmeans|sample\.int|jitter'` returns

* **`methods/outlier_detection/`** (RU-MCCDs.R, SU-MCCDs.R, UN-MCCD.R,
  SUN-MCCD.R) — **no matches at all**.
* **`R/ccds/`** — matches only in
  * `Kest.R` / `NN_Dist_Est.R`: `rpoisball.unit`, `rpoisbox.unit`, `mvrnorm`.
    These are the **Monte-Carlo null-distribution generators**, reachable only via
    `Kest.simpois.edge.quantile()` / `NNDest.simpois.lower.quant()`, which
    `ccd.Kest.edge.quantile()` and `nnccd.radi()` call **only when `simul` is
    NULL**. Every wrapper in `wp0_mccd_methods.R` passes a non-NULL `simul`
    loaded from `R/{RK,NN}-test_quantile/*.RData` by `get_simul()`. The branch is
    dead at inference time — the Monte Carlo draws are frozen in those files.
  * `K-S_CCD_Rui_trash.R` and the `R/*-test_quantile/*.R` table generators, none
    of which is sourced by the detector chain.

The specific candidates named in the brief are all deterministic:

| Candidate | Finding |
|---|---|
| k-means / clustering initialisation | **There is no k-means anywhere.** Macro-clusters come from a greedy approximate dominating set, not a randomly initialised partition. |
| Silhouette choice of number of clusters | `rccd.silhouette` (`R/ccds/RK_CCD_New.R:108`) and `nnccd.silhouette` (`R/ccds/UN_CCD.R:99`) are forward scans over `i = 2..lenD` keeping a strict `datasi > maxsi` maximum — ties resolve to the first index, deterministically. |
| Approximate MDS / `dom.method = "greedy2"` | `dominate.mat.greedy2` (`R/ccds/ccdfunctions.R:63`) and `dominate.mat.ks` (`:108`) select with `which.max()`. "Approximate" here means *not provably minimum*, **not randomised**. |
| Tie-breaking in greedy selection | `which.max()` — lowest index wins. Deterministic. |
| `connected.ksccd.m` | `R/ccds/mKNN_CCD_functions.R:50` — a bracketing search on a scalar. `igraph::components` and `cluster::silhouette` are deterministic. |

### Empirical proof — two independent checks, 20 cells

5 datasets spanning size and dimension (hepatitis 74/19, glass 213/9, stamps
340/9, WDBC 367/30, pima 555/8) × 4 methods. Evidence in
**`wp2b_determinism.csv`**.

**(a) Outcome check** — same data, `set.seed(1)` vs `set.seed(20260809)`,
compared with `identical()` (bit-for-bit, not to 3 dp):

| Column | Result |
|---|---|
| `score_identical` | TRUE, 20/20 |
| `radii_identical` | TRUE, 20/20 |
| `cluster_identical` | TRUE, 20/20 |
| `delta_identical` | TRUE, 20/20 |
| `metrics_identical` | TRUE, 20/20 |

**(b) RNG-stream check** — this is the check that separates *"the RNG was never
consumed"* from *"the RNG was consumed but happened not to change the answer on
these data"*. Two claims, and only (b) can distinguish them.

| Column | Result |
|---|---|
| `random_seed_touched` — `.Random.seed` captured immediately before and after the detector call, compared with `identical()` | **FALSE, 20/20** — the stream is untouched |
| `runif_identical` — `runif(1)` from `set.seed(777)` with **no** detector call vs `runif(1)` from `set.seed(777)` **with** the detector call in between | **TRUE, 20/20** — identical draws, so not one variate was consumed |

Columns `runif_seed777_no_call` and `runif_seed777_after_call` hold the two raw
draws.

**Consequence.** A seed grid is meaningless here: 50 seeds would return 50
bit-identical results. Step 2 was therefore **not run**, per the task brief.
Randomness is not where R3.1's concern lives — Step 3 is.

---

## Step 4 — Anchor check: **PASS, 16/16**

The unperturbed cell (`ptype = "none"`, rep 0) reproduces
`results/tr1/final_comparison.csv` to 3 decimal places on **all four** of BA, F2,
TPR and TNR, for every one of the 4 datasets × 4 methods.
Supporting columns in `wp2b_stability_summary.csv`: `BA_base` / `anchor_BA`,
`F2_base` / `anchor_F2`, `TPR_base` / `anchor_TPR`, `TNR_base` / `anchor_TNR`,
verdict in **`anchor_pass`** (TRUE in all 32 summary rows).

Two harness guards also pass, in `build_summary()`'s "TAP GUARD" block:

* `radii_len == n_eff` on **1,616 / 1,616** rows — the radius tap really did
  capture a per-point vector in original row order.
* `n_clusters == 0` on **0** rows — the `connected.ksccd.m` tap fired on every
  cell, so cluster counts are measured, not assumed.

---

## Step 3 — Perturbation study (the substantive answer to R3.1)

4 datasets × 4 methods × 2 perturbation types × 50 replicates = 1,600 runs,
plus 16 unperturbed reference cells. `S_min = 0.0625`. Perturbation seed
`set.seed(1000 + rep)`, recorded in the **`seed`** column of every row.

* **`jitter`** — `X[, j] += rnorm(n, 0, 0.01 * sd(X[, j]))`, i.i.d. per feature.
* **`dropout`** — leave-p-out: drop `round(0.01 * n)` rows at random.

Comparisons are against the **unperturbed baseline** run (well defined for both
perturbation types; for `dropout`, restricted to retained rows). Modal-partition
and modal-flagged-set comparisons are also recorded
(`ari_vs_modal_mean`, `jaccard_vs_modal_mean`, `n_runs_flagset_eq_modal`,
`modal_eq_base`).

### Headline, averaged over the four datasets

| Perturbation | Method | ARI | Jaccard | flip rate | BA sd | F2 sd |
|---|---|---|---|---|---|---|
| jitter σ=0.01 | U-MCCD | 0.336 | 0.428 | 0.238 | 0.0863 | 0.1313 |
| | SU-MCCD | 1.000\* | 0.752 | 0.154 | 0.0709 | 0.0761 |
| | UN-MCCD | 0.606 | 0.740 | 0.086 | 0.0603 | 0.0893 |
| | **SUN-MCCD** | 0.950\* | **0.914** | **0.040** | **0.0403** | **0.0434** |
| dropout 1% | U-MCCD | 0.838 | 0.822 | 0.057 | 0.0395 | 0.0654 |
| | SU-MCCD | 1.000\* | 0.931 | 0.041 | 0.0457 | 0.0519 |
| | UN-MCCD | 0.873 | 0.906 | 0.036 | 0.0543 | 0.0711 |
| | **SUN-MCCD** | 0.970\* | **0.968** | **0.014** | **0.0253** | **0.0299** |

**SUN-MCCD — the manuscript's headline recommendation — is the most stable of
the four on Jaccard, per-point flip rate, BA sd and F2 sd, under both
perturbations. U-MCCD is the least stable on every one of them.**

\* **Important caveat, do not quote the SU-/SUN-MCCD ARI as a stability win.**
With `S_min = 0.0625` those two detectors return a **single macro-cluster** on
all four of these data sets (`n_clusters_base = 1`; see `n_clusters_table`,
which reads `1:50` for almost every SU-/SUN-MCCD cell). ARI between two
one-cluster partitions is 1 by definition. Their ARI is therefore *degenerate*,
not evidence. The informative stability measures for SU-/SUN-MCCD are Jaccard
and the flip rate, which do move, because the MCG component split *within* the
single macro-cluster still changes.

### R3.1's discontinuity claim is confirmed

The reviewer is right, and the effect is not small. Under a jitter of just 1% of
one per-feature SD:

* **stamps × U-MCCD**: Jaccard vs. baseline falls to **0.138** (mean), min
  **0.039**; per-point flip rate **0.389**, and **90.3%** of points change label
  at least once across the 50 replicates. Flagged-set size swings from 15
  (baseline) to a range of **[9, 278]**. BA sd 0.114, F2 sd 0.204.
* **glass × U-MCCD**: ARI **0.239**, Jaccard **0.491**, flip rate 0.310, one
  point flipping in **96%** of replicates (`flip_rate_max = 0.960`).
* Cluster counts are the discontinuity made visible. **stamps × U-MCCD** finds
  **79** macro-clusters unperturbed; under jitter the count collapses to **2** in
  30 of 50 replicates and otherwise scatters over 74–84
  (`n_clusters_table = 2:30;74:2;75:5;78:4;79:4;80:2;81:1;82:1;84:1`).
  **hepatitis × U-MCCD** goes from 12 to a spread of 1–18.
* Even the calmest cell moves: **stamps × SUN-MCCD** under dropout has Jaccard
  0.990 and a flip rate of 0.0016, but **6.8%** of points still flip at least
  once.

Dropout (1% of rows) is uniformly gentler than jitter for the RK-based
detectors, and roughly comparable for the NND-based ones.

### Radii

The radius tap records the per-point radius vector in original row order.

**Consistency check that the tap is doing what it claims**: U-MCCD and SU-MCCD
report *identical* radius statistics on the same dataset, as do UN-MCCD and
SUN-MCCD — correct, because each pair shares one radius routine
(`ccd.Kest.edge.quantile` for the RK pair, `nnccd.radi` for the NND pair) on the
same input. Radius variability is a property of the **(dataset, radius family)**,
not of the detector.

| Dataset | RK family (U-, SU-MCCD) mean CV | NND family (UN-, SUN-MCCD) mean CV |
|---|---|---|
| | jitter / dropout | jitter / dropout |
| hepatitis | 0.017 / 0.016 | 0.173 / 1.422 |
| glass | 0.073 / 0.026 | 0.598 / 0.919 |
| stamps | 0.073 / 0.026 | 0.527 / 0.981 |
| pima | 0.055 / 0.023 | 0.273 / 0.588 |

Columns: `radius_mean`, `radius_mean_pointsd`, `radius_mean_pointcv`,
`radius_max_pointcv`, `radius_frac_points_changed`.

**The NND-family CV is inflated by exact zeros and must not be read as
"NND radii are 20–50× noisier".** With `scores = FALSE`, `nnccd.radi`'s
`j == low.num` early-exit branch sets `R[i] = 0`, and the zero-avoidance
fallback (`R[i] = sort(ddx[i,])[2]`) is applied only on the `scores = TRUE`
path. Measured on the unperturbed cells:

| Dataset | RK zero radii | NND zero radii |
|---|---|---|
| hepatitis | 0 / 74 (0.0%) | 26 / 74 (**35.1%**) |
| glass | 0 / 213 (0.0%) | 66 / 213 (**31.0%**) |
| stamps | 0 / 340 (0.0%) | 103 / 340 (**30.3%**) |
| pima | 0 / 555 (0.0%) | 83 / 555 (**15.0%**) |

A point whose radius is 0 in 49 of 50 replicates and non-zero in one has
CV = √49 = 7.071 — which is exactly the `radius_max_pointcv` value recorded for
every NND-family cell. So the NND radius CV is dominated by points toggling
across the zero boundary, not by broad radius jitter. The RK radii, which are
never zero here, have CV ≤ 0.073 throughout — i.e. **radii themselves are
stable; it is the connectivity built on them that is not.** That is precisely
R3.1's point, and it localises the instability in the graph step rather than the
radius step.

---

## What this supports in the manuscript

1. The pipeline is **deterministic** — reruns reproduce published numbers
   bit-for-bit. Proved two ways (outcome + RNG stream), not asserted.
2. **R3.1's discontinuity concern is real and is now quantified** rather than
   waved away. A 1%-of-one-SD jitter can move U-MCCD's flagged set to Jaccard
   0.14 against its own unperturbed output.
3. **SUN-MCCD is the most perturbation-stable of the four detectors** on the
   label-level and metric-level measures, which strengthens the headline
   recommendation on a dimension the paper had not previously measured.
4. The instability sits in **cluster coverage / graph connectivity**, not in
   radius estimation — consistent with the over-flagging failure mode already
   reported for high dimension.

## Caveats to carry into the response letter

* SU-/SUN-MCCD ARI ≈ 1 is **degenerate** (single macro-cluster at
  `S_min = 0.0625` on these four sets), not a stability result. Quote Jaccard
  and flip rate for those two.
* Four datasets, n ≤ 555. Not extended to the large sets on cost grounds.
* One perturbation magnitude only (σ = 0.01, drop 1%). No sensitivity to the
  perturbation size itself.
* `radius_max_pointcv = 7.071` is a saturation artefact of the zero-radius
  toggle described above, not a meaningful dispersion figure.

## Provenance note

The run was interrupted once when the external USB-NVMe enclosure holding `G:`
dropped off the bus (Windows `disk` event 51 and `Ntfs` event 50, "Delayed Write
Failed", 2026-08-09 22:08). The pima jitter arm lost cells 39–50 for three
methods. Because every cell is appended to the CSV as soon as it is computed,
the loss was bounded and the missing cells were recomputed after the volume
re-enumerated. Post-recovery integrity checks: 1,616 rows parse with 30 columns,
0 blank timestamps, 0 empty `radii_vec`, `radii_aligned` TRUE on every row,
`status == "ok"` on every row.
