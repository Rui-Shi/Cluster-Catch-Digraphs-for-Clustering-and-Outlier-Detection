# WP5 protocol — real data above d = 21

**Declared 2026-09-05, before any download or run.** Answers AE.4, R1.8, R3.8,
R5.6. R1.8 names the gap explicitly: every real data set in Section 6 tops
out at d = 21, but the paper claims degradation at d >= 50 using only
synthetic evidence. This protocol fixes the data, the preprocessing, the
methods, the n/a rule, and what gets reported, before any of it is run. Any
deviation forced by an implementation failure is recorded as an appended,
dated note rather than an edit to the declarations below it — the same
discipline `WP4_PROTOCOL.md` uses.

Factual base: `revision_experiments/tr1/WP5_INVENTORY.md` (read-only
reconnaissance, written before this protocol). Nothing in that file is
re-litigated here except where this protocol has to make a call the
inventory left open (musk's n, and the mnist subsample seed) — both are
recorded in §6.

## 1. Data sets

| Dataset | n | d | n_outliers | contamination |
|---|---|---|---|---|
| letter | 1600 | 32 | 100 | 6.25% |
| mnist | 1000 (subsample) | 100 | measured at fetch time (~9.2% source rate carried into the subsample) | ~9% |
| musk | 1000 (subsample, reused) | 166 | 32 | 3.2% |
| arrhythmia | 452 | 274 | 66 | 14.6% |

### Sources

- **letter**: `https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/20_letter.npz`
  — ADBench mirror of the ODDS "Letter Recognition" outlier set. Verified
  reachable 2026-09-05 (HTTP 200, 23,946 bytes). ADBench repo commit
  `f7fed68dea6901fe9a81ed51b251b1d2456790cd` (last touch of this path, per the
  GitHub commits API, 2023-07-17) — the same mirror and vintage already used
  for Musk/Speech/InternetAds (`revision_experiments/results/datasets_csv/manifest.csv`).
  Format: `.npz` with arrays `X` (float, n x d) and `y` (int, source polarity
  **1 = outlier, 0 = regular** — confirmed directly against the sibling
  `musk.npz` already in this repo: `(X.shape, y.shape) = ((3062,166),(3062,))`,
  `y` unique counts `{0: 2965, 1: 97}`, matching the manifest's 97-outlier
  figure exactly).
- **mnist**: `https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/24_mnist.npz`
  — ADBench mirror of the ODDS MNIST outlier set. Verified reachable
  2026-09-05 (HTTP 200, 455,359 bytes). Same repo, same commit, same `.npz`
  layout and label polarity as letter. The literature figure for this set is
  n=7603, d=100, 700 outliers (9.21%); `84a_wp5_fetch_convert.py` prints the
  values it actually reads and treats a d mismatch as fatal (d=100 is the
  point of including this set) but does not hard-gate on the literature n /
  outlier count, since no file in this repo pre-declares them authoritatively
  — see §6 for why this is a judgement call, not a discovered fact.
- **musk**: reused as-is from `revision_experiments/results/datasets_csv/Musk_sub1000.csv`
  (TR2's n=1000 contamination-preserving subsample of the full n=3062 Musk
  set, seed 20260716, single draw — `revision_experiments/tr2/07_wp5_subsample_ccd.R`
  header and `FINDINGS.md:105`). Not redrawn. See §6 for why the full n=3062
  set is not used here even though this paper's own detectors could likely
  afford it at that n (WP5_INVENTORY.md §5).
- **arrhythmia**: reused as-is from `revision_experiments/results/datasets_csv/Arrhythmia.csv`
  (full n=452, community ODDS mirror, already exported and gated against
  (n,d,n_outliers) by `tr2/02_load_data.R`).

## 2. Subsampling rule

**mnist only.** letter and arrhythmia are used in full; musk is already a
subsample and is not redrawn (§1, §6).

Contamination-preserving, single draw, fixed seed, no replacement:

1. Let `n_out_full`, `n_reg_full` be the outlier/regular counts in the full
   downloaded set (source polarity, before the label flip in §3) and
   `cont = n_out_full / (n_out_full + n_reg_full)`.
2. `n_out_sub = round(cont * 1000)`, `n_reg_sub = 1000 - n_out_sub`.
3. Draw `n_out_sub` indices from the outlier pool and `n_reg_sub` indices from
   the regular pool, independently, without replacement, via
   `numpy.random.default_rng(seed)`.
4. Concatenate; row order is not otherwise meaningful (`evaluate()` in
   `shared/harness.R` jointly reorders `(Y, score)` before scoring, so no
   downstream step depends on regulars-first ordering).

**Seed: 20260905** (the date this protocol is declared) — chosen fresh for
this draw and recorded here before the fetch script exists to run it. This is
a different seed from musk's inherited 20260716; the two draws are unrelated
and neither is retroactively reconciled to match the other.

## 3. Preprocessing

**Label convention in every CSV under `results/tr1/wp5/data/`: 1 = regular,
0 = outlier** (this repo's convention throughout, opposite of the ADBench/
ODDS source files). letter and mnist are flipped from source polarity
(`y_repo = 1 - y_source`) by `84a_wp5_fetch_convert.py`, the same idiom
`RealData_Collection.R` and `tr2/02_load_data.R` use for their own raw
`.mat`/`.npz` sources (`ifelse(y == 1, 0, 1)`).

**letter and mnist get the robust MADN standardization the study uses for
its larger high-d sets** — the exact function is `scale_R_safe()`, defined at
`revision_experiments/tr2/02_load_data.R:142-153`:

```r
scale_R_safe <- function(x) {
  M <- median(x); madn <- mad(x)
  if (madn > 0) return((x - M) / madn)
  sdx <- sd(x)
  if (sdx > 0) return((x - M) / sdx)
  return(x - M)  # constant column -> all zeros, no NaN
}
```

applied column-wise: `RealData_Collection.R`'s own `scale_R()` is
`(x - median(x)) / mad(x)` with no fallback, which is exactly what
`scale_R_safe()` reduces to whenever `mad(x) > 0`; the SD/constant fallback
only fires on a column where the robust scale collapses to zero, which is
untested territory for the paper's original low-d sets but did fire on Musk
(3/166 columns) and Arrhythmia (149/274 columns) per the existing manifest.
`84a_wp5_fetch_convert.py` reimplements this identically in `numpy`
(median / `1.4826 * median(|x - median(x)|)`, matching R's `mad()` default
constant exactly, then the same two-level fallback) and reports how many
columns of letter and mnist hit each branch — declared now, not decided after
seeing the counts.

**musk and arrhythmia keep whatever the existing manifest records for
them** — both are copied byte-for-byte (feature values unchanged) from
`results/datasets_csv/{Musk_sub1000,Arrhythmia}.csv` into
`results/tr1/wp5/data/`, i.e. robust median/MADN per column with the same
SD/constant fallback, already applied before the CSVs were written. Nothing
in `84a_wp5_fetch_convert.py` re-scales them; it only re-verifies (n, d,
n_outliers) against §1's table and renames the file.

## 4. Methods and settings

### 4.1 The four proposed detectors

`U-MCCD`, `SU-MCCD`, `UN-MCCD`, `SUN-MCCD` via the registry wrappers in
`tr1/wp0_mccd_methods.R` (`umccd_method`, `sumccd_method`, `unmccd_method`,
`sunmccd_method`), unchanged from every other WP that has used them.

- **S_min = 0.05** passed as `min.cls = 0.05` to `SU-MCCD` and `SUN-MCCD`
  only (`U-MCCD`/`UN-MCCD` take no `min.cls` argument — they are the
  uniform-coverage reference constructions, no shape-adaptive step to gate).
  This is the single constant CLAUDE.md records for "every set and for the
  simulations" — no WP5-specific value.
- **alpha** resolved by the paper's own three resolvers in `shared/harness.R`
  (`rk_quant_label_paper`, `nn_quant_label_paper_UN`,
  `nn_quant_label_paper_SUN`), called with no override. At all four of
  letter/mnist/musk/arrhythmia (d = 32, 100, 166, 274, all >= 20) every
  resolver returns `"999"` (alpha = 0.1%) — there is no dimension-dependent
  alpha schedule left to exercise in this WP; every running cell is at the
  same significance level. This is worth stating plainly in the manuscript
  text next to the table, since a reader could otherwise assume the four
  data sets probe different alpha regimes the way the d = 2..21 real-data
  table does.
- Threshold 0.5 for all four (`REAL_DATA_THRESHOLDS`, unchanged).

### 4.2 The five original baselines

`LOF`, `DBSCAN`, `MST`, `ODIN`, `iForest`, via `METHOD_REGISTRY` in
`shared/harness.R`, defaults unchanged (`LOF` MinPts 11:30 max; `DBSCAN`
k=4, oracle contamination; `MST` cont=0.02, thresh=1.2; `ODIN` default
k=round(sqrt(n)); `iForest` ntrees=1000, sample_size=min(256,n), seed=1).

### 4.3 The eight WP4 competitors

`ECOD`, `COPOD`, `DIF` (5 seeds), `LUNAR` (5 seeds), `HDBSCAN`/`GLOSH`,
`OPTICS`, mutual-kNN (k in {5,10,15,20,30}, oracle-k on F2), `SNN`
(k in {5,10,15,20,30}, oracle-k on F2) — identical settings, thresholding
rule (T1: pyod-style `contamination = 0.1`, strict `>` comparison), and
seeding convention to `tr1/WP4_PROTOCOL.md`, run via the same
`81_wp4_baselines.py` driver against the WP5 data/scores folders (the only
change to that script — see §7 — is where it reads data from and writes
scores to; every method, hyperparameter, and threshold rule is untouched).
No T2 (oracle-contamination) sensitivity is computed for WP5 unless the
main-table T1 result raises a question that needs it, to keep this WP inside
its 2-day budget; T1 is the reported regime, exactly as it is the reported
regime in WP4.

## 5. The n/a rule

**U-MCCD and SU-MCCD (the RK-based detectors) are reported `n/a` at musk
(d=166) and arrhythmia (d=274)**, with the reason
`"RK envelope 100% zero quantiles"` recorded explicitly in the output row —
not skipped, not omitted. Mechanism: no production RK quantile table exists
at d=166 or d=274 in `R/RK-test_quantile/` at all (confirmed by directory
listing), and a niter=20 probe table built specifically to measure this
(`revision_experiments/results/tr2/probes/RK-test-simul_{166,274}d_999%.RData`)
shows **100.0% of quantile cells are zero** at both dimensions
(WP5_INVENTORY.md §3). A covering ball radius drawn from an all-zero
quantile envelope collapses every point to isolation by construction, so
running the detector would not produce a meaningful score — it would produce
noise dressed as a result. `84_wp5_highd.R` checks for the RK table's
existence before attempting the call and writes the n/a row directly; it does
not attempt-then-catch a failure, because the failure mode here is not an
error, it is a well-understood structural fact about the RK covering-sphere
construction in high dimensions (the same mechanism R1.2 asks the paper to
explain).

**At letter (d=32) and mnist (d=100), U-MCCD and SU-MCCD are run, not
n/a'd** — RK tables exist at both dimensions — **but both rows carry a
degeneracy caveat** in the output: 51.276% of the d=32 RK envelope and 91.4%
of the d=100 envelope are zero quantiles (measured directly on the
production tables, WP5_INVENTORY.md §3, matching HANDOFF_FROM_TR2.md's
independently-derived figures exactly). This is the same mechanism as the
n/a rows, one and two dimensionality tiers earlier and not yet total. The
manuscript should present d=32/100/166/274 as one degradation curve
(51% -> 91% -> 100% -> 100% zero), not as two running cells and two
unrelated n/a's.

`84_wp5_highd.R` also writes a small companion table,
`results/tr1/wp5/rk_degeneracy.csv`, with one row per d in {32, 100, 166,
274}: the zero-quantile fraction, its source (production table vs. niter=20
probe table), and the table file path — this is the number the n/a rows and
the caveated rows both cite, and it is what answers R1.2 with a measurement
rather than an assertion.

## 6. Decisions this protocol makes that the revision plan and inventory left open

- **Musk runs at n=1000, not the full n=3062.** The revision plan's WP5 table
  lists musk's n as "1000 (subsample)"; WP5_INVENTORY.md §6 flags that this
  paper's own detectors are cheap enough (WP0's measured runtime, ~n^2-ish,
  not TR2's ~n^4 UN-CCD cost) that full-data musk would very likely be
  affordable, and that a spliced n=3062 NN table already exists from TR2's
  own work. This protocol does **not** take that option: `get_simul("NN",
  166, quant = nn_quant_label_paper_SUN(166), n = 3062)` resolves the
  filename `NN-test-simul_166d_999%.RData` (the paper's own resolver, not the
  spliced file's non-standard name `..._n3062_spliced.RData`), and that file's
  extent is 1000 rows — `get_simul()`'s own extent check
  (`shared/harness.R:162-186`) would refuse it as "too short for the data"
  before any compute happened. Reaching for the spliced table by name would
  mean hand-picking a table outside the paper's own resolver logic for one
  data set only, which is exactly the kind of ad hoc exception this revision
  has spent effort eliminating elsewhere (`wp0_mccd_methods.R`'s header on the
  RK/NN resolver unification). The n=1000 subsample the resolver already
  supports cleanly is used instead, and the full-data run is named here as a
  possible future addition, not attempted.
- **mnist's subsample seed (20260905) is newly chosen, not inherited from
  musk's 20260716.** The two draws are independent judgement calls at
  different times for different data sets; nothing requires them to match and
  reusing musk's seed for an unrelated draw would suggest a coupling that
  does not exist.
- **mnist's full-set n/d/outlier counts are not hard-gated**, unlike every
  other data set in this WP (see §1). No file in this repo pre-declares an
  authoritative ODDS/ADBench mnist outlier-set (n, n_outliers) the way
  `CCD_OutlierDetection_Neurocomputing.tex`'s Table `tab:Real_Data` does for
  the sixteen Section 6 sets — the only number this WP truly depends on is
  d=100, which is checked and is fatal if wrong. `84a_wp5_fetch_convert.py`
  prints the values it reads (expected, from the literature: n=7603,
  n_outliers=700, 9.21%) so a large disagreement is visible, but a data set
  whose sole purpose is "the d=100 anchor point" does not need its exact
  outlier count locked in advance the way a set feeding a headline comparison
  table would.

## 7. The one permitted edit below script 84

`tr1/81_wp4_baselines.py` gains two optional CLI arguments, `--data-dir` and
`--out-dir`, whose *defaults reproduce current behaviour exactly*
(`results/tr1/wp4/data` and `results/tr1/wp4` respectively, computed the same
way the hardcoded constants were). No method, hyperparameter, seed, threshold
rule, or file-naming convention inside the script changes. `84_wp5_highd.R`'s
companion, `84b_wp5_metrics.R`, calls the modified script (via a documented
invocation, not by re-implementing its logic) with
`--data-dir results/tr1/wp5/data --out-dir results/tr1/wp5` so the eight WP4
competitors run against the WP5 data folder and write into the WP5 results
tree instead of WP4's. The exact diff is recorded in the WP5 hand-off report,
not duplicated here.

## 8. Outputs

```
results/tr1/wp5/data/{letter,mnist,musk,arrhythmia}.csv     features (V1..Vd) + label (1=regular, 0=outlier)
results/tr1/wp5/data/manifest.csv                           n, d, n_outliers, contamination, preprocessing, source, subsample_seed
results/tr1/wp5/wp5_highd_results.csv                       one row per (dataset, method[, min_cls]) -- MCCD (4) + baselines (5); n/a rows explicit
results/tr1/wp5/rk_degeneracy.csv                           zero-quantile fraction per d in {32,100,166,274}, source, table path
results/tr1/wp5/scores/                                     raw WP4-competitor scores (81's output, redirected via --out-dir)
results/tr1/wp5/fit_log.csv, fit_errors.log, versions.txt   81's own per-cell log, redirected via --out-dir
results/tr1/wp5/wp5_metrics_main.csv                        84b's merged per-(dataset,method) TPR/TNR/BA/F2, all 17 methods
results/tr1/wp5/WP5_FINDINGS.md                             84b's narrative summary
results/tr1/wp5/smoke/                                      --smoke outputs for 84 and 81 (proof-of-wiring only, not results)
```

`results/` is gitignored in this nested repo; only this protocol file and the
scripts that implement it are committed.

## 9. What is reported

- **Per-set table**: one row per (letter, mnist, musk, arrhythmia) x
  (4 proposed + 5 baseline + 8 competitor) cell, TPR/TNR/BA/F2, with the
  U-MCCD/SU-MCCD n/a rows at musk/arrhythmia carrying the reason string from
  §5 in place of numbers.
- **The RK-degeneracy fraction per d** (§5's `rk_degeneracy.csv`) as the
  mechanism behind the n/a rows and the caveated letter/mnist rows — this is
  the evidence R1.2 asks for, extended two dimensionality tiers past the
  91.4%-at-d=100 figure already measured for WP2.
- Whatever the R3.3 (mutual reciprocity vs. mutual catch) and R1.3 (density
  clustering with varying-density support) comparisons already computed for
  WP4 look like on these four higher-d sets, following the same structure as
  `82_wp4_metrics.R`'s sections 6-7, since the same eight competitors are
  present and the question does not stop mattering above d=21.
