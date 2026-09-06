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

> **Appended note, 2026-09-05 (WP5 R6, verifier pass):** "copied byte-for-byte"
> above overstates it. `handle_reused()` loads each source CSV with
> `pandas.read_csv` and writes it back out with `DataFrame.to_csv` — the
> feature and label *values* are unchanged, but the bytes on disk are not
> guaranteed identical (pandas' float formatting/serialization can differ
> from whatever produced the original file, and line endings are not
> verified either). The correct description is **value-identical,
> re-serialized by pandas**, not byte-for-byte. Nothing about the gate logic
> or the (n, d, n_outliers) verification changes; only this sentence's
> precision does.

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

> **Appended note, 2026-09-05 (deviation forced by an implementation
> failure, discovered during the run).** letter (d=32) does not, in fact,
> produce a caveated U-MCCD/SU-MCCD row as declared above. Both cells error
> deterministically: `Error in integrate(integrand, 0, acos(t/2)) : a limit
> is NA or NaN`, raised from `Kest.f.edge()`'s Ripley's-K edge-correction
> integral, reached via
> `RUMCCD_outlier -> RKCCD_correct_quant -> rccd.clustering_correct_quantile
> -> ccd.Kest.edge.quantile -> Kest.f.edge -> sapply/lapply -> integrate`
> (full traceback captured, not merely observed as a crash). mnist (d=100)
> does not hit this failure — both of its RK-based rows complete and carry
> the declared caveat. The mechanism is consistent with §5's own explanation
> (51.276% zero quantiles at d=32 is enough to produce a degenerate geometry
> feeding the edge-correction integral a `t` outside `acos`'s domain), but
> it manifests as a hard error rather than a caveated number at this
> dimension specifically. This script does not attempt-then-catch-and-round
> the failure into a caveated result — the run's own `84_wp5_highd.R`
> `tryCatch` records it faithfully as `status = "error"`. No fix is applied
> here: `Kest.f.edge()` lives under `R/` and `RUMCCD_outlier()` under
> `methods/`, both out of scope for this WP. See `WP5_FINDINGS.md` for the
> reported consequence — the d=32 RK-family row in R1.8's literal band is
> `error`, not a number, revising the "one degradation curve" framing above:
> the actual d=32/100/166/274 status sequence is error/caveated-ok/n-a/n-a
> for U-MCCD and SU-MCCD, not caveated-ok/caveated-ok/n-a/n-a.

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

> **Appended note, 2026-09-05 (WP5 R6, verifier pass):** "`84b_wp5_metrics.R`
> calls the modified script" above is imprecise — `84b_wp5_metrics.R` does
> not invoke `81_wp4_baselines.py` itself; it only *reads* the score files
> `81` writes under `results/tr1/wp5/scores/`. `81` must be run separately,
> to completion (`--status` reporting 0 remaining), before `84b` is run —
> `84b` hard-errors on a missing score file rather than launching `81` to
> produce it.

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

> **Appended note, 2026-09-05 (WP5 R6, verifier pass):** "all 17 methods" in
> the `wp5_metrics_main.csv` row above (and the identical phrase at
> `84b_wp5_metrics.R:35`) undercounts the row count. The 4 proposed + 5
> original baselines + 8 WP4 competitors are 17 *methods*, but three of the
> eight competitors (HDBSCAN/GLOSH, OPTICS, mutual-kNN) are each reported
> under **two** rows per data set — a T1-thresholded row and a separate
> native-label row (`HDBSCAN-noise`, `OPTICS-noise`, `MutualKNN-m0`) — so the
> merged table carries **20 rows per data set**, not 17. Nothing about which
> cells are computed changes; only the row count this note corrects.

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

> **Appended note, 2026-09-05 (WP5 R3, verifier pass):** a precision that
> matters for how this WP's band-coverage claim should be read. R1.8 names a
> literal band, "20 < d < 50" (the gap between the current Section 6 ceiling
> of d=21 and the paper's claimed degradation onset at d>=50). Of this WP's
> four data sets, **only letter (d=32) falls inside that literal band.**
> mnist (d=100), musk (d=166), and arrhythmia (d=274) all sit above it and
> extend the degradation curve past d=50, not within it. So this WP answers
> R1.8 with **one point inside the named band plus a three-point degradation
> trajectory beyond it**, not with four points spanning the band — the
> manuscript text must not imply four in-band observations. The protocol
> pre-commits, now, to reporting whatever letter's cell actually shows,
> whichever direction it goes, rather than choosing after the run which
> framing to use.

> **Appended note, 2026-09-05 (WP5 follow-up, deviation decided POST-HOC —
> after seeing the S5 crash, not before it).** The letter (d=32) RK
> integration error recorded in S5's appended note above is not a fluke of
> `Kest.f.edge()`'s edge-correction integral in general — it has a specific,
> confirmed cause in the input data. `letter.csv` (as written by
> `84a_wp5_fetch_convert.py` from the ADBench mirror, before any change
> described here) contains **exactly 2 duplicated feature rows** —
> `sum(duplicated(X)) == 2` — i.e. 2 pairs of exact-duplicate points (row
> indices 373/476 and 94/384 in the written CSV), for 4 duplicate-involved
> rows total. **Both duplicated pairs are labelled regular** (repo
> convention: `label == 1`), not outlier — no outlier point is involved.
> `letter.csv`'s minimum non-zero pairwise Euclidean distance is 0.2248303;
> the duplicate pairs sit at distance exactly 0. mnist, musk, and arrhythmia
> were checked the same way and contain **zero** duplicated feature rows —
> they do not exhibit this crash and are not affected by anything in this
> note.
>
> Mechanism, confirmed (not merely inferred from the traceback already
> logged under S5): a zero-distance pair hands `Kest.f.edge()` a candidate
> covering-ball radius of exactly 0. That zero radius makes the ratio `M/sc`
> evaluate to `NaN` (0/0-shaped), and `NaN` propagates into the integration
> bound `integrate(integrand, 0, acos(t/2))`, which is exactly the "a limit
> is NA or NaN" error S5 recorded. This is a strictly different mechanism
> from the acos-domain concern raised when S5's crash was first discovered:
> no traced call had `M/sc > 2` — the failure is a zero radius, not an
> out-of-domain acos argument.
>
> **This paper's own real-data collection is duplicate-free by convention
> already.** The 16 data sets used in Section 6 are drawn from the DAMI
> "withoutdupl" file family (ODDS/ELKI), which deduplicates by construction.
> letter is the one WP5 data set sourced from ADBench instead, and ADBench's
> mirror does not deduplicate. Applying the collection's own existing
> duplicate-free convention to letter — drop exact duplicate feature rows,
> keep the first occurrence of each duplicate pair — is therefore not a new
> rule invented to route around a crash; it is the rule every other data set
> in this study already satisfies, applied to the one set that happens not
> to.
>
> **This is disclosed as a POST-HOC deviation, not a pre-declared one.** The
> decision to deduplicate letter was made after seeing the `integrate()`
> crash, specifically to explain and resolve it — S1-S9 above did not
> anticipate or declare deduplication as a preprocessing step for any WP5
> data set. The remedy applied is: remove the 2 duplicate rows (keep first
> occurrence), giving `n = 1598`, `n_outliers = 100` unchanged (both removed
> rows were regular points), contamination rising marginally from 6.25% to
> 100/1598 = 6.258%. A one-off diagnostic run of `U-MCCD` (registry wrapper,
> study settings, `min.cls` not applicable to U-MCCD) on the deduplicated
> n=1598 matrix completed without error in ~45s: TPR=0.320, TNR=0.899,
> BA=0.610, F2=0.274 — confirming the fix resolves the crash before the full
> rerun (§ below) is launched.
>
> **Both the original-run rows (n=1600, with the crash) and the
> deduplicated rerun (n=1598) are kept on disk.** The original letter rows
> from `wp5_highd_results.csv`/`wp5_highd_done.csv`/`fit_log.csv` and the
> original `results/tr1/wp5/scores/letter_*` files are preserved under
> `results/tr1/wp5/letter_with_duplicates/` before the live files are
> stripped of their letter rows and letter is rerun end to end
> (`84_wp5_highd.R`, `81_wp4_baselines.py`, `84b_wp5_metrics.R`) against the
> deduplicated data. `84a_wp5_fetch_convert.py` is amended to perform this
> deduplication as a declared step for letter only, recording `n_raw` and
> `n_duplicates_removed` in `manifest.csv` for every data set (0 for the
> three sets unaffected). See `WP5_FINDINGS.md` for the full before/after
> comparison across all 18 non-RK methods.
