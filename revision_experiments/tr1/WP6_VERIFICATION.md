# WP6 verification (Opus, 2026-09-06) — verdict: PASS WITH REQUIRED CHANGES

Reviewed commits caecfe8 (protocol) and bd0681d (91_wp6_runtime.R, 91b_wp6_runtime_py.py).
Confirmed PASS: the four wp0 MCCD wrappers call the detector once (no double build; the harness's
generic wrappers that build twice are never called); BLAS is the single-threaded reference build;
no MCCD code path is parallel; iForest nthreads = 1; torch threads = 1; table coverage for every
cell (d = 5, 10: 5000 rows; d = 50, 100: 1000 rows); get_simul(n = nrow(X)) everywhere; the eight
competitors line-by-line identical to 81_wp4_baselines.py except that only fit() is timed; per-cell
append with keys-only done files; the pandas NaN fix; idle check refuses correctly; complexity
targets match the manuscript (tex 921-927, 947: O(n^3 log n) for U/SU-MCCD, O(n^3) for UN/SUN-MCCD;
space O(n^2) at 933-943).

## Required changes (all must land before the run)

**A1 (blocker) — the RK table load is inside the driver's timed region.** umccd_method/sumccd_method
call get_simul() (wp0_mccd_methods.R:276, 293), which load()s from disk every call with no cache
(harness.R:230-236), before the wrapper's own t0 but inside the driver's proc.time() at
91_wp6_runtime.R:406-410. Measured: U-MCCD outer 3.96 s vs wrapper 1.28 s at n = 250 (2.68 s of
load); the smoke row shows time_s 2.95 with table_load_s 2.66. A constant 2.7 s offset turns the
log-log slope from ~1.0-1.6 into 0.36. Fix: memoize get_simul in the driver, inserted right after
91_wp6_runtime.R:81:

```r
local({
  .orig <- get_simul
  .cache <- new.env(parent = emptyenv())
  get_simul <<- function(variant = c("RK", "NN"), d, quant = NULL, n = NULL) {
    variant <- match.arg(variant)
    key <- paste(variant, d, quant, sep = "|")
    if (!exists(key, envir = .cache, inherits = FALSE)) {
      rm(list = ls(.cache, all.names = TRUE), envir = .cache)  # one entry: ~763 MB each
      assign(key, .orig(variant, d, quant, n = NULL), envir = .cache)
    }
    tab <- get(key, envir = .cache, inherits = FALSE)
    check_simul_extent(tab$simul, variant, d, n, tab$file)   # keep the loud short-table abort
    tab
  }
})
```
(Adapt the field names to what get_simul actually returns; the point is one load per
(variant, d, quant) outside every timed region, with the extent check still run per call.)

**A2 (blocker)** — move `tls <- table_load_seconds(...)` (91:402) out of the rep loop to just inside
`for (method in R_METHOD_ORDER) {` (after :393) and reuse it on all ten rows; that call both
measures the cold load and warms the cache.

**B1 (blocker) — memory is an absolute peak, not a delta.** gc(reset=TRUE) then sum(gc()[,6])
reports baseline + high-water of new allocation; the 862.7 MB at n = 100 is the wrapper's own
in-region table load (every RK file carries an unused 763 MB Kest.m matrix; the detector reads only
$r and $quan). Fix at 91:404 / :240 / :412:
```r
gc_before <- gc(reset = TRUE)
...
mem_delta_mb <- sum(gc_after[, 6]) - sum(gc_before[, 2])
```
Add mem_delta_mb to ROW_COLS (:261-262); keep mem_peak_mb for provenance; report the delta.
Python tracemalloc is already a delta (started after X and model construction); disclose that the
R and Python memory columns are deltas of different allocators and not comparable across languages.

**B2** — add a memory slope: lm(log(mem_delta_mb) ~ log(n)) at d = 10 on all reps, per MCCD
method, into 91_wp6_slope.csv (91:329-345), compared against the manuscript's space claim O(n^2).

**C (major) — hdbscan is multi-threaded by default** (core_dist_n_jobs = 4, joblib, not OMP).
91b_wp6_runtime_py.py:293 → `hdbscan.HDBSCAN(min_cluster_size=5, min_samples=None, core_dist_n_jobs=1)`.
Amend protocol §4's sentence claiming hdbscan has no separate thread knob; disclose the deviation
from 81_wp4_baselines.py (timing only, scores unchanged).

**D (blocker) — the three d-sweep exports collide on one filename.** CELLS$cell_value (91:173) is
"500" for cells 6-8 and export_path_for tags all three "d", so only d = 5 is ever written and the
Python d-sweep collapses to one point. Fix 91:173:
```r
cell_value = c("100", "250", "500", "1000", "2000", "5", "50", "100"),
```
and make export_rep1 (91:270-277) write cell 3 (n = 500, d = 10) under BOTH n_500_rep1.csv and
d_10_rep1.csv so the Python d-sweep has four points. results/tr1/wp6/data/ is empty now; nothing
stale to clear.

**E (major) — first-call warm-up inside the Python timed region.** `from sklearn.neighbors import
NearestNeighbors` (91b:227) and `from scipy import sparse` (91b:242) are inside functions called
from the timed region; lazy init shows as ECOD 0.51 s vs COPOD 0.002 s, and LUNAR seed1 5.03 s /
65 MB vs seed2 1.93 s / 0.19 MB. Hoist both imports to module scope and run one untimed warm-up fit
per (method, variant) at the smallest cell before the measured reps.

**Summariser** — add mean_time_s next to median (CV = sd/mean needs its own centre); add the
log n-folded fit for U/SU-MCCD (fit log(time_s/log(n)) ~ log(n), compare to 3); summarise the
Python raw CSV too and PRE-DECLARE the collapse rules: MutualKNN/SNN reported at k = 10 (median over
reps) with the full k range in the supplement; DIF/LUNAR median over all 50 seed x rep fits.

## Disclosures for the caption
DBSCAN eps uses the true contamination; MST is run at cont = 0.05 (generator's rate) not the
registry's 0.02; ODIN's k = round(sqrt(n)) grows with n so its slope is not fixed-parameter; the
idle checks fail open if tasklist is unavailable and do not see non-R/Python load.

## Cost ruling
The WP8 single-population slowdown (125 s at n = 200, d = 3) does not apply: WP6 data are two
uniform clusters 3 apart with a clean connectivity break, d >= 5. Estimate n = 2000 MCCD at 40-150 s
per rep; R side 1-2.5 h, Python side 1-2.5 h (DIF/LUNAR at n = 2000 dominate). Total 2.5-5 h idle.
No rep cap needed. Run --grid=n and --grid=d as separate invocations.

## Launch sequence (after all changes; only when the WP8 grids have finished; never --force)
```
Rscript revision_experiments/tr1/91_wp6_runtime.R --smoke           # expect U-MCCD time_s ~0.3 s, mem_delta low hundreds MB
.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --smoke
Rscript revision_experiments/tr1/91_wp6_runtime.R --idle-check
.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --idle-check
Rscript revision_experiments/tr1/91_wp6_runtime.R --grid=n --reps=10
Rscript revision_experiments/tr1/91_wp6_runtime.R --grid=d --reps=10
ls revision_experiments/results/tr1/wp6/data/                      # 9 export files
.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --idle-check
.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --reps 10
Rscript revision_experiments/tr1/91_wp6_runtime.R --summarize
```
Sanity checks before the numbers go near the manuscript: U-MCCD at n = 100 has table_load_s ~2.7
and time_s ~0.3 (not the reverse); mem_delta_mb for U-MCCD and SUN-MCCD within one order of
magnitude at equal n; LUNAR rep 1 no longer an outlier; the d-sweep CSV carries four distinct d.
