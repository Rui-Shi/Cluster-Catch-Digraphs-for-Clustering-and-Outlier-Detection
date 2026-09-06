# WP6 protocol — runtime and memory

**Declared 2026-09-05, before any timing cell is run.** Answers AE.5,
R3.10's runtime/memory/scalability half, R5.5, and R5.6's scalability
angle. R3.10 and R5.5 both name the same complaint: the manuscript proves
time complexity of `O(n^3)` or `O(n^3 log n)` (space `O(n^2)`) but reports
no measured runtime, memory, or scalability. This file fixes the grid, the
generator, the methods, what is timed, the memory method, and the
complexity-reconciliation plan in advance. Nothing here changes after a
result is seen; a deviation forced by an implementation failure is recorded
as an appended, dated note, not a silent edit.

Written by the WP6 author sub-session. **Not run for real by that
session** — an idle machine is required and four background R jobs were
running at write time. Only a tiny smoke (`--smoke`, one cell, one rep,
seconds of compute) was executed, to prove the scripts parse and work
end to end. The real grid is the verifier's / a later session's job, on an
idle machine, per the idle-check below.

## 0. The claim being reconciled

`CCD_OutlierDetection_Neurocomputing.tex`, Section "Time and Space
Complexity" (`sec:time_space_complexity`, lines 842–947), Table
`tab:space_time`:

| Algorithm | Time complexity |
|---|---|
| KS-CCDs | `O(n^3+n^2(d+\log n))` |
| RK-CCDs | `O(n^3(\log n+N)+n^2 d)` |
| U-MCCD | `O(n^3(\log n+N)+n^2(d+\log n))` |
| SU-MCCD | `O(n^3(\log n+N)+n^2(d+\log n))` |
| UN-MCCD | `O((N+d+\log n)n^2+n^3)` |
| SUN-MCCD | `O((N+d+\log n)n^2+n^3)` |

Two distinct claims hide inside reviewers' shorthand "`O(n^3)` or
`O(n^3\log n)`": the **RK-based** pair (U-MCCD, SU-MCCD) has a genuine
`n^3\log n` leading term (the `n^3(\log n+N)` factor); the **NND-based**
pair (UN-MCCD, SUN-MCCD) has leading term `n^3` alone — the `\log n` there
multiplies only the lower-order `n^2` term. All four are within-cluster
bounds and take the number of Monte Carlo replicates `N` (`niter`, fixed at
1000 by every wrapper) and `d` as constants for a fixed-`(N,d)` runtime
sweep, which is what the `n`-grid below measures.

**TR2's independent measurement is the reason this needs reconciling, not
just proving.** TR2 (Pattern Recognition revision, shared codebase) timed
its own **generic** `UNCCD-OOS`/`UNCCD-IOS` construction (a 22-core
parallel port of `nnccd.radi`, `revision_experiments/tr2/04_wp4_runtime.R`)
at `n ∈ {100,250,500,1000,2000}`, `d=10`, Gaussian two-cluster data, and
found empirical ratios implying **an exponent climbing 1.6 → 4.1**, i.e.
`~n^4` and steepening (`tr2/FINDINGS.md` §6/§10, `results/tr2/
wp4_runtime2_n.csv`) — steeper than either stated bound, on a *parallel*
implementation of the same nearest-neighbour radius search this paper's
UN-MCCD/SUN-MCCD also use (serially).

Against that, this paper's own four detectors (wired by
`tr1/wp0_mccd_methods.R`, single monolithic call per detector, no
parallelism) were already timed once, incidentally, on **real data** up to
`n=4819` during the WP0 reproduction gate (`results/tr1/wp0_gate*.csv`,
summarized in `tr1/WP5_INVENTORY.md` §5):

| Dataset | n | d | U-MCCD | SU-MCCD | UN-MCCD | SUN-MCCD |
|---|---|---|---|---|---|---|
| stamps | 340 | 9 | 2.4 s | 0.4 s | 0.5 s | 0.3 s |
| vowels | 1452 | 12 | 26.1 s | 18.0 s | 24.7 s | 16.6 s |
| waveform | 3443 | 21 | 381.3 s | 266–278 s | 424.5 s | 255–270 s |
| wilt | 4819 | 5 | 1130.6 s | 948–1015 s | 1051.1 s | 948–1015 s |

Fitting even this uncontrolled, mixed-`d`, four-point real-data series
(vowels → wilt, `n` ratio 3.32, time ratio ≈55.8×) gives a log-log slope
near **3.3–3.4**, not TR2's ~4, and not a constant either — consistent
with "close to the stated cubic bound but not exactly it," which is a
different and much less alarming answer than TR2's. **This is exactly why
a controlled synthetic sweep (fixed `d`, fixed generator, varying only
`n`) is needed rather than trusting either the real-data anecdote or the
TR2 cross-codebase analogy**: real data confounds `n` and `d`, and TR2's
number is for a different (parallelized, generic) construction, not for
this paper's own wired detectors.

**Pre-declared reporting rule:** the fitted slope (± SE) for each of the
four MCCD methods will be reported in Section 4.7 exactly as fitted,
compared explicitly against 3 (the stated NND-based bound) or against "3
with a log factor" (RK-based). **A fitted slope above the stated exponent
will be reported as such, in Section 4.7, not smoothed over or omitted.**
If UN-MCCD/SUN-MCCD's fitted slope over `n ∈ {100,...,2000}` also
steepens toward TR2's ~4 once `n` reaches 2000, that is the headline
finding and will be stated as an implementation-vs-bound gap, following
`REVISION_PLAN`'s own instruction: "say which, plainly, in Section 4.7."

## 1. Grid

Two sweeps, sharing the `(n=500, d=10)` cell (computed once, used by both
summaries):

- **Sweep A ("n")**: `n ∈ {100, 250, 500, 1000, 2000}` at `d = 10`.
- **Sweep B ("d")**: `d ∈ {5, 10, 50, 100}` at `n = 500`.

10 reps per (cell, method). Single thread (§4). Idle machine (§6).
Median and CV (`sd/mean`) reported per cell; the log-log slope fit (§7)
uses every rep, not just the medians.

### Quantile-table coverage (confirmed before declaring the grid final)

Every MCCD cell above resolves to a table that already exists and is long
enough, checked directly (`get_simul()`'s own `n` argument would abort
otherwise):

| d | RK label used | RK table rows | NN(UN) label | NN(UN) rows | NN(SUN) label | NN(SUN) rows |
|---|---|---|---|---|---|---|
| 5 | 99% | 5000 | 95% | 5000 | 95% | 5000 |
| 10 | 999% | 5000 | 99% | 5000 | 999% | 5000 |
| 50 | 999% | 1000 | 999% | 1000 | 999% | 1000 |
| 100 | 999% | 1000 | 999% | 1000 | 999% | 1000 |

(RK label from `rk_quant_label_paper(d)`; NN(UN) from
`nn_quant_label_paper_UN(d)`; NN(SUN) from `nn_quant_label_paper_SUN(d)`,
all three defined once in `shared/harness.R` §2.) Confirms the task
brief's two questions directly: **the `d=10` tables (both RK and NN) have
5000 rows, comfortably covering `n=2000`**, and **the `d=50`/`d=100`
tables have exactly 1000 rows, covering `n=500` with margin** — verified
by loading each of the eight `.RData` files and reading
`length(simul$average)`/`nrow(simul$quan[[1]])` directly (a few seconds of
disk I/O, not a timed computation; done once while writing this protocol).
No cell in the grid is `SKIPPED_NO_TABLE` by design.

## 2. Synthetic data generator

Reuses the uniform two-cluster geometry from `55_wp2c_simulation_arm.R`'s
`gen_uniform` (also the source WP8's `85_wp8_boundary_fp.R` copied from,
per that file's own header note) — **not** TR2's Gaussian generator in
`04_wp4_runtime.R`, per the task brief's explicit "uniform two-cluster"
instruction. `55_wp2c_simulation_arm.R` is below script 84 and is
read-only for this work package; the generator is reproduced here (as WP8
already did once), not sourced, and is fixed at declaration time:

```r
gen_uniform_wp6 <- function(seed, n, d, cont = 0.05) {
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- (mu1 + mu2) / 2
  n1 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)
  n2 <- n - n1 - n0                      # exact total for every n (TR2's
                                          # generalization of the reference
                                          # script's n=500-specific "-1")
  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) +
             matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) +
             matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis && sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1
    }
  }
  rownames(outlier) <- NULL
  X <- rbind(data1, data2, outlier)
  colnames(X) <- paste0("V", seq_len(d))
  Y <- c(rep(1, n1 + n2), rep(0, n0))    # 1 = regular, 0 = outlier
  list(X = X, Y = Y)
}
```

`rpoisball.unit` (`R/ccds/Kest.R`) is already in scope once
`shared/harness.R` is sourced (it sources
`methods/outlyingness_scores/RKCCD_OOS_IOS.R` /
`UNCCD_OOS_IOS.R`, which chain to it). 5% contamination throughout, fixed.
**Seeding**: `seed = 60000 + cell_index*100 + rep`, `cell_index` from the
fixed enumeration in §3 — one fresh draw per rep, same idiom as
`04_wp4_runtime.R` and `85_wp8_boundary_fp.R`. Data generation is **not**
timed (matches both prior harnesses' convention).

## 3. Cells and methods

Cell table (8 distinct `(n,d)` pairs; `(500,10)` serves both sweeps and is
computed once):

| cell_index | n | d | in sweep(s) |
|---|---|---|---|
| 1 | 100 | 10 | n |
| 2 | 250 | 10 | n |
| 3 | 500 | 10 | n, d |
| 4 | 1000 | 10 | n |
| 5 | 2000 | 10 | n |
| 6 | 500 | 5 | d |
| 7 | 500 | 50 | d |
| 8 | 500 | 100 | d |

**Methods, 17 total, identical across cells:**

- **4 proposed (R, `tr1/wp0_mccd_methods.R`)**: U-MCCD, SU-MCCD, UN-MCCD,
  SUN-MCCD, via `umccd_method`/`sumccd_method`/`unmccd_method`/
  `sunmccd_method`. `S_min = min.cls = 0.05` passed to SU-MCCD and
  SUN-MCCD (the only two that take it); `method = "ascend"` for UN-MCCD
  and SUN-MCCD at every `d` — this is the paper's **own** stock
  configuration (confirmed by WP9's ablation: "both statistics, centre
  removed, ascend... is retained... every published result uses it"),
  **not** TR2's generic `d`-dependent ascend/descend switch, which applies
  to a different construction and would be the wrong analogy here. `quant`
  is left `NULL` in every call so each wrapper resolves alpha through the
  paper's own resolver (`rk_quant_label_paper`/`nn_quant_label_paper_UN`/
  `nn_quant_label_paper_SUN`, all in `shared/harness.R` §2) — runtime is
  measured under the manuscript's actual operating alpha at each `d`, not
  a fixed level.
- **5 old baselines (R, `shared/harness.R` registry)**: LOF, DBSCAN, MST,
  ODIN, iForest — `METHOD_REGISTRY[["LOF"|"DBSCAN"|"MST"|"ODIN"]]` called
  unchanged. **iForest is the one exception**: the registry's
  `iforest_method` does not pin `nthreads`, and `isotree` defaults to all
  cores, which would break single-worker comparability. This driver
  therefore does **not** call `METHOD_REGISTRY[["iForest"]]`; it defines
  its own `iforest_1thread()` (same hyperparameters — `ntrees=1000,
  sample_size=min(256,n)` — plus `nthreads=1`), the same pattern TR2 used
  in `04_wp4_runtime.R`, and calls that instead. This is an addition inside
  the new WP6 script, not an edit to `harness.R`.
- **8 new competitors (Python, same settings as `tr1/81_wp4_baselines.py`,
  timed by `91b_wp6_runtime_py.py`)**: ECOD, COPOD, HDBSCAN, OPTICS,
  MutualKNN, SNN, DIF, LUNAR. `MutualKNN`/`SNN` are timed once per `k` in
  `{5,10,15,20,30}` (own code, `sklearn.neighbors.NearestNeighbors`);
  `DIF`/`LUNAR` once per seed in `{1,...,5}` — same replication scheme as
  `81_wp4_baselines.py`, reused for consistency with the metrics side of
  this project rather than TR2's own PyOD runtime script (which times a
  different, three-method set: ECOD/LUNAR/AutoEncoder). `--reps` applies
  per (cell, method, k-or-seed) exactly as for the R side.

**Deliberately excluded**: `shared/harness.R`'s generic `RKCCD-OOS/IOS`,
`UNCCD-OOS/IOS` wrappers. They are not this paper's four named methods and
not one of its baselines; they exist for TR2's own study. They are also
the wrappers the harness audit flagged (see §5) — excluding them sidesteps
that pitfall entirely rather than working around it.

## 4. What is timed, and how "single-threaded" is enforced

**The double-construction pitfall, and why it does not apply to the four
proposed methods.** `shared/harness.R`'s `rkccd_oos_method`/
`rkccd_ios_method`/`unccd_oos_method`/`unccd_ios_method` each call the
underlying construction function **twice** — once (wrapped in
`invisible()`) purely to produce a `t_construct` number, and again inside
the actual scoring call that produces `t_total` — so a `t_total` read off
those four wrappers is not the method's true one-pass cost; a real
single-pass timing would have to call the detector once and use that. This
driver never touches those four wrappers (§3), so the pitfall is moot for
them — but the *general lesson applies to how this driver times
everything*: **every method here is timed as exactly one call**, wrapped
in `proc.time()`, nothing before or after it inside the timed region. For
the four MCCD methods this is additionally already guaranteed by
`wp0_mccd_methods.R` itself, which documents that its four wrappers are
"single monolithic calls... `t_construct` and `t_total` are the SAME
single measured wall-clock time" — this driver's own outer `proc.time()`
call and that inner measurement time the identical span, so no double
counting is possible even by construction.

**Timed region, precisely, per rep:**

```r
invisible(gc(reset = TRUE))
t0 <- proc.time()
res <- do.call(METHOD_REGISTRY_OR_LOCAL[[method]], list(X = X, d = d, Y = Y, ...))
elapsed <- (proc.time() - t0)[["elapsed"]]
gc_after <- gc(reset = FALSE)
```

`elapsed` (wall clock, `proc.time()`'s `elapsed` field, matching both
`04_wp4_runtime.R`'s and `85_wp8`'s convention) is `time_s`. Excluded from
the timed region: data generation (§2), quantile-table loading for the
MCCD methods (below), and result post-processing / CSV writes.

**Quantile-table load timing, reported separately.** The four MCCD
wrappers in `wp0_mccd_methods.R` (off-limits to edit) each call
`get_simul()` themselves, before their own internal `t0` — so the load is
already excluded from the number the wrapper hands back. To surface a
`table_load_s` figure without editing that file, this driver does its own
`get_simul(variant, d, quant = <same resolver>, n = nrow(X))` call
immediately before invoking the wrapper, timed with its own
`proc.time()`, and reports that as `table_load_s`. **This means the table
is loaded twice per rep in practice** (once by this driver for the
measurement, once again inside the wrapper) — both loads sit outside any
region counted toward `time_s`, so it does not inflate the reported
runtime; it is flagged here as a known, harmless redundancy forced by the
append-only constraint on `wp0_mccd_methods.R`. `table_load_s` is `NA` for
the 5 R baselines and all 8 Python competitors.

**Single-threading, and how it is verified:**

- R: `Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")` set at the top of
  `91_wp6_runtime.R`, before `shared/harness.R` is sourced (same ordering
  TR2 used). Verified two ways: (a) `extSoftVersion()["BLAS"]` /
  `La_library()` return empty strings on this R 4.6.1 install, i.e. R is
  linked against the single-threaded reference BLAS, not an
  OpenMP-parallel one, so the env vars are a defensive no-op rather than
  load-bearing; (b) a source-level check (`grep`) of
  `R/ccds/UN_CCD.R`, `R/ccds/RK_CCD_New.R`, and the four
  `methods/outlier_detection/*.R` files for `parallel`/`foreach`/
  `mclapply`/`parSapply`/`makeCluster` returns nothing — none of the four
  detectors' own code paths ever parallelizes, unlike TR2's deliberate
  22-core override of `nnccd.radi` (§0), which this driver does not use.
  **iForest is the one method that needs an explicit override, not just an
  env var**: `isotree::isolation.forest(..., nthreads = 1)`, passed
  directly (§3), since `isotree`'s OpenMP thread count is controlled by
  its own argument, not reliably by the environment variable alone.
  `threads` is recorded as `1` on every row for provenance.
- Python: the same four env vars plus `VECLIB_MAXIMUM_THREADS` are set
  **before** `numpy`/`torch` are imported (order matters — BLAS thread
  pools are fixed at import time), and `torch.set_num_threads(1)` /
  `torch.set_num_interop_threads(1)` are called explicitly, identical to
  TR2's `05_wp4_runtime_pyod.py`. `hdbscan` and scikit-learn's `OPTICS`
  have no separate thread knob beyond BLAS/OpenMP, which the env vars
  cover. `threads` is recorded as `1` on every Python row too.

## 5. Memory

**R**: `gc(reset = TRUE)` immediately before the timed call (§4);
`gc(reset = FALSE)` immediately after. `gc()`'s return matrix has columns
`used, (Mb), gc trigger, (Mb), max used, (Mb)` for its two rows (`Ncells`,
`Vcells`); column 6 (`max used (Mb)`) is the high-water mark reached
between the two `gc()` calls, for each of R's two heaps. `mem_peak_mb =
sum(gc_after[, 6])`, the standard R idiom for the peak-since-reset
figure. This is a coarse, R-heap-only number — it does not see memory a
C-level allocation outside R's GC ever mapped and freed without R
noticing, but it is "nearly free" (`REVISION_PLAN`'s own phrase) and is
what R3.10/R5.5 ask for.

**Python**: `tracemalloc` peak, not `psutil` RSS delta. Rationale: RSS
delta is contaminated by the interpreter's own long-lived allocations
(numpy/torch/sklearn import footprint, which happens once per process, not
once per fit) and by OS-level allocator behaviour (freed memory is not
always returned to the OS, so a later cell's RSS delta can read
artificially low or even negative). `tracemalloc.start()` /
`tracemalloc.get_traced_memory()` (returns `(current, peak)` bytes since
the last `start()`/`clear_traces()`) brackets exactly the `model.fit(X)`
call, `tracemalloc.stop()` after, mirroring the R side's reset-before /
read-after structure. **Known limitation, stated in advance**: PyTorch
tensor allocations inside `DIF`/`LUNAR` that live on the C++ side (not
through Python's allocator) are invisible to `tracemalloc` — those two
methods' `mem_peak_mb` will under-count relative to their true footprint,
and this is disclosed in the summary output, not hidden.

## 6. Idle-machine requirement

Both timing scripts refuse to start the **real** grid (`--smoke` and
`--summarize` are exempt — the former is seconds of compute, the latter
touches no timed region at all) unless the machine looks idle:

- R (`91_wp6_runtime.R`): lists all `tasklist /FO CSV /NH` processes,
  filters to `Rscript.exe`/`Rterm.exe`/`R.exe`/`python.exe`/`pythonw.exe`,
  excludes this process's own PID (`Sys.getpid()`), and refuses to start
  (printing the offending PIDs and image names) unless `--force` is
  given. `--idle-check` alone runs just this check and exits, without
  starting anything — useful for a verifier to confirm before launching
  the real grid.
- Python (`91b_wp6_runtime_py.py`): same check via `psutil.process_iter()`
  (or, if `psutil` is unavailable in the pinned venv, a `tasklist`
  subprocess identical to the R side), same `--force` / `--idle-check`
  contract.

**Why this matters here specifically**: `HANDOFF_FROM_TR2.md` and
`tr2/FINDINGS.md` both record that TR2 lost a full WP4 runtime grid to
undetected co-scheduling and had to re-time the contention-flagged cells
from scratch (`FINDINGS.md` §10 — "d=10 and n=2000 RKCCD cells re-run on
an idle machine... corrections to §6's headline numbers"). The idle-check
exists to catch that failure mode before it happens rather than after.

## 7. Complexity reconciliation: the log-log slope fit

For each of the four MCCD methods, fit
`log(time_s) ~ log(n)` by ordinary least squares on **every individual
rep** in sweep A (5 cells × 10 reps = 50 points per method, `d=10`
fixed) — not on the 5 per-cell medians — so the residual degrees of
freedom (48) come from genuine repeated measurement rather than from
treating 5 cell-medians as if they were 5 independent draws. The fitted
slope is the estimated exponent; `summary(lm(...))$coefficients["log(n)",
c("Estimate","Std. Error")]` gives the slope and its SE directly.

**Reported, per method**: fitted slope ± SE, `R^2`, and the explicit
comparison sentence — "stated bound: `n^3`(`× log n` for RK-based
methods); fitted: `<slope> ± <SE>`" — with **no editorial softening** if
the fitted value exceeds the stated one. A slope whose 95% CI excludes 3
(for UN-MCCD/SUN-MCCD) or excludes "consistent with `n^3\log n`" (for
U-MCCD/SU-MCCD, judged by comparing the fit to a null model with the
`\log n` factor folded in) is reported as a discrepancy in Section 4.7,
per `REVISION_PLAN`'s instruction to "say which, plainly." `--summarize`
writes this table to `results/tr1/wp6/91_wp6_slope.csv`.

Sweep B (the `d`-grid) is descriptive only — four points is too few to fit
a `d`-exponent with a meaningful SE, and the complexity bounds are stated
as fixed-`d` sweeps in `n` in the first place. Its role is confirming (or
not) that runtime tracks `n` and not `d`, echoing the WP0 gate's
real-data observation (waveform `d=21` and wilt `d=5` cost the same
order; §0) under a controlled synthetic design.

## 8. Estimated cost of the real grid

**Cost model, from the WP0 gate's real-data numbers (§0), not TR2's
parallel `~n^4` cross-codebase figure**, since these are this paper's own
detectors, unparallelized, exactly as the real grid will run them. Using
the vowels→wilt anchor (slope ≈ 3.3) and vowels' own point (`n=1452`,
≈17–26 s per proposed method) as the scaling base:

- `n=2000, d=10`: `(2000/1452)^3.3 ≈ 3.0×` vowels' time ⇒ roughly
  **50–80 s per MCCD method per rep**, i.e. ≈500–800 s per method across
  10 reps, ≈35–55 min for all 4 MCCD methods at this one cell.
- `n=1000, d=10`: ≈`0.35×` vowels ⇒ a few seconds per rep.
  `n≤500`: well under a second to a few seconds per rep.
- The `d`-grid (`n=500` fixed) is cheap throughout — real-data evidence
  (§0) says cost tracks `n`, not `d`, at these dimensions.
- 5 R baselines: LOF is the only one with a super-linear real cost driver
  (`MinPts` sweep 11–30 on an `n×n` neighbour computation); expect at most
  a few seconds per rep even at `n=2000`. DBSCAN/MST/ODIN/iForest are
  sub-second to a few seconds regardless of cell.
- 8 Python competitors: ECOD/COPOD/MutualKNN/SNN are near-linear and fast
  throughout; HDBSCAN/OPTICS are the `O(n^2)`-ish pair and may take up to
  low tens of seconds at `n=2000`; DIF/LUNAR (5 seeds each) are the
  slowest, plausibly tens of seconds per fit at `n=2000` given their
  neural-network training loops — call it up to ~10 minutes total for
  DIF+LUNAR's 100 timed fits (10 reps × 5 seeds × 2 methods) at the most
  expensive cell.

**Total estimate: on the order of 2–4 hours of wall clock for the entire
grid** (both sweeps, all 17 methods, 10 reps), dominated by the `n=2000`
MCCD cell — a very different picture from a naive extrapolation off TR2's
parallel `~n^4` UN-CCD number (which alone would suggest multi-hour-to-day
costs for a single `n=2000` rep). This is exactly the ambiguity §0 and §7
exist to resolve with a real, controlled measurement rather than an
analogy. As a defensive measure only (not expected to fire), each rep
still runs under a **30-minute per-rep `setTimeLimit`**, matching the
spirit of TR2's timeout discipline at a much shorter cap appropriate to
this paper's measured (not TR2's) cost model; a timeout is recorded as
`status = "FLAGGED_TIMEOUT"` and the remaining reps of that
`(cell, method)` are skipped, never silently dropped.

## 9. Files

- `revision_experiments/tr1/91_wp6_runtime.R` — R side (4 MCCD + 5
  baselines).
- `revision_experiments/tr1/91b_wp6_runtime_py.py` — Python side (8 WP4
  competitors), reading the CSV exports `91_wp6_runtime.R` writes.
- `results/tr1/wp6/91_wp6_runtime_raw.csv` — one row per (method, n, d,
  rep): `method, n, d, rep, seed, time_s, mem_peak_mb, table_load_s,
  threads, status`.
- `results/tr1/wp6/91_wp6_runtime_done.csv` — **keys-only** checkpoint
  file (`method, n, d, rep`), separate from the raw file on purpose:
  `harness.R`'s `has_result()` treats any `NA` in a non-key payload column
  as an incomplete/partial row, and `table_load_s` is legitimately `NA`
  for every baseline by design (not because of truncation) — checkpointing
  against the raw file directly would make every baseline rep look
  perpetually unfinished and rerun forever. The done file carries no
  payload columns, so that check never fires against it.
- `results/tr1/wp6/data/<grid>_<cell_value>_rep1.csv` — rep-1 dataset
  export per cell (features `V1..Vd` + `label`), for the Python script.
  Same asymmetry as `04_wp4_runtime.R`/`05_wp4_runtime_pyod.py`: the R
  side draws a fresh dataset per rep, the Python side times all its reps
  against the single rep-1 export (its reps vary model seed / system
  noise only). Declared here, not hidden.
- `results/tr1/wp6/91b_wp6_runtime_py_raw.csv`,
  `results/tr1/wp6/91b_wp6_runtime_py.csv` — Python raw and aggregate.
- `results/tr1/wp6/91_wp6_runtime_n.csv`, `_d.csv` — `--summarize` output
  (median, CV per cell/method).
- `results/tr1/wp6/91_wp6_slope.csv` — `--summarize` output (§7).
- `results/tr1/wp6/smoke/` — `--smoke` output (R and Python), same file
  names as the corresponding production paths, isolated by directory so a
  smoke run can never contaminate or be mistaken for a production row.

## 10. Handoff

This author sub-session's job stops at: protocol committed, scripts
written, `parse()`/`py_compile` clean, `--smoke` proven (including the
done-file skip on a second `--smoke` call), everything committed. **The
real grid, on an idle machine, is the verifier's or a later session's
task.** Before running it: confirm the four background R jobs mentioned in
the task brief have finished, run `--idle-check` on both scripts, and only
then run without `--force`.

## 11. Fixer notes, 2026-09-06 (post-verification changes)

Opus verified this package on 2026-09-06 (`tr1/WP6_VERIFICATION.md`):
**PASS WITH REQUIRED CHANGES**. All required changes (A1, A2, B1, B2, C, D,
E, the summariser items, the caption disclosures) have now landed in
`91_wp6_runtime.R` / `91b_wp6_runtime_py.py`. This section records what
changed and amends two sentences in the declarations above that the
verification found to be wrong. **Nothing above this section was edited.**

**A1 (blocker, fixed)** — `get_simul()` is now memoized inside
`91_wp6_runtime.R` only (one cache entry, evicted on every new
`(variant, d, quant)` key; `check_simul_extent` still runs on every call).
`shared/harness.R` is untouched. Before the fix, every rep of every MCCD
method reloaded its quantile table from disk *inside* the timed region
(the wrapper's own internal `t_total` excludes the load, but this driver's
outer `proc.time()` around the whole wrapper call does not) — measured as
U-MCCD outer 3.96s vs wrapper-internal 1.28s at n=250. Smoke evidence after
the fix: U-MCCD `time_s` dropped from 2.95s to 0.36s at n=100; `SU-MCCD`,
`UN-MCCD`, `SUN-MCCD` each show `table_load_s ≈ 0` immediately after
U-MCCD's cold load populates the shared cache key.

**A2 (blocker, fixed)** — `table_load_s` is now measured once per
`(cell, method)`, hoisted above the rep loop, and shared by every rep of
that method at that cell (§4's timed-region description is otherwise
unchanged — table loading was never inside `time_s` itself, only inside
this driver's separate, always-untimed `table_load_seconds()` measurement,
which was simply being repeated wastefully once per rep before this fix).

**B1 (blocker, fixed)** — `mem_delta_mb` (peak minus pre-call baseline,
`sum(gc_after[,6]) - sum(gc_before[,2])`) is now reported alongside the
pre-existing `mem_peak_mb` (kept for provenance, per the verification's
instruction). Smoke evidence: U-MCCD's `mem_peak_mb` is still ≈864 MB
(dominated by the RK table's own unused 763 MB `Kest.m` matrix, resident
before the call and therefore part of the *baseline*, not the delta), but
`mem_delta_mb` is 55.3 MB — two orders of magnitude smaller, and now
comparable across methods regardless of which quantile table each one
loads.

**B2 (fixed)** — a memory log-log slope, `lm(log(mem_delta_mb) ~ log(n))`
at `d=10`, per MCCD method, is now written to `91_wp6_slope.csv` alongside
the time slope, compared explicitly against 2 (the manuscript's `O(n^2)`
space claim) rather than 3.

**C (major, fixed)** — `hdbscan.HDBSCAN(...)` now passes
`core_dist_n_jobs=1`. **Amendment to §4's Python single-threading
paragraph**: the sentence "`hdbscan` and scikit-learn's `OPTICS` have no
separate thread knob beyond BLAS/OpenMP" is **wrong for `hdbscan`** —
`core_dist_n_jobs` defaults to 4 and is controlled via `joblib`, not
OMP/BLAS, so the pinned env vars do not reach it. It remains true for
`OPTICS`. This is a documented deviation from `81_wp4_baselines.py` (which
leaves `hdbscan` at its default thread count) — timing only; detection
scores are unaffected by this knob.

**D (blocker, fixed)** — `CELLS$cell_value` for the three d-only cells
(6, 7, 8) now holds the cell's own `d` value ("5", "50", "100") instead of
n's value ("500", "500", "500"), which had made all three d-only export
paths collide on `d_500_rep1.csv`. `export_rep1()` additionally writes
cell 3 (n=500, d=10, shared by both sweeps) under **both**
`n_500_rep1.csv` and `d_10_rep1.csv`. A new `--check-exports` dry-check mode
(no idle-check, no data generation, no file writes — pure path-string
resolution over the `CELLS` table) confirms the fix directly: run and
verified 2026-09-06, "8 CELLS rows resolve to 9 export names (9 unique)."
`--smoke` itself only ever touches cell 1, so this dry check is the actual
confirmation that all 9 names are distinct across the full grid, not just
smoke's one cell.

**E (major, fixed)** — `from sklearn.neighbors import NearestNeighbors` and
`from scipy import sparse` are hoisted to module scope in
`91b_wp6_runtime_py.py` (were previously imported lazily inside
`knn_index()`/`knn_adjacency()`, i.e. inside the timed region on first
call). Additionally, one untimed warm-up fit per (method, variant) now runs
at the smallest available cell before any measured rep, for every method
including the previously-unaffected `ECOD`/`COPOD`/`HDBSCAN`/`OPTICS`/
`DIF`/`LUNAR` (first-call costs are not unique to the two hand-rolled
neighbour methods). **Warm-up fits are discarded**: not written to the raw
CSV, not written to the done file, logged only as `[warm-up] <method>
<variant> ok` lines to stdout — chosen over logging them as rows because a
warm-up fit is deliberately not comparable to a measured rep (different
JIT/cache state on either side of it) and keeping it out of the raw CSV
means no downstream aggregation code has to know to filter it out. Smoke
evidence: `LUNAR seed1` was 5.03s/65MB vs `seed2` 1.93s/0.19MB before the
fix; after the fix, `seed1` 1.96s/0.16MB vs `seed2` 2.00s/0.16MB — no longer
distinguishable from sampling noise.

**Summariser (fixed)** — `91_wp6_runtime_n.csv`/`_d.csv` now carry
`mean_time_s` beside `median_time_s` (CV = sd/mean needs its own centre)
and `median_mem_delta_mb`/`cv_mem_delta` beside the legacy
`mem_peak_mb` columns. `91_wp6_slope.csv` gains, per MCCD method: a `mem`
row (§B2, compared to `n^2`) and, for the two RK-based methods only, a
`time_folded` row fitting `log(time_s/log(n)) ~ log(n)` and comparing the
result directly to 3 (rather than "3 with a log factor", which cannot be
compared to a single fitted number without folding). `--summarize` now also
reads `91b_wp6_runtime_py_raw.csv`, when present, and writes
`91b_wp6_runtime_py_summary.csv` applying the collapse rules pre-declared
here: **MutualKNN/SNN summarized at `k=10` only** (median over reps; the
full `k ∈ {5,10,15,20,30}` sweep stays in the raw CSV for the supplement),
**DIF/LUNAR summarized as the median over all seed × rep fits pooled
together** (not per-seed), and the four deterministic methods
(ECOD/COPOD/HDBSCAN/OPTICS) as the median over reps of their one variant.
Each output row records which rule produced it, in a `collapse_rule`
column. Dry-run tested 2026-09-06 against a temporary copy of the smoke raw
CSVs (copied into the production paths, summarized, then deleted again —
the production `results/tr1/wp6/` tree is empty again, as it was before
this session and as the real grid still expects); with only one `n` value
present every slope fit correctly falls back to `NA` (fewer than 3 points),
confirming the fallback path rather than the fitted path, which is the
correct behaviour for a single-cell smoke sample.

**Caption disclosures (added, for whichever section/table in the
manuscript reports this grid's results)**:
- DBSCAN's `eps` is tuned using the **true** contamination rate of the
  synthetic generator, not an oracle-free rule — a methodological
  convenience for this runtime grid, not the paper's real-data protocol.
- MST is run at `cont = 0.05` (this generator's own contamination rate),
  not the registry's default `cont = 0.02` — set explicitly in
  `call_method()` to match `gen_uniform_wp6()`.
- ODIN's `k = round(sqrt(n))` grows with `n` by construction, so ODIN's
  measured runtime slope is **not** a fixed-parameter comparison the way
  the four MCCD methods' slopes are (their `N`, `d` are held fixed by
  design; ODIN's own hyperparameter is not).
- R and Python `mem_*_mb`/`mem_delta_mb` columns are **not comparable
  across languages**: R's figures come from `gc()`'s heap accounting
  (§5), Python's from `tracemalloc` (§5) — different allocators,
  different blind spots (R misses non-GC C-level allocation; Python's
  `tracemalloc` misses PyTorch's C++-side tensor allocations inside
  DIF/LUNAR, as already disclosed in §5).
- The idle-machine check (§6) **fails open**: it proceeds as idle if
  `tasklist` itself is unavailable, and it only ever sees
  `Rscript.exe`/`Rterm.exe`/`R.exe`/`python.exe`/`pythonw.exe` processes —
  any other process consuming CPU/memory (another user's job, a scheduled
  task, a non-R/Python compute process) is invisible to it.

**Cost ruling (reaffirmed, unchanged)** — the estimate in §8 stands: no rep
cap is needed for the real grid, and `--grid=n` and `--grid=d` should still
be run as two separate invocations (as already specified in the launch
sequence in §10-adjacent material above), not combined into one `--grid=all`
call, so a failure or interruption in one sweep does not require re-running
the other.

Both scripts re-verified 2026-09-06 after all of the above:
`parse()`/`py_compile` clean; `--check-exports` confirms 9 distinct export
names; `--smoke` run twice on each side confirms the checkpoint skip fires
on the second call; the summarizer was dry-run against a temporary copy of
the smoke CSVs (not committed, deleted after inspection) and produced all
expected columns with no errors. The real grid was **not** run by this
fixer session, per the task brief — WP8 background R jobs were still the
stated reason to stay off the machine, and this session's job was the
required-changes list, not the grid itself.

## Dated note (2026-09-06 04:40, coordinator)
During the production run, one manuscript-writing agent was active on the same machine: a single short Python table-generation script at its start and one pdflatex/bibtex build sequence (single-threaded, under 15 minutes total) at its end. No R process and no other compute ran. The idle check passed at launch; this overlap is disclosed as the only known deviation from the idle-machine requirement.
