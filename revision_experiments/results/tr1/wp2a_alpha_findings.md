# WP2(a) — α (quantile level) sensitivity sweep

Driver: `revision_experiments/42_wp2a_alpha_sweep.R`
Per-cell results: `revision_experiments/results/tr1/wp2a_alpha_sweep.csv`
Per-group summary: `revision_experiments/results/tr1/wp2a_alpha_summary.csv`

Answers reviewer points R1.4, R3.1, R3.2, R5.2 and AE.2 — all of which object that α
in the submitted manuscript varies with both dimension and method with no stated
justification.

α enters the four detectors only through *which Monte-Carlo quantile table the
spatial-randomness test consults*. U-MCCD and SU-MCCD read the RK tables
(`R/RK-test_quantile/`); UN-MCCD and SUN-MCCD read the NN tables
(`R/NN-test_quantile/`). Nothing else in the pipeline is touched.

`S_min = 0.0625` (proportion of n) throughout, the manuscript value.

---

## Step 1 — which α levels exist on disk

Filename tokens present, for each of the 12 distinct dimensions carried by the 16
real data sets. Token `999` = quantile 0.999, i.e. α = 0.001.

| d | RK tokens | RK n | NN tokens | NN n |
|---|---|---|---|---|
| 5 | 99 999 | 2 | 90 95 | 2 |
| 6 | 99 999 | 2 | 95 99 | 2 |
| 7 | 99 999 | 2 | 95 99 | 2 |
| 8 | 99 999 | 2 | 95 99 | 2 |
| 9 | 99 999 | 2 | 95 99 | 2 |
| 10 | 99 999 | 2 | 99 999 | 2 |
| 12 | 999 | 1 | 99 999 | 2 |
| 16 | 999 | 1 | 99 999 | 2 |
| 18 | 999 | 1 | 99 999 | 2 |
| 19 | 999 | 1 | 99 999 | 2 |
| 21 | 999 | 1 | 99 999 | 2 |
| 30 | 999 | 1 | 999 | 1 |

Plainly: **no dimension carries more than two α levels on disk, and six of the twelve
carry only one on the RK side.** Taken at face value this makes the sweep a
two-point comparison at best. Step 2 changes that for RK.

Reproduce with `Rscript revision_experiments/42_wp2a_alpha_sweep.R --inventory`.

---

## Step 2 — can α be changed without regenerating anything?

**The answer is different for the two variants, and it is decisive for the cost of
this whole work package.**

### RK — yes, every α is free

`R/RK-test_quantile/RK-test-simul_9d_99%.RData` contains one object, `simul`, a
3-element list:

```
$Kest.m   matrix    2000 x 50000     (= niter x (m * rn), m = 5000, rn = 10)
$quan     list      1 key: "0.99"
$r        numeric   length 10
```

`$Kest.m` is the **full matrix of raw per-iteration Monte-Carlo draws**, retained
verbatim (`Kest.R:412` returns `list(Kest.m = Kest.m, quan = Kest.quan, r = r)`).
The reduced envelope in `$quan` is nothing but

```r
temp <- apply(Kest.m, 2, quantile, probs = q);  matrix(temp, nrow = m)
```

(`Kest.R:406-411`). So **any α is recoverable from a file already on disk at zero
simulation cost** — one `apply()` over the stored draws, about 3 seconds per level.
That is also why the RK tables are 265–711 MB each and the NN tables are 70 KB.

This was verified, not assumed. Before any derived cell was allowed to run, the
driver recomputes the base file's **own** stored `$quan` entry from its `$Kest.m` and
aborts unless the recomputation is `all.equal(..., tolerance = 0)` identical to the
stored matrix. The guard passed on every base file loaded (`[rk base] derivation
guard OK` lines in the run log). The RK half of this sweep is therefore **dense
(α ∈ {0.10, 0.05, 0.01, 0.001}) at every dimension, including the six where only one
file exists.**

### NN — no, regeneration is required

`R/NN-test_quantile/NN-test-simul_9d_95%.RData` contains one object, `simul`, a
2-element list:

```
$average  numeric   length 5000
$median   numeric   length 5000
```

These are already-reduced quantile vectors — exactly what `nnccd.radi()` slices as
`simul$average[1:n]` / `simul$median[1:n]` (`harness.R` ~line 161,
`R/ccds/UN_CCD.R:241`). The generator discards the draws inside its own body:
`NNDestP.simpois.lower.quant` returns `list(average = ..., median = ...)` and nothing
else (`R/ccds/NN_Dist_Est.R:111`). **A new α on the NN side needs a new Monte-Carlo
run.** The NN half of the sweep is limited to the tokens on disk.

---

## Step 2b — the NN tables at nine of eleven dimensions are duplicate files

This was not part of the brief; it fell out of Step 2 and it changes the
interpretation of the headline high-dimensional result, so it is reported first.

Comparing the two NN tables that exist at each d
(`Rscript revision_experiments/42_wp2a_alpha_sweep.R --nncmp`, corroborated by MD5):

| d | tokens compared | `$average` identical? | max abs diff | MD5 |
|---|---|---|---|---|
| 5 | 90 vs 95 | **no** | 7.6e-2 | differ |
| 6 | 95 vs 99 | yes | 0 | **same** |
| 7 | 95 vs 99 | yes | 0 | **same** |
| 8 | 95 vs 99 | yes | 0 | **same** |
| 9 | 95 vs 99 | yes | 0 | **same** |
| 10 | 99 vs 999 | **no** | 1.3e-1 | differ |
| 12 | 99 vs 999 | yes | 0 | **same** |
| 16 | 99 vs 999 | yes | 0 | **same** |
| 18 | 99 vs 999 | yes | 0 | **same** |
| 19 | 99 vs 999 | yes | 0 | **same** |
| 21 | 99 vs 999 | yes | 0 | **same** |
| 30 | — | only one table | — | — |

At d ∈ {6,7,8,9,12,16,18,19,21} the two files are **byte-identical** — same MD5, same
byte length, same modification timestamp to the minute. One Monte-Carlo run was
saved twice under two names. Only d = 5 and d = 10 hold two genuinely distinct NN
tables.

Consequence: **the NN detectors returning identical output at "99%" and "99.9%" is
not evidence that α has stopped binding.** It is guaranteed by the input files being
the same numbers. Any conclusion of the form "the MC-SRT saturates at d = 21" that
rests on comparing those two tables is unsupported. See the Saturation section.

The RK tables show no such problem — all twelve pairs checked have distinct MD5s and
modification times about a day apart, consistent with genuinely separate runs.

---

## Step 3 — coverage actually achieved

`wp2a_alpha_sweep.csv`, column `status`:

| status | cells |
|---|---|
| `ok` | 144 |
| `timeout` | 1 |
| `not_run` | 65 |

**11 of the 16 data sets are fully covered**: hepatitis, lymphography, glass, WBC,
vertebral, ecoli, stamps, WDBC (the eight small sets), plus pima, shuffle and vowels.

Per data set the plan is: 4 derived RK levels × 2 RK methods, plus one "disk" RK
replicate per method wherever a second RK file exists at that d, plus every NN token
on disk × 2 NN methods.

**5 data sets were not run**: pageblocks, wilt, PenDigits, waveform, thyroid. The
justification is measured, not assumed. `UN-MCCD` on pageblocks (n = 4795, d = 10)
**exceeded the 900 s cell timeout** — that single cell is in the CSV with
`status = timeout`, `elapsed_sec = 900.1`. The detectors are roughly O(n³): U-MCCD on
vowels (n = 1452) took 58 s (`elapsed_sec`), which projects to 600–2100 s per cell at
n = 3200–4819, times 12–14 cells per data set, i.e. 2–8 h per data set. That is
outside the WP2(a) budget. The 65 unrun cells carry `status = not_run` and the same
reason string in `note`.

Cost of what was run: ~25 min wall clock for the 144 cells, on a contended machine
(three other agents were running R jobs concurrently), so `elapsed_sec` is an upper
bound, not a clean benchmark.

Losses from the truncation, stated plainly: **pageblocks (d = 10) and wilt (d = 5)
are the only two data sets whose NN tables at two α levels are genuinely different**
(Step 2b). They are also the two largest. The consequence is that **this sweep
contains no informative measurement of NN-side α sensitivity at all** — every NN
group that ran compares two duplicate tables. That gap is a direct consequence of the
duplicate-table problem, not of the compute budget alone: even with unlimited compute,
only 2 of the 16 data sets could have contributed an informative NN comparison from
the tables on disk.

---

## Step 4 — how much does α move the answer?

### Headline

Of the 42 (dataset × method × source) groups with at least two α levels, **22 change
the flagged set and 20 do not** (`wp2a_alpha_summary.csv`, column
`labels_identical_across_alpha`). But that 20 is entirely the 20 uninformative NN
groups built on duplicate tables. Restricting to the groups where the underlying
tables genuinely differ:

> **Every single informative group — all 22 RK groups, at every dimension from 6 to 30
> — changes its flagged set when α changes. Not one is stable.**

### Magnitude, RK derived groups (22 groups, α ∈ {0.10, 0.05, 0.01, 0.001})

Columns `BA_range`, `F2_range`, `min_pairwise_jaccard`, `n_flagged_min/max` in
`wp2a_alpha_summary.csv`.

| statistic across the 22 groups | value |
|---|---|
| median BA range | 0.160 |
| mean BA range | 0.206 |
| median F2 range | 0.158 |
| median minimum pairwise Jaccard | 0.216 |

A median Jaccard of 0.22 between the flagged sets at two α levels means the two
"answers" overlap by about a fifth. This is not a second-order effect.

Largest moves:

- **BA: lymphography × U-MCCD, d = 18 — 0.468 → 0.951, range 0.482.**
  BA by α: `900=0.6350; 950=0.4683; 990=0.5411; 999=0.9507`. Jaccard floor 0.074.
- **F₂: the same cell — range 0.682** (0.000 at α = 0.05 to 0.682 at α = 0.001).
- Next largest BA: shuffle × U-MCCD, d = 9, range 0.428, Jaccard floor **0.017** —
  the two flagged sets are effectively disjoint.
- Third: lymphography × SU-MCCD, d = 18, range 0.419.

Smallest moves: vowels × SU-MCCD (BA range 0.054, Jaccard 0.846), WBC × SU-MCCD
(0.070), ecoli × U-MCCD (0.067). Even these change the label set.

### Is there a plateau?

**No usable one.** BA is **non-monotone in α in 19 of the 22 RK groups** — it rises,
falls and rises again along the ladder (e.g. vertebral × U-MCCD:
`900=0.507; 950=0.393; 990=0.536; 999=0.564`; stamps × U-MCCD:
`900=0.530; 950=0.837; 990=0.511; 999=0.653`). Only 3 groups are monotone, and two of
them are the same data set (WBC × U-MCCD and × SU-MCCD, both increasing; WDBC ×
SU-MCCD, decreasing). Short flats do occur, always at the permissive end: 4 of 22
groups are flat between α = 0.10 and α = 0.05, and WBC × SU-MCCD is flat at
BA = 0.7465 across α ∈ {0.10, 0.05, 0.01} before jumping to 0.817 at α = 0.001. There
is no interval of α over which the four
detectors are jointly insensitive, and no α that is best on more than a handful of
data sets. The manuscript's per-d, per-method α schedule cannot be defended as
"anywhere in this range works"; on this evidence it is a choice the results depend on.

### Monte-Carlo replicate noise is a confound of comparable size

Where two genuinely different RK tables exist at the same d, the sweep runs the
detector at the *same nominal α = 0.001* off both — once derived from the base file's
draws, once loaded from the second file (`alpha_source` = `derived` vs `disk`). These
should agree up to Monte-Carlo error alone.

7 of 14 such pairs are label-identical. The other 7 differ, sometimes grossly:
shuffle × U-MCCD flags 269 points from one table and 39 from the other at the same α
(Jaccard 0.069); ecoli × SU-MCCD flags 198 vs 88 (Jaccard 0.444); vertebral × U-MCCD
69 vs 91. **Re-running the quantile simulation at a fixed α can change the answer as
much as changing α does.** Any α recommendation the paper makes should be accompanied
by a statement of how many Monte-Carlo iterations the table behind it used.

### Saturation at high d — refuted as stated, and relocated

The prior claim on this codebase was that vowels (d = 12) and waveform (d = 21)
SUN-MCCD return byte-identical output at the 99 % and 99.9 % NN tables, read as α
ceasing to bind at high dimension.

The output is indeed byte-identical — reproduced here for vowels: UN-MCCD and
SUN-MCCD give `BA = 0.7706` / `0.7704` at both tokens, `flagged_idx` identical.
**But the cause is not saturation.** `NN-test-simul_12d_99%.RData` and
`NN-test-simul_12d_999%.RData` are the same file byte-for-byte (Step 2b). The
detector was handed the same numbers twice. The same holds at d = 18, 19 and 21, so
the waveform observation has the same explanation.

**There is therefore no measured NN saturation dimension.** The only two dimensions
where the question could be asked from the tables on disk are d = 5 and d = 10, and
those are wilt and pageblocks, both unrun (above).

On the RK side, where every comparison is informative, **α binds at every dimension
tested — 6, 7, 8, 9, 12, 18, 19 and 30 — with no saturation up to d = 30.** The
Jaccard floor does not trend upward with d: it is 0.207 at d = 6, 0.017 at d = 9,
0.074 at d = 12 and 18, 0.112 at d = 30. If anything the RK test becomes *more*
α-sensitive at moderate d, not less.

What does trend with d is over-flagging, and it is visible here: at d = 30 (WDBC)
SU-MCCD flags 194–341 of 367 points depending on α, and at d = 12 (vowels) 981–1160
of 1452. The failure mode is real; it is just not attributable to α saturating.

### α versus S_min

Reported for the α side only, per the brief. α changes the flagged set in 100 % of
informative groups, with a median BA swing of 0.160 and a median Jaccard of 0.216.
Whoever holds the S_min measurement can compare against those three numbers directly
(`wp2a_alpha_summary.csv`: `labels_identical_across_alpha`, `BA_range`,
`min_pairwise_jaccard`, filtered to `variant == "RK" & alpha_source == "derived"`).

---

## Step 5 — cost of a dense grid

### Measured

`Rscript revision_experiments/42_wp2a_alpha_sweep.R --cost 5 200` ran the production
generator `01_gen_quantile_table.R` for **NN, d = 5, n = 1000, niter = 200, cores = 4**,
writing to `results/tr1/cost_probe/` (never to `R/NN-test_quantile/`).

```
Elapsed: 470.6 s for niter = 200  ->  2.355 s per iteration at 4 cores
```

The machine was contended (three other agents running R), so this is an upper bound.

### Projected

Production setting is n = 1000, niter = 10000 (the "high d" tier reverse-engineered in
`01_gen_quantile_table.R`'s header, and confirmed by the one surviving driver
`R/NN-test_quantile/50100d_999%.R`).

| quantity | projection |
|---|---|
| one NN table, d = 5, niter = 10000, 4 cores | **23 554 s ≈ 6.5 h** |
| one α level across the 12 real-data dimensions | ~78 h (cost per iteration grows with d) |
| a 4-level dense NN grid over 12 dimensions | **~310 h, lower bound** |

The d = 5 figure is the cheap end; per-iteration cost rises with d through the
`dist()` and sampling steps, so the multi-dimension figures are lower bounds. **No
production run was launched.**

### Recommendation: one run should emit every α

**Feasible, and cheap to implement.** The NN generator already materialises the full
`niter × n` matrices `NN.dist.ave.mat` and `NN.dist.med.mat` internally and then
collapses them with a single `quantile(NN.dist.ave.mat[, x], 1 - quant)` per column
(`R/ccds/NN_Dist_Est.R:102-111`). `quant` is a scalar there purely by convention — the
RK generator, in the same codebase, already takes a *vector* `quan` and returns a
named sub-list `Kest.quan[[as.character(cur_quan)]]`, one entry per level
(`Kest.R:406-411`). The NN side can adopt the identical pattern:

1. accept a vector `quant`;
2. replace the two `sapply(..., quantile, 1 - quant)` reductions with a loop over the
   requested levels, storing `list(average = ..., median = ...)` per level;
3. return the named list; teach `get_simul()` to index it.

That turns **one** Monte-Carlo run into the whole α grid. Cost is unchanged (the draws
are already being computed and thrown away), the change is ~15 lines in a new override
file — the same "override after sourcing, never touch the original" pattern the
revision experiments already use — and it needs no change to the detectors.

The alternative, which is strictly better still: **also persist the raw
`niter × n` draw matrices**, exactly as the RK generator does. At n = 1000,
niter = 10000 that is ~80 MB per component in double precision (the RK tables are
265–711 MB, so this is well within the precedent already set on disk), and it makes
every future α, at every d, free forever — which is precisely what makes the RK half
of this sweep dense while the NN half is stuck at two points.

**Both are recommendations only. No generator was written and no regeneration run was
started.**

---

## Reproducing

```bash
cd Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection
Rscript revision_experiments/42_wp2a_alpha_sweep.R --inventory   # Steps 1, 2, 2b
Rscript revision_experiments/42_wp2a_alpha_sweep.R --nncmp       # Step 2b alone
Rscript revision_experiments/42_wp2a_alpha_sweep.R "hepatitis,glass"   # grid subset
Rscript revision_experiments/42_wp2a_alpha_sweep.R --summary     # Step 4
Rscript revision_experiments/42_wp2a_alpha_sweep.R --cost 5 200  # Step 5
```

The grid is checkpointed on (`dataset`, `method`, `alpha_token`, `alpha_source`), so a
re-invocation resumes rather than recomputes.

### CSV columns backing each claim

| claim | file | column(s) |
|---|---|---|
| 144 cells ran, 1 timed out, 65 not run | `wp2a_alpha_sweep.csv` | `status` |
| which α level a row used | `wp2a_alpha_sweep.csv` | `alpha_token`, `alpha_q`, `alpha` |
| derived vs on-disk table | `wp2a_alpha_sweep.csv` | `alpha_source`, `base_token` |
| the exact flagged set | `wp2a_alpha_sweep.csv` | `flagged_idx` (`;`-joined row indices) |
| label set changes with α | `wp2a_alpha_summary.csv` | `labels_identical_across_alpha`, `min_pairwise_jaccard` |
| size of the BA / F₂ move | `wp2a_alpha_summary.csv` | `BA_range`, `F2_range`, `BA_by_alpha` |
| over-flagging at high d | `wp2a_alpha_sweep.csv` | `n_flagged`, `n`, `d` |
| cost justifying the truncation | `wp2a_alpha_sweep.csv` | `elapsed_sec`, `note` |
