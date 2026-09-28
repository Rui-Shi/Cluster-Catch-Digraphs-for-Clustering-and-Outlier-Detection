# WP2a alpha-provenance: which alpha generated the shipped NN tables?

Covers **all nine duplicated dimensions that carry a real data set** —
d = 6, 7, 8, 9 (tier 1, 95/99 pairs), d = 18, 19 (tier 1, 99/999), and
d = 12, 16, 21 (tier 2, 99/999) — serving twelve data sets. Tier 2 was run
after tier 1's verdicts were verified independently by the coordinator
(`67_verify_provenance.R`).

Scripts: `revision_experiments/64_nn_provenance_gen.R` (generation),
`revision_experiments/65_nn_provenance_test.R` (tests).
CSVs: `wp2a_provenance_signtest.csv`, `wp2a_provenance_control.csv`,
`wp2a_provenance_behavioural.csv`, `wp2a_provenance_behavioural_summary.csv`,
`wp2a_provenance_gen_manifest.csv`.
Fresh tables: `R/NN-test_quantile_provenance/rep01..rep08/`. `R/NN-test_quantile`
was not modified.

## The question

21 pairs of shipped NN quantile tables are identical in content — one Monte
Carlo run saved under two level names (`62_audit_nn_duplicates.R`). Where the
pair is duplicated, the alpha that produced a published result is not readable
off the filename, and Table `tab:alpha_real` asserts one anyway.

**The level is a property of the TABLE, not of the data set.** Identifying the
alpha a shipped table holds requires only that the fresh replicates cover a
usable range of subsample sizes k — not that they reach the data set's n. That
is why d = 6 at n = 555 settles `thyroid` (n = 3656), d = 9 settles `shuffle`,
and d = 21 at n = 555 settles `waveform` (n = 3443), with no large-n generation
anywhere. It also bounds what Test 2 can be asked to do (see its truncation
note below).

## Design

For each dimension, R = 8 independent replicate Monte Carlo runs at niter = 250,
reduced at every candidate level from the *same* draws (`01i_nn_multiquant_table.R`,
`engine = "orig_list"`, a verbatim transcription of `R/ccds/NN_Dist_Est.R`'s loop
body). The shipped table is then compared **within its own dimension** against
that replicate cloud. This is the fix for `46_verify_wp2a_alpha.R`, which
compared against a *neighbouring dimension* and failed because the
per-unit-dimension effect (0.01734) is 9.4x the alpha effect (0.00185).

Three statistics, in decreasing order of trust:

- **Test 1 (commissioned primary)** — pointwise sign test over subsample size k:
  is the shipped `S_k` nearer the level-A cloud median or the level-B cloud
  median? The generator draws a fresh sample for every k, so the k are
  independent trials. k where the cloud medians coincide are excluded from the
  denominator.
- **Test 1L (added)** — empirical level. Pool the raw draws across replicates
  (M = 2000 per (d,k)) and report `p_k = Fhat_k(S_k)`; `E[p_k]` is the alpha
  that generated `S_k`. Added because Test 1 compares a shipped quantile
  estimate against fresh estimates at niter = 250, and a lower-tail quantile
  estimate is biased upward at small niter in a level-dependent way
  (`quantile(x, 0.001)` from 250 draws is essentially the sample minimum, whose
  expected rank is 1/251 = 0.004, not 0.001). The shipped tables' niter is not
  recorded, so that bias cannot be matched away; Test 1L does not depend on it.
- **Test 2 (secondary)** — behavioural: rerun UN-MCCD and SUN-MCCD with each
  fresh replicate and ask which level reproduces the published metric.

Generation n: 555 for every dimension except d = 9, generated at n = 1013 to
cover `shuffle`. Tier 2 (d = 12, 16, 21) is n = 555 throughout, which is far
below those data sets' n (1452 / 3200 / 3443) and does not need to be above it:
the level is a property of the table. This deviates from the commissioned per-dimension n (74 at
d = 19, 148 at d = 18, …) upward: the number of usable k *is* the sign test's
sample size, and at 24 cores the extra cost was minutes. The commissioned n are
all ≤ 555, so nothing was lost. Measured cost was ~4x the brief's model
(33 s wall per (d, replicate) at n = 555, 200 s at n = 1013, 24 cores).

## Positive control — the discriminator has power (12/12)

**Re-demonstrated in the same execution as the tier-2 verdicts**, not inherited
from the tier-1 run: `run_test1()` re-derives every control row from the
replicate clouds on each invocation, so the calibration below and the d = 12 /
16 / 21 verdicts above came out of one process.

Three dimensions carry genuinely distinct shipped pairs, so the true label is
known from the filename. Both files at each were run through the identical
machinery. d = 10 and d = 20 control the extreme regime (where d = 18, 19 sit);
d = 5 (90%/95%) was added to control the less-extreme regime (where d = 6..9
sit) — without it the 95-vs-99 verdicts would have rested on a discriminator
validated only at 99-vs-999.

| d | true label | sign-test frac_A | verdict (sign) | level estimate | verdict (level) | correct |
|---|---|---|---|---|---|---|
| 5 | 90 | 0.998 | 90 | 0.1005 | 90 | yes |
| 5 | 95 | 0.000–0.002 | 95 | 0.0505 | 95 | yes |
| 10 | 99 | 0.886–0.912 | 99 | 0.0105–0.0107 | 99 | yes |
| 10 | 999 | 0.0036 | 999 | 0.00149–0.00153 | 999 | yes |
| 20 | 99 | 0.899 | 99 | 0.0103–0.0104 | 99 | yes |
| 20 | 999 | 0.000–0.007 | 999 | 0.00150–0.00151 | 999 | yes |

12 of 12 known-answer cases recovered correctly by **both** statistics, spanning
alpha = 10%, 5%, 1%, 0.1%, on both the `average` and `median` vectors. The level
estimator's readings (0.1005, 0.0505, 0.0104, 0.0015) track the nominal alphas
(0.10, 0.05, 0.01, 0.001) closely; the 0.001 reading runs ~1.5x high, which is
the expected residual of estimating a 1-in-1000 tail from a 2000-draw pool and
is far smaller than the 10x gap to the competing hypothesis.

The control also calibrates the sign test's scale. A **true-99** file produces
frac_A = 0.886–0.912 (not 1.0, because the niter = 250 999-cloud is biased
toward the 99-cloud); a **true-999** file produces frac_A = 0.000–0.007. These
two signatures are far apart, and the unknown dimensions fall squarely on one
of them.

## Per-dimension verdict

| d | candidates | n informative k | frac_A | sign p | level estimate | **verdict** |
|---|---|---|---|---|---|---|
| 6 | 95 / 99 | 554 | 1.000 | 3.4e-167 | 0.0504–0.0507 | **95% (alpha = 5%)** |
| 7 | 95 / 99 | 554 | 1.000 | 3.4e-167 | 0.0504–0.0505 | **95% (alpha = 5%)** |
| 8 | 95 / 99 | 554 | 1.000 | 3.4e-167 | 0.0497–0.0498 | **95% (alpha = 5%)** |
| 9 | 95 / 99 | 1012 | 0.999–1.000 | 4.6e-302 | 0.0503–0.0506 | **95% (alpha = 5%)** |
| 18 | 99 / 999 | 554 | 0.908–0.917 | 1.6e-99 | 0.0103–0.0108 | **99% (alpha = 1%)** |
| 19 | 99 / 999 | 554 | 0.904–0.912 | 1.7e-92 | 0.0105 | **99% (alpha = 1%)** |
| 12 | 99 / 999 | 554 | 0.890–0.904 | 5.6e-85 | 0.0102–0.0105 | **99% (alpha = 1%)** |
| 16 | 99 / 999 | 554 | 0.904–0.908 | 1.7e-92 | 0.0103 | **99% (alpha = 1%)** |
| 21 | 99 / 999 | 554 | 0.892 | 6.9e-86 | 0.0103–0.0104 | **99% (alpha = 1%)** |

Both statistics agree at every dimension, on both the `average` and `median`
vectors.

### The level estimate with spread and ratio to both candidates

This is what makes the verdicts identifying rather than directional. `level` is
the mean of `Fhat_k(S_k)` over independent k; the CI is +/- 1.96 SE of that
mean. The per-k IQR is also given, but it describes the *grain of a single
column* — at alpha = 0.001 with M = 2000 pooled draws `phat_k` can only land on
multiples of 1/2000 — not the precision of the estimate.

| d | level (avg / med) | 95% CI (avg) | per-k IQR | x alpha_A | x alpha_B | verdict |
|---|---|---|---|---|---|---|
| **controls** | | | | | | |
| 5 (true 90) | 0.10047 / 0.10065 | [0.0997, 0.1012] | 0.094–0.107 | **1.00x** of 0.10 | 2.01x of 0.05 | 90 correct |
| 5 (true 95) | 0.05047 / 0.05080 | [0.0499, 0.0510] | 0.046–0.055 | 0.50x of 0.10 | **1.01x** of 0.05 | 95 correct |
| 10 (true 99) | 0.01052 / 0.01072 | [0.0102, 0.0108] | 0.0085–0.0125 | **1.05x** of 0.01 | 10.52x of 0.001 | 99 correct |
| 10 (true 999) | 0.00149 / 0.00153 | [0.00138, 0.00159] | 0.0005–0.0020 | 0.15x of 0.01 | **1.49x** of 0.001 | 999 correct |
| 20 (true 99) | 0.01037 / 0.01032 | [0.0101, 0.0106] | 0.0080–0.0125 | **1.04x** of 0.01 | 10.37x of 0.001 | 99 correct |
| 20 (true 999) | 0.00151 / 0.00150 | [0.00140, 0.00162] | 0.0005–0.0020 | 0.15x of 0.01 | **1.51x** of 0.001 | 999 correct |
| **unknown** | | | | | | |
| 6 | 0.05071 / 0.05039 | [0.0501, 0.0513] | 0.046–0.056 | **1.01x** of 0.05 | 5.07x of 0.01 | **95** |
| 7 | 0.05043 / 0.05048 | [0.0499, 0.0510] | 0.046–0.055 | **1.01x** of 0.05 | 5.04x of 0.01 | **95** |
| 8 | 0.04971 / 0.04975 | [0.0491, 0.0503] | 0.045–0.055 | **0.99x** of 0.05 | 4.97x of 0.01 | **95** |
| 9 | 0.05026 / 0.05059 | [0.0498, 0.0507] | 0.046–0.055 | **1.01x** of 0.05 | 5.03x of 0.01 | **95** |
| 12 | 0.01021 / 0.01054 | [0.0099, 0.0105] | 0.0080–0.0120 | **1.02x** of 0.01 | 10.21x of 0.001 | **99** |
| 16 | 0.01025 / 0.01031 | [0.0100, 0.0105] | 0.0080–0.0120 | **1.03x** of 0.01 | 10.25x of 0.001 | **99** |
| 18 | 0.01076 / 0.01028 | [0.0105, 0.0110] | 0.0085–0.0125 | **1.08x** of 0.01 | 10.76x of 0.001 | **99** |
| 19 | 0.01053 / 0.01051 | [0.0103, 0.0108] | 0.0085–0.0125 | **1.05x** of 0.01 | 10.53x of 0.001 | **99** |
| 21 | 0.01028 / 0.01043 | [0.0100, 0.0106] | 0.0080–0.0120 | **1.03x** of 0.01 | 10.28x of 0.001 | **99** |

Every unknown dimension sits within 1–8% of one candidate and a factor of 5 or
10 from the other, with a CI far too tight to reach the alternative. The
controls establish that this reading is accurate to the same tolerance at all
four levels actually in use.

### d = 21 corroborated on an independent Monte Carlo

d = 21 got a second, independent check, because it is the dimension with the
most history attached to it. `57_waveform_d21_rerun.R`'s generation left a raw
draw matrix at **n = 5000** in
`R/NN-test_quantile_d21_regen/samedraws/NN-draws_21d_n5000.RData` — a different
script, different seeds, and the `fast_stream` engine rather than `orig_list`.
Splitting it into blocks reproduces this study's design on entirely different
data, and it reaches k = 5000 instead of k = 555.

| source | engine | k used | level (avg / med) | 95% CI (avg) | x 0.01 | x 0.001 |
|---|---|---|---|---|---|---|
| this study, rep01–08, n=555 | orig_list | 554 | 0.01028 / 0.01043 | [0.01001, 0.01056] | 1.03x | 10.28x |
| d21_regen samedraws, n=5000 | fast_stream | 4999 | 0.01039 / 0.01052 | [0.01020, 0.01057] | 1.04x | 10.39x |

Two independent Monte Carlos, two engines, two values of n, a 9x difference in
the number of k — same answer. Only the auxiliary *level* estimate is
interpretable; its per-block reduction is niter = 50, which makes a 0.999
quantile the sample minimum, so the auxiliary sign test's B-cloud is badly
biased and its `frac_A` is not comparable to the main runs'. That is recorded
in the CSV rather than quietly dropped. The d = 18/19 sign-test fractions (0.90–0.92) match the control's
**true-99** signature (0.886–0.912) and are nowhere near its true-999 signature
(0.000–0.007). The containment counts tell the same story: at d = 18/19,
`n_inside_A_only` = 233–258 against `n_inside_B_only` = 25–30, the same pattern
as the known-99 controls (231–251 vs 28–37) and the opposite of the known-999
controls (0–1 vs 397–419).

### The "less extreme survives" pattern held, 9 of 9

After tier 1 I flagged, as directional and untested, that in all six resolved
dimensions the surviving content was the **less extreme** of the two levels.
Tier 2 was the test, and it was a real one: d = 12, 16 and 21 could each have
come back 0.001. **All three came back 0.01.** The pattern now holds in 9 of 9
duplicated dimensions — the duplicate always carries the less extreme level's
numbers under the more extreme filename.

This is consistent with the duplication being a save-under-two-names of a
single run, but it remains an empirical regularity across the nine dimensions
that carry data sets, not a proven property of the generation process. The
duplicated dimensions with no data set attached (d = 11, 13, 14, 15, 17, 22–28)
were not tested and are not claimed.

## Per data set: actual alpha vs the paper's claim

Table `tab:alpha_real` claims: d 5–9 → UN 5%, SUN 5%; d 10–19 → UN 1%,
SUN 0.1%; d ≥ 20 → all 0.1%. Both UN-MCCD and SUN-MCCD load the same file at
d ≤ 9 (`nn_quant_for_d` → "95"); at d = 10–19 UN loads `_99%` and SUN loads
`_999%`, but those two files are identical in content; at d ≥ 20 both load
`_999%`.

| data set | d | file content is | UN-MCCD actual | claim | SUN-MCCD actual | claim | |
|---|---|---|---|---|---|---|---|
| vertebral | 6 | 5% | 5% | 5% | 5% | 5% | OK |
| ecoli | 7 | 5% | 5% | 5% | 5% | 5% | OK |
| pima | 8 | 5% | 5% | 5% | 5% | 5% | OK |
| glass | 9 | 5% | 5% | 5% | 5% | 5% | OK |
| WBC | 9 | 5% | 5% | 5% | 5% | 5% | OK |
| stamps | 9 | 5% | 5% | 5% | 5% | 5% | OK |
| shuffle | 9 | 5% | 5% | 5% | 5% | 5% | OK |
| vowels | 12 | 1% | 1% | 1% | **1%** | 0.1% | **MISMATCH** |
| PenDigits | 16 | 1% | 1% | 1% | **1%** | 0.1% | **MISMATCH** |
| lymphography | 18 | 1% | 1% | 1% | **1%** | 0.1% | **MISMATCH** |
| hepatitis | 19 | 1% | 1% | 1% | **1%** | 0.1% | **MISMATCH** |
| waveform | 21 | 1% | **1%** | 0.1% | **1%** | 0.1% | **MISMATCH, both** |

**Seven of the twelve data sets match the paper's claim on both methods**
(`vertebral`, `ecoli`, `pima`, `glass`, `WBC`, `stamps`, `shuffle` — all at
alpha = 5%, as stated).

Five do not:

- `vowels`, `PenDigits`, `lymphography`, `hepatitis` — **SUN-MCCD** ran at
  alpha = 1%, not the 0.1% the table states. UN-MCCD is correct on all four,
  because at d = 10–19 UN loads the `_99%` filename, which is the one whose
  content matches its label.
- `waveform` (d = 21) — **both methods are misstated.** At d >= 20 both
  `nn_quant_label_paper_UN` and `nn_quant_label_paper_SUN` resolve to `"999"`,
  so both load `NN-test-simul_21d_999%.RData`, and that file holds the 1%
  table. The published waveform rows for UN-MCCD and SUN-MCCD were both
  produced at alpha = 1%, against a claimed 0.1%. This is the only data set
  where the UN-MCCD claim also fails.

In every case this is a labelling error, not a different experiment: the number
that was published is the number the 1% table produces. Nothing needs
recomputing; Table `tab:alpha_real` needs correcting.

### d = 21 specifically: this disagrees with the earlier suggestion

`57_waveform_d21_rerun.R` left the d = 21 level "suggestive, hedge it" at
**0.001**, inferred behaviourally from F2 agreement to 3 dp between the shipped
run and a fresh 999 run. **That suggestion is wrong, and this study contradicts
it directly.** The shipped d = 21 table holds the **1%** quantile: 0.01028
(average) and 0.01043 (median), CI [0.01001, 0.01056], which is 1.03x alpha =
0.01 and 10.3x alpha = 0.001 — and the independent n = 5000 Monte Carlo agrees
at 0.01039. The positive control at d = 20, the adjacent dimension, recovers
both of its known labels correctly with the same machinery in the same
execution.

The earlier inference failed for a reason worth recording: **F2 agreement to
3 dp does not imply the same confusion matrix.** F2 is a coarse summary, two
different (TP, FP, FN) triples can round to the same F2 at 3 dp, and the
behavioural test in this study reproduces exactly that failure mode at scale —
52 of 72 tier-1 cells came out "consistent with both". Behavioural agreement
was never going to identify the level; reading the level off the table is.

## Test 2 (behavioural) — mostly no power, as expected

### Truncation, and the second baseline it forced

The tier-2 data sets are much larger than the fresh tables (vowels 1452,
PenDigits 3200, waveform 3443, against 555). This does not touch Test 1, but it
does distort Test 2: `nnccd.radi`'s `pmin(1:n, length(...))` clamp
(`harness.R:154-158`) silently reuses the last table entry past the end, so a
fresh run on waveform uses 555 entries where the shipped run uses 5000.
Comparing those two directly would measure table **length**, not alpha — the
same category of error that sank `46_verify_wp2a_alpha.R`.

So for every data set whose n exceeds the fresh table length, a second baseline
is run: the shipped table **truncated in memory to the fresh length**. Fresh-A,
fresh-B and `shipped_trunc` then sit under an identical clamp and are
comparable; the untruncated `shipped` row is kept purely as the reproduction
anchor against the published number. For those data sets the summary's
comparison target is `shipped_trunc`, not the published value — recorded in the
`target_is` column so nothing is silently swapped.

### Results

Detector runs: 9 tier-1 data sets and `vowels` at {shipped, 8 fresh level-A,
8 fresh level-B}; `PenDigits` and `waveform` at 4 replicates per level instead
of 8, because a cell there costs ~5x what it does on vowels and the full grid
would have run ~3 h for a test that is secondary by construction. All with
`min.cls = 0.05` for SUN-MCCD and metrics via the harness's `evaluate()`.

- **Reproduction anchor: the shipped-table run reproduces the published value in
  all 80 (data set x method x metric) cells** to 3 dp, `vowels` included. The
  pipeline is sound.
- Of 80 verdict cells: 60 "consistent with both", 11 "NO POWER (A and B give
  identical values)", 6 pointing at 95, 3 pointing at 99. The 6-vs-3 split is
  not coherent — `vertebral` UN-MCCD F2 points at 99 while `vertebral` SUN-MCCD
  BA points at 95 on the same data set — so Test 2 is noise here and contributes
  nothing beyond the anchor. This was the expected outcome: Test 2 can only ask
  which fresh level reproduces the shipped table's *behaviour*, which Test 1
  asks of the table directly and far more sharply.
- **`vowels`: all 8 cells "consistent with both".** No power, exactly as
  predicted. The B-range strictly contains the A-range on 6 of 8 metrics
  (e.g. SUN-MCCD BA: A = [0.658, 0.770], B = [0.658, 0.941], target 0.7704),
  which is the geometry that makes containment tests uninformative here.
- **`PenDigits` and `waveform` Test 2 were still running when this was
  written.** Measured cost is worse than the vowels extrapolation suggested:
  the first PenDigits UN-MCCD cell alone ran ~6 min (single-threaded and fully
  CPU-bound — 298 CPU-seconds in 300 s of wall time), which puts the pair at
  roughly 3-4 h at 4 replicates per level. Rows land per cell in
  `wp2a_provenance_behavioural.csv` and `--summarize2` can be re-run against
  whatever has landed. **They do not bear on any verdict** — those come from
  Test 1, which is complete. Given that `vowels` returned "consistent with
  both" on all 8 cells, the expected outcome here is the same.

## Explicit NO POWER statements

- **Test 2 has no power at all on:** `lymphography` SUN-MCCD (all four metrics),
  `WBC` UN-MCCD (all four metrics), and TPR alone on `glass` SUN-MCCD, `WBC`
  SUN-MCCD and `shuffle` SUN-MCCD — 11 cells where the level-A and level-B
  replicate runs give identical metrics, so nothing about alpha can be read from
  them. `lymphography` SUN-MCCD and `WBC` UN-MCCD are alpha-invariant on this
  data outright; the three TPR-only cells are at TPR = 1 under every table.
- **Test 2 is effectively uninformative on the remaining 60 cells** ("consistent
  with both"), including all 8 `vowels` cells: the comparison target sits inside
  both replicate ranges.
- **Test 1 has power everywhere in this tier.** `n_informative` is the full k
  range at every dimension (554, or 1012 at d = 9); no dimension collapsed.
- The stricter "clean separation" subset (clouds not overlapping at all) is
  465–492 of 554 k at d = 6–8 and 872–886 of 1012 at d = 9, but only 43–52 of
  554 at d = 18/19 — at the 99/999 pair the niter = 250 clouds overlap at most
  k. The d = 18/19 verdicts therefore rest on the cloud-median comparison and
  the level estimator, not on non-overlap. Both are validated by the control at
  exactly those levels.

## Caveats, assumptions, and what surprised me

1. **The shipped tables' niter is unknown.** The generation scripts left in
   `R/NN-test_quantile/*.R` use `n = 1000, iteN = 10000`, but every shipped
   table is length 5000 — so those scripts are *not* what produced the shipped
   d = 6..28 tables, and the real niter is unrecorded. This is why Test 1L
   exists; it is the statistic that does not depend on matching niter.
2. **Assumption, stated in the brief and verified in code:** a table entry at
   subsample size k does not depend on the n the table was generated at.
   `R/ccds/NN_Dist_Est.R:24` draws a fresh sample of size x for every x, and
   column x is reduced from that sample alone. This is what licenses comparing
   n = 555 fresh tables against n = 5000 shipped ones.
3. **The alpha effect is concentrated at small k.** At d = 10 the 99/999 gap is
   0.132 at k = 2, 0.011 at k = 100, and 0.0018 at k = 1000. Most of the
   discriminating power lives in the first ~200 columns; the tail contributes
   little. That is also, in hindsight, exactly why attempt 46 failed.
4. **Consistent direction of the duplication — RESOLVED, and the tier-1
   extrapolation was correct.** After tier 1 I predicted from six cases that
   d = 12, 16 and 21 would also carry the less extreme level, and flagged it as
   directional. Tier 2 confirmed all three. Recorded here because the
   prediction preceded the test, not because it makes the untested dimensions
   (d = 11, 13, 14, 15, 17, 22–28) safe to assume — they remain unclaimed.
5. **Not settled:** which of the two identical files was the original and which
   the copy. The content identifies the alpha; it does not identify the write
   order. Nothing downstream depends on this.
6. **The d = 5 control is at 90%/95%, not 95%/99%.** It is the only known-label
   case available below d = 10. It validates that the estimator reads the right
   alpha in the less-extreme regime, but it does not literally exercise the
   95-vs-99 discrimination. Given that the level estimates at d = 6–9 come out
   at 0.0497–0.0507 against a 0.05-vs-0.01 choice — a 10x gap, and the d = 5
   control returned 0.0505 for a known 95% file — the residual risk is small.


> **Superseded 2026-09-05.** The "tab:alpha_real needs correcting" conclusion above predates the table repair (`71`/`72`/`75`, 2026-08-11 to 08-13). The tables were regenerated at the stated levels and the six cells rerun, so the table is correct as printed. Re-verified row by row on 2026-09-05: 32/32 NND cells and 16/16 RK cells match the measured level of the file each published number came from. Record of truth: `revision_experiments/INSTALLED_REGEN999_TABLES.md`.
