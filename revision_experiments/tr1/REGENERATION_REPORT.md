# Real-data regeneration — final report

Completed 2026-08-09. Sixteen data sets × nine methods, one pipeline, one
documented configuration. Supersedes the published Section 6 tables.

Configuration frozen in `BENCHMARK_EXPANSION_RULE.md`; assembled by
`34_final_summary.R` into `results/tr1/final_comparison.csv`.

---

## 1. Defects found in the published tables

Each is evidenced, not inferred.

**1.1 The loader mis-sorts glass and ecoli.** `RealData_Collection.R:34` sorts
glass by column 9 and `:112` sorts ecoli by column 6 — both feature columns,
where the labels sit in columns 10 and 8. The other fourteen sets are correct.

The published counting code (`count_scores2`, `count_DBSCAN`, `count_MST2`,
`count_ODIN`) slices positionally — it calls the last n0 rows the outliers and
never reads Y. On glass the true outliers sit at rows 182–190 of 213; on ecoli
at rows 26, 60, 116, 119, 123, 125, 141, 163. **In both cases the last n0 rows
contain zero true outliers.**

This is *not* a TPR ceiling of 0 — an earlier draft of this report said so and
was wrong, and the published glass and ecoli TPRs of 0.778 and 0.500 refute it
directly. A positional count reports the flag rate over the last n0 positions,
and a detector can flag those rows freely. The actual consequence: the
published "TPR" for these two sets is the flag rate on a block of **regular**
points, and the published "TNR" is computed over a block that still contains
**every** true outlier. Neither is the rate its heading names.

Two independent proofs:

- ecoli: scoring our own scores positionally returns TPR 0.750 on all four
  detectors — exactly the published value — while correct counting returns
  0.875. ecoli UN-MCCD reproduces the published row *exactly* under positional
  counting: 0.750 / 0.668 / 0.709 / 0.204.
- LOF, recomputed at its unchanged published configuration with only the
  counting changed (`38_lof_original8.R`): 3 of 8 rows reproduce exactly
  (hepatitis, stamps, waveform), and the two largest moves are glass
  (TNR 0.618 → 0.794) and ecoli (TPR 0.500 → 0.750) — the mis-sorted pair.

`harness.R`'s `evaluate()` reorders `(Y, score)` jointly before counting, so no
number in this report is affected.

**1.2 Two published rows are duplicates.** glass U-MCCD and SU-MCCD are
byte-identical in both supplement tables (1.000/0.363/0.681/0.257) though the
detectors give TNR 0.441 and 0.647. On stamps, our SUN-MCCD output matches the
published *SU-MCCD* row on all four metrics, and no S_min in [0.01, 0.20]
produces the published SU-MCCD row.

**1.3 DBSCAN had two stacked label leaks.** `Real_Data_DBSCAN.R:24` passed
`df` — features plus the ground-truth column — and read `eps` at the true
contamination rate. Reproducing that exact configuration matches the published
row on 4 of 8 sets; the two failures are glass and ecoli, the mis-sorted pair.

**1.4 MST was tuned per data set**, `thresh` ranging 1.05–1.4. Reproducing the
per-set values matches 7 of 8.

**1.5 Two published rows are internally inconsistent.** The wilt DBSCAN row is
arithmetically impossible: BA 0.673 and F₂ 0.381 cannot follow from TPR 0.000
and TNR 0.959. The vertebral LOF row reports TPR 0.033, TNR 0.938, BA 0.488,
but those rates give 0.485. Recomputation resolves it in favour of the BA —
the true TNR is 0.943 (198/210, not 197/210), from which 0.488 follows exactly.

`37_arith_audit.R` checks BA = (TPR+TNR)/2 across all 144 cells of the current
table; it holds everywhere to rounding.

**1.6 No real-data driver exists for the four proposed methods.** Every
baseline has one under `Algo_Compare_OutlierDetection/*/Simulation/Real_Data/`.
The four MCCD trees contain only `Simulation/`. Nothing in the linked
repository produces the published Section 6 rows for U-MCCD, SU-MCCD, UN-MCCD
or SUN-MCCD.

## 2. S_min recovered

The manuscript specifies S_min only for the simulations (main text L605,
"half the contamination level"); it is never specified for the real data.
Inverting an eight-point sweep against the published rows
(`25_smin_match.R`) gives **S_min = 0.0625** as the unique value in
[0.01, 0.20] reproducing every reproducible published cell.

vertebral pins it exactly: SU-MCCD requires the cluster floor ≥ 15 points and
SUN-MCCD requires ≤ 15, so `round(S_min × 240) = 15`, i.e. S_min = 15/240.

A single fixed constant is also what WP2c needs — it uses no label
information, answering the fairness objection three reviewers raised against
deriving S_min from the contamination rate.

Note the simulation drivers actually use `min.cls = max(0.04, cont/2)`; the
0.04 floor is absent from the manuscript's description.

## 3. Reproduction status, original eight

At S_min = 0.0625: **12 cells reproduce exactly, 5 more within one
observation** (max diff 0.004).

| | U-MCCD | SU-MCCD | UN-MCCD | SUN-MCCD |
|---|---|---|---|---|
| hepatitis | 0.029 | exact | exact | exact |
| glass | 0.078 | dup row | exact | 0.001 |
| vertebral | 0.034 | exact | exact | 0.001 |
| ecoli | pub. wrong | pub. wrong | pub. wrong | pub. wrong |
| stamps | exact | dup row | exact | 0.002 |
| vowels | exact | 0.056 | exact | 0.113 |
| waveform | 0.100 | 0.360 | 0.320 | 0.300 |
| wilt | exact | 0.004 | 0.001 | exact |

Unexplained after exhausting every hypothesis: waveform (all four), vowels
SU/SUN, and U-MCCD on hepatitis, glass, vertebral.

## 4. waveform — three hypotheses tested, all refuted

- **Counting**: waveform sorts correctly, so 1.1 does not apply.
- **α**: UN-MCCD and SUN-MCCD give **bit-identical** results at the 99% and
  999% NN tables. At d = 21 the MC-SRT is *saturated* — the significance level
  no longer moves a single decision. This is a substantive finding for AE.4:
  the manuscript dates high-dimensional degradation to d ≥ 50.
- **Preprocessing**: MADN (the loader default) is already the closest of the
  three. raw and z-score are worse for U-MCCD (0.420 and 0.320 vs 0.100).

Given 1.6 — no driver exists — the published waveform row is likely
unrecoverable.

## 5. Final comparison, 16 data sets

F₂, best per row in **bold**.

| data set | d | U-MCCD | SU-MCCD | UN-MCCD | SUN-MCCD | LOF | DBSCAN | MST | ODIN | iForest |
|---|---|---|---|---|---|---|---|---|---|---|
| wilt | 5 | **0.336** | 0.182 | 0.118 | 0.206 | 0.022 | 0.000 | 0.232 | 0.069 | 0.004 |
| vertebral | 6 | **0.311** | 0.140 | 0.036 | 0.109 | 0.038 | 0.000 | 0.282 | 0.159 | 0.000 |
| thyroid | 6 | 0.293 | 0.327 | 0.261 | 0.389 | 0.172 | 0.412 | 0.138 | 0.093 | **0.664** |
| ecoli | 7 | 0.235 | 0.321 | 0.238 | 0.333 | 0.492 | 0.152 | 0.186 | 0.294 | **0.638** |
| pima | 8 | 0.399 | **0.504** | 0.232 | 0.497 | 0.098 | 0.084 | 0.327 | 0.144 | 0.083 |
| glass | 9 | 0.283 | 0.385 | 0.116 | 0.324 | **0.412** | 0.122 | 0.286 | 0.205 | 0.100 |
| stamps | 9 | 0.072 | 0.376 | 0.381 | **0.455** | 0.162 | 0.144 | 0.373 | 0.262 | 0.108 |
| WBC | 9 | 0.316 | 0.316 | 0.318 | 0.350 | 0.543 | 0.652 | 0.354 | 0.342 | **0.656** |
| shuffle | 9 | 0.224 | 0.179 | 0.183 | 0.162 | **0.389** | 0.317 | 0.065 | 0.243 | 0.174 |
| pageblocks | 10 | **0.570** | 0.525 | 0.100 | 0.526 | 0.218 | 0.545 | 0.378 | 0.078 | 0.343 |
| vowels | 12 | 0.196 | 0.193 | 0.263 | 0.267 | 0.364 | 0.249 | 0.153 | **0.427** | 0.027 |
| PenDigits | 16 | 0.046 | 0.033 | 0.034 | 0.033 | **0.053** | 0.000 | 0.043 | 0.027 | 0.026 |
| lymphography | 18 | 0.682 | 0.732 | 0.600 | 0.638 | 0.667 | **0.862** | 0.224 | 0.532 | 0.750 |
| hepatitis | 19 | 0.278 | 0.286 | **0.446** | **0.446** | 0.000 | 0.000 | 0.375 | 0.313 | 0.122 |
| waveform | 21 | **0.252** | 0.219 | 0.186 | 0.228 | 0.000 | 0.107 | 0.130 | 0.193 | 0.000 |
| WDBC | 30 | 0.170 | 0.131 | 0.146 | 0.172 | 0.441 | 0.111 | 0.242 | 0.286 | **0.472** |

**Proposed methods take 7 of 16.** (Was 8 before LOF was recomputed; the
corrected LOF now takes glass at 0.412, ahead of SU-MCCD's 0.385.)

### Mean across all sixteen

| method | mean F₂ | mean BA |
|---|---|---|
| **SUN-MCCD** | **0.321** | **0.722** |
| SU-MCCD | 0.303 | 0.689 |
| U-MCCD | 0.291 | 0.709 |
| iForest | 0.260 | 0.632 |
| LOF | 0.254 | 0.654 |
| MST | 0.237 | 0.606 |
| DBSCAN | 0.234 | 0.599 |
| UN-MCCD | 0.229 | 0.658 |
| ODIN | 0.229 | 0.619 |

**SUN-MCCD ranks first on both mean F₂ and mean BA**, against honestly
configured baselines, on a benchmark whose membership was fixed before any
result existed. Three of the four proposed methods occupy the top three
positions.

### Dimension

| range | data sets | proposed wins |
|---|---|---|
| d ≤ 10 | 10 | 5 |
| d ≥ 12 | 6 | 2 |

WDBC (d = 30) is the weakest relative showing — F₂ 0.172 against iForest's
0.472, with TNR collapsing to 0.073 for SU-MCCD. Combined with the α
saturation at d = 21, this gives a measured account of where and why the
methods degrade, replacing the vague "d ≥ 50" caveat.

The failure signature is consistent: **TNR collapses while TPR stays at
1.000.** The methods over-flag rather than miss — a cluster-coverage problem,
not a spatial-randomness-test problem.

## 6. Open items

- waveform, vowels SU/SUN, U-MCCD on hepatitis/glass/vertebral: unexplained.
- MST has no good fixed threshold: at 2 it flags nothing on glass and
  waveform; at 1.2 it flags nearly everything once n > 1000 (pageblocks TNR
  0.023). Uniform 1.2 is retained as the honest choice, with the failure
  documented.
- hepatitis n0 is 7, not the 6 recorded in earlier notes. Table 5 should be
  checked.
- LOF was initially retained from the published tables, on the correct grounds
  that its *configuration* was never in doubt (MinPts 11–30, threshold 1.5,
  matching main text L1007). That left it the only transcribed block in the
  table, and transcription was the problem. `38_lof_original8.R` recomputes it
  at the same configuration; 3 of 8 rows reproduce exactly, glass and ecoli
  move for the mis-sort reason, and vertebral/vowels/wilt move slightly.
  Nothing about how LOF is run changed — only how its output is counted.
