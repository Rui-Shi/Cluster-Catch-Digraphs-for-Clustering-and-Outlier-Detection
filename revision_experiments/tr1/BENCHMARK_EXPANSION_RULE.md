# Benchmark expansion rule — declared 2026-08-09, before any result was computed

The revision expands the Section 6 real-data study beyond the original eight
data sets, to answer AE.3 (narrow evaluation) and AE.4 (weak high-dimensional
validation).

To keep the expansion defensible, the inclusion rule is fixed **in advance** and
depends only on structural properties. It was written and committed before any
detector was run on any new data set. The structural inventory it draws on
(`32_dataset_inventory.R` → `results/tr1/dataset_inventory.csv`) reports sample
size, dimension, contamination and sort correctness, and nothing else — no
method was executed to produce it.

## The rule

Include every outlier-detection benchmark already constructed by
`data/outlier_detection/RealData_Collection.R`, subject to:

1. **n ≤ 5000.** CCD construction is O(n²); larger sets are out of scope for
   this revision's compute budget and are handled separately in the scalability
   evaluation.
2. **Contamination ≤ 15%.** Above that the minority class is not meaningfully
   an outlier population.
3. **The required RK and NN quantile tables exist at the data set's
   dimension.** Verified before declaring the rule; all candidates pass.

**Every data set satisfying the rule is reported. No data set may be dropped
after its results are known, and no data set may be added later on the basis of
how it scores.**

## What the rule admits

All sixteen candidates pass all three criteria — the rule performs no selection
at all, which is the strongest position available.

| data set | n | d | n0 | contam % | in original study |
|---|---|---|---|---|---|
| lymphography | 148 | 18 | 6 | 4.1 | |
| glass | 213 | 9 | 9 | 4.2 | yes |
| WBC | 223 | 9 | 10 | 4.5 | |
| vertebral | 240 | 6 | 30 | 12.5 | yes |
| ecoli | 336 | 7 | 8 | 2.4 | yes |
| stamps | 340 | 9 | 31 | 9.1 | yes |
| WDBC | 367 | 30 | 10 | 2.7 | |
| pima | 555 | 8 | 55 | 9.9 | |
| shuffle | 1013 | 9 | 13 | 1.3 | |
| vowels | 1452 | 12 | 46 | 3.2 | yes |
| PenDigits | 3200 | 16 | 20 | 0.6 | |
| waveform | 3443 | 21 | 100 | 2.9 | yes |
| thyroid | 3656 | 6 | 93 | 2.5 | |
| pageblocks | 4795 | 10 | 510 | 10.6 | |
| wilt | 4819 | 5 | 257 | 5.3 | yes |
| hepatitis | 74 | 19 | 7 | 9.5 | yes |

Sixteen data sets, n from 74 to 4819, d from 5 to 30, contamination 0.6% to
12.5%.

**WDBC raises the dimensional ceiling from 21 to 30**, which is the point of the
exercise for AE.4. waveform is retained: it is the case where the MC-SRT was
shown to saturate (identical output at the 99% and 999% levels at d = 21), and
removing the second-highest-dimension set while claiming improved
high-dimensional validation would be indefensible.

## Configuration, also fixed in advance

| component | setting |
|---|---|
| proposed detectors | S_min = 0.0625, MC-SRT α from the paper's per-dimension rule |
| LOF | MinPts 11–30, threshold 1.5 — unchanged; published rows retained for the original eight, run fresh on the new eight |
| DBSCAN | MinPts = 4, eps from the 4-distance curve at the known contamination rate — the published protocol, with the `df`→`X` defect corrected so the label column is no longer an input feature |
| MST | cont = 0.02, thresh = 1.2 applied uniformly to every data set, no per-set override |
| ODIN | defaults, k = round(√n), threshold = round(n^(1/3)) |
| iForest | ntrees = 1000, sample_size = min(256, n), threshold 0.55, seed 1 |

## Known data caveat

`RealData_Collection.R` mis-sorts glass (line 34, sorts by column 9) and ecoli
(line 112, sorts by column 6) — feature columns, where the labels are in
columns 10 and 8. All fourteen other sets sort correctly. This does not affect
any number reported here, because `evaluate()` reorders `(Y, score)` jointly
before counting, but it does invalidate any previously published figure for
those two sets that was counted positionally.
