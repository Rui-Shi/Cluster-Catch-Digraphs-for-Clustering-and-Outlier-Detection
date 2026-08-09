# REGEN_SPEC.md — real-data table regeneration contract

Shared contract for the parallel agents regenerating the Section 6 real-data
tables of NEUCOM-D-26-15191. Read this in full before running anything.

Created 2026-08-08. Supersedes nothing; `13_wp0_gate.R` and `wp0_mccd_methods.R`
are used as-is.

---

## 0. Hard rules

1. **Do not modify** any of: `harness.R`, `wp0_mccd_methods.R`, `13_wp0_gate.R`,
   `published_realdata_truth.csv`, or anything under `methods/`, `R/`, `data/`,
   `simulations/`. They are shared by every agent running right now; an edit
   would corrupt another agent's run. If you believe one of them is wrong,
   **report it, do not fix it.**
2. **Do not touch** any file in the manuscript repo root (`*.tex`, `*.bib`,
   `*.md`, `*.png`). Manuscript edits are the main session's job.
3. Create **new numbered scripts only**, in `revision_experiments/`, using the
   number assigned to you in your brief. Never reuse another agent's number.
4. Write results **only** to the CSV shard named in your brief. Two writers on
   one CSV interleave partial lines and destroy it.
5. **Never `git push`.** Committing inside the nested experiment repo is fine.

## 1. Environment

- R: `C:\Program Files\R\R-4.6.1\bin\Rscript.exe`
- **Invoke R from the PowerShell tool, not Bash** — R segfaults under this
  Bash environment.
- Repo root: `G:/Submissions/TR1/TR1_Neurocomputing_resubmit/Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection`
- Run with the repo root as the working directory; scripts resolve paths via
  `here::here()`.
- 24 cores, 63 GB RAM, but five other agents are running. Keep to **one R
  process at a time** unless your brief says otherwise.

### Running jobs — read this, it has bitten us

**Run every R job synchronously and report what it printed.** Do not launch it
with `run_in_background` and then idle-wait: previous agents in this project
stalled indefinitely doing that and had to be killed.

The PowerShell tool caps at 600 s per call (`timeout` parameter, max 600000 ms).
If a single cell needs longer, **split the work so each invocation fits** — one
dataset per call, or one method per call. Everything in these briefs has been
sized to fit except the `wilt` proposed-method cells, which the main session is
running itself.

If a job genuinely cannot be split under 600 s, stop and report that rather
than backgrounding it.

## 2. Data

Load **only** via `load_real_dataset(name)` from `harness.R`. Do **not** read
`results/datasets_csv/` — those are TR2 exports on a different preprocessing
path and will not match the published tables.

```r
suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
dat <- load_real_dataset("glass")   # list(X, Y, d, n)
```

The eight datasets, as the loader returns them:

| dataset | n | d | n0 (outliers) |
|---|---|---|---|
| hepatitis | 74 | 19 | 6 |
| glass | 213 | 9 | 9 |
| vertebral | 240 | 6 | 30 |
| ecoli | 336 | 7 | 8 |
| stamps | 340 | 9 | 31 |
| vowels | 1452 | 12 | 46 |
| waveform | 3443 | 21 | 100 |
| wilt | 4819 | 5 | 257 |

Confirm n/d/n0 at the top of your run and report any disagreement.

**Label convention: `Y == 1` is a regular point, `Y == 0` is an outlier.**
This is inverted relative to most outlier-detection code. Score polarity is
*larger = more outlying*.

Always score with `evaluate(Y, score, threshold)` from `harness.R`. It jointly
reorders `(Y, score)` regulars-first before calling `count_scores2`, which is
what guards the `glass` loader's mis-sort. Never call `count_scores2`,
`count_DBSCAN`, `count_MST2` or `count_ODIN` directly — they assume positional
ordering and will silently miscount `glass`.

## 3. Output schema

Append with `append_result(csv_path, row)` from `harness.R`. Every row:

| column | meaning |
|---|---|
| `dataset` | one of the eight |
| `method` | `U-MCCD`, `SU-MCCD`, `UN-MCCD`, `SUN-MCCD`, `DBSCAN`, `MST`, `ODIN`, `iForest` |
| `variant` | short slug for the parameter setting, e.g. `default`, `knee-4`, `published-exact` |
| `params` | human-readable full parameter string, e.g. `minPts=4, eps=knee(4-dist)=0.831` |
| `n`, `d` | as loaded |
| `TPR`,`TNR`,`BA`,`F2` | from `evaluate()` |
| `published_TPR`,`published_TNR`,`published_BA`,`published_F2` | from `published_realdata_truth.csv`, or `NA` |
| `max_abs_diff` | max abs difference over the four metrics, or `NA` |
| `match_3dp` | `TRUE` if all four agree to 3 d.p. |
| `t_total` | wall-clock seconds |
| `status` | `ok`, `error`, or `timeout` |
| `note` | anything worth knowing; `NA` otherwise |
| `timestamp` | `format(Sys.time())` |

Wrap each cell in `tryCatch` and write an `error` row rather than aborting the
run. Use `has_result()` keyed on `(dataset, method, variant)` so a restart
skips completed cells.

## 4. Paper-faithful settings for the four proposed methods

Already implemented in `wp0_mccd_methods.R`; **do not re-derive them.** Drive
them through `13_wp0_gate.R`, which is fully parameterized:

```
Rscript 13_wp0_gate.R "<datasets>" "<methods>" "<min_cls_reading>"
```

with two environment variables:

- `WP0_GATE_OUT_CSV` — your CSV shard (**always set this**; the default is
  another agent's file)
- `WP0_GATE_TIMEOUT_SEC` — per-cell timeout, default 600

`min_cls_reading` must be `half_contam` (the manuscript's literal rule,
`0.5 * n0/n`) or `full_contam` (`n0/n`, comparison only). `min.cls` is a
**proportion of n, not a count** — the detector applies `round(min.cls * n)`
internally. Passing a count collapses the shape-adaptive step to nothing.
`min.cls` is ignored for `U-MCCD` and `UN-MCCD`, which use uniform coverage.

The MC-SRT significance level is resolved per dimension inside the wrappers:
`rk_quant_label_paper(d)` = 99% for d<10 else 999%; `nn_quant_label_paper_UN(d)`
= 85/90/95/99/999 at d ≤ 2/4/9/19/else; `nn_quant_label_paper_SUN(d)` = the NN
rule for d<10, forced 999% for d≥10. Threshold is 0.5 against the binarized
mutual-catch-graph output.

## 5. Baseline protocol — what changed and why

Three reviewer points (R3.9, R5.3, AE.3) attack the real-data baseline protocol
as unfair. Two concrete defects were confirmed in the published drivers:

1. **`Real_Data_DBSCAN.R:24` passes `df`, not `X`** — the full data frame
   including the ground-truth label column — into `DBSCAN()`. The published
   DBSCAN row therefore had the label as an input feature. Every other baseline
   driver correctly passes `X`.
2. **DBSCAN's `eps` is read off the 4-distance curve at the true contamination
   rate**, and **MST's `thresh` is tuned per data set** (1.2 for hepatitis,
   glass, vertebral, ecoli, stamps; 1.3 vowels; 1.05 waveform; 1.4 wilt).

The regenerated tables use **label-free settings only**. LOF, ODIN and iForest
were already label-free and are unchanged in configuration.

**LOF is not re-run.** `LOF/LOF.R` uses `L_MinPts=11, U_MinPts=30` and
`Real_Data_LOF_1.5.R` calls it with exactly those bounds and `Thresh = 1.5`,
matching the manuscript. The published LOF rows stand.

## 6. Reference variants

For DBSCAN and MST, also compute the published configuration as a
**reference row** (`variant = "published-exact"`). This is not a candidate for
the paper — it exists to prove that the published numbers came from the leaked
configuration, which is what lets us describe the change honestly in the
response letter.

## 7. Reporting

Finish with a markdown table of every cell you ran: dataset, method, variant,
the four metrics, the published values, and the delta. Then, in prose:

- which cells reproduce and which do not,
- anything that errored or timed out, and the message,
- any assumption you had to make,
- anything you noticed that contradicts this spec.

**Report numbers, not adjectives.** Do not tune anything to make a number look
better; if a setting produces a bad result, that is the result. Do not
summarize by saying a run "succeeded" — paste the table.
