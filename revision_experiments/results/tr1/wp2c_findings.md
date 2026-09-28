# WP2(c) — a label-free rule for `S_min`, and what it costs

Reviewer points **R3.2** and **R5.2**. Scripts `revision_experiments/54_wp2c_smin_rule.R`
(real data, Steps 1–2) and `revision_experiments/55_wp2c_simulation_arm.R` (simulations, Step 3).

---

## 0. What the objection is, stated accurately

`S_min` appears in the main text three times: its definition, the simulation sentence at
L605 ("evaluates SU-MCCD under the same **simulation** settings … with `S_min` set to half the
contamination level"), and a future-work mention. **It is never specified for the real data.**

So the accurate charge is two-part, and only the first part is a methodological defect:

1. the **simulation** rule is oracle-derived — it reads the true contamination, i.e. the labels;
2. the **real-data** value was never documented. It was already label-free (a fixed constant),
   but the paper does not say so.

The paper does not misdescribe itself and this document does not claim it does.

## 1. The structural fact everything below rests on

`min.cls` (= `S_min`) reaches the detectors through exactly one expression:

```
R/ccds/RK_CCD_New.R:113   lenD = length(which(graph$catch >= round(min.cls*n)))   # SU-MCCD
R/ccds/UN_CCD.R:102       lenD = length(which(graph$catch >= round(min.cls*n)))   # SUN-MCCD
```

(`SU-MCCDs.R:11` → `RKCCD_correct_quant` → `rccd.silhouette`; `SUN-MCCD.R:14` → `nnccd.silhouette`.
Verified by grep across the repo: no other consumer is on either call path. U-MCCD and UN-MCCD
take no `min.cls` argument at all, so no `S_min` rule can move them.)

Therefore **`k = round(S_min · n)` is a sufficient statistic for `S_min`**. Two `S_min` values that
round to the same `k` give bit-identical output, by construction and with no experiment needed.
This is used throughout to (a) define plateaus exactly, (b) avoid re-running query points that
share a `k` with a measured grid point, (c) run each replicate's detector once per distinct `k`
and share the result across rules. `load_grid()` in `54_wp2c_smin_rule.R` asserts the property
empirically and would abort on a violation; it never fired.

`S_min` is a **proportion of n**, never a count (`RK_CCD_New.R:110`, "the minimum percentage
accepted as a cluster").

---

## 2. Step 1 — plateau structure on the real data

`results/tr1/wp2c_plateau_table.csv` — one row per (dataset, method, plateau); a plateau is a
maximal contiguous run of measured `S_min` grid points with a **bit-identical flagged set**
(`flagged_idx`). Built from `wp2a_smin_sweep16.csv` (227 ok cells, read-only) +
`wp2a_smin_bigfour.csv` (read-only, partial) + `wp2c_extra_cells.csv` (32 new cells measured
here, at the query points whose `k` was not already on the grid).

Coverage: 25 (dataset, method) pairs over 13 datasets.

| question | answer | column |
|---|---|---|
| all three of 0.0625, `n0/n`, `0.5·n0/n` in the **same** plateau | **13 of 25 pairs** | `flag_fixed`/`flag_contam`/`flag_half` |
| … for **both** methods of a dataset | **5 of 13** — hepatitis, pima, stamps, vowels, WDBC | — |
| does not coincide | ecoli, glass, lymphography, shuffle (SU only), vertebral (SUN only), WBC (SU only); plus pageblocks and PenDigits where the contamination `k` is not affordable to measure | — |
| single macro-cluster at `S_min = 0.0625` | **23 of 25 pairs** | `n_clusters` on the `flag_fixed` plateau |
| single macro-cluster at half-contamination | **12 of 22 measured pairs** | `n_clusters` on the `flag_half` plateau |

The two non-degenerate pairs at 0.0625 are vertebral/SUN-MCCD and WBC/SUN-MCCD (2 clusters each).

## 3. Step 2 — candidate rules on the real data

`results/tr1/wp2c_rule_candidates.csv`. A rule **hits** a (dataset, method) pair when the
flagged set it produces is identical to the reference's. Flagged-set identity is used rather than
plateau bookkeeping because it is exact and grid-independent; `plateau_id` is carried alongside so
the two can be cross-checked. Reference = the oracle `0.5·n0/n` (`hit_vs_half`).

| rank | rule | `rule_param` | `hit_vs_half` | `hit_vs_fixed` | mean BA | mean ΔBA vs oracle | frac 1 cluster |
|---|---|---|---|---|---|---|---|
| 1 | fixed constant | 0.0300 | **17 / 25** | 13 | 0.7008 | −0.0040 | 0.500 |
| 2 | fixed constant | 0.0200 | 16 / 25 | 12 | 0.6810 | −0.0109 | 0.417 |
| 3 | fixed constant | **0.0625** | 14 / 25 | 25 | 0.7094 | +0.0191 | 0.920 |
| 3 | fixed constant | 0.0100 | 14 / 25 | 9 | 0.6873 | −0.0195 | 0.304 |
| 5 | fixed constant | 0.0500 | 13 / 25 | 21 | **0.7287** | **+0.0239** | 0.864 |
| 5 | fixed constant | 0.1000 | 13 / 25 | 21 | 0.7224 | +0.0176 | 0.955 |
| 7 | fixed constant | 0.1500 | 12 / 25 | 21 | 0.7226 | +0.0174 | 1.000 |
| 7 | adaptive, fixed point | ρ=1.00 | 12 / 22 | 20 | 0.7222 | +0.0174 | 1.000 |
| 7 | fixed constant | 0.2000 / 0.3000 / 0.0050 | 12 / 25 | 20 / 21 / 6 | — | — | — |
| 12 | adaptive, one step | ρ=1.00 | 11 / 22 | 19 | 0.7197 | +0.0149 | 0.955 |
| 13 | adaptive, one step / fixed pt | ρ=0.50 | 10 / 22 | 18 | 0.7167 | +0.0119 | 0.909 |
| 15 | adaptive, one step / fixed pt | ρ=0.25 | 9 / 22 | 17 | 0.7174 | +0.0126 | 0.864 |

Oracle reference over the same 25 pairs (22 measured): mean BA 0.7048, mean F2 0.3086,
single cluster in 12 of 22.

### Candidate rule 3 (fraction of the median CCD cluster size) — and why it is unusable

The rule is circular: the partition depends on `S_min`. Both resolutions the brief allows were
implemented and both are recorded in `results/tr1/wp2c_adaptive_real.csv`:

* **seed** — `S_min⁽⁰⁾ = 0`, the detectors' own shipped default (`min.cls = 0`, no size floor).
  It uses no label and no tuning constant.
* **one step** — `S_min = ρ · m₀ / n` where `m₀` is the median macro-cluster size at the seed. Stop.
  (`iter == 1` in the trace.)
* **fixed point** — iterate `S_min⁽ᵗ⁺¹⁾ = ρ · mₜ / n`, testing convergence on `k` (exact, because
  `k` is the sufficient statistic), cycle-detected, capped at 8 iterations.

**The iteration converges in 2–4 steps on all 22 real-data pairs, and in 21 of 22 it converges to
the degenerate point `S_min = ρ`.** The failure is self-reinforcing and structural: raising `S_min`
raises the cluster-size floor, which merges macro-clusters, which raises the median cluster size to
`n`, which sets `S_min = ρ·n/n = ρ` — an absorbing state independent of the data. A typical trace
(`k_path` column, WBC / SU-MCCD, ρ=0.5) is `0 > 56 > 112` and stops: 2 clusters → 1 cluster → 1
cluster forever. The one exception is shuffle / SUN-MCCD at ρ=0.25, which converges to `k = 4`
with 57 clusters.

So at its fixed point the "adaptive" rule is **not adaptive at all** — it is the fixed constant ρ
in disguise, and a very large one. The one-step variant is genuinely adaptive but, as §4 shows,
performs catastrophically in the simulations. **The rule is rejected.** It is not ranked above a
fixed constant here on any average, and the averages in the table above that look competitive
(+0.0119 to +0.0174 ΔBA) are an artefact of it collapsing to a single macro-cluster in 86–100% of
pairs — the same collapse that costs it 0.20 BA in the simulation arm.

---

## 4. Step 3 — the simulation arm (the deliverable)

`results/tr1/wp2c_simulation_arm.csv` (10,608 rows, **all `status == "ok"`**) and
`results/tr1/wp2c_simulation_summary.csv`.

**Design.** Generators lifted verbatim from `09_wp3_synthetic.R:124-195`, which copies them from
`simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/{Uniform,Gaussian}/10d/10d_2cls_n500_cont5%.R`.
Only `n`, `d` and the contamination level are promoted to arguments; every other constant
(2 clusters, `cls_dis = 3`, `otl_dis = 2`, `r_min/r_max = 0.7/1.3`, `noise_level = 0.01`, outlier
radius 5) is byte-identical. Seeds `123 + rep`, recorded in the `seed` column.

| axis | values | maps to |
|---|---|---|
| generator | uniform, gaussian | Section 5's uniform and Gaussian cluster studies |
| n | 200 (main), 500 (anchor) | both in the manuscript's `n ∈ {50,100,200,500,1000}` |
| d | 10 | the drivers' own d |
| contamination | 0.02, 0.05, 0.10, 0.20 | Section 5's contamination robustness check |
| replicates | 100 (84 for the n=500 anchor) | reduced from the published 1000 per the revision plan |

9 settings × 2 methods = **18 cells**. `uniform_d10_n500_c05` is the drivers' exact configuration,
retained as a fidelity anchor so the reduction to n=200 can be checked rather than assumed.

Metrics via `evaluate(Y, score, 0.5)`; `Y == 1` regular, `Y == 0` outlier. Never bypassed.

**Not run.** d = 5 and d = 20 at n = 200 (16 further settings) were designed but not executed.
Measured cost: d=5 ≈ 17 min/setting (SU-MCCD is 9.0 s per call at n=200, d=5, against 4.1 s at
d=10 and 2.8 s at d=20), d=20 ≈ 5.5 min/setting → about **3 h** for the pair at 6 cores.
Also not run: the gaussian n=500 anchor (≈ 30 min).

### 4.1 Paired label-free minus oracle — the number the response letter needs

Per replicate, rule minus `oracle_half_contamination`; mean over paired replicates, sd across
replicates. Pooled over the 18 cells:

| rule | mean ΔBA | worst cell ΔBA | best cell ΔBA | cells worse | cells better | mean ΔF₂ | mean BA | frac 1 cluster | k̂ correct |
|---|---|---|---|---|---|---|---|---|---|
| fixed 0.0300 | **+0.0240** | **0.0000** | +0.2409 | **0 / 18** | 4 | +0.0255 | 0.9677 | 0.000 | 0.951 |
| `min.cls = 0` | +0.0238 | −0.0032 | +0.2409 | 1 / 18 | 4 | +0.0252 | 0.9676 | 0.000 | 0.950 |
| fixed 0.0500 | +0.0144 | −0.1446 | +0.2409 | 3 / 18 | 3 | +0.0024 | 0.9581 | 0.040 | 0.922 |
| fixed 0.0625 | −0.0057 | −0.2166 | +0.2131 | 7 / 18 | 3 | −0.0293 | 0.9380 | 0.140 | 0.838 |
| adaptive one-step ρ=0.5 | **−0.2018** | −0.2495 | 0.0000 | **17 / 18** | 0 | −0.4182 | 0.7419 | 0.999 | 0.001 |

Oracle reference: BA 0.9428, F₂ 0.8395, single cluster in 11.7% of runs, correct k̂ in 84.3%.

Worst cell per rule (`setting_id` / `method`):

| rule | worst cell | ΔBA | ΔF₂ | BA there |
|---|---|---|---|---|
| `min.cls = 0` | `uniform_d10_n200_c20` / SUN-MCCD | −0.0032 | −0.0058 | 0.9957 |
| fixed 0.0300 | `uniform_d10_n500_c05` / SU-MCCD | **0.0000** | 0.0000 | 0.9935 |
| fixed 0.0500 | `uniform_d10_n500_c05` / SU-MCCD | −0.1446 (sd 0.1208, se 0.0132, 84 reps) | −0.3790 | 0.8488 |
| fixed 0.0625 | `uniform_d10_n500_c05` / SU-MCCD | −0.2166 (sd 0.0834, se 0.0091, 84 reps) | −0.5625 | 0.7768 |
| adaptive one-step | `uniform_d10_n500_c05` / SUN-MCCD | −0.2495 (sd 0.0007) | −0.6524 | 0.7501 |

**Headline: not using the labels costs nothing.** Two label-free rules (fixed 0.03, `min.cls = 0`)
are on average **better** than the oracle, by +0.0240 and +0.0238 BA respectively, and fixed 0.03
is never worse in any of the 18 cells. The oracle's own weakness is at high contamination: at 20%
it sets `S_min = 0.10`, which collapses SU-MCCD to one macro-cluster in 99–100% of runs and costs
it 0.24 BA against any small constant (`uniform_d10_n200_c20`, SU-MCCD: `min.cls = 0` scores
ΔBA = **+0.2409**, sd 0.0351, 100 reps). Reading the true contamination is therefore not merely
unnecessary here, it is actively harmful once contamination is high.

The adaptive rule is the opposite of a close call: −0.2018 BA mean, worse in 17 of 18 cells,
never better in any.

### 4.2 The cluster-count question

True cluster count is **2 in every synthetic row** — both generators build `data1`, `data2` and
outliers. There is **no ground-truth-k column** in any output file and none is needed for the
synthetic arm; k̂ accuracy is `mean(n_clusters == 2)`, reported above and in §4.1. The real-data
files carry no ground-truth cluster count, so a WP7 real-data k̂ study would need its own labels.

Also, to prevent a recurrence of a misreading: **`k` is `round(s_min · n)`**, the integer
cluster-size floor — *not* a cluster count. 16 of 10,608 rows fail a naive re-derivation
`k == round(s_min*n)`; all 16 are `adaptive_one_step_rho0.5` rows whose product lands exactly on a
`.5` boundary (31.5, 26.5, 32.5, 36.5, 37.5, 28.5), where R's half-to-even rounding plus a 1-ULP
shift in the re-read `s_min` flips the result. The recorded `k` is the one the detector used.

Collapse fraction by rule, over 1,768 runs each (`frac_ncls_1` in the report):

| rule | mean `S_min` | mean `k` | frac `n_clusters == 1` | frac == 2 (correct) |
|---|---|---|---|---|
| fixed 0.0300 | 0.0300 | 6.9 | **0.000** | 0.952 |
| `min.cls = 0` | 0.0000 | 0.0 | **0.000** | 0.951 |
| fixed 0.0500 | 0.0500 | 11.4 | 0.036 | 0.926 |
| oracle | 0.0442 | 9.5 | 0.117 | 0.843 |
| fixed 0.0625 | 0.0625 | 13.8 | 0.135 | 0.843 |
| adaptive one-step | 0.2450 | 55.9 | 0.999 | 0.001 |

**The collapse is caused by the `S_min` value, not by the method.** At `S_min = 0.03` or
`min.cls = 0`, zero of 1,768 runs collapse, for either detector. This matters for WP7: the
degeneracy is a parameter choice, not evidence that shape-adaptive coverage is broken.

By method, at the published `S_min = 0.0625` (884 runs each):

| method | frac `n_clusters == 1` | frac == 2 |
|---|---|---|
| SU-MCCD | **0.269** | 0.687 |
| SUN-MCCD | **0.000** | 1.000 |

And by setting, `fixed_0.0625` only:

| setting | `k` | SU-MCCD frac 1 cluster | SUN-MCCD frac 1 cluster |
|---|---|---|---|
| `uniform_d10_n200_c02 … c20` | 12 | 0.02 – 0.08 | 0.00 |
| `gaussian_d10_n200_c02 … c10` | 12 | 0.30 – 0.32 | 0.00 |
| `gaussian_d10_n200_c20` | 12 | 0.54 | 0.00 |
| `uniform_d10_n500_c05` | 31 | **0.89** | 0.00 |

**So the hypothesis that the collapse is confined to the real data is half right, and the half
that is wrong is the important half:**

* **SUN-MCCD — confirmed.** 0 of 884 synthetic runs collapse at 0.0625, against 11 of 13 real
  datasets. For the headline method the collapse really is a real-data phenomenon.
* **SU-MCCD — refuted.** It collapses in 26.9% of synthetic runs at 0.0625 overall, and in 89% at
  n = 500 — the configuration closest to the paper's own driver. The trend is monotone in
  `k = round(0.0625·n)`: 12 → 2–8% (uniform), 12 → 30–54% (gaussian, harder shapes), 31 → 89%.
  Since `k` grows with `n`, this predicts exactly what is seen on the larger real datasets.

---

## 5. Step 4 — recommendation, and what it costs

**Recommend a single fixed label-free constant, `S_min = 0.05`, for both the simulation and the
real-data study.** Justification, with the number behind each clause:

| claim | evidence |
|---|---|
| label-free | a constant; reads nothing from `Y` |
| **no published SUN-MCCD number moves** | identical BA *and* F₂ to 0.0625 on **11 of 11** real datasets; paired ΔBA = 0.0000 |
| one published SU-MCCD number moves, upward | 10 of 11 identical; vertebral changes 0.3881 → **0.4929** BA (2 clusters instead of 1) |
| beats the oracle in the simulations | mean paired ΔBA **+0.0144** over 18 cells, mean ΔF₂ +0.0024 |
| strictly better than 0.0625 on both arms | real data mean BA 0.7287 vs 0.7094; simulations +0.0144 vs −0.0057; collapse 4.0% vs 14.0% |
| less degenerate | correct k̂ in 92.2% of synthetic runs vs 83.8% at 0.0625 |

**Where it loses — flagged, not buried.** `S_min = 0.05` is worse than the oracle in 3 of 18
cells. The worst is `uniform_d10_n500_c05` / **SU-MCCD**: ΔBA = **−0.1446** (sd 0.1208, se 0.0132,
84 paired reps), ΔF₂ = −0.3790, caused by SU-MCCD collapsing to one macro-cluster in 57% of those
runs. That is the largest-n setting and the one nearest the paper's own driver configuration, so
it cannot be dismissed as a corner. Two consequences:

* The loss is confined to **SU-MCCD**, which the manuscript already frames as a subordinate
  construction, not the recommendation. SUN-MCCD's ΔBA in that same cell is exactly 0.0000.
* If a rule that never loses is wanted, it is **fixed 0.03** (ΔBA ≥ 0 in all 18 simulation cells,
  0 collapses in 1,768 runs). It is **not** recommended, because it costs the headline method
  −0.0505 mean BA on the real data and changes the flagged set on 4 of 11 SUN-MCCD datasets —
  including glass (0.7696 → 0.4935) and pima (0.7215 → 0.5333). Buying a clean simulation table
  with a 0.28 BA loss on a published real-data row is a bad trade, and it would move published
  numbers.

**Would any published number move under the recommendation?** One: vertebral / SU-MCCD, from
BA 0.3881 to 0.4929. Every SUN-MCCD real-data number is bit-identical.

**What the response letter can say.** The simulation `S_min` is replaced by a fixed constant that
reads nothing from the labels; across 9 synthetic settings × 2 detectors × 100 replicates the
label-free constant is on average **better** than the oracle rule it replaces (ΔBA = +0.0144), and
the oracle is in fact harmful at 20% contamination (−0.24 BA for SU-MCCD against any small
constant). The real-data study was already label-free; it is now documented as such, and the
constant is the same one in both arms.

---

## 6. Files

| file | rows | contents |
|---|---|---|
| `54_wp2c_smin_rule.R` | — | plateaus, query resolution, extra cells, adaptive traces, candidate ranking |
| `55_wp2c_simulation_arm.R` | — | generators, the 9-setting arm, summary and report modes |
| `wp2c_plateau_table.csv` | 47 | one row per (dataset, method, plateau) + `flag_fixed`/`flag_contam`/`flag_half` |
| `wp2c_rule_candidates.csv` | 432 | one row per (dataset, method, rule); `hit_vs_half`, `hit_vs_fixed`, `BA_minus_oracle` |
| `wp2c_extra_cells.csv` | 32 | query-point cells measured here; schema matches `wp2a_smin_sweep16.csv` |
| `wp2c_adaptive_real.csv` | 188 | fixed-point traces, `k_path`, `converged`, `cycled` |
| `wp2c_simulation_arm.csv` | 10,608 | one row per (setting, rep, method, rule); all `status == "ok"` |
| `wp2c_simulation_summary.csv` | 108 | per (setting, method, rule) with `dBA_mean`, `dBA_sd`, `dBA_se` |

Read-only inputs, never written: `wp2a_smin_sweep16.csv`, `wp2a_smin_bigfour.csv`,
`regen_smin_grid.csv`, `dataset_inventory.csv`, `harness.R`, `wp0_mccd_methods.R`,
`R/ccds/*`, `methods/outlier_detection/*`. No detector code was modified; every rule is applied
by choosing the `min.cls` value passed in.
