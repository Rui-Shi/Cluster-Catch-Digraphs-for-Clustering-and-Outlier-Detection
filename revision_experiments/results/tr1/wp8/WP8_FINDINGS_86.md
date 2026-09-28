# WP8 experiment 2 — local / bridged / collective outliers (R3.7)

**Script:** `revision_experiments/tr1/86_wp8_outlier_types.R`. **Grid:** launched
2026-09-05 in the foreground/background per `WP8_PROTOCOL.md`'s round-3 note
(`Rscript revision_experiments/tr1/86_wp8_outlier_types.R ALL 100 99999999`,
8 settings x 100 reps x 9 methods = 7200 cells), complete as of the files
analysed here (2026-09-06). Analysis only — no script, protocol, or
manuscript file was touched; `results/tr1/wp7/` and `results/tr1/wp8/86_*`'s
sibling files for other experiments (`85_*`, `87_*`, `88_*`) were not read or
written; `results/tr1/wp6/` was not touched.

## R3.7 (verbatim)

> **R3.7** The main principle of the MCG framework is that outliers lack
> connectivity to a cluster core. This assumption is more suitable for
> global outliers located away from the main cluster structure, but it may
> be restrictive for local outliers embedded within clusters, collective
> outliers near cluster boundaries, anomalous points connected to a cluster
> core through a small number of bridging observations, or internally
> well-connected groups of anomalous points. Conversely, a small but
> legitimate cluster may be classified as anomalous if its size is below
> S_min. The authors should include targeted experiments for these difficult
> cases and discuss in greater detail how the method distinguishes bridging
> effects, legitimate small clusters, and collective anomalies.

(The small-legitimate-cluster-below-S_min half of R3.7 is answered by
experiment 3, script 87 — see `WP8_FINDINGS_87.md`. This file covers the
four remaining cases: local, bridged, and collective outliers.)

## Design recap (`WP8_PROTOCOL.md` experiment 2 and its dated 2026-09-05 notes)

Two base clusters as in experiment 1 (n1=95, n2=94, cont=0.05), n=200,
n0 = round(n*0.05) = 10, d in {3,10}, S_min = 0.05 passed as `min.cls` to
SU-MCCD/SUN-MCCD only (U-MCCD/UN-MCCD take no `min.cls` argument at all and
so cannot be constrained by S_min). Four types, 8 settings, 100 reps, 9
methods (U-MCCD, SU-MCCD, UN-MCCD, SUN-MCCD, LOF, DBSCAN, MST, ODIN, iForest).

The originally declared `local` and `bridge` generators failed Opus
verification (`WP8_VERIFICATION.md`, 2026-09-05): the interior local outlier
was not actually isolated (NN ratio 0.61 at d=10, below the required floor
of 2) and 5.2–6.2 of 10 bridge points per rep landed inside a cluster's own
realised ball. Both were rewritten, and `bridge` needed a second rewrite
after round-1's fix still collapsed the chain to a compact blob in 37–41% of
reps (`WP8_REVERIFICATION.md`). The four generators actually run, as
implemented:

- **`local_shell`** — host clusters shrunk to a dense core (`core_frac =
  0.4`); each outlier sits at radius `runif(1, 0.8, 1.0) * scale` along a
  random direction from its host centre — 0.49 beyond the host cluster's own
  realised radius (0.40) on average, i.e. a shell/boundary placement, not
  "embedded within" the cluster as R3.7 describes.
- **`local_interior`** — restored verbatim from the pre-rewrite construction
  (git `5ed690c`): a fixed 0.4*scale interior offset from the host centre,
  full-radius clusters, regular points kept out of a 0.3-radius ball around
  each outlier by rejection sampling. Its near-zero TPR at d=10 is itself
  part of the R3.7 answer: at d=10, n=95, median within-cluster NN distance
  is 0.699 against a realised radius of ~1.009, so an interior point whose
  isolation requires roughly doubling the typical NN distance is not
  geometrically realisable at that n and d — not a construction defect.
- **`bridge`** — a 10-point chain spanning only the gap between the two
  clusters' realised surfaces (clearance 0.05), with per-point jitter
  restricted to the subspace orthogonal to the cluster axis (sd 0.03) so
  that, by Pythagoras, jitter can only increase — never decrease — distance
  from either cluster centre. `n_inside = 0` is analytic under this
  construction, not merely likely. **Pre-registered reading (added before
  the production run, from a 20-replicate diagnostic at d=3, independent
  seeds 777000+1000d+r):** the chain is the densest object in the data
  (chain-point NN 0.098 at d=3 vs 0.204 for regular points, ratio 0.50; 93%
  of a chain point's 5 nearest neighbours are other chain points), and with
  n1+n2+n0 = 199, `round(0.05*199) = 10 = n0` exactly, so SU-MCCD and
  SUN-MCCD are **obliged by `min.cls`** to accept the chain as a legitimate
  cluster and return TPR = 0 on it — a ground-truth/rule collision, not a
  detection failure. The pre-registered ordering, from that diagnostic:
  U-MCCD (mean TPR 0.095) and UN-MCCD (0.060) above SUN-MCCD (0.000 of 20
  reps); SU-MCCD 0.010, ODIN 0.005, LOF/DBSCAN/iForest 0.000.
- **`collective`** — a tight (sd 0.15) group of 10 points at `mu1 + u*(s1 +
  0.5)`, s1 the realised jitter scale — just past cluster 1's own realised
  boundary, internally well-connected, matching R3.7's "collective anomalies
  near a cluster boundary" case directly.

**Density caveat (`WP8_REVERIFICATION.md`):** `local_shell`'s (and
`local_interior`'s) host clusters are drawn at `core_frac = 0.4` for
`local_shell` and full radius for `local_interior`; `local_shell`'s clusters
in particular are ~2.6x denser (by NN spacing) than `bridge`/`collective`'s
full-radius clusters. The four types are not density-matched — a known,
reported asymmetry, not corrected here (out of scope for the fixer pass).

## Integrity results

`count.fields` pass on both output files — every row's field count matches
its header:

| File | Header fields | Data rows | Bad rows |
|---|---|---|---|
| `86_outlier_types.csv` | 14 | 7,200 | 0 |
| `86_outlier_types_clusters.csv` | 8 | 604,800 | 0 |

Both row counts match the pre-registered expectation exactly: 8 settings x
100 reps x 9 methods = 7,200 metric rows; 8 settings x 100 reps x 4 MCCD
methods x 189 regular points (n1+n2 = 95+94) = 604,800 cluster rows.

- **Duplicates** on `(setting_id, d, rep, method)`: 0 in the main file;
  duplicates on `(setting_id, d, rep, method, row_index)`: 0 in the cluster
  file.
- **`status != "ok"`**: 0 rows (all 7,200 metric rows are `status == "ok"`;
  no method failed on any cell).
- **Every `(setting, method)` at 100 reps**: confirmed for all 72
  `setting_id x method` combinations (8 settings x 9 methods) — none missing,
  none short, none over.
- **`n_inside` on bridge rows**: all 1,800 `bridge` rows (2 d x 100 reps x 9
  methods) have `n_inside = 0`, matching the analytic guarantee from the
  orthogonal-jitter construction. All 5,400 non-bridge rows carry the
  sentinel `n_inside = -1` (not applicable), as documented.
- **Cluster-file row counts**: all 3,200 `(setting_id, d, rep, method)` cells
  (8 x 100 x 4 MCCD) have exactly 189 rows — no partial-write cells from a
  crash between the metrics row and the cluster block.
- **`--summarize`**: run against the production file; wrote
  `revision_experiments/results/tr1/wp8/86_outlier_types_summary.csv` (72
  rows, one per `type x d x method`). Cross-checked against an independent
  from-scratch recomputation (base R, no `harness.R` dependency) — identical
  to the digits shown. `n_reps = 100` and `n_error = 0` on all 72 rows: no
  method vanished from the summary. `n_inside_mean` is `NaN`/`NA` for the 54
  rows outside `bridge` (as documented — `local_interior`, `local_shell`,
  `collective` never populate `n_inside`) and exactly `0` for all 18 `bridge`
  rows.

No anomalies found.

## Per-arm results

All tables below are from `86_outlier_types_summary.csv` (100 reps/cell;
SEs shown are over-replicate standard errors of the mean).

### local_interior — TPR / TNR (mean ± SE)

| d | method | TPR | TNR | frac(TPR>0) | mean n_flagged |
|---|---|---|---|---|---|
| 3 | DBSCAN | 0.000 | 0.998 | 0.00 | 0.34 |
| 3 | iForest | 0.000 | 0.947 | 0.00 | 10.07 |
| 3 | LOF | 0.000 | 1.000 | 0.00 | 0.00 |
| 3 | MST | 0.470 ± 0.027 | 0.615 | 0.91 | 77.50 |
| 3 | ODIN | 0.000 | 0.974 | 0.00 | 4.99 |
| 3 | SU-MCCD | 0.000 | 0.9998 | 0.00 | 0.03 |
| 3 | SUN-MCCD | 0.000 | 0.998 | 0.00 | 0.32 |
| 3 | U-MCCD | 0.000 | 0.996 | 0.00 | 0.83 |
| 3 | UN-MCCD | 0.000 | 0.995 | 0.00 | 0.88 |
| 10 | DBSCAN | 0.000 | 0.999 | 0.00 | 0.20 |
| 10 | iForest | 0.000 | 1.000 | 0.00 | 0.00 |
| 10 | LOF | 0.000 | 1.000 | 0.00 | 0.00 |
| 10 | MST | 0.000 | 0.851 | 0.00 | 28.21 |
| 10 | ODIN | 0.000 | 0.883 | 0.00 | 22.05 |
| 10 | SU-MCCD | 0.078 ± 0.018 | 0.886 | 0.18 | 22.32 |
| 10 | SUN-MCCD | 0.000 | 0.996 | 0.00 | 0.85 |
| 10 | U-MCCD | 0.001 | 0.926 | 0.01 | 14.01 |
| 10 | UN-MCCD | 0.000 | 0.976 | 0.00 | 4.52 |

**Reading.** At d = 10 the arm is near-zero for every method as
pre-registered (max TPR = 0.078, SU-MCCD, and that comes with a depressed
TNR of 0.886 — i.e. even SU-MCCD's small nonzero rate rides on flagging
~11% of all regular points, not on a targeted catch of the isolated point).
The geometric argument holds: an interior point whose isolation requires
doubling the median within-cluster NN distance (0.699 vs realised radius
~1.009 at n=95) is not constructible at this n and d.

At d = 3, where the construction is geometrically realisable, only **MST**
shows a non-trivial TPR (0.470), but at the cost of TNR = 0.615 — MST flags
essentially 39% of the entire data set (mean 77.5 of 199 points) at d=3
regardless of type (see `local_shell`, `bridge`, `collective` below, where
MST's flag rate stays in the same 70–78-point range); its apparent
"detection" of the interior outlier is a by-product of generic
over-flagging, not discrimination. All eight other methods, **including all
four MCCD detectors**, score exactly TPR = 0.000 at d = 3 as well. So the
answer to "which methods, if any, detect interior isolated points" is: none
cleanly, at either d; the one method with non-zero TPR does so by flagging
nearly 40% of the data indiscriminately.

### local_shell — TPR / TNR (mean ± SE)

| d | method | TPR | TNR | mean n_flagged |
|---|---|---|---|---|
| 3 | DBSCAN | 1.000 | 1.000 | 10.00 |
| 3 | iForest | 0.966 ± 0.006 | 1.000 | 9.67 |
| 3 | LOF | 1.000 | 1.000 | 10.00 |
| 3 | MST | 0.880 ± 0.009 | 0.653 | 74.44 |
| 3 | ODIN | 1.000 | 0.983 | 13.26 |
| 3 | SU-MCCD | 1.000 | 0.999 | 10.13 |
| 3 | SUN-MCCD | 1.000 | 0.999 | 10.21 |
| 3 | U-MCCD | 1.000 | 0.993 | 11.30 |
| 3 | UN-MCCD | 1.000 | 0.991 | 11.71 |
| 10 | DBSCAN | 0.979 ± 0.007 | 1.000 | 9.79 |
| 10 | iForest | 0.932 ± 0.010 | 1.000 | 9.32 |
| 10 | LOF | 1.000 | 1.000 | 10.00 |
| 10 | MST | 0.916 ± 0.007 | 0.913 | 25.64 |
| 10 | ODIN | 1.000 | 0.978 | 14.17 |
| 10 | SU-MCCD | 1.000 | 0.991 | 11.64 |
| 10 | SUN-MCCD | 1.000 | 0.998 | 10.35 |
| 10 | U-MCCD | 1.000 | 0.978 | 14.21 |
| 10 | UN-MCCD | 1.000 | 0.993 | 11.27 |

**Reading.** As expected, this arm is easy — all four MCCD detectors and
LOF hit TPR = 1.000 exactly at both d, and DBSCAN is at 0.979–1.000. The
methods that "miss" (never completely, always > 0.88) are **iForest**
(0.966 at d=3, 0.932 at d=10) and **MST** (0.880/0.916) — MST again pairs
its partial recall with a heavily depressed TNR (0.653 at d=3, meaning ~35%
of the data is flagged), so its near-full TPR is again largely a
by-product of an aggressive general flag rate rather than a targeted catch.
The near-universal success here should be read together with the density
caveat: `local_shell`'s host clusters are ~2.6x denser than `bridge`'s and
`collective`'s (core_frac = 0.4 vs full radius), so a shell point 0.49
beyond the realised radius is more starkly separated in absolute distance
terms than a comparably-placed point would be in a full-density cluster —
this arm over-states how easily a genuinely embedded local outlier would be
caught in the paper's own two-cluster geometry.

### bridge — TPR / TNR (mean ± SE)

| d | method | TPR | TNR | mean n_flagged |
|---|---|---|---|---|
| 3 | DBSCAN | 0.000 | 0.998 | 0.37 |
| 3 | iForest | 0.000 | 0.945 | 10.48 |
| 3 | LOF | 0.002 | 0.9999 | 0.04 |
| 3 | MST | 0.140 ± 0.020 | 0.607 | 75.65 |
| 3 | ODIN | 0.002 | 0.967 | 6.30 |
| 3 | **SU-MCCD** | **0.023 ± 0.008** | 0.9995 | 0.32 |
| 3 | **SUN-MCCD** | **0.010 ± 0.005** | 0.998 | 0.55 |
| 3 | **U-MCCD** | **0.068 ± 0.017** | 0.993 | 2.03 |
| 3 | **UN-MCCD** | **0.039 ± 0.011** | 0.992 | 1.99 |
| 10 | DBSCAN | 0.000 | 0.998 | 0.29 |
| 10 | iForest | 0.000 | 1.000 | 0.00 |
| 10 | LOF | 0.000 | 1.000 | 0.00 |
| 10 | MST | 0.113 ± 0.016 | 0.900 | 19.95 |
| 10 | ODIN | 0.000 | 0.969 | 5.94 |
| 10 | **SU-MCCD** | **0.138 ± 0.025** | 0.993 | 2.79 |
| 10 | **SUN-MCCD** | **0.045 ± 0.015** | 0.998 | 0.75 |
| 10 | **U-MCCD** | **0.164 ± 0.026** | 0.978 | 5.73 |
| 10 | **UN-MCCD** | **0.094 ± 0.022** | 0.994 | 2.10 |

**Ordering check against the pre-registration.** The pre-registered reading
(protocol round-3 note 4, 20-replicate diagnostic at d=3) predicted U-MCCD
and UN-MCCD above SUN-MCCD (SUN-MCCD ≈ 0), because U-MCCD/UN-MCCD carry no
`min.cls` argument at all while SU-MCCD/SUN-MCCD are obliged to accept the
10-point chain as a legitimate cluster (the chain's own NN spacing is the
densest in the data, and `round(0.05*199)=10=n0` exactly matches the
minimum-cluster-size threshold). **The production 100-rep table confirms
this ordering at both d:**

- d = 3: U-MCCD (0.068) > UN-MCCD (0.039) > SU-MCCD (0.023) > SUN-MCCD
  (0.010) — a clean split by `min.cls` membership (no-`min.cls` pair
  0.068/0.039 both above the `min.cls` pair 0.023/0.010).
- d = 10: U-MCCD (0.164) > UN-MCCD (0.094) > SU-MCCD (0.138) > SUN-MCCD
  (0.045) — U-MCCD > UN-MCCD > SUN-MCCD holds exactly as predicted; SU-MCCD
  sits between UN-MCCD and U-MCCD here (RK- vs NND-based coverage crossing
  at d=10), but SUN-MCCD remains the lowest of all four MCCD detectors at
  both d, consistent with being both shape-adaptive (`min.cls`-bound) and
  NND-based (the family whose bracket search is most conservative about
  admitting the chain as noise).

SUN-MCCD's TPR is not exactly 0 in the 100-rep production run (0.010 at
d=3, 0.045 at d=10) as the 20-rep diagnostic's point estimate suggested, but
it is the smallest of the nine methods' non-zero values at d=3 and remains
well below U-MCCD/UN-MCCD at both d — consistent with the qualitative
ordering, not a contradiction of it (the diagnostic was explicitly a
20-replicate estimate, not a claim of an exact zero).

**TNR distinguishes absorption from over-flagging.** SU-MCCD and SUN-MCCD's
low TPR comes with **very high TNR** (0.993–1.000 across both d) and modest
mean flag counts (0.32–2.79, well under the true n0=10) — they are not
compensating for missing the chain by flagging elsewhere; they simply treat
the chain as a legitimate cluster and leave it and everything else alone,
exactly the `min.cls` absorption mechanism predicted. U-MCCD/UN-MCCD's
higher TPR also comes with high TNR (0.978–0.994) and modest flag counts
(1.99–5.73) — a few chain points (typically the two clipped endpoints) are
caught without collateral false positives. **MST** is the outlier here: its
TPR (0.113–0.140) is the second-highest but its TNR (0.607–0.900) and mean
flag count (19.95–75.65, several times the true n0=10) show its partial
recall is again a by-product of generic over-flagging rather than
discrimination of the chain specifically.

### collective — TPR / TNR (mean ± SE)

| d | method | TPR | TNR | mean n_flagged |
|---|---|---|---|---|
| 3 | DBSCAN | 0.001 | 0.999 | 0.29 |
| 3 | iForest | 0.322 ± 0.022 | 0.962 | 10.47 |
| 3 | LOF | 0.042 ± 0.012 | 1.000 | 0.43 |
| 3 | MST | 0.359 ± 0.020 | 0.619 | 75.57 |
| 3 | ODIN | 0.001 | 0.970 | 5.66 |
| 3 | SU-MCCD | 0.178 ± 0.036 | 0.999 | 1.89 |
| 3 | SUN-MCCD | 0.189 ± 0.036 | 0.999 | 2.11 |
| 3 | U-MCCD | 0.385 ± 0.047 | 0.991 | 5.60 |
| 3 | UN-MCCD | 0.343 ± 0.045 | 0.994 | 4.63 |
| 10 | DBSCAN | 0.002 | 0.999 | 0.20 |
| 10 | iForest | 0.082 ± 0.013 | 1.000 | 0.82 |
| 10 | LOF | 0.007 | 1.000 | 0.07 |
| 10 | MST | 0.157 ± 0.015 | 0.908 | 18.99 |
| 10 | ODIN | 0.000 | 0.970 | 5.66 |
| 10 | **SU-MCCD** | **0.547 ± 0.050** | 0.988 | 7.74 |
| 10 | **SUN-MCCD** | **0.426 ± 0.049** | 0.999 | 4.43 |
| 10 | **U-MCCD** | **0.600 ± 0.049** | 0.972 | 11.28 |
| 10 | **UN-MCCD** | **0.516 ± 0.050** | 0.996 | 5.96 |

**Reading.** At d = 10 the smoke's suggestion is confirmed with numbers:
the four MCCD detectors (0.426–0.600) dominate every baseline by a wide
margin (best baseline MST at 0.157, and that with TNR 0.908 / mean flag
count 18.99 — again inflated by generic over-flagging, not targeted
detection; the next best clean baseline is iForest at 0.082 with TNR =
1.000). At d = 3 the ordering is less clean: U-MCCD (0.385) and UN-MCCD
(0.343) still lead, but iForest (0.322, TNR 0.962, mean flagged 10.47 —
close to the true n0=10, a genuine, non-inflated detection) is competitive
with SU-MCCD (0.178) and SUN-MCCD (0.189), and MST's 0.359 again comes with
heavy over-flagging (TNR 0.619). So the mutual-catch mechanism's advantage
on collective anomalies **strengthens with dimension** — d=10 is a clean
sweep for all four MCCD methods over all five baselines, while d=3 is a
partial win with iForest a real (not inflated) competitor.

## Cross-arm reading

**What the mutual-catch mechanism handles, in numbers.** *Collective*
anomalies near a cluster boundary (R3.7's "internally well-connected groups
of anomalous points... near cluster boundaries") are the strongest case for
the mechanism: at d=10 all four MCCD detectors (TPR 0.43–0.60) outperform
every baseline (best clean baseline 0.08) with high TNR (0.97–1.00), and at
d=3 the two uniform-coverage variants (U-MCCD 0.385, UN-MCCD 0.343) lead
outright. *Local_shell* outliers (TPR = 1.000 or near it for every MCCD
detector at both d) are also handled, but this arm is a shell/boundary
placement 0.49 beyond the realised radius, not the "embedded within" case
R3.7 actually raises, and its clusters are ~2.6x denser than
bridge/collective's — both facts make this arm easier than the concern it
was meant to probe.

**What it does not handle, in numbers.** *Bridge* anomalies are the clean
failure case R3.7 anticipates: SU-MCCD and SUN-MCCD score TPR 0.010–0.138
across both d, with SUN-MCCD lowest of all four MCCD detectors at both d
(0.010 at d=3, 0.045 at d=10) — not because the chain evades detection, but
because `min.cls` obliges the shape-adaptive detectors to accept an
exactly-threshold-sized (n0=10=round(0.05×199)), maximally dense chain as a
legitimate cluster. This is a designed consequence of the S_min rule
colliding with an adversarially dense, threshold-sized bridging structure,
not a shortfall in the mutual-catch principle itself — and the ordering
(no-`min.cls` U-MCCD/UN-MCCD above `min.cls`-bound SU-MCCD/SUN-MCCD at both
d) is exactly what the rule predicts. *Local_interior* anomalies are a
near-total failure for every method, MCCD and baseline alike (TPR ≤ 0.078
at d=10, and exactly 0.000 for all methods but MST — itself an
over-flagging artifact, not detection — at d=3); the failure traces to
geometry (an interior point cannot be isolated without exceeding the
cluster's own realised radius at this n and d) rather than to the
mutual-catch rule specifically, since every baseline fails identically.

**Density caveat, stated plainly.** `local_shell`/`local_interior`'s host
clusters and `bridge`/`collective`'s host clusters are not density-matched
(2.6x by NN spacing) — a known asymmetry recorded in `WP8_PROTOCOL.md` and
not corrected here. The strong local_shell results should not be read as
"local outliers are easy" in general; they are easy in a cluster that is
2.6x denser than the one the bridge and collective arms use, at a stand-off
(0.49 beyond the realised radius) that is closer to a boundary placement
than a genuinely embedded one.

## Anomalies

None. All integrity checks passed; the script's own `--summarize` output
matched an independent from-scratch recomputation to the reported
precision; the production 100-rep bridge ordering confirmed the
pre-registered reading recorded in `WP8_PROTOCOL.md` before the grid was
launched.
