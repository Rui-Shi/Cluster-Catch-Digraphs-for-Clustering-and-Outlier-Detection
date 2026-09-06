# WP5 preparation inventory — real data above d = 21

Read-only reconnaissance for WP5 (AE.4, R1.8, R3.8, R5.6). No experiments were run;
no manuscript, harness, or data files were modified. All facts below are backed by
file paths, sizes, and small R snippets loaded in the foreground
(`C:/Program Files/R/R-4.6.1/bin/Rscript.exe` via PowerShell — Rscript segfaults
under the Bash tool per HANDOFF_FROM_TR2.md).

Candidate sets (REVISION_PLAN_NEUCOM-D-26-15191.md, WP5 table): letter (n=1600,
d=32), mnist (n=1000 subsample, d=100), musk (n=1000 subsample, d=166), arrhythmia
(n=452, d=274). Cut-list item 2 (§ REVISION_PLAN cut list) proposes letter +
arrhythmia as the minimal pair.

---

## 1. Data availability

| Dataset | Present in CCD/? | Evidence |
|---|---|---|
| **letter** | **No.** | No file anywhere under `CCD/` matches `*letter*` except unrelated matches (`.venv` Tk demo images, an R comment "letter." in `43_wp2a_direction_sweep.R:63`). Not in `revision_experiments/results/datasets_csv/`, not in `data/outlier_detection/`. |
| **mnist** | **No.** | Same search, zero hits (only an unrelated PyTorch header path inside `.venv`). |
| **musk** | **Yes**, full data + n=1000 subsample. | `data/outlier_detection/musk_raw.csv` (2,045,145 B), `musk.npz` (736,255 B); `revision_experiments/results/datasets_csv/Musk.csv` (3063 rows incl. header, 167 cols incl. label → n=3062, d=166) and `Musk_sub1000.csv` (1001 rows → n=1000, d=166). |
| **arrhythmia** | **Yes.** | `data/outlier_detection/arrhythmia_raw.csv` (362,490 B), `arrhythmia.mat` (105,540 B); `revision_experiments/results/datasets_csv/Arrhythmia.csv` (453 rows incl. header, 275 cols incl. label → n=452, d=274). |

### Present sets — n, d, outliers, preprocessing, label convention

Both loaded via `revision_experiments/tr2/02_load_data.R` (not
`RealData_Collection.R` / `load_real_dataset()` — see §4 for why that matters),
manifest at `revision_experiments/results/datasets_csv/manifest.csv`:

| Dataset | n | d | n_outliers | Source | Normalization |
|---|---|---|---|---|---|
| Musk | 3062 | 166 | 97 (3.2%) | `https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/25_musk.npz` (ADBench mirror of ODDS musk.mat, commit `3dac822`, fetched 2026-07-10) | robust median/MADN per column, SD/constant fallback on 3/166 zero-MADN columns |
| Arrhythmia | 452 | 274 | 66 (14.6%) | `https://raw.githubusercontent.com/BELLoney/Outlier-detection/master/ODDS/original/arrhythmiaori.mat` (community ODDS mirror, commit `a90424a`, fetched 2026-07-10; odds.cs.stonybrook.edu TLS-unreachable at acquisition time; Arrhythmia is absent from ADBench's Classical set) | robust median/MADN per column, SD/constant fallback on 149/274 zero-MADN columns |

Label convention verified directly, not just by manifest claim: CSV tail column is
named `label`; tally of Arrhythmia.csv gives 66 rows of `0` / 386 of `1`, Musk.csv
gives 97 of `0` / 2965 of `1`, Musk_sub1000.csv gives 32 of `0` / 968 of `1` — all
match the manifest's outlier counts exactly. `02_load_data.R:14-16,102-106,200-203`
documents the convention explicitly: the raw ADBench/ODDS `.npz`/`.mat` source uses
**1 = outlier**; the loader flips polarity so that, in every CSV under
`datasets_csv/`, **1 = regular, 0 = outlier**, matching this paper's harness
contract (`Y`: 1 = regular, 0 = outlier) and matching `RealData_Collection.R`'s own
convention ("1 represent regular observations, 0 are outliers", line 1) — no
adapter needed. Musk_sub1000's 32 outliers out of 1000 matches
REVISION_PLAN's candidate table exactly (32).

### Absent sets — canonical source, for download approval

Both are in ADBench's "Classical" collection on GitHub, following the identical
naming/fetch convention already used for Musk/Speech/InternetAds (confirmed via
GitHub API listing of `Minqi824/ADBench/adbench/datasets/Classical/`):

- **letter**: `https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/20_letter.npz` — ODDS "Letter Recognition" outlier set, canonically n=1600, d=32, 100 outliers (6.25%), matching REVISION_PLAN's table exactly.
- **mnist**: `https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/24_mnist.npz` — ODDS MNIST outlier set (full set is larger; REVISION_PLAN specifies a 1000-row contamination-preserving subsample, d=100, ~70 outliers).

Format: `.npz` (NumPy), same as Musk/Speech/InternetAds — `02_load_data.R`'s
existing `.npz` reader and label-flip/dedup/normalization pipeline would apply
unchanged; no new loader code needed, only new entries in that script's dataset
table and manifest, mirroring the Musk/Arrhythmia pattern. **Nothing was
downloaded** — this is flagged for the user's approval per the read-only
constraint.

---

## 2. Quantile-table inventory

Checked with `get_simul`-equivalent loads (`nrow(simul$quan[[1]])` for RK,
`length(simul$average)` / `length(simul$median)` for NN). Alpha-level tokens per
`revision_experiments/shared/harness.R`'s resolvers: `rk_quant_label_paper(d)` →
"999" for all four d (all ≥10); `nn_quant_label_paper_UN(d)` → "999" for all four
(all ≥20); `nn_quant_label_paper_SUN(d)` → "999" for all four (all ≥10). So every
cell in this WP needs the "999" (α = 0.1%) file, both families.

| d (dataset) | RK-test-simul_\<d\>d_999%.RData | Extent (nrow) | NN-test-simul_\<d\>d_999%.RData | Extent (avg/median len) |
|---|---|---|---|---|
| 32 (letter) | **EXISTS** `R/RK-test_quantile/` (264,952,612 B, built 2024-05-16) | 5000 | **EXISTS** `R/NN-test_quantile/` (70,300 B, built 2024-09-18) | 5000 / 5000 |
| 100 (mnist) | **EXISTS** (2,694,823 B, built 2023-10-01) | 1000 | **EXISTS** (14,077 B, built 2023-10-02) | 1000 / 1000 |
| 166 (musk) | **MISSING** from `R/RK-test_quantile/` (never promoted to production — deliberate scope decision, FINDINGS.md §3c) | — | **EXISTS** (14,033 B, built 2026-08-09) | 1000 / 1000; **plus** a spliced full-n table `NN-test-simul_166d_999%_n3062_spliced.RData` (42,053 B) covering **3062** |
| 274 (arrhythmia) | **MISSING** from production dir | — | **EXISTS** (14,031 B, built 2026-08-09) | 1000 / 1000 |

Table-vs-n coverage:

| Dataset | Planned n | RK covers it? | NN covers it? |
|---|---|---|---|
| letter | 1600 | Yes (table extent 5000 ≥ 1600) | Yes (5000 ≥ 1600) |
| mnist | 1000 | Yes (extent 1000 = n exactly, boundary-valid) | Yes (1000 = n exactly) |
| musk | 1000 subsample or 3062 full | **n/a — no RK table at d=166 at all** | Yes, both n=1000 (regular table) and n=3062 (the spliced table TR2 already built for its own full-data run) |
| arrhythmia | 452 | **n/a — no RK table at d=274** | Yes (1000 ≥ 452) |

**Correction to the revision plan's own assumption:** the plan's WP5 table marks
letter's table as "**generate** (~2 h, 20 cores)". That is stale — a production
RK table at d=32/999% and a production NN table at d=32/999% both already exist,
both with 5000-row extent, comfortably covering n=1600. **No new table generation
is needed for letter under either method family.** (These d=32/d=100 tables
predate this revision cycle — 2023–2024 timestamps — and are part of the original
51/58-file production inventory, not something TR2 built for WP5.)

For RK at d=166/274: not a gap to fill by generation. FINDINGS.md §3/§3c records
a **deliberate decision** not to build these, because the tables are numerically
generatable but statistically vacuous (see §3 below) — building them would not
change the "n/a" outcome, only waste 3–4 core-hours each (probe-measured:
~3.1 h and ~3.6 h at niter=10000, `probe_report.csv`).

**Cost estimate if letter/mnist RK or NN tables ever needed generating** (they
don't, per above), using the measured cadence model from
`revision_experiments/71_regen_nnd_alpha001.R` header (anchor: 396 s/iter at
d=16, n=3200; NN cost ≈ n³·d/3, so `sec_per_iter ≈ 396 × (n/3200)³ × (d/16)`,
niter=10000 for a 0.1% table): at n=1600, d=32 this gives ≈ 396 × (1600/3200)³ ×
(32/16) = 99 s/iter × 10000 iters ≈ 275 h serial, ÷ 20 cores ≈ **13.8 h**. That
would have been the generation cost had the table not already existed — moot here,
but relevant if any future WP5 extension needs a table at an uncovered n.

---

## 3. RK degeneracy — measured directly, not just cited

Computed `mean(simul$quan[[1]] == 0)` on the production/probe tables:

| d | Fraction zero (measured here) | Source | HANDOFF's figure |
|---|---|---|---|
| 32 (letter) | **51.276%** (25,638 / 50,000 cells, 5000×10 matrix) | production table | (not listed; consistent with CLAUDE.md's "51% zeros by d=35") |
| 100 (mnist) | **91.4%** (9,140 / 10,000 cells, 1000×10 matrix) | production table | "91.4% at d=100" — exact match |
| 166 (musk) | **100.0%** (probe table, 1000×10, niter=20) | `revision_experiments/results/tr2/probes/RK-test-simul_166d_999%.RData` | "100% at d=166" — exact match |
| 274 (arrhythmia) | **100.0%** (probe table, 1000×10, niter=20) | `.../probes/RK-test-simul_274d_999%.RData` | "100% at d=274" — exact match |

**U-MCCD and SU-MCCD (the RK-based methods) would be reported `n/a`** at **musk
(d=166) and arrhythmia (d=274)** — no RK table exists there at all, consistent
with the 100% structural degeneracy, and REVISION_PLAN already states this rule
("U-MCCD and SU-MCCD rows at d ≥ 166 are n/a"). At **letter (d=32) and mnist
(d=100)** the RK tables exist and are technically usable, but 51%/91% of the CSR
envelope is zero — over half (letter) to nearly all (mnist) of the radius
searches would short-circuit at the first candidate. REVISION_PLAN's own n/a
threshold is d≥166, so it does not call these two n/a, but any reported
U-MCCD/SU-MCCD row at d=100 (mnist) should carry the same degeneracy caveat as
prose, since 91.4% zero is barely different in kind from the 100% case just one
tier up.

---

## 4. Reuse of TR2's existing WP5 work — mostly not directly applicable

TR2 already ran high-d real-data experiments (`06_wp5_highdim.R`,
`07_wp5_subsample_ccd.R`, `07b_wp5_fulldata_ccd.R`) on Arrhythmia, Musk (both
n=1000 subsample and full n=3062), and Speech — but on a **different method
family**: `UNCCD-OOS` / `UNCCD-IOS`, TR2's own simple outlyingness-score
thresholding built directly on the UN-CCD digraph (`harness.R`'s
`METHOD_REGISTRY`, `unccd_oos_method` / `unccd_ios_method`, threshold=2 from
`REAL_DATA_THRESHOLDS`). This paper's four detectors (U-MCCD, SU-MCCD, UN-MCCD,
SUN-MCCD) are wired into a **separate** registry
(`revision_experiments/tr1/wp0_mccd_methods.R`, entry points `RUMCCD_outlier`,
`SUMCCD_outlier`, `UNMCCD_outlier`, `SUNMCCD_outlier`) that adds the
mutual-catch-graph layer and cluster-connectivity outlier rule on top of the same
digraph construction — a materially different outlier decision, not a relabeling
of the same score. **No TR2 WP5 cell can be read as a result for this paper's
methods; UN-CCD (the shared construction) is a component, but the outlier call
differs.**

What IS reusable, narrowly:

- **Data**: Musk.csv, Musk_sub1000.csv, Arrhythmia.csv and their manifest rows
  (§1) — directly usable by this paper's loaders without re-preprocessing.
- **Quantile tables**: the NN d=166/274 production tables and the d=166 spliced
  n=3062 table (§2) — built by TR2 for its own UNCCD-OOS/IOS runs, but the table
  content (Monte Carlo NND quantiles) is method-agnostic; this paper's
  `UNMCCD_outlier`/`SUNMCCD_outlier` consume the identical file format via the
  identical `get_simul("NN", d, quant=...)` call. These ARE directly reusable.
- **S_min / min.cls**: TR2's UNCCD-OOS/IOS calls have no `S_min` concept at all
  (a simple score threshold, not a cluster-coverage rule) — nothing to compare
  or mismatch against this paper's S_min=0.05.
- **Alpha**: TR2 used `nn_quant_label_paper_UN(d)` → "999" for d=166/274/400,
  the same resolver and token this paper's UN-MCCD would use at those d. SUN-MCCD
  also resolves to "999" at these d, so no mismatch there either. (RK: not
  applicable, TR2 never built or used an RK table at d≥166.)
- **Runtime/failure-mode findings** (FINDINGS.md, not reusable as numbers but
  as prior art for what to expect): Musk's UNCCD-OOS/IOS AUC is near-degenerate
  (score IQR exactly 0 before the std_MADN fix, 0.7357 after); MST baseline
  degenerates on Musk (TPR=1, TNR=0); the full-Musk radius search cost 53,006 s
  (14.7 h) at 20 cores, full-Speech cost 112,968 s (31.4 h) — both far above the
  a-priori "~2 h" forecast recorded in REVISION_PLAN's own WP5 table, which was
  written before TR2's measurement. These are cautionary numbers for scoping,
  not results to cite for this paper's four methods.

**Bottom line: nothing computationally reusable for the paper's own four MCCD
detectors' scores** — every (dataset, method) cell for U-MCCD/SU-MCCD/UN-MCCD/
SUN-MCCD on musk/arrhythmia/letter/mnist still needs to be run fresh. What is
reused is infrastructure: data files, quantile tables, and the alpha-token
resolution — all confirmed compatible above.

---

## 5. Runtime — measured on this paper's own four methods, not TR2's

HANDOFF_FROM_TR2.md's WP4 runtime numbers are for TR2's UN-CCD/baselines, not
this paper's methods, and REVISION_PLAN's WP5 section budgets "≈540 s per method
at n=1000 regardless of d" — but that figure appears to originate from the same
TR2-era ~n⁴ UN-CCD construction cost, which WP0's outcome explicitly states does
**not** apply to this paper's wired detectors ("Runtime is ~n², not the ~n⁴
inherited from TR2"). The WP0 gate run (`revision_experiments/results/tr1/
wp0_gate*.csv`) already measured this paper's four detectors directly on real
data, at n up to 4819:

| Dataset | n | d | U-MCCD | SU-MCCD | UN-MCCD | SUN-MCCD |
|---|---|---|---|---|---|---|
| stamps | 340 | 9 | 2.4 s | 0.4 s | 0.5 s | 0.3 s |
| vowels | 1452 | 12 | 26.1 s | 18.0 s | 24.7 s | 16.6 s |
| waveform | 3443 | 21 | 381.3 s | 266–278 s | 424.5 s | 255–270 s |
| wilt | 4819 | 5 | 1130.6 s | 948–1015 s | 1051.1 s | 948–1015 s |

(`t_construct` and `t_total` are identical by design — see `wp0_mccd_methods.R`'s
comment: all four detectors are single monolithic calls with no exposed
sub-timing hook.)

These points show runtime tracking **n**, not d (waveform at d=21 and wilt at
d=5 cost the same order of magnitude; stamps at d=9 and vowels at d=12 likewise
track n, not d). **Extrapolating from vowels (n=1452, ~17–26 s/method) down to
n=1000 and up to n=1600, expect roughly 10–30 s per method per data set** — one
to two orders of magnitude below the plan's 540 s/method budget. This is good
news for scope: the plan's WP5 compute budget is very likely a large
overestimate once one substitutes this paper's own measured cost model for
TR2's inherited one.

At n=3062 (musk full) or n=4819 (wilt, the largest measured point), cost is in
the 950–1130 s range per method; a full-data musk run (n=3062, between waveform's
3443 and vowels' 1452) would land in the same few-hundred-to-thousand-second
range per method, cheap enough to run in full rather than subsampling — worth
reconsidering whether the plan's "1000-row subsample" for musk is even necessary
for this paper's methods (it was necessary for TR2's ~n⁴ UN-CCD construction,
which does not apply here).

---

## 6. Recommended scope

**Minimal scope (cut-list item 2): letter + arrhythmia only.**

- Data: letter needs downloading (`20_letter.npz`, §1) — the only new-data step.
  Arrhythmia is already present.
- Tables: **zero generation needed.** letter's RK/NN d=32/999% tables already
  exist at sufficient extent (5000 ≥ 1600); arrhythmia's NN d=274/999% table
  already exists (1000 ≥ 452); arrhythmia's RK table doesn't exist and won't be
  built (100% degenerate, §3), so U-MCCD/SU-MCCD report `n/a` there by design,
  not by omission.
- Compute: ~4 detector calls × 2 datasets, each ≲30 s at these n (§5) — call it
  under 5 minutes of CCD compute, plus baselines (LOF/DBSCAN/MST/ODIN/iForest,
  cheap) and quantile-table loads.
- **Answers**: letter sits at d=32, inside R1.8's named 20<d<50 band, with all
  four methods reporting (U-MCCD/SU-MCCD caveated at 51% RK degeneracy, not
  n/a). Arrhythmia at d=274 gives the extreme high-d anchor with UN-MCCD/
  SUN-MCCD only (U-MCCD/SU-MCCD n/a, itself a documented finding answering
  R1.2's "why do RK-CCD balls collapse in high-d" question directly).
- **Does not answer**: the d=100 midpoint (nothing in 50–166 gets tested), so
  the trajectory of degradation between the d=32 partial-degeneracy point and
  the d=166+ total-degeneracy point is not shown for the RK family; only two
  points anchor R1.8's band rather than three or four.

**Alternative needing no new table generation: arrhythmia + musk.**

- Data: both already present in full (§1) — literally zero download or table
  work of any kind.
- Tables: arrhythmia NN d=274 exists (1000≥452); musk NN d=166 exists at both
  n=1000 (subsample) and n=3062 (full, via the spliced table) — full data is
  usable without subsampling given this paper's actual (~n²-ish, not TR2's ~n⁴)
  cost profile (§5).
- Compute: at n=452/1000/3062, expect on the order of tens of seconds
  (arrhythmia, musk-sub1000) up to roughly 15–20 minutes (musk full, by
  interpolation between waveform's 3443/381s and wilt's 4819/1130s) per
  UN-MCCD/SUN-MCCD call; U-MCCD/SU-MCCD are `n/a` at both (both d≥166).
- **Answers**: two points deep in the high-d degenerate regime (d=166, d=274),
  both with a clean UN-MCCD/SUN-MCCD number and a documented RK n/a, and — if
  full-data musk is run instead of the subsample — one genuinely large-n
  high-d cell (n=3062) that stress-tests the calibrated-cutoff transfer TR2
  already found breaks down for its own scores (FINDINGS.md §9).
- **Does not answer**: R1.8's own named band (20<d<50) is untouched entirely —
  both sets sit at d≥166, well outside "20 < d < 50" — so this alternative
  answers "does the method still run/degrade further at very high d" but not
  the reviewer's literal request for evidence in the gap between the current
  d=21 ceiling and the claimed d≥50 degradation onset. It also gives no new
  RK-family evidence beyond confirming n/a again.

**Assessment**: the minimal scope is worth strictly more per unit of new work —
it requires one download (letter) against zero new table generation either way,
and it is the only option of the two that puts a data point inside R1.8's
literally-quoted band. The "no new tables" framing turns out to be true of
*both* options, since none of the four candidate d's need table generation
regardless of which pair is chosen (§2's correction). The real cost/benefit
choice is therefore just "one dataset download" (letter) vs. "an extra
15–20 minute full-data musk run" (musk) — both cheap; letter buys the band
coverage the AE explicitly asked for and musk does not.
