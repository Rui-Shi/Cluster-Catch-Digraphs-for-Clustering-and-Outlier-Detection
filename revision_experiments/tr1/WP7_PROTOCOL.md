# WP7 protocol — clustering quality (R3.6)

**Declared 2026-09-05, before `90_wp7_clustering_quality.R` was run on anything
but the WP8 smoke outputs used to develop it.** WP7 is scoring only: it runs
no detector and generates no new score. It reads the cluster-assignment files
WP8's four scripts already write and scores them. Per the cut list
(`REVISION_PLAN_NEUCOM-D-26-15191.md`, item 3), WP7 is **synthetic settings
only, no contamination sweep** — the four WP8 experiments (85 boundary FP,
86 outlier types, 87 small cluster, 88 CSR violation) are the entire input.

## 1. Inputs

Read-only. WP7 never writes into `results/tr1/wp8/`.

| script | production file | smoke file | setting identity |
|---|---|---|---|
| 85 | `85_boundary_fp_clusters.csv` | `smoke/85_wp8_boundary_fp_clusters.csv` | `setting_id` (`uniform_d3` etc.), generator recovered by stripping `_d<d>$` |
| 86 | `86_outlier_types_clusters.csv` | `smoke/86_wp8_outlier_types_clusters.csv` | `setting_id`, type recovered the same way |
| 87 | `87_small_cluster_clusters.csv` | `smoke/87_wp8_small_cluster_clusters.csv` | `(m, d)`, `setting_key <- sprintf("m%d_d%d", m, d)` |
| 88 | `88_csr_violation_clusters.csv` | `smoke/88_wp8_csr_violation_clusters.csv` | `setting_id`, generator recovered the same way as 85/86 |

All four share the schema `(..., d, rep, seed, method, row_index, true_cluster,
detected_cluster)` per `WP8_PROTOCOL.md`'s outputs summary and the task
brief: `detected_cluster = NA` means unassigned by the MCCD method;
`true_cluster` is `{1,2}` for 85/86, `{1,2,3}` for 87, always `1` for 88 (a
single homogeneous population — see §3 below on why 88's ARI/NMI/AMI are
individually uninformative there).

`--smoke` points at `results/tr1/wp8/smoke/`; as of this writing that
directory holds only 86's and 87's cluster files (85 and 88 passed their own
smoke gate during WP8 development and were launched straight to production,
per `WP8_REVERIFICATION.md`), so a `--smoke` run processes 86 and 87 and logs
"no cluster file found" for 85 and 88 rather than failing.

## 2. De-duplication and completeness (per WP8_REVERIFICATION.md's
"Analysis-time items for all four")

- **De-dup key**, keeping the **last** occurrence: `(setting/m-key, d, rep,
  method, row_index)`. This is the natural key WP8_REVERIFICATION.md names —
  a retried cell can leave duplicate rows on the USB-enclosure-interrupted
  writes this repo is exposed to.
- **Completeness assertion**, per `(setting-key, d, rep, method)` cell: the
  number of distinct `row_index` values must equal the script's fixed regular
  -point count — **189 for 85/86** (`n1 = round(200*0.95*0.5) = 95`,
  `n2 = 95 - 1 = 94`, both constant across every setting and `d` in those two
  scripts because they depend only on the fixed `n = 200, cont = 0.05`) and
  **200 for 87/88** (no contamination outliers in either generator, so every
  point is "regular"). A cell that fails this check is **not silently
  dropped**: it is written to the output with `status = "incomplete"` and
  sentinel metric values (see §5), so its absence is visible in the summary
  rather than invisible.
- **A done cell may lack cluster rows after a crash** (WP8_PROTOCOL.md notes
  this for 86/88 specifically, where the cluster block is written *after*
  the metrics row). WP7 treats a `(setting, d, rep, method)` combination that
  never appears in the cluster file at all as simply absent from WP7's
  input — it contributes to neither a completeness failure nor a computed
  row, and the smoke/full run logs the count of methods actually seen per
  script for visibility.

## 3. Metrics

**ARI**: `mclust::adjustedRandIndex` (installed, confirmed —
`mclust 6.1.3`). Checked against hand-built cases (identical partitions → 1;
one degenerate all-one-cluster partition against a 2-cluster one → 0; both
degenerate → 1; a random pairing → a small negative number) before use.

**NMI and AMI: `aricode` is NOT installed** in this R installation (checked
2026-09-05, `requireNamespace("aricode", quietly=TRUE)` returns `FALSE`), so
both are implemented directly from the contingency table in
`90_wp7_clustering_quality.R`:

- Contingency table `cm[i,j]` = count of points with true label `i` and
  predicted label `j`; row/column sums `a`, `b`; `N = sum(cm)`.
- **A real bug found and fixed while writing this**, worth recording because
  it produces silently wrong (not erroring) output: `table()` on factors
  carries dimnames, and single-bracket indexing (`cm[i,j]`, `a[i]`, `b[j]`)
  propagates those dimnames onto the resulting scalar as a `names`
  attribute. `c(nmi = nmi, ami = ami)` then does not produce a vector named
  `c("nmi","ami")` — it produces `c("nmi.<level>", "ami.<level>")`, because
  `c()` concatenates the outer name with the already-present inner name.
  `r["nmi"]` on that vector returns `NA`, silently, with no warning. Fixed by
  `cm <- unname(as.matrix(cm))` immediately after building the table, before
  any indexing.
- **Mutual information**: `MI = sum_{i,j} (n_ij/N) log(N n_ij / (a_i b_j))`
  over `n_ij > 0`.
- **Entropy**: `H(x) = -sum (x_k/N) log(x_k/N)` over `x_k > 0`.
- **NMI** (arithmetic-mean normalization, matching sklearn's default
  `average_method="arithmetic"`): `NMI = 2*MI/(H(U)+H(V))`, with the
  standard degenerate-case convention — both partitions single-cluster
  (`H(U)=H(V)=0`) → `NMI = 1`; exactly one degenerate → `NMI = 0`.
- **AMI**: the exact hypergeometric expected-MI correction of Vinh, Epps &
  Bailey (2010) — `E[MI] = sum_{i,j} sum_{n_ij=max(1,a_i+b_j-N)}^{min(a_i,b_j)}
  (n_ij/N) log(N n_ij/(a_i b_j)) * P(n_ij)`, where `P(n_ij)` is the
  hypergeometric probability of that overlap given the margins, computed in
  log-space via `lfactorial()` for numerical stability. `AMI = (MI - E[MI]) /
  (mean(H(U),H(V)) - E[MI])`, again matching sklearn's arithmetic-mean
  default, with the same both-degenerate → 1 convention.
- **Validated against `sklearn.metrics` (scikit-learn 1.9.0, the same
  `revision_experiments/.venv`)** on 7 hand-built cases — identical labels,
  a permutation of the same labels, one degenerate partition, both
  degenerate, two random-label pairs, and an all-singletons partition — with
  exact agreement to 6 decimal places on NMI, AMI, and ARI in every case
  (see the task report for the full table). This is the closest available
  substitute for a published-package cross-check given `aricode` is absent.

## 4. Unassigned handling — **both treatments reported, per method**

A point with `detected_cluster = NA` (unassigned by `mccd_translate()`, or
DBSCAN/HDBSCAN noise, see §6) is handled two ways, both cheap from the same
contingency table machinery, so both are computed and written:

- **`excluded`** (primary/headline): unassigned points are dropped before
  scoring. ARI/NMI/AMI then describe the quality of the partition **among
  the points the method actually committed to a cluster** — the convention
  used when comparing HDBSCAN-style noise-labelling methods against ground
  truth in the density-clustering literature. Chosen as primary because
  R3.6's ask ("quality of the recovered cluster structure") is about the
  clusters themselves, and a method that abstains on hard points but
  clusters the rest well is a different failure mode from one that clusters
  everyone poorly.
- **`singleton`** (robustness check): every unassigned point gets its own
  unique singleton label (`paste0("UNASSN_", row_index)`), so it is never
  treated as agreeing with anything. This is the conservative reading — a
  method cannot improve its score by simply declining to commit — and is
  reported alongside `excluded` rather than replacing it.

`unassigned_frac = n_unassigned / n_reg` is reported once per cell
(identical under both treatments, since it is defined before either is
applied) and is the number that should carry the headline "how often does
this method abstain" claim; the two ARI/NMI/AMI treatments are the "and what
does that mean for what it does commit to" complement.

**Degenerate cases** (documented, not treated as errors): if `excluded`
leaves fewer than 2 points, or fewer than 2 points and only one now-visible
true or predicted label, metrics are set to `0` with
`status = "degenerate"` and a `note` explaining why — never `NA`, per the
NA-avoidance convention below. 88's single true cluster (`true_k = 1`
throughout) makes ARI/NMI/AMI individually uninformative there by
construction, exactly as `WP8_PROTOCOL.md` already anticipated ("the WP7 use
is k-hat accuracy... more than ARI/NMI, which are degenerate with one true
cluster") — WP7 still computes them for completeness (mclust's ARI handles
the constant-partition case cleanly, returning 0 or 1 rather than `NaN`,
confirmed empirically) but the k-hat distribution is what carries the actual
88 finding.

## 5. k-hat convention

`k_hat = length(unique(detected_cluster[!is.na(detected_cluster)]))` over
the cell's regular points — the count of distinct **assigned** macro-cluster
labels, independent of the unassigned treatment (the unassigned bucket is
not counted as a cluster in either treatment; `singleton` fragments it into
many size-1 clusters for the ARI/NMI/AMI computation only, not for k-hat).
`true_k` is 2 (85/86), 3 (87), 1 (88).

Reported as **a count table**, per the task brief ("as counts per value, not
a mean"): for each `(script, setting-key, d, method)`, how many replicates
landed at each observed `k_hat` value, alongside `true_k` for reference.

**Cross-check**: 87's own driver already computes an equivalent quantity —
`87_small_cluster.csv`'s `n_clusters` column is `length(unique(res$cluster[
!is.na(res$cluster)]))`, computed at run time in `87_wp8_small_cluster.R`
(line ~152) — the identical definition WP7 uses independently from the
cluster file. `90_wp7_clustering_quality.R`'s smoke run cross-checks WP7's
own `k_hat` against that column for every smoke cell where the main metrics
file is available, and reports any mismatch as an error (would indicate a
translation-layer bug, not a WP7 bug, since both compute the same quantity
from the same `res$cluster` object at different removes).

## 6. DBSCAN / HDBSCAN comparison design

R3.6 asks for a direct comparison against "algorithms that jointly identify
clusters and noise points." WP4 already runs DBSCAN and HDBSCAN on the real
data; WP7 needs them on the **same synthetic cells** the four MCCD methods
were scored on, which means regenerating each cell's data from its recorded
seed rather than re-deriving it independently.

**Regeneration.** Each of `85_wp8_boundary_fp.R`, `86_wp8_outlier_types.R`,
`87_wp8_small_cluster.R`, `88_wp8_csr_violation.R` guards its own
CLI-dispatch tail with `if (!isTRUE(getOption("wp8.no_main", FALSE)))`
(`WP8_PROTOCOL.md`'s "Schema and chunking interface" note), specifically so
its generator functions can be `source()`d for read-only reuse without
triggering a live run. `90_wp7_clustering_quality.R` sets
`options(wp8.no_main = TRUE)`, `sys.source()`s each of the four scripts into
its **own** environment (all four define `GENERATORS`, `CANON`, `BASE_SEED`,
etc. under the same names, so each must live in its own environment or they
clobber each other), and restores the option afterward. The seed used for
regeneration is **read directly from the cluster file's own `seed` column**
— WP8's seed formula (`BASE_SEED + 100000*s_i + rep`, `s_i` = the setting's
row index in the script's fixed `CANON` table) is never recomputed, so
regeneration cannot silently drift if `CANON`'s row order ever changed.
Calling the WP8 script's own `GENERATORS[[key]](seed, n, d, ...)` (or, for
87, the one `gen_small_cluster(seed, n_total, d, m)`) reuses the exact code
path WP8 ran, not a reimplementation — the only way to guarantee bit-for-bit
reproduction rather than an approximation of it.

**Verification (not just "a few cells" — every cell WP7 touches)**: after
regenerating, WP7 compares the regenerated `true_cluster` vector and row
count against the cluster file's own `true_cluster`/row count for that cell
(using whichever MCCD method's rows are present as the reference, with an
additional cross-method consistency check when more than one method's rows
exist for the same cell). A mismatch is logged as `status = "error"` for
that cell's DBSCAN/HDBSCAN rows rather than silently trusted or fatal to the
whole run.

**DBSCAN**: `harness.R`'s own `DBSCAN(X, k=4, quant=cont)` (from
`simulations/outlier_detection/Algo_Compare_OutlierDetection/DBSCAN/DBSCAN.R`,
already sourced transitively by every WP8 script) — same `k=4` (MinPts) as
every other use of this wrapper in the revision. Its native cluster labels
(`dbscan::dbscan(...)$cluster`, `0` = noise) are used directly, not the
binarized outlier score the WP8 driver computes from them.
`cont` = oracle contamination, **exactly as the WP8 driver would compute it
if it called DBSCAN on that cell**:

| script | cont |
|---|---|
| 85, 86 | `n0 / (n0 + n_reg)` from the regenerated cell (matches `dbscan_method()`'s own `sum(Y==0)/length(Y)`, since the WP8 driver's `Y` vector is exactly `c(rep(1,n_reg), rep(0,n0))`) |
| 87 | `0` — 87 never calls DBSCAN or any baseline (only the four MCCD methods, per `WP8_PROTOCOL.md`'s method-list restriction for experiment 3); there is no precedent to match, and cluster 3 is a **legitimate cluster**, not an outlier population, so a nonzero oracle contamination would misrepresent what DBSCAN is being asked to find. `0` is the closest analogue to what the MCCD methods themselves receive (`Y = NULL` — no contamination hint at all, only S_min) |
| 88 | `0`, by design — no outliers exist in any of the four generators, and `WP8_PROTOCOL.md` already documents DBSCAN's `cont=0` behaviour there ("flag rate 0 by construction") for the binarized-score use; the **native cluster labels** are still informative here even though the binarized outlier score is not, since `cont=0` sets `eps` to the widest (most merging) quantile of the 4th-NN distances, and whether that still fragments the population is itself a finding |

**Limitation, stated plainly**: `DBSCAN()`'s `eps` rule was designed for the
real-data outlier-detection use, where an oracle contamination level always
exists and is meaningful. 87 and 88 have no true outliers, so `cont=0`
there is not "no information" in the way it is for 88's honest CSR-null
design — for 87 specifically it means DBSCAN's `eps` is derived from a
quantity (contamination) that has no referent in that generator at all. This
is reused rather than redesigned because building a DBSCAN variant with a
different eps rule for 87 would be a new baseline, not the paper's own
DBSCAN reused for a fair comparison — the mismatch is reported as a named
limitation of the comparison, not silently absorbed.

**HDBSCAN**: **python**, not R (`hdbscan.HDBSCAN` 0.8.44, the same install
WP4 already verified in `revision_experiments/.venv`,
`revision_experiments/.venv/python.exe`) — chosen for consistency with WP4's
own HDBSCAN comparison on the real data, and because `WP4_PROTOCOL.md`
already records that sklearn 1.9's own `HDBSCAN` exposes no GLOSH score
(irrelevant to WP7, which only needs native labels, but the python package
was the one already vetted and installed for this exact purpose). Same
defaults as WP4: `min_cluster_size=5, min_samples=None`. Called via a small
helper, `90b_wp7_hdbscan.py`, over a per-cell CSV round trip (R writes the
regenerated `X` to a temp CSV, calls the helper via `system2()`, reads the
label CSV back, deletes both temp files) — resumable in the same sense as
the rest of WP7 (a crash mid-cell leaves no half-written row in the output,
since the round trip completes or errors before any `append_result()`
call). `labels_ == -1` → unassigned; the helper's raw non-negative labels are
used as-is (ARI/NMI/AMI/k-hat are invariant to a relabelling, so no
1-based/0-based reconciliation is needed).

**Limits of the DBSCAN/HDBSCAN comparison overall**: neither method receives
`S_min` or any hint about the intended cluster count; DBSCAN's `k=4`
(MinPts) and HDBSCAN's `min_cluster_size=5` are both defaults inherited from
elsewhere in this revision, not tuned for `n=200` two/three-cluster
synthetic geometries — the same criticism the manuscript itself makes of
naïve baselines, so it is stated here for fairness. A single HDBSCAN fit per
cell (deterministic, no seed sensitivity beyond the data's own seed) and a
single DBSCAN fit per cell (also deterministic given `cont`).

## 7. Outputs

```
results/tr1/wp7/wp7_metrics_long.csv
  script, setting_key, d, rep, seed, method, unassigned_treatment,
  n_reg, n_unassigned, unassigned_frac, true_k, k_hat,
  ari, nmi, ami, status, note
  -- one block of 2 rows (unassigned_treatment in {excluded, singleton})
  per (script, setting_key, d, rep, method); method is one of the 4 MCCD
  methods present in that script's cluster file, or DBSCAN / HDBSCAN.
  Resumable: has_result() gated on
  (script, setting_key, d, rep, method, unassigned_treatment = "excluded")
  as the marker that a cell's 2-row block was fully written (the block is a
  single append_result() call, so a torn write leaves a ragged trailing line
  that has_result()'s own count.fields() check already treats as absent).

results/tr1/wp7/wp7_summary_metrics.csv   (--summarize)
  script, setting_key, d, method, unassigned_treatment,
  n_reps, ari_mean, ari_se, nmi_mean, nmi_se, ami_mean, ami_se,
  unassigned_frac_mean

results/tr1/wp7/wp7_khat_table.csv        (--summarize)
  script, setting_key, d, method, true_k, k_hat, n_reps_at_this_khat

results/tr1/wp7/smoke/...                 --smoke output, same schemas,
                                           never the production files
```

**NA-avoidance convention** (same discipline as `WP8_PROTOCOL.md`'s
"NA convention" note, and for the same reason — `has_result()` treats any
`NA` in a non-key column as a truncated/partial row and would re-run an
already-complete cell): `note` is always a non-empty string (`"-"` on
success); numeric payload columns are never `NA` on a written row — a
degenerate or incomplete cell gets sentinel `0` values plus a `status`
column that is not `"ok"` and a `note` explaining why, never `NA`.

## 8. Cost

See the task report for the measured smoke timing and the full-grid
estimate (cells × seconds for the metrics-only pass, and the added cost of
the DBSCAN/HDBSCAN regeneration pass, dominated by per-cell python
subprocess startup for HDBSCAN).

## 9. Open questions / choices not fully specified by the revision plan

Recorded here per the declare-before-look discipline used throughout this
revision's WP8/WP4 protocols:

1. `excluded` is designated primary and `singleton` secondary for the
   ARI/NMI/AMI headline (§4); the revision plan does not choose between
   them, so this is the author's call, stated in advance rather than picked
   after seeing which one looks better.
2. DBSCAN's `cont` for 87 (§6) is set to `0` as the closest analogue to "no
   contamination hint," not because `0` is a good eps rule for that
   generator — it is not, and the limitation is stated rather than fixed by
   inventing a new DBSCAN variant.
3. HDBSCAN runs via python, not the R package, because WP4 already
   vetted the python 0.8.44 install for exactly this purpose and no
   R `hdbscan` package install has been attempted or verified in this
   environment.
4. k-hat's cross-check against 87's own `n_clusters` column is opportunistic
   (87 is the only WP8 script whose main metrics file already carries an
   equivalent quantity) — 85/86/88 have no such column to cross-check
   against, so their k-hat correctness rests on the shared, tested
   `clustering_metrics()`/counting code alone.
