# WP4 protocol — modern baselines on real data

**Declared 2026-09-05, before any competitor was run on any data set.**

WP4 answers AE.3, R1.3, R1.5, R3.3 and R6.5: the five published baselines (LOF,
DBSCAN, MST, ODIN, iForest) all predate 2010 except iForest, none of them is a
density-clustering method with varying-density support, none is a modern
tabular deep detector, and — the sharpest of the objections, R3.3 — none is a
plain mutual-kNN detector, so the paper cannot currently show that its gain
comes from the mutual *catch* mechanism rather than from reciprocity as such.

This file fixes the data, the competitors, their settings, and the
thresholding rule in advance, in the same spirit as
`BENCHMARK_EXPANSION_RULE.md`. Nothing here may be changed after a result is
seen. Any deviation forced by an implementation failure is recorded as an
appended, dated note rather than an edit to the declarations above it.

## 1. Data

The same sixteen data sets as the Section 6 real-data study, admitted by the
pre-declared inclusion rule in `BENCHMARK_EXPANSION_RULE.md`:

hepatitis, lymphography, glass, WBC, vertebral, ecoli, stamps, WDBC, pima,
Shuttle, vowels, PenDigits, waveform, thyroid, pageblocks, wilt.

**Preprocessing is exactly what the R loader
`data/outlier_detection/RealData_Collection.R` applies**, with no additional
scaling, imputation, or feature selection. The competitors therefore see the
identical feature matrices the four proposed detectors saw. Per data set the
loader applies:

| regime | data sets |
|---|---|
| classical z-score (`scale`) | glass, hepatitis, Shuttle (`shuffle`) |
| robust median/MADN (`scale_R`) | vertebral, vowels, thyroid, PenDigits, waveform, wilt, pageblocks |
| none (source `.arff` already scaled to [0,1]) | ecoli, stamps, WBC, WDBC, lymphography, pima |

Deduplication (`distinct`) and the two subsampling steps (PenDigits to 3180
regulars, pageblocks to 4285 regulars, both at `set.seed(123)`) are the
loader's own and are inherited unchanged. The loader's R object names are
used verbatim; note **Shuttle is the object `shuffle`**.

Export is by `tr1/80_export_datasets_wp4.R`, which sources
`shared/harness.R` (read-only) and calls `load_real_dataset()`. Each data set
is written to `results/tr1/wp4/data/<dataset>.csv` with the preprocessed
feature columns followed by a final `label` column. The export is gated on an
exact match of (n, d, n_outliers) against Table `tab:Real_Data` of
`CCD_OutlierDetection_Neurocomputing.tex`; a single disagreement aborts the
run.

### Label polarity — two opposite conventions, kept explicit

| where | convention |
|---|---|
| this repo, `RealData_Collection.R`, `harness.R`, every exported `label` column, `evaluate()` | **1 = regular, 0 = outlier** |
| PyOD `labels_`, sklearn/hdbscan noise label `-1`, every `is_outlier` column written by WP4 | **1 = outlier, 0 = regular** |

Nothing in WP4 silently converts between them. Score files carry no labels at
all; native-label files carry a column named `is_outlier`, in PyOD polarity.
The conversion happens once, in R, where the metrics are computed.

### A known loader defect, and why it does not bite here

`RealData_Collection.R` sorts glass by column 9 and ecoli by column 6 — both
feature columns, where the labels sit in columns 10 and 8. So the rows of
those two sets are not ordered regulars-then-outliers. WP4 never relies on row
order: the exported `label` column travels row-aligned with the features, every
score file preserves data row order, and `evaluate()` reorders `(Y, score)`
jointly before counting. Positional scoring is not used anywhere in WP4.

## 2. Methods and settings

Every method is fit on the feature matrix only. Labels are never passed to any
`fit`.

**Score polarity for every score-emitting method: higher = more outlying.**
Where a method's natural quantity is a density, the score is its negation.

| method | family | package | settings | replication |
|---|---|---|---|---|
| ECOD | parameter-free, 2022 | `pyod.models.ecod.ECOD` | pyod defaults | deterministic, one fit (seed slot 0) |
| COPOD | copula-based, 2020 | `pyod.models.copod.COPOD` | pyod defaults | deterministic, one fit (seed slot 0) |
| DIF | deep isolation forest, 2023 | `pyod.models.dif.DIF` | pyod defaults | seeds 1–5 |
| LUNAR | modern deep tabular, graph-based, 2022 | `pyod.models.lunar.LUNAR` | pyod defaults | seeds 1–5 |
| HDBSCAN + GLOSH | density clustering, varying densities | `hdbscan.HDBSCAN` | defaults (`min_cluster_size=5`, `min_samples=None`) | deterministic, one fit |
| OPTICS | density clustering, ordering-based | `sklearn.cluster.OPTICS` | defaults (`min_samples=5`, `xi=0.05`) | deterministic, one fit |
| mutual-kNN | reciprocal neighbourhood | own code | k ∈ {5, 10, 15, 20, 30} | deterministic, one fit per k |
| SNN | shared-nearest-neighbour density | own code | k ∈ {5, 10, 15, 20, 30} | deterministic, one fit per k |

Details that are not defaults, or that need stating:

- **Seeding.** For DIF and LUNAR, `random_state=seed` is passed to the
  constructor *and* numpy, `random`, and `torch` global RNGs are seeded before
  each fit, following `tr2/08_wp6_pyod_baselines.py`: pyod's `random_state` is
  not known to cover every source of randomness in every model in this version.
- **`contamination`.** Left at pyod's default 0.1 for all four pyod models.
  It does not influence fitting or `decision_scores_` — pyod consumes it only
  downstream, in `BaseDetector._process_decision_scores`, to set
  `threshold_`/`labels_`. WP4 saves raw `decision_scores_` and derives labels
  separately (§3), so the constructor value is inert here.
- **GLOSH.** The HDBSCAN score is `outlier_scores_`, hdbscan's own GLOSH
  implementation (Campello et al. 2015). Higher = more outlying already; used
  as-is. Non-finite entries, if any, are recorded in the fit log and left in
  place for R to handle rather than silently patched.
- **HDBSCAN native labels.** `labels_ == -1` is the threshold-free variant:
  `is_outlier = 1` for noise points, 0 otherwise.
- **OPTICS score.** `reachability_`, which is `inf` for the first point in the
  ordering and for any point never reached. **Declared rule: every non-finite
  reachability is replaced by (max finite reachability) × 1.01**, so those
  points are the most outlying and the ordering among the finite values is
  untouched. The count of replaced entries is logged per data set.
- **OPTICS native labels.** `labels_ == -1` is the threshold-free variant.
- **mutual-kNN.** For each i, kNN(i) is the k nearest neighbours in Euclidean
  distance excluding i itself. The mutual degree is
  m_i = #{ j : j ∈ kNN(i) and i ∈ kNN(j) }, and the score is **k − m_i**, so a
  point no neighbour reciprocates scores k and a fully reciprocated point
  scores 0. Threshold-free variant: `is_outlier = 1` iff m_i = 0.
- **SNN.** Ertöz–Steinbach–Kumar style. For each i,
  density_i = Σ_{j ∈ kNN(i)} |kNN(i) ∩ kNN(j)|, and the score is **−density_i**.
  No threshold-free variant is defined for SNN.
- **Ties in kNN.** Broken by index order, which is what
  `sklearn.neighbors.NearestNeighbors` does. Both neighbour-based detectors use
  one shared exact kNN query per (data set, k), so they see identical
  neighbourhoods.

### The oracle concession on k, stated in advance

For mutual-kNN and SNN the **best k per data set is selected on F₂ using the
true labels, and the selected cell is what gets reported**. This is a
deliberately generous setting for the competitor: no unsupervised rule picks k,
so these two baselines are given ground-truth information that none of the
proposed methods receives. It is declared here so that it is reported as an
oracle in the manuscript and the response letter, not presented as a fair
tuning-free comparison. If a proposed method still beats an oracle-k mutual-kNN
detector, R3.3's question — mutual catch, or mere reciprocity — has an answer
that does not depend on how k was chosen.

## 3. Thresholding and metrics

WP4 computes **no metrics**. It emits raw scores and, where applicable, native
labels. Labels and metrics are derived downstream in R, by the same
`evaluate()` used for the nine methods already in the study, so that
TPR/TNR/BA/F₂ are computed identically for competitors and proposed methods.

Two thresholds, both declared now:

- **T1, contamination = 0.1.** PyOD's library default, identical for every
  data set, label-free. **This is the main-table setting.**
- **T2, contamination = the data set's true outlier rate.** An oracle,
  reported only as a sensitivity, never in the headline comparison.

In both cases the label rule is pyod's own:
`threshold = percentile(score, 100 × (1 − contamination))` and
`is_outlier = score > threshold`.

The threshold-free variants (HDBSCAN noise, OPTICS noise, mutual-kNN m_i = 0)
take no contamination argument and are reported as they come.

## 4. Outputs

```
results/tr1/wp4/data/<dataset>.csv                       features + label (1 = regular, 0 = outlier)
results/tr1/wp4/data/manifest.csv                        n, d, n_outliers, scaling regime, table check
results/tr1/wp4/scores/<dataset>_<method>[_k<k>][_seed<s>].csv    one column `score`, higher = more outlying
results/tr1/wp4/scores/<dataset>_<method>[...]_labels.csv         one column `is_outlier` (1 = outlier), native variants only
results/tr1/wp4/fit_log.csv                              one appended row per finished cell
results/tr1/wp4/fit_errors.log                           tracebacks
```

Row order in every score and label file is the data row order of the
corresponding `<dataset>.csv`.

`fit_log.csv` columns: `dataset, method, k, seed, n, d, fit_seconds, status,
error`. **One line is appended per finished cell**, never buffered to the end,
because the volume this repo sits on drops off the bus under sustained write
load. A cell already logged with `status == "ok"` is skipped on restart, so the
driver is resumable and can be called repeatedly under a short wall-clock cap.

**Fit times are indicative only.** Another R job was running on this machine
throughout, so `fit_seconds` is not a clean timing measurement and must not be
used for the WP6 scalability claims.

## 5. Environment

Recorded 2026-09-05 from `revision_experiments/.venv` (python at the env root,
not `Scripts/`), before any WP4 run:

| package | version |
|---|---|
| python | 3.11.15 (conda-forge, MSC v.1944, 64-bit) |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| scikit-learn | 1.9.0 |
| torch | 2.13.0+cpu |
| pyod | 3.6.1 |
| hdbscan | 0.8.44 |

`hdbscan` was installed for WP4 on 2026-09-05 with
`.venv/python.exe -m pip install hdbscan --only-binary :all:`, which resolved
the prebuilt wheel `hdbscan-0.8.44-cp311-cp311-win_amd64.whl`; no compilation
was needed and no other package changed version. This closes the gap recorded
in the revision plan, where GLOSH was unavailable because sklearn 1.9's own
`HDBSCAN` exposes no `outlier_scores_`.

pyod constructor defaults as reported by `get_params()` in this environment,
for the record:

- `ECOD`: `contamination=0.1, n_jobs=1`
- `COPOD`: `contamination=0.1, n_jobs=1`
- `DIF`: `batch_size=1000, contamination=0.1, device=cpu, hidden_activation='tanh', hidden_neurons=[500,100], max_samples=256, n_ensemble=50, n_estimators=6, random_state=None, representation_dim=20, skip_connection=False`
- `LUNAR`: `algorithm='auto', contamination=0.1, epsilon=0.1, leaf_size=30, lr=0.001, metric='minkowski', model_type='WEIGHT', n_epochs=200, n_neighbours=5, negative_sampling='MIXED', p=2, proportion=1.0, random_state=None, scaler=None, val_size=0.1, wd=0.1`
- `hdbscan.HDBSCAN`: `min_cluster_size=5, min_samples=None, metric='euclidean', cluster_selection_method='eom', alpha=1.0, algorithm='best'`
- `sklearn.cluster.OPTICS`: `min_samples=5, xi=0.05, cluster_method='xi', metric='minkowski', p=2, max_eps=inf`

## 6. Cell count

Per data set: ECOD 1 + COPOD 1 + DIF 5 + LUNAR 5 + HDBSCAN 1 + OPTICS 1 +
mutual-kNN 5 + SNN 5 = **24 cells**, so **384 cells** over sixteen data sets.
Work is ordered smallest n first so that partial results are usable early. A
method that fails on a data set is logged with its error and the sweep
continues; a failure is reported, never dropped.
