# Cluster Catch Digraphs (CCDs) for Outlier Detection and Clustering

Source code, simulation scripts, and benchmark datasets accompanying the
following work on Cluster Catch Digraphs:

- Shi, R., Ceyhan, E., Billor, N. **Outlyingness Scores with Cluster Catch
  Digraphs for Identifying High-Dimensional and Collective Outliers.**
  *Pattern Recognition* (in revision, 2026).
- Shi, R., Ceyhan, E., Billor, N. **Outlier Detection with Cluster Catch
  Digraphs.** arXiv:2409.11596, 2024.
  <https://arxiv.org/abs/2409.11596>
- Shi, R., Ceyhan, E. **Clustering with Uniformity- and Neighbor-Based
  Random Geometric Graphs.** arXiv:2501.06268, 2025.
  <https://arxiv.org/abs/2501.06268>

The repository is structured so that any external researcher can locate
and re-run the exact scripts and data used in the manuscript. Every
`source()`, `load()`, and `setwd()` call has been rewritten to use
[`here::here()`](https://here.r-lib.org/), so the scripts run as-is from
a fresh clone after `install.packages("here")`.

---

## Repository layout

```
.
├── R/                              shared library code (sourced by everything else)
│   ├── ccds/                       core CCD function libraries
│   │   ├── NN-test_quantile/       per-dimension NN test quantile R scripts
│   │   └── RK-test_quantile/       per-dimension RK test quantile R scripts
│   ├── general_functions/          generation, counting, ratio helpers
│   ├── NN-test_quantile/           cluster-tree NN test quantile R scripts
│   └── RK-test_quantile/           cluster-tree RK test quantile R scripts
│
├── methods/                        method definitions / driver wrappers
│   ├── clustering/                 RK-, UN-, KS-CCDs clustering wrappers
│   ├── outlier_detection/          RU-, SU-, UN-, SUN-MCCDs (arXiv:2409.11596)
│   └── outlyingness_scores/        OOS / IOS scoring wrappers (the manuscript's contribution)
│
├── simulations/                    Monte-Carlo experiments grouped by task
│   ├── clustering/                 RK_CCDs, UN_CCDs, KS_CCDs + baselines
│   ├── outlier_detection/          RU-, SU-, UN-, SUN-MCCDs, D-MCGs + baselines
│   └── outlyingness_scores/        OOS/IOS simulations + cutoff-calibration runs
│
├── data/                           benchmark real-world datasets
│   ├── clustering/                 16 datasets used for clustering benchmarks
│   └── outlier_detection/          17 datasets used for outlier-detection benchmarks
│
└── README.md                       this file
```

---

## Mapping: scripts → algorithms

### Algorithms proposed in this body of work

| Algorithm | Description | Method file (definition) | Simulation tree |
|---|---|---|---|
| **RK-CCDs** | Clustering via Ripley's *K*–function spatial-randomness test | `methods/clustering/RK-CCDs.R` | `simulations/clustering/RK_CCDs/` |
| **UN-CCDs** | Clustering via Nearest-Neighbour Distance (NND) spatial-randomness test | `methods/clustering/UN-CCDs.R` | `simulations/clustering/UN_CCDs/` |
| **KS-CCDs** | Clustering via Kolmogorov–Smirnov-type density statistic | `methods/clustering/KS-CCDs.R` | `simulations/clustering/KS_CCDs/` |
| **D-MCGs** | Density-based Mutual Catch Graphs for single-cluster outlier detection | (defined inside the simulation tree) | `simulations/outlier_detection/D-MCGs/` |
| **RU-MCCDs** | Rapid Uniformity-Based CCDs with Mutual Catch Graphs (RK-CCD + MCG) | `methods/outlier_detection/RU-MCCDs.R` | `simulations/outlier_detection/RU-MCCDs/` |
| **UN-MCCDs** | Uniformity- and Neighbour-based CCDs adapted for outlier detection | `methods/outlier_detection/UN-MCCD.R` | `simulations/outlier_detection/UN-MCCDs/` |
| **SU-MCCDs** | Shape-adaptive variant of RU-MCCDs | `methods/outlier_detection/SU-MCCDs.R` | `simulations/outlier_detection/SU-MCCDs/` |
| **SUN-MCCDs** | Shape-adaptive variant of UN-MCCDs | `methods/outlier_detection/SUN-MCCD.R` | `simulations/outlier_detection/SUN-MCCDs/` |
| **OOS** | Outbound Outlyingness Score: ratio of mean neighbour-density to a point's vicinity density | `methods/outlyingness_scores/Outlyingness_Score.R` | `simulations/outlyingness_scores/RKCCD_OOS_IOS/`, `simulations/outlyingness_scores/UNCCD_OOS_IOS/` |
| **IOS** | Inbound Outlyingness Score: inverse of cumulative inbound-neighbour influence; robust to masking | `methods/outlyingness_scores/Outlyingness_Score.R` (shared with OOS) | same as OOS |
| **RKCCD-OOS / RKCCD-IOS** | OOS / IOS computed on the RK-CCD covering-ball digraph | `methods/outlyingness_scores/RKCCD_OOS_IOS.R` | `simulations/outlyingness_scores/RKCCD_OOS_IOS/` |
| **UNCCD-OOS / UNCCD-IOS** | OOS / IOS computed on the UN-CCD covering-ball digraph | `methods/outlyingness_scores/UNCCD_OOS_IOS.R` | `simulations/outlyingness_scores/UNCCD_OOS_IOS/` |

### Baseline methods benchmarked against the proposals

| Baseline | Used for | Implementation |
|---|---|---|
| DBSCAN | clustering & outlier detection | `simulations/clustering/Algo_Compare_Clustering/DBSCAN/`, `simulations/outlier_detection/Algo_Compare_OutlierDetection/DBSCAN/` |
| K-means++ | clustering | `simulations/clustering/Algo_Compare_Clustering/K-means++/` |
| Louvain clustering | clustering | `simulations/clustering/Algo_Compare_Clustering/Louvain/` |
| Spectral clustering | clustering | `simulations/clustering/Algo_Compare_Clustering/Spectral_Clustering/` |
| Minimum Spanning Tree (MST) | clustering & outlier detection | `simulations/clustering/Algo_Compare_Clustering/MST/`, `simulations/outlier_detection/Algo_Compare_OutlierDetection/MST/` |
| Local Outlier Factor (LOF) | outlier detection | `simulations/outlier_detection/Algo_Compare_OutlierDetection/LOF/` |
| Isolation Forest (iForest) | outlier detection | `simulations/outlier_detection/Algo_Compare_OutlierDetection/Isolation Forest/` |
| ODIN (in-degree number) | outlier detection | `simulations/outlier_detection/Algo_Compare_OutlierDetection/ODIN/` |

### Cutoff calibration (OOS / IOS thresholds)

| Cluster-process model | Folder |
|---|---|
| Uniform clusters | `simulations/outlyingness_scores/Uniform_cutoffs/` |
| Gaussian clusters | `simulations/outlyingness_scores/Gaussian_cutoffs/` |

Each folder contains threshold-calibration scripts for FCCD, FNNCCD,
NNCCD, RKCCD, and UNCCD variants.

---

## Mapping: datasets → benchmark experiments

### Outlier-detection real datasets — `data/outlier_detection/`

The 10 datasets used in the OOS/IOS manuscript are listed first; the
remaining datasets are archived alongside them and used in extended
benchmarks for the MCCD-family methods (arXiv:2409.11596).

**Used in the OOS/IOS manuscript:**

| Dataset | File(s) |
|---|---|
| Hepatitis | `Hepatitis_withoutdupl_10_v01.arff` |
| Lymphography | `Lymphography_withoutdupl_norm_idf.arff` |
| Glass | `glass.mat`, `glass_desc.pdf` |
| WBC (Wisconsin Breast Cancer) | `WBC_withoutdupl_norm_v02.arff` |
| Vertebral | `vertebral.mat`, `Vertebral_desc.pdf` |
| Stamps | `Stamps_withoutdupl_norm_09.arff` |
| WDBC (Diagnostic Breast Cancer) | `WDBC_withoutdupl_norm_v02.arff` |
| Vowels | `vowels.mat`, `Vowels_desc.pdf` |
| Thyroid | `thyroid.mat`, `Thyroid_desc.pdf` |
| Wilt | `Wilt_withoutdupl_05.arff` |

**Additional datasets (extended MCCD benchmarks):**

| Dataset | File(s) |
|---|---|
| PageBlocks | `PageBlocks_withoutdupl_09.arff` |
| PenDigits | `PenDigits_withoutdupl_norm_v01.arff` |
| Pima | `Pima_withoutdupl_norm_10_v01.arff` |
| Shuttle | `Shuttle_withoutdupl_norm_v01.arff` |
| Waveform | `Waveform_withoutdupl_v01.arff` |
| Ecoli | `ecoli.csv`, `ecoli.mat` |
| Yeast | `yeast.data` |

`data/outlier_detection/RealData_Collection.R` and `LoadData.R` are
loader scripts that read each dataset into a uniform format. Per-dataset
descriptive PDFs sit alongside the data files.

### Clustering real datasets — `data/clustering/`

| Dataset | File(s) |
|---|---|
| Iris | `iris.data`, `Iris.pdf` |
| Wisconsin Breast Cancer (clustering) | `breast-cancer-wisconsin.data`, `Breast Cancer Wisconsin.pdf` |
| Lymphography | `lymphography.data`, `Lymphography .pdf` |
| Hepatitis | `hepatitis.data` |
| Glass | `glass.txt`, `glass_desc.pdf` |
| Yeast | `yeast.data`, `Yeast.pdf` |
| Seeds | `seeds_dataset.txt`, `seeds.pdf` |
| Wholesale Customers | `Wholesale customers data.csv`, `Wholesale customers.pdf` |
| Wine Quality | `winequality-red.csv`, `winequality-white.csv`, `wine_quality.pdf` |
| Blood Transfusion | `transfusion.data`, `Blood Transfusion Service Center.pdf` |
| Travel Reviews | `Travel Reviews.pdf` |
| User Knowledge Modeling | `user_knowlege.xls`, `User Knowledge Modeling.pdf` |
| Website Phishing | `PhishingData.arff`, `Website Phishing.pdf` |
| Waveform | `Waveform_withoutdupl_v01.arff`, `Waveform_desc.pdf` |
| Synthetic 2-D shapes | `D31.txt`, `R15.txt`, `aggregation.txt`, `compound.txt`, `jain.txt`, `asymmetric.txt`, `asymmetric_label.txt` |

`data/clustering/Real_Data_Collection.R` and
`Real_Data_Collection_full.R` load these into a uniform format.

---

## Reproducing the simulations

### Prerequisites

Install R (>= 4.1) and the `here` package, which is used by every script
to resolve paths relative to the repository root:

```r
install.packages("here")
```

The other R packages used across the simulations include `mclust`,
`igraph`, `parallel`, `doParallel`, `MASS`, `dbscan`, `isotree`,
`solitude`, `cluster`, `mvtnorm`, `spatstat`, and `farff`. Install
on demand as scripts complain.

### Running a single simulation

From any working directory, after cloning the repo:

```bash
cd Cluster-Catch-Digraphs-for-Clustering-and-Outlier-Detection
Rscript simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/3d/3d_n500_eps05.R
```

`here::here(...)` resolves paths against the repository root (anchored
by the `.git` directory), so this is independent of the current working
directory at invocation time.

### Running an entire experiment family

Each simulation folder contains the per-dimension and per-sample-size
scripts together with their SLURM `.out` job logs from the original HPC
runs. The folder layout typically reads:

```
<method>/Simulations/<distribution>/<dim>d/<dim>d_<setting>_n<size>.R
```

so a depth-first walk in `xargs Rscript` will reproduce the entire grid.
Expect long run times (hours to days) for the high-dimensional settings.

### Reading the precomputed quantile tables

`R/NN-test_quantile/`, `R/RK-test_quantile/`, and the analogous
subfolders under `R/ccds/` hold the simulation drivers that produce the
quantile-table `.RData` files referenced by the simulations. The
`.RData` files themselves are large and were not committed; rerun the
relevant `*.R` driver to regenerate them, then `load()` calls in the
simulation scripts will resolve.

---

## Notes on bundled vs missing artifacts

The MCCD simulation trees originally `source()`d files named `M-FCCDs.R`,
`M-FNNCCDs.R`, and `M-CCDs.R`. These are earlier names for the
outlier-detection method wrappers and have been re-pointed at their
present-day equivalents:

| Old name | Present-day file |
|---|---|
| `M-CCDs.R` | `methods/outlier_detection/RU-MCCDs.R` |
| `M-FCCDs.R` | `methods/outlier_detection/SU-MCCDs.R` |
| `M-FNNCCDs.R` | `methods/outlier_detection/SUN-MCCD.R` |

The only `source()` references that remain unmappable are NNCCD
threshold scripts referenced from deprecated `_old` simulation folders:

| Missing target | Where it would live | Affected simulations |
|---|---|---|
| NNCCD threshold scripts | `simulations/outlyingness_scores/NNCCD_OOS_IOS/Simulation/...` | only deprecated `_old/` runs under `simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/` |

The OOS/IOS results reported in the manuscript do not depend on these
missing files; the active RKCCD- and UNCCD-based simulation trees are
fully reproducible.

---

## Citation

If you use this code, please cite the manuscript and the relevant arXiv
preprints listed at the top of this file.
