# Cluster Catch Digraphs (CCDs) for Outlier Detection and Clustering

Source code, simulation scripts, and benchmark datasets accompanying the
following work on Cluster Catch Digraphs:

- Shi, R., Ceyhan, E., Billor, N. **Outlyingness Scores with Cluster Catch
  Digraphs for Identifying High-Dimensional and Collective Outliers.**
  *Pattern Recognition* (in revision, 2026).
- Shi, R. **Novel Outlier Detection and Clustering Methods Based on
  Cluster Catch Digraphs.** PhD dissertation, Auburn University.

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
| **RU-MCCDs** | Rapid Uniformity-Based CCDs with Mutual Catch Graphs (RK-CCD + MCG) | `simulations/outlier_detection/RU-MCCDs/RU-MCCDs.R` | `simulations/outlier_detection/RU-MCCDs/` |
| **UN-MCCDs** | Uniformity- and Neighbour-based CCDs adapted for outlier detection | `simulations/outlier_detection/UN-MCCDs/UN-MCCD.R` | `simulations/outlier_detection/UN-MCCDs/` |
| **SU-MCCDs** | Shape-adaptive variant of RU-MCCDs | `simulations/outlier_detection/SU-MCCDs/SU-MCCDs.R` | `simulations/outlier_detection/SU-MCCDs/` |
| **SUN-MCCDs** | Shape-adaptive variant of UN-MCCDs | `simulations/outlier_detection/SUN-MCCDs/SUN-MCCD.R` | `simulations/outlier_detection/SUN-MCCDs/` |
| **OOS** *(manuscript)* | Outbound Outlyingness Score: ratio of mean neighbour-density to a point's vicinity density | `methods/outlyingness_scores/Outlyingness_Score.R` | `simulations/outlyingness_scores/RKCCD_OOS_IOS/`, `simulations/outlyingness_scores/UNCCD_OOS_IOS/` |
| **IOS** *(manuscript)* | Inbound Outlyingness Score: inverse of cumulative inbound-neighbour influence; robust to masking | `methods/outlyingness_scores/Outlyingness_Score.R` (shared with OOS) | same as OOS |
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

The 10 datasets cited in the OOS/IOS manuscript are listed first; the
remaining datasets in this folder were collected during the dissertation
work and are kept for reproducibility of related experiments.

| Dataset | File(s) | Used by |
|---|---|---|
| Hepatitis | `Hepatitis_withoutdupl_10_v01.arff` | OOS/IOS, MCCD family |
| Lymphography | `Lymphography_withoutdupl_norm_idf.arff` | OOS/IOS, MCCD family |
| Glass | `glass.mat`, `glass_desc.pdf` | OOS/IOS, MCCD family |
| WBC (Wisconsin Breast Cancer) | `WBC_withoutdupl_norm_v02.arff` | OOS/IOS, MCCD family |
| Vertebral | `vertebral.mat`, `Vertebral_desc.pdf` | OOS/IOS, MCCD family |
| Stamps | `Stamps_withoutdupl_norm_09.arff` | OOS/IOS, MCCD family |
| WDBC (Diagnostic Breast Cancer) | `WDBC_withoutdupl_norm_v02.arff` | OOS/IOS, MCCD family |
| Vowels | `vowels.mat`, `Vowels_desc.pdf` | OOS/IOS, MCCD family |
| Thyroid | `thyroid.mat`, `Thyroid_desc.pdf` | OOS/IOS, MCCD family |
| Wilt | `Wilt_withoutdupl_05.arff` | OOS/IOS, MCCD family |
| PageBlocks | `PageBlocks_withoutdupl_09.arff` | dissertation |
| PenDigits | `PenDigits_withoutdupl_norm_v01.arff` | dissertation |
| Pima | `Pima_withoutdupl_norm_10_v01.arff` | dissertation |
| Shuttle | `Shuttle_withoutdupl_norm_v01.arff` | dissertation |
| Waveform | `Waveform_withoutdupl_v01.arff` | dissertation |
| Ecoli | `ecoli.csv`, `ecoli.mat` | dissertation |
| Yeast | `yeast.data` | dissertation |

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

A small number of `source()` references map to companion methods or
threshold scripts that were not bundled in this repository. They are now
written as `here::here("methods/outlyingness_scores/M-FCCDs/M-FCCDs.R")`
(etc.) — i.e., the path tells you where the file should live if you
obtain it. The affected references are:

| Missing target | Where it would live | Affected simulations |
|---|---|---|
| M-FCCDs.R | `methods/outlyingness_scores/M-FCCDs/M-FCCDs.R` | F-CCD-based OOS/IOS variants |
| M-FNNCCDs.R | `methods/outlyingness_scores/M-FNNCCDs/M-FNNCCDs.R` | FNN-CCD-based OOS/IOS variants |
| NNCCD threshold scripts | `simulations/outlyingness_scores/NNCCD_OOS_IOS/Simulation/...` | NNCCD threshold calibration |

The OOS/IOS results reported in the manuscript do not depend on these
missing files (RKCCD and UNCCD variants are fully reproducible).

---

## Manuscript correspondence (Pattern Recognition / SLDS)

The manuscript cited above benchmarks **OOS** and **IOS** against the
MCCD family (RU, SU, UN, SUN) and against the established detectors
LOF, DBSCAN, MST, ODIN, and Isolation Forest. The exact correspondence
is:

| Manuscript element | Repository location |
|---|---|
| OOS / IOS definition (Section 2) | `methods/outlyingness_scores/Outlyingness_Score.R` |
| Synthetic-data Monte Carlo (Section 3) | `simulations/outlyingness_scores/RKCCD_OOS_IOS/`, `simulations/outlyingness_scores/UNCCD_OOS_IOS/` |
| Cutoff calibration (Section 3.1, Tables 2-3) | `simulations/outlyingness_scores/Uniform_cutoffs/`, `simulations/outlyingness_scores/Gaussian_cutoffs/` |
| Real-data benchmarks (Section 4) | `data/outlier_detection/` + the 10 datasets listed above |
| Baseline detectors (Section 4) | `simulations/outlier_detection/Algo_Compare_OutlierDetection/` |

---

## Citation

If you use this code, please cite both the manuscript and the
dissertation listed at the top of this file.
