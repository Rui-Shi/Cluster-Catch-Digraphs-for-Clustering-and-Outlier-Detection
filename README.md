# Cluster Catch Digraphs (CCDs) for Outlier Detection and Clustering

This repository contains the source code and simulation scripts for the dissertation, "Novel Outlier Detection and Clustering Methods Based on Cluster Catch Digraphs."

## Algorithms

The repository implements a family of methods for clustering, outlier detection, and scoring based on Cluster Catch Digraphs (CCDs).

### Clustering Method

* **RK-CCDs (CCDs based on Ripley's $K$ function):** A clustering method that uses the Ripley's $K$ function for for spatial randomness test. It is designed to be more effective in high-dimensional spaces compared to previous CCD methods.

* **UN-CCDs (Uniformity- and Neighbor-based CCDs):** A clustering method that uses the Nearest-Neighbor Distance (NND) for spatial randomness test. It is designed to be more effective in high-dimensional spaces compared to previous CCD methods.

* **KS-CCDs (The CCDs based on a Kolmogorov-Smirnov (KS)-type statistic):** A clustering method that uses a density-based metric to determine the volume of each covering ball, versatile across different types of datasets. However, it need an input density parameter $\delta$.

### Outlier Detection Methods
* **RU-MCCD (Rapid Uniformity-Based CCDs with Mutual catch graph):** This method first uses RK-CCDs to form initial clusters and then applies Mutual Catch Graphs (MCGs) to find outliers in low-density regions.
* **UN-MCCD (Uniformity- and Neighbor-based CCDs with Mutual catch graph):** An adaptation of the UN-CCD clustering algorithm specifically tailored for outlier detection tasks.
* **SU-MCCD & SUN-MCCD:** Shape-adaptive versions of the RU-MCCD and UN-MCCD methods. They are designed to handle datasets with arbitrarily shaped clusters by using a more flexible cluster coverage approach.

### Outlyingness Scores
* **OOS (Outbound Outlyingness Score):** A score that measures how much a point deviates from its neighbors based on a comparison of local vicinity densities.
* **IOS (Inbound Outlyingness Score):** A score that quantifies a point's degree of being an outlier by measuring the inverse of the "cumulative influence" it receives from its neighbors. This score is highly robust to the masking problem, where groups of outliers obscure one another.

## The Structure of the Repository 

We introduce the details of the structure of this repository.

### CCDs_Clustering

* **KS_CCDs**: The simulations study for `KS-CCDs` on clustering.
* **RK_CCDs**: The simulations study for `RK-CCDs` on clustering.
* **UN_CCDs**: The simulations study for `UN-CCDs` on clustering.
* **Algo_Compare_Clustering**: The simulation study for existing clustering method, including `DBSCAN`, `K-means++`, `Louvain clustering`, `Minimal Spanning Tree`, and `Spectral clustering`.

### CCDs_Outlier_Detection

* **D-MCGs**: The simulation study for `D-MCGs` on outlier detection with single cluster.
* **RU-MCCDs**: The simulation study for `D-MCGs` on outlier detection.
* **SU-MCCDs**: The simulation study for `SU-MCCDs` on outlier detection.
* **SUN-MCCDs**: The simulation study for `SUN-MCCDs` on outlier detection.
* **UN-MCCDs**: The simulation study for `UN-MCCDs` on outlier detection.

### CCDs_Outlyingness_Score
* **RKCCD_OOS_IOS**: The simulation study for `OOS` and `IOS` based on `RK-CCDs` on outlier detection.
* **UNCCD_OOS_IOS**: The simulation study for `OOS` and `IOS` based on `UN-CCDs` on outlier detection.
* **Gaussian_cutoffs** and **Uniform_cutoffs**: The experiments to find the optimal threshold of `OOS` and `IOS` for outlier detection under different types of datasets and dimensions.

### Necessary_Functions
The essential scripts to conduct CCD-based methods.

### Real_Datasets_Clustering, Real_Datasets_Outlier_Detection
The benchmark real datasets for clustering and outlier detection.



