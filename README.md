# Cluster Catch Digraphs for Outlier Detection and Clustering

This repository contains the source code and simulation scripts for the dissertation, "Novel Outlier Detection and Clustering Methods Based on Cluster Catch Digraphs."

## Algorithms

The repository implements a family of methods for clustering, outlier detection, and scoring based on Cluster Catch Digraphs (CCDs).

### Clustering Method
* **UN-CCDs (Uniformity- and Neighbor-based CCDs):** A clustering method that uses the Nearest-Neighbor Distance (NND) to test for spatial randomness. It is designed to be more effective in high-dimensional spaces compared to previous CCD methods.

### Outlier Detection Methods
* **RU-MCCD (Rapid Uniformity-Based CCDs with Mutual catch graph):** This method first uses RK-CCDs to form initial clusters and then applies Mutual Catch Graphs (MCGs) to find outliers in low-density regions.
* **UN-MCCD (Uniformity- and Neighbor-based CCDs with Mutual catch graph):** An adaptation of the UN-CCD clustering algorithm specifically tailored for outlier detection tasks.
* **SU-MCCD & SUN-MCCD:** Shape-adaptive versions of the RU-MCCD and UN-MCCD methods. They are designed to handle datasets with arbitrarily shaped clusters by using a more flexible cluster coverage approach.

### Outlyingness Scores
* **OOS (Outbound Outlyingness Score):** A score that measures how much a point deviates from its neighbors based on a comparison of local vicinity densities.
* **IOS (Inbound Outlyingness Score):** A score that quantifies a point's degree of being an outlier by measuring the inverse of the "cumulative influence" it receives from its neighbors. This score is highly robust to the masking problem, where groups of outliers obscure one another.