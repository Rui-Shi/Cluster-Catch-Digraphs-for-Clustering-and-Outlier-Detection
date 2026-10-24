# EXPLORATORY, NOT ADOPTED. Tried on 2026-10-23 and not used in the paper.
# Replaces the fixed threshold 2 of Real_Data_MST_fix_cutoff.R by a label-free per-dataset rule that
# was specified before any run: the threshold that maximizes the average silhouette index of the MST
# partition over seq(1.00, 3.00, by = 0.01), first maximum on ties, Sil = -Inf for a partition with
# fewer than 2 clusters (or k = n). It over-segments: the chosen threshold is 1.03, 1.03, 1.04, 1.00,
# 1.00, 1.00, 1.04, 1.00, 1.00 on iris, seeds, user_knowledge, wholesale, breast_cancer, aggregation,
# asymmetric, R15, D31 (the lower grid edge on 5 of 9) and the partitions have 74, 101, 130, 227, 226,
# 404, 465, 282, 1495 clusters. The paper keeps thresh = 2 (Real_Data_MST_fix_cutoff.R), Zahn's ratio
# criterion. Results: revision_experiments/results/tr3/real_baselines_fixed_mstrule_MST_run1.csv
# (and _run2, _params), README.md of revision_experiments/tr3, "MST threshold rule (2026-10-23)".
# Run with the driver at HEAD: the driver prints the block's line (threshold included) but does not
# write the chosen threshold to a params file.
source(here::here("data/clustering/Real_Data_Collection.R"))
source(here::here("simulations/clustering/Algo_Compare_Clustering/MST/MST_Clustering.R"))
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)
library(foreach)
library(doParallel)

cores = detectCores()
cl = makeCluster(cores)
registerDoParallel(cl)

# Threshold rule (authors' decision, 2026-10-23, fixed before any run). On each data set, the
# threshold is the grid value that maximizes the average silhouette index of the MST partition,
# mean(silhouette(label, dist(X))[,3]). Ties go to the first maximum (which.max). Class labels are
# not used. A partition with fewer than 2 clusters (or with as many clusters as points) has no
# silhouette; it gets -Inf, so it is chosen only if every grid value gives such a partition.
# Until 2026-10-23 the script used thresh = 2, "chosen after a trial-and-error process".
thresh_grid = seq(1.00, 3.00, by = 0.01)  # the ratio to cut a edge when comparing its adjacent edges

# mst_prep() and mst_cut() split MST_Clustering() in two: the MST and the adjacent-edge averages
# (the slow part, computed once per data set) and the pruning at one threshold. The lines are
# copied from MST_Clustering(); mst_cut(mst_prep(X), t) is identical to MST_Clustering(X, t).
mst_prep = function(data){
  n = dim(data)[1]

  # construct a full graph
  g = graph_from_adjacency_matrix(as.matrix(dist(data)), mode = "undirected", weighted = TRUE)

  # Calculate MST
  mst = mst(g)

  # the end points of each edge
  end_pts = ends(mst, E(mst))

  # Calculate edge lengths
  edge_lengths = E(mst)$weight

  # Calculate average length of adjacent edges for each edge
  adj_avg = sapply(1:length(edge_lengths), function(i){
    end_temp = end_pts[i,]
    # the two end points of current edge
    end1 = end_temp[1]
    end2 = end_temp[2]
    adj_ind = apply(end_pts, 1, function(row){
      sum(row == end1 | row == end2) == 1})
    return(mean(E(mst)$weight[adj_ind]))
  })
  return(list(mst = mst, adj_avg = adj_avg))
}

mst_cut = function(p, thresh){
  # Remove edges significantly longer than adjacent average
  remove_ind <- which(E(p$mst)$weight > thresh*p$adj_avg)
  mst_pruned <- delete_edges(p$mst, remove_ind)

  # Extract clusters
  cls = components(mst_pruned)$membership
  return(cls)
}

# the average silhouette index of a partition; -Inf where silhouette() is undefined
mst_sil = function(label, d_X){
  k = length(unique(label))
  if (k < 2 || k >= length(label)) return(-Inf)
  return(mean(silhouette(label, d_X)[,3]))
}
##### iris data #####
set.seed(1)
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("iris x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### seeds data #####
set.seed(1)
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("seeds x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### user_knowledge data #####
set.seed(1)
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("user_knowledge x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### wholesale data #####
set.seed(1)
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("wholesale x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### breast_cancer #####
set.seed(1)
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("breast_cancer x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### aggregation data #####
set.seed(1)
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("aggregation x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### asymmetric data #####
set.seed(1)
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("asymmetric x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### R15 data #####
set.seed(1)
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("R15 x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


##### D31 data #####
set.seed(1)
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# choose the threshold on this data set: the first grid value with the largest average silhouette index
d_X = dist(X)
p_X = mst_prep(X)
sil_grid = sapply(thresh_grid, function(th) mst_sil(mst_cut(p_X, th), d_X))
thresh_op = thresh_grid[which.max(sil_grid)]
# the split MST code gives the labels of MST_Clustering()
stopifnot(identical(mst_cut(p_X, thresh_op), MST_Clustering(X, thresh_op)))

# construct MST for clustering with the chosen threshold
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mst_sil(MST_label, dist(X))

print(paste("D31 x MST: the threshold is", thresh_op, "(grid index", which.max(sil_grid), "of", length(thresh_grid), ", grid values without silhouette:", sum(sil_grid == -Inf), ") the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


stopCluster(cl)