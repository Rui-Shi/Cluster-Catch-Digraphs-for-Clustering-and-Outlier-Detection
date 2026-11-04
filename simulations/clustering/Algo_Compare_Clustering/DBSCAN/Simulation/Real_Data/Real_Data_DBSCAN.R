source(here::here("data/clustering/Real_Data_Collection.R"))
library(dbscan)
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)

k = 4 # fix the minpts to 5 (4 neighors)

# Eps at the elbow of the sorted k-NN distance curve, found without looking at the plot:
# scale the sorted distances d_(1) <= ... <= d_(n) to [0, 1] on both axes, x_i = (i-1)/(n-1)
# and y_i = (d_(i) - d_(1)) / (d_(n) - d_(1)), and take the point farthest below the
# diagonal, i* = argmax_i (x_i - y_i) (the first index on ties). Eps = d_(i*) * (1 + 1e-9).
# The factor puts the point that defines d_(i*) inside the radius (dbscan counts distances <= Eps),
# so that rounding in the distance cannot decide whether that point is a core point.
knee_eps = function(X, k){
  d = sort(kNNdist(X, k))
  n = length(d)
  x = (seq_len(n) - 1) / (n - 1)
  y = (d - d[1]) / (d[n] - d[1])
  d[which.max(x - y)] * (1 + 1e-9)
}

##### iris data #####
set.seed(1)
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("iris x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
set.seed(1)
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("seeds x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
set.seed(1)
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("user_knowledge x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
set.seed(1)
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("wholesale x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
set.seed(1)
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("breast_cancer x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
set.seed(1)
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("aggregation x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
set.seed(1)
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("asymmetric x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))
      
      

##### R15 data #####
set.seed(1)
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("R15 x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
set.seed(1)
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = knee_eps(X, k)
abline(h = eps, col = "red")  # Visualize a potential eps value at the "elbow"

# conduct clustering
db_result = dbscan(X, eps, minPts = k+1)
db_label = db_result$cluster

# number of clusters:
n_clusters = length(unique(db_label))-1

# ARI
ari = adjustedRandIndex(db_label, Y)

# the average silhouette index
sil = mean(silhouette(db_label, dist(X))[,3])

print(paste("D31 x DBSCAN: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))