library(ClusterR)
source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)
library(foreach)
library(doParallel)

##### iris data #####
k=3
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("iris x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
k=3
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("seeds x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
k=4
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("user_knowledge x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
k=6
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("wholesale x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
k=2
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("breast_cancer x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
k=7
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("aggregation x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
k=5
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("asymmetric x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### R15 data #####
k=15
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("R15 x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
k=31
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# construct kmean++ for clustering
kmeanpp_label = KMeans_rcpp(X, clusters = k, initializer = "kmeans++", verbose = F)$clusters

# number of clusters:
n_clusters = length(unique(kmeanpp_label))

# ARI
ari = adjustedRandIndex(kmeanpp_label, Y)

# the average silhouette index
sil = mean(silhouette(kmeanpp_label, dist(X))[,3])

print(paste("D31 x kmeanpp: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))