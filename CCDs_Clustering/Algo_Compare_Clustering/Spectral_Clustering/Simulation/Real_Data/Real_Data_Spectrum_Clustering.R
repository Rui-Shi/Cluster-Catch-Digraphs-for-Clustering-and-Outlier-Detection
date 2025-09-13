library(kernlab)
library(igraph)
library(ClusterR)
source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(mclust)
library(cluster)
library(foreach)
library(doParallel)

sigma=1

##### iris data #####
k=3 # number of clusters
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# Gaussian kernel similarity matrix.
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("iris x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
k=3 # number of clusters
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# Gaussian kernel similarity
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("seeds x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
k=4 # number of clusters
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# Gaussian kernel similarity
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("user_knowledge x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
k=3 # number of clusters
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# Gaussian kernel similarity
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("wholesale x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
k=2 # number of clusters
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# Gaussian kernel similarity
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("breast_cancer x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
k=7 # number of clusters
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# Gaussian kernel similarity matrix.
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("aggregation x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
k=5 # number of clusters
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# Gaussian kernel similarity matrix.
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("asymmetric x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### R15 data #####
k=15 # number of clusters
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# Gaussian kernel similarity matrix.
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("R15 x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
k=31 # number of clusters
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# Gaussian kernel similarity matrix.
similarity_matrix = as.matrix(dist(X))
# similarity_matrix = exp(-similarity_matrix^2 / (2 * sigma^2))

# Spectral clustering
clusters = specc(similarity_matrix, centers = k)
spectrum_label = clusters@.Data

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("D31 x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))