library(kernlab)
library(igraph)
library(ClusterR)
source(here::here("data/clustering/Real_Data_Collection.R"))
library(mclust)
library(cluster)
library(foreach)
library(doParallel)

##### iris data #####
set.seed(1)
k=3 # number of clusters
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("iris x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### seeds data #####
set.seed(1)
k=3 # number of clusters
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("seeds x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### user_knowledge data #####
set.seed(1)
k=4 # number of clusters
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("user_knowledge x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### wholesale data #####
set.seed(1)
k=3 # number of clusters
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("wholesale x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### breast_cancer data #####
set.seed(1)
k=2 # number of clusters
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("breast_cancer x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### aggregation data #####
set.seed(1)
k=7 # number of clusters
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("aggregation x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### asymmetric data #####
set.seed(1)
k=5 # number of clusters
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("asymmetric x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### R15 data #####
set.seed(1)
k=15 # number of clusters
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("R15 x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))



##### D31 data #####
set.seed(1)
k=31 # number of clusters
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# Spectral clustering with the Gaussian (rbfdot) kernel on the data itself;
# specc chooses the kernel width (kpar = "automatic")
clusters = specc(as.matrix(X), centers = k)
spectrum_label = clusters@.Data

# the kernel width chosen by specc: k(x, y) = exp(-sigma_used * ||x - y||^2)
sigma_used = kernlab::kpar(kernlab::kernelf(clusters))$sigma

# number of clusters:
n_clusters = length(unique(spectrum_label))

# ARI
ari = adjustedRandIndex(spectrum_label, Y)

# the average silhouette index
sil = mean(silhouette(spectrum_label, dist(X))[,3])

print(paste("D31 x spectrum clustering: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters, "specc sigma", sigma_used))