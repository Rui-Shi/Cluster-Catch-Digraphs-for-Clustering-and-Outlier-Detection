t1 = Sys.time()

source("/media/rui/exNVME/code_working_folder/CCDs_Clustering/KS-CCDs.R")
source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(mclust)
library(parallel)
library(doParallel)
library(MASS)
library(cluster)
library(igraph)

cores = detectCores()
cores=20 # 14700K
cl = makeCluster(cores)
registerDoParallel(cl)

delta_range = c(seq(0.01, 1, 0.01), seq(1.1, 10, 0.1), seq(11, 20, 1))

##### iris data #####
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("iris x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("seeds x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("user_knowledge x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("wholesale x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("breast_cancer x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("aggregation x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("asymmetric x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### R15 data #####
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("R15 x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

Sil_seq = foreach(delta=delta_range,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(X, m=delta)$si
Sil_seq = unlist(Sil_seq)

# identify the best density parameter delta
delta_best = delta_range[which.max(Sil_seq)]
print(delta_best)

# clustering with the best delta
KSCCD_label = KSCCD_clustering(X, m=delta_best)$label

# number of clusters:
n_clusters = length(unique(KSCCD_label))

# ARI
ari = adjustedRandIndex(KSCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(KSCCD_label, dist(X))[,3])

print(paste("D31 x KS-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



stopCluster(cl)

t2 = Sys.time()
t2-t1