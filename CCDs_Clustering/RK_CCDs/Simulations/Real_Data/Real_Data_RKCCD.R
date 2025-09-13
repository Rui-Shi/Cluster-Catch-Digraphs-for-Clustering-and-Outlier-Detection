t1 = Sys.time()

source("/media/rui/exNVME/code_working_folder/CCDs_Clustering/RK-CCDs.R")
source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(mclust)
library(parallel)
library(doParallel)
library(MASS)
library(cluster)
library(igraph)

##### iris data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_4d_99%.RData")
quant = 0.99

n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("iris x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_7d_99%.RData")
quant = 0.99

n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("seeds x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_6d_99%.RData")
quant = 0.99

n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("user_knowledge x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_8d_99%.RData")
quant = 0.99

n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("wholesale x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_9d_999%.RData")
quant = 0.999

n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("breast_cancer x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_2d_99%.RData")
quant = 0.99

n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("aggregation x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_2d_99%.RData")
quant = 0.99

n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("asymmetric x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### R15 data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_2d_99%.RData")
quant = 0.99

n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("R15 x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_2d_99%.RData")
quant = 0.99

n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# clustering with RK-CCD
RKCCD_label = RKCCD_clustering(datax=X, simul=simul, quant=quant)$label

# number of clusters:
n_clusters = length(unique(RKCCD_label))

# ARI
ari = adjustedRandIndex(RKCCD_label, Y)

# the average silhouette index
sil = mean(silhouette(RKCCD_label, dist(X))[,3])

print(paste("D31 x RK-CCDs: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))


t2 = Sys.time()
t2-t1