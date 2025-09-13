source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(dbscan)
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)

k = 4 # fix the minpts to 5 (4 neighors)

##### iris data #####
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 0.8
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
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 0.9
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
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 1.5
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
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 1
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
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 2
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
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 1.2
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
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 0.15
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
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 0.35
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
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

kNNdistplot(X, k) # set minPts = 5, 4 neighorhood distance

eps = 0.6
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