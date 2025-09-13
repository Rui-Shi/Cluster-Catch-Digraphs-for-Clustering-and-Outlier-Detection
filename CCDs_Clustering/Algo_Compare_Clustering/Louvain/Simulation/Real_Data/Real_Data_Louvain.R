source("/media/rui/exNVME/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)
library(igraph)

sigma = 2

##### iris data #####
n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("iris x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####
n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("seeds x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### user_knowledge data #####
n_col = ncol(user_knowledge)
X = user_knowledge[, -n_col] # remove the class column for clustering
Y = user_knowledge[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("user_knowledge x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####
n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("wholesale x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer data #####
n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("breast_cancer x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### aggregation data #####
n_col = ncol(aggregation)
X = aggregation[, -n_col] # remove the class column for clustering
Y = aggregation[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("aggregation x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### asymmetric data #####
n_col = ncol(asymmetric)
X = asymmetric[, -n_col] # remove the class column for clustering
Y = asymmetric[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("asymmetric x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### R15 data #####
n_col = ncol(R15)
X = R15[, -n_col] # remove the class column for clustering
Y = R15[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("R15 x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### D31 data #####
n_col = ncol(D31)
X = D31[, -n_col] # remove the class column for clustering
Y = D31[, n_col] # true cluster labels

# distance matrix
distM = as.matrix(dist(X))

# Apply Gaussian kernel
similarity_matrix = exp(-distM^2 / (2 * sigma^2))

# created a weighted graph
graph = graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = T)

# clustering with louvain
louvain_clustering = cluster_louvain(graph)

louvain_label = membership(louvain_clustering)

# number of clusters:
n_clusters = length(unique(louvain_label))

# ARI
ari = adjustedRandIndex(louvain_label, Y)

# the average silhouette index
sil = mean(silhouette(louvain_label, dist(X))[,3])

print(paste("D31 x louvain: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))