source("G:/code_working_folder/Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R")
source("G:/code_working_folder/Algo_Compare_Clustering/MST/MST_Clustering.R")
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)
library(foreach)
library(doParallel)

cores = detectCores()
cl = makeCluster(cores)
registerDoParallel(cl)

thresh_list=seq(1.01, 2, 0.01)  # the range for threshes: the ratio to cut a edge when comparing its adjacent edges

##### iris data #####

n_col = ncol(iris)
X = iris[, -n_col] # remove the class column for clustering
Y = iris[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]
print(thresh_op)

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("iris x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### seeds data #####

n_col = ncol(seeds)
X = seeds[, -n_col] # remove the class column for clustering
Y = seeds[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]
print(thresh_op)

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("seeds x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### glass data #####

n_col = ncol(glass)
X = glass[, -n_col] # remove the class column for clustering
Y = glass[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("glass x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wholesale data #####

n_col = ncol(wholesale)
X = wholesale[, -n_col] # remove the class column for clustering
Y = wholesale[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("wholesale x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### breast_cancer #####

n_col = ncol(breast_cancer)
X = breast_cancer[, -n_col] # remove the class column for clustering
Y = breast_cancer[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("breast_cancer x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### transfusion #####

n_col = ncol(transfusion)
X = transfusion[, -n_col] # remove the class column for clustering
Y = transfusion[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("transfusion x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### phishing #####

n_col = ncol(phishing)
X = phishing[, -n_col] # remove the class column for clustering
Y = phishing[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  if(length(unique(MST_label))==1) sil=0
  else sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("phishing x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wine_red #####

n_col = ncol(wine_red)
X = wine_red[, -n_col] # remove the class column for clustering
Y = wine_red[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  if(length(unique(MST_label))==1) sil=0
  else sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("wine_red x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### yeast #####

n_col = ncol(yeast)
X = yeast[, -n_col] # remove the class column for clustering
Y = yeast[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  if(length(unique(MST_label))==1) sil=0
  else sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("yeast x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



##### wine_white #####

n_col = ncol(wine_white)
X = wine_white[, -n_col] # remove the class column for clustering
Y = wine_white[, n_col] # true cluster labels

# search for the best cutoff value by maximizing the silhouette index
sil_list = foreach(i=thresh_list,.packages = c("MASS","cluster","igraph")) %dopar% {
  MST_label= MST_Clustering(X, i)
  # the average silhouette index
  if(length(unique(MST_label))==1) sil=0
  else sil = mean(silhouette(MST_label, dist(X))[,3])
  return(sil)
}

sil_array = unlist(sil_list)

thresh_op = thresh_list[which.max(sil_array)]

# construct MST for clustering with the best cutoff value
MST_label= MST_Clustering(X, thresh_op)
print(thresh_op)

# number of clusters:
n_clusters = length(unique(MST_label))

# ARI
ari = adjustedRandIndex(MST_label, Y)

# the average silhouette index
sil = mean(silhouette(MST_label, dist(X))[,3])

print(paste("wine_white x MST: the ARI is", ari, "the average silhouette index is", sil, "number of clusters detected", n_clusters))



stopCluster(cl)