#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/CCDs_Clustering/KS-CCDs.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
setwd("/mmfs1/home/rzs0112/code_working_folder/CCDs_Clustering/KS_CCDs/Simulations/Uniform/5d")
library(mclust)
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

set.seed(2024)
t1 = Sys.time()
cores = detectCores()

n = 50
d = 5

delta_range = c(seq(0.01, 1, 0.01), seq(1.1, 10, 0.1), seq(11, 20, 1))

iteN = 500 # the number of simulated data set.
cls_dis = 3 # the distances between each cluster center

# the min and max of the radii of clusters
r_min = 0.8
r_max = 1.2



##### 2 clusters #####
# the radius of clusters are random numbers between 0.7-1.3
# N: number of clusters
# n1, n2...: the number of each clusters
# label_T: the ground truth labels
N = 2
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
MU = rbind(mu1,mu2)
mu = apply(MU,2,mean)
n1 = round(n/N)
n2 = round(n/N)
n_cls=c(n1, n2)
label_T = c(rep(0, n1), rep(1, n2))

# generate data sets
data_list = lapply(1:iteN, uniform_cls, n_cls=n_cls, 
                                       d=d, 
                                       centers=MU,
                                       noise=F,
                                       rad_range=c(r_min, r_max))

cl <- makeCluster(cores)
registerDoParallel(cl)
# find the best density parameter delta
Sil_seq_2cls=c()
for(delta in delta_range){
  clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta)
  Sil_seq_2cls = c(Sil_seq_2cls, mean(sapply(clustering_result, function(x) x$si)))
}

# identify the best density parameter delta
delta_best_2cls = delta_range[which.max(Sil_seq_2cls)]

clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta_best_2cls)


# store the ARI and the number of times identifying correct number of clusters
ARI_2cls = sapply(clustering_result, function(x) adjustedRandIndex(x$label, label_T))
num_clusters_detected_2cls = sapply(clustering_result, function(x) x$si.ind)
sil_2cls = sapply(clustering_result, function(x) x$si)

print(paste("KS-CCDs: The mean ARI is", mean(ARI_2cls),",",
            "The mean average silhouette index is", mean(sil_2cls),",",
            "and, the percentatage of times identifying the correct number of clusters is", mean(num_clusters_detected_2cls==N)*100,"%,", 
            "best delta:", delta_best_2cls))



##### 3 clusters #####
# the radius of clusters are random numbers between 0.7-1.3
# N: number of clusters
# n1, n2...: the number of each clusters
# label_T: the ground truth labels
N = 3
mu1 = rep(3,d)
mu2 = c(3+cls_dis, rep(3,d-1))
mu3 = c(3+cls_dis, 3+cls_dis, rep(3,d-2))
MU = rbind(mu1, mu2, mu3)
mu = apply(MU,2,mean)
n1 = round(n/N)-1
n2 = round(n/N)
n3 = round(n/N)
n_cls=c(n1, n2, n3)
label_T = c(rep(0, n1), rep(1, n2), rep(2, n3))

# generate data sets
data_list = lapply(1:iteN, uniform_cls, n_cls=n_cls, 
                   d=d, 
                   centers=MU,
                   noise=F,
                   rad_range=c(r_min, r_max))

# find the best density parameter delta
Sil_seq_3cls=c()
for(delta in delta_range){
  clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta)
  Sil_seq_3cls = c(Sil_seq_3cls, mean(sapply(clustering_result, function(x) x$si)))
}

# identify the best density parameter delta
delta_best_3cls = delta_range[which.max(Sil_seq_3cls)]

clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta_best_3cls)


# store the ARI and the number of times identifying correct number of clusters
ARI_3cls = sapply(clustering_result, function(x) adjustedRandIndex(x$label, label_T))
num_clusters_detected_3cls = sapply(clustering_result, function(x) x$si.ind)
sil_3cls = sapply(clustering_result, function(x) x$si)

print(paste("KS-CCDs: The mean ARI is", mean(ARI_3cls),",",
            "The mean average silhouette index is", mean(sil_3cls),",",
            "and, the percentatage of times identifying the correct number of clusters is", mean(num_clusters_detected_3cls==N)*100,"%,", 
            "best delta:", delta_best_3cls))



##### 5 clusters #####
# the radius of clusters are random numbers between 0.7-1.3
# N: number of clusters
# n1, n2...: the number of each clusters
# label_T: the ground truth labels
N = 5
mu1 = rep(3,d)
mu2 = c(3+cls_dis, rep(3,d-1))
mu3 = c(3+cls_dis, 3+cls_dis, rep(3,d-2))
mu4 = c(3, 3+cls_dis-0.5, rep(3,d-2))
mu5 = c(3+2*cls_dis-0.5, rep(3,d-1))
MU = rbind(mu1, mu2, mu3, mu4, mu5)
mu = apply(MU,2,mean)
n1 = round(n/N)
n2 = round(n/N)
n3 = round(n/N)
n4 = round(n/N)
n5 = round(n/N)
n_cls=c(n1, n2, n3, n4, n5)
label_T = c(rep(0, n1), rep(1, n2), rep(2, n3), rep(3, n4), rep(4, n5))

# generate data sets
data_list = lapply(1:iteN, uniform_cls, n_cls=n_cls, 
                   d=d, 
                   centers=MU,
                   noise=F,
                   rad_range=c(r_min, r_max))

# find the best density parameter delta
Sil_seq_5cls=c()
for(delta in delta_range){
  clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta)
  Sil_seq_5cls = c(Sil_seq_5cls, mean(sapply(clustering_result, function(x) x$si)))
}

# identify the best density parameter delta
delta_best_5cls = delta_range[which.max(Sil_seq_5cls)]

clustering_result = foreach(x=data_list,.packages = c("MASS","cluster")) %dopar% KSCCD_clustering(datax=x, m=delta_best_5cls)

# store the ARI and the number of times identifying correct number of clusters
ARI_5cls = sapply(clustering_result, function(x) adjustedRandIndex(x$label, label_T))
num_clusters_detected_5cls = sapply(clustering_result, function(x) x$si.ind)
sil_5cls = sapply(clustering_result, function(x) x$si)

print(paste("KS-CCDs: The mean ARI is", mean(ARI_5cls),",",
            "The mean average silhouette index is", mean(sil_5cls),",",
            "and, the percentatage of times identifying the correct number of clusters is", mean(num_clusters_detected_5cls==N)*100,"%,", 
            "best delta:", delta_best_5cls))

t2 = Sys.time()
t2-t1

stopCluster(cl)
save.image("5d_235clusters_n50.RData")