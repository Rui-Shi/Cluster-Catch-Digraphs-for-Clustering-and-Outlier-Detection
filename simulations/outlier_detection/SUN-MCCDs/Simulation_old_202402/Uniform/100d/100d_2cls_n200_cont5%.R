#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_100d_999%.RData")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

set.seed(123)
t1 = Sys.time()
cores = detectCores()
# cores = 24 # for 13900K


# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiments

n = 200
d = 100
cont = 0.05

if(d==2){quant=0.85 # the level of K-test
} else if(d==3) {
  quant=0.90
} else if(d==5) {
  quant=0.95
} else if(d==10) {
  quant=0.999
} else {
  quant=0.999
}

M = 10000 # the number of simulated K-function values
iteN = 1000 # the number of simulated data set.
cls_dis = 3 # the distances between each cluster center
otl_dis = 2 # the minimal distances of outliers to cluster centers

# the min and max of the radii of clusters
r_min = 0.7
r_max = 1.3


# Begin
# Simulate the lower tail quantile of NN distance
set.seed(321)
# simul = NNDest.simpois.lower.quant(n, d, quant=quant, niter=M, shape="sphere")


MFNNCCD.outlier = function(datax,simul){
  set.seed(321)
  # find the MCG of NNCCD
  Mgraph = nnccd.clustering.mutual.connected(datax, low.num=3, quantile="lower", 
                                             method="decend", dom.method="greedy2", simul=simul, niter=100)
  
  # find the optimal number of clusters by maximizing silhouette index
  result.cls = nnccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,k=20, min.cls = max(0.05,cont/2))
  
  #the connected points of each components
  connect.info = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(Mgraph$D.member==Mgraph$member[x])],])
  })
  
  #find the clustering info
  nnccd.clusters = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(result.cls$label==Mgraph$member[x])],])
  })
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(connect.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:result.cls$sil.ind, function(x){
    return(ksccd.connected(nnccd.clusters[[x]],delta.max[x],sequential=FALSE)$member)})
  return(c(clusters=list(nnccd.clusters),label=list(cluster.component)))
}



# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)*0.5)
n2 = round(n*(1-cont)*0.5)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    if(r1 > otl_dis & r2 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFNNCCD.outlier(x,simul)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

t2 = Sys.time()
t2-t1