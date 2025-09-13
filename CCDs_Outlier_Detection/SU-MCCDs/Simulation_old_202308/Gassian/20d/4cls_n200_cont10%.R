#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiment
n = 200
d = 20
cont = 0.10
iteN = 1000
quant=0.99
cores = detectCores()
# cores = 24 # 13900K
print("5d_4cls_n200_cont10%")

# Simulate the lower tail quantile of NN distance
set.seed(123)
simul = Kest.simpois.edge.quantile(n, d, rn=10, quan = quant, niter=1000)

MFCCD.outlier = function(datax,simul){
  set.seed(321)
  Mgraph = rccd.clustering.mutual.connected_correct(datax, low.num=2, r.seq=10, dom.method="greedy2", quan=quant, simul=simul, niter=1000)
  
  # find the optimal number of clusters by maximizing silhouette index
  result.cls = rccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,k=5, min.cls = max(0.05,cont/2))
  
  #the connected points of each components
  connect.info = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(Mgraph$D.member==Mgraph$member[x])],])
  })
  
  #find the clustering info
  rccd.clusters = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(result.cls$label==Mgraph$member[x])],])
  })
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(connect.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:result.cls$sil.ind, function(x){
    return(ksccd.connected(rccd.clusters[[x]],delta.max[x],sequential=FALSE)$member)})
  return(c(clusters=list(rccd.clusters),label=list(cluster.component)))
}

# simulate two clusters of equal size within two unit balls centered at (3,3) and (6,3)
mu1 = rep(3,d)
mu2 = c(9,rep(3,d-1))
mu3 = c(3,9,rep(3,d-2))
mu4 = c(3,3,9,rep(3,d-3))
mu = rbind(mu1,mu2,mu3,mu4)
mu = apply(mu, 2, mean)

n1 = round(n*(1-cont)*1/4)
n2 = round(n*(1-cont)*1/4)
n3 = round(n*(1-cont)*1/4)
n4 = round(n*(1-cont)*1/4)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1,mu1,diag(d))
  data2 = mvrnorm(n2,mu2,diag(d))
  data3 = mvrnorm(n3,mu3,diag(d))
  data4 = mvrnorm(n4,mu4,diag(d))
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*12 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    r4 = sqrt(sum((temp-mu4)^2))
    if(r1 > 4 & r2 > 4 & r3 > 4 & r4 > 4){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,data4,outlier))
})


t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFCCD.outlier(x,simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# try different density for the 4 clusters: 2/10, 3/10, 3/10, 2/10
n1 = round(n*(1-cont)*2/10)
n2 = round(n*(1-cont)*3/10)
n3 = round(n*(1-cont)*3/10)
n4 = round(n*(1-cont)*2/10)
n0 = round(n*cont)

set.seed(321)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1,mu1,diag(d))
  data2 = mvrnorm(n2,mu2,diag(d))
  data3 = mvrnorm(n3,mu3,diag(d))
  data4 = mvrnorm(n4,mu4,diag(d))
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*12 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    r4 = sqrt(sum((temp-mu4)^2))
    if(r1 > 4 & r2 > 4 & r3 > 4 & r4 > 4){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,data4,outlier))
})

t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFCCD.outlier(x,simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# #draw plots for report
# library(plotrix)
# plot(datax, pch=16,asp=1)
# points(nnccd.clusters[[2]],pch=16,col="blue")
# points(nnccd.clusters[[1]],pch=16,col="red")
# for(i in 1:n){
#   draw.circle(datax[Mgraph$D[i],1],datax[Mgraph$D[i],2],Mgraph$R[i])
# }