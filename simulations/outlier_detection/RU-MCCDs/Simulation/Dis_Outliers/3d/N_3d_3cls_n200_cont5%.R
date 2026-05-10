#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/RK-test_quantile/RK-test-simul_3d_99%.RData")
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
d = 3
cont = 0.05

if(d<10){quant=0.99 # the level of K-test
} else {
  quant=0.999
}

M = 10000 # the number of simulated K-function values
iteN = 1000 # the number of simulated data set.
cls_dis = 3 # the distances between each cluster center
otl_dis = 1.25 # the minimal distances of outliers to cluster centers

# the min and max of the radii of clusters
r_min = 0.7
r_max = 1.3


# Begin
# Simulate the lower tail quantile of K-function
# simul = Kest.simpois.edge.quantile(n, d, rn=10, quan = quant, niter=M)

## Main outlier detection function
ksccd.rkccd.outlier = function(datax){
  RK.result = RKCCD_correct_quant(datax,r.seq=10, dom.method="greedy2", 
                                  quan=quant, simul=simul, cls=NULL,ind=NULL, niter=M)
  # a list of clusters of the dataset
  rccd.clusters = lapply(1:RK.result$si.ind, function(x){return(datax[RK.result$label==x,])})
  
  #the cover info of each Dominated graph
  catch.num = RK.result$catch
  ddatax = as.matrix(dist(datax))
  catch.info= lapply(1:RK.result$si.ind, function(x){return(datax[ddatax[RK.result$Int.D[x],]<=RK.result$Int.R[x],])})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(catch.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:RK.result$si.ind, function(x){
    return(ksccd.connected(rccd.clusters[[x]],delta.max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(rccd.clusters),label=list(cluster.component)))
}


mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3,3+cls_dis,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d)*runif(1,r_min,r_max) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > otl_dis & r2 > otl_dis & r3 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.rkccd.outlier(x)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# distance of outliers to cluster centers = 1.75
otl_dis = 1.5

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3,3+cls_dis,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d)*runif(1,r_min,r_max) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > otl_dis & r2 > otl_dis & r3 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.rkccd.outlier(x)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# distance of outliers to cluster centers = 2
otl_dis = 1.75

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3,3+cls_dis,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d)*runif(1,r_min,r_max) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > otl_dis & r2 > otl_dis & r3 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.rkccd.outlier(x)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# distance of outliers to cluster centers = 2.25
otl_dis = 2

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3,3+cls_dis,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d)*runif(1,r_min,r_max) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > otl_dis & r2 > otl_dis & r3 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.rkccd.outlier(x)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# distance of outliers to cluster centers = 2.5
otl_dis = 2.25

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3,3+cls_dis,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d)*runif(1,r_min,r_max) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > otl_dis & r2 > otl_dis & r3 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.rkccd.outlier(x)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

t2 = Sys.time()
t2-t1