source("/media/rui/exNVME/code_working_folder/ccds/RK_CCD_New.R")
source("/media/rui/exNVME/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/media/rui/exNVME/code_working_folder/general functions/count.R")
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_10d_999%.RData")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

set.seed(123)
t1 = Sys.time()
cores = detectCores()
cores = 24 # for 13900K

# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiments

n = 200
d = 10

if(d<10){quant=0.99 # the level of K-test
} else {
  quant=0.999
}

iteN = 1000 # the number of simulated data set.
cls_dis = 6 # the distances between each cluster center
cont = 0.05 # % of outliers = 0.05

# the min and max of the radii of clusters
r_min = 0.7
r_max = 1.3

# Begin
# Simulate the lower tail quantile of K-function
# simul = Kest.simpois.edge.quantile(n, d, rn=10, quan = quant, niter=M)

## Main outlier detection function
MFCCD.outlier = function(datax,simul){
  set.seed(321)
  Mgraph = rccd.clustering.mutual.connected_correct(datax, low.num=2, r.seq=10, 
                                                    dom.method="greedy2", quan=quant, simul=simul, niter=1000)
  
  # find the optimal number of clusters by maximizing silhouette index
  result.cls = rccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,
                                      k=max(d,5), min.cls = max(0.05,cont/2))
  
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



# the distances between cluster centers and outlier center
otl_dis = 1.5 

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = mu1 + c(otl_dis,rep(0,d-1))

n1 = round(n*(1-cont)/2)
n2 = round(n*(1-cont)/2)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  
  outlier = t(sapply(1:n0,function(x){
    outlier_temp = rpoisball.unit(1,d) + mu
  }))
  
  rownames(outlier) = NULL
  return(rbind(data1,data2,outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% MFCCD.outlier(x, simul)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# the distances between cluster centers and outlier center
otl_dis = 2 

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = mu1 + c(otl_dis,rep(0,d-1))

n1 = round(n*(1-cont)/2)
n2 = round(n*(1-cont)/2)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  
  outlier = t(sapply(1:n0,function(x){
    outlier_temp = rpoisball.unit(1,d) + mu
  }))
  
  rownames(outlier) = NULL
  return(rbind(data1,data2,outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% MFCCD.outlier(x, simul)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# the distances between cluster centers and outlier center
otl_dis = 2.5 

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = mu1 + c(otl_dis,rep(0,d-1))

n1 = round(n*(1-cont)/2)
n2 = round(n*(1-cont)/2)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  
  outlier = t(sapply(1:n0,function(x){
    outlier_temp = rpoisball.unit(1,d) + mu
  }))
  
  rownames(outlier) = NULL
  return(rbind(data1,data2,outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% MFCCD.outlier(x, simul)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# the distances between cluster centers and outlier center
otl_dis = 3

mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = mu1 + c(otl_dis,rep(0,d-1))

n1 = round(n*(1-cont)/2)
n2 = round(n*(1-cont)/2)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d)*runif(1,r_min,r_max) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d)*runif(1,r_min,r_max) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  
  outlier = t(sapply(1:n0,function(x){
    outlier_temp = rpoisball.unit(1,d) + mu
  }))
  
  rownames(outlier) = NULL
  return(rbind(data1,data2,outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% MFCCD.outlier(x, simul)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

t2 = Sys.time()
t2-t1