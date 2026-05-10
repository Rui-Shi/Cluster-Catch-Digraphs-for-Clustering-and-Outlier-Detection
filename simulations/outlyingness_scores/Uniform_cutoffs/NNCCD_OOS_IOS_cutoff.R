source("/media/rui/exNVME/code_working_folder/Outlyingness_Score/NNCCD_OOS_IOS.R")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)
library(cluster)
t1 = Sys.time()
cutoff_OOS = function(x, method="ascend", simul,d,n,n0){
  scores = sort(NNCCD_OOS(datax = x, method=method, simul=simul, d=d))
  return(scores[n-n0-3])}

cutoff_IOS = function(x, method="ascend", simul, d, n, n0){
  scores = sort(NNCCD_IOS(datax = x, method=method, simul=simul, d=d))
  return(scores[n-n0-3])}

set.seed(123)
cores = detectCores()
cores = 24 # for 13900K


# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiments
n = 200
d = 2
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

n1 = round(n*(1-cont)*0.5)
n2 = round(n*(1-cont)*0.5)
n0 = round(n*cont)

##d=2##
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_2d_85%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_2d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_2d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_2d))
median(unlist(NNCCD_OOS_2d))



##d=3##
d=3
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_3d_90%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_3d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_3d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_3d))
median(unlist(NNCCD_OOS_3d))



##d=5##
d=5
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_5d_95%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_5d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_5d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_5d))
median(unlist(NNCCD_OOS_5d))



##d=10##
d=10
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_10d_99%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_10d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_10d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_10d))
median(unlist(NNCCD_OOS_10d))




##d=20##
d=20
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_20d_999%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_20d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_20d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_20d))
median(unlist(NNCCD_OOS_20d))



##d=50##
d=50
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_50d_999%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_50d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_50d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

median(unlist(NNCCD_IOS_50d))
median(unlist(NNCCD_OOS_50d))

##d=100##
d=100
load("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile/NN-test-simul_100d_999%.RData")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

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
NNCCD_OOS_100d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_OOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)

cl <- makeCluster(cores)
registerDoParallel(cl)
NNCCD_IOS_100d = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% cutoff_IOS(x,method="descend",simul=simul,d=d,n=n,n0=n0)
stopCluster(cl)


median(unlist(NNCCD_IOS_100d))
median(unlist(NNCCD_OOS_100d))

save.image("/media/rui/exNVME/code_working_folder/Outlyingness_Score/NNCCD_OSS_IOS_cutoff.RData")
t2 = Sys.time()
t2-t1