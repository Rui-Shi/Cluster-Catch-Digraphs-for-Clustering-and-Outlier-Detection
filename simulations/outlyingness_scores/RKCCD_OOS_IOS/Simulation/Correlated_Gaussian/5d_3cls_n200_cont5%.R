t1 = Sys.time()
source("/media/rui/exNVME/code_working_folder/ccds/RK_CCD_New.R")
source("/media/rui/exNVME/code_working_folder/Outlyingness_Score/RKCCD_OOS_IOS.R")
source("/media/rui/exNVME/code_working_folder/general functions/count.R")
load("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile/RK-test-simul_5d_99%.RData")
source("/media/rui/exNVME/code_working_folder/RKCCD_OOS_IOS/Simulation/Gaussian/Threshold.R")
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
# method: ascend or descend order when finding the radius
n = 200
d = 5
cont = 0.05

if(d==2){quant=0.99;threshold_OOS=Threshold_OOS[1];threshold_IOS=Threshold_IOS[1]# the level of K-test
} else if(d==3) {
  quant=0.99;threshold_OOS=Threshold_OOS[2];threshold_IOS=Threshold_IOS[2]
} else if(d==5) {
  quant=0.99;threshold_OOS=Threshold_OOS[3];threshold_IOS=Threshold_IOS[3]
} else if(d==10) {
  quant=0.999;threshold_OOS=Threshold_OOS[4];threshold_IOS=Threshold_IOS[4]
} else if(d==20){
  quant=0.999;threshold_OOS=Threshold_OOS[5];threshold_IOS=Threshold_IOS[5]
} else if(d==50){
  quant=0.999;threshold_OOS=Threshold_OOS[6];threshold_IOS=Threshold_IOS[6]
} else if(d==100){
  quant=0.999;threshold_OOS=Threshold_OOS[7];threshold_IOS=Threshold_IOS[7]
}

M = 10000 # the number of simulated K-function values
iteN = 1000 # the number of simulated data set.
cls_dis = 3 # the distances between each cluster center
otl_dis = 2 # the minimal distances of outliers to cluster centers

# the min and max of the radii of clusters
r_min = 0.7
r_max = 1.3

# Begin
# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu3 = c(3+cls_dis, 3+cls_dis, rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)/3)
n2 = round(n*(1-cont)/3)
n3 = round(n*(1-cont)/3)+1
n0 = round(n*cont)


sigma = 0.25

# corrlation = 0
rho = 0

# covariance matrix
covM1 = diag(1,d)
covM2 = diag(1,d)
covM3 = diag(1,d)


set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  data2 = mvrnorm(n2, mu2, covM2*(sigma*runif(1,r_min,r_max))^2)
  data3 = mvrnorm(n3, mu3, covM3*(sigma*runif(1,r_min,r_max))^2)
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
  return(rbind(data1, data2, data3, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
OOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_OOS(datax=x,simul=simul,d=d,quant=quant)
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_IOS(datax=x,simul=simul,d=d,quant=quant)
stopCluster(cl)

count.result_OOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=OOS_scores,threshold=threshold_OOS, n=n, n0=n0)
mean_OOS = c(mean(count.result_OOS[,1]), mean(count.result_OOS[,2]), mean(count.result_OOS[,3]), mean(count.result_OOS[,4]))
print(paste("OOS: the mean TPR is", mean_OOS[1],",","the mean TNR", mean_OOS[2], ",", "the mean BA", mean_OOS[3], "the mean F2", mean_OOS[4]))

count.result_IOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=IOS_scores,threshold=threshold_IOS, n=n, n0=n0)
mean_IOS = c(mean(count.result_IOS[,1]), mean(count.result_IOS[,2]), mean(count.result_IOS[,3]), mean(count.result_IOS[,4]))
print(paste("IOS: the mean TPR is", mean_IOS[1],",","the mean TNR", mean_IOS[2], ",", "the mean BA", mean_IOS[3], "the mean F2", mean_IOS[4]))


# corrlation = 0.1
rho = 0.1

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

covM2 = diag(1,d)
covM2[lower.tri(covM2)] = -(rho-0.1)
covM2[upper.tri(covM2)] = -(rho-0.1)

covM3 = diag(1,d)
covM3[lower.tri(covM3)] = rho+0.1
covM3[upper.tri(covM3)] = rho+0.1


set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  data2 = mvrnorm(n2, mu2, covM2*(sigma*runif(1,r_min,r_max))^2)
  data3 = mvrnorm(n3, mu3, covM3*(sigma*runif(1,r_min,r_max))^2)
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
  return(rbind(data1, data2, data3, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
OOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_OOS(datax=x,simul=simul,d=d,quant=quant)
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_IOS(datax=x,simul=simul,d=d,quant=quant)
stopCluster(cl)

count.result_OOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=OOS_scores,threshold=threshold_OOS, n=n, n0=n0)
mean_OOS = c(mean(count.result_OOS[,1]), mean(count.result_OOS[,2]), mean(count.result_OOS[,3]), mean(count.result_OOS[,4]))
print(paste("OOS: the mean TPR is", mean_OOS[1],",","the mean TNR", mean_OOS[2], ",", "the mean BA", mean_OOS[3], "the mean F2", mean_OOS[4]))

count.result_IOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=IOS_scores,threshold=threshold_IOS, n=n, n0=n0)
mean_IOS = c(mean(count.result_IOS[,1]), mean(count.result_IOS[,2]), mean(count.result_IOS[,3]), mean(count.result_IOS[,4]))
print(paste("IOS: the mean TPR is", mean_IOS[1],",","the mean TNR", mean_IOS[2], ",", "the mean BA", mean_IOS[3], "the mean F2", mean_IOS[4]))


# corrlation = 0.2
rho = 0.2

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

covM2 = diag(1,d)
covM2[lower.tri(covM2)] = -(rho-0.1)
covM2[upper.tri(covM2)] = -(rho-0.1)

covM3 = diag(1,d)
covM3[lower.tri(covM3)] = rho+0.1
covM3[upper.tri(covM3)] = rho+0.1


set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  data2 = mvrnorm(n2, mu2, covM2*(sigma*runif(1,r_min,r_max))^2)
  data3 = mvrnorm(n3, mu3, covM3*(sigma*runif(1,r_min,r_max))^2)
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
  return(rbind(data1, data2, data3, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
OOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_OOS(datax=x,simul=simul,d=d,quant=quant)
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_IOS(datax=x,simul=simul,d=d,quant=quant)
stopCluster(cl)

count.result_OOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=OOS_scores,threshold=threshold_OOS, n=n, n0=n0)
mean_OOS = c(mean(count.result_OOS[,1]), mean(count.result_OOS[,2]), mean(count.result_OOS[,3]), mean(count.result_OOS[,4]))
print(paste("OOS: the mean TPR is", mean_OOS[1],",","the mean TNR", mean_OOS[2], ",", "the mean BA", mean_OOS[3], "the mean F2", mean_OOS[4]))

count.result_IOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=IOS_scores,threshold=threshold_IOS, n=n, n0=n0)
mean_IOS = c(mean(count.result_IOS[,1]), mean(count.result_IOS[,2]), mean(count.result_IOS[,3]), mean(count.result_IOS[,4]))
print(paste("IOS: the mean TPR is", mean_IOS[1],",","the mean TNR", mean_IOS[2], ",", "the mean BA", mean_IOS[3], "the mean F2", mean_IOS[4]))


# corrlation = 0.5
rho = 0.5

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

covM2 = diag(1,d)
covM2[lower.tri(covM2)] = rho-0.1
covM2[upper.tri(covM2)] = rho-0.1

covM3 = diag(1,d)
covM3[lower.tri(covM3)] = rho+0.1
covM3[upper.tri(covM3)] = rho+0.1


set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  data2 = mvrnorm(n2, mu2, covM2*(sigma*runif(1,r_min,r_max))^2)
  data3 = mvrnorm(n3, mu3, covM3*(sigma*runif(1,r_min,r_max))^2)
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
  return(rbind(data1, data2, data3, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
OOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_OOS(datax=x,simul=simul,d=d,quant=quant)
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_IOS(datax=x,simul=simul,d=d,quant=quant)
stopCluster(cl)

count.result_OOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=OOS_scores,threshold=threshold_OOS, n=n, n0=n0)
mean_OOS = c(mean(count.result_OOS[,1]), mean(count.result_OOS[,2]), mean(count.result_OOS[,3]), mean(count.result_OOS[,4]))
print(paste("OOS: the mean TPR is", mean_OOS[1],",","the mean TNR", mean_OOS[2], ",", "the mean BA", mean_OOS[3], "the mean F2", mean_OOS[4]))

count.result_IOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=IOS_scores,threshold=threshold_IOS, n=n, n0=n0)
mean_IOS = c(mean(count.result_IOS[,1]), mean(count.result_IOS[,2]), mean(count.result_IOS[,3]), mean(count.result_IOS[,4]))
print(paste("IOS: the mean TPR is", mean_IOS[1],",","the mean TNR", mean_IOS[2], ",", "the mean BA", mean_IOS[3], "the mean F2", mean_IOS[4]))


t2 = Sys.time()
t2-t1

