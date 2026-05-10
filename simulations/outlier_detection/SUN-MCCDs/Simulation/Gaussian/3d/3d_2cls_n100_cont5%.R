#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/M-FNNCCDs/M-FNNCCDs.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_3d_90%.RData")
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

n = 100
d = 3
cont = 0.05

method = "ascend"
if(d<=5){
  min.cls = max(0.04, cont/2)
} else {
  min.cls=0
}

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


# simulate two clusters of equal size within two unit balls centered at (3,3) and (3+cls_dis,3)
# the radius of clusters are random numbers between 0.7-1.3
mu1 = rep(3,d)
mu2 = c(3+cls_dis,rep(3,d-1))
mu = rbind(mu1,mu2)
mu = apply(mu,2,mean)

n1 = round(n*(1-cont)*0.5)
n2 = round(n*(1-cont)*0.5)-1
n0 = round(n*cont)

# The level of noise (the number of points of a Gaussian cluster outside the above radius)
noise_level = 0.01
sigma = 1/sqrt(qchisq(1-noise_level,d))

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1,mu1,diag(d)*(sigma*runif(1,r_min,r_max))^2)
  data2 = mvrnorm(n2,mu2,diag(d)*(sigma*runif(1,r_min,r_max))^2)
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
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFNNCCD_outlier(x, simul, min.cls=min.cls, method=method)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean True positive rate is", 1-mean[2]))

t2 = Sys.time()
t2-t1

save.image("/mmfs1/home/rzs0112/code_working_folder/M-FNNCCDs/Simulation/Gaussian/3d/3d_2cls_n100_cont5%.RData")