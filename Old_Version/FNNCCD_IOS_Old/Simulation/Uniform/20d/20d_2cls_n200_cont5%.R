#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/Outlyingness_Score/FNNCCD_IOS.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_20d_999%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/FNNCCD_IOS/Simulation/Uniform/Threshold.R")
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
# method: ascend or descend order when finding the radius
n = 200
d = 20
cont = 0.05
method="descend"

if(d==2){quant=0.85;threshold=Threshold[1]# the level of K-test
} else if(d==3) {
  quant=0.9;threshold=Threshold[2]
} else if(d==5) {
  quant=0.95;threshold=Threshold[3]
} else if(d==10) {
  quant=0.99;threshold=Threshold[4]
} else if(d==20){
  quant=0.999;threshold=Threshold[5]
} else if(d==50){
  quant=0.999;threshold=Threshold[6]
} else if(d==100){
  quant=0.999;threshold=Threshold[7]
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
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% FNNCCD_IOS(datax=x,simul=simul,method=method,d=d)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x=x,scores=IOS_scores,threshold=threshold, n=n, n0=n0)

mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("IOS: the mean TPR is", mean[1],",","and, the mean TNR", mean[2]))

t2 = Sys.time()
t2-t1

save.image("/mmfs1/home/rzs0112/code_working_folder/FNNCCD_IOS/Simulation/Uniform/20d/20d_2cls_n200_cont5%.RData")