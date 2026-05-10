t1 = Sys.time()
library(parallel)
library(doParallel)
library(MASS)
library(igraph)
source("G:/code_working_folder/Algo_Compare/LOF/LOF.R")
source("G:/code_working_folder/general functions/count.R")
source("G:/code_working_folder/general functions/Strauss.R")

set.seed(123)
t1 = Sys.time()
cores = detectCores()
cores = 20 # for 14700

LB = 11 # Lower bound for MinPts
UB = 30 # Upper bound for MinPts
Thresh = 1.5 # threshhold for outliers

# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiments
# method: ascend or descend order when finding the radius
n = 200
d = 3
cont = 0.05

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
mu = mu1

n1 = round(n*(1-cont))
n0 = round(n*cont)


sigma = 0.35

# correlation = 0
rho = 0

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*4 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    if(r1 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
scores = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% LOF(x,LB,UB)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x, scores, Thresh, n, n0)
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("LOF: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))

# correlation = 0.1
rho = 0.1

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*4 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    if(r1 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
scores = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% LOF(x,LB,UB)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x, scores, Thresh, n, n0)
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("LOF: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))


# corrlation = 0.2
rho = 0.2

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*4 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    if(r1 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
scores = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% LOF(x,LB,UB)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x, scores, Thresh, n, n0)
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("LOF: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))

# corrlation = 0.5
rho = 0.5

# covariance matrix
covM1 = diag(1,d)
covM1[lower.tri(covM1)] = rho
covM1[upper.tri(covM1)] = rho

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = mvrnorm(n1, mu1, covM1*(sigma*runif(1,r_min,r_max))^2)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*4 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    if(r1 > otl_dis){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1, outlier))
})

cl <- makeCluster(cores)
registerDoParallel(cl)
scores = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% LOF(x,LB,UB)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores(x, scores, Thresh, n, n0)
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("LOF: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))


t2 = Sys.time()
t2-t1