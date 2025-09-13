library(parallel)
library(doParallel)
library(MASS)
library(igraph)
source("G:/code_working_folder/Algo_Compare/DBSCAN/DBSCAN.R")
source("G:/code_working_folder/general functions/Strauss.R")

set.seed(123)
t1 = Sys.time()
cores = detectCores()
# cores = 24 # for 13900K

MinPts = 4 # The minimum number of points required within the eps radius to form a dense region.

# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiments

n = 100
d = 10
cont = 0.05

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
labels_list = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% DBSCAN(x, MinPts, cont)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_DBSCAN(data.list[[x]], cont, labels_list[[x]])
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("DBSCAN: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))

t2 = Sys.time()
t2-t1

save.image("G:/code_working_folder/Algo_Compare/DBSCAN/Simulation/Uniform/10d/10d_2cls_n100_cont5%.RData")