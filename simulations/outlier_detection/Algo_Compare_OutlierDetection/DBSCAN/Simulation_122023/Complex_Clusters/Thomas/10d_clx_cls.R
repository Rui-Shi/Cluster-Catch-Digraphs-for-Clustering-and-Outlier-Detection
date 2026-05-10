tt1 = Sys.time()
source("G:/code_working_folder/general functions/Uni-Gau_cls.R")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)
library(dbscan)
ratio = c(1.000000, 1.064383, 1.193162, 1.604715, 2.859147)

d = 10
k = 4 # the value of MinPTS for DBSCAN
quant = 0.09 # noise level
iteN = 1000
cores = detectCores()
# cores = 24 # 13900K

# simulation settings
kappa1 = 0
mu1 = 37.1*ratio[1]
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = 37.1*ratio[1]
expand2 = 0
slen = 1
kappa_O = 20


# simulate clusters of random sizes and positions
set.seed(1234)
data.listNum = lapply(1:iteN, function(x){
  data_simu =  Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu$Matérn_children
  cls2 = data_simu$Thomas_children
  cls3 = data_simu$noise
  outlier = data_simu$Outlier
  return(c(data = list(rbind(cls1, cls2, cls3, outlier)), num = list(data_simu$num)))
})

data.list = lapply(1:iteN, function(x){
  return(data.listNum[[x]]$data)
})

data.num = lapply(1:iteN, function(x){
  return(data.listNum[[x]]$num)
})

perf = sapply(1:iteN, function(x){
  datax = data.list[[x]]
  dist_M = as.matrix(dist(datax))
  
  k_dist = apply(dist_M, 1, function(row){  # compute the k-distance for each point
    return(sort(row)[k+1])
  })
  
  eps = quantile(k_dist, 1-quant) # find the eps of DBSCAN based on noise level
  
  result = dbscan(datax, eps, k) # conduct DBSCAN
  
  clusters = result$cluster  # count SR and FP
  FP = sum(clusters[1:data.num[[x]][1]]==0)/data.num[[x]][1]
  SR = sum(clusters[-c(1:data.num[[x]][1])]==0)/data.num[[x]][2]
  return(c(SR,FP))
})

result_mean = apply(perf, 1, mean)

print(result_mean)

tt2 = Sys.time()
tt2 - tt1