tt1 = Sys.time()
source("G:/code_working_folder/general functions/Uni-Gau_cls.R")
source("G:/code_working_folder/general functions/ratio2.R")
source("G:/code_working_folder/Algo_Compare/DBSCAN/DBSCAN.R")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)
library(dbscan)

d = 20
iteN = 1000
cores = detectCores()
cores = 24 # 13900K

# simulation settings
kappa1 = 0
mu1 = ratio[5]
expand1 = 0
r = 0.1
kappa2 = 6
scale = 0.005
mu2 = ratio[5]
expand2 = 0
slen = 1
kappa_O = 20

MinPts = 4 # The minimum number of points required within the eps radius to form a dense region.
cont = 0.09 # the average contamination level

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

cl <- makeCluster(cores)
registerDoParallel(cl)
labels_list = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% DBSCAN(x, MinPts, cont)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_DBSCAN2(data.num[[x]][1], data.num[[x]][2], labels_list[[x]])
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("DBSCAN: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))

tt2 = Sys.time()
tt2 - tt1