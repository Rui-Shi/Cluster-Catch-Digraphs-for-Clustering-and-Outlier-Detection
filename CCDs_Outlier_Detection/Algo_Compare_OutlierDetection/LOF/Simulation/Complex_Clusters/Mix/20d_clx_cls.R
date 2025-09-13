tt1 = Sys.time()
library(parallel)
library(doParallel)
library(MASS)
library(igraph)
source("G:/code_working_folder/Algo_Compare/LOF/LOF.R")
source("G:/code_working_folder/general functions/count.R")
source("G:/code_working_folder/general functions/ratio3.R")
source("G:/code_working_folder/general functions/Uni-Gau_cls.R")

set.seed(123)
t1 = Sys.time()
cores = detectCores()
cores = 24 # for 13900K

LB = 11 # Lower bound for MinPts
UB = 30 # Upper bound for MinPts
Thresh = 1.5 # threshhold for outliers

iteN = 1000 # the number of simulated data set.

d=20
kappa1 = 3
mu1 = ratio[5]
expand1 = 0
r = 0.1
kappa2 = 3
scale = 0.005
mu2 = ratio[5]
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

cl <- makeCluster(cores)
registerDoParallel(cl)
scores = foreach(x=data.list,.packages = c("MASS","cluster","dbscan")) %dopar% LOF(x,LB,UB)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores1(x=x,scores=scores,threshold=Thresh)
mean = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]),mean(count.result[,4]))

print(paste("LOF: the mean TPR is", mean[1],",","the mean TNR", mean[2],",","the mean BA", mean[3],",","the mean F2", mean[4]))


tt2 = Sys.time()
tt2 - tt1