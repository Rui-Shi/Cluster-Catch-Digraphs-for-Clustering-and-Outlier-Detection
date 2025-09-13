#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/Outlyingness_Score/FNNCCD_IOS.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_3d_90%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/FNNCCD_IOS/Simulation/Uniform/Threshold.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/ratio3.R")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

t1 = Sys.time()
cores = detectCores()
# cores = 24 # for 13900K

d = 3
cont = 0.1
method="ascend"

# simulation settings
kappa1 = 3
mu1 = ratio[2]
expand1 = 0
r = 0.1
kappa2 = 3
scale = 0.005
mu2 = ratio[2]
expand2 = 0
slen = 1
kappa_O = 20

if(d==2){quant=0.99;threshold=Threshold[1]# the level of K-test
} else if(d==3) {
  quant=0.99;threshold=Threshold[2]
} else if(d==5) {
  quant=0.99;threshold=Threshold[3]
} else if(d==10) {
  quant=0.999;threshold=Threshold[4]
} else if(d==20){
  quant=0.999;threshold=Threshold[5]
} else if(d==50){
  quant=0.999;threshold=Threshold[6]
} else if(d==100){
  quant=0.999;threshold=Threshold[7]
}

M = 10000 # the number of simulated K-function values
iteN = 1000 # the number of simulated data set.

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
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% FNNCCD_IOS(datax=x,simul=simul,method=method,d=d)
stopCluster(cl)

count.result = foreach(x=1:iteN,.combine = rbind) %do% count_scores1(x=x,scores=IOS_scores,threshold=threshold)
mean_IOS = c(mean(count.result[,1]),mean(count.result[,2]),mean(count.result[,3]))
print(paste("IOS: the mean TPR is", mean_IOS[1], "and, the mean TNR", mean_IOS[2], "and, the mean F2-score is", mean_IOS[3]))

save.image("/mmfs1/home/rzs0112/code_working_folder/FNNCCD_IOS/Simulation/Complex_Clusters/Mix/3d_clx_cls.RData")

t2 = Sys.time()
t2-t1