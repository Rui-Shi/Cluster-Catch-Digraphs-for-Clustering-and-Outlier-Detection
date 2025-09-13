#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/Outlyingness_Score/RKCCD_OOS_IOS.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/RK-test_quantile/RK-test-simul_100d_999%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/RKCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/Threshold.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/ratio1.R")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

t1 = Sys.time()
cores = detectCores()
# cores = 24 # for 13900K

d = 100
min.cls=0.04

# simulation settings
kappa1 = 6
mu1 = ratio[7]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = ratio[7]
expand2 = 0
slen = 1
kappa_O = 20

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
OOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_OOS(datax=x,simul=simul,d=d,quant=quant)
IOS_scores = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% RKCCD_IOS(datax=x,simul=simul,d=d,quant=quant,min.cls=min.cls)
stopCluster(cl)

count.result_OOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores1(x=x,scores=OOS_scores,threshold=threshold_OOS)
mean_OOS = c(mean(count.result_OOS[,1]),mean(count.result_OOS[,2]),mean(count.result_OOS[,3]), mean(count.result_OOS[,4]))
print(paste("OOS: the mean TPR is", mean_OOS[1], "and, the mean TNR", mean_OOS[2], "and, BA is", mean_OOS[3], "and, the mean F2-score is", mean_OOS[4]))

count.result_IOS = foreach(x=1:iteN,.combine = rbind) %do% count_scores1(x=x,scores=IOS_scores,threshold=threshold_IOS)
mean_IOS = c(mean(count.result_IOS[,1]),mean(count.result_IOS[,2]),mean(count.result_IOS[,3]), mean(count.result_IOS[,4]))
print(paste("IOS: the mean TPR is", mean_IOS[1], "and, the mean TNR", mean_IOS[2], "and, BA is", mean_IOS[3], "and, the mean F2-score is", mean_IOS[4]))


t2 = Sys.time()
t2-t1

save.image("/mmfs1/home/rzs0112/code_working_folder/RKCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/100d_clx_cls.RData")