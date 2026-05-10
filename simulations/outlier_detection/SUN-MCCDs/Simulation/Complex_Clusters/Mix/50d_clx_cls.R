#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/M-FNNCCDs/M-FNNCCDs.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_50d_999%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/ratio3.R")


library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 50
iteN = 1000
cont = 0.10
min.cls = 0
method = "descend"
cores = detectCores()
# cores = 24 # 13900K

if(d==2){quant=0.85 # the level of K-test
} else if(d==3) {
  quant=0.90
} else if(d==5) {
  quant=0.95
} else if(d==10) {
  quant=0.99
} else {
  quant=0.999
}


# simulation settings
kappa1 = 3
mu1 = ratio[6]
expand1 = 0
r = 0.1
kappa2 = 3
scale = 0.005
mu2 = ratio[6]
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


t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFNNCCD_outlier(x, simul, min.cls=min.cls, method=method)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count1(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]), mean(count.result[,3]))
print(paste("The mean success rate is", mean[1],",","the mean true negative rate is", 1-mean[2], ",", "and the mean F2-score is", mean[3]))

# outlier.result[[1]]$clusters
# x = 10
# len = length(outlier.result[[x]]$clusters)
# for(i in 1:len){
#   p = 2*i
#   plot(x=outlier.result[[x]]$cluster[[i]][,1],y=outlier.result[[x]]$cluster[[i]][,2], asp=1,ylim=c(0,1),xlim=c(0,1), pch = p)
# }

save.image("/mmfs1/home/rzs0112/code_working_folder/M-FNNCCDs/Simulation/Complex_Clusters/Mix/50d_clx_cls.RData")
tt2 = Sys.time()
tt2-tt1