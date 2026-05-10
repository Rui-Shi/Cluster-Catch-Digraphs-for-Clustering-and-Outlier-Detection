#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_2d_85%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/ratio1.R")


library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 2
iteN = 1000
method = "ascend"
cont = 0.13
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
kappa1 = 6
mu1 = ratio[1]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = ratio[1]
expand2 = 0
slen = 1
kappa_O = 20


MFNNCCD.outlier = function(datax,simul){
  # find the MCG of NNCCD
  Mgraph = nnccd.clustering.mutual.connected(datax, low.num=3, quantile="lower", 
                                             method=method, dom.method="greedy2", simul=simul, niter=100)
  
  # find the optimal number of clusters by maximizing silhouette index
  result.cls = nnccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,k=20, min.cls = max(0.05,cont/2))
  
  #the connected points of each components
  connect.info = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(Mgraph$D.member==Mgraph$member[x])],])
  })
  
  #find the clustering info
  nnccd.clusters = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(result.cls$label==Mgraph$member[x])],])
  })
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(connect.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:result.cls$sil.ind, function(x){
    return(ksccd.connected(nnccd.clusters[[x]],delta.max[x],sequential=FALSE)$member)})
  return(c(clusters=list(nnccd.clusters),label=list(cluster.component)))
}


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
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFNNCCD.outlier(x, simul)
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

save.image("/mmfs1/home/rzs0112/code_working_folder/M-FNNCCDs/Simulation/Complex_Clusters/Matern/2d_clx_cls.RData")
tt2 = Sys.time()
tt2-tt1