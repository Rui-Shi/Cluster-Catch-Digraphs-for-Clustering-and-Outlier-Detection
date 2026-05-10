#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/RK-test_quantile/RK-test-simul_5d_99%.RData")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/ratio3.R")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 5
iteN = 1000
quant = 0.99
cont = 0.13
cores = detectCores()
# cores = 24 # 13900K



# simulation settings
kappa1 = 3
mu1 = ratio[3]
expand1 = 0
r = 0.1
kappa2 = 3
scale = 0.005
mu2 = ratio[3]
expand2 = 0
slen = 1
kappa_O = 20

# Simulate the lower tail quantile of NN distance
# set.seed(123)
# simul = Kest.simpois.edge.quantile(n, d, rn=10, quan = quant, niter=1000)

MFCCD.outlier = function(datax,simul){
  Mgraph = rccd.clustering.mutual.connected_correct(datax, low.num=2, r.seq=10, 
                                                    dom.method="greedy2", quan=quant, simul=simul, niter=1000)
  
  # find the optimal number of clusters by maximizing silhouette index
  result.cls = rccd.silhouette_mutual(Mgraph,datax,ind=NULL, lenClimit=Inf,
                                      k=max(d,5), min.cls = max(0.05,cont/2))
  
  #the connected points of each components
  connect.info = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(Mgraph$D.member==Mgraph$member[x])],])
  })
  
  #find the clustering info
  rccd.clusters = lapply(1:result.cls$sil.ind, function(x){
    return(datax[Mgraph$D[which(result.cls$label==Mgraph$member[x])],])
  })
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(connect.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:result.cls$sil.ind, function(x){
    return(ksccd.connected(rccd.clusters[[x]],delta.max[x],sequential=FALSE)$member)})
  return(c(clusters=list(rccd.clusters),label=list(cluster.component)))
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
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MFCCD.outlier(x, simul)
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

save.image("/mmfs1/home/rzs0112/code_working_folder/M-FCCDs/Simulation/Complex_Clusters/Mix/5d_clx_cls.RData")
tt2 = Sys.time()
tt2-tt1