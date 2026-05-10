#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
load("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile/NN-test-simul_20d_999%.RData")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)
ratio = c(1.000000, 1.024079, 1.068926, 1.163903, 1.291596)

d = 20
iteN = 1000
quant = 0.999
cores = detectCores()
# cores = 24 # 13900K



# simulation settings
kappa1 = 6
mu1 = 36.4*ratio[5]
expand1 = 0
r = 0.1
kappa2 = 0
scale = 0.005
mu2 = 36.4*ratio[5]
expand2 = 0
slen = 1
kappa_O = 20

# main function for outlier function
# datax: the data set
# simul: the simulated windows or envelope 
nnccd.rkccd.outlier = function(datax, simul){

  set.seed(321)
  graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method="decend", dom.method="greedy2", 
                                     simul=simul, niter=100)
  NN.result = nnccd.silhouette(graph,datax,cls=NULL,ind=NULL, lenDlimit=Inf)
  NN.result = c(graph,NN.result)
  
  # a list of clusters of the dataset
  nnccd.clusters = lapply(1:NN.result$si.ind, function(x){return(datax[NN.result$label==x,])})
  
  #the cover info of each Dominated graph
  catch.num = NN.result$catch
  ddatax = as.matrix(dist(datax))
  catch.info = lapply(1:NN.result$si.ind, function(x){return(datax[ddatax[NN.result$Int.D[x],]<=NN.result$Int.R[x],])})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(catch.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:NN.result$si.ind, function(x){
    return(ksccd.connected(nnccd.clusters[[x]],delta.max[x],sequential=FALSE,alpha=0.05)$member)})
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
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% nnccd.rkccd.outlier(x, simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count1(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

# outlier.result[[1]]$clusters
# x = 10
# len = length(outlier.result[[x]]$clusters)
# for(i in 1:len){
#   p = 2*i
#   plot(x=outlier.result[[x]]$cluster[[i]][,1],y=outlier.result[[x]]$cluster[[i]][,2], asp=1,ylim=c(0,1),xlim=c(0,1), pch = p)
# }


tt2 = Sys.time()
tt2-tt1