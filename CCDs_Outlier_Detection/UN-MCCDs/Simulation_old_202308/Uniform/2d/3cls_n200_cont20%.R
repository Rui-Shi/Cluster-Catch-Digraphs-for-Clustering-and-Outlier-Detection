#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
library(parallel)
library(doParallel)
library(MASS)


# main function for outlier function
# datax: the data set
# simul: the simulated windows or envelope 
nnccd.rkccd.outlier = function(datax, simul){
  set.seed(321)
  graph <- nnccd_clustering_quantile(datax, low.num=3, quantile="lower", 
                                     method="accend", dom.method="greedy2", 
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


# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiment
n = 200
d = 2
cont = 0.20
iteN = 1000
quant=0.80
cores = detectCores()
# cores = 24 # 13900K

# Simulate the lower tail quantile of NN distance
simul = NNDest.simpois.lower.quant(n, d, quant=quant, niter=1000, shape="sphere")
print("3cls_n200_cont20%")

# simulate two clusters of equal size within two unit balls centered at (3,3) and (6,3)
mu1 = rep(3,d)
mu2 = c(6,rep(3,d-1))
mu3 = c(3,6,rep(3,d-2))
mu = rbind(mu1,mu2,mu3)
mu = apply(mu, 2, mean)

n1 = round(n*(1-cont)*1/3)
n2 = round(n*(1-cont)*1/3)
n3 = round(n*(1-cont)*1/3)+1
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > 1.5 & r2 > 1.5 & r3 > 1.5){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})

t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% nnccd.rkccd.outlier(x, simul=simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

# count the success rate and false positive rate
count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# try different density for the two clusters: 1/6, 2/6, 3/6
n1 = round(n*(1-cont)*1/6)
n2 = round(n*(1-cont)*2/6)
n3 = round(n*(1-cont)*3/6)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*5 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    if(r1 > 1.5 & r2 > 1.5 & r3 > 1.5){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,outlier))
})


t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list, .packages = c("MASS","cluster")) %dopar% nnccd.rkccd.outlier(x, simul=simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

# count the success rate and false positive rate
count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))