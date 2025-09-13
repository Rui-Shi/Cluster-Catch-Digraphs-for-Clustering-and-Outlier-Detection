#!/usr/bin/env Rscript
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

# d: dimension
# cont: contamination level
# n0: number of outliers
# n: size of the dataset
# n1,n2,...: the size of each clusters
# iteN: number of experiment
n = 40
d = 2
cont = 0.10
iteN = 1000
quant=0.99
cores = detectCores()
# cores = 24 # 13900K
print("2d_4cls_n40_cont10%")


# Simulate the lower tail quantile of NN distance
set.seed(123)
simul = Kest.simpois.edge.quantile(n, d, rn=10, quan = quant, niter=1000)

MCCD.outlier = function(datax,simul){
  set.seed(321)
  RK.result = RKCCD_correct_quant(datax,r.seq=10, dom.method="greedy2", 
                                  quan=0.99, simul=simul, cls=NULL,ind=NULL, niter=100)
  # a list of clusters of the dataset
  rccd.clusters = lapply(1:RK.result$si.ind, function(x){return(datax[RK.result$label==x,])})
  
  #the cover info of each Dominated graph
  catch.num = RK.result$catch
  ddatax = as.matrix(dist(datax))
  catch.info= lapply(1:RK.result$si.ind, function(x){return(datax[ddatax[RK.result$Int.D[x],]<=RK.result$Int.R[x],])})
  
  #find the max density delta such that the KS-CCD for each Dominated graph is connected
  delta.max = sapply(catch.info, connected.ksccd.m)
  
  # find the connected components for each cluster with the delta obtained above.
  cluster.component = lapply(1:RK.result$si.ind, function(x){
    return(ksccd.connected(rccd.clusters[[x]],delta.max[x],sequential=FALSE,alpha=0.05)$member)})
  return(c(clusters=list(rccd.clusters),label=list(cluster.component)))
}

# simulate two clusters of equal size within two unit balls centered at (3,3) and (6,3)
mu1 = c(3,3)
mu2 = c(6,3)
mu3 = c(3,6)
mu4 = c(6,6)
mu = rbind(mu1,mu2,mu3,mu4)
mu = apply(mu, 2, mean)

n1 = round(n*(1-cont)*1/4)
n2 = round(n*(1-cont)*1/4)
n3 = round(n*(1-cont)*1/4)
n4 = round(n*(1-cont)*1/4)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  data4 = rpoisball.unit(n4,d) + matrix(rep(mu4,n4),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*6 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    r4 = sqrt(sum((temp-mu4)^2))
    if(r1 > 1.5 & r2 > 1.5 & r3 > 1.5 & r4>1.5){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,data4,outlier))
})


t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MCCD.outlier(x,simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# try different density for the 4 clusters: 2/10, 3/10, 3/10, 2/10
n1 = round(n*(1-cont)*2/10)
n2 = round(n*(1-cont)*3/10)
n3 = round(n*(1-cont)*3/10)
n4 = round(n*(1-cont)*2/10)
n0 = round(n*cont)

set.seed(123)
data.list = lapply(1:iteN, function(x){
  data1 = rpoisball.unit(n1,d) + matrix(rep(mu1,n1),ncol=d,byrow=T)
  data2 = rpoisball.unit(n2,d) + matrix(rep(mu2,n2),ncol=d,byrow=T)
  data3 = rpoisball.unit(n3,d) + matrix(rep(mu3,n3),ncol=d,byrow=T)
  data4 = rpoisball.unit(n4,d) + matrix(rep(mu4,n4),ncol=d,byrow=T)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*6 + mu
    r1 = sqrt(sum((temp-mu1)^2))
    r2 = sqrt(sum((temp-mu2)^2))
    r3 = sqrt(sum((temp-mu3)^2))
    r4 = sqrt(sum((temp-mu4)^2))
    if(r1 > 1.5 & r2 > 1.5 & r3 > 1.5 & r4>1.5){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data1,data2,data3,data4,outlier))
})

t1 = Sys.time()
# cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)
outlier.result = foreach(x=data.list,.packages = c("MASS","cluster","igraph")) %dopar% MCCD.outlier(x,simul)
stopCluster(cl)
t2 = Sys.time()
t2-t1

count.result = foreach(x=1:iteN,.combine = rbind) %do% count(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))


# #draw plots for report
# library(plotrix)
# plot(datax, pch=16,asp=1)
# points(nnccd.clusters[[2]],pch=16,col="blue")
# points(nnccd.clusters[[1]],pch=16,col="red")
# for(i in 1:n){
#   draw.circle(datax[Mgraph$D[i],1],datax[Mgraph$D[i],2],Mgraph$R[i])
# }