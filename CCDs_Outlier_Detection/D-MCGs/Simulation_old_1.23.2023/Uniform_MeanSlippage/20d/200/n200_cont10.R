library(parallel)
library(plotrix)
library(dplyr)
library(MASS)
library(doParallel)
source("/media/shirui001/WD SN350/code_working folder/ccds/ccd_ks_NEW.R")
source("/media/shirui001/WD SN350/code_working folder/ccds/functions.R")
source("/media/shirui001/WD SN350/code_working folder/Heuristic 1/functions/mKNN_CCD_functions.R")
source("/media/shirui001/WD SN350/code_working folder/general functions/count.R")

t1 = Sys.time()

set.seed(1)
len = 1000
n = 200
d = 20
cont = 0.10 #contaminate levels
i = 1
range = NULL
while(i <= d){
  range.temp = c(min=0,max=1)
  range = rbind(range,range.temp)
  i = i + 1
}
rownames(range) = c(1:d)

data.list = NULL
for(i in 1:len){
  data = apply(range,1,function(t){
    return(runif(n,t[1],t[2]))
  })
  colnames(data) = c(1:d)
  data.list = c(data.list,list(data))
}

#####
#parallel computation
cl = makeCluster(24)

m.max.seq = parSapply(cl, data.list, connected.ksccd.m)
stopCluster(cl)

m.max.0.05 = quantile(m.max.seq,0.05) # the 5% quantile of the 1000 generated m ##
print(m.max.0.05)

######################### datasets with outliers ###############################
#####
## synthetic datasets that contain 6% outliers(3 outliers out of 50)
# n: number of points in each data set.
# n0: number of outliers (the last 6 objects are outliers)
# d: dimensionality
# min and max: ranges in each dimension (normal data)
# data.otl: datasets that contains outliers
# len: number of datasets
# sft: the distance to shift (5,10,15,20)

#####
##shift by 15
set.seed(2)
len = 100
n1 = round(n*(1-cont))
n0 = n-n1
sft = c(2,3,4)
data.list = NULL
for(i in 1:len){
  normal = apply(range,1,function(t){
    return(runif(n1,t[1],t[2]))
  })
  
  range.outlier = range
  range.outlier[1,] = range.outlier[1,] + sft[1]
  outlier = apply(range.outlier,1,function(t){
    return(runif(n0,t[1],t[2]))   # mean shift by 15
  })
  data.otl = rbind(normal,outlier)
  colnames(data.otl) = c(1:d)
  #plot(data.otl)
  data.list = c(data.list,list(data.otl))
}

#parallel computation to find the m for the data set that contains outliers
cores = detectCores()
cl = makeCluster(24)
m.max.otl.seq = parLapply(cl, data.list, connected.ksccd.m)
stopCluster(cl)
m.max.otl.seq = unlist(m.max.otl.seq)
length(which(m.max.otl.seq<m.max.0.05))/len # the percentage of the data.otl where mutiple components (outliers or clusters) could be detected.

cl <- makeCluster(24)
registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.connected(x,m.max.0.05,sequential=FALSE)
stopCluster(cl)

count.result = foreach(x=1:len,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))






#####
##shift by 20
set.seed(3)
data.list = NULL
for(i in 1:len){
  normal = apply(range,1,function(t){
    return(runif(n1,t[1],t[2]))
  })
  
  range.outlier = range
  range.outlier[1,] = range.outlier[1,] + sft[2]
  outlier = apply(range.outlier,1,function(t){
    return(runif(n0,t[1],t[2]))   # mean shift by 20
  })
  data.otl = rbind(normal,outlier)
  colnames(data.otl) = c(1:d)
  #plot(data.otl)
  data.list = c(data.list,list(data.otl))
}

#parallel computation to find the m for the data set that contains outliers
cores = detectCores()
cl = makeCluster(24)
m.max.otl.seq = parLapply(cl, data.list, connected.ksccd.m)
stopCluster(cl)
m.max.otl.seq = unlist(m.max.otl.seq)
length(which(m.max.otl.seq<m.max.0.05))/len # the percentage of the data.otl where mutiple components (outliers or clusters) could be detected.

cl <- makeCluster(24)
registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.connected(x,m.max.0.05,sequential=FALSE)
stopCluster(cl)

count.result = foreach(x=1:len,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))




#####
##shift by 25
set.seed(4)
data.list = NULL
for(i in 1:len){
  normal = apply(range,1,function(t){
    return(runif(n1,t[1],t[2]))
  })
  
  range.outlier = range
  range.outlier[1,] = range.outlier[1,] + sft[3]
  outlier = apply(range.outlier,1,function(t){
    return(runif(n0,t[1],t[2]))   # mean shift by 25
  })
  data.otl = rbind(normal,outlier)
  colnames(data.otl) = c(1:d)
  #plot(data.otl)
  data.list = c(data.list,list(data.otl))
}

#parallel computation to find the m for the data set that contains outliers
cores = detectCores()
cl = makeCluster(24)
m.max.otl.seq = parLapply(cl, data.list, connected.ksccd.m)
stopCluster(cl)
m.max.otl.seq = unlist(m.max.otl.seq)
length(which(m.max.otl.seq<m.max.0.05))/len # the percentage of the data.otl where mutiple components (outliers or clusters) could be detected.

cl <- makeCluster(24)
registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% ksccd.connected(x,m.max.0.05,sequential=FALSE)
stopCluster(cl)

count.result = foreach(x=1:len,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

t2 = Sys.time()
t2-t1