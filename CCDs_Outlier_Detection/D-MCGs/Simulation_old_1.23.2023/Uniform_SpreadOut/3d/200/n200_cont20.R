library(parallel)
library(plotrix)
library(dplyr)
library(MASS)
library(doParallel)
source("/media/shirui001/WD SN350/code_working folder/ccds/ccd_ks_NEW.R")
source("/media/shirui001/WD SN350/code_working folder/ccds/functions.R")
source("/media/shirui001/WD SN350/code_working folder/Heuristic 1/functions/mKNN_CCD_functions.R")
source("/media/shirui001/WD SN350/code_working folder/general functions/count.R")
source("/media/shirui001/WD SN350/code_working folder/Algorithm 3/functions/RK_CCD_New.R")

t1 = Sys.time()

set.seed(1)
len = 1000
n = 200
d = 3
cont = 0.20 #contaminate levels

data.list = NULL
for(i in 1:len){
  data = rpoisball.unit(n,d)
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
## synthetic datasets that contain outliers
# n: number of points in each data set.
# n0: number of outliers (the last 6 objects are outliers)
# d: dimensionality
# data.otl: datasets that contains outliers
# len: number of datasets


set.seed(2)
len = 100
n1 = round(n*(1-cont))
n0 = n-n1

data.list = lapply(1:len, function(x){
  data = rpoisball.unit(n1,d)
  i = 0
  outlier = NULL
  while(i < n0){
    temp = rpoisball.unit(1,d)*3
    r = sqrt(temp[1]^2+temp[2]^2)
    if(r > 1.5){
      outlier = rbind(outlier,temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  return(rbind(data,outlier))
})

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