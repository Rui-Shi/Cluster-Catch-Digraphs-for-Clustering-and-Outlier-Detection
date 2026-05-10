library(parallel)
library(plotrix)
library(dplyr)
library(MASS)
library(doParallel)
source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccd_ks_NEW.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/ccdfunctions.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/Kest.R")
source("/mmfs1/home/rzs0112/code_working_folder/KS-MCGs/simulations/Uniform_MeanSlippage/svdd_sphere.R")

t1 = Sys.time()

cores = detectCores()

######################### datasets with outliers ###############################
#####
## synthetic datasets that contain 6% outliers(3 outliers out of 50)
# n: number of points in each data set.
# n0: number of outliers
# d: dimensionality
# data.otl: datasets that contains outliers
# len1: number of datasets
# len2: number of datasets to simulate
# sft: the distance to shift
# cont: comtamination level

#####
##shift by 3
set.seed(123)

d = 20

n = 100
cont = 0.06 #contaminate levels
n1 = round(n*(1-cont))
n0 = n-n1

len1 = 1000
len2 = 200
sft = 3


C = 0.05 # penalty parameter for svdd

data.list = lapply(1:len1, function(x){
  data = rpoisball.unit(n,d)
  # outliers
  data[-c(1:n1),1] = data[-c(1:n1),1] + sft
  return(data)
})

cl = makeCluster(cores)

# the KS-MCG algorithm
KS_MCG = function(x){
  svdd_res = svdd_sphere(x, C) # perform svdd with sphere boundary
  radius = as.numeric(svdd_res$radius)
  center = as.numeric(svdd_res$center)
  
  data.simul = NULL
  for(i in 1:len2){
    data = rpoisball.unit(n,d)*radius
    data = t(t(data)+ center)
    data.simul = c(data.simul,list(data))
  }
  m.max.seq = sapply(data.simul, connected.ksccd.m)
  m.max.0.05 = quantile(m.max.seq,0.05) # the 5% quantile of the 1000 generated m #
  res = ksccd.connected(x,m.max.0.05,sequential=FALSE)
  return(res)
}

registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% KS_MCG(x)
stopCluster(cl)

count.result = foreach(x=1:len1,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# contamination level = 10%
cont = 0.1 #contaminate levels
n1 = round(n*(1-cont))
n0 = n-n1

data.list = lapply(1:len1, function(x){
  data = rpoisball.unit(n,d)
  # outliers
  data[-c(1:n1),1] = data[-c(1:n1),1] + sft
  return(data)
})

cl = makeCluster(cores)

registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% KS_MCG(x)
stopCluster(cl)

count.result = foreach(x=1:len1,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))



# contamination level = 20%
cont = 0.2 #contaminate levels
n1 = round(n*(1-cont))
n0 = n-n1

data.list = lapply(1:len1, function(x){
  data = rpoisball.unit(n,d)
  # outliers
  data[-c(1:n1),1] = data[-c(1:n1),1] + sft
  return(data)
})

cl = makeCluster(cores)

registerDoParallel(cl)
components_result = foreach(x=data.list,.packages = c("MASS","cluster")) %dopar% KS_MCG(x)
stopCluster(cl)

count.result = foreach(x=1:len1,.combine = rbind) %do% count2(x)
mean = c(mean(count.result[,1]),mean(count.result[,2]))
print(paste("The mean success rate is", mean[1],",","and, the mean false positive rate is", mean[2]))

t2 = Sys.time()
t2-t1