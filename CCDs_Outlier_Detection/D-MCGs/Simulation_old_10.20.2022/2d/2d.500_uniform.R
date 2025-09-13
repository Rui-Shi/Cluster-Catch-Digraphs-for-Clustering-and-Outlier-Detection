library(parallel)

###code for KS-CCDs###
source("G:/OneDrive - Auburn University/Research Outliers Detection/code/ccds/ccd_ks_NEW.R")
source("G:/OneDrive - Auburn University/Research Outliers Detection/code/ccds/functions.R")

#####
## simulate n uniform data
# n: number of points in each data set.
# d: dimensionality
# min and max: ranges in each dimension
# len: number of data set

len = 1000
n = 500
d = 2
range1 = c(min=0,max=10)
range2 = c(min=0,max=10)
range = rbind(range1,range2)
data.list = NULL
for(i in 1:len){
  data = apply(range,1,function(t){
    return(runif(n,t[1],t[2]))
  })
  colnames(data) = c(1:d)
  #plot(data)
  data.list = c(data.list,list(data))
}


source("G:/OneDrive - Auburn University/Research Outliers Detection/code/KS-CCDs & mKNN/functions/mKNN_CCD_functions.R")


#####
#parallel computation
cores = detectCores()
cl = makeCluster(20)
t1 = Sys.time()
m.max.seq = parLapply(cl, data.list, connected.ksccd.m)
stopCluster(cl)
m.max.seq = unlist(m.max.seq)
t2 = Sys.time()
t2-t1

m.max.0.05 = quantile(m.max.seq,0.05) # the 5% quantile of the 1000 generated m ##
print(m.max.0.05)

mean(m.max.seq)

mean(m.max.seq)
den = density(m.max.seq)
plot(den,main="the empirical distribution of m (d=2, n=500, uniform)")



#####
## simulate n uniform data and n0 outliers
# n: number of points in each data set.
# n0: number of outliers
# d: dimensionality
# min and max: ranges in each dimension (normal data)
# data.otl: data set that contains outliers
# len: number of data set
len = 100
n = 450
n0 = 50
d = 2
range1 = c(min=0,max=10)
range2 = c(min=0,max=10)
range = rbind(range1,range2)
data.otl.list = NULL
for(i in 1:len){
  normal = apply(range,1,function(t){
    return(runif(n,t[1],t[2]))
  })
  
  outlier = apply(range,1,function(t){
    return(runif(n0,t[1]+12,t[2]+12))
  })
  data.otl = rbind(normal,outlier)
  colnames(data.otl) = c(1:d)
  #plot(data.otl)
  data.otl.list = c(data.otl.list,list(data.otl))
}



#####
#parallel computation to find the m for the data set that contains outliers
cores = detectCores()
cl = makeCluster(20)
t1 = Sys.time()
m.max.otl.seq = parLapply(cl, data.otl.list, connected.ksccd.m)
stopCluster(cl)
m.max.otl.seq = unlist(m.max.otl.seq)
t2 = Sys.time()
t2-t1

print(m.max.otl.seq)
length(which(m.max.otl.seq<m.max.0.05))/len # the percentage of the data.otl where mutiple components (outliers or clusters) could be detected.
