tt1 = Sys.time()
source("/media/rui/exNVME/code_working_folder/ccds/RK_CCD_New.R")
source("/media/rui/exNVME/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/media/rui/exNVME/code_working_folder/general functions/count.R")
source("/media/rui/exNVME/code_working_folder/general functions/Uni-Gau_cls.R")
setwd("/media/rui/exNVME/code_working_folder/general functions/RK-test_quantile")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 20
rn = 10
quan = 0.999
niter = 10000
n = 1000

# set.seed(1234)
# cores = detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# simul.list = foreach(x=c(1:n),.packages = c("MASS","cluster","igraph")) %dopar% 
# Kest.simpois.edge.quantile(x, d, rn, quan, niter)
# stopCluster(cl)

simul = KestP.simpois.edge.quantile(n, d, rn, quan, niter)

save(simul, file = "RK-test-simul_20d_999%.RData")

tt2 = Sys.time()
tt2-tt1