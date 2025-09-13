#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/ccds/RK_CCD_New.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
setwd("/mmfs1/home/rzs0112/code_working_folder/general functions/RK-test_quantile")

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 26
rn = 10
quan = 0.999
niter = 2000
n = 5000

# set.seed(1234)
# cores = detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# simul.list = foreach(x=c(1:n),.packages = c("MASS","cluster","igraph")) %dopar% 
#   Kest.simpois.edge.quantile(x, d, rn, quan, niter)
# stopCluster(cl)

simul = KestP.simpois.edge.quantile(n, d, rn, quan, niter)
save(simul, file = "RK-test-simul_26d_999%.RData")


tt2 = Sys.time()
tt2-tt1