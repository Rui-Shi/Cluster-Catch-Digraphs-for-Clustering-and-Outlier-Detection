tt1 = Sys.time()
source("/media/rui/exNVME/code_working_folder/ccds/NN_CCD.R")
source("/media/rui/exNVME/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/media/rui/exNVME/code_working_folder/general functions/count.R")
source("/media/rui/exNVME/code_working_folder/general functions/Uni-Gau_cls.R")
setwd("/media/rui/exNVME/code_working_folder/general functions/NN-test_quantile")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 2
iteN = 10000
quant = 0.85
n = 1000
simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")
save(simul, file = "NN-test-simul_2d_85%.RData")

quant = 0.90
simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")
save(simul, file = "NN-test-simul_2d_90%.RData")

tt2 = Sys.time()
tt2-tt1