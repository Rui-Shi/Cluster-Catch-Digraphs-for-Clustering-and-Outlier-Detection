#!/usr/bin/env Rscript
tt1 = Sys.time()
source("/mmfs1/home/rzs0112/code_working_folder/ccds/NN_CCD.R")
source("/mmfs1/home/rzs0112/code_working_folder/ccds/mKNN_CCD_functions.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/count.R")
source("/mmfs1/home/rzs0112/code_working_folder/general functions/Uni-Gau_cls.R")
setwd("/mmfs1/home/rzs0112/code_working_folder/general functions/NN-test_quantile")
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 11
iteN = 2000
quant = 0.99
n = 5000
# simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")
# save(simul, file = "NN-test-simul_11d_99%.RData")

quant = 0.999
simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")
save(simul, file = "NN-test-simul_11d_999%.RData")

tt2 = Sys.time()
tt2-tt1