#!/usr/bin/env Rscript
tt1 = Sys.time()
source(here::here("R/ccds/UN_CCD.R"))
source(here::here("R/ccds/mKNN_CCD_functions.R"))
source(here::here("R/general_functions/count.R"))
source(here::here("R/general_functions/Uni-Gau_cls.R"))
setwd(here::here("R/NN-test_quantile"))
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 9
iteN = 2000
quant = 0.95
n = 5000
simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")
save(simul, file = "NN-test-simul_9d_95%.RData")

quant = 0.99
save(simul, file = "NN-test-simul_9d_99%.RData")

tt2 = Sys.time()
tt2-tt1