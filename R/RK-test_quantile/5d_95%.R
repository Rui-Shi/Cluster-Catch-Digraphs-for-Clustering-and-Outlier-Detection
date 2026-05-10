tt1 = Sys.time()
source(here::here("R/ccds/NN_CCD.R"))
source(here::here("R/ccds/mKNN_CCD_functions.R"))
source(here::here("R/general_functions/count.R"))
source(here::here("R/general_functions/Uni-Gau_cls.R"))
setwd(here::here("R/general_functions/NN-test_quantile"))
library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 5
iteN = 10000
quant = 0.95
n = 1000
simul = NNDestP.simpois.lower.quant(n, d, quant=quant, niter=iteN, shape="sphere")

save(simul, file = "NN-test-simul_5d_90%.RData")


tt2 = Sys.time()
tt2-tt1