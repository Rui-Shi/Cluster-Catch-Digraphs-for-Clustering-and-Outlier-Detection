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

# Define constants 
iteN = 10000
quant = 0.999
n = 1000

cat(sprintf("Starting sequential loop. Inner function will handle parallelization.\n"))

# Standard for loop
for (d in 35:99) {
  
  cat(sprintf("Running simulation for d = %d at %s\n", d, Sys.time()))
  
  # The function will use its own internal parallelization here
  simul = NNDestP.simpois.lower.quant(n, d, quant = quant, niter = iteN, shape = "sphere")
  
  # Dynamically generate the file name and save
  file_name = paste0("NN-test-simul_", d, "d_999%.RData")
  save(simul, file = file_name)
  
}

tt2 = Sys.time()
print(tt2 - tt1)