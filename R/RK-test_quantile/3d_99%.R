tt1 = Sys.time()
source(here::here("R/ccds/RK_CCD_New.R"))
source(here::here("R/ccds/mKNN_CCD_functions.R"))
source(here::here("R/general_functions/count.R"))
source(here::here("R/general_functions/Uni-Gau_cls.R"))
setwd(here::here("R/general_functions/RK-test_quantile"))

library(parallel)
library(doParallel)
library(MASS)
library(igraph)

d = 3
rn = 10
quan = 0.99
niter = 5000
n = 1000

# set.seed(1234)
# cores = detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# simul.list = foreach(x=c(1:n),.packages = c("MASS","cluster","igraph")) %dopar% 
#   Kest.simpois.edge.quantile(x, d, rn, quan, niter)
# stopCluster(cl)

simul = KestP.simpois.edge.quantile(n, d, rn, quan, niter)

save(simul, file = "RK-test-simul_3d_99%.RData")

tt2 = Sys.time()
tt2-tt1