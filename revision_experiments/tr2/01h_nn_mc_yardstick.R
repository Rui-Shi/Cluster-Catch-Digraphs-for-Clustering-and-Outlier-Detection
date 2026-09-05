#!/usr/bin/env Rscript
# 01h_nn_mc_yardstick.R  (T3 Phase B-prep)
# Monte-Carlo-error yardstick for the equivalence tolerance used in
# 01f_validate_nn_fast.R: runs the ORIGINAL NNDest.simpois.lower.quant
# twice with DIFFERENT seeds at the same scales used there
# (d=10 n=200 niter=1000; d=50 n=200 niter=500) and reports the same
# max-rel-diff / correlation statistics. The orig-vs-FAST differences are
# acceptable iff they are comparable to this orig-vs-ORIG baseline --
# i.e. indistinguishable from pure Monte-Carlo seed noise.
#
# Usage: Rscript "revision_experiments/tr2/01h_nn_mc_yardstick.R"

suppressPackageStartupMessages({
  library(here)
  library(MASS)
})
source(here::here("R/ccds/NN_Dist_Est.R"))

yardstick <- function(d, n, niter, quant = 0.999, seed1 = 101, seed2 = 202) {
  cat(sprintf("\n--- MC yardstick: ORIGINAL vs ORIGINAL, d=%d n=%d niter=%d, seeds %d vs %d ---\n",
              d, n, niter, seed1, seed2))
  set.seed(seed1)
  a <- NNDest.simpois.lower.quant(n, d, quant = quant, niter = niter, shape = "sphere")
  set.seed(seed2)
  b <- NNDest.simpois.lower.quant(n, d, quant = quant, niter = niter, shape = "sphere")
  rel_ave <- abs(a$average - b$average) / pmax(abs(a$average), 1e-8)
  rel_med <- abs(a$median - b$median) / pmax(abs(a$median), 1e-8)
  cat(sprintf("average: max rel diff %.4f | mean rel diff %.4f | Pearson cor %.6f\n",
              max(rel_ave), mean(rel_ave), cor(a$average, b$average)))
  cat(sprintf("median : max rel diff %.4f | mean rel diff %.4f | Pearson cor %.6f\n",
              max(rel_med), mean(rel_med), cor(a$median, b$median)))
}

yardstick(d = 10, n = 200, niter = 1000)
yardstick(d = 50, n = 200, niter = 500)

cat("\nMC YARDSTICK DONE\n")
