#!/usr/bin/env Rscript
# 01d_nn_serial_iteration_probe.R  (T3 Phase A support)
# Real-code-path timing check for NN at very high d, where a full-size probe
# (n=1000, niter=20, cores=20) is infeasible within the probe budget:
# one iteration at n=1000, d=1555 costs ~1 h (eigen-dominated) and 20 workers
# would need ~125 GB RAM (~6.2 GB/worker for the per-iteration dataset list).
#
# Instead this runs the ORIGINAL serial function NNDest.simpois.lower.quant
# (R/ccds/NN_Dist_Est.R, unmodified) with niter=1 at a REDUCED n, times it,
# and compares against the component model of 01b_nn_component_probe.R
# evaluated at the same reduced n. Agreement validates the model, which is
# then used to extrapolate the production (n=1000) per-iteration cost.
#
# Usage: Rscript "revision_experiments/tr2/01d_nn_serial_iteration_probe.R" <d> [n=100]

suppressPackageStartupMessages({
  library(here)
  library(MASS)
})
source(here::here("R/ccds/NN_Dist_Est.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: 01d_nn_serial_iteration_probe.R <d> [n]")
d <- as.integer(args[[1]])
n <- if (length(args) >= 2) as.integer(args[[2]]) else 100L

cat(sprintf("=== NN serial single-iteration probe: ORIGINAL NNDest.simpois.lower.quant, d=%d n=%d niter=1 ===\n", d, n))
set.seed(123)
t0 <- Sys.time()
res <- NNDest.simpois.lower.quant(n, d, quant = 0.999, niter = 1, shape = "sphere")
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("names: %s | average length %d | median length %d | any NA: %s\n",
            paste(names(res), collapse = ","), length(res$average),
            length(res$median), any(is.na(res$average)) || any(is.na(res$median))))
cat(sprintf("MEASURED_ONE_ITERATION_SECONDS=%.3f (at n=%d, d=%d)\n", elapsed, n, d))
cat("Compare with 01b component model evaluated at this n (eigen scales ~n calls,\n")
cat("gen ~ n(n+1)/2 draws, dist ~ sum x^2). Model validation at production n=1000\n")
cat("then justifies extrapolating iteration cost without running the full hour.\n")
