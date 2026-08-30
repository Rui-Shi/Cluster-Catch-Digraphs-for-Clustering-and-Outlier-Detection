#!/usr/bin/env Rscript
# 01c_validate_rk_stable.R  (T3 Phase A support)
# Validates the log-space-stable RK weight computation used by
# 01_gen_quantile_table.R for d >= 342 against the ORIGINAL implementation
# (Kest.simpois.edge.quantile in R/ccds/Kest.R) at dimensions where the
# original works (d = 20, d = 100).
#
# Method: serial (non-parallel) versions of both, run with the SAME RNG seed
# so they see literally identical simulated data; only the edge-correction
# weight computation differs (direct gamma()/integrate(sin^d) vs log-space).
# integrate()'s default rel.tol is ~1.2e-4, so agreement is asserted at
# max relative difference < 1e-3 (observed differences should be far smaller).
#
# Also smoke-checks the stable path at d = 344 and d = 1555 (original is NaN
# there): output must be finite, NaN-free, and quantile() must succeed.
#
# Usage: Rscript "revision_experiments/01c_validate_rk_stable.R"

suppressPackageStartupMessages({
  library(here)
  library(MASS)
})

source(here::here("R/ccds/Kest.R"))  # original Kest.simpois.edge.quantile + rpoisball.unit

# Serial twin of the stable computation in 01_gen_quantile_table.R
# (identical math and identical RNG-consumption order to the serial original).
Kest.simpois.edge.quantile.stable.serial <- function(m, d, rn, quan, niter) {
  r <- seq(1 / rn, 1, 1 / rn)
  Kest.m <- NULL
  for (i in 1:niter) {
    temp <- rpoisball.unit(m, d)          # same global generator as the original
    temp.dist <- as.matrix(dist(temp))

    log_cons <- log(sqrt(pi)) + lgamma((d + 1) / 2) - log(2) - lgamma(d / 2 + 1)
    log_ftemp_one <- function(t) {
      b <- acos(min(max(t / 2, -1), 1))
      if (b <= 0) return(Inf)
      lsb <- log(sin(b))
      I_scaled <- integrate(function(u) exp(d * (log(sin(u)) - lsb)), 0, b)$value
      logI <- d * lsb + log(I_scaled)
      return(log_cons - logI)
    }
    ftemp <- sapply(temp.dist, log_ftemp_one, simplify = TRUE)
    ftemp <- matrix(exp(ftemp), nrow = nrow(temp.dist), byrow = FALSE)

    diag(temp.dist) <- Inf
    result <- sapply(r, function(x) {
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      W <- ftemp
      W[Mtemp == 0] <- 0
      sumM <- cumsum(2 * colSums(W))
      return(sumM / ((1:m) * (1:m)))
    }, simplify = TRUE)
    Kest.m <- rbind(Kest.m, as.vector(result))
  }

  Kest.quan <- list()
  for (cur_quan in quan) {
    temp <- apply(Kest.m, 2, quantile, probs = as.numeric(cur_quan))
    Kest.quan[[as.character(cur_quan)]] <- matrix(temp, nrow = m)
  }
  return(list(Kest.m = Kest.m, quan = Kest.quan, r = r))
}

compare_at <- function(d, m = 200, rn = 10, quan = 0.999, niter = 2, seed = 42) {
  cat(sprintf("\n--- Equivalence check at d=%d (m=%d, niter=%d, seed=%d) ---\n",
              d, m, niter, seed))
  set.seed(seed)
  t0 <- Sys.time()
  orig <- Kest.simpois.edge.quantile(m, d, rn, quan, niter)
  t_orig <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  set.seed(seed)
  t0 <- Sys.time()
  stab <- Kest.simpois.edge.quantile.stable.serial(m, d, rn, quan, niter)
  t_stab <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  stopifnot(identical(dim(orig$Kest.m), dim(stab$Kest.m)))
  nz <- which(orig$Kest.m != 0 | stab$Kest.m != 0)
  abs_diff <- abs(orig$Kest.m - stab$Kest.m)
  rel_diff <- if (length(nz)) {
    max(abs_diff[nz] / pmax(abs(orig$Kest.m[nz]), .Machine$double.xmin))
  } else 0
  cat(sprintf("nonzero cells: %d of %d | max abs diff: %.3e | max rel diff: %.3e\n",
              length(nz), length(orig$Kest.m), max(abs_diff), rel_diff))
  cat(sprintf("quan matrices max abs diff: %.3e\n",
              max(abs(orig$quan[[1]] - stab$quan[[1]]))))
  cat(sprintf("timing serial: original %.2f s | stable %.2f s (per %d iters)\n",
              t_orig, t_stab, niter))
  stopifnot("stable/original mismatch beyond integrate() tolerance" = rel_diff < 1e-3)
  cat(sprintf("PASS: stable == original at d=%d (rel tol 1e-3)\n", d))
  invisible(list(rel_diff = rel_diff))
}

smoke_at <- function(d, m = 200, rn = 10, quan = 0.999, niter = 1, seed = 7) {
  cat(sprintf("\n--- Stable-only smoke check at d=%d (m=%d, niter=%d) ---\n", d, m, niter))
  # Show the original's constant is indeed broken here:
  cons <- (sqrt(pi) * gamma((d + 1) / 2)) / (2 * gamma(d / 2 + 1))
  cat(sprintf("original cons at d=%d: %s (log-space cons: %.6f)\n",
              d, format(cons),
              exp(log(sqrt(pi)) + lgamma((d + 1) / 2) - log(2) - lgamma(d / 2 + 1))))
  set.seed(seed)
  res <- Kest.simpois.edge.quantile.stable.serial(m, d, rn, quan, niter)
  stopifnot(
    "stable Kest.m has NA/NaN"  = !any(is.na(res$Kest.m)),
    "stable Kest.m has Inf"     = all(is.finite(res$Kest.m)),
    "stable quan has NA/NaN"    = !any(is.na(res$quan[[1]]))
  )
  cat(sprintf("PASS: finite, NaN-free output at d=%d | Kest.m frac zero: %.4f | max: %s\n",
              d, mean(res$Kest.m == 0), format(max(res$Kest.m))))
}

compare_at(20)
compare_at(100)
smoke_at(344)
smoke_at(1555)

cat("\nALL RK STABLE-OVERRIDE VALIDATION CHECKS PASSED\n")
