#!/usr/bin/env Rscript
# 01f_validate_nn_fast.R  (T3 Phase B-prep)
# Validates the two NN generator-side optimizations in
# revision_experiments/tr2/01e_nn_fast.R against the ORIGINAL
# NNDestP.simpois.lower.quant (R/ccds/NN_Dist_Est.R, unmodified):
#   (1) rpoisball.unit.fast(): matrix(rnorm(n*d),n,d) in place of
#       mvrnorm(n, rep(0,d), diag(d))
#   (2) NNDestP.simpois.lower.quant.fast(): streaming SimuOnce (no
#       pre-materialized data.simu.list)
#
# Checks (all validation-scale, minutes total):
#   A. rnorm-swap: bitwise or merely statistical equivalence? (instant)
#   B. Statistical equivalence of the full pipeline, original vs fast,
#      niter=1000 @ d=10 and niter=500 @ d=50 (serial, small n) (~8 min).
#      Tolerance yardstick: run 01h_nn_mc_yardstick.R (original vs original,
#      different seeds, same scales) -- the orig-vs-fast differences must be
#      comparable to that pure-Monte-Carlo baseline.
#   C. Timing at d=166: HEAD-TO-HEAD original vs fast, back-to-back in this
#      process (n=1000, niter=20, cores=20). Back-to-back cancels ambient
#      machine load, unlike comparing against the older clean probe log
#      (results/probes/NN_166_probe.log: 71.244 s, 2026-07-10), which is
#      still printed for reference. (~2 x 1-3 min depending on load)
#   D. RAM + timing at d=400: ONE serial SimuOnce-equivalent call (n=1000 --
#      exactly what one %dopar% worker executes per assigned iteration) for
#      original vs fast, each in its OWN CHILD Rscript process so Windows
#      peak working set (ps::ps_memory_info 'peak_wset'; monotonic within a
#      process) reflects that variant alone. (~2 x 2-6 min)
#
# Usage: Rscript "revision_experiments/tr2/01f_validate_nn_fast.R" [sections]
#   sections: subset of "abcd" (default "abcd"), e.g. "cd" reruns only the
#   timing/RAM checks. (Internal child modes: "_d_orig", "_d_fast".)

suppressPackageStartupMessages({
  library(here)
  library(parallel)
  library(doParallel)
  library(MASS)
})
# `ps` is only needed for the RAM measurements (section D / child modes);
# sections A-C run without it (e.g. under an R install lacking ps).
HAVE_PS <- requireNamespace("ps", quietly = TRUE)

cli <- commandArgs(trailingOnly = TRUE)
SECTIONS <- if (length(cli) >= 1) tolower(cli[[1]]) else "abcd"
run_sec <- function(s) grepl(s, SECTIONS, fixed = TRUE)

source(here::here("R/ccds/NN_Dist_Est.R"))               # originals, unmodified
source(here::here("revision_experiments/tr2/01e_nn_fast.R")) # fast overrides

# ---------------------------------------------------------------------
# Child modes for section D (run in a fresh process so peak_wset is clean)
# ---------------------------------------------------------------------
mem_line <- function(tag) {
  if (!HAVE_PS) {
    cat(sprintf("%s: [ps package not installed -- RAM not measured]\n", tag))
    return(invisible(NA_real_))
  }
  mi <- ps::ps_memory_info(ps::ps_handle())
  peak <- if ("peak_wset" %in% names(mi)) as.numeric(mi[["peak_wset"]]) / 1024^2 else NA_real_
  rss <- as.numeric(mi[["rss"]]) / 1024^2
  cat(sprintf("%s: rss %.1f MB | peak_wset %.1f MB\n", tag, rss, peak))
  invisible(peak)
}

if (SECTIONS %in% c("_d_orig", "_d_fast")) {
  n <- 1000L; d <- 400L
  set.seed(99)
  mem_line("BASELINE (post-load, pre-run)")
  t0 <- Sys.time()
  if (SECTIONS == "_d_orig") {
    # Original SimuOnce body (NNDestP.simpois.lower.quant, verbatim shape):
    # full data.simu.list materialized first.
    data.simu.list <- lapply(1:n, rpoisball.unit, d = d)
    mem_line("AFTER list materialization")
    NN.dist.temp <- sapply(2:n, function(x) {
      data.temp <- data.simu.list[[x]]
      data.dist <- as.matrix(dist(data.temp))
      diag(data.dist) <- Inf
      NN.dist.ttemp <- apply(data.dist, 1, min)
      c(mean(NN.dist.ttemp), median(NN.dist.ttemp))
    })
  } else {
    # Streaming SimuOnce body (01e_nn_fast.R, verbatim shape).
    NN.dist.temp <- sapply(2:n, function(x) {
      data.temp <- rpoisball.unit.fast(x, d)
      data.dist <- as.matrix(dist(data.temp))
      diag(data.dist) <- Inf
      NN.dist.ttemp <- apply(data.dist, 1, min)
      c(mean(NN.dist.ttemp), median(NN.dist.ttemp))
    })
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("CHILD_ELAPSED_SECONDS=%.3f\n", elapsed))
  cat(sprintf("CHILD_RESULT_DIM=%s | any NA: %s\n",
              paste(dim(NN.dist.temp), collapse = "x"), any(is.na(NN.dist.temp))))
  mem_line("AFTER SimuOnce")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------
# A. rnorm-swap: bitwise identical under the same seed?
# ---------------------------------------------------------------------
if (run_sec("a")) {
  cat("=== A. rnorm-swap: bitwise vs statistical (same seed, d=10) ===\n")
  set.seed(42); a <- mvrnorm(50, rep(0, 10), diag(10))
  set.seed(42); b <- matrix(rnorm(500), 50, 10)
  cat(sprintf(
    "identical(): %s | max abs diff: %.4f | same multiset of draws: %s\n",
    identical(a, b), max(abs(a - b)),
    isTRUE(all.equal(sort(as.vector(a)), sort(as.vector(b))))
  ))
  eS <- eigen(diag(10), symmetric = TRUE)
  cat(sprintf(
    "eigen(diag(10))$vectors == identity: %s  (rotation != identity is WHY draws\n  differ bitwise; QZ ~ N(0,I) for ANY orthogonal Q when Z ~ N(0,I), so the\n  distribution is unaffected -- rotational invariance of the isotropic Gaussian)\n",
    isTRUE(all.equal(eS$vectors, diag(10)))
  ))
}

# ---------------------------------------------------------------------
# B. Statistical equivalence, full pipeline, original vs fast (serial)
# ---------------------------------------------------------------------
NNDest.simpois.lower.quant.fast.serial <- function(n, d, quant = 0.99, niter = 100) {
  NN.dist.ave.mat <- NULL
  NN.dist.med.mat <- NULL
  for (i in 1:niter) {
    NN.dist.temp <- sapply(2:n, function(x) {
      data.temp <- rpoisball.unit.fast(x, d)
      data.dist <- as.matrix(dist(data.temp))
      diag(data.dist) <- Inf
      NN.dist.ttemp <- apply(data.dist, 1, min)
      c(mean(NN.dist.ttemp), median(NN.dist.ttemp))
    })
    NN.dist.ave.mat <- rbind(NN.dist.ave.mat, c(0, NN.dist.temp[1, ]))
    NN.dist.med.mat <- rbind(NN.dist.med.mat, c(0, NN.dist.temp[2, ]))
  }
  list(
    average = sapply(1:n, function(x) unname(quantile(NN.dist.ave.mat[, x], 1 - quant))),
    median  = sapply(1:n, function(x) unname(quantile(NN.dist.med.mat[, x], 1 - quant)))
  )
}

compare_stat <- function(d, n, niter, quant = 0.999, seed = 1) {
  cat(sprintf("\n--- B. Statistical equivalence at d=%d n=%d niter=%d quant=%s ---\n", d, n, niter, quant))
  set.seed(seed)
  t0 <- Sys.time()
  orig <- NNDest.simpois.lower.quant(n, d, quant = quant, niter = niter, shape = "sphere")
  t_orig <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Deliberately a DIFFERENT seed: this is a distributional-agreement check
  # (does the fast path sample the same population?), not a paired-draw
  # check -- the RNG call sequences differ between the two code paths
  # regardless of seed (see check A), so pairing on seed would not mean
  # pairing on draws anyway.
  set.seed(seed + 1)
  t0 <- Sys.time()
  fast <- NNDest.simpois.lower.quant.fast.serial(n, d, quant = quant, niter = niter)
  t_fast <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rel_ave <- abs(orig$average - fast$average) / pmax(abs(orig$average), 1e-8)
  rel_med <- abs(orig$median - fast$median) / pmax(abs(orig$median), 1e-8)
  cat(sprintf("average: max rel diff %.4f | mean rel diff %.4f | Pearson cor %.6f\n",
              max(rel_ave), mean(rel_ave), cor(orig$average, fast$average)))
  cat(sprintf("median : max rel diff %.4f | mean rel diff %.4f | Pearson cor %.6f\n",
              max(rel_med), mean(rel_med), cor(orig$median, fast$median)))
  cat(sprintf("timing serial (niter=%d): original %.2f s | fast %.2f s | fast/orig speedup %.2fx\n",
              niter, t_orig, t_fast, t_orig / t_fast))
  invisible(list(rel_ave = rel_ave, rel_med = rel_med))
}

if (run_sec("b")) {
  compare_stat(d = 10, n = 200, niter = 1000)
  compare_stat(d = 50, n = 200, niter = 500)
}

# ---------------------------------------------------------------------
# C. Timing at d=166: HEAD-TO-HEAD original vs fast, same load
# ---------------------------------------------------------------------
if (run_sec("c")) {
  cat("\n=== C. Head-to-head at d=166 (n=1000, niter=20, cores=20) ===\n")
  cat("Clean reference for context (results/probes/NN_166_probe.log,\n")
  cat("2026-07-10, original, niter=20 cores=20, lightly loaded): 71.244 s\n")

  requested_cores <- 20L
  detectCores <- function(...) requested_cores  # nolint: intentional shadow, mirrors 01_gen_quantile_table.R

  t0 <- Sys.time()
  o166 <- NNDestP.simpois.lower.quant(1000, 166, quant = 0.999, niter = 20, shape = "sphere")
  t_orig_166 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("ORIGINAL: %.3f s | any NA: %s\n", t_orig_166,
              any(is.na(o166$average)) || any(is.na(o166$median))))

  t0 <- Sys.time()
  f166 <- NNDestP.simpois.lower.quant.fast(1000, 166, quant = 0.999, niter = 20, shape = "sphere")
  t_fast_166 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("FAST    : %.3f s | any NA: %s\n", t_fast_166,
              any(is.na(f166$average)) || any(is.na(f166$median))))

  cat(sprintf("HEADTOHEAD_SPEEDUP_166=%.4f\n", t_orig_166 / t_fast_166))
  cat(sprintf("HEADTOHEAD_T_ORIG_166=%.3f\nHEADTOHEAD_T_FAST_166=%.3f\n", t_orig_166, t_fast_166))
}

# ---------------------------------------------------------------------
# D. RAM + timing at d=400 via fresh child processes (clean peak_wset)
# ---------------------------------------------------------------------
if (run_sec("d")) {
  cat("\n=== D. Per-worker RAM + timing at d=400, n=1000 (one SimuOnce each, fresh child process each) ===\n")
  me <- here::here("revision_experiments/tr2/01f_validate_nn_fast.R")
  rscript <- file.path(R.home("bin"), "Rscript")
  for (mode in c("_d_orig", "_d_fast")) {
    cat(sprintf("--- child %s ---\n", mode))
    out <- system2(rscript, args = c(shQuote(me), mode), stdout = TRUE, stderr = TRUE)
    cat(paste(out, collapse = "\n"), "\n")
  }
  cat("Interpretation: compare each child's 'AFTER SimuOnce' peak_wset against\n")
  cat("its own BASELINE -- the difference is that variant's per-worker footprint.\n")
}

cat("\nALL NN FAST-OVERRIDE VALIDATION CHECKS COMPLETE\n")
