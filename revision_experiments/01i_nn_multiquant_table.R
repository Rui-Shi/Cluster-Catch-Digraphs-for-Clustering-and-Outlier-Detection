#!/usr/bin/env Rscript
# 01i_nn_multiquant_table.R -- one Monte-Carlo run, ALL alpha levels.
#
# WHY THIS FILE EXISTS
# --------------------
# R/ccds/NN_Dist_Est.R's generators (NNDest.simpois.lower.quant, line 18, and
# its parallel twin NNDestP.simpois.lower.quant, line 56) accumulate the
# per-iteration draws into NN.dist.ave.mat / NN.dist.med.mat, reduce them to a
# SINGLE lower-tail quantile at lines 101-110, and then throw the draw matrices
# away when the function returns list(average=, median=). Every new alpha
# therefore costs a whole new Monte-Carlo run, which is the reason the on-disk
# NN tables are 70 KB while the RK tables (which DO store their raw Kest.m) are
# 280-620 MB and can be re-reduced to any alpha for free.
#
# This file removes that asymmetry for NN:
#   * nn_mc_draws()      runs the Monte Carlo ONCE and RETURNS the draw matrices
#   * nn_reduce_quant()  the reduction, verbatim from NN_Dist_Est.R:101-110
#   * nn_multiquant()    reduces one set of draws at as many alphas as you like
#   * write_nn_tables()  emits (a) one drop-in list(average=, median=) `simul`
#                        file per alpha, named exactly like the production
#                        tables so an existing loader can read them, and
#                        (b) ONE raw-draws file per (d, n), so every FUTURE
#                        alpha at that d is free -- the same contract the RK
#                        generator already offers through $Kest.m.
#
# R/ccds/NN_Dist_Est.R, 01e_nn_fast.R and 01_gen_quantile_table.R are NOT
# modified. This is a new, additional generator.
#
# TWO ENGINES
# -----------
# engine = "orig_list" (DEFAULT, and what the tables shipped with this study
#     were built with): the per-iteration body is a verbatim transcription of
#     the ORIGINAL NNDest.simpois.lower.quant loop body (NN_Dist_Est.R:24-35) --
#     the same `lapply(1:n, rpoisball.unit, d=d)` prefetch, calling the same
#     top-level mvrnorm-based rpoisball.unit() from NN_Dist_Est.R. Run serially
#     (cores = 1) under set.seed(seed) it consumes the RNG in exactly the same
#     order as the original, so its draw matrices are BIT-IDENTICAL to the
#     original's. That is what makes the exactness test in
#     48_wp2a_nn_saturation.R --validate possible at all.
#
# engine = "fast_stream": the optimized body from 01e_nn_fast.R (rnorm instead
#     of mvrnorm -- distributionally identical, proven in that file's header --
#     and streaming instead of prefetching, to cut per-worker RAM). Cheaper at
#     large n/d, but its RNG stream differs from the original by construction,
#     so it can only ever be validated statistically (01f_validate_nn_fast.R).
#     Provided here for large-n use; NOT used for this study's tables, where
#     n = 200 makes the exact engine affordable.
#
# PARALLELISM AND REPRODUCIBILITY
# -------------------------------
# cores == 1  -> plain serial for-loop after set.seed(seed). Reproducible, and
#                RNG-compatible with the original serial generator.
# cores >  1  -> parallel::parLapply over a PSOCK cluster with
#                clusterSetRNGStream(cl, seed). parLapply chunks the task list
#                statically and deterministically, so the result is
#                reproducible for a FIXED (niter, cores, seed) triple -- all
#                three are recorded in the saved metadata. It is NOT the same
#                stream as the serial path (L'Ecuyer-CMRG vs Mersenne-Twister);
#                the estimator is identical, the particular draws are not.
#
# FILENAME TOKEN
# --------------
# Production files are named e.g. NN-test-simul_10d_99%.RData for quant = 0.99
# and ..._999%.RData for 0.999, i.e. the token is the PERCENTAGE with the dot
# removed: 0.90 -> "90", 0.95 -> "95", 0.99 -> "99", 0.999 -> "999". Note that
# 01_gen_quantile_table.R:109 computes the token as sub("^0\\.", "", quant),
# which agrees for 0.99/0.999 but would emit "9" for 0.90; we use the
# percentage rule so the names match what is actually on disk.
#
# CLI
#   Rscript revision_experiments/01i_nn_multiquant_table.R \
#       <d> <n> <niter> <cores> <outdir> [seed] [quants,csv] [engine]
#   defaults: seed = 20260809 + d, quants = 0.90,0.95,0.99,0.999,
#             engine = orig_list

suppressPackageStartupMessages({
  library(here)
  library(parallel)
  library(MASS)
})

# rpoisball.unit() (mvrnorm-based) comes from here, unmodified.
source(here::here("R/ccds/NN_Dist_Est.R"))

# ---------------------------------------------------------------------------
# Filename token
# ---------------------------------------------------------------------------
nn_token_of_q <- function(q) {
  stopifnot(is.numeric(q), all(q > 0), all(q < 1))
  gsub("\\.", "", formatC(q * 100, format = "f", digits = 6, drop0trailing = TRUE))
}
nn_q_of_token <- function(tok) as.numeric(paste0("0.", tok))

# ---------------------------------------------------------------------------
# Per-iteration bodies
# ---------------------------------------------------------------------------

# ENGINE "orig_list": verbatim transcription of NNDest.simpois.lower.quant's
# loop body (R/ccds/NN_Dist_Est.R lines 24-35), lifted into a function so the
# draws can be returned instead of accumulated-and-discarded. rpoisball.unit is
# resolved from NN_Dist_Est.R (the mvrnorm version) -- redefined inside so a
# PSOCK worker, which receives this closure but not the master's .GlobalEnv,
# can resolve it. Keep textually identical to NN_Dist_Est.R:8-15.
nn_simu_once_orig <- function(n, d) {
  rpoisball.unit <- function(n, d) {
    r1 <- runif(n, 0, 1)^(1 / d)
    norm.data <- matrix(mvrnorm(n, rep(0, d), diag(d)), ncol = d, byrow = T)
    data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
    data1 <- apply(data1, 1, function(x) x * r1)
    return(data1)
  }
  data.simu.list <- lapply(1:n, rpoisball.unit, d = d)
  NN.dist.temp <- sapply(2:n, function(x) {
    data.temp <- data.simu.list[[x]]
    data.dist <- as.matrix(dist(data.temp))
    diag(data.dist) <- Inf
    NN.dist.ttemp <- apply(data.dist, 1, min)
    c(mean(NN.dist.ttemp), median(NN.dist.ttemp))
  })
  list(ave = c(0, NN.dist.temp[1, ]), med = c(0, NN.dist.temp[2, ]))
}

# ENGINE "fast_stream": the 01e_nn_fast.R body (rnorm-equivalence + streaming).
nn_simu_once_fast <- function(n, d) {
  rpoisball.unit.fast <- function(n, d) {
    r1 <- runif(n, 0, 1)^(1 / d)
    norm.data <- matrix(rnorm(n * d), nrow = n, ncol = d)
    data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
    data1 <- apply(data1, 1, function(x) x * r1)
    return(data1)
  }
  NN.dist.temp <- sapply(2:n, function(x) {
    data.temp <- rpoisball.unit.fast(x, d)
    data.dist <- as.matrix(dist(data.temp))
    diag(data.dist) <- Inf
    NN.dist.ttemp <- apply(data.dist, 1, min)
    c(mean(NN.dist.ttemp), median(NN.dist.ttemp))
  })
  list(ave = c(0, NN.dist.temp[1, ]), med = c(0, NN.dist.temp[2, ]))
}

# ---------------------------------------------------------------------------
# The Monte-Carlo run: returns the DRAWS, not a quantile.
# ---------------------------------------------------------------------------
#' @return list(ave = niter x n matrix, med = niter x n matrix, meta = list(...))
nn_mc_draws <- function(n, d, niter, cores = 1L, seed = 1L, engine = c("orig_list", "fast_stream")) {
  engine <- match.arg(engine)
  simu_once <- if (engine == "orig_list") nn_simu_once_orig else nn_simu_once_fast

  t0 <- Sys.time()
  if (cores <= 1L) {
    set.seed(seed)
    out <- vector("list", niter)
    for (i in seq_len(niter)) out[[i]] <- simu_once(n, d)
  } else {
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl), add = TRUE)
    clusterSetRNGStream(cl, seed)
    clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS)))
    # Closure over a MINIMAL environment: serializing the nn_mc_draws frame
    # itself would drag the cluster handle and the result list to every worker.
    fenv <- new.env(parent = globalenv())
    fenv$.n <- n; fenv$.d <- d; fenv$.f <- simu_once
    fn <- function(i) .f(.n, .d)
    environment(fn) <- fenv
    out <- parLapply(cl, seq_len(niter), fn)
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ave <- do.call(rbind, lapply(out, `[[`, "ave"))
  med <- do.call(rbind, lapply(out, `[[`, "med"))
  stopifnot(dim(ave) == c(niter, n), dim(med) == c(niter, n))

  list(ave = ave, med = med,
       meta = list(n = n, d = d, niter = niter, cores = cores, seed = seed,
                   engine = engine, elapsed_sec = elapsed,
                   sec_per_iter = elapsed / niter,
                   generated = format(Sys.time()),
                   R_version = R.version.string))
}

# ---------------------------------------------------------------------------
# The reduction. VERBATIM from R/ccds/NN_Dist_Est.R lines 101-110 (identical in
# the serial version at lines 41-49): the lower-tail quantile at probability
# 1 - quant, taken column by column, with names() stripped.
# ---------------------------------------------------------------------------
nn_reduce_quant <- function(ave.mat, med.mat, quant) {
  n <- ncol(ave.mat)
  quant.ave.lower <- sapply(1:n, function(x) {
    quant <- quantile(ave.mat[, x], 1 - quant)
  })
  names(quant.ave.lower) <- NULL
  quant.med.lower <- sapply(1:n, function(x) {
    quant <- quantile(med.mat[, x], 1 - quant)
  })
  names(quant.med.lower) <- NULL
  list(average = quant.ave.lower, median = quant.med.lower)
}

#' Reduce ONE set of draws at every requested quant.
#' @return named list, keys = filename tokens ("90","95","99","999"),
#'   values = list(average=, median=) -- the exact `simul` object shape the
#'   NN detectors expect.
nn_multiquant <- function(draws, quants) {
  out <- lapply(quants, function(q) nn_reduce_quant(draws$ave, draws$med, q))
  names(out) <- nn_token_of_q(quants)
  out
}

# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
#' Emit one drop-in per-alpha table file plus one raw-draws file.
#' NEVER overwrites: aborts if a target file already exists unless overwrite=TRUE.
write_nn_tables <- function(draws, quants, outdir, overwrite = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  n <- draws$meta$n; d <- draws$meta$d
  tabs <- nn_multiquant(draws, quants)

  written <- character(0)
  for (tok in names(tabs)) {
    f <- file.path(outdir, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))
    if (file.exists(f) && !overwrite) stop("refusing to overwrite existing table: ", f)
    simul <- tabs[[tok]]
    attr(simul, "nn_multiquant_meta") <- c(draws$meta, list(quant = nn_q_of_token(tok)))
    save(simul, file = f)
    written <- c(written, f)
  }

  fdraw <- file.path(outdir, sprintf("NN-draws_%dd_n%d.RData", d, n))
  if (file.exists(fdraw) && !overwrite) stop("refusing to overwrite existing draws file: ", fdraw)
  NN.dist.ave.mat <- draws$ave
  NN.dist.med.mat <- draws$med
  nn_meta <- draws$meta
  save(NN.dist.ave.mat, NN.dist.med.mat, nn_meta, file = fdraw)
  written <- c(written, fdraw)

  # zero-quantile diagnostic, per alpha
  zf <- do.call(rbind, lapply(names(tabs), function(tok) data.frame(
    d = d, n = n, niter = draws$meta$niter, token = tok, quant = nn_q_of_token(tok),
    alpha = 1 - nn_q_of_token(tok),
    frac_zero_average = mean(tabs[[tok]]$average == 0),
    frac_zero_median  = mean(tabs[[tok]]$median  == 0),
    # position 1 is the hand-added structural 0 in the original generator
    frac_zero_average_excl_pos1 = mean(tabs[[tok]]$average[-1] == 0),
    frac_zero_median_excl_pos1  = mean(tabs[[tok]]$median[-1]  == 0),
    min_nonzero_average = suppressWarnings(min(tabs[[tok]]$average[tabs[[tok]]$average > 0])),
    stringsAsFactors = FALSE)))

  list(files = written, zero_frac = zf, tables = tabs)
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
.nn_mq_is_main <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  length(f) == 1L && basename(f) == "01i_nn_multiquant_table.R"
}

if (.nn_mq_is_main()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: Rscript revision_experiments/01i_nn_multiquant_table.R <d> <n> <niter> <cores> <outdir> [seed] [quants,csv] [engine]",
         call. = FALSE)
  }
  d      <- as.integer(args[[1]])
  n      <- as.integer(args[[2]])
  niter  <- as.integer(args[[3]])
  cores  <- as.integer(args[[4]])
  outdir <- args[[5]]
  seed   <- if (length(args) >= 6) as.integer(args[[6]]) else (20260809L + d)
  quants <- if (length(args) >= 7) as.numeric(strsplit(args[[7]], ",")[[1]]) else c(0.90, 0.95, 0.99, 0.999)
  engine <- if (length(args) >= 8) args[[8]] else "orig_list"

  cat(sprintf("=== 01i_nn_multiquant_table.R ===\nd=%d n=%d niter=%d cores=%d seed=%d engine=%s\nquants=%s -> tokens %s\noutdir=%s\n",
              d, n, niter, cores, seed, engine,
              paste(quants, collapse = ","), paste(nn_token_of_q(quants), collapse = ","), outdir))

  draws <- nn_mc_draws(n, d, niter, cores, seed, engine)
  cat(sprintf("MC done: %.1f s total, %.5f s/iteration (cores=%d, contended)\n",
              draws$meta$elapsed_sec, draws$meta$sec_per_iter, cores))

  res <- write_nn_tables(draws, quants, outdir)
  cat("Wrote:\n"); cat(paste0("  ", res$files, collapse = "\n"), "\n")
  cat("\nZero-quantile fractions:\n")
  print(res$zero_frac, row.names = FALSE, digits = 4)
  cat(sprintf("TOTAL_ELAPSED_SECONDS=%.3f\n", draws$meta$elapsed_sec))
  cat("=== DONE ===\n")
}
