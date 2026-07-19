#!/usr/bin/env Rscript
# revision_experiments/04b_validate_par_radi.R -- bit-exactness check for the
# parallel nnccd.radi override ported into 04_wp4_runtime.R (task: redesign
# WP4 grids to eliminate empty cells / add 22-core parallel UNCCD radius
# search). Confirms the parallel port (arithmetic ported from
# 07b_wp5_fulldata_ccd.R, WITHOUT its chunk disk cache) reproduces
# harness.R's clamped SERIAL nnccd.radi exactly on a small synthetic dataset
# using 04_wp4_runtime.R's own generator: n=300, d=10.
#
# USAGE (PowerShell; Rscript under Bash segfaults):
#   Rscript "revision_experiments/04b_validate_par_radi.R" [cores]
#   (cores defaults to 4 -- validation runs are capped at 4 cores so the
#   concurrently-running 20-core full-data Musk/Speech chain is not starved)
#
# Prints "[PASS]"/"[FAIL]" for radii and OOS scores, then "VALIDATE PASS" /
# "VALIDATE FAIL" (exit status 1 on failure). Writes nothing to disk.

suppressPackageStartupMessages({ library(parallel); library(here) })

args  <- commandArgs(trailingOnly = TRUE)
CORES <- { ci <- suppressWarnings(as.integer(args[[1]])); if (length(args) >= 1 && !is.na(ci)) ci else 4L }

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

source(here::here("revision_experiments/harness.R"))

cat(sprintf("=== 04b validate === synthetic n=300,d=10, serial-vs-parallel, cores=%d\n", CORES))

# ---------------------------------------------------------------------------
# Capture harness.R's clamped SERIAL nnccd.radi, then define the parallel
# override (verbatim port from 04_wp4_runtime.R / 07b_wp5_fulldata_ccd.R,
# minus 07b's chunk disk cache).
# ---------------------------------------------------------------------------
nnccd.radi.serial <- nnccd.radi

nnccd.radi.parallel <- function(dx, quantile = "lower", method = "ascend", low.num, quant,
                                simul = NULL, niter, scores = F) {
  ddx <- as.matrix(dist(dx))
  n <- nrow(dx)
  d <- ncol(dx)

  if (quantile != "lower") stop("parallel nnccd.radi: only quantile='lower' supported (as used by all callers)")

  if (!is.null(simul)) {
    avg_len <- length(simul$average)
    med_len <- length(simul$median)
    NN.envelop <- list(
      average = simul$average[pmin(1:n, avg_len)],
      median  = simul$median[pmin(1:n, med_len)]
    )
  } else {
    NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)
  }

  radi_one <- function(i) {
    R_i <- 0
    if (method == "ascend") {
      o.d <- order(ddx[i, ])
      for (j in low.num:n) {
        r <- ddx[i, o.d[j]]
        NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j], o.d[2:j]], r)
        lower.bound.ave <- NN.envelop$average[j - 1]
        lower.bound.med <- NN.envelop$median[j - 1]
        if (NN.dist.obs$averge < lower.bound.ave | NN.dist.obs$median < lower.bound.med) {
          if (j == low.num) R_i <- 0
          else R_i <- ddx[i, o.d[j - 1]]
          break
        }
      }
    }
    if (method == "descend") {
      o.d <- order(ddx[i, ], decreasing = T)
      for (j in 1:(n - low.num)) {
        r <- ddx[i, o.d[j]]
        NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n - 1)], o.d[j:(n - 1)]], r)
        lower.bound.ave <- rev(NN.envelop$average)[j + 2]
        lower.bound.med <- rev(NN.envelop$median)[j + 2]
        if (NN.dist.obs$averge > lower.bound.ave & NN.dist.obs$median > lower.bound.med) {
          R_i <- r
          break
        }
      }
    }
    if (scores && R_i == 0) R_i <- sort(ddx[i, ])[2]
    R_i
  }

  cl <- makeCluster(CORES)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("ddx", "n", "low.num", "method", "scores", "NN.envelop", "NNDest.dist.f"),
                envir = environment())
  R <- unname(parSapplyLB(cl, 1:n, radi_one))
  return(list(R = R, KS = NULL))
}

# ---------------------------------------------------------------------------
# Synthetic dataset: 04_wp4_runtime.R's own generator, n=300, d=10.
# ---------------------------------------------------------------------------
gen_gaussian_2cls_runtime <- function(n, d, cont = 0.05, seed,
                                       cls_dis = 3, otl_dis = 2,
                                       r_min = 0.7, r_max = 1.3,
                                       noise_level = 0.01) {
  n1 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)
  n2 <- n - n1 - n0

  mu1 <- rep(3, d)
  mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- colMeans(rbind(mu1, mu2))

  sigma <- 1 / sqrt(qchisq(1 - noise_level, d))

  set.seed(seed)
  data1 <- MASS::mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 <- MASS::mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)

  i <- 0
  outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    r1 <- sqrt(sum((temp - mu1)^2))
    r2 <- sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis && r2 > otl_dis) {
      outlier <- rbind(outlier, temp)
      i <- i + 1
    }
  }
  rownames(outlier) <- NULL

  X <- rbind(data1, data2, outlier)
  colnames(X) <- paste0("V", seq_len(d))
  Y <- c(rep(1, n1 + n2), rep(0, n0))
  list(X = X, Y = Y)
}

dat <- gen_gaussian_2cls_runtime(300, 10, seed = 314159L)
tab <- get_simul("NN", 10)
dir <- if (10 <= 5) "ascend" else "descend"   # matches 04_wp4_runtime.R's unccd_dir_for_d(10)

# ---------------------------------------------------------------------------
# 1. Radii bit-exactness
# ---------------------------------------------------------------------------
cat("serial radii...\n"); t0 <- Sys.time()
Rs <- nnccd.radi.serial(dat$X, quantile = "lower", method = dir, low.num = 3,
                        quant = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE)
cat(sprintf("  serial wall %.2f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat(sprintf("parallel radii (cores=%d)...\n", CORES)); t0 <- Sys.time()
Rp <- nnccd.radi.parallel(dat$X, quantile = "lower", method = dir, low.num = 3,
                          quant = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE)
cat(sprintf("  parallel wall %.2f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

ok_r <- identical(Rs$R, Rp$R)
cat(sprintf("[%s] radii identical (n=%d)\n", if (ok_r) "PASS" else "FAIL", nrow(dat$X)))

# ---------------------------------------------------------------------------
# 2. Full OOS score bit-exactness (parallel construction vs serial
#    construction, using the SAME override-swap trick as 07b's validate mode)
# ---------------------------------------------------------------------------
nnccd.radi <- nnccd.radi.parallel
cat("full OOS score via parallel construction...\n")
sc_p <- NNCCD_OOS(datax = dat$X, simul = tab$simul, method = dir, d = 10)

nnccd.radi <- nnccd.radi.serial
cat("full OOS score via serial construction...\n")
sc_s <- NNCCD_OOS(datax = dat$X, simul = tab$simul, method = dir, d = 10)

ok_s <- identical(sc_s, sc_p)
cat(sprintf("[%s] OOS scores identical\n", if (ok_s) "PASS" else "FAIL"))

if (ok_r && ok_s) {
  cat("VALIDATE PASS\n")
} else {
  cat("VALIDATE FAIL\n")
  quit(status = 1)
}
