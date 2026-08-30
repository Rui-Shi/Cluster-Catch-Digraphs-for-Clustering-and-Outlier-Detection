#!/usr/bin/env Rscript
# revision_experiments/tr2/20_wp_collinearity.R
#
# Ad hoc collinearity-sensitivity sweep, requested to back the manuscript's
# claim that OOS/IOS are "empirically stable under the collinearity settings
# tested" (currently unsupported -- no such experiment exists in the repo).
# NOT part of the WP1-WP8 revision-plan task numbering; this is a standalone
# check. Follows the pattern of 09_wp3_synthetic.R (gaussian setting, same
# generator, same calibrated cutoffs) and 04_wp4_runtime.R (parallel
# nnccd.radi override for UNCCD, ported verbatim, needed because serial
# UNCCD-OOS/IOS construction at n=500,d=10 measured ~700s/call -- infeasible
# for a multi-rho x multi-rep sweep; the parallel override cuts this to
# ~35s/call on 20 cores, verified by a timing probe before this script was
# written).
#
# DESIGN
#   Data: two-Gaussian-cluster + uniform-outlier generator, IDENTICAL to
#     gen_gaussian() in 09_wp3_synthetic.R (n=500, d=10, cont=0.05,
#     cls_dis=3, otl_dis=2, r_min/r_max=0.7/1.3, noise_level=0.01), except
#     each cluster's covariance is changed from isotropic
#     diag(d)*(sigma*runif(1,r_min,r_max))^2 to the equicorrelated form
#     (sigma*runif(1,r_min,r_max))^2 * ((1-rho)*I_d + rho*J_d) (unit
#     variances, all off-diagonal correlations = rho). rho=0 reproduces the
#     isotropic generator exactly. Both clusters draw an independent
#     runif(1,r_min,r_max) radius scale, as in the original.
#   rho grid: 0, 0.3, 0.5, 0.7, 0.9 (equicorrelation matrix is PSD for
#     rho > -1/(d-1); comfortably satisfied up to 0.9 at d=10).
#   d = 10 only (matches the manuscript's synthetic sensitivity study).
#   Methods: RKCCD-OOS, RKCCD-IOS, UNCCD-OOS, UNCCD-IOS.
#   Cutoffs: the SAME calibrated gaussian-setting cutoffs 09_wp3_synthetic.R
#     reads from Threshold.R (d=10 -> index 4): RKCCD-OOS=3.5,
#     RKCCD-IOS=6.5, UNCCD-OOS=3.5, UNCCD-IOS=6.5. Read live from the same
#     Threshold.R files (not hard-coded) with the same stopifnot
#     cross-check 09_wp3_synthetic.R uses.
#   Replications: one dataset per (rho, rep), all 4 methods scored on the
#     SAME draw (paired design, matches 09_wp3_synthetic.R). REPS from CLI,
#     default 20 (dropped from the nominal 50 -- see header timing note).
#   Metrics: TPR/TNR/BA/F2 via harness::evaluate() at the single calibrated
#     cutoff (no multiplier sweep -- this is a stability check across rho,
#     not a cutoff-sensitivity study).
#
# CLI:  Rscript "revision_experiments/tr2/20_wp_collinearity.R" <n_reps> <cores>
#       defaults: n_reps = 20, cores = 20 (cores only affects the UNCCD
#       parallel radius search; RKCCD is single-threaded and fast).
#
# Checkpointing: rows appended per (rho, rep) to
# results/tr2/wp_collinearity_raw.csv; a (rho, rep) with all 4 method rows
# present is skipped on rerun. Aggregation (mean, sd, se over reps per
# (rho, method)) written to results/tr2/wp_collinearity_agg.csv at the end
# of every invocation.

args <- commandArgs(trailingOnly = TRUE)
N_REPS <- if (length(args) >= 1) as.integer(args[1]) else 20L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 20L
MAX_SECONDS <- if (length(args) >= 3) as.numeric(args[3]) else Inf  # wall-clock budget: finish
  # the current (rho,rep) block, checkpoint, aggregate, and exit(0) cleanly once
  # exceeded, rather than running unboundedly -- lets each invocation fit inside a
  # single foreground call with no background promotion needed.
stopifnot(!is.na(N_REPS), N_REPS >= 1, !is.na(CORES), CORES >= 1)

suppressPackageStartupMessages({
  library(here)
  library(parallel)
})

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

source(here::here("revision_experiments/shared/harness.R"))

REPO_ROOT <- here::here()
RESULTS_DIR <- file.path(REPO_ROOT, "revision_experiments/results/tr2")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
RAW_CSV <- file.path(RESULTS_DIR, "wp_collinearity_raw.csv")
AGG_CSV <- file.path(RESULTS_DIR, "wp_collinearity_agg.csv")

# ---------------------------------------------------------------------------
# Single-instance lock. Multiple concurrent processes appending to RAW_CSV
# via write.table(append=TRUE) can interleave and corrupt the CSV, so refuse
# to run if another instance already holds the lock. A lock older than 3
# hours is treated as stale (crashed holder) and reclaimed.
# ---------------------------------------------------------------------------
LOCK_FILE <- file.path(RESULTS_DIR, "wp_collinearity.lock")
if (file.exists(LOCK_FILE)) {
  age_hr <- as.numeric(difftime(Sys.time(), file.info(LOCK_FILE)$mtime, units = "hours"))
  if (age_hr < 3) {
    cat(sprintf("[LOCK] %s exists (age %.2f h) -- another instance appears to be running. Exiting without doing any work.\n",
                LOCK_FILE, age_hr))
    quit(save = "no", status = 0)
  } else {
    cat(sprintf("[LOCK] stale lock (age %.2f h) -- reclaiming.\n", age_hr))
  }
}
writeLines(as.character(Sys.getpid()), LOCK_FILE)
on.exit(unlink(LOCK_FILE), add = TRUE)

D <- 10L
D_INDEX <- 4L  # d=10 -> index 4 (vectors ordered d = 2,3,5,10,20,50,100), per 09_wp3_synthetic.R
RHOS <- c(0, 0.3, 0.5, 0.7, 0.9)
OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
UNCCD_METHODS <- c("UNCCD-OOS", "UNCCD-IOS")

cat(sprintf("20_wp_collinearity.R: n_reps=%d, cores=%d, rhos=%s\n",
            N_REPS, CORES, paste(RHOS, collapse = ",")))

# ---------------------------------------------------------------------------
# Parallel nnccd.radi override (ported verbatim from 04_wp4_runtime.R;
# UNCCD-OOS/IOS only. RKCCD is untouched and stays single-threaded.)
# ---------------------------------------------------------------------------
PAR_RADI_CORES <- CORES

nnccd.radi <- function(dx, quantile = "lower", method = "ascend", low.num, quant,
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

  cl <- makeCluster(PAR_RADI_CORES)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("ddx", "n", "low.num", "method", "scores", "NN.envelop", "NNDest.dist.f"),
                envir = environment())
  R <- unname(parSapplyLB(cl, 1:n, radi_one))
  return(list(R = R, KS = NULL))
}

# ---------------------------------------------------------------------------
# Calibrated cutoffs (gaussian setting, d=10), read live from Threshold.R,
# same files and same cross-check as 09_wp3_synthetic.R
# ---------------------------------------------------------------------------
read_thresholds <- function(rel_path) {
  e <- new.env()
  sys.source(here::here(rel_path), envir = e)
  list(OOS = get("Threshold_OOS", envir = e),
       IOS = get("Threshold_IOS", envir = e))
}
th_rk_gau <- read_thresholds("simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/Threshold.R")
th_nn_gau <- read_thresholds("simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/Gaussian/Threshold.R")

CUTOFFS <- c("RKCCD-OOS" = th_rk_gau$OOS[D_INDEX], "RKCCD-IOS" = th_rk_gau$IOS[D_INDEX],
             "UNCCD-OOS" = th_nn_gau$OOS[D_INDEX], "UNCCD-IOS" = th_nn_gau$IOS[D_INDEX])
stopifnot(identical(unname(CUTOFFS), c(3.5, 6.5, 3.5, 6.5)))
cat("Calibrated cutoffs (d=10, gaussian setting):",
    paste(sprintf("%s=%.1f", names(CUTOFFS), CUTOFFS), collapse = "  "), "\n")

# NOTE: harness.R's default rk_quant_for_d(10) currently resolves to "99",
# but 09_wp3_synthetic.R (the closest existing precedent for this exact
# gaussian/d=10 setting, whose calibrated cutoffs we reuse below) explicitly
# documents "RKCCD: quant = 0.999 (table RK-test-simul_10d_999%.RData),
# matched to the original scripts" and asserts quant_label == "999" for its
# RK table. We match that precedent exactly (explicit override) rather than
# the harness default, since the cutoffs 3.5/6.5 we are reusing were
# calibrated against the 999% table's convention. NN default ("99" for
# d in [10,19]) already agrees with 09_wp3_synthetic.R's own assertion, so
# no NN override is needed.
RK_TAB <- get_simul("RK", D, quant = "999")
NN_TAB <- get_simul("NN", D)
stopifnot(RK_TAB$quant_label == "999", NN_TAB$quant_label == "99")

# ---------------------------------------------------------------------------
# Data generator: gen_gaussian() from 09_wp3_synthetic.R, generalized to an
# equicorrelated within-cluster covariance parameterized by rho. rho = 0
# reproduces the original isotropic generator exactly.
# ---------------------------------------------------------------------------
gen_gaussian_rho <- function(seed, rho) {
  n = 200; d = 10; cont = 0.05
  cls_dis = 3; otl_dis = 2
  r_min = 0.7; r_max = 1.3

  mu1 = rep(3, d)
  mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = apply(rbind(mu1, mu2), 2, mean)

  n1 = round(n * (1 - cont) * 0.5)
  n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)

  noise_level = 0.01
  sigma = 1 / sqrt(qchisq(1 - noise_level, d))

  eqcorr <- function(scale2) scale2 * ((1 - rho) * diag(d) + rho * matrix(1, d, d))

  set.seed(seed)
  data1 = mvrnorm(n1, mu1, eqcorr((sigma * runif(1, r_min, r_max))^2))
  data2 = mvrnorm(n2, mu2, eqcorr((sigma * runif(1, r_min, r_max))^2))
  i = 0
  outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    r1 = sqrt(sum((temp - mu1)^2))
    r2 = sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) {
      outlier = rbind(outlier, temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

# ---------------------------------------------------------------------------
# Checkpointing
# ---------------------------------------------------------------------------
ROW_COLS <- c("rho", "rep", "seed", "method", "n", "n0", "cutoff",
              "TPR", "TNR", "BA", "F2", "t_score", "status", "timestamp")

rep_done <- function(rho, rep) {
  if (!file.exists(RAW_CSV)) return(FALSE)
  df <- tryCatch(read.csv(RAW_CSV, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  sub <- df[abs(df$rho - rho) < 1e-9 & df$rep == rep & df$status == "OK", ]
  length(unique(sub$method)) == length(OS_METHODS)
}

append_row <- function(row) {
  stopifnot(identical(names(row), ROW_COLS))
  df <- as.data.frame(row, stringsAsFactors = FALSE)
  if (!file.exists(RAW_CSV)) {
    write.csv(df, RAW_CSV, row.names = FALSE)
  } else {
    write.table(df, RAW_CSV, sep = ",", col.names = FALSE, row.names = FALSE,
                append = TRUE, qmethod = "double")
  }
}

make_row <- function(rho, rep, seed, method, n, n0, cutoff, v, t_score, status) {
  num_or_na <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else round(as.numeric(x), 4)
  list(rho = rho, rep = rep, seed = seed, method = method, n = n, n0 = n0,
       cutoff = cutoff,
       TPR = num_or_na(v[["TPR"]]), TNR = num_or_na(v[["TNR"]]),
       BA = num_or_na(v[["BA"]]), F2 = num_or_na(v[["F2"]]),
       t_score = round(t_score, 3), status = status,
       timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
BASE_SEED <- 4200L  # arbitrary, disjoint from 09_wp3_synthetic.R's 123/1234
run_start <- Sys.time()

BUDGET_HIT <- FALSE
for (rho in RHOS) {
  if (BUDGET_HIT) break
  for (rep in seq_len(N_REPS)) {
    if (rep_done(rho, rep)) {
      cat(sprintf("[rho=%.1f rep=%d/%d] already complete -- skip\n", rho, rep, N_REPS))
      next
    }
    if (as.numeric(difftime(Sys.time(), run_start, units = "secs")) > MAX_SECONDS) {
      cat(sprintf("[budget] %.0fs elapsed > MAX_SECONDS=%.0f -- stopping before next block, checkpoint is clean.\n",
                  as.numeric(difftime(Sys.time(), run_start, units = "secs")), MAX_SECONDS))
      BUDGET_HIT <- TRUE
      break
    }
    seed <- BASE_SEED + rep * 10L + which(RHOS == rho)
    dat <- gen_gaussian_rho(seed, rho)
    X <- dat$X; n <- dat$n; n0 <- dat$n0
    Y <- c(rep(1, n - n0), rep(0, n0))

    score_calls <- list(
      "RKCCD-OOS" = function() RKCCD_OOS(datax = X, simul = RK_TAB$simul, d = D, quant = RK_TAB$quant),
      "RKCCD-IOS" = function() RKCCD_IOS(datax = X, simul = RK_TAB$simul, d = D, quant = RK_TAB$quant, min.cls = 0),
      "UNCCD-OOS" = function() NNCCD_OOS(datax = X, simul = NN_TAB$simul, method = "descend", d = D),
      "UNCCD-IOS" = function() NNCCD_IOS(datax = X, simul = NN_TAB$simul, method = "descend", d = D, min.cls = 0)
    )

    for (m in OS_METHODS) {
      t0 <- Sys.time()
      res <- tryCatch(list(score = score_calls[[m]]()), error = function(e) e)
      t_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (inherits(res, "error")) {
        status <- paste0("ERROR: ", substr(gsub("[\r\n,]+", " ", conditionMessage(res)), 1, 160))
        append_row(make_row(rho, rep, seed, m, n, n0, CUTOFFS[[m]],
                            c(TPR = NA, TNR = NA, BA = NA, F2 = NA), t_score, status))
        cat(sprintf("[rho=%.1f rep=%d %-10s] %s\n", rho, rep, m, status))
        next
      }
      score <- res$score
      if (length(score) != n || any(is.na(score))) {
        status <- sprintf("ERROR: bad score vector (len %d vs n %d, NA %s)", length(score), n, any(is.na(score)))
        append_row(make_row(rho, rep, seed, m, n, n0, CUTOFFS[[m]],
                            c(TPR = NA, TNR = NA, BA = NA, F2 = NA), t_score, status))
        cat(sprintf("[rho=%.1f rep=%d %-10s] %s\n", rho, rep, m, status))
        next
      }
      v <- evaluate(Y, score, CUTOFFS[[m]])
      append_row(make_row(rho, rep, seed, m, n, n0, CUTOFFS[[m]], v, t_score, "OK"))
      cat(sprintf("[rho=%.1f rep=%d/%d %-10s] BA=%.3f F2=%.3f t=%.1fs (elapsed %.1f min)\n",
                  rho, rep, N_REPS, m, v["BA"], v["F2"], t_score,
                  as.numeric(difftime(Sys.time(), run_start, units = "mins"))))
    }
  }
}

cat(sprintf("\n---- sweep done in %.1f min; aggregating ----\n",
            as.numeric(difftime(Sys.time(), run_start, units = "mins"))))

# ---------------------------------------------------------------------------
# Aggregation: mean / sd / se over reps, per (rho, method)
# ---------------------------------------------------------------------------
raw <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
ok <- raw[raw$status == "OK", ]

agg <- do.call(rbind, lapply(split(ok, list(ok$rho, ok$method), drop = TRUE), function(g) {
  se <- function(x) if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  data.frame(
    rho = g$rho[1], method = g$method[1], n_reps = nrow(g),
    TPR_mean = mean(g$TPR), TNR_mean = mean(g$TNR),
    BA_mean = mean(g$BA), F2_mean = mean(g$F2),
    TPR_sd = sd(g$TPR), TNR_sd = sd(g$TNR), BA_sd = sd(g$BA), F2_sd = sd(g$F2),
    BA_se = se(g$BA), F2_se = se(g$F2),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$method, agg$rho), ]
rownames(agg) <- NULL
write.csv(agg, AGG_CSV, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", AGG_CSV, nrow(agg)))

cat("\nBA_mean / F2_mean by method x rho:\n")
for (m in OS_METHODS) {
  sub <- agg[agg$method == m, ]
  cat(sprintf("  %-10s ", m))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("rho=%.1f:BA=%.3f/F2=%.3f  ", sub$rho[i], sub$BA_mean[i], sub$F2_mean[i]))
  }
  cat("\n")
}

cat("\n20_wp_collinearity.R DONE\n")
