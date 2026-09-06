#!/usr/bin/env Rscript
# revision_experiments/tr1/85_wp8_boundary_fp.R
#
# WP8 experiment 1 -- boundary false positives (R1.1). See WP8_PROTOCOL.md
# for the pre-declared design; this file implements it and must not diverge
# from it silently.
#
# R1.1 asks a yes/no question that deserves a number: mutual capture (each
# side of an edge must catch the other) is a stricter condition than one-sided
# coverage, so points near a cluster's own boundary are structurally more
# likely to be excluded from the mutual-catch-graph's majority component than
# points near its centre. This script measures the false-positive rate of
# REGULAR points as a function of their normalized distance from their own
# cluster's centre, in 10 bins of r/r_max, for the two-cluster uniform and
# Gaussian settings already used elsewhere in this revision
# (55_wp2c_simulation_arm.R's gen_uniform/gen_gaussian, reproduced here with
# the extra bookkeeping -- per-point cluster id and centre -- that the bin
# computation needs and the original drivers never returned).
#
# Settings: generator in {uniform, gaussian} x d in {3, 10}, n = 200,
# contamination = 0.05. 4 settings, 100 reps, 9 methods (4 MCCD + 5
# baselines) per rep -- see WP8_PROTOCOL.md's shared conventions.
#
# ATOMICITY. Ten bin-rows share one detector call. A companion "done" file
# (keys: setting_id, d, rep, method) is appended only after all 10 bin rows
# AND (for the 4 MCCD methods) the WP7 cluster rows for that cell have been
# written; has_result() against the done file, not the main file, is what
# gates skip/restart, so a rep is never half-recorded across a restart even
# though its result is not one physical row.
#
# Usage:
#   Rscript 85_wp8_boundary_fp.R --smoke
#     One rep of one setting (uniform, d=3), all 9 default methods, written
#     to results/tr1/wp8/smoke/85_wp8_boundary_fp.csv. Safe to run twice to
#     see the done-file skip fire on the second call.
#   Rscript 85_wp8_boundary_fp.R [reps] [dims] [methods]
#     reps:    integer, default 100
#     dims:    comma list, default "3,10"
#     methods: comma list, default the 9-method WP8 default (see below)

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

ALL_THRESHOLDS <- REAL_DATA_THRESHOLDS  # merged by wp0_mccd_methods.R's append
DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
                      "LOF", "DBSCAN", "MST", "ODIN", "iForest")
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

BASE_SEED <- 8501L
CONT <- 0.05
N_NOMINAL <- 200L

RESDIR    <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV  <- file.path(RESDIR, "85_boundary_fp.csv")
DONE_CSV  <- file.path(RESDIR, "85_boundary_fp_done.csv")
CLUST_CSV <- file.path(RESDIR, "85_boundary_fp_clusters.csv")
SMOKE_DIR <- file.path(RESDIR, "smoke")
SMOKE_MAIN  <- file.path(SMOKE_DIR, "85_wp8_boundary_fp.csv")
SMOKE_DONE  <- file.path(SMOKE_DIR, "85_wp8_boundary_fp_done.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "85_wp8_boundary_fp_clusters.csv")

# ---------------------------------------------------------------------------
# Generators (WP8_PROTOCOL.md experiment 1). Identical geometry constants to
# 55_wp2c_simulation_arm.R's gen_uniform/gen_gaussian; extra return values
# (per-point cluster id, cluster centres) added for the boundary-bin
# computation, which those drivers never needed.
# ---------------------------------------------------------------------------
gen_uniform_boundary <- function(seed, n, d, cont) {
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- (mu1 + mu2) / 2
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis && sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1
    }
  }
  rownames(outlier) <- NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0,
       n1 = n1, n2 = n2, true_cluster = c(rep(1L, n1), rep(2L, n2), rep(NA_integer_, n0)),
       mu1 = mu1, mu2 = mu2)
}

gen_gaussian_boundary <- function(seed, n, d, cont) {
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- (mu1 + mu2) / 2
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  noise_level <- 0.01
  sigma <- 1 / sqrt(qchisq(1 - noise_level, d))
  set.seed(seed)
  data1 <- mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 <- mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis && sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1
    }
  }
  rownames(outlier) <- NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0,
       n1 = n1, n2 = n2, true_cluster = c(rep(1L, n1), rep(2L, n2), rep(NA_integer_, n0)),
       mu1 = mu1, mu2 = mu2)
}

GENERATORS <- list(uniform = gen_uniform_boundary, gaussian = gen_gaussian_boundary)

# ---------------------------------------------------------------------------
# Settings table
# ---------------------------------------------------------------------------
build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(generator = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$generator, s$d)
  s
}

# ---------------------------------------------------------------------------
# Per-point normalized-radius bin
# ---------------------------------------------------------------------------
compute_bins <- function(X, true_cluster, mu1, mu2, n_reg) {
  bin <- rep(NA_integer_, n_reg)
  for (k in c(1L, 2L)) {
    idx <- which(true_cluster[seq_len(n_reg)] == k)
    if (!length(idx)) next
    mu_k <- if (k == 1L) mu1 else mu2
    r <- sqrt(rowSums(sweep(X[idx, , drop = FALSE], 2, mu_k)^2))
    r_max <- max(r)
    rn <- if (r_max > 0) r / r_max else rep(0, length(r))
    bin[idx] <- pmin(pmax(ceiling(10 * rn), 1L), 10L)
  }
  bin
}

# ---------------------------------------------------------------------------
# One (setting, d, rep, method) cell
# ---------------------------------------------------------------------------
run_cell <- function(setting_id, generator, d, rep_id, seed, method,
                      main_csv, done_csv, clust_csv) {
  keys <- c(setting_id = setting_id, d = d, rep = rep_id, method = method)
  if (isTRUE(has_result(done_csv, keys))) {
    cat(sprintf("[skip] %s d=%d rep=%d method=%s\n", setting_id, d, rep_id, method))
    return(invisible("skip"))
  }

  dat <- GENERATORS[[generator]](seed, N_NOMINAL, d, CONT)
  X <- dat$X; n <- dat$n; n0 <- dat$n0; n_reg <- dat$n1 + dat$n2
  Y <- c(rep(1L, n_reg), rep(0L, n0))

  extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN) else list()
  res <- tryCatch(do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra)),
                   error = function(e) e)
  if (inherits(res, "error")) {
    cat(sprintf("CELL_FAIL %s d=%d rep=%d method=%s note=%s\n",
                setting_id, d, rep_id, method, substr(conditionMessage(res), 1, 150)))
    return(invisible("error"))
  }

  bin <- compute_bins(X[seq_len(n_reg), , drop = FALSE], dat$true_cluster, dat$mu1, dat$mu2, n_reg)
  flagged <- res$score[seq_len(n_reg)] >= ALL_THRESHOLDS[[method]]

  for (b in 1:10) {
    idx_b <- which(bin == b)
    row <- list(setting_id = setting_id, generator = generator, d = d, rep = rep_id,
                seed = seed, method = method, bin = b,
                n_regular_in_bin = length(idx_b),
                n_flagged_in_bin = sum(flagged[idx_b]))
    append_result(main_csv, row)
  }

  if (method %in% MCCD_METHODS && !is.null(res$cluster)) {
    for (i in seq_len(n_reg)) {
      append_result(clust_csv, list(
        setting_id = setting_id, d = d, rep = rep_id, seed = seed, method = method,
        row_index = i, true_cluster = dat$true_cluster[i],
        detected_cluster = if (is.na(res$cluster[i])) NA_integer_ else res$cluster[i]))
    }
  }

  append_result(done_csv, as.list(keys))
  cat(sprintf("CELL_DONE %s d=%d rep=%d method=%s n_flagged_total=%d\n",
              setting_id, d, rep_id, method, sum(flagged)))
  flush.console()
  invisible("ok")
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
args <- args[args != "--smoke"]

if (MODE_SMOKE) {
  dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
  st <- build_settings(3L)  # uniform_d3 only (first generator)
  s <- st[st$generator == "uniform", ][1, ]
  cat(sprintf("==== 85_wp8_boundary_fp --smoke: %s, 1 rep, %d methods ====\n",
              s$setting_id, length(DEFAULT_METHODS)))
  seed <- BASE_SEED + 100000L * 1L + 9001L
  for (m in DEFAULT_METHODS) {
    run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_DONE, SMOKE_CLUST)
  }
  cat("---- re-invoking the same cells to demonstrate done-file skip ----\n")
  for (m in DEFAULT_METHODS) {
    run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_DONE, SMOKE_CLUST)
  }
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

N_REPS  <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 100L
DIMS    <- if (length(args) >= 2 && nzchar(args[2])) as.integer(strsplit(args[2], ",")[[1]]) else c(3L, 10L)
METHODS <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else DEFAULT_METHODS

SETTINGS <- build_settings(DIMS)
cat(sprintf("85_wp8_boundary_fp: %d settings x %d reps x %d methods\n",
            nrow(SETTINGS), N_REPS, length(METHODS)))

for (s_i in seq_len(nrow(SETTINGS))) {
  s <- SETTINGS[s_i, ]
  for (rep_id in seq_len(N_REPS)) {
    seed <- BASE_SEED + 100000L * s_i + rep_id
    for (m in METHODS) {
      run_cell(s$setting_id, s$generator, s$d, rep_id, seed, m, MAIN_CSV, DONE_CSV, CLUST_CSV)
    }
  }
}
cat("85_wp8_boundary_fp: done.\n")
