#!/usr/bin/env Rscript
# revision_experiments/tr1/86_wp8_outlier_types.R
#
# WP8 experiment 2 -- local / bridged / collective outliers (R3.7). See
# WP8_PROTOCOL.md for the pre-declared design; this file implements it.
#
# R3.7 argues the mutual-catch-graph rule (outliers lack connectivity to a
# cluster core) is tuned for GLOBAL outliers and may be restrictive for:
#   - local outliers embedded within a cluster's own extent,
#   - anomalies bridged to a cluster core through a chain of connecting
#     points,
#   - internally well-connected collective groups of anomalies near a
#     cluster boundary.
# Three generators, one per case, each replacing only the outlier placement
# in the same two-cluster base (n1=95, n2=94 at cont=0.05, same geometry
# constants as 85_wp8_boundary_fp.R / 55_wp2c_simulation_arm.R). n0 = 10
# throughout, so evaluate()/count_scores2 is safe (n0 > 0), unlike experiment
# 4 (88), which has no outliers at all.
#
# Settings: type in {local, bridge, collective} x d in {3, 10}, n = 200,
# n0 = 10. 6 settings, 100 reps, 9 methods (4 MCCD + 5 baselines) per rep.
#
# Single row per cell -- has_result() on the main CSV gates skip/restart
# directly (no done-marker needed, unlike 85).
#
# Usage:
#   Rscript 86_wp8_outlier_types.R --smoke
#     One rep of one setting (local, d=3), all 9 default methods, written to
#     results/tr1/wp8/smoke/86_wp8_outlier_types.csv.
#   Rscript 86_wp8_outlier_types.R [reps] [dims] [methods]

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

ALL_THRESHOLDS <- REAL_DATA_THRESHOLDS
DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
                      "LOF", "DBSCAN", "MST", "ODIN", "iForest")
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

BASE_SEED <- 8601L
CONT <- 0.05
N_NOMINAL <- 200L
CLS_DIS <- 3; R_MIN <- 0.7; R_MAX <- 1.3

RESDIR    <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV  <- file.path(RESDIR, "86_outlier_types.csv")
CLUST_CSV <- file.path(RESDIR, "86_outlier_types_clusters.csv")
SMOKE_DIR <- file.path(RESDIR, "smoke")
SMOKE_MAIN  <- file.path(SMOKE_DIR, "86_wp8_outlier_types.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "86_wp8_outlier_types_clusters.csv")

# ---------------------------------------------------------------------------
# Base two-cluster regular points, shared by all three generators. `avoid`,
# if given, is a list of (center, radius) pairs; candidate regular points
# falling within `radius` of any `center` are redrawn (up to 200 attempts,
# then kept anyway and the shortfall logged) -- used only by the "local"
# generator to carve out the density hole its outliers sit in.
# ---------------------------------------------------------------------------
draw_cluster <- function(n_k, d, mu_k, avoid = NULL) {
  scale <- runif(1, R_MIN, R_MAX)
  pts <- matrix(NA_real_, nrow = n_k, ncol = d)
  n_forced <- 0L
  for (i in seq_len(n_k)) {
    ok <- FALSE
    for (attempt in 1:200) {
      cand <- as.numeric(rpoisball.unit(1, d)) * scale + mu_k
      clash <- FALSE
      if (!is.null(avoid)) {
        for (a in avoid) if (sqrt(sum((cand - a$center)^2)) < a$radius) { clash <- TRUE; break }
      }
      if (!clash) { ok <- TRUE; break }
    }
    if (!ok) n_forced <- n_forced + 1L
    pts[i, ] <- cand
  }
  attr(pts, "n_forced") <- n_forced
  pts
}

# ---------------------------------------------------------------------------
# Three outlier-type generators (WP8_PROTOCOL.md experiment 2)
# ---------------------------------------------------------------------------
gen_local <- function(seed, n, d, cont) {
  cls_dis <- CLS_DIS
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  n_k1 <- ceiling(n0 / 2); n_k2 <- n0 - n_k1
  gaps <- vector("list", n0)
  for (i in seq_len(n0)) {
    host <- if (i <= n_k1) 1L else 2L
    mu_k <- if (host == 1L) mu1 else mu2
    g <- mu_k + as.numeric(rpoisball.unit(1, d)) * 0.4
    gaps[[i]] <- list(center = g, radius = 0.3, host = host)
  }
  avoid1 <- Filter(function(a) a$host == 1L, gaps)
  avoid2 <- Filter(function(a) a$host == 2L, gaps)
  data1 <- draw_cluster(n1, d, mu1, avoid1)
  data2 <- draw_cluster(n2, d, mu2, avoid2)
  outlier <- do.call(rbind, lapply(gaps, function(a) a$center))
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)))
}

gen_bridge <- function(seed, n, d, cont) {
  cls_dis <- CLS_DIS
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- draw_cluster(n1, d, mu1)
  data2 <- draw_cluster(n2, d, mu2)
  outlier <- do.call(rbind, lapply(seq_len(n0), function(i) {
    t_i <- i / (n0 + 1)
    mu1 + t_i * (mu2 - mu1) + rnorm(d, 0, 0.15)
  }))
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)))
}

gen_collective <- function(seed, n, d, cont) {
  cls_dis <- CLS_DIS
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- draw_cluster(n1, d, mu1)
  data2 <- draw_cluster(n2, d, mu2)
  u <- rnorm(d); u <- u / sqrt(sum(u^2))
  center <- mu1 + u * (R_MAX + 0.5)
  outlier <- matrix(rep(center, n0), nrow = n0, byrow = TRUE) + matrix(rnorm(n0 * d, 0, 0.15), nrow = n0)
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)))
}

GENERATORS <- list(local = gen_local, bridge = gen_bridge, collective = gen_collective)

build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(type = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$type, s$d)
  s
}

# ---------------------------------------------------------------------------
# One (setting, d, rep, method) cell
# ---------------------------------------------------------------------------
run_cell <- function(setting_id, type, d, rep_id, seed, method, main_csv, clust_csv) {
  keys <- c(setting_id = setting_id, d = d, rep = rep_id, method = method)
  if (isTRUE(has_result(main_csv, keys))) {
    cat(sprintf("[skip] %s d=%d rep=%d method=%s\n", setting_id, d, rep_id, method))
    return(invisible("skip"))
  }

  dat <- GENERATORS[[type]](seed, N_NOMINAL, d, CONT)
  X <- dat$X; n <- dat$n; n0 <- dat$n0; n_reg <- dat$n1 + dat$n2
  Y <- c(rep(1L, n_reg), rep(0L, n0))

  extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN) else list()
  out <- tryCatch({
    res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra))
    stopifnot(length(res$score) == n, !anyNA(res$score))
    m <- evaluate(Y, res$score, ALL_THRESHOLDS[[method]])
    # note = "-" (never "" or NA) on success. has_result() treats ANY NA
    # in a non-key column as an incomplete/partial row (its own
    # truncation guard) -- and an empty string doesn't dodge this: a CSV
    # column that is "" on EVERY row round-trips through read.csv's
    # automatic type.convert() as logical NA (an all-empty character
    # column has no type to infer, so it collapses to `logical NA`),
    # which is exactly the NA has_result is watching for. A non-empty
    # placeholder avoids the type.convert collapse.
    list(m = m, res = res, status = "ok", note = "-")
  }, error = function(e) list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
                              res = NULL, status = "error",
                              note = substr(conditionMessage(e), 1, 200)))

  row <- list(setting_id = setting_id, type = type, d = d, rep = rep_id, seed = seed,
              method = method,
              TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
              BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
              n_flagged = if (is.null(out$res)) NA_integer_ else sum(out$res$score >= ALL_THRESHOLDS[[method]]),
              status = out$status, note = out$note)
  append_result(main_csv, row)

  if (identical(out$status, "ok") && method %in% MCCD_METHODS && !is.null(out$res$cluster)) {
    for (i in seq_len(n_reg)) {
      append_result(clust_csv, list(
        setting_id = setting_id, d = d, rep = rep_id, seed = seed, method = method,
        row_index = i, true_cluster = dat$true_cluster[i],
        detected_cluster = if (is.na(out$res$cluster[i])) NA_integer_ else out$res$cluster[i]))
    }
  }

  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE %s d=%d rep=%d method=%s TPR=%.3f TNR=%.3f\n",
                setting_id, d, rep_id, method, out$m[["TPR"]], out$m[["TNR"]]))
  } else {
    cat(sprintf("CELL_FAIL %s d=%d rep=%d method=%s note=%s\n",
                setting_id, d, rep_id, method, out$note))
  }
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
args <- args[args != "--smoke"]

if (MODE_SMOKE) {
  dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
  st <- build_settings(3L)
  s <- st[st$type == "local", ][1, ]
  cat(sprintf("==== 86_wp8_outlier_types --smoke: %s, 1 rep, %d methods ====\n",
              s$setting_id, length(DEFAULT_METHODS)))
  seed <- BASE_SEED + 100000L * 1L + 9001L
  for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$type, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
  cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
  for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$type, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

N_REPS  <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 100L
DIMS    <- if (length(args) >= 2 && nzchar(args[2])) as.integer(strsplit(args[2], ",")[[1]]) else c(3L, 10L)
METHODS <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else DEFAULT_METHODS

SETTINGS <- build_settings(DIMS)
cat(sprintf("86_wp8_outlier_types: %d settings x %d reps x %d methods\n",
            nrow(SETTINGS), N_REPS, length(METHODS)))

for (s_i in seq_len(nrow(SETTINGS))) {
  s <- SETTINGS[s_i, ]
  for (rep_id in seq_len(N_REPS)) {
    seed <- BASE_SEED + 100000L * s_i + rep_id
    for (m in METHODS) run_cell(s$setting_id, s$type, s$d, rep_id, seed, m, MAIN_CSV, CLUST_CSV)
  }
}
cat("86_wp8_outlier_types: done.\n")
