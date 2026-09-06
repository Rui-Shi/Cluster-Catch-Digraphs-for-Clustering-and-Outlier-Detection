#!/usr/bin/env Rscript
# revision_experiments/tr1/87_wp8_small_cluster.R
#
# WP8 experiment 3 -- a small legitimate cluster below S_min (R3.7). See
# WP8_PROTOCOL.md for the pre-declared design.
#
# R3.7: "a small but legitimate cluster may be classified as anomalous if
# its size is below S_min." Two background clusters plus a third genuine
# cluster whose size m is swept through {2,4,6,8,10,12,15,20} out of n=200,
# with S_min=0.05 (round(0.05*200)=10) as the size at which the flip is
# expected. No contamination outliers -- the only question is what happens
# to cluster 3 itself.
#
# Methods: SU-MCCD, SUN-MCCD (min.cls = S_min = 0.05, the focal methods --
# they are the only two whose wrapper takes a min.cls argument at all) and
# U-MCCD, UN-MCCD as controls (called identically at every m, since nothing
# in their call changes with m).
#
# Settings: m in {2,4,6,8,10,12,15,20} x d in {3,10}. 16 settings, 100 reps,
# 4 methods per rep. Single row per cell -- has_result() on the main CSV
# gates skip/restart directly.
#
# Usage:
#   Rscript 87_wp8_small_cluster.R --smoke
#     One rep of one setting (m=10, d=3), all 4 methods, written to
#     results/tr1/wp8/smoke/87_wp8_small_cluster.csv.
#   Rscript 87_wp8_small_cluster.R [reps] [dims] [methods]

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
S_MIN <- 0.05

BASE_SEED <- 8701L
N_TOTAL <- 200L
SIZES <- c(2L, 4L, 6L, 8L, 10L, 12L, 15L, 20L)
CLS_DIS <- 3; R_MIN <- 0.7; R_MAX <- 1.3

RESDIR   <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV <- file.path(RESDIR, "87_small_cluster.csv")
SMOKE_DIR  <- file.path(RESDIR, "smoke")
SMOKE_MAIN <- file.path(SMOKE_DIR, "87_wp8_small_cluster.csv")

# ---------------------------------------------------------------------------
# Generator: two background clusters plus a third genuine cluster of size m.
# mu1 = rep(3,d), mu2 = mu1 + cls_dis*e1, mu3 = mu1 + cls_dis*e2 (first two
# standard basis vectors -- needs d >= 2, satisfied at d = 3 and d = 10).
# Each cluster gets its own independent radius-jitter draw (three separate
# runif(1, R_MIN, R_MAX) calls), matching the "one shared jitter per cluster"
# idiom of every other WP8/WP2 generator.
# ---------------------------------------------------------------------------
gen_small_cluster <- function(seed, n_total, d, m) {
  stopifnot(d >= 2)
  e1 <- c(1, rep(0, d - 1)); e2 <- c(0, 1, rep(0, d - 2))
  mu1 <- rep(3, d); mu2 <- mu1 + CLS_DIS * e1; mu3 <- mu1 + CLS_DIS * e2
  n_bg <- n_total - m
  n1 <- floor(n_bg / 2); n2 <- n_bg - n1
  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, R_MIN, R_MAX) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, R_MIN, R_MAX) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  data3 <- rpoisball.unit(m,  d) * runif(1, R_MIN, R_MAX) + matrix(rep(mu3, m),  ncol = d, byrow = TRUE)
  list(X = rbind(data1, data2, data3), n = n1 + n2 + m, n1 = n1, n2 = n2, m = m)
}

# ---------------------------------------------------------------------------
# One (m, d, rep, method) cell
# ---------------------------------------------------------------------------
run_cell <- function(m, d, rep_id, seed, method, out_csv) {
  keys <- c(m = m, d = d, rep = rep_id, method = method)
  if (isTRUE(has_result(out_csv, keys))) {
    cat(sprintf("[skip] m=%d d=%d rep=%d method=%s\n", m, d, rep_id, method))
    return(invisible("skip"))
  }

  dat <- gen_small_cluster(seed, N_TOTAL, d, m)
  X <- dat$X; n <- dat$n
  third_idx <- (dat$n1 + dat$n2 + 1):n

  extra <- if (method %in% MIN_CLS_METHODS) list(min.cls = S_MIN) else list()
  out <- tryCatch({
    res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = NULL), extra))
    stopifnot(length(res$score) == n, !anyNA(res$score))
    ncls <- length(unique(res$cluster[!is.na(res$cluster)]))
    # note = "-" (never "" or NA) on success: has_result() treats ANY NA
    # in a non-key column as an incomplete row, and an all-"" column
    # round-trips through read.csv's type.convert() as logical NA (see
    # 86_wp8_outlier_types.R for the full explanation), so a non-empty
    # placeholder is required, not just a non-NA one.
    list(frac_flagged = mean(res$score[third_idx] == 1), n_clusters = ncls, status = "ok",
         note = "-")
  }, error = function(e) list(frac_flagged = NA_real_, n_clusters = NA_integer_,
                              status = "error", note = substr(conditionMessage(e), 1, 200)))

  row <- list(m = m, d = d, rep = rep_id, seed = seed, method = method,
              frac_flagged = out$frac_flagged, n_clusters = out$n_clusters,
              status = out$status, note = out$note)
  append_result(out_csv, row)

  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE m=%d d=%d rep=%d method=%s frac_flagged=%.3f n_clusters=%d\n",
                m, d, rep_id, method, out$frac_flagged, out$n_clusters))
  } else {
    cat(sprintf("CELL_FAIL m=%d d=%d rep=%d method=%s note=%s\n", m, d, rep_id, method, out$note))
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
  cat(sprintf("==== 87_wp8_small_cluster --smoke: m=10, d=3, 1 rep, %d methods ====\n",
              length(DEFAULT_METHODS)))
  seed <- BASE_SEED + 100000L * 1L + 9001L
  for (mth in DEFAULT_METHODS) run_cell(10L, 3L, 9001L, seed, mth, SMOKE_MAIN)
  cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
  for (mth in DEFAULT_METHODS) run_cell(10L, 3L, 9001L, seed, mth, SMOKE_MAIN)
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

N_REPS  <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 100L
DIMS    <- if (length(args) >= 2 && nzchar(args[2])) as.integer(strsplit(args[2], ",")[[1]]) else c(3L, 10L)
METHODS <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else DEFAULT_METHODS

SETTINGS <- do.call(rbind, lapply(SIZES, function(m)
  do.call(rbind, lapply(DIMS, function(d) data.frame(m = m, d = d)))))
cat(sprintf("87_wp8_small_cluster: %d settings x %d reps x %d methods\n",
            nrow(SETTINGS), N_REPS, length(METHODS)))

for (s_i in seq_len(nrow(SETTINGS))) {
  s <- SETTINGS[s_i, ]
  for (rep_id in seq_len(N_REPS)) {
    seed <- BASE_SEED + 100000L * s_i + rep_id
    for (mth in METHODS) run_cell(s$m, s$d, rep_id, seed, mth, MAIN_CSV)
  }
}
cat("87_wp8_small_cluster: done.\n")
