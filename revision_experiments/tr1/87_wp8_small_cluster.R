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
# Settings: m in {2,4,6,8,10,12,15,20} x d in {3,10}. 16 settings. Reps
# differ BY METHOD FAMILY (2026-09-05 cost decision, see the dated note
# below): RK methods (U-MCCD, SU-MCCD) at 50 reps, NND methods (UN-MCCD,
# SUN-MCCD) at 100 reps -- U/SU-MCCD cost ~11 s/cell at d=3 and are 96% of
# the experiment's total cost. Single row per cell -- has_result() on the
# main CSV gates skip/restart directly.
#
# Usage:
#   Rscript 87_wp8_small_cluster.R --smoke
#     One rep of one setting (m=10, d=3), all 4 methods, written to
#     results/tr1/wp8/smoke/87_wp8_small_cluster.csv.
#   Rscript 87_wp8_small_cluster.R --summarize
#     Mean/SE of frac_flagged and n_clusters per (m, d, method), plus the
#     fraction of reps with n_clusters == 3. Writes
#     results/tr1/wp8/87_small_cluster_summary.csv.
#   Rscript 87_wp8_small_cluster.R [settings] [reps] [budget] [msel]
#     settings: comma list of setting_id (from the fixed canonical grid
#               m x d in {3,10}, e.g. "m10_d3,m20_d10"), or "ALL" (default).
#               The canonical grid never changes with this selection, so a
#               setting's seed never depends on which subset of settings a
#               chunk happens to cover.
#     reps:     if given, overrides the per-method-family rep count
#               UNIFORMLY for all 4 methods (useful for quick manual
#               chunking/testing); if omitted, RK methods run 50 reps and
#               NND methods run 100 reps (the 2026-09-05 cost decision).
#     budget:   stop starting new reps after this many seconds (default 480)
#     msel:     comma list of m values (subset of SIZES) to further
#               restrict `settings` -- e.g. "2,4" runs only those sizes
#               across whichever d's `settings` already selected.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
RK_METHODS  <- c("U-MCCD", "SU-MCCD")
NND_METHODS <- c("UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

# Cost decision (2026-09-05, WP8_VERIFICATION.md item 87.7): RK methods at
# 50 reps, NND methods at 100 reps -- see the dated note below.
REPS_BY_FAMILY <- c("U-MCCD" = 50L, "SU-MCCD" = 50L, "UN-MCCD" = 100L, "SUN-MCCD" = 100L)

BASE_SEED <- 8701L
N_TOTAL <- 200L
SIZES <- c(2L, 4L, 6L, 8L, 10L, 12L, 15L, 20L)
FULL_DIMS <- c(3L, 10L)
CLS_DIS <- 3; R_MIN <- 0.7; R_MAX <- 1.3

RESDIR   <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV <- file.path(RESDIR, "87_small_cluster.csv")
CLUST_CSV <- file.path(RESDIR, "87_small_cluster_clusters.csv")
SUMM_CSV <- file.path(RESDIR, "87_small_cluster_summary.csv")
SMOKE_DIR  <- file.path(RESDIR, "smoke")
SMOKE_MAIN <- file.path(SMOKE_DIR, "87_wp8_small_cluster.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "87_wp8_small_cluster_clusters.csv")

# ---------------------------------------------------------------------------
# Generator: two background clusters plus a third genuine cluster of size m.
# mu1 = rep(3,d), mu2 = mu1 + cls_dis*e1, mu3 = mu1 + cls_dis*e2 (first two
# standard basis vectors -- needs d >= 2, satisfied at d = 3 and d = 10).
# Background clusters each get their own independent radius-jitter draw
# (two separate runif(1, R_MIN, R_MAX) calls).
#
# FIX (2026-09-05, WP8_VERIFICATION.md item 87.1): the original generator
# gave cluster 3 the SAME radius draw as clusters 1-2, so shrinking m also
# shrank cluster 3's DENSITY, not just its size (measured NN-spacing ratio
# vs background at d=3: 4.42 at m=2, 2.26 at m=10, 1.71 at m=20 -- an
# increasingly sparse, not just increasingly small, cluster). Cluster 3's
# radius is now scaled by (m/n1)^(1/d) relative to an independent
# runif(1,R_MIN,R_MAX) draw, so that its point DENSITY (points per unit
# volume) matches cluster 1's at every m (volume ~ radius^d, so radius ~
# (m/n1)^(1/d) holds density m/volume constant relative to cluster 1's
# n1-point cluster).
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
  data3 <- rpoisball.unit(m,  d) * (runif(1, R_MIN, R_MAX) * (m / n1)^(1 / d)) +
    matrix(rep(mu3, m), ncol = d, byrow = TRUE)
  list(X = rbind(data1, data2, data3), n = n1 + n2 + m, n1 = n1, n2 = n2, m = m,
       true_cluster = c(rep(1L, n1), rep(2L, n2), rep(3L, m)))
}

build_settings <- function(dims) {
  s <- do.call(rbind, lapply(SIZES, function(m)
    do.call(rbind, lapply(dims, function(d) data.frame(m = m, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("m%d_d%d", s$m, s$d)
  s
}
CANON <- build_settings(FULL_DIMS)

resolve_settings <- function(sel) {
  if (is.null(sel) || (length(sel) == 1 && toupper(sel) == "ALL")) return(CANON)
  out <- CANON[match(sel, CANON$setting_id), ]
  if (anyNA(out$setting_id)) stop(sprintf("unknown setting_id(s): %s",
                                          paste(sel[is.na(match(sel, CANON$setting_id))], collapse = ", ")))
  out
}

# ---------------------------------------------------------------------------
# One (m, d, rep, method) cell
# ---------------------------------------------------------------------------
run_cell <- function(m, d, rep_id, seed, method, out_csv, clust_csv) {
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
    list(frac_flagged = mean(res$score[third_idx] == 1), n_clusters = ncls,
         n_unassigned_third = sum(is.na(res$cluster[third_idx])),
         singleton_lost = if (!is.null(res$singleton_lost_rows)) as.integer(res$singleton_lost_rows) else -1L,
         res = res, status = "ok", note = "-")
  }, error = function(e) list(frac_flagged = NA_real_, n_clusters = NA_integer_,
                              n_unassigned_third = NA_integer_, singleton_lost = NA_integer_,
                              res = NULL, status = "error", note = substr(conditionMessage(e), 1, 200)))

  row <- list(m = m, d = d, rep = rep_id, seed = seed, method = method,
              frac_flagged = out$frac_flagged, n_clusters = out$n_clusters,
              n_unassigned_third = out$n_unassigned_third, singleton_lost = out$singleton_lost,
              status = out$status, note = out$note)
  append_result(out_csv, row)

  if (identical(out$status, "ok") && !is.null(out$res$cluster)) {
    # Block write: all n cluster rows for this cell in ONE append_result
    # call (WP7 file; new for 87 -- the original protocol declared "no WP7
    # cluster file for this experiment", superseded by
    # WP8_VERIFICATION.md item 87.4, see the dated note below).
    idx <- seq_len(n)
    dc <- out$res$cluster[idx]
    append_result(clust_csv, list(
      m = rep(m, n), d = rep(d, n), rep = rep(rep_id, n), seed = rep(seed, n),
      method = rep(method, n), row_index = idx, true_cluster = dat$true_cluster[idx],
      detected_cluster = ifelse(is.na(dc), NA_integer_, dc)))
  }

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
# --summarize: mean/SE of frac_flagged and n_clusters per (m, d, method),
# plus the fraction of reps with n_clusters == 3. status=="ok" only,
# de-duplicated on (m, d, rep, method) keeping the LAST occurrence.
# ---------------------------------------------------------------------------
do_summarize <- function() {
  if (!file.exists(MAIN_CSV)) { cat("no results yet\n"); return(invisible(NULL)) }
  df_all <- read.csv(MAIN_CSV, stringsAsFactors = FALSE)
  key <- paste(df_all$m, df_all$d, df_all$rep, df_all$method, sep = "\r")
  df_all <- df_all[!duplicated(key, fromLast = TRUE), ]
  df <- df_all[df_all$status == "ok", ]

  rows <- list()
  for (mm in unique(df_all$m)) for (dd in unique(df_all$d[df_all$m == mm])) for (meth in DEFAULT_METHODS) {
    sub <- df[df$m == mm & df$d == dd & df$method == meth, ]
    n_err <- sum(df_all$m == mm & df_all$d == dd & df_all$method == meth & df_all$status != "ok")
    if (!nrow(sub)) next
    rows[[length(rows) + 1L]] <- data.frame(
      m = mm, d = dd, method = meth, n_reps = nrow(sub),
      frac_flagged_mean = mean(sub$frac_flagged), frac_flagged_se = sd(sub$frac_flagged) / sqrt(nrow(sub)),
      n_clusters_mean = mean(sub$n_clusters), n_clusters_se = sd(sub$n_clusters) / sqrt(nrow(sub)),
      frac_reps_ncls3 = mean(sub$n_clusters == 3),
      n_unassigned_third_mean = mean(sub$n_unassigned_third),
      n_error = n_err, stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  S <- S[order(S$method, S$d, S$m), ]
  write.csv(S, SUMM_CSV, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMM_CSV, nrow(S)))
  print(S, row.names = FALSE)
  invisible(S)
}

# ---------------------------------------------------------------------------
# --run: full (subset) grid, checkpointed by budget. Reps are PER METHOD
# (REPS_BY_FAMILY) unless overridden uniformly by the `reps` argument.
# ---------------------------------------------------------------------------
do_run <- function(sel, msel, reps_override, budget, methods) {
  st <- resolve_settings(sel)
  if (!is.null(msel)) st <- st[st$m %in% as.integer(msel), , drop = FALSE]
  stopifnot(nrow(st) > 0)
  cat(sprintf("87_wp8_small_cluster: %d settings x %d methods (budget %.0fs)\n",
              nrow(st), length(methods), budget))
  t_start <- Sys.time()
  for (m_id in methods) {
    n_reps <- if (!is.null(reps_override)) reps_override else REPS_BY_FAMILY[[m_id]]
    for (row in seq_len(nrow(st))) {
      s <- st[row, ]
      s_i <- match(s$setting_id, CANON$setting_id)
      for (rep_id in seq_len(n_reps)) {
        if (as.numeric(difftime(Sys.time(), t_start, units = "secs")) > budget) {
          cat(sprintf("[budget %gs exhausted -- rerun to continue] next up: %s x %s rep %d\n",
                      budget, s$setting_id, m_id, rep_id))
          return(invisible("budget_stop"))
        }
        seed <- BASE_SEED + 100000L * s_i + rep_id
        run_cell(s$m, s$d, rep_id, seed, m_id, MAIN_CSV, CLUST_CSV)
      }
    }
  }
  cat("87_wp8_small_cluster: run complete.\n")
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
MODE_SUMM  <- "--summarize" %in% args
args <- args[!args %in% c("--smoke", "--summarize")]

if (!isTRUE(getOption("wp8.no_main", FALSE))) {
  if (MODE_SMOKE) {
    dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("==== 87_wp8_small_cluster --smoke: m=10, d=3, 1 rep, %d methods ====\n",
                length(DEFAULT_METHODS)))
    s_i <- match("m10_d3", CANON$setting_id)
    seed <- BASE_SEED + 100000L * s_i + 9001L
    for (mth in DEFAULT_METHODS) run_cell(10L, 3L, 9001L, seed, mth, SMOKE_MAIN, SMOKE_CLUST)
    cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
    for (mth in DEFAULT_METHODS) run_cell(10L, 3L, 9001L, seed, mth, SMOKE_MAIN, SMOKE_CLUST)
    cat("ALL_CELLS_COMPLETE\n")
    quit(save = "no", status = 0)
  } else if (MODE_SUMM) {
    do_summarize()
    quit(save = "no", status = 0)
  } else {
    sel     <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
    reps_ov <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else NULL
    budget  <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 480
    msel    <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else NULL
    do_run(sel, msel, reps_ov, budget, DEFAULT_METHODS)
  }
}
