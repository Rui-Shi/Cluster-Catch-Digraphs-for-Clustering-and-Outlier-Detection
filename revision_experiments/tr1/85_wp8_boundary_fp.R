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
# cluster's centre, in 10 bins, for the two-cluster uniform and Gaussian
# settings already used elsewhere in this revision
# (55_wp2c_simulation_arm.R's gen_uniform/gen_gaussian, reproduced here with
# the extra bookkeeping -- per-point cluster id and centre -- that the bin
# computation needs and the original drivers never returned).
#
# Settings: generator in {uniform, gaussian} x d in {3, 10}, n = 200,
# contamination = 0.05. 4 settings, 100 reps, 9 methods (4 MCCD + 5
# baselines) per rep -- see WP8_PROTOCOL.md's shared conventions.
#
# ATOMICITY. 20 bin-rows (10 equal-width + 10 equal-mass, `bin_type` column;
# see the 2026-09-05 dated note in WP8_PROTOCOL.md) and, for the 4 MCCD
# methods, up to 189 cluster rows share one detector call. Both are written
# as ONE block each (append_result() accepts a list of equal-length vectors
# and builds a multi-row data.frame from it -- 6.4 ms per 189-row block vs
# 773 ms row by row). A companion "done" file (keys: setting_id, d, rep,
# method) is appended only after the bin-row block AND the cluster-row block
# (when applicable) have both landed; has_result() against the done file,
# not the main file, is what gates skip/restart, so a rep is never
# half-recorded across a restart even though its result is not one physical
# row.
#
# Usage:
#   Rscript 85_wp8_boundary_fp.R --smoke
#     One rep of one setting (uniform, d=3), all 9 default methods, written
#     to results/tr1/wp8/smoke/85_wp8_boundary_fp.csv. Safe to run twice to
#     see the done-file skip fire on the second call.
#   Rscript 85_wp8_boundary_fp.R --summarize
#     Per (setting_id, d, method, bin_type, bin): sum(n_flagged)/sum(n_regular)
#     as the point estimate (never a mean of per-rep ratios), plus mean/SE of
#     per-rep ratios over non-empty reps. Writes
#     results/tr1/wp8/85_boundary_fp_summary.csv.
#   Rscript 85_wp8_boundary_fp.R [settings] [reps] [budget] [methods]
#     settings: comma list of setting_id (from the fixed canonical grid
#               generator x d in {3,10}, e.g. "uniform_d3,gaussian_d10"), or
#               "ALL" (default). The canonical grid never changes with this
#               selection, so a setting's seed never depends on which subset
#               of settings a chunk happens to cover.
#     reps:     target replicate count (default 100)
#     budget:   stop starting new reps after this many seconds (default 480)
#     methods:  comma list, default the 9-method WP8 default (see below)

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
SUMM_CSV  <- file.path(RESDIR, "85_boundary_fp_summary.csv")
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
# Settings table -- fixed canonical grid (generator x d in {3,10}, nothing
# else declared). CANON never changes with how a run is chunked, which is
# what fixes the seed bug: s_i is always this table's row index, never the
# row index of a re-built, possibly-DIMS-restricted table.
# ---------------------------------------------------------------------------
build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(generator = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$generator, s$d)
  s
}
FULL_DIMS <- c(3L, 10L)
CANON <- build_settings(FULL_DIMS)

resolve_settings <- function(sel) {
  if (is.null(sel) || (length(sel) == 1 && toupper(sel) == "ALL")) return(CANON)
  out <- CANON[match(sel, CANON$setting_id), ]
  if (anyNA(out$setting_id)) stop(sprintf("unknown setting_id(s): %s",
                                          paste(sel[is.na(match(sel, CANON$setting_id))], collapse = ", ")))
  out
}

# ---------------------------------------------------------------------------
# Per-point normalized-radius bins. Two schemes (2026-09-05 dated note):
#   "width" -- ceiling(10 * r/r_max), the original scheme. At d=10, 66% of
#              regular points land in bin 10 and bins 1-5 hold ~0.1%
#              combined -- no interior to contrast. Also ill-defined for the
#              Gaussian arm, whose support has no hard edge (r_max is just
#              the largest OBSERVED radius, noisy from rep to rep).
#   "mass"  -- ceiling(10 * rank(r, ties.method="first") / length(r)), an
#              equal-COUNT scheme computed per cluster (same pooling as
#              "width": both clusters' regular points land in the same 10
#              bins). Immune to both defects: every bin gets ~n_k/10 points
#              by construction, and there is no r_max to be unbounded.
# Returns a list(width = <bin vector>, mass = <bin vector>), same length and
# indexing as the input.
# ---------------------------------------------------------------------------
compute_bins <- function(X, true_cluster, mu1, mu2, n_reg) {
  bin_width <- rep(NA_integer_, n_reg)
  bin_mass  <- rep(NA_integer_, n_reg)
  for (k in c(1L, 2L)) {
    idx <- which(true_cluster[seq_len(n_reg)] == k)
    if (!length(idx)) next
    mu_k <- if (k == 1L) mu1 else mu2
    r <- sqrt(rowSums(sweep(X[idx, , drop = FALSE], 2, mu_k)^2))
    r_max <- max(r)
    rn <- if (r_max > 0) r / r_max else rep(0, length(r))
    bin_width[idx] <- pmin(pmax(ceiling(10 * rn), 1L), 10L)
    bin_mass[idx]  <- pmin(pmax(ceiling(10 * rank(r, ties.method = "first") / length(r)), 1L), 10L)
  }
  list(width = bin_width, mass = bin_mass)
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

  out <- tryCatch({
    extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN) else list()
    Y <- c(rep(1L, n_reg), rep(0L, n0))
    res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra))
    stopifnot(length(res$score) == dat$n, !anyNA(res$score))
    list(res = res, status = "ok", note = "-")
  }, error = function(e) list(res = NULL, status = "error",
                              note = substr(conditionMessage(e), 1, 200)))

  if (identical(out$status, "error")) {
    cat(sprintf("CELL_FAIL %s d=%d rep=%d method=%s note=%s\n",
                setting_id, d, rep_id, method, out$note))
    flush.console()
    return(invisible("error"))
  }
  res <- out$res

  bins <- compute_bins(X[seq_len(n_reg), , drop = FALSE], dat$true_cluster, dat$mu1, dat$mu2, n_reg)
  flagged <- res$score[seq_len(n_reg)] >= ALL_THRESHOLDS[[method]]

  # Block write: 20 bin-rows (10 width + 10 mass) in ONE append_result call.
  n_reg_width <- sapply(1:10, function(b) sum(bins$width == b, na.rm = TRUE))
  n_flg_width <- sapply(1:10, function(b) sum(flagged[which(bins$width == b)]))
  n_reg_mass  <- sapply(1:10, function(b) sum(bins$mass == b, na.rm = TRUE))
  n_flg_mass  <- sapply(1:10, function(b) sum(flagged[which(bins$mass == b)]))
  append_result(main_csv, list(
    setting_id = rep(setting_id, 20), generator = rep(generator, 20), d = rep(d, 20),
    rep = rep(rep_id, 20), seed = rep(seed, 20), method = rep(method, 20),
    bin_type = rep(c("width", "mass"), each = 10), bin = rep(1:10, 2),
    n_regular_in_bin = c(n_reg_width, n_reg_mass),
    n_flagged_in_bin = c(n_flg_width, n_flg_mass)))

  # Block write: up to n_reg cluster rows (MCCD methods only) in ONE call.
  # unassigned_rows (from mccd_translate(), via wp0_mccd_methods.R -- not
  # edited here) has no baseline-method analogue; sentinel -1L (never
  # NA_integer_) is used for baselines so has_result()'s NA-based
  # partial-row guard is never tripped by a legitimately complete row (see
  # the NA-convention note in WP8_PROTOCOL.md).
  unassigned_rows <- if (!is.null(res$unassigned_rows)) as.integer(res$unassigned_rows) else -1L
  if (method %in% MCCD_METHODS && !is.null(res$cluster)) {
    idx <- seq_len(n_reg)
    dc <- res$cluster[idx]
    append_result(clust_csv, list(
      setting_id = rep(setting_id, n_reg), d = rep(d, n_reg), rep = rep(rep_id, n_reg),
      seed = rep(seed, n_reg), method = rep(method, n_reg), row_index = idx,
      true_cluster = dat$true_cluster[idx],
      detected_cluster = ifelse(is.na(dc), NA_integer_, dc)))
  }

  append_result(done_csv, c(as.list(keys), list(unassigned_rows = unassigned_rows)))
  cat(sprintf("CELL_DONE %s d=%d rep=%d method=%s n_flagged_total=%d\n",
              setting_id, d, rep_id, method, sum(flagged)))
  flush.console()
  invisible("ok")
}

# ---------------------------------------------------------------------------
# --summarize: per (setting_id, d, method, bin_type, bin), the point
# estimate is sum(n_flagged)/sum(n_regular) over reps (never a mean of
# per-rep ratios, which would over-weight reps with few points in that bin),
# plus mean and SE of the PER-REP ratio over reps with a non-empty bin and
# the count of such reps. De-duplicated on (setting_id,d,rep,method) keeping
# the last occurrence; only cells with a done marker are included (a bin row
# whose cell errored out never gets a done row, so it is naturally excluded
# once we inner-join against DONE_CSV).
# ---------------------------------------------------------------------------
do_summarize <- function() {
  if (!file.exists(MAIN_CSV) || !file.exists(DONE_CSV)) { cat("no results yet\n"); return(invisible(NULL)) }
  df   <- read.csv(MAIN_CSV, stringsAsFactors = FALSE)
  done <- read.csv(DONE_CSV, stringsAsFactors = FALSE)

  # De-dup key includes bin_type/bin (unlike the done file's cell key): 85's
  # main CSV is 20 rows per cell, not one, so a key of just
  # (setting_id,d,rep,method) would collapse all 20 legitimate bin rows of a
  # cell down to 1, keeping only whichever bin happened to be written last.
  key_df <- paste(df$setting_id, df$d, df$rep, df$method, df$bin_type, df$bin, sep = "\r")
  df <- df[!duplicated(key_df, fromLast = TRUE), ]
  key_done <- paste(done$setting_id, done$d, done$rep, done$method, sep = "\r")
  done <- done[!duplicated(key_done, fromLast = TRUE), ]

  df$cell_key <- paste(df$setting_id, df$d, df$rep, df$method, sep = "\r")
  df <- df[df$cell_key %in% key_done, ]

  rows <- list()
  for (grp in split(df, list(df$setting_id, df$d, df$method, df$bin_type, df$bin), drop = TRUE)) {
    per_rep_ratio <- ifelse(grp$n_regular_in_bin > 0, grp$n_flagged_in_bin / grp$n_regular_in_bin, NA_real_)
    nonempty <- !is.na(per_rep_ratio)
    rows[[length(rows) + 1L]] <- data.frame(
      setting_id = grp$setting_id[1], d = grp$d[1], method = grp$method[1],
      bin_type = grp$bin_type[1], bin = grp$bin[1],
      n_reps = length(unique(grp$rep)), n_reps_nonempty = sum(nonempty),
      sum_n_regular = sum(grp$n_regular_in_bin), sum_n_flagged = sum(grp$n_flagged_in_bin),
      pooled_ratio = if (sum(grp$n_regular_in_bin) > 0) sum(grp$n_flagged_in_bin) / sum(grp$n_regular_in_bin) else NA_real_,
      mean_per_rep_ratio = if (any(nonempty)) mean(per_rep_ratio[nonempty]) else NA_real_,
      se_per_rep_ratio = if (sum(nonempty) > 1) sd(per_rep_ratio[nonempty]) / sqrt(sum(nonempty)) else NA_real_,
      stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  S <- S[order(S$setting_id, S$method, S$bin_type, S$bin), ]
  write.csv(S, SUMM_CSV, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMM_CSV, nrow(S)))
  invisible(S)
}

# ---------------------------------------------------------------------------
# --run: full (subset) grid, checkpointed by budget
# ---------------------------------------------------------------------------
do_run <- function(sel, n_reps, budget, methods) {
  st <- resolve_settings(sel)
  cat(sprintf("85_wp8_boundary_fp: %d settings x %d reps x %d methods (budget %.0fs)\n",
              nrow(st), n_reps, length(methods), budget))
  t_start <- Sys.time()
  for (row in seq_len(nrow(st))) {
    s <- st[row, ]
    s_i <- match(s$setting_id, CANON$setting_id)
    for (rep_id in seq_len(n_reps)) {
      if (as.numeric(difftime(Sys.time(), t_start, units = "secs")) > budget) {
        cat(sprintf("[budget %gs exhausted -- rerun to continue] next up: %s rep %d\n",
                    budget, s$setting_id, rep_id))
        return(invisible("budget_stop"))
      }
      seed <- BASE_SEED + 100000L * s_i + rep_id
      for (m in methods) run_cell(s$setting_id, s$generator, s$d, rep_id, seed, m, MAIN_CSV, DONE_CSV, CLUST_CSV)
    }
  }
  cat("85_wp8_boundary_fp: run complete.\n")
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
    s <- CANON[CANON$generator == "uniform" & CANON$d == 3L, ][1, ]
    cat(sprintf("==== 85_wp8_boundary_fp --smoke: %s, 1 rep, %d methods ====\n",
                s$setting_id, length(DEFAULT_METHODS)))
    s_i <- match(s$setting_id, CANON$setting_id)
    seed <- BASE_SEED + 100000L * s_i + 9001L
    for (m in DEFAULT_METHODS) {
      run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_DONE, SMOKE_CLUST)
    }
    cat("---- re-invoking the same cells to demonstrate done-file skip ----\n")
    for (m in DEFAULT_METHODS) {
      run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_DONE, SMOKE_CLUST)
    }
    cat("ALL_CELLS_COMPLETE\n")
    quit(save = "no", status = 0)
  } else if (MODE_SUMM) {
    do_summarize()
    quit(save = "no", status = 0)
  } else {
    sel     <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
    n_reps  <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 100L
    budget  <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 480
    methods <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else DEFAULT_METHODS
    do_run(sel, n_reps, budget, methods)
  }
}
