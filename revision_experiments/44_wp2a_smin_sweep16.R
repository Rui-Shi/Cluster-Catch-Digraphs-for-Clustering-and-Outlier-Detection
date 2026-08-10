#!/usr/bin/env Rscript
# revision_experiments/44_wp2a_smin_sweep16.R
#
# WP2(a) S_min sensitivity sweep, extended to all 16 real data sets.
# Builds on 21_regen_smin_grid.R / 15_wp0_constant_smin.R / 16_wp0_su_highsmin.R
# (read, not modified). Those covered 5 datasets (hepatitis, glass, vertebral,
# ecoli, stamps) with an 8-point grid; this covers all 16 with a 10-point
# grid. Their existing output, results/tr1/regen_smin_grid.csv, is READ-ONLY
# and is never written to by this script -- a separate output file is used.
#
# Only SU-MCCD and SUN-MCCD take min.cls (a PROPORTION of n -- see
# R/ccds/RK_CCD_New.R:110-113, R/ccds/UN_CCD.R -- passing a raw count
# collapses lenD to 0). U-MCCD / UN-MCCD are out of scope (no min.cls arg).
#
# USAGE (chunked; each invocation must finish inside the calling shell's
# 600s hard timeout, so a per-cell internal elapsed-time limit below that is
# always in force):
#
#   Rscript 44_wp2a_smin_sweep16.R <datasets> <methods> <smins> <budget_sec>
#
#   <datasets>   comma list of dataset names, or "ALL" (all 16, small-first
#                order), or one of the named groups: SMALL8, MEDIUM2, LARGE6
#   <methods>    comma list, default "SU-MCCD,SUN-MCCD"
#   <smins>      comma list of numeric S_min values, or "ALL" for the full
#                10-point grid
#   <budget_sec> stop STARTING new cells once this many seconds have elapsed
#                in THIS invocation (a cell already running is allowed to
#                finish or hit its own internal timeout); default 520,
#                comfortably under the 600s hard external limit.
#
#   Rscript 44_wp2a_smin_sweep16.R --fill-not-run
#     Sweeps the full 16 x 2 x 10 = 320-cell space; for every (dataset,
#     method, s_min) combination not already present in the output CSV,
#     appends a status="not_run" row with a note citing the measured
#     per-cell cost (mean elapsed_sec already observed for that dataset in
#     the CSV, if any) that justified not attempting it.
#
# Checkpointed on (dataset, method, s_min) via a float-safe cell_done() (NOT
# harness.R's has_result(), which does a string compare and is unsafe for
# floats whose printed form drops trailing zeros, e.g. 0.10 -> "0.1").

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

REPO_ROOT <- here::here()
OUT_CSV   <- file.path(REPO_ROOT, "revision_experiments/results/tr1/wp2a_smin_sweep16.csv")

METHODS_ALL <- c("SU-MCCD", "SUN-MCCD")
stopifnot(all(METHODS_ALL %in% names(METHOD_REGISTRY)))

SMIN_GRID <- c(0.005, 0.01, 0.02, 0.03, 0.05, 0.0625, 0.10, 0.15, 0.20, 0.30)

SMALL8  <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli", "stamps", "WDBC")
MEDIUM2 <- c("pima", "shuffle")
LARGE6  <- c("vowels", "PenDigits", "waveform", "thyroid", "pageblocks", "wilt")
DATASETS_ALL <- c(SMALL8, MEDIUM2, LARGE6)   # small-first order, per the task brief

CELL_TIMEOUT_SEC <- 550   # internal cap; leaves buffer under the 600s hard external limit

ROW_COLS <- c("dataset", "method", "s_min", "n", "d",
              "TPR", "TNR", "BA", "F2",
              "n_flagged", "n_clusters", "median_cluster_size", "cluster_sizes",
              "flagged_idx", "elapsed_sec", "status", "note", "timestamp")

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
collapse_num <- function(v) if (length(v) == 0) "" else paste(v, collapse = ";")

#' Float-safe existence check: is (dataset, method, s_min) already recorded?
cell_done <- function(csv_path, dataset, method, s_min) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  any(df$dataset == dataset & df$method == method & abs(df$s_min - s_min) < 1e-9)
}

run_cell <- function(dataset, method, s_min) {
  if (cell_done(OUT_CSV, dataset, method, s_min)) {
    cat(sprintf("[skip, already recorded] %-12s x %-9s s_min=%.4f\n", dataset, method, s_min))
    return(invisible("skip"))
  }

  dat <- load_real_dataset(dataset)
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n

  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT_SEC, transient = TRUE)
    res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = s_min)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

    if (length(res$score) != n || anyNA(res$score)) {
      stop(sprintf("score sanity check failed: length=%d (n=%d), any NA=%s",
                    length(res$score), n, anyNA(res$score)))
    }
    thr <- REAL_DATA_THRESHOLDS[[method]]
    m <- evaluate(Y, res$score, thr)

    cl <- res$cluster
    cl <- cl[!is.na(cl)]
    tab <- table(cl)
    csizes <- as.integer(tab)

    list(m = m, score = res$score, n_clusters = length(tab),
         median_cluster_size = if (length(csizes)) median(csizes) else NA_real_,
         cluster_sizes = csizes,
         status = "ok", note = NA_character_)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    is_timeout <- grepl("elapsed time limit|reached CPU time limit", msg)
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
         score = NULL, n_clusters = NA_integer_, median_cluster_size = NA_real_,
         cluster_sizes = integer(0),
         status = if (is_timeout) sprintf("timeout(>%ds)", CELL_TIMEOUT_SEC) else "error",
         note = substr(msg, 1, 300))
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  flagged <- if (is.null(out$score)) integer(0) else which(out$score == 1)

  row <- as.data.frame(list(
    dataset = dataset, method = method, s_min = s_min, n = n, d = d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
    n_flagged = if (is.null(out$score)) NA_integer_ else length(flagged),
    n_clusters = out$n_clusters, median_cluster_size = out$median_cluster_size,
    cluster_sizes = collapse_num(out$cluster_sizes),
    flagged_idx = collapse_idx(flagged),
    elapsed_sec = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ), stringsAsFactors = FALSE)[ROW_COLS]

  append_result(OUT_CSV, row)

  cat(sprintf("  %-12s x %-9s s_min=%-7.4f n=%-5d d=%-3d BA=%s F2=%s nflag=%-5s ncls=%-3s %-16s wall=%.2fs\n",
              dataset, method, s_min, n, d,
              if (is.na(out$m[["BA"]])) "  NA  " else sprintf("%.3f", out$m[["BA"]]),
              if (is.na(out$m[["F2"]])) "  NA  " else sprintf("%.3f", out$m[["F2"]]),
              if (is.null(out$score)) "NA" else length(flagged),
              ifelse(is.na(out$n_clusters), "NA", out$n_clusters),
              out$status, wall))
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# --fill-not-run mode
# ---------------------------------------------------------------------------
fill_not_run <- function() {
  existing <- if (file.exists(OUT_CSV)) read.csv(OUT_CSV, stringsAsFactors = FALSE) else NULL
  mean_cost <- function(dataset) {
    if (is.null(existing)) return(NA_real_)
    sub <- existing[existing$dataset == dataset & existing$status == "ok", ]
    if (nrow(sub) == 0) return(NA_real_)
    mean(sub$elapsed_sec)
  }
  n_added <- 0L
  for (dataset in DATASETS_ALL) {
    mc <- mean_cost(dataset)
    for (method in METHODS_ALL) {
      for (s_min in SMIN_GRID) {
        if (cell_done(OUT_CSV, dataset, method, s_min)) next
        note <- if (is.na(mc)) {
          "not_run: dataset never attempted in this sweep (no measured cost)"
        } else {
          sprintf("not_run: skipped to stay within the ~40min large-set compute budget; measured mean cost for this dataset (ok cells) = %.1fs/cell", mc)
        }
        dat_nd <- tryCatch(load_real_dataset(dataset), error = function(e) NULL)
        row <- as.data.frame(list(
          dataset = dataset, method = method, s_min = s_min,
          n = if (!is.null(dat_nd)) dat_nd$n else NA_integer_,
          d = if (!is.null(dat_nd)) dat_nd$d else NA_integer_,
          TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
          n_flagged = NA_integer_, n_clusters = NA_integer_, median_cluster_size = NA_real_,
          cluster_sizes = "", flagged_idx = "", elapsed_sec = NA_real_,
          status = "not_run", note = note,
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        ), stringsAsFactors = FALSE)[ROW_COLS]
        append_result(OUT_CSV, row)
        n_added <- n_added + 1L
      }
    }
  }
  cat(sprintf("44_wp2a_smin_sweep16.R --fill-not-run: appended %d not_run rows.\n", n_added))
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1 && args[1] == "--fill-not-run") {
  fill_not_run()
  quit(save = "no", status = 0)
}

resolve_datasets <- function(s) {
  if (is.null(s) || !nzchar(s) || toupper(s) == "ALL") return(DATASETS_ALL)
  if (toupper(s) == "SMALL8")  return(SMALL8)
  if (toupper(s) == "MEDIUM2") return(MEDIUM2)
  if (toupper(s) == "LARGE6")  return(LARGE6)
  strsplit(s, ",")[[1]]
}
resolve_methods <- function(s) if (is.null(s) || !nzchar(s) || toupper(s) == "ALL") METHODS_ALL else strsplit(s, ",")[[1]]
resolve_smins   <- function(s) if (is.null(s) || !nzchar(s) || toupper(s) == "ALL") SMIN_GRID else as.numeric(strsplit(s, ",")[[1]])

DATASETS   <- resolve_datasets(if (length(args) >= 1) args[1] else NULL)
METHODS    <- resolve_methods(if (length(args) >= 2) args[2] else NULL)
SMINS      <- resolve_smins(if (length(args) >= 3) args[3] else NULL)
BUDGET_SEC <- if (length(args) >= 4) as.numeric(args[4]) else 520

cat(sprintf("44_wp2a_smin_sweep16.R: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("44_wp2a_smin_sweep16.R: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("44_wp2a_smin_sweep16.R: s_mins   = %s\n", paste(SMINS, collapse = ", ")))
cat(sprintf("44_wp2a_smin_sweep16.R: budget   = %.0fs this invocation; output = %s\n", BUDGET_SEC, OUT_CSV))

t_start <- Sys.time()
stopped_early <- FALSE
for (dataset in DATASETS) {
  for (method in METHODS) {
    for (s_min in SMINS) {
      elapsed_so_far <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      if (elapsed_so_far > BUDGET_SEC) {
        cat(sprintf("[budget exhausted at %.0fs > %.0fs -- stopping this invocation]\n", elapsed_so_far, BUDGET_SEC))
        stopped_early <- TRUE
        break
      }
      run_cell(dataset, method, s_min)
    }
    if (stopped_early) break
  }
  if (stopped_early) break
}

cat("\n44_wp2a_smin_sweep16.R: invocation done.\n")
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  cat(sprintf("%s now has %d rows.\n", OUT_CSV, nrow(final)))
  cat(sprintf("status counts: %s\n", paste(names(table(final$status)), table(final$status), sep = "=", collapse = ", ")))
}
