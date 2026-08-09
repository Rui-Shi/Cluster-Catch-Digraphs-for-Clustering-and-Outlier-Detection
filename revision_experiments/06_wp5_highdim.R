#!/usr/bin/env Rscript
# revision_experiments/06_wp5_highdim.R
#
# WP5: high-dimensional real-data experiments (Task T5). Evaluates baseline
# outlier detectors (LOF, DBSCAN, MST, ODIN, iForest) on all four new
# high-dimensional real datasets, and the two UN-CCD outlyingness scores
# (UNCCD-OOS, UNCCD-IOS) on the three of those datasets that are in scope for
# CCD experiments, via the harness's METHOD_REGISTRY. RKCCD is out of scope
# here: FINDINGS.md item 3 documents that the RK-CCD CSR envelope is already
# >91% zero-quantile at d=100 and 100% zero at d>=166 -- the RK machinery is
# numerically generatable at these dimensions but statistically vacuous, so
# the user decided (FINDINGS.md 3c) to drop RKCCD rows at d>=166 entirely.
#
# USAGE
#   Rscript "revision_experiments/06_wp5_highdim.R" [methods: all|baselines|ccd]
#     all        (default) baselines on every dataset + UN-CCD on its 3
#                in-scope datasets
#     baselines  LOF/DBSCAN/MST/ODIN/iForest on all 4 datasets (20 rows)
#     ccd        UNCCD-OOS/UNCCD-IOS on Arrhythmia/Musk/Speech only
#
# SCOPE (user decisions 2026-07-11, FINDINGS.md section 3c)
#   Datasets  : Arrhythmia (452x274), Musk (3062x166), InternetAds
#               (1966x1555), Speech (3686x400) -- results/datasets_csv/.
#   Baselines : ALL FOUR datasets.
#   UN-CCD OOS/IOS: Arrhythmia, Musk, Speech only. InternetAds is excluded by
#               deliberate scope decision (not merely a missing-table
#               accident) and is therefore never visited by the CCD method
#               loop at all -- it produces no row of any kind under mode=ccd.
#   No RKCCD  : dropped at these dimensions per the degeneracy finding above.
#
# DATASET LOADING -- deviation from the task's assumption: Musk/Arrhythmia/
# Speech/InternetAds are NOT loadable via the harness's load_real_dataset().
# That function sources data/outlier_detection/RealData_Collection.R and
# looks up an R object matching the dataset name; RealData_Collection.R only
# defines objects for the original 10 datasets (hepatitis, glass, WBC, ...)
# and was never touched by the high-dim data-acquisition task (02_load_data.R
# loads the 4 new datasets from separate raw CSVs and writes them straight to
# results/datasets_csv/*.csv without registering any RealData_Collection.R
# object). Calling load_real_dataset("Musk") would source RealData_
# Collection.R (a no-op find) and then fail on get("Musk", envir=.GlobalEnv).
# Instead we read directly from results/datasets_csv/<Name>.csv, mirroring
# the CSV-fallback idiom already used for cache-only datasets in
# 10_wp3_real.R (lines 99-106 there). harness.R itself is not touched.
#
# REAL-DATA PROTOCOL (matches the paper / harness registry, unchanged):
#   - OS threshold = 2 for UNCCD-OOS/IOS (REAL_DATA_THRESHOLDS in harness.R).
#   - Baselines use the registry's built-in real-data defaults (LOF MinPts
#     11:30 Thresh=1.5; DBSCAN k=4 with ORACLE contamination via Y, matching
#     Real_Data_DBSCAN.R; MST cont=0.02/thresh=1.2, the Real_Data_MST.R
#     defaults; ODIN k=round(sqrt(n)), indegree_threshold=round(n^(1/3)),
#     the Real_Data_ODIN.R defaults; iForest ntrees=1000, sample_size=256,
#     threshold=0.55) -- called with NO extra overrides, unlike
#     04_wp4_runtime.R's synthetic-data grid which overrides MST's cont to
#     match its own generator's known contamination. Metrics via the
#     harness's evaluate() (TPR/TNR/BA/F2; includes the regulars-first /
#     outliers-last reorder guard needed for any dataset whose native row
#     order violates count_scores2's positional contract).
#   - UN-CCD radius-search direction: "ascend" for d<=5, "descend" for d>=10
#     (FINDINGS.md / 04_wp4_runtime.R's unccd_dir_for_d convention). All
#     three WP5 CCD dimensions (166/274/400) are >=10 -> "descend".
#
# MISSING QUANTILE TABLES: no NN(UN-CCD) quantile table exists yet at
# d in {166, 274, 400} (FINDINGS.md 3a/3c; Phase B table generation is a
# separate, later task). UNCCD-OOS/IOS cells at these d get ONE marker row
# (status = SKIPPED_NO_TABLE) and are skipped -- mirrors 04_wp4_runtime.R's
# SKIPPED_NO_TABLE pattern exactly. Marker rows never block a rerun: once a
# table lands, rerunning (any mode covering that dataset/method) detects the
# table via has_nn_table() and fills in the real row; the old marker becomes
# stale and is dropped in favor of the real row at aggregation time.
#
# CHECKPOINTING: every completed (dataset, method) evaluation appends
# immediately to results/wp5_highdim_raw.csv, keyed by (dataset, method).
# Skip-if-present: any existing row with status != SKIPPED_NO_TABLE blocks a
# rerun of that key (delete the row to force a retry, same philosophy as
# 04_wp4_runtime.R's FLAGGED_TIMEOUT rows). Score vectors are cached
# separately to results/scores_cache/<dataset>_<method>.rds (skip recompute
# if cached) so a metrics-only rerun (e.g. after a threshold change) is
# free. Aggregation (results/wp5_highdim_metrics.csv) is rewritten in full
# at the end of every invocation.
#
# TIMEOUTS: none imposed on baselines (per task instruction -- LOF/ODIN kNN
# and the MST O(n^2) edge-pruning pass can be slow-ish on Speech/InternetAds;
# that is expected and reported, not aborted). Every method call is wrapped
# in tryCatch so one failure records an ERROR row and moves on rather than
# killing the run. A console NOTE is printed (not enforced) if any single
# (dataset, method) exceeds 45 minutes of wall clock.
#
# SCORE SANITY: mirrors 03_smoke.R's assert_scores_sane/assert_metrics_sane
# convention -- NA/NaN scores are always an error; -Inf is always an error;
# +Inf is allowed ONLY for UNCCD-OOS (the manuscript's documented convention
# for a point with an empty outbound-neighbor set, OOS(x_i) := +infinity).
# Metrics must land in [0,1]. Any violation records an ERROR row (message
# truncated to 200 chars) instead of crashing the run.
#
# Run (from the CLONE root, via PowerShell -- Rscript under Bash segfaults):
#   Rscript "revision_experiments/06_wp5_highdim.R" baselines
#   Rscript "revision_experiments/06_wp5_highdim.R" ccd
#   Rscript "revision_experiments/06_wp5_highdim.R" all

suppressPackageStartupMessages({
  library(here)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[[1]] else "all"
if (!(MODE %in% c("all", "baselines", "ccd"))) {
  stop("Usage: Rscript \"revision_experiments/06_wp5_highdim.R\" [methods: all|baselines|ccd]",
       call. = FALSE)
}

# Pin math-library threading (comparability with the rest of the pipeline;
# also plays nicely with any other concurrent single-core R job).
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

source(here::here("revision_experiments/harness.R"))

SHARED_DIR       <- here::here("revision_experiments/results")       # shared infra: not project-specific
RESULTS_DIR      <- file.path(SHARED_DIR, "tr2")
DATASETS_CSV_DIR <- file.path(SHARED_DIR, "datasets_csv")
CACHE_DIR        <- file.path(SHARED_DIR, "scores_cache")
RAW_CSV          <- file.path(RESULTS_DIR, "wp5_highdim_raw.csv")
METRICS_CSV      <- file.path(RESULTS_DIR, "wp5_highdim_metrics.csv")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== 06_wp5_highdim.R === mode=%s\n", MODE))

# ---------------------------------------------------------------------------
# Scope (FIXED enumeration + FIXED execution order: cheap first)
# ---------------------------------------------------------------------------
ALL_DATASETS     <- c("Arrhythmia", "Musk", "InternetAds", "Speech")
CCD_DATASETS     <- c("Arrhythmia")   # InternetAds excluded by decision; Musk + Speech
                                      # excluded by decision 2026-07-16: high-d (d>=50)
                                      # RK/NN quantile tables are simulated for dataset
                                      # size n<=1000 only (harness clamps beyond -- see
                                      # FINDINGS §7); Musk n=3062, Speech n=3686
BASELINE_METHODS <- c("LOF", "DBSCAN", "MST", "ODIN", "iForest")
CCD_METHODS      <- c("UNCCD-OOS", "UNCCD-IOS")

methods_for_mode <- function(mode, dataset) {
  m <- character(0)
  if (mode %in% c("all", "baselines")) m <- c(m, BASELINE_METHODS)
  if (mode %in% c("all", "ccd") && dataset %in% CCD_DATASETS) m <- c(m, CCD_METHODS)
  m
}

# +Inf allowance: only UNCCD-OOS carries the documented isolation convention
# in this driver's method set (RKCCD is out of scope for WP5).
is_oos_method <- function(m) identical(m, "UNCCD-OOS")

# ---------------------------------------------------------------------------
# Dataset loader (direct CSV; see header for why load_real_dataset() does
# not apply to these 4 datasets)
# ---------------------------------------------------------------------------
load_wp5_dataset <- function(name) {
  path <- file.path(DATASETS_CSV_DIR, paste0(name, ".csv"))
  stopifnot("dataset CSV not found" = file.exists(path))
  df <- read.csv(path)
  stopifnot("expected final column named 'label'" = colnames(df)[ncol(df)] == "label")
  d <- ncol(df) - 1
  Y <- df[[ncol(df)]]
  list(X = as.matrix(df[, seq_len(d), drop = FALSE]), Y = Y, d = d, n = nrow(df),
       n_outliers = sum(Y == 0))
}

# ---------------------------------------------------------------------------
# UN-CCD dispatch direction + quantile-table availability
# ---------------------------------------------------------------------------
unccd_dir_for_d <- function(d) if (d <= 5) "ascend" else "descend"

has_nn_table <- function(d) {
  file.exists(file.path(NN_QUANT_TABLE_DIR,
                        sprintf("NN-test-simul_%dd_%s%%.RData", d, nn_quant_for_d(d))))
}

# ---------------------------------------------------------------------------
# Score cache
# ---------------------------------------------------------------------------
cache_path <- function(dataset, method) file.path(CACHE_DIR, sprintf("%s_%s.rds", dataset, method))

# ---------------------------------------------------------------------------
# Raw CSV checkpoint history (status-aware; generalizes 04_wp4_runtime.R's
# HIST/rep_recorded/skip_marker_exists pattern from (cell, method, rep) keys
# down to (dataset, method) keys -- WP5 has exactly one evaluation per pair)
# ---------------------------------------------------------------------------
ROW_COLS <- c("dataset", "method", "n", "d", "n_outliers", "t_construct", "t_total",
              "wall_seconds", "TPR", "TNR", "BA", "F2", "status", "timestamp")

HIST <- new.env()
HIST$df <- if (file.exists(RAW_CSV)) read.csv(RAW_CSV, stringsAsFactors = FALSE) else NULL

hist_subset <- function(dataset, method) {
  df <- HIST$df
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df[df$dataset == dataset & df$method == method, , drop = FALSE]
}

row_recorded <- function(dataset, method) {
  sub <- hist_subset(dataset, method)
  if (is.null(sub) || nrow(sub) == 0) return(FALSE)
  any(sub$status != "SKIPPED_NO_TABLE")   # a stale skip marker never blocks a rerun
}

skip_marker_exists <- function(dataset, method) {
  sub <- hist_subset(dataset, method)
  if (is.null(sub) || nrow(sub) == 0) return(FALSE)
  any(sub$status == "SKIPPED_NO_TABLE")
}

append_and_track <- function(row) {
  stopifnot(identical(names(row), ROW_COLS))
  append_result(RAW_CSV, row)
  new_df <- as.data.frame(row, stringsAsFactors = FALSE)
  HIST$df <- if (is.null(HIST$df)) new_df else rbind(HIST$df, new_df)
}

NA_METRICS <- c(TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_)

make_row <- function(dataset, method, ds, t_construct, t_total, wall_seconds, metrics, status) {
  num_or_na <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else round(as.numeric(x), 4)
  list(
    dataset = dataset, method = method, n = ds$n, d = ds$d, n_outliers = ds$n_outliers,
    t_construct = num_or_na(t_construct), t_total = num_or_na(t_total),
    wall_seconds = num_or_na(wall_seconds),
    TPR = num_or_na(metrics[["TPR"]]), TNR = num_or_na(metrics[["TNR"]]),
    BA = num_or_na(metrics[["BA"]]), F2 = num_or_na(metrics[["F2"]]),
    status = status, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# ---------------------------------------------------------------------------
# Score sanity (mirrors 03_smoke.R's assert_scores_sane, made non-fatal --
# returns a problem string instead of stop()-ing, so the caller can record an
# ERROR row and continue)
# ---------------------------------------------------------------------------
score_problem <- function(score, n, method) {
  if (!is.numeric(score)) return(sprintf("score is not numeric (class %s)", class(score)[1]))
  if (length(score) != n) return(sprintf("score length %d != n %d", length(score), n))
  if (any(is.na(score))) return(sprintf("score contains %d NA/NaN value(s)", sum(is.na(score))))
  bad_inf <- is.infinite(score) & !(is_oos_method(method) & score > 0)
  if (any(bad_inf)) {
    return(sprintf(
      "score contains %d unexpected non-finite value(s) (-Inf always disallowed; +Inf only allowed for UNCCD-OOS)",
      sum(bad_inf)))
  }
  NULL
}

# ---------------------------------------------------------------------------
# Per (dataset, method) processing
# ---------------------------------------------------------------------------
process_one <- function(dataset, method, ds) {
  if (row_recorded(dataset, method)) {
    cat(sprintf("  [%-10s] already recorded -- skip\n", method))
    return(invisible(NULL))
  }

  if (method %in% CCD_METHODS && !has_nn_table(ds$d)) {
    if (!skip_marker_exists(dataset, method)) {
      append_and_track(make_row(dataset, method, ds, NA, NA, NA, NA_METRICS, "SKIPPED_NO_TABLE"))
      cat(sprintf("  [%-10s] SKIPPED_NO_TABLE (no NN quantile table at d=%d; rerun once it is generated)\n",
                  method, ds$d))
    } else {
      cat(sprintf("  [%-10s] SKIPPED_NO_TABLE (marker already recorded)\n", method))
    }
    return(invisible(NULL))
  }

  extra <- list()
  if (method %in% CCD_METHODS) extra$method <- unccd_dir_for_d(ds$d)

  cp <- cache_path(dataset, method)
  from_cache <- file.exists(cp)
  wall <- NA_real_

  if (from_cache) {
    rec <- readRDS(cp)
    cat(sprintf("  [%-10s] score loaded from cache (%s)\n", method, basename(cp)))
  } else {
    t0 <- Sys.time()
    res <- tryCatch(
      do.call(METHOD_REGISTRY[[method]], c(list(X = ds$X, d = ds$d, Y = ds$Y), extra)),
      error = function(e) e
    )
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(res, "error")) {
      msg <- conditionMessage(res)
      status <- paste0("ERROR: ", substr(gsub("[\r\n,]+", " ", msg), 1, 200))
      append_and_track(make_row(dataset, method, ds, NA, NA, wall, NA_METRICS, status))
      cat(sprintf("  [%-10s] %s after %.1fs\n", method, status, wall))
      return(invisible(NULL))
    }
    rec <- list(score = res$score, t_construct = res$t_construct, t_total = res$t_total,
                t_wall = wall, dataset = dataset, method = method, n = ds$n, d = ds$d)
  }

  prob <- score_problem(rec$score, ds$n, method)
  if (!is.null(prob)) {
    status <- paste0("ERROR: ", prob)
    append_and_track(make_row(dataset, method, ds, rec$t_construct, rec$t_total, wall, NA_METRICS, status))
    cat(sprintf("  [%-10s] %s\n", method, status))
    return(invisible(NULL))
  }

  if (!from_cache) saveRDS(rec, cp)   # cache only sane, successful scores

  thr <- REAL_DATA_THRESHOLDS[[method]]
  v <- evaluate(ds$Y, rec$score, thr)
  if (any(is.na(v)) || any(v < -1e-9 | v > 1 + 1e-9)) {
    status <- sprintf("ERROR: metrics out of range/NA: %s", paste(round(v, 4), collapse = ","))
    append_and_track(make_row(dataset, method, ds, rec$t_construct, rec$t_total, wall, NA_METRICS, status))
    cat(sprintf("  [%-10s] %s\n", method, status))
    return(invisible(NULL))
  }

  append_and_track(make_row(dataset, method, ds, rec$t_construct, rec$t_total, wall, v, "OK"))
  cat(sprintf("  [%-10s] OK TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f t_total=%.2fs wall=%s\n",
              method, v["TPR"], v["TNR"], v["BA"], v["F2"], rec$t_total,
              if (from_cache) "cached" else sprintf("%.2fs", wall)))
  if (!from_cache && is.finite(wall) && wall > 45 * 60) {
    cat(sprintf("  [NOTE] %s / %s exceeded 45 min wall clock (%.1f min)\n", dataset, method, wall / 60))
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
aggregate_wp5 <- function() {
  if (!file.exists(RAW_CSV)) { cat("[aggregate] no raw file yet\n"); return(invisible(NULL)) }
  df <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))

  # Dedup per (dataset, method): prefer the most recent non-SKIPPED row if one
  # exists (handles a table arriving after an earlier SKIPPED_NO_TABLE marker
  # was recorded), else the most recent skip marker.
  keys <- interaction(df$dataset, df$method, drop = TRUE)
  agg <- do.call(rbind, lapply(split(seq_len(nrow(df)), keys), function(idx) {
    sub <- df[idx, , drop = FALSE]
    real <- sub[sub$status != "SKIPPED_NO_TABLE", , drop = FALSE]
    if (nrow(real) > 0) real[nrow(real), , drop = FALSE] else sub[nrow(sub), , drop = FALSE]
  }))
  rownames(agg) <- NULL

  ds_order <- ALL_DATASETS
  m_order  <- c(BASELINE_METHODS, CCD_METHODS)
  agg <- agg[order(match(agg$dataset, ds_order), match(agg$method, m_order)), ]

  write.csv(agg, METRICS_CSV, row.names = FALSE)
  cat(sprintf("[aggregate] wrote %s (%d rows)\n", METRICS_CSV, nrow(agg)))
  invisible(agg)
}

# ---------------------------------------------------------------------------
# Main loop (dataset-major, FIXED order: Arrhythmia -> Musk -> InternetAds ->
# Speech; a dataset with no in-scope method for MODE is not even loaded)
# ---------------------------------------------------------------------------
run_start <- Sys.time()

for (dataset in ALL_DATASETS) {
  methods <- methods_for_mode(MODE, dataset)
  if (length(methods) == 0) {
    cat(sprintf("\n==== %s: no methods in scope for mode=%s -- not loaded ====\n", dataset, MODE))
    next
  }
  cat(sprintf("\n==== %s (methods: %s) ====\n", dataset, paste(methods, collapse = ", ")))
  ds <- load_wp5_dataset(dataset)
  cat(sprintf("  loaded: n=%d, d=%d, n_outliers=%d (%.1f%%)\n",
              ds$n, ds$d, ds$n_outliers, 100 * ds$n_outliers / ds$n))
  for (method in methods) {
    process_one(dataset, method, ds)
  }
}

cat(sprintf("\n---- run done in %.1f min; aggregating ----\n",
            as.numeric(difftime(Sys.time(), run_start, units = "mins"))))
aggregate_wp5()
cat("=== 06_wp5_highdim.R DONE ===\n")
