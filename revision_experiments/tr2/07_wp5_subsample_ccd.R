#!/usr/bin/env Rscript
# revision_experiments/tr2/07_wp5_subsample_ccd.R
#
# WP5 follow-up (task #3 of the revision-experiment plan): UNCCD-OOS / UNCCD-IOS
# on n=1000 subsamples of Musk and Speech.
#
# WHY THIS SCRIPT EXISTS (see FINDINGS.md item 7, "Scope decision 2026-07-16"):
# 06_wp5_highdim.R deliberately excludes Musk (3062x166) and Speech (3686x400)
# from CCD scope because the RK/NN quantile tables are simulated for dataset
# size n<=1000 only; beyond that the harness's defensive clamp (harness.R
# nnccd.radi(), ~line 153) silently reuses the table's last entry, i.e. wrong
# critical values rather than a real error. FINDINGS.md item 7 names
# subsampling to n<=1000 as the practical mitigation. This script is that
# mitigation: it runs UNCCD-OOS/IOS on seeded n=1000 subsamples
# (Musk_sub1000.csv / Speech_sub1000.csv, results/datasets_csv/, contamination
# rate preserved, seed 20260716, single stratified draw) using the NN quantile
# tables at d=166 / d=400 that were generated for exactly this purpose
# (FINDINGS.md item 7: both tables complete as of 2026-07-16).
#
# USAGE (one single-threaded process per (dataset, method) cell; this is the
# unit the orchestrator launches as a detached job):
#   Rscript "revision_experiments/tr2/07_wp5_subsample_ccd.R" <dataset> <method>
#     <dataset>  Musk_sub1000 | Speech_sub1000
#     <method>   UNCCD-OOS | UNCCD-IOS
#
# DATASET LOADING: identical convention to 06_wp5_highdim.R's
# load_wp5_dataset() -- read directly from results/datasets_csv/<dataset>.csv
# (feature columns V1..Vd, final column named "label", 1=regular/0=outlier).
#
# UN-CCD DIRECTION: unccd_dir_for_d(d) = "descend" for d>=10 (matches
# 04_wp4_runtime.R's convention, reused verbatim in 06_wp5_highdim.R). Both
# Musk (d=166) and Speech (d=400) get "descend".
#
# QUANTILE TABLE: resolved by harness.R's get_simul("NN", d) exactly as the
# registry's unccd_oos_method/unccd_ios_method already do internally
# (nn_quant_for_d(166) = nn_quant_for_d(400) = "999", i.e.
# R/NN-test_quantile/NN-test-simul_{166,400}d_999%.RData). No override passed
# here; this script does not touch that resolution logic.
#
# REAL-DATA PROTOCOL: same as 06_wp5_highdim.R -- OS threshold = 2
# (REAL_DATA_THRESHOLDS in harness.R), metrics via harness's evaluate()
# (TPR/TNR/BA/F2, regulars-first reorder guard included).
#
# SCORE SANITY: identical convention to 06_wp5_highdim.R's score_problem() --
# NA/NaN always an error, -Inf always an error, +Inf allowed ONLY for
# UNCCD-OOS (manuscript's documented isolation convention, OOS(x_i):=+Inf
# when N_O(x_i)=empty).
#
# CHECKPOINTING: appends to results/wp5_subsample_raw.csv, keyed by
# (dataset, method); an existing OK row for that key blocks a rerun (delete
# the row to force a retry). Score vectors cached to
# results/scores_cache/<dataset>_<method>.rds (skip recompute if cached, same
# as 06). This is a NEW raw CSV -- 06's wp5_highdim_raw.csv is never touched.
#
# Run (from the CLONE root, via PowerShell -- Rscript under Bash segfaults):
#   Rscript "revision_experiments/tr2/07_wp5_subsample_ccd.R" "Musk_sub1000" "UNCCD-OOS"
#   Rscript "revision_experiments/tr2/07_wp5_subsample_ccd.R" "Musk_sub1000" "UNCCD-IOS"
#   Rscript "revision_experiments/tr2/07_wp5_subsample_ccd.R" "Speech_sub1000" "UNCCD-OOS"
#   Rscript "revision_experiments/tr2/07_wp5_subsample_ccd.R" "Speech_sub1000" "UNCCD-IOS"

suppressPackageStartupMessages({
  library(here)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript \"revision_experiments/tr2/07_wp5_subsample_ccd.R\" <dataset: Musk_sub1000|Speech_sub1000> <method: UNCCD-OOS|UNCCD-IOS>",
       call. = FALSE)
}
DATASET <- args[[1]]
METHOD  <- args[[2]]

# Dataset name is validated by load_dataset() below (file must exist under
# results/datasets_csv/<dataset>.csv with a trailing "label" column) rather
# than a hard whitelist here -- production cells are Musk_sub1000 /
# Speech_sub1000, but this also lets a *_smoke CSV be smoke-tested through the
# exact same code path before launching the production jobs.
if (!(METHOD %in% c("UNCCD-OOS", "UNCCD-IOS"))) {
  stop(sprintf("Unknown method '%s'; expected UNCCD-OOS or UNCCD-IOS", METHOD), call. = FALSE)
}

SUBSAMPLE_SEED <- 20260716L  # documents the seed used to draw this CSV (gen_subsamples_wp5.R)

# Pin math-library threading (matches 06_wp5_highdim.R / 04_wp4_runtime.R; each
# production job is one single-threaded process, four of them run concurrently).
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

cat(sprintf("=== 07_wp5_subsample_ccd.R === dataset=%s method=%s pid=%d\n",
            DATASET, METHOD, Sys.getpid()))

source(here::here("revision_experiments/shared/harness.R"))

SHARED_DIR       <- here::here("revision_experiments/results")       # shared infra: not project-specific
RESULTS_DIR      <- file.path(SHARED_DIR, "tr2")
DATASETS_CSV_DIR <- file.path(SHARED_DIR, "datasets_csv")
CACHE_DIR        <- file.path(SHARED_DIR, "scores_cache")
RAW_CSV          <- file.path(RESULTS_DIR, "wp5_subsample_raw.csv")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Dataset loader (identical convention to 06_wp5_highdim.R's load_wp5_dataset())
# ---------------------------------------------------------------------------
load_dataset <- function(name) {
  path <- file.path(DATASETS_CSV_DIR, paste0(name, ".csv"))
  stopifnot("dataset CSV not found" = file.exists(path))
  df <- read.csv(path)
  stopifnot("expected final column named 'label'" = colnames(df)[ncol(df)] == "label")
  d <- ncol(df) - 1
  Y <- df[[ncol(df)]]
  list(X = as.matrix(df[, seq_len(d), drop = FALSE]), Y = Y, d = d, n = nrow(df),
       n_outliers = sum(Y == 0))
}

# UN-CCD radius-search direction (04_wp4_runtime.R / 06_wp5_highdim.R convention)
unccd_dir_for_d <- function(d) if (d <= 5) "ascend" else "descend"

is_oos_method <- function(m) identical(m, "UNCCD-OOS")

cache_path <- function(dataset, method) file.path(CACHE_DIR, sprintf("%s_%s.rds", dataset, method))

# ---------------------------------------------------------------------------
# Checkpoint guard
# ---------------------------------------------------------------------------
ROW_COLS <- c("dataset", "method", "n", "d", "n_outliers", "seed", "t_construct", "t_total",
              "wall_seconds", "TPR", "TNR", "BA", "F2", "status", "timestamp")

row_recorded_ok <- function(dataset, method) {
  if (!file.exists(RAW_CSV)) return(FALSE)
  df <- tryCatch(read.csv(RAW_CSV, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  any(df$dataset == dataset & df$method == method & df$status == "OK")
}

make_row <- function(dataset, method, ds, t_construct, t_total, wall_seconds, metrics, status) {
  num_or_na <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) NA_real_ else round(as.numeric(x), 4)
  list(
    dataset = dataset, method = method, n = ds$n, d = ds$d, n_outliers = ds$n_outliers,
    seed = SUBSAMPLE_SEED,
    t_construct = num_or_na(t_construct), t_total = num_or_na(t_total),
    wall_seconds = num_or_na(wall_seconds),
    TPR = num_or_na(metrics[["TPR"]]), TNR = num_or_na(metrics[["TNR"]]),
    BA = num_or_na(metrics[["BA"]]), F2 = num_or_na(metrics[["F2"]]),
    status = status, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

append_row <- function(row) {
  stopifnot(identical(names(row), ROW_COLS))
  append_result(RAW_CSV, row)
}

NA_METRICS <- c(TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_)

# ---------------------------------------------------------------------------
# Score sanity (mirrors 06_wp5_highdim.R's score_problem())
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
# Main (single (dataset, method) cell; tryCatch wraps everything so the
# process always prints a DONE/ERROR line rather than crashing silently)
# ---------------------------------------------------------------------------
main <- function() {
  if (row_recorded_ok(DATASET, METHOD)) {
    cat("SKIP (already recorded)\n")
    return(invisible(NULL))
  }

  ds <- load_dataset(DATASET)
  cat(sprintf("  loaded: n=%d, d=%d, n_outliers=%d (%.2f%%)\n",
              ds$n, ds$d, ds$n_outliers, 100 * ds$n_outliers / ds$n))

  cp <- cache_path(DATASET, METHOD)
  from_cache <- file.exists(cp)
  wall <- NA_real_
  extra <- list(method = unccd_dir_for_d(ds$d))

  if (from_cache) {
    rec <- readRDS(cp)
    cat(sprintf("  score loaded from cache (%s)\n", basename(cp)))
  } else {
    t0 <- Sys.time()
    res <- do.call(METHOD_REGISTRY[[METHOD]], c(list(X = ds$X, d = ds$d, Y = ds$Y), extra))
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rec <- list(score = res$score, t_construct = res$t_construct, t_total = res$t_total,
                t_wall = wall, dataset = DATASET, method = METHOD, n = ds$n, d = ds$d)
  }

  prob <- score_problem(rec$score, ds$n, METHOD)
  if (!is.null(prob)) stop(prob, call. = FALSE)

  if (!from_cache) saveRDS(rec, cp)   # cache only sane, successful scores

  thr <- REAL_DATA_THRESHOLDS[[METHOD]]
  v <- evaluate(ds$Y, rec$score, thr)
  if (any(is.na(v)) || any(v < -1e-9 | v > 1 + 1e-9)) {
    stop(sprintf("metrics out of range/NA: %s", paste(round(v, 4), collapse = ",")), call. = FALSE)
  }

  append_row(make_row(DATASET, METHOD, ds, rec$t_construct, rec$t_total, wall, v, "OK"))
  cat(sprintf("  OK TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f t_total=%.2fs wall=%s\n",
              v["TPR"], v["TNR"], v["BA"], v["F2"], rec$t_total,
              if (from_cache) "cached" else sprintf("%.2fs", wall)))
  cat(sprintf("DONE %s %s\n", DATASET, METHOD))
}

result <- tryCatch({ main(); NULL }, error = function(e) e)
if (!is.null(result)) {
  msg <- gsub("[\r\n]+", " ", conditionMessage(result))
  cat(sprintf("ERROR %s %s: %s\n", DATASET, METHOD, msg))
  # Record an ERROR row too, unless we never even got to `ds` (load failure).
  if (exists("ds", inherits = FALSE)) {
    tryCatch(
      append_row(make_row(DATASET, METHOD, ds, NA, NA, NA, NA_METRICS,
                           paste0("ERROR: ", substr(msg, 1, 200)))),
      error = function(e2) cat(sprintf("  (failed to append ERROR row: %s)\n", conditionMessage(e2)))
    )
  }
}
