#!/usr/bin/env Rscript
# revision_experiments/23_regen_baselines.R
#
# Agent C (parallel regen, 5-way split) -- BASELINES shard for
# NEUCOM-D-26-15191. Regenerates label-free DBSCAN / MST / ODIN / iForest
# rows for the eight real datasets, per REGEN_SPEC.md. LOF is NOT re-run
# (LOF/LOF.R with L_MinPts=11, U_MinPts=30, Thresh=1.5 was already label-free
# and matches the manuscript).
#
# Two confirmed defects in the published drivers, both reproduced here as
# REFERENCE rows (variant = "published-exact") rather than silently fixed:
#   1. Real_Data_DBSCAN.R:24 calls DBSCAN(df, MinPts, cont) -- the FULL data
#      frame, including the ground-truth Y column, not X. Every other
#      baseline driver passes X.
#   2. DBSCAN's eps was read off the k-distance curve at the TRUE
#      (oracle) contamination rate n0/n; MST's thresh was hand-tuned per
#      dataset (1.2 for hepatitis/glass/vertebral/ecoli/stamps, 1.3 vowels,
#      1.05 waveform, 1.4 wilt; cont=0.02 throughout).
#
# This script does NOT modify harness.R, MST_Outlier.R, DBSCAN.R, ODIN.R, or
# ISO.R -- it only calls the already-sourced functions (via harness.R) or
# dbscan::dbscan() directly, per REGEN_SPEC.md's instruction that label-free
# DBSCAN variants call dbscan::dbscan(data, eps, minPts) directly rather than
# going through the source DBSCAN() wrapper.
#
# Usage:
#   Rscript revision_experiments/23_regen_baselines.R <dataset>
#   <dataset> one of: hepatitis, glass, vertebral, ecoli, stamps, vowels,
#   waveform, wilt. One dataset per invocation, run synchronously (never
#   run_in_background) -- wilt in particular needs its own invocation, since
#   it builds a 4819x4819 distance matrix and a complete-graph MST.
#
# Output: revision_experiments/results/regen_baselines.csv (this script's
# OWN shard -- never written to by any other agent). Resumable: checkpoints
# via has_result()/append_result() keyed on (dataset, method, variant).

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))

OUT_CSV   <- here::here("revision_experiments", "results", "regen_baselines.csv")
TRUTH_CSV <- here::here("revision_experiments", "published_realdata_truth.csv")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1])) {
  stop("Usage: Rscript 23_regen_baselines.R <dataset>  (one of hepatitis, glass, vertebral, ecoli, stamps, vowels, waveform, wilt)")
}
DATASET <- args[1]

# REGEN_SPEC.md section 2 expected (n, d, n0) -- confirmed at runtime below.
EXPECTED <- list(
  hepatitis = list(n = 74,   d = 19, n0 = 6),
  glass     = list(n = 213,  d = 9,  n0 = 9),
  vertebral = list(n = 240,  d = 6,  n0 = 30),
  ecoli     = list(n = 336,  d = 7,  n0 = 8),
  stamps    = list(n = 340,  d = 9,  n0 = 31),
  vowels    = list(n = 1452, d = 12, n0 = 46),
  waveform  = list(n = 3443, d = 21, n0 = 100),
  wilt      = list(n = 4819, d = 5,  n0 = 257)
)
if (!(DATASET %in% names(EXPECTED))) {
  stop(sprintf("Unknown dataset '%s'. Must be one of: %s", DATASET, paste(names(EXPECTED), collapse = ", ")))
}

truth <- read.csv(TRUTH_CSV, stringsAsFactors = FALSE)
get_published <- function(dataset, method) {
  sub <- truth[tolower(truth$dataset) == tolower(dataset) & truth$method == method, ]
  if (nrow(sub) == 0) return(NULL)
  v <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
  for (i in seq_len(nrow(sub))) if (sub$metric[i] %in% names(v)) v[sub$metric[i]] <- sub$value[i]
  v
}

# Normalized-Kneedle knee finder -- deterministic, parameter-free, exactly
# per REGEN_SPEC.md.
knee_eps <- function(kd) {
  y <- sort(kd); m <- length(y)
  if (m < 3 || y[m] <= y[1]) return(stats::median(kd))
  xn <- (seq_len(m) - 1) / (m - 1)
  yn <- (y - y[1]) / (y[m] - y[1])
  y[which.max(abs(yn - xn))]
}

# k-distance vector exactly as DBSCAN.R computes it:
#   dist_M = as.matrix(dist(data)); k_dist[i] = sort(dist_M[i,])[k+1]
kdist_vec <- function(data_mat, k) {
  distM <- as.matrix(dist(data_mat))
  apply(distM, 1, function(row) sort(row)[k + 1])
}

ROW_COLS <- c("dataset", "method", "variant", "params", "n", "d",
              "TPR", "TNR", "BA", "F2",
              "published_TPR", "published_TNR", "published_BA", "published_F2",
              "max_abs_diff", "match_3dp", "t_total", "status", "note", "timestamp")

build_row <- function(dataset, method, variant, params, n, d, metrics, t_total, status, note) {
  pub <- get_published(dataset, method)
  has_metrics <- !is.null(metrics) && !anyNA(metrics)
  if (is.null(pub) || !has_metrics) {
    published_vals <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
    max_diff <- NA_real_; match3 <- NA
  } else {
    published_vals <- pub
    diffs <- abs(metrics[c("TPR", "TNR", "BA", "F2")] - pub[c("TPR", "TNR", "BA", "F2")])
    max_diff <- max(diffs)
    match3 <- all(round(metrics[c("TPR", "TNR", "BA", "F2")], 3) == round(pub[c("TPR", "TNR", "BA", "F2")], 3))
  }
  if (is.null(metrics)) metrics <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
  data.frame(
    dataset = dataset, method = method, variant = variant, params = params,
    n = n, d = d,
    TPR = unname(metrics["TPR"]), TNR = unname(metrics["TNR"]),
    BA = unname(metrics["BA"]), F2 = unname(metrics["F2"]),
    published_TPR = unname(published_vals["TPR"]), published_TNR = unname(published_vals["TNR"]),
    published_BA = unname(published_vals["BA"]), published_F2 = unname(published_vals["F2"]),
    max_abs_diff = max_diff, match_3dp = match3,
    t_total = t_total, status = status,
    note = if (is.null(note) || is.na(note)) NA_character_ else note,
    timestamp = format(Sys.time()),
    stringsAsFactors = FALSE
  )[ROW_COLS]
}

run_cell <- function(dataset, method, variant, fn) {
  keys <- c(dataset = dataset, method = method, variant = variant)
  if (has_result(OUT_CSV, keys)) {
    cat(sprintf("[skip, already recorded] %s x %s x %s\n", dataset, method, variant))
    return(invisible(NULL))
  }
  t0 <- Sys.time()
  out <- tryCatch(list(ok = fn()), error = function(e) list(err = conditionMessage(e)))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(out$err)) {
    row <- build_row(dataset, method, variant, params = NA_character_, n = NA_integer_, d = NA_integer_,
                      metrics = NULL, t_total = wall, status = "error", note = substr(out$err, 1, 300))
    cat(sprintf("  %-10s x %-8s x %-16s ERROR: %s\n", dataset, method, variant, out$err))
  } else {
    r <- out$ok
    row <- build_row(dataset, method, variant, params = r$params, n = r$n, d = r$d,
                      metrics = r$metrics, t_total = wall, status = "ok", note = r$note)
    cat(sprintf("  %-10s x %-8s x %-16s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | %s | %.2fs\n",
                dataset, method, variant, r$metrics["TPR"], r$metrics["TNR"], r$metrics["BA"], r$metrics["F2"],
                r$params, wall))
  }
  append_result(OUT_CSV, row)
  invisible(row)
}

# ---------------------------------------------------------------------------
cat(sprintf("=== 23_regen_baselines.R : %s ===\n", DATASET))
dat <- load_real_dataset(DATASET)
X <- dat$X; Y <- dat$Y; n <- dat$n; d <- dat$d
n0 <- sum(Y == 0)
exp_ <- EXPECTED[[DATASET]]
if (n != exp_$n || d != exp_$d || n0 != exp_$n0) {
  cat(sprintf("*** DISAGREEMENT with REGEN_SPEC.md table: loaded n=%d d=%d n0=%d, expected n=%d d=%d n0=%d ***\n",
              n, d, n0, exp_$n, exp_$d, exp_$n0))
} else {
  cat(sprintf("confirmed n=%d d=%d n0=%d (matches REGEN_SPEC.md)\n", n, d, n0))
}
Xmat <- as.matrix(X)

DBSCAN_THR  <- REAL_DATA_THRESHOLDS[["DBSCAN"]]
MST_THR     <- REAL_DATA_THRESHOLDS[["MST"]]
ODIN_THR    <- REAL_DATA_THRESHOLDS[["ODIN"]]
IFOREST_THR <- REAL_DATA_THRESHOLDS[["iForest"]]

# ---------------------------------------------------------------------------
# DBSCAN -- 5 variants. All go through dbscan::dbscan(data, eps, minPts)
# directly, using the SAME eps-from-k-distance formula DBSCAN.R uses
# (dist_M -> k_dist -> quantile); only (data, k, eps-source) differ. Keeping
# every variant, including the two references, on one code path means the
# resolved eps is always available to report.
# ---------------------------------------------------------------------------

kd4_X <- kdist_vec(Xmat, 4)
kd5_X <- kdist_vec(Xmat, 5)

run_cell(DATASET, "DBSCAN", "knee-4", function() {
  eps <- knee_eps(kd4_X)
  labels <- dbscan::dbscan(Xmat, eps = eps, minPts = 4)$cluster
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, DBSCAN_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("minPts=4, eps=%.4f (knee of 4-dist curve, input=X)", eps),
       note = NA_character_)
})

run_cell(DATASET, "DBSCAN", "knee-5", function() {
  eps <- knee_eps(kd5_X)
  labels <- dbscan::dbscan(Xmat, eps = eps, minPts = 5)$cluster
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, DBSCAN_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("minPts=5, eps=%.4f (knee of 5-dist curve, input=X)", eps),
       note = NA_character_)
})

run_cell(DATASET, "DBSCAN", "fixed-q05", function() {
  eps <- as.numeric(quantile(kd4_X, 0.95))
  labels <- dbscan::dbscan(Xmat, eps = eps, minPts = 4)$cluster
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, DBSCAN_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("minPts=4, eps=%.4f (quantile(4-dist,0.95), input=X)", eps),
       note = NA_character_)
})

run_cell(DATASET, "DBSCAN", "oracle-X", function() {
  q <- 1 - n0 / n
  eps <- as.numeric(quantile(kd4_X, q))
  labels <- dbscan::dbscan(Xmat, eps = eps, minPts = 4)$cluster
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, DBSCAN_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("minPts=4, eps=%.4f (quantile(4-dist, 1-n0/n=%.4f), input=X)", eps, q),
       note = "REFERENCE: isolates the oracle eps from the df->X label leak")
})

run_cell(DATASET, "DBSCAN", "published-exact", function() {
  XY <- cbind(Xmat, Y)
  kd4_XY <- kdist_vec(XY, 4)
  q <- 1 - n0 / n
  eps <- as.numeric(quantile(kd4_XY, q))
  labels <- dbscan::dbscan(XY, eps = eps, minPts = 4)$cluster
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, DBSCAN_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("minPts=4, eps=%.4f (quantile(4-dist, 1-n0/n=%.4f), input=cbind(X,Y))", eps, q),
       note = "REFERENCE: reproduces the published driver's label leak (Real_Data_DBSCAN.R:24)")
})

# ---------------------------------------------------------------------------
# MST -- 3 variants, all via the sourced MST_Outlier(data, cont, thresh).
# ---------------------------------------------------------------------------

MST_PUBLISHED_THRESH <- c(hepatitis = 1.2, glass = 1.2, vertebral = 1.2, ecoli = 1.2,
                           stamps = 1.2, vowels = 1.3, waveform = 1.05, wilt = 1.4)

run_cell(DATASET, "MST", "default", function() {
  labels <- MST_Outlier(Xmat, cont = 0.04, thresh = 2)
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, MST_THR)
  list(metrics = m, n = n, d = d,
       params = "cont=0.04, thresh=2 (MST_Outlier() function's own defaults)", note = NA_character_)
})

run_cell(DATASET, "MST", "uniform-1.2", function() {
  labels <- MST_Outlier(Xmat, cont = 0.02, thresh = 1.2)
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, MST_THR)
  list(metrics = m, n = n, d = d,
       params = "cont=0.02, thresh=1.2 (uniform across all 8 datasets, no per-dataset override)", note = NA_character_)
})

run_cell(DATASET, "MST", "published-exact", function() {
  thresh <- unname(MST_PUBLISHED_THRESH[DATASET])
  labels <- MST_Outlier(Xmat, cont = 0.02, thresh = thresh)
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, MST_THR)
  list(metrics = m, n = n, d = d,
       params = sprintf("cont=0.02, thresh=%.2f (per-dataset published value)", thresh),
       note = "REFERENCE: reproduces the published per-dataset thresh")
})

# ---------------------------------------------------------------------------
# ODIN -- 1 variant, default only (already label-free).
# ---------------------------------------------------------------------------

run_cell(DATASET, "ODIN", "default", function() {
  labels <- ODIN(Xmat)
  score <- ifelse(labels == 0, 1, 0)
  m <- evaluate(Y, score, ODIN_THR)
  k_resolved <- round(sqrt(n))
  indeg_resolved <- round(n^(1 / 3))
  list(metrics = m, n = n, d = d,
       params = sprintf("k=round(sqrt(n))=%d, indegree_threshold=round(n^(1/3))=%d", k_resolved, indeg_resolved),
       note = NA_character_)
})

# ---------------------------------------------------------------------------
# iForest -- 2 variants: published single-seed, and 10-seed mean +- sd(BA).
#
# NOTE: harness.R's iforest_method() only calls R's set.seed(seed) and then
# isotree::isolation.forest(X, ntrees=1000, sample_size=...) WITHOUT passing
# a `seed=` argument through. isotree::isolation.forest() has its own
# internal `seed = 1` default (confirmed via args(isotree::isolation.forest),
# isotree 0.6.1.5) that is independent of R's global RNG -- so routing
# through iforest_method() with different R-level set.seed() calls silently
# produces byte-identical models every time (sd(BA) = 0 by construction, not
# by finding). Per REGEN_SPEC.md's instruction to call
# isotree::isolation.forest() directly with resolved hyperparameters (same
# pattern as the direct dbscan::dbscan() calls above), both variants below
# call isolation.forest() directly and pass `seed = seed` explicitly, so the
# 10 seeds in "mean10" are genuinely different models. harness.R itself is
# untouched.
# ---------------------------------------------------------------------------

run_cell(DATASET, "iForest", "seed1", function() {
  sample_size <- min(256, n)
  model <- isotree::isolation.forest(Xmat, ntrees = 1000, sample_size = sample_size, seed = 1)
  score <- as.numeric(predict(model, Xmat))
  m <- evaluate(Y, score, IFOREST_THR)
  list(metrics = m, n = n, d = d,
       params = "ntrees=1000, sample_size=min(256,n), seed=1 (published single-seed config, isotree seed= param)",
       note = NA_character_)
})

run_cell(DATASET, "iForest", "mean10", function() {
  sample_size <- min(256, n)
  mat <- sapply(1:10, function(sd_) {
    model <- isotree::isolation.forest(Xmat, ntrees = 1000, sample_size = sample_size, seed = sd_)
    score <- as.numeric(predict(model, Xmat))
    evaluate(Y, score, IFOREST_THR)
  })
  means <- rowMeans(mat)
  sd_BA <- sd(mat["BA", ])
  list(metrics = means, n = n, d = d,
       params = "ntrees=1000, sample_size=min(256,n), seeds=1..10 (isotree seed= param), metrics=mean across seeds",
       note = sprintf("sd(BA) across 10 seeds = %.4f", sd_BA))
})

cat(sprintf("\n=== 23_regen_baselines.R : %s done ===\n", DATASET))
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  final_ds <- final[final$dataset == DATASET, ]
  cat(sprintf("%s now has %d rows total, %d for %s (status=ok: %d)\n",
              OUT_CSV, nrow(final), nrow(final_ds), DATASET, sum(final_ds$status == "ok")))
}
