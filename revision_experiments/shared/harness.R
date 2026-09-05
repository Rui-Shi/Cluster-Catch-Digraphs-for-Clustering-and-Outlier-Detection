#!/usr/bin/env Rscript
# revision_experiments/shared/harness.R
#
# Sourceable experiment harness for the revision-experiments pipeline,
# shared by the TR1 (Neurocomputing) and TR2 (Pattern Recognition) drivers.
# Provides:
#   - rk_quant_label_paper(), nn_quant_label_paper_UN(),
#     nn_quant_label_paper_SUN()             the manuscript's alpha schedule;
#                                            the ONLY definitions in the repo
#   - get_simul(variant, d, quant, n)        quantile-table loader, size-checked
#   - METHOD_REGISTRY: 9 named scoring functions (4 CCD-based OS + 5 baselines)
#   - evaluate(Y, score, threshold)          TPR/TNR/BA/F2 wrapper, guarded
#   - append_result()/has_result()           CSV checkpointing for restarts
#
# This file only READS the original R/, methods/, simulations/, data/ trees.
# Nothing in those directories is modified, and -- since 2026-09-05 -- nothing
# in them is overridden here either.
#
# REMOVED 2026-09-05: the "defensive index-clamp overrides" of
# ccd.Kest.edge.quantile() (R/ccds/RK_CCD_New.R) and nnccd.radi()
# (R/ccds/UN_CCD.R) that used to sit in section 2 of this file. Two reasons,
# both verified:
#
#   1. They were DEAD CODE. tr1/wp0_mccd_methods.R sources the four detectors
#      in methods/outlier_detection/{RU-MCCDs,SU-MCCDs,UN-MCCD,SUN-MCCD}.R,
#      and line 1 of each of those re-sources R/ccds/RK_CCD_New.R or
#      R/ccds/UN_CCD.R into .GlobalEnv -- silently reverting both overrides.
#      Every run that loaded wp0 (i.e. every MCCD run) used the originals.
#
#   2. The clamp was a hazard, not a safeguard. Several regenerated tables in
#      R/NN-test_quantile are exactly DATASET-SIZED while carrying a generic
#      filename: d=19 -> 74 entries (hepatitis), d=18 -> 148 (lymphography),
#      d=12 -> 1452 (vowels), d=16 -> 3200 (PenDigits), d=21 -> 3443
#      (waveform). Under a clamp, running any LARGER data set at one of those
#      dimensions would silently reuse the table's last critical value for
#      every point past its end. Unclamped, the original yields NA and fails
#      loudly. Loud is what we want.
#
# The size question is now handled once, up front, by get_simul()'s `n`
# argument (section 2) instead of by patching the numerics downstream.

suppressPackageStartupMessages({
  library(here)
  library(MASS)
  library(igraph)
  library(cluster)
  library(dbscan)
  library(isotree)
  library(FNN)
})

# ---------------------------------------------------------------------------
# 1. Source original (read-only) source files
# ---------------------------------------------------------------------------

# CCD-based outlyingness scores. Each of these sources its own dependency
# chain (Outlyingness_Score.R, RK_CCD_New.R / UN_CCD.R, which in turn source
# ccdfunctions.R, Kest.R / NN_Dist_Est.R) via here::here(), landing all
# definitions in .GlobalEnv.
source(here::here("methods/outlyingness_scores/RKCCD_OOS_IOS.R"))
source(here::here("methods/outlyingness_scores/UNCCD_OOS_IOS.R"))

# Baseline outlier-detection methods, exactly as used in the paper's real-data
# benchmarking scripts.
source(here::here("simulations/outlier_detection/Algo_Compare_OutlierDetection/LOF/LOF.R"))
source(here::here("simulations/outlier_detection/Algo_Compare_OutlierDetection/DBSCAN/DBSCAN.R"))
source(here::here("simulations/outlier_detection/Algo_Compare_OutlierDetection/MST/MST_Outlier.R"))
source(here::here("simulations/outlier_detection/Algo_Compare_OutlierDetection/ODIN/ODIN.R"))
source(here::here("simulations/outlier_detection/Algo_Compare_OutlierDetection/Isolation Forest/ISO.R"))

# Metrics.
source(here::here("R/general_functions/count.R"))

# ---------------------------------------------------------------------------
# 2. The paper's alpha schedule, and get_simul(): quantile-table loader
# ---------------------------------------------------------------------------
#
# SINGLE SOURCE OF TRUTH for the MC-SRT significance level per dimension.
#
# What was here before, and why it is gone. This file used to carry its own
# rk_quant_for_d()/nn_quant_for_d() buckets, reverse-engineered from the
# RKCCD_OOS_IOS / UNCCD_OOS_IOS simulation drivers, while
# tr1/wp0_mccd_methods.R carried a second, paper-derived set of resolvers.
# The two DISAGREED:
#   * rk_quant_for_d() switched at d <= 10; the manuscript switches at d < 10,
#     so d = 10 got alpha = 1% here and 0.1% in the paper.
#   * nn_quant_for_d() gave every NND method alpha = 1% for d in 10..19; the
#     manuscript gives that only to UN-MCCD, while SUN-MCCD is already at
#     0.1% from d = 10.
# Both buckets are deleted. The three resolvers below are now the only
# definitions anywhere in revision_experiments/, and every caller names the
# one it means.
#
# Authority: CCD_OutlierDetection_Neurocomputing.tex, Section "Uniform Cluster
# Settings" -- RK-based methods (U-MCCD, SU-MCCD) use alpha = 1% for d < 10
# and 0.1% for d >= 10; UN-MCCD uses alpha = 15%, 10%, 5%, 1%, 0.1% at
# d = 2, 3, 5, 10, {20,50,100}; SUN-MCCD matches UN-MCCD but reaches 0.1%
# already at d = 10 -- and Table tab:alpha_real, which tabulates the same
# schedule for the sixteen real data sets (RK switches once at d=10; UN-MCCD
# at d=10 and again at d=20; SUN-MCCD once at d=10). Levels are anchored at
# the tabulated dimensions and applied as a step function: each level holds
# from its own anchor up to, but not including, the next.
#
# The value returned is the FILE-LABEL STRING as it appears in the .RData
# filename ("85"/"90"/"95"/"99"/"999"), not a probability; that is what
# get_simul()'s `quant` argument takes.

#' RK-based methods (U-MCCD, SU-MCCD, and the RKCCD-OOS/IOS wrappers below).
#' alpha = 1% for d < 10; alpha = 0.1% for d >= 10.
rk_quant_label_paper <- function(d) if (d < 10) "99" else "999"

#' NND-based UN-MCCD. alpha = 15% (d = 2), 10% (d = 3,4), 5% (d = 5..9),
#' 1% (d = 10..19), 0.1% (d >= 20).
nn_quant_label_paper_UN <- function(d) {
  if (d <= 2) "85"
  else if (d <= 4) "90"
  else if (d <= 9) "95"
  else if (d <= 19) "99"
  else "999"
}

#' NND-based SUN-MCCD. Identical to UN-MCCD below d = 10; forces alpha = 0.1%
#' ("999") already at d = 10, one step earlier than UN-MCCD.
nn_quant_label_paper_SUN <- function(d) {
  if (d < 10) nn_quant_label_paper_UN(d) else "999"
}

RK_QUANT_TABLE_DIR <- here::here("R/RK-test_quantile")
NN_QUANT_TABLE_DIR <- here::here("R/NN-test_quantile")

#' Assert that a loaded quantile table covers a data set of `n` points.
#'
#' The tables are indexed by running point count, so a table shorter than the
#' data set is not a rounding issue -- it is missing critical values. What the
#' consumers do with the overrun is silent in both directions:
#'   * nnccd.radi() (R/ccds/UN_CCD.R) slices simul$average[1:n] /
#'     simul$median[1:n]; R pads an over-long index with NA rather than
#'     erroring, and the NAs then propagate into every radius and score.
#'   * ccd.Kest.edge.quantile() (R/ccds/RK_CCD_New.R) indexes
#'     simul$quan[[...]][j, ] and does throw "subscript out of bounds", but
#'     only partway through a run.
#' Neither is a good place to find out, so we check here, before any compute.
#'
#' Table extents are NOT uniform and must not be assumed. The RK tables are
#' generated by KestP.simpois.edge.quantile(m, d, rn, quan, niter), which
#' returns matrix(..., nrow = m): the row count is the maximum supported n,
#' and it is 1000 in some shipped files (d = 2,3,4,10@99%,20@99%,50,100) and
#' 5000 in others -- `niter` is a separate argument (the Monte-Carlo
#' replication count) and is not the row count. The NN tables are worse: the
#' regenerated 0.1% files at d = 12,16,18,19,21 are exactly dataset-sized
#' (1452, 3200, 148, 74, 3443) under an otherwise generic filename.
#'
#' Exported (not internal) so that the get_simul() shadows in the numbered
#' driver scripts, which load tables from their own directories, can apply
#' the identical check.
#'
#' @param simul the loaded table object
#' @param variant "RK" or "NN"
#' @param d dimensionality (for the message only)
#' @param n required extent; NULL skips the check
#' @param path file the table came from (for the message only)
check_simul_extent <- function(simul, variant, d, n, path = NA_character_) {
  if (is.null(n) || is.na(n)) return(invisible(NULL))
  n <- as.integer(n)
  have <- if (variant == "NN") {
    c(average = length(simul$average), median = length(simul$median))
  } else {
    c(quan = nrow(simul$quan[[1]]))
  }
  short <- have[have < n]
  if (length(short) > 0) {
    stop(sprintf(paste0(
      "get_simul(): quantile table is too short for the data.\n",
      "  variant = %s, d = %d, required n = %d\n",
      "  table extent: %s\n",
      "  file: %s\n",
      "This table cannot serve a data set of this size. Do NOT clamp the index ",
      "-- reusing the last critical value for every point past the end is a ",
      "silent error. Regenerate the table at n >= %d, or use a table that ",
      "already covers it."),
      variant, d, n,
      paste(sprintf("%s = %d", names(have), have), collapse = ", "),
      path, n))
  }
  invisible(NULL)
}

#' Resolve and load the RK-CCD or NN(UN)-CCD quantile lookup table.
#'
#' @param variant "RK" or "NN"
#' @param d dimensionality
#' @param quant REQUIRED file-label string, e.g. "99", "999", "95" (as it
#'   appears in the filename, without the trailing "%"). There is no default:
#'   the level depends on which method is asking, not only on d, so the caller
#'   must name a resolver -- rk_quant_label_paper(d) for U-MCCD/SU-MCCD and
#'   the RKCCD wrappers, nn_quant_label_paper_UN(d) for UN-MCCD and the
#'   generic UNCCD wrappers, nn_quant_label_paper_SUN(d) for SUN-MCCD.
#' @param n number of rows in the data set the table will serve. Required
#'   whenever it is known; the loaded table is checked against it and the call
#'   stops if the table is shorter. NULL only for callers that are inspecting
#'   a table rather than running a data set through it.
#' @return list(simul = <loaded object>, quant = <numeric quantile, e.g.
#'   0.999>, quant_label = <string used in filename>, file = <path used>)
get_simul <- function(variant = c("RK", "NN"), d, quant = NULL, n = NULL) {
  variant <- match.arg(variant)
  if (is.null(quant)) {
    stop(paste0(
      "get_simul(): `quant` is required and has no default.\n",
      "  The significance level is method-specific, not a function of d alone.\n",
      "  Call one of: rk_quant_label_paper(d)      (U-MCCD, SU-MCCD, RKCCD-OOS/IOS)\n",
      "               nn_quant_label_paper_UN(d)   (UN-MCCD, UNCCD-OOS/IOS)\n",
      "               nn_quant_label_paper_SUN(d)  (SUN-MCCD)\n",
      "  (The former rk_quant_for_d()/nn_quant_for_d() defaults were removed on\n",
      "   2026-09-05: they disagreed with the manuscript at d = 10.)"))
  }
  q <- quant
  if (variant == "RK") {
    fname <- sprintf("RK-test-simul_%dd_%s%%.RData", d, q)
    path <- file.path(RK_QUANT_TABLE_DIR, fname)
  } else {
    fname <- sprintf("NN-test-simul_%dd_%s%%.RData", d, q)
    path <- file.path(NN_QUANT_TABLE_DIR, fname)
  }
  if (!file.exists(path)) {
    stop(sprintf(
      "get_simul(): missing quantile table for variant=%s, d=%d, quant=%s.\nExpected file: %s",
      variant, d, q, path
    ))
  }
  e <- new.env()
  load(path, envir = e)
  if (!exists("simul", envir = e)) {
    stop(sprintf("get_simul(): file %s does not contain an object named 'simul'", path))
  }
  simul <- get("simul", envir = e)
  check_simul_extent(simul, variant, d, n, path)
  list(
    simul = simul,
    quant = as.numeric(paste0("0.", q)),
    quant_label = q,
    file = path
  )
}

# ---------------------------------------------------------------------------
# 3. Method registry
# ---------------------------------------------------------------------------
#
# Every entry is function(X, d, Y = NULL, ...) returning
#   list(score = numeric length n, larger = more outlying,
#        t_construct = seconds or NA,
#        t_total = seconds)
#
# Score-polarity notes (documented per task instructions):
#   - RKCCD-OOS/IOS, UNCCD-OOS/IOS: wrapper output is already median/MADN
#     standardized with larger = more outlying (no inversion needed).
#   - LOF: max LOF across MinPts 11:30, larger = more outlying (no inversion;
#     matches simulations/.../LOF/LOF.R and Real_Data_LOF.R, Thresh = 1.5).
#   - DBSCAN, MST, ODIN: the original baseline code only produces a binary
#     0/1 cluster-membership label (0 = outlier), not a continuous score. We
#     binarize: score = 1 if flagged outlier by the original algorithm, else
#     0; matched threshold = 0.5. This reproduces the exact TPR/TNR/BA/F2 the
#     original count_DBSCAN/count_MST2/count_ODIN functions would compute.
#   - iForest: rather than call the ISO.R wrapper (which internally
#     binarizes at threshold=0.55 and discards the raw score), we call
#     isotree::isolation.forest/predict directly with the SAME
#     hyperparameters (ntrees=1000, sample_size=256) to expose a genuine
#     continuous anomaly score (larger = more outlying), and apply the same
#     threshold = 0.55 externally at evaluate() time.
#
# t_construct is populated ONLY for the two CCD-based methods (it times the
# digraph-construction call in isolation: RKCCD_correct_quant() /
# nnccd_clustering_quantile()); it is NA for the five baselines, which have
# no analogous "construction" phase.
#
# ALPHA RESOLUTION IN THESE FOUR WRAPPERS (changed 2026-09-05). They used to
# call get_simul(variant, d) and take whatever the harness-local bucket
# returned. `quant` now has no default, so each names a resolver:
#   * RKCCD-OOS/IOS -> rk_quant_label_paper(d). This CHANGES BEHAVIOUR AT
#     d = 10 ONLY: the old rk_quant_for_d() bucket switched at d <= 10 and
#     returned "99" there, the manuscript switches at d < 10 and gives "999".
#     Every other dimension is unaffected. Any cached d = 10 RKCCD result
#     predates the change.
#   * UNCCD-OOS/IOS -> nn_quant_label_paper_UN(d). These wrappers are generic
#     UN-CCD constructions, not one of the paper's four named MCCD detectors,
#     so there is no SUN-specific schedule to apply; the UN-MCCD schedule is
#     the right one and is named explicitly here so the choice is visible.
#     It matches the old nn_quant_for_d() bucket exactly, so no behaviour
#     changes at any d.
# The `n = nrow(X)` argument is what makes a table too short for the data a
# hard error instead of a silent NA fill.

rkccd_oos_method <- function(X, d, Y = NULL, ...) {
  X <- as.matrix(X)
  tab <- get_simul("RK", d, quant = rk_quant_label_paper(d), n = nrow(X))
  t0 <- Sys.time()
  invisible(RKCCD_correct_quant(X, r.seq = 10, dom.method = "greedy2",
                                 quan = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE))
  t_construct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0b <- Sys.time()
  score <- RKCCD_OOS(datax = X, simul = tab$simul, d = d, quant = tab$quant)
  t_total <- as.numeric(difftime(Sys.time(), t0b, units = "secs"))

  list(score = score, t_construct = t_construct, t_total = t_total)
}

rkccd_ios_method <- function(X, d, Y = NULL, min.cls = 0, ...) {
  X <- as.matrix(X)
  tab <- get_simul("RK", d, quant = rk_quant_label_paper(d), n = nrow(X))
  t0 <- Sys.time()
  invisible(RKCCD_correct_quant(X, r.seq = 10, dom.method = "greedy2",
                                 quan = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE,
                                 min.cls = min.cls))
  t_construct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0b <- Sys.time()
  score <- RKCCD_IOS(datax = X, simul = tab$simul, d = d, quant = tab$quant, min.cls = min.cls)
  t_total <- as.numeric(difftime(Sys.time(), t0b, units = "secs"))

  list(score = score, t_construct = t_construct, t_total = t_total)
}

unccd_oos_method <- function(X, d, Y = NULL, method = "ascend", ...) {
  X <- as.matrix(X)
  tab <- get_simul("NN", d, quant = nn_quant_label_paper_UN(d), n = nrow(X))
  t0 <- Sys.time()
  invisible(nnccd_clustering_quantile(X, low.num = 3, quantile = "lower", method = method,
                                       dom.method = "greedy2", simul = tab$simul, niter = 1000, scores = TRUE))
  t_construct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0b <- Sys.time()
  score <- NNCCD_OOS(datax = X, simul = tab$simul, method = method, d = d)
  t_total <- as.numeric(difftime(Sys.time(), t0b, units = "secs"))

  list(score = score, t_construct = t_construct, t_total = t_total)
}

unccd_ios_method <- function(X, d, Y = NULL, method = "ascend", min.cls = 0, ...) {
  X <- as.matrix(X)
  tab <- get_simul("NN", d, quant = nn_quant_label_paper_UN(d), n = nrow(X))
  t0 <- Sys.time()
  invisible(nnccd_clustering_quantile(X, low.num = 3, quantile = "lower", method = method,
                                       dom.method = "greedy2", simul = tab$simul, niter = 1000, scores = TRUE))
  t_construct <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t0b <- Sys.time()
  score <- NNCCD_IOS(datax = X, simul = tab$simul, method = method, d = d, min.cls = min.cls)
  t_total <- as.numeric(difftime(Sys.time(), t0b, units = "secs"))

  list(score = score, t_construct = t_construct, t_total = t_total)
}

lof_method <- function(X, d, Y = NULL, ...) {
  X <- as.matrix(X)
  t0 <- Sys.time()
  score <- LOF(X, L_MinPts = 11, U_MinPts = 30)
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

dbscan_method <- function(X, d, Y = NULL, ...) {
  if (is.null(Y)) stop("dbscan_method(): requires Y to compute the oracle contamination level (matches Real_Data_DBSCAN.R)")
  X <- as.matrix(X)
  cont <- sum(Y == 0) / length(Y)
  t0 <- Sys.time()
  labels <- DBSCAN(X, k = 4, quant = cont)  # MinPts = 4, per Real_Data_DBSCAN.R
  score <- ifelse(labels == 0, 1, 0)
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

mst_method <- function(X, d, Y = NULL, cont = 0.02, thresh = 1.2, ...) {
  X <- as.matrix(X)
  t0 <- Sys.time()
  labels <- MST_Outlier(X, cont = cont, thresh = thresh)  # defaults per Real_Data_MST.R
  score <- ifelse(labels == 0, 1, 0)
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

odin_method <- function(X, d, Y = NULL, ...) {
  X <- as.matrix(X)
  t0 <- Sys.time()
  labels <- ODIN(X)  # default k = round(sqrt(n)), indegree_threshold = round(n^(1/3)), per Real_Data_ODIN.R
  score <- ifelse(labels == 0, 1, 0)
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

iforest_method <- function(X, d, Y = NULL, seed = 1, ...) {
  X <- as.matrix(X)
  t0 <- Sys.time()
  set.seed(seed)
  sample_size <- min(256, nrow(X))
  model <- isotree::isolation.forest(X, ntrees = 1000, sample_size = sample_size)
  score <- as.numeric(predict(model, X))
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

METHOD_REGISTRY <- list(
  "RKCCD-OOS" = rkccd_oos_method,
  "RKCCD-IOS" = rkccd_ios_method,
  "UNCCD-OOS" = unccd_oos_method,
  "UNCCD-IOS" = unccd_ios_method,
  "LOF"       = lof_method,
  "DBSCAN"    = dbscan_method,
  "MST"       = mst_method,
  "ODIN"      = odin_method,
  "iForest"   = iforest_method
)

# Real-data thresholds. OS methods: 2, per manuscript line ~1123
# ("The four OSs use a threshold of 2 ..."). LOF: 1.5 (default Thresh in
# Real_Data_LOF.R, used for WBC). DBSCAN/MST/ODIN: 0.5 against the binarized
# {0,1} score defined above. iForest: 0.55, per ISO.R's default `threshold`.
REAL_DATA_THRESHOLDS <- list(
  "RKCCD-OOS" = 2, "RKCCD-IOS" = 2, "UNCCD-OOS" = 2, "UNCCD-IOS" = 2,
  "LOF" = 1.5, "DBSCAN" = 0.5, "MST" = 0.5, "ODIN" = 0.5, "iForest" = 0.55
)

# ---------------------------------------------------------------------------
# 4. evaluate(): TPR/TNR/BA/F2 wrapper
# ---------------------------------------------------------------------------
#
# Y convention (R/general_functions/count.R, count_scores2): 1 = regular,
# 0 = outlier; and count_scores2 assumes the data is POSITIONALLY ordered
# with all regular points first and all outliers last (it slices
# label_pred[1:(n-n0)] and label_pred[(n-n0+1):n], it does not look up Y at
# each index beyond using it to count n0). RealData_Collection.R intends to
# sort every real dataset outliers-last ("# move outliers to the end"), and
# the synthetic generator idiom (rbind(cluster1, cluster2, outliers)) does so
# by construction -- BUT the loader's glass block sorts by glass[,9] (the 9th
# FEATURE, Fe; glass has 10 columns with the label in column 10), so glass
# comes back with its outliers in the middle rather than last, and
# count_scores2 would silently miscount it. Guard: jointly reorder (Y, score)
# regulars-first (stable sort, so relative order within each class is
# preserved) before calling count_scores2. A no-op for already-sorted data.
#
# SCORE GUARDS (added 2026-09-05). count_scores2 initialises
#   label_pred <- rep(0, n)     # 0 == OUTLIER
# and then overwrites positions by `which(score >= threshold)` and
# `which(score < threshold)`. `which()` drops NA, so any NA/NaN score keeps
# its initial value and is silently predicted an OUTLIER; and if `score` is
# shorter than `Y`, every trailing position keeps its initial value and is
# also silently predicted an outlier. Both inflate TPR without any warning,
# which is the worst possible failure for a detector benchmark. Neither is
# hypothetical: a quantile table shorter than its data set (see
# check_simul_extent above) produces exactly the NA case. Refuse both.
evaluate <- function(Y, score, threshold) {
  if (length(score) != length(Y)) {
    stop(sprintf(paste0(
      "evaluate(): score/label length mismatch -- length(score) = %d, length(Y) = %d.\n",
      "count_scores2() initialises predictions to 0 (= outlier), so the %d ",
      "unmatched position(s) would be silently counted as detected outliers."),
      length(score), length(Y), abs(length(Y) - length(score))))
  }
  bad <- !is.finite(score)
  if (any(bad)) {
    stop(sprintf(paste0(
      "evaluate(): %d of %d score(s) are not finite (%d NA, %d NaN, %d Inf) ",
      "at position(s) %s%s.\n",
      "count_scores2() drops these from both which() calls, leaving them at ",
      "the initial value 0 (= outlier), so they would be silently counted as ",
      "detected outliers."),
      sum(bad), length(score),
      sum(is.na(score) & !is.nan(score)), sum(is.nan(score)),
      sum(is.infinite(score)),
      paste(utils::head(which(bad), 10), collapse = ", "),
      if (sum(bad) > 10) ", ..." else ""))
  }
  ord <- order(Y, decreasing = TRUE)  # Y: 1 = regular first, 0 = outlier last
  v <- count_scores2(Y[ord], score[ord], threshold)
  names(v) <- c("TPR", "TNR", "BA", "F2")
  v
}

# ---------------------------------------------------------------------------
# 5. CSV checkpointing
# ---------------------------------------------------------------------------

# Per-cell append is deliberate and must stay. G: is a USB-NVMe enclosure that
# has dropped off the bus under sustained write load; buffering a long run in
# memory and writing once at the end turns a drop into total loss, whereas
# appending per cell bounds it to the cells not yet written. The cost of that
# choice is that a drop can leave a half-written final line, which is what the
# has_result() truncation handling below exists for.

#' Append one result row to a CSV, creating it with a header if it doesn't
#' exist yet. `row` should be a named list/vector of scalars.
#'
#' HEADER CHECK (added 2026-09-05). The append branch writes with
#' col.names = FALSE, i.e. POSITIONALLY: whatever order `row`'s names happen
#' to be in is the order the values land in, with no reference to the header
#' already on disk. A caller that builds its row with fields in a different
#' order, or with a field added or dropped, silently produces a file whose
#' columns no longer mean what the header says. Compare names against the
#' existing header first and refuse the write if they differ.
append_result <- function(csv_path, row) {
  df <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(csv_path)) {
    write.csv(df, csv_path, row.names = FALSE)
    return(invisible(csv_path))
  }
  hdr_line <- readLines(csv_path, n = 1, warn = FALSE)
  if (length(hdr_line) == 1 && nzchar(hdr_line)) {
    hdr <- names(utils::read.csv(text = hdr_line, stringsAsFactors = FALSE,
                                 check.names = TRUE))
    got <- names(df)
    if (!identical(hdr, got)) {
      stop(sprintf(paste0(
        "append_result(): row names do not match the existing header of %s.\n",
        "  header on disk : %s\n",
        "  row being added: %s\n",
        "Rows are appended positionally (col.names = FALSE), so writing this ",
        "row would put values under the wrong column names."),
        csv_path, paste(hdr, collapse = ", "), paste(got, collapse = ", ")))
    }
  }
  utils::write.table(df, csv_path, sep = ",", col.names = FALSE, row.names = FALSE,
                     append = TRUE, qmethod = "double")
  invisible(csv_path)
}

#' Check whether a result matching `keys` (named list/vector, e.g.
#' c(dataset="WBC", method="LOF", seed="1")) already exists in csv_path, for
#' skip-if-present restart logic.
#'
#' TRUNCATION HANDLING (added 2026-09-05). read.csv() defaults to fill = TRUE,
#' so a line cut short by a drive drop mid-append comes back as a COMPLETE row
#' whose missing tail is NA. If the key columns survived the cut, that row
#' matches, has_result() returns TRUE, and the restart skips a cell whose
#' result was never actually written. And if a key column itself is the NA,
#' the comparison yields NA, `any()` returns NA, and `if (has_result(...))`
#' errors out with "missing value where TRUE/FALSE needed".
#'
#' So: count fields per line first and treat any line whose field count
#' differs from the header's as absent (warning, with line numbers); treat a
#' key-matching row that is NA in any non-key column as absent too, since its
#' payload is what the caller wanted; and wrap the answer in isTRUE() so this
#' function can only ever return TRUE or FALSE.
has_result <- function(csv_path, keys) {
  if (!file.exists(csv_path)) return(FALSE)

  nf <- tryCatch(utils::count.fields(csv_path, sep = ",", quote = "\""),
                 error = function(e) NULL)
  if (is.null(nf) || length(nf) < 2) return(FALSE)
  n_hdr <- nf[1]
  nf_data <- nf[-1]
  # NA means count.fields could not parse the line at all (e.g. an unbalanced
  # quote from a cut-off write); that is a truncation too.
  ragged <- which(is.na(nf_data) | nf_data != n_hdr)  # 1 = first data row
  if (length(ragged) > 0) {
    warning(sprintf(paste0(
      "has_result(): %s has %d line(s) with the wrong number of fields ",
      "(expected %d) at file line(s) %s -- most likely truncated by an ",
      "interrupted append. Treating them as absent, so the affected cells ",
      "will be recomputed."),
      csv_path, length(ragged), n_hdr,
      paste(ragged + 1L, collapse = ", ")), call. = FALSE)
  }

  df <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  for (k in names(keys)) {
    if (!(k %in% names(df))) return(FALSE)
  }

  ok <- rep(TRUE, nrow(df))
  if (length(ragged) > 0) ok[ragged[ragged <= nrow(df)]] <- FALSE

  match_all <- ok
  for (k in names(keys)) {
    col <- as.character(df[[k]])
    match_all <- match_all & !is.na(col) & (col == as.character(keys[[k]]))
  }
  match_all[is.na(match_all)] <- FALSE

  # A row whose keys match but whose payload columns are NA is a partial
  # write, not a result.
  payload <- setdiff(names(df), names(keys))
  if (length(payload) > 0 && any(match_all)) {
    incomplete <- apply(df[match_all, payload, drop = FALSE], 1,
                        function(r) any(is.na(r)))
    match_all[which(match_all)[incomplete]] <- FALSE
  }

  isTRUE(any(match_all))
}

# ---------------------------------------------------------------------------
# 6. Real-data loader (WBC and friends)
# ---------------------------------------------------------------------------
#
# data/outlier_detection/RealData_Collection.R does setwd() to
# data/outlier_detection at its top; we capture and restore the working
# directory so later here::here()-based / relative-path code in callers is
# unaffected. All objects it creates (WBC, glass, hepatitis, ...) land in
# .GlobalEnv per source()'s default `local = FALSE`.
load_real_dataset <- function(name) {
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  if (!exists(name, envir = .GlobalEnv)) {
    source(here::here("data/outlier_detection/RealData_Collection.R"))
  }
  df <- get(name, envir = .GlobalEnv)
  d <- ncol(df) - 1
  X <- as.data.frame(df[, 1:d, drop = FALSE])
  Y <- df[, d + 1]
  list(X = X, Y = Y, d = d, n = nrow(df))
}
