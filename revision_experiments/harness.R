#!/usr/bin/env Rscript
# revision_experiments/harness.R
#
# Sourceable experiment harness for the PR-D-26-05767 revision-experiments
# pipeline (Task T2). Provides:
#   - get_simul(variant, d, quant)          quantile-table loader
#   - defensive index-clamp overrides for ccd.Kest.edge.quantile / nnccd.radi
#   - METHOD_REGISTRY: 9 named scoring functions (4 CCD-based OS + 5 baselines)
#   - evaluate(Y, score, threshold)          TPR/TNR/BA/F2 wrapper
#   - append_result()/has_result()           CSV checkpointing for restarts
#
# This file only READS the original R/, methods/, simulations/, data/ trees.
# Nothing in those directories is modified. Where a defensive tweak is
# required (index clamps against n > 5000), the function is redefined here,
# in the global environment, AFTER the originals are sourced, so the
# redefinition wins for every subsequent call (R resolves free variables by
# looking up the *current* binding in the enclosing environment, which for
# all of these functions is .GlobalEnv, since methods/outlyingness_scores/*
# source their dependencies with `source(here::here(...))` at top level).

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
# 2. Defensive index-clamp overrides (redefined AFTER sourcing originals)
# ---------------------------------------------------------------------------
#
# ccd.Kest.edge.quantile() (R/ccds/RK_CCD_New.R, ~line 184) indexes
# Kest.slopes$quan[[as.character(quan)]][j, ] by the running point-count j,
# which ranges up to n = nrow(dx). The recovered RK quantile tables have
# exactly 5000 rows, so any dataset with n > 5000 would index past the end of
# the matrix and R would throw "subscript out of bounds". None of our smoke
# datasets are anywhere near that size, but per the task spec we clamp
# defensively: j_idx = min(j, nrow(table)), reusing the table's last row for
# any j beyond its extent. Everything else is byte-for-byte identical to the
# original (including the pre-existing, unchanged "avoid 0 radius" logic in
# the scores=TRUE branch) -- we are not "fixing" that logic, only guarding
# the index.
ccd.Kest.edge.quantile <- function(dx, ddx, low.num, method = "non-dynamic", r.seq, quan, simul = NULL, niter, scores = F) {

  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0, n)

  if (!is.null(simul)) {
    Kest.slopes <- simul
  } else if (method == "non-dynamic") {
    Kest.slopes <- Kest.simpois.edge.quantile(n, d, r.seq, quan, niter)
  } else {
    Kest.slopes <- Kest.simpois.edge.quantile.dynamic(n, d, r.seq, quan, niter)
  }
  rr <- Kest.slopes$r

  quan_tab <- Kest.slopes$quan[[as.character(quan)]]
  tab_rows <- nrow(quan_tab)

  if (!scores) {
    for (i in 1:n) {
      o.d <- order(ddx[i, ])
      for (j in low.num:n) {
        r <- ddx[i, o.d[j]] * rr
        sc <- ddx[i, o.d[j]]
        Kest.obs <- Kest.f.edge(ddx[o.d[1:j], o.d[1:j]], r, sc, d)

        j_idx <- min(j, tab_rows)  # <-- defensive clamp (only change vs. original)
        lo.hi <- quan_tab[j_idx, ]
        flag <- (Kest.obs > lo.hi)
        if (any(flag)) {
          if (j == low.num) R[i] <- 0
          else R[i] <- ddx[i, o.d[j - 1]]
          break
        } else if (j == n) {
          R[i] <- 0
          break
        }
      }
    }
  } else {
    for (i in 1:n) {
      o.d <- order(ddx[i, ])
      for (j in low.num:n) {
        r <- ddx[i, o.d[j]] * rr
        sc <- ddx[i, o.d[j]]
        Kest.obs <- Kest.f.edge(ddx[o.d[1:j], o.d[1:j]], r, sc, d)

        j_idx <- min(j, tab_rows)  # <-- defensive clamp (only change vs. original)
        lo.hi <- quan_tab[j_idx, ]
        flag <- (Kest.obs > lo.hi)
        if (any(flag)) {
          if (j == low.num) R[i] <- 0
          else R[i] <- ddx[i, o.d[j - 1]]
          break
        } else if (j == n) {
          R[i] <- 0
          break
        }
        if (R[i] == 0) { R[i] = sort(ddx[i, ])[2] }  # avoid 0 radius (necessary for outlyingness scores!) -- unchanged
      }
    }
  }
  return(list(R = R, KS = NULL))
}

# nnccd.radi() (R/ccds/UN_CCD.R, ~line 241) slices simul$average[1:n] /
# simul$median[1:n] and then indexes those slices by j-1 (or j) up to n.
# The recovered NN tables have length-5000 vectors; simul$average[1:n] with
# n > 5000 would silently pad with NA (R does not error on this, it just
# returns NA for out-of-range positions), which would then propagate NAs
# into every downstream radius/score. We clamp the index into the raw
# simul$average / simul$median vectors at length(simul$average) so an
# oversized n reuses the table's last entry instead of introducing NAs.
# Everything else is byte-for-byte identical to the original.
nnccd.radi <- function(dx, quantile = "lower", method = "ascend", low.num, quant, simul = NULL, niter, scores = F) {

  ddx <- as.matrix(dist(dx))
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0, n)

  if (quantile == "lower") {
    if (!is.null(simul)) {
      # <-- defensive clamp (only change vs. original): pmin(1:n, length(...))
      #     instead of a bare 1:n slice, so indices beyond the table length
      #     reuse the last available table entry rather than becoming NA.
      avg_len <- length(simul$average)
      med_len <- length(simul$median)
      NN.envelop <- list(
        average = simul$average[pmin(1:n, avg_len)],
        median  = simul$median[pmin(1:n, med_len)]
      )
    } else {
      NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)
    }
    if (!scores) {
      for (i in 1:n) {
        if (method == "ascend") {
          o.d <- order(ddx[i, ])
          for (j in low.num:n) {
            r <- ddx[i, o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j], o.d[2:j]], r)

            lower.bound.ave = NN.envelop$average[j - 1]
            lower.bound.med = NN.envelop$median[j - 1]
            if (NN.dist.obs$averge < lower.bound.ave | NN.dist.obs$median < lower.bound.med) {
              if (j == low.num) R[i] = 0
              else R[i] = ddx[i, o.d[j - 1]]
              break
            }
          }
        }
        if (method == "descend") {
          o.d <- order(ddx[i, ], decreasing = T)
          for (j in 1:(n - low.num)) {
            r <- ddx[i, o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n - 1)], o.d[j:(n - 1)]], r)

            lower.bound.ave = rev(NN.envelop$average)[j + 2]
            lower.bound.med = rev(NN.envelop$median)[j + 2]
            if (NN.dist.obs$averge > lower.bound.ave & NN.dist.obs$median > lower.bound.med) {
              R[i] = r
              break
            }
          }
        }
      }
    } else {
      for (i in 1:n) {
        if (method == "ascend") {
          o.d <- order(ddx[i, ])
          for (j in low.num:n) {
            r <- ddx[i, o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j], o.d[2:j]], r)

            lower.bound.ave = NN.envelop$average[j - 1]
            lower.bound.med = NN.envelop$median[j - 1]
            if (NN.dist.obs$averge < lower.bound.ave | NN.dist.obs$median < lower.bound.med) {
              if (j == low.num) R[i] = 0
              else R[i] = ddx[i, o.d[j - 1]]
              break
            }
          }
        }
        if (method == "descend") {
          o.d <- order(ddx[i, ], decreasing = T)
          for (j in 1:(n - low.num)) {
            r <- ddx[i, o.d[j]]
            NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n - 1)], o.d[j:(n - 1)]], r)

            lower.bound.ave = rev(NN.envelop$average)[j + 2]
            lower.bound.med = rev(NN.envelop$median)[j + 2]
            if (NN.dist.obs$averge > lower.bound.ave & NN.dist.obs$median > lower.bound.med) {
              R[i] = r
              break
            }
          }
        }
        if (R[i] == 0) { R[i] = sort(ddx[i, ])[2] }  # avoid 0 radius (necessary for outlyingness scores!) -- unchanged
      }
    }
  }
  return(list(R = R, KS = NULL))
}

# ---------------------------------------------------------------------------
# 3. get_simul(): quantile-table loader
# ---------------------------------------------------------------------------
#
# Quantile-per-dimension convention settled by reading the actual (non-cutoff)
# simulation scripts that generated the published Section 3.1 tables:
#
#   RK-CCD (RKCCD_OOS_IOS/Simulation/Gaussian/{5,10}d/*.R and
#   Uniform_cutoffs/RKCCD_OOS_IOS_cutoff.R, identical rule at every d block):
#       if (d <= 5) quant = 0.99 else quant = 0.999
#   This is a genuine two-bucket boolean rule (not merely enumerated cases),
#   so it generalizes cleanly to untested d, including WBC's d = 9 -> 0.999.
#
#   UN-CCD / NN-CCD (UNCCD_OOS_IOS/Simulation/Gaussian/{5,10,20}d/*.R):
#       d=2 -> table "85%", d=3 -> "90%", d=5 -> "95%", d=10 -> "99%",
#       d=20 -> "999%"
#   There is no single boolean rule here (the schedule has 5 anchor points,
#   and the auxiliary Uniform_cutoffs/NNCCD_OOS_IOS_cutoff.R script actually
#   disagrees with the real Gaussian d=10 script -- cutoff says 0.999, the
#   real simulation script says 0.99 with table NN-test-simul_10d_99%.RData;
#   we trust the real simulation script since that is what produced the
#   published numbers). For untested d we use nearest-lower-anchor bucketing:
#       d==2 -> 85%; d in [3,4] -> 90%; d in [5,9] -> 95%;
#       d in [10,19] -> 99%; d>=20 -> 999%
#   This puts WBC's d=9 at 95%. Both the RK and NN conventions are validated
#   by the Part-C reproduction gate in 03_smoke.R.

RK_QUANT_TABLE_DIR <- here::here("R/RK-test_quantile")
NN_QUANT_TABLE_DIR <- here::here("R/NN-test_quantile")

rk_quant_for_d <- function(d) if (d <= 5) "99" else "999"

nn_quant_for_d <- function(d) {
  if (d <= 2) "85"
  else if (d <= 4) "90"
  else if (d <= 9) "95"
  else if (d <= 19) "99"
  else "999"
}

#' Resolve and load the RK-CCD or NN(UN)-CCD quantile lookup table for a
#' given dimension.
#'
#' @param variant "RK" or "NN"
#' @param d dimensionality
#' @param quant optional override, e.g. "99", "999", "95" (as it appears in
#'   the filename, without the trailing "%"). If NULL, resolved via
#'   rk_quant_for_d()/nn_quant_for_d().
#' @return list(simul = <loaded object>, quant = <numeric quantile, e.g.
#'   0.999>, quant_label = <string used in filename>, file = <path used>)
get_simul <- function(variant = c("RK", "NN"), d, quant = NULL) {
  variant <- match.arg(variant)
  if (variant == "RK") {
    q <- if (is.null(quant)) rk_quant_for_d(d) else quant
    fname <- sprintf("RK-test-simul_%dd_%s%%.RData", d, q)
    path <- file.path(RK_QUANT_TABLE_DIR, fname)
  } else {
    q <- if (is.null(quant)) nn_quant_for_d(d) else quant
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
  list(
    simul = get("simul", envir = e),
    quant = as.numeric(paste0("0.", q)),
    quant_label = q,
    file = path
  )
}

# ---------------------------------------------------------------------------
# 4. Method registry
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

rkccd_oos_method <- function(X, d, Y = NULL, ...) {
  X <- as.matrix(X)
  tab <- get_simul("RK", d)
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
  tab <- get_simul("RK", d)
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
  tab <- get_simul("NN", d)
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
  tab <- get_simul("NN", d)
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
# 5. evaluate(): TPR/TNR/BA/F2 wrapper
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
evaluate <- function(Y, score, threshold) {
  ord <- order(Y, decreasing = TRUE)  # Y: 1 = regular first, 0 = outlier last
  v <- count_scores2(Y[ord], score[ord], threshold)
  names(v) <- c("TPR", "TNR", "BA", "F2")
  v
}

# ---------------------------------------------------------------------------
# 6. CSV checkpointing
# ---------------------------------------------------------------------------

#' Append one result row to a CSV, creating it with a header if it doesn't
#' exist yet. `row` should be a named list/vector of scalars.
append_result <- function(csv_path, row) {
  df <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(csv_path)) {
    write.csv(df, csv_path, row.names = FALSE)
  } else {
    write.table(df, csv_path, sep = ",", col.names = FALSE, row.names = FALSE,
                append = TRUE, qmethod = "double")
  }
}

#' Check whether a result matching `keys` (named list/vector, e.g.
#' c(dataset="WBC", method="LOF", seed="1")) already exists in csv_path, for
#' skip-if-present restart logic.
has_result <- function(csv_path, keys) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  for (k in names(keys)) {
    if (!(k %in% names(df))) return(FALSE)
  }
  match_all <- rep(TRUE, nrow(df))
  for (k in names(keys)) {
    match_all <- match_all & (as.character(df[[k]]) == as.character(keys[[k]]))
  }
  any(match_all)
}

# ---------------------------------------------------------------------------
# 7. Real-data loader (WBC and friends)
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
