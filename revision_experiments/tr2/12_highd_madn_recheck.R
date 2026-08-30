#!/usr/bin/env Rscript
# revision_experiments/12_highd_madn_recheck.R
#
# WP8 audit follow-up (FINDINGS.md line 208, "OPEN QUESTION worth one cheap
# experiment"). The std_MADN bug -- mad(x)==0 jumped straight to an all-zero
# vector with NO SD-fallback step -- suppressed OOS *precisely* in the
# degenerate regime (all-zero -> the tie-break dies on an out-of-bounds
# neighbour -> OOS collapses to rank-AUC 0.5). The paper's high-dimensional
# *simulation* cells (d = 50, 100) were computed on that same buggy path and
# were NEVER cached (only the printed TPR/TNR survive in slurm-*.out), so the
# fingerprint scan (scan_madn0.R, which needs cached .rds scores) does not
# cover them. If any of those cells hit MADN=0 for OOS, the published
# "IOS > OOS at high d" comparison could be partly an artifact of OOS being
# zeroed.
#
# THIS SCRIPT settles it by MEASUREMENT, not inference. For each requested
# high-d cell it regenerates the IDENTICAL Monte Carlo data (seed 123, gen
# code copied verbatim from the first-cycle drivers), builds each digraph
# ONCE, and applies BOTH std_MADN variants -- the current fixed one and a
# byte-faithful copy of the old buggy one -- to the same graph. The digraph
# and covering radii are independent of std_MADN (it only rescales the final
# scores), so a single build yields OOS and IOS under both paths, and the
# fixed-vs-buggy difference is EXACT per replicate. It records, per cell:
#   * whether mad(rawOOS)==0 fired (the whole-vector OOS standardization), and
#     for IOS how many per-cluster mad==0 events fired and how many points
#     they cover;
#   * rank-AUC (threshold-free) and the paper's TPR/TNR/BA/F2 at the published
#     cutoff, for OOS and IOS under fixed and buggy std_MADN.
# If fixed == buggy for every rep, the fix is a no-op at that cell and the
# published comparison is real (not a standardization artifact).
#
# FIDELITY GUARD: before the production loop, an equivalence assertion checks
# that the copied post-graph logic (oos_apply / ios_apply with the fixed
# madn) reproduces the real RKCCD_OOS/IOS or NNCCD_OOS/IOS to < 1e-8 on the
# first dataset. If the copy has drifted from the shipped functions the run
# hard-stops -- no silent infidelity.
#
# CLI:
#   Rscript 12_highd_madn_recheck.R <family> <dist> <d> <n> [reps] [cores]
#     family : RKCCD | UNCCD
#     dist   : Gaussian | Uniform
#   e.g. Rscript 12_highd_madn_recheck.R RKCCD Gaussian 100 500 100 22
#
# Output (appended): results/highd_madn_recheck.csv  (one row per cell)
#         plus a per-rep dump results/highd_madn_recheck_<cell>_perrep.csv
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(here)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(MASS)
  library(igraph)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: <family> <dist> <d> <n> [reps] [cores]")
FAMILY <- args[[1]]                       # RKCCD | UNCCD
DIST   <- args[[2]]                       # Gaussian | Uniform
D      <- as.integer(args[[3]])
N      <- as.integer(args[[4]])
REPS   <- if (length(args) >= 5) as.integer(args[[5]]) else 100L
CORES  <- if (length(args) >= 6) as.integer(args[[6]]) else 22L
# arg 7: skip the fidelity guard (only after a family has already been
# validated once -- the copied post-graph logic is family-independent and the
# build_graph reproduction is proven per family). Saves ~3 serial builds/cell,
# which matters for UNCCD (~minutes/build at n=500).
SKIP_FIDELITY <- if (length(args) >= 7) as.logical(args[[7]]) else FALSE
stopifnot(FAMILY %in% c("RKCCD","UNCCD"), DIST %in% c("Gaussian","Uniform"))

CELL <- sprintf("%s_%s_d%d_n%d", FAMILY, DIST, D, N)
message(sprintf("[recheck] cell=%s reps=%d cores=%d  %s",
                CELL, REPS, CORES, format(Sys.time())))

## ---- sources --------------------------------------------------------------
source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))  # fixed std_MADN, OOS, IOS, Vic_Den
source(here::here("R/general_functions/count.R"))                        # count_scores
if (FAMILY == "RKCCD") {
  source(here::here("R/ccds/RK_CCD_New.R"))
  source(here::here("methods/outlyingness_scores/RKCCD_OOS_IOS.R"))
} else {
  source(here::here("R/ccds/UN_CCD.R"))
  source(here::here("methods/outlyingness_scores/UNCCD_OOS_IOS.R"))
}

## ---- the two std_MADN variants -------------------------------------------
# fixed: MADN -> SD -> 0 (the shipped Outlyingness_Score.R std_MADN)
std_MADN_fixed <- std_MADN
# buggy: the pre-fix code -- mad==0 jumps straight to all-zeros, no SD step
std_MADN_buggy <- function(x){
  if (mad(x) != 0) {
    s <- (x - median(x)) / mad(x)
  } else {
    s <- rep(0, length(x))
  }
  return(s)
}

## ---- table + thresholds + method params for this cell ---------------------
# quant / table / Threshold vectors, matched to the first-cycle drivers.
# Dimension index: d = 2,3,5,10,20,50,100 -> 1..7.
d_idx <- match(D, c(2,3,5,10,20,50,100))
if (is.na(d_idx)) stop("unsupported d (expect one of 2,3,5,10,20,50,100)")

if (FAMILY == "RKCCD") {
  load(here::here(sprintf("R/RK-test_quantile/RK-test-simul_%dd_999%%.RData", D)))  # -> simul
  QUANT  <- if (D >= 10) 0.999 else 0.99   # matches the driver if/else chain
  METHOD <- NA_character_
  thr_src <- sprintf("simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/%s/Threshold.R", DIST)
} else {
  load(here::here(sprintf("R/NN-test_quantile/NN-test-simul_%dd_999%%.RData", D)))  # -> simul
  QUANT  <- c(0.85,0.9,0.95,0.99,0.999,0.999,0.999)[d_idx]  # UNCCD driver chain
  METHOD <- "descend"
  thr_src <- sprintf("simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/%s/Threshold.R", DIST)
}
source(here::here(thr_src))                # -> Threshold_OOS, Threshold_IOS
THR_OOS <- Threshold_OOS[d_idx]
THR_IOS <- Threshold_IOS[d_idx]
MIN_CLS <- 0                               # both families' drivers use min.cls=0 (Gaussian/Uniform)

message(sprintf("[recheck] quant=%s method=%s thrOOS=%s thrIOS=%s",
                QUANT, ifelse(is.na(METHOD),"-",METHOD), THR_OOS, THR_IOS))

## ---- data generation (verbatim from the first-cycle drivers) -------------
gen_data <- function(dist, d, n, cont = 0.05, iteN = 100L, seed = 123L,
                     cls_dis = 3, otl_dis = 2, r_min = 0.7, r_max = 1.3) {
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- apply(rbind(mu1, mu2), 2, mean)
  n1  <- round(n * (1 - cont) * 0.5)
  n2  <- round(n * (1 - cont) * 0.5) - 1
  n0  <- round(n * cont)
  set.seed(seed)
  if (dist == "Gaussian") {
    noise_level <- 0.01
    sigma <- 1 / sqrt(qchisq(1 - noise_level, d))
    dl <- lapply(1:iteN, function(x) {
      data1 <- mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
      data2 <- mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
      i <- 0; outlier <- NULL
      while (i < n0) {
        temp <- rpoisball.unit(1, d) * 5 + mu
        if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
          outlier <- rbind(outlier, temp); i <- i + 1
        }
      }
      rownames(outlier) <- NULL
      rbind(data1, data2, outlier)
    })
  } else { # Uniform
    dl <- lapply(1:iteN, function(x) {
      data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
      data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
      i <- 0; outlier <- NULL
      while (i < n0) {
        temp <- rpoisball.unit(1, d) * 5 + mu
        if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
          outlier <- rbind(outlier, temp); i <- i + 1
        }
      }
      rownames(outlier) <- NULL
      rbind(data1, data2, outlier)
    })
  }
  list(data = dl, n0 = n0)
}

## ---- faithful post-graph copies, parameterized by the madn function -------
# OOS from a prebuilt radius vector R.  Byte-for-byte the post-graph block of
# RKCCD_OOS / NNCCD_OOS, with std_MADN(...) -> madn(...).
oos_apply <- function(datax, R, d, madn) {
  scores <- OOS(datax, R, d)
  scores <- madn(scores)
  vd <- Vic_Den(datax, R, d)
  score_vd <- cbind(scores, vd)
  order_index <- order(scores)
  score_vd <- score_vd[order_index, ]
  frequency <- table(scores)
  repeated_values <- as.numeric(names(frequency[frequency > 1]))
  for (x in repeated_values) {
    if (x == 0) index <- which(score_vd[, 1] == x)
    else        index <- which(abs(score_vd[, 1] - x) / x < 0.0001)
    if (length(index) < 2) break
    index_max <- max(index); index_min <- min(index)
    if (index_min == 1) {
      s <- score_vd[, 1][1]; e <- score_vd[, 1][index_max + 1]; diff <- e - s
    } else if (index_max == length(datax[, 1])) {
      s <- score_vd[, 1][index_min - 1]; e <- score_vd[, 1][index_max]; diff <- e - s
    } else {
      s <- score_vd[, 1][index_min - 1]; e <- score_vd[, 1][index_max + 1]; diff <- e - s
    }
    if (is.na(diff)) break
    den_sum <- sum(score_vd[, 2][index])
    score_vd[, 1][index] <- sapply(index, function(x) s + diff * score_vd[, 2][x] / den_sum)
    scores <- score_vd[order(order_index), 1]
  }
  scores
}

# IOS from a prebuilt radius vector R + cluster labels.  Byte-for-byte the
# post-graph block of RKCCD_IOS / NNCCD_IOS.
ios_apply <- function(datax, R, label, n_cls, d, madn) {
  vd <- Vic_Den(datax, R, d)
  member <- unique(label)
  scores <- lapply(member, function(x) {
    index_cls <- which(label == x)
    cluster <- datax[index_cls, ]
    R_cls <- R[index_cls]
    madn(IOS(cluster, R_cls, d))
  })
  scores_whole <- rep(0, length(datax[, 1]))
  for (i in 1:n_cls) scores_whole[which(label == member[i])] <- scores[[i]]
  score_vd <- cbind(scores_whole, vd)
  order_index <- order(scores_whole)
  score_vd <- score_vd[order_index, ]
  frequency <- table(scores_whole)
  repeated_values <- as.numeric(names(frequency[frequency > 1]))
  for (x in repeated_values) {
    if (x == 0) index <- which(score_vd[, 1] == x)
    else        index <- which(abs(score_vd[, 1] - x) / x < 0.0001)
    if (length(index) < 2) break
    index_max <- max(index); index_min <- min(index)
    if (index_min == 1) {
      s <- score_vd[, 1][1]; e <- score_vd[, 1][index_max + 1]; diff <- e - s
    } else if (index_max == length(datax[, 1])) {
      s <- score_vd[, 1][index_min - 1]; e <- score_vd[, 1][index_max]; diff <- e - s
    } else {
      s <- score_vd[, 1][index_min - 1]; e <- score_vd[, 1][index_max + 1]; diff <- e - s
    }
    if (is.na(diff)) break
    den_sum <- sum(score_vd[, 2][index])
    score_vd[, 1][index] <- sapply(index, function(x) s + diff * score_vd[, 2][x] / den_sum)
    scores_whole <- score_vd[order(order_index), 1]
  }
  scores_whole
}

# build a digraph ONCE and return R, cluster labels, si.ind
build_graph <- function(datax) {
  if (FAMILY == "RKCCD") {
    g <- RKCCD_correct_quant(datax, r.seq = 10, dom.method = "greedy2",
                             quan = QUANT, simul = simul, niter = 1000,
                             scores = TRUE, min.cls = MIN_CLS)
    list(R = g$R[order(g$D)], label = g$label, n_cls = g$si.ind)
  } else {
    g <- nnccd_clustering_quantile(datax, low.num = 3, quantile = "lower",
                                   method = METHOD, dom.method = "greedy2",
                                   simul = simul, niter = 1000, scores = TRUE)
    sil <- nnccd.silhouette(g, datax, cls = NULL, min.cls = MIN_CLS,
                            ind = NULL, lenDlimit = Inf)
    list(R = g$R[order(g$D)], label = sil$label, n_cls = sil$si.ind)
  }
}

## ---- rank-AUC (threshold-free), outlier = last n0, higher score = outlier -
auc_out <- function(scores, n, n0) {
  if (any(is.na(scores))) return(NA_real_)   # NaN/NA unrankable; +Inf is fine
  y <- c(rep(0L, n - n0), rep(1L, n0))
  r <- rank(scores, ties.method = "average")
  (sum(r[y == 1]) - n0 * (n0 + 1) / 2) / (as.numeric(n0) * (n - n0))
}

metrics <- function(scores, n, n0, thr) {
  cs <- count_scores(1, list(scores), thr, n, n0)   # TPR,TNR,BA,F2
  c(AUC = auc_out(scores, n, n0), cs)
}

## ---- generate the data ----------------------------------------------------
gd <- gen_data(DIST, D, N, cont = 0.05, iteN = REPS, seed = 123L)
data.list <- gd$data
n0 <- gd$n0
n  <- N
message(sprintf("[recheck] generated %d datasets, n=%d n0=%d d=%d (dist=%s)",
                length(data.list), n, n0, D, DIST))

## ---- FIDELITY GUARD: copied logic must reproduce the shipped functions ----
if (SKIP_FIDELITY) {
  message("[recheck] fidelity check SKIPPED (family previously validated)")
} else {
  message("[recheck] fidelity check against shipped functions on rep 1 ...")
  d1 <- data.list[[1]]
  g1 <- build_graph(d1)
  mine_oos <- oos_apply(d1, g1$R, D, std_MADN_fixed)
  mine_ios <- ios_apply(d1, g1$R, g1$label, g1$n_cls, D, std_MADN_fixed)
  if (FAMILY == "RKCCD") {
    real_oos <- RKCCD_OOS(d1, simul = simul, d = D, quant = QUANT)
    real_ios <- RKCCD_IOS(d1, simul = simul, d = D, quant = QUANT, min.cls = MIN_CLS)
  } else {
    real_oos <- NNCCD_OOS(d1, simul = simul, method = METHOD, d = D)
    real_ios <- NNCCD_IOS(d1, simul = simul, method = METHOD, d = D, min.cls = MIN_CLS)
  }
  do <- max(abs(mine_oos - real_oos)); di <- max(abs(mine_ios - real_ios))
  message(sprintf("[recheck]   OOS max|mine-real| = %.3e ; IOS = %.3e", do, di))
  if (!(is.finite(do) && do < 1e-8)) stop("FIDELITY FAIL: oos_apply != real OOS")
  if (!(is.finite(di) && di < 1e-8)) stop("FIDELITY FAIL: ios_apply != real IOS")
  message("[recheck]   fidelity OK -- copied post-graph logic matches shipped functions")
}

## ---- production loop (parallel over reps) ---------------------------------
cl <- makeCluster(CORES)
registerDoParallel(cl)
clusterExport(cl, c("FAMILY","DIST","D","N","QUANT","METHOD","MIN_CLS","simul",
                    "THR_OOS","THR_IOS",
                    "std_MADN_fixed","std_MADN_buggy","oos_apply","ios_apply",
                    "build_graph","auc_out","metrics","n0","n"),
              envir = environment())
invisible(clusterEvalQ(cl, {
  suppressWarnings(suppressMessages({library(here); library(MASS); library(igraph)}))
  source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))
  source(here::here("R/general_functions/count.R"))
  if (FAMILY == "RKCCD") {
    source(here::here("R/ccds/RK_CCD_New.R"))
    source(here::here("methods/outlyingness_scores/RKCCD_OOS_IOS.R"))
  } else {
    source(here::here("R/ccds/UN_CCD.R"))
    source(here::here("methods/outlyingness_scores/UNCCD_OOS_IOS.R"))
  }
  NULL
}))

per <- foreach(rep = seq_along(data.list), .combine = rbind,
               .packages = c("MASS","igraph","here")) %dopar% {
  out <- tryCatch({
    datax <- data.list[[rep]]
    # Use the ACTUAL row count, not the nominal N: the generator makes
    # n1+n2+n0 points and 2*round(N*0.475)-1+round(N*0.05) equals N only for
    # some N (e.g. 500) -- at N=1000 it is 999. The outliers are always the
    # last n0 rows, so labelling by (n_rep - n0) is correct; passing nominal N
    # to auc_out mismatched the score-vector length and NA'd every rep.
    n_rep <- nrow(datax)
    g <- build_graph(datax)
    R <- g$R
    # --- MADN=0 fingerprints ---
    raw_oos <- OOS(datax, R, D)
    madn0_oos <- as.integer(isTRUE(mad(raw_oos) == 0))
    # per-cluster IOS mad==0 events
    member <- unique(g$label)
    ncl_mad0 <- 0L; npts_mad0 <- 0L
    for (m in member) {
      idx <- which(g$label == m)
      raw_cls <- IOS(matrix(datax[idx, ], ncol = D), R[idx], D)
      if (isTRUE(mad(raw_cls) == 0)) { ncl_mad0 <- ncl_mad0 + 1L; npts_mad0 <- npts_mad0 + length(idx) }
    }
    # --- scores under both std_MADN paths (same graph) ---
    oos_fx <- oos_apply(datax, R, D, std_MADN_fixed)
    oos_bg <- oos_apply(datax, R, D, std_MADN_buggy)
    ios_fx <- ios_apply(datax, R, g$label, g$n_cls, D, std_MADN_fixed)
    ios_bg <- ios_apply(datax, R, g$label, g$n_cls, D, std_MADN_buggy)
    m_oos_fx <- metrics(oos_fx, n_rep, n0, THR_OOS)
    m_oos_bg <- metrics(oos_bg, n_rep, n0, THR_OOS)
    m_ios_fx <- metrics(ios_fx, n_rep, n0, THR_IOS)
    m_ios_bg <- metrics(ios_bg, n_rep, n0, THR_IOS)
    identical_oos <- as.integer(isTRUE(all.equal(oos_fx, oos_bg, tolerance = 1e-9)))
    identical_ios <- as.integer(isTRUE(all.equal(ios_fx, ios_bg, tolerance = 1e-9)))
    data.frame(
      rep = rep, n_cls = g$n_cls, madn0_oos = madn0_oos,
      ios_ncl_mad0 = ncl_mad0, ios_npts_mad0 = npts_mad0,
      oos_identical = identical_oos, ios_identical = identical_ios,
      AUC_OOS_fix = m_oos_fx["AUC"], AUC_OOS_bug = m_oos_bg["AUC"],
      AUC_IOS_fix = m_ios_fx["AUC"], AUC_IOS_bug = m_ios_bg["AUC"],
      TPR_OOS_fix = m_oos_fx["TPR"], TNR_OOS_fix = m_oos_fx["TNR"],
      BA_OOS_fix  = m_oos_fx["BA"],  F2_OOS_fix  = m_oos_fx["F2"],
      TPR_OOS_bug = m_oos_bg["TPR"], TNR_OOS_bug = m_oos_bg["TNR"],
      TPR_IOS_fix = m_ios_fx["TPR"], TNR_IOS_fix = m_ios_fx["TNR"],
      BA_IOS_fix  = m_ios_fx["BA"],  F2_IOS_fix  = m_ios_fx["F2"],
      TPR_IOS_bug = m_ios_bg["TPR"], TNR_IOS_bug = m_ios_bg["TNR"],
      stringsAsFactors = FALSE)
  }, error = function(e) data.frame(rep = rep, err = conditionMessage(e)))
  out
}
stopCluster(cl)

if (!"AUC_OOS_fix" %in% names(per)) {
  print(utils::head(per))
  stop("worker rows malformed -- every rep errored. See the per-rep frame above.")
}

ok <- per[!is.na(per$AUC_OOS_fix), , drop = FALSE]
nrep_ok <- nrow(ok); nrep_err <- REPS - nrep_ok

## ---- aggregate + report ---------------------------------------------------
agg <- data.frame(
  cell = CELL, family = FAMILY, dist = DIST, d = D, n = N,
  reps = REPS, reps_ok = nrep_ok, reps_err = nrep_err,
  mean_n_cls          = mean(ok$n_cls),
  reps_madn0_oos      = sum(ok$madn0_oos),
  reps_ios_any_mad0   = sum(ok$ios_ncl_mad0 > 0),
  mean_ios_ncl_mad0   = mean(ok$ios_ncl_mad0),
  reps_oos_identical  = sum(ok$oos_identical),
  reps_ios_identical  = sum(ok$ios_identical),
  AUC_OOS_fix = mean(ok$AUC_OOS_fix), AUC_OOS_bug = mean(ok$AUC_OOS_bug),
  AUC_IOS_fix = mean(ok$AUC_IOS_fix), AUC_IOS_bug = mean(ok$AUC_IOS_bug),
  TPR_OOS_fix = mean(ok$TPR_OOS_fix), TNR_OOS_fix = mean(ok$TNR_OOS_fix),
  BA_OOS_fix  = mean(ok$BA_OOS_fix),  F2_OOS_fix  = mean(ok$F2_OOS_fix),
  TPR_OOS_bug = mean(ok$TPR_OOS_bug), TNR_OOS_bug = mean(ok$TNR_OOS_bug),
  TPR_IOS_fix = mean(ok$TPR_IOS_fix), TNR_IOS_fix = mean(ok$TNR_IOS_fix),
  BA_IOS_fix  = mean(ok$BA_IOS_fix),  F2_IOS_fix  = mean(ok$F2_IOS_fix),
  TPR_IOS_bug = mean(ok$TPR_IOS_bug), TNR_IOS_bug = mean(ok$TNR_IOS_bug),
  stringsAsFactors = FALSE)

res_dir <- here::here("revision_experiments/results/tr2")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
perrep_file <- file.path(res_dir, sprintf("highd_madn_recheck_%s_perrep.csv", CELL))
write.csv(ok, perrep_file, row.names = FALSE)
agg_file <- file.path(res_dir, "highd_madn_recheck.csv")
write.table(agg, agg_file, sep = ",", row.names = FALSE,
            col.names = !file.exists(agg_file), append = file.exists(agg_file))

message("\n========== ", CELL, " ==========")
message(sprintf("reps ok/err            : %d / %d", nrep_ok, nrep_err))
message(sprintf("mean #clusters         : %.2f", agg$mean_n_cls))
message(sprintf("reps with MADN0(OOS)   : %d / %d", agg$reps_madn0_oos, nrep_ok))
message(sprintf("reps with any IOS mad0 : %d / %d (mean clusters affected %.2f)",
                agg$reps_ios_any_mad0, nrep_ok, agg$mean_ios_ncl_mad0))
message(sprintf("reps OOS fixed==buggy  : %d / %d", agg$reps_oos_identical, nrep_ok))
message(sprintf("reps IOS fixed==buggy  : %d / %d", agg$reps_ios_identical, nrep_ok))
message(sprintf("OOS  AUC fix=%.4f bug=%.4f | TPR fix=%.4f bug=%.4f (TNR %.4f)",
                agg$AUC_OOS_fix, agg$AUC_OOS_bug, agg$TPR_OOS_fix, agg$TPR_OOS_bug, agg$TNR_OOS_fix))
message(sprintf("IOS  AUC fix=%.4f bug=%.4f | TPR fix=%.4f bug=%.4f (TNR %.4f)",
                agg$AUC_IOS_fix, agg$AUC_IOS_bug, agg$TPR_IOS_fix, agg$TPR_IOS_bug, agg$TNR_IOS_fix))
message("wrote: ", agg_file)
message("wrote: ", perrep_file)
message("[recheck] DONE ", CELL, "  ", format(Sys.time()))
