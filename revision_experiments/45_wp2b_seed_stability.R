#!/usr/bin/env Rscript
# 45_wp2b_seed_stability.R -- WP2(b): stability of the four MCCD detectors
# under (i) the R random-number stream and (ii) small perturbations of the
# input data.  Answers reviewer point R3.1.
#
#   Rscript 45_wp2b_seed_stability.R --determinism            # Step 1
#   Rscript 45_wp2b_seed_stability.R [datasets] [methods] [ptypes] [repspec]
#   Rscript 45_wp2b_seed_stability.R --summary                # Steps 2-4 + anchor
#
# ---------------------------------------------------------------------------
# STEP 1 -- WHERE DOES RANDOMNESS ENTER?  (static reading, verified by --determinism)
# ---------------------------------------------------------------------------
# Static trace of every code path reachable from the four detectors:
#
#   U-MCCD   RUMCCD_outlier  -> RKCCD_correct_quant -> rccd.clustering_correct_quantile
#                              -> ccd.Kest.edge.quantile (harness override)
#                              -> dominate.mat.greedy2, dominate.mat.ks
#                              -> rccd.silhouette -> cluster::silhouette
#                              -> connected.ksccd.m -> ksccd.connected
#   SU-MCCD  SUMCCD_outlier   same, plus min.cls in rccd.silhouette
#   UN-MCCD  UNMCCD_outlier  -> nnccd_clustering_quantile -> nnccd.radi (harness override)
#                              -> dominate.mat.greedy2, dominate.mat.ks
#                              -> nnccd.silhouette -> cluster::silhouette
#                              -> connected.ksccd.m -> ksccd.connected
#   SUN-MCCD SUNMCCD_outlier  same, plus min.cls in nnccd.silhouette
#
# `grep -nE 'sample\(|runif|rnorm|set\.seed|kmeans|sample\.int|jitter'` over
# methods/outlier_detection/ returns NOTHING, and over R/ccds/ returns hits only in
#   - Kest.R / NN_Dist_Est.R  : rpoisball.unit / rpoisbox.unit / mvrnorm, i.e. the
#     Monte-Carlo NULL-DISTRIBUTION generators.  These are reached only through
#     Kest.simpois.edge.quantile() / NNDest.simpois.lower.quant(), which
#     ccd.Kest.edge.quantile() and nnccd.radi() call ONLY when `simul` is NULL.
#     Every wrapper in wp0_mccd_methods.R passes a non-NULL `simul` loaded from
#     R/{RK,NN}-test_quantile/*.RData by get_simul(), so that branch is dead at
#     inference time.  The Monte Carlo draws are frozen in those .RData files.
#   - K-S_CCD_Rui_trash.R and the R/*-test_quantile/*.R table generators, none of
#     which is sourced by the detector chain.
#
# The remaining candidates named in the brief are all deterministic:
#   * There is no k-means anywhere; macro-clusters come from a greedy
#     approximate dominating set, not from a randomly initialised partition.
#   * dominate.mat.greedy2 (ccdfunctions.R:63) and dominate.mat.ks (:108) break
#     ties with which.max(), i.e. lowest index wins -- deterministic, not random.
#   * The silhouette-based choice of the number of clusters (rccd.silhouette
#     RK_CCD_New.R:108 / nnccd.silhouette UN_CCD.R:99) is a plain forward scan
#     over i = 2..lenD keeping a strict `datasi > maxsi` maximum -- again ties
#     resolve to the FIRST index, deterministically.
#   * connected.ksccd.m (mKNN_CCD_functions.R:50) is a bracketing search on a
#     scalar; igraph::components and cluster::silhouette are deterministic.
#   * dom.method = "greedy2" is an APPROXIMATE minimum-dominating-set heuristic,
#     but "approximate" here means "not provably minimum", not "randomised".
#
# Predicted verdict: given fixed data and a fixed quantile table the four
# detectors are exactly deterministic and never touch the RNG stream.  Two
# independent empirical checks are run by --determinism and neither is allowed
# to substitute for the other:
#   (a) OUTCOME check -- same data, different set.seed(), compare the score
#       vector, per-point radii, macro-cluster vector and delta vector with
#       identical().  Bit-for-bit, not to 3 dp.
#   (b) STREAM check -- capture .Random.seed immediately before and after the
#       detector call and compare with identical(); and, independently, compare
#       runif(1) drawn from seed 777 with NO detector call in between against
#       runif(1) drawn from seed 777 WITH the detector call in between.  If the
#       detector consumed even one variate the two draws must differ.  This
#       separates "the RNG was never consumed" from "the RNG was consumed but
#       happened not to change the answer on these data", which are different
#       claims.  (Note the two draws must come from the SAME seed; comparing
#       draws made after two different set.seed() values proves nothing.)
#
# ---------------------------------------------------------------------------
# STEP 3 -- PERTURBATION STUDY (runs regardless of the Step 1 verdict)
# ---------------------------------------------------------------------------
# R3.1's substantive point is that graph connectivity is a discontinuous
# function of the data, so an arbitrarily small perturbation can split a
# component and flip a label.  That is testable without any RNG inside the
# method: perturb the INPUT.
#
#   ptype = "none"    the unperturbed reference cell (rep 0). Also the Step-4
#                     anchor against results/tr1/final_comparison.csv.
#   ptype = "jitter"  X[, j] <- X[, j] + rnorm(n, 0, SIGMA * sd(X[, j]))
#                     with SIGMA = 0.01, i.i.d. per cell, per feature.
#   ptype = "dropout" leave-p-out bootstrap: drop round(0.01 * n) rows at random.
#
# Perturbation seed is set.seed(1000 + rep) and is recorded in the `seed`
# column of every row, so any cell can be reproduced in isolation.
#
# ---------------------------------------------------------------------------
# TAPS
# ---------------------------------------------------------------------------
# The detectors return radii = list(R) where R = r[D] -- the per-point radius
# vector PERMUTED by the dominating-set order D, and D is not returned (see
# SU-MCCDs.R:14, which itself has to undo the permutation with order(D)).
# Point identity is therefore unrecoverable from the returned object, and the
# radius comparison R3.1 asks for needs point identity.  We install two
# read-only taps -- wrappers that call the untouched original and record what
# it produced:
#   RADTAP  around ccd.Kest.edge.quantile / nnccd.radi -> per-point radii in
#           ORIGINAL row order (both are called with dx = datax).
#   DTAP    around connected.ksccd.m -> one delta per macro-cluster, so
#           length(DTAP$deltas) == si.ind == the number of macro-clusters
#           (same tap idiom as 40_wp2a_tolerance_sweep.R).
# Both are installed AFTER every source() has run, because
# methods/outlier_detection/*.R re-source their dependencies and would clobber
# an earlier binding.  connected.ksccd.m internally re-sources ccd_ks_NEW.R and
# ccdfunctions.R; verified that neither of those files defines
# ccd.Kest.edge.quantile, nnccd.radi or connected.ksccd.m, so no tap is lost
# mid-run.  A guard in build_summary() asserts the taps actually fired.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- taps: install only now, after every source() has run --------------------

TAP <- new.env(parent = emptyenv())
TAP$radii   <- NULL          # per-point radii, original row order, first call only
TAP$radii_n <- 0L            # how many times a radius function fired this cell
TAP$deltas  <- numeric(0)    # one delta per macro-cluster

tap_reset <- function() { TAP$radii <- NULL; TAP$radii_n <- 0L; TAP$deltas <- numeric(0) }

local({
  orig_kest  <- get("ccd.Kest.edge.quantile", envir = globalenv())
  orig_nnrad <- get("nnccd.radi",             envir = globalenv())
  orig_cksm  <- get("connected.ksccd.m",      envir = globalenv())

  assign("ccd.Kest.edge.quantile", function(...) {
    v <- orig_kest(...)
    TAP$radii_n <- TAP$radii_n + 1L
    if (is.null(TAP$radii)) TAP$radii <- v$R
    v
  }, envir = globalenv())

  assign("nnccd.radi", function(...) {
    v <- orig_nnrad(...)
    TAP$radii_n <- TAP$radii_n + 1L
    if (is.null(TAP$radii)) TAP$radii <- v$R
    v
  }, envir = globalenv())

  assign("connected.ksccd.m", function(t) {
    v <- orig_cksm(t)
    TAP$deltas <- c(TAP$deltas, v)
    v
  }, envir = globalenv())
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

S_MIN           <- 0.0625        # proportion of n (manuscript value)
SIGMA           <- 0.01          # jitter scale, in per-feature SD units
DROP_FRAC       <- 0.01          # leave-p-out fraction
N_REPS          <- 50
PERT_SEED_BASE  <- 1000
CELL_TIMEOUT    <- 420           # seconds
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
ALL_METHODS     <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")

OUT_CSV     <- here::here("revision_experiments/results/tr1/wp2b_seed_stability.csv")
SUMMARY_CSV <- here::here("revision_experiments/results/tr1/wp2b_stability_summary.csv")
DET_CSV     <- here::here("revision_experiments/results/tr1/wp2b_determinism.csv")
ANCHOR_CSV  <- here::here("revision_experiments/results/tr1/final_comparison.csv")

args <- commandArgs(trailingOnly = TRUE)
MODE_DET <- any(args == "--determinism")
MODE_SUM <- any(args == "--summary")
args <- args[!(args %in% c("--determinism", "--summary"))]

DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else
              c("hepatitis", "glass", "stamps", "pima")
METHODS  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else ALL_METHODS
PTYPES   <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else
              c("none", "jitter", "dropout")
REPS     <- if (length(args) >= 4 && nzchar(args[4])) {
              p <- as.integer(strsplit(args[4], ":")[[1]]); if (length(p) == 1) p else p[1]:p[2]
            } else 1:N_REPS

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

collapse_num <- function(v) if (length(v) == 0) "" else paste(signif(v, 10), collapse = ";")
collapse_int <- function(v) if (length(v) == 0) "" else paste(as.integer(v), collapse = ";")
parse_num <- function(s) if (is.na(s) || !nzchar(s)) numeric(0) else as.numeric(strsplit(s, ";")[[1]])
parse_int <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else as.integer(strsplit(s, ";")[[1]])

#' Adjusted Rand index of two label vectors of equal length (NA allowed; pairs
#' with an NA on either side are dropped).  Hand-rolled so the script has no
#' dependency beyond what harness.R already attaches.
adj_rand <- function(a, b) {
  keep <- !is.na(a) & !is.na(b)
  a <- a[keep]; b <- b[keep]
  n <- length(a)
  if (n < 2) return(NA_real_)
  if (length(unique(a)) == 1 && length(unique(b)) == 1) return(1)
  tab <- table(a, b)
  cc  <- function(x) sum(x * (x - 1) / 2)
  sij <- cc(as.vector(tab)); sa <- cc(rowSums(tab)); sb <- cc(colSums(tab))
  tot <- n * (n - 1) / 2
  expected <- sa * sb / tot
  maxi     <- (sa + sb) / 2
  if (isTRUE(all.equal(maxi, expected))) return(1)
  (sij - expected) / (maxi - expected)
}

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(1)
  length(intersect(a, b)) / u
}

#' Run one detector on one (X, Y) and return everything WP2(b) needs.
run_cell <- function(meth, X, Y, d) {
  tap_reset()
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
    res <- if (meth %in% MIN_CLS_METHODS) {
      METHOD_REGISTRY[[meth]](X = X, d = d, Y = Y, min.cls = S_MIN)
    } else {
      METHOD_REGISTRY[[meth]](X = X, d = d, Y = Y)
    }
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    list(res = res, status = "ok", note = NA_character_)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    list(res = NULL,
         status = if (grepl("elapsed time limit", msg)) "timeout" else "error",
         note = msg)
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(out$res)) {
    return(list(status = out$status, note = out$note, elapsed = wall,
                score = NULL, metrics = setNames(rep(NA_real_, 4), c("TPR","TNR","BA","F2")),
                radii = numeric(0), deltas = numeric(0), cluster = integer(0),
                unassigned = NA_integer_))
  }
  m <- evaluate(Y, out$res$score, REAL_DATA_THRESHOLDS[[meth]])
  list(status = "ok", note = NA_character_, elapsed = wall,
       score = out$res$score, metrics = m,
       radii = if (is.null(TAP$radii)) numeric(0) else TAP$radii,
       deltas = TAP$deltas, cluster = out$res$cluster,
       unassigned = out$res$unassigned_rows)
}

# ---------------------------------------------------------------------------
# STEP 1 -- determinism probe
# ---------------------------------------------------------------------------

run_determinism <- function() {
  ds_list <- if (length(args) >= 1 && nzchar(args[1])) DATASETS else
               c("hepatitis", "glass", "WDBC", "pima")   # 74/19, 213/9, 367/30, 555/8
  seeds <- c(1L, 20260809L)
  rows <- list()
  cat("45_wp2b STEP 1: determinism probe\n")
  cat(sprintf("  datasets = %s\n  methods  = %s\n  seeds    = %s\n",
              paste(ds_list, collapse = ", "), paste(METHODS, collapse = ", "),
              paste(seeds, collapse = ", ")))
  for (ds in ds_list) {
    dat <- load_real_dataset(ds)
    for (meth in METHODS) {
      keep <- list()
      rs_touched <- NA
      for (k in seq_along(seeds)) {
        set.seed(seeds[k])
        rs_before <- get(".Random.seed", envir = globalenv())
        r <- run_cell(meth, dat$X, dat$Y, dat$d)
        rs_after <- get(".Random.seed", envir = globalenv())
        if (k == 1) rs_touched <- !identical(rs_before, rs_after)
        keep[[k]] <- r
      }
      # STREAM check: both draws start from the SAME seed; they can differ only
      # if the detector call between set.seed() and runif() consumed variates.
      set.seed(777L); u_noskip <- runif(1)
      set.seed(777L); invisible(run_cell(meth, dat$X, dat$Y, dat$d)); u_after <- runif(1)
      runif_after <- c(u_noskip, u_after)
      a <- keep[[1]]; b <- keep[[2]]
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, method = meth, n = dat$n, d = dat$d,
        seed_a = seeds[1], seed_b = seeds[2],
        status_a = a$status, status_b = b$status,
        score_identical   = identical(a$score,   b$score),
        radii_identical   = identical(a$radii,   b$radii),
        cluster_identical = identical(a$cluster, b$cluster),
        delta_identical   = identical(a$deltas,  b$deltas),
        metrics_identical = identical(a$metrics, b$metrics),
        random_seed_touched = rs_touched,
        runif_seed777_no_call   = runif_after[1],
        runif_seed777_after_call = runif_after[2],
        runif_identical = identical(runif_after[1], runif_after[2]),
        n_flagged_a = if (is.null(a$score)) NA_integer_ else sum(a$score > 0.5),
        n_flagged_b = if (is.null(b$score)) NA_integer_ else sum(b$score > 0.5),
        n_clusters_a = length(a$deltas), n_clusters_b = length(b$deltas),
        BA_a = unname(a$metrics[["BA"]]), BA_b = unname(b$metrics[["BA"]]),
        radii_tap_calls = TAP$radii_n,
        elapsed_a = a$elapsed, elapsed_b = b$elapsed,
        stringsAsFactors = FALSE)
      rr <- rows[[length(rows)]]
      cat(sprintf("  %-10s %-9s score=%-5s radii=%-5s cluster=%-5s delta=%-5s | RNG touched=%-5s runif same=%-5s | %.1fs\n",
                  ds, meth, rr$score_identical, rr$radii_identical, rr$cluster_identical,
                  rr$delta_identical, rr$random_seed_touched, rr$runif_identical, rr$elapsed_a))
      flush.console()
    }
  }
  out <- do.call(rbind, rows)
  dir.create(dirname(DET_CSV), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(DET_CSV)) {
    old <- read.csv(DET_CSV, stringsAsFactors = FALSE)
    key_old <- paste(old$dataset, old$method); key_new <- paste(out$dataset, out$method)
    out <- rbind(old[!(key_old %in% key_new), , drop = FALSE], out)
  }
  out <- out[order(out$dataset, out$method), ]
  write.csv(out, DET_CSV, row.names = FALSE)
  cat(sprintf("\n45_wp2b: wrote %s (%d rows)\n", DET_CSV, nrow(out)))
  ok <- out$status_a == "ok" & out$status_b == "ok"
  cat(sprintf("  outcome identical across seeds : %d / %d cells\n",
              sum(ok & out$score_identical & out$radii_identical &
                  out$cluster_identical & out$delta_identical), sum(ok)))
  cat(sprintf("  .Random.seed UNCHANGED by call : %d / %d cells\n",
              sum(ok & !out$random_seed_touched), sum(ok)))
  cat(sprintf("  runif(1) after call identical  : %d / %d cells\n",
              sum(ok & out$runif_identical), sum(ok)))
  invisible(out)
}

# ---------------------------------------------------------------------------
# STEP 3 -- perturbation grid
# ---------------------------------------------------------------------------

make_perturbed <- function(ptype, X, Y, rep) {
  n <- nrow(X)
  if (ptype == "none") {
    return(list(X = X, Y = Y, retained = seq_len(n), seed = NA_integer_))
  }
  seed <- PERT_SEED_BASE + rep
  set.seed(seed)
  if (ptype == "jitter") {
    Xm <- as.matrix(X)
    sds <- apply(Xm, 2, stats::sd)
    sds[!is.finite(sds)] <- 0
    for (j in seq_len(ncol(Xm))) {
      Xm[, j] <- Xm[, j] + stats::rnorm(n, 0, SIGMA * sds[j])
    }
    return(list(X = as.data.frame(Xm), Y = Y, retained = seq_len(n), seed = seed))
  }
  if (ptype == "dropout") {
    n_drop <- max(1L, round(DROP_FRAC * n))
    drop <- sort(sample.int(n, n_drop))
    keep <- setdiff(seq_len(n), drop)
    return(list(X = X[keep, , drop = FALSE], Y = Y[keep], retained = keep, seed = seed))
  }
  stop("unknown ptype: ", ptype)
}

load_done <- function() {
  if (!file.exists(OUT_CSV)) return(character(0))
  df <- tryCatch(read.csv(OUT_CSV, stringsAsFactors = FALSE, colClasses = "character"),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(character(0))
  paste(df$dataset, df$method, df$ptype, df$rep, sep = "|")
}

run_grid <- function() {
  cat(sprintf("45_wp2b STEP 3: datasets=%s  methods=%s  ptypes=%s  reps=%s\n",
              paste(DATASETS, collapse = ","), paste(METHODS, collapse = ","),
              paste(PTYPES, collapse = ","),
              if (length(REPS) > 1) sprintf("%d:%d", min(REPS), max(REPS)) else as.character(REPS)))
  cat(sprintf("  sigma=%g  drop_frac=%g  S_min=%g  out=%s\n", SIGMA, DROP_FRAC, S_MIN, OUT_CSV))
  done <- load_done()
  for (ds in DATASETS) {
    dat <- load_real_dataset(ds)
    for (meth in METHODS) {
      for (pt in PTYPES) {
        reps <- if (pt == "none") 0L else REPS
        for (rp in reps) {
          key <- paste(ds, meth, pt, rp, sep = "|")
          if (key %in% done) { cat(sprintf("[skip] %s\n", key)); next }
          pr <- make_perturbed(pt, dat$X, dat$Y, rp)
          r  <- run_cell(meth, pr$X, pr$Y, dat$d)

          flagged_local <- if (is.null(r$score)) integer(0) else which(r$score > 0.5)
          flagged_orig  <- pr$retained[flagged_local]
          radii_ok      <- length(r$radii) == length(pr$retained)

          append_result(OUT_CSV, list(
            dataset = ds, method = meth, ptype = pt, rep = rp, seed = pr$seed,
            n_full = dat$n, n_eff = nrow(pr$X), d = dat$d,
            n_outliers_eff = sum(pr$Y == 0),
            TPR = unname(r$metrics[["TPR"]]), TNR = unname(r$metrics[["TNR"]]),
            BA  = unname(r$metrics[["BA"]]),  F2  = unname(r$metrics[["F2"]]),
            n_flagged = if (is.null(r$score)) NA_integer_ else length(flagged_local),
            n_clusters = length(r$deltas),
            unassigned_rows = r$unassigned,
            radii_len = length(r$radii),
            radii_aligned = radii_ok,
            radii_mean = if (length(r$radii)) mean(r$radii) else NA_real_,
            radii_sd   = if (length(r$radii) > 1) stats::sd(r$radii) else NA_real_,
            delta_values  = collapse_num(r$deltas),
            retained_idx  = collapse_int(pr$retained),
            flagged_idx   = collapse_int(flagged_orig),
            cluster_vec   = collapse_int(r$cluster),
            radii_vec     = collapse_num(r$radii),
            elapsed_sec = r$elapsed, s_min = if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_,
            status = r$status, note = r$note, timestamp = format(Sys.time())
          ))
          cat(sprintf("  %-10s %-9s %-8s rep=%-3d n=%-5d BA=%s F2=%s nflag=%-4s ncls=%-2d %s %.1fs\n",
                      ds, meth, pt, rp, nrow(pr$X),
                      if (is.na(r$metrics[["BA"]])) "  NA " else sprintf("%.3f", r$metrics[["BA"]]),
                      if (is.na(r$metrics[["F2"]])) "  NA " else sprintf("%.3f", r$metrics[["F2"]]),
                      if (is.null(r$score)) "NA" else length(flagged_local),
                      length(r$deltas), r$status, r$elapsed))
          flush.console()
        }
      }
    }
  }
  cat("45_wp2b: grid chunk done\n")
}

# ---------------------------------------------------------------------------
# STEPS 2/3/4 -- summary, per-point flip rates, anchor check
# ---------------------------------------------------------------------------

build_summary <- function() {
  chr <- c(delta_values = "character", retained_idx = "character",
           flagged_idx = "character", cluster_vec = "character",
           radii_vec = "character")
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE, colClasses = chr)
  for (k in names(chr)) df[[k]][is.na(df[[k]])] <- ""

  anchor <- read.csv(ANCHOR_CSV, stringsAsFactors = FALSE)

  rows <- list()
  for (ds in unique(df$dataset)) {
    for (meth in unique(df$method[df$dataset == ds])) {
      sub  <- df[df$dataset == ds & df$method == meth, ]
      base <- sub[sub$ptype == "none", ]
      if (nrow(base) == 0) next
      base <- base[nrow(base), ]
      base_flag  <- parse_int(base$flagged_idx)
      base_clust <- parse_int(base$cluster_vec)
      base_rad   <- parse_num(base$radii_vec)
      n_full     <- base$n_full

      # ---- Step 4 anchor -----------------------------------------------------
      arow <- anchor[anchor$dataset == ds & anchor$method == meth, ]
      anch_pass <- NA; anch_BA <- NA_real_; anch_F2 <- NA_real_
      anch_TPR <- NA_real_; anch_TNR <- NA_real_
      if (nrow(arow) == 1) {
        anch_BA <- arow$BA; anch_F2 <- arow$F2; anch_TPR <- arow$TPR; anch_TNR <- arow$TNR
        anch_pass <- all(abs(round(c(base$BA, base$F2, base$TPR, base$TNR), 3) -
                             c(anch_BA, anch_F2, anch_TPR, anch_TNR)) < 5e-4)
      }

      for (pt in setdiff(unique(sub$ptype), "none")) {
        s <- sub[sub$ptype == pt & sub$status == "ok", ]
        if (nrow(s) == 0) next

        flags <- lapply(s$flagged_idx,  parse_int)
        clus  <- lapply(s$cluster_vec,  parse_int)
        rads  <- lapply(s$radii_vec,    parse_num)
        rets  <- lapply(s$retained_idx, parse_int)

        # modal flagged set / modal partition (exact string mode)
        fmode_key <- names(sort(table(s$flagged_idx), decreasing = TRUE))[1]
        modal_flag <- parse_int(fmode_key)
        modal_flag_count <- sum(s$flagged_idx == fmode_key)
        cmode_key <- names(sort(table(s$cluster_vec), decreasing = TRUE))[1]
        modal_clust <- parse_int(cmode_key)
        modal_clust_count <- sum(s$cluster_vec == cmode_key)

        # ARI / Jaccard vs. the unperturbed baseline, restricted to retained rows
        ari_b <- numeric(nrow(s)); jac_b <- numeric(nrow(s))
        ari_m <- numeric(nrow(s)); jac_m <- numeric(nrow(s))
        for (i in seq_len(nrow(s))) {
          keep <- rets[[i]]
          cv <- rep(NA_integer_, n_full); cv[keep] <- clus[[i]]
          ari_b[i] <- adj_rand(cv, base_clust)
          jac_b[i] <- jaccard(flags[[i]], intersect(base_flag, keep))
          mv <- rep(NA_integer_, n_full)
          mk <- rets[[which(s$cluster_vec == cmode_key)[1]]]
          mv[mk] <- modal_clust
          ari_m[i] <- adj_rand(cv, mv)
          jac_m[i] <- jaccard(flags[[i]], modal_flag)
        }

        # per-point flip rate vs. baseline label, conditional on being retained
        present <- integer(n_full); flipped <- integer(n_full)
        base_is_flag <- logical(n_full); base_is_flag[base_flag] <- TRUE
        for (i in seq_len(nrow(s))) {
          keep <- rets[[i]]
          present[keep] <- present[keep] + 1L
          cur <- logical(n_full); cur[flags[[i]]] <- TRUE
          flipped[keep] <- flipped[keep] + as.integer(cur[keep] != base_is_flag[keep])
        }
        rate <- ifelse(present > 0, flipped / present, NA_real_)

        # per-point radius variability across reps (aligned by original row id)
        radmat <- matrix(NA_real_, nrow = nrow(s), ncol = n_full)
        for (i in seq_len(nrow(s))) {
          if (length(rads[[i]]) == length(rets[[i]])) radmat[i, rets[[i]]] <- rads[[i]]
        }
        pt_mean <- apply(radmat, 2, mean, na.rm = TRUE)
        pt_sd   <- apply(radmat, 2, stats::sd,  na.rm = TRUE)
        pt_cv   <- ifelse(pt_mean > 0, pt_sd / pt_mean, NA_real_)

        nc_tab <- table(s$n_clusters)
        rows[[length(rows) + 1]] <- data.frame(
          dataset = ds, method = meth, ptype = pt, n_reps = nrow(s),
          n_full = n_full, d = base$d, sigma = if (pt == "jitter") SIGMA else NA_real_,
          drop_frac = if (pt == "dropout") DROP_FRAC else NA_real_,

          radius_mean          = mean(pt_mean, na.rm = TRUE),
          radius_mean_pointsd  = mean(pt_sd,   na.rm = TRUE),
          radius_mean_pointcv  = mean(pt_cv,   na.rm = TRUE),
          radius_max_pointcv   = suppressWarnings(max(pt_cv, na.rm = TRUE)),
          radius_frac_points_changed = mean(pt_sd > 0, na.rm = TRUE),

          n_clusters_base = base$n_clusters,
          n_clusters_mode = as.integer(names(nc_tab)[which.max(nc_tab)]),
          n_clusters_min  = min(s$n_clusters), n_clusters_max = max(s$n_clusters),
          n_clusters_mean = mean(s$n_clusters), n_clusters_sd = stats::sd(s$n_clusters),
          n_clusters_table = paste(sprintf("%s:%d", names(nc_tab), as.integer(nc_tab)),
                                   collapse = ";"),

          ari_vs_base_mean = mean(ari_b), ari_vs_base_min = min(ari_b),
          ari_vs_base_sd = stats::sd(ari_b),
          ari_vs_modal_mean = mean(ari_m), ari_vs_modal_min = min(ari_m),
          n_runs_partition_eq_modal = modal_clust_count,

          jaccard_vs_base_mean = mean(jac_b), jaccard_vs_base_min = min(jac_b),
          jaccard_vs_base_sd = stats::sd(jac_b),
          jaccard_vs_modal_mean = mean(jac_m), jaccard_vs_modal_min = min(jac_m),
          n_runs_flagset_eq_modal = modal_flag_count,
          modal_eq_base = identical(sort(modal_flag), sort(base_flag)),

          flip_rate_mean = mean(rate, na.rm = TRUE),
          flip_rate_max  = suppressWarnings(max(rate, na.rm = TRUE)),
          n_points_ever_flipped = sum(rate > 0, na.rm = TRUE),
          frac_points_ever_flipped = mean(rate > 0, na.rm = TRUE),

          n_flagged_base = base$n_flagged,
          n_flagged_mean = mean(s$n_flagged), n_flagged_sd = stats::sd(s$n_flagged),
          n_flagged_min = min(s$n_flagged), n_flagged_max = max(s$n_flagged),

          BA_base = base$BA, BA_mean = mean(s$BA), BA_sd = stats::sd(s$BA),
          BA_min = min(s$BA), BA_max = max(s$BA), BA_range = max(s$BA) - min(s$BA),
          F2_base = base$F2, F2_mean = mean(s$F2), F2_sd = stats::sd(s$F2),
          F2_min = min(s$F2), F2_max = max(s$F2), F2_range = max(s$F2) - min(s$F2),
          TPR_base = base$TPR, TPR_sd = stats::sd(s$TPR),
          TNR_base = base$TNR, TNR_sd = stats::sd(s$TNR),

          anchor_BA = anch_BA, anchor_F2 = anch_F2,
          anchor_TPR = anch_TPR, anchor_TNR = anch_TNR, anchor_pass = anch_pass,

          n_status_not_ok = sum(sub$ptype == pt & sub$status != "ok"),
          radii_tap_ok = all(sub$radii_aligned[sub$ptype == pt & sub$status == "ok"]),
          elapsed_mean = mean(s$elapsed_sec),
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$dataset, out$method, out$ptype), ]
  dir.create(dirname(SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("45_wp2b: wrote %s (%d rows)\n", SUMMARY_CSV, nrow(out)))

  # ---- guards --------------------------------------------------------------
  cat("\n==== TAP GUARD (radii vector must be per-point, original order) ====\n")
  okrows <- df[df$status == "ok", ]
  cat(sprintf("  rows with radii_len == n_eff: %d / %d\n",
              sum(okrows$radii_aligned), nrow(okrows)))
  cat(sprintf("  rows with n_clusters == 0   : %d  (tap on connected.ksccd.m never fired)\n",
              sum(okrows$n_clusters == 0)))

  cat("\n==== STEP 4 ANCHOR: unperturbed cell vs final_comparison.csv (3 dp) ====\n")
  an <- unique(out[, c("dataset", "method", "BA_base", "F2_base", "TPR_base", "TNR_base",
                       "anchor_BA", "anchor_F2", "anchor_TPR", "anchor_TNR", "anchor_pass")])
  for (i in seq_len(nrow(an))) {
    cat(sprintf("  %-10s %-9s BA %.3f/%.3f  F2 %.3f/%.3f  TPR %.3f/%.3f  TNR %.3f/%.3f  -> %s\n",
                an$dataset[i], an$method[i], an$BA_base[i], an$anchor_BA[i],
                an$F2_base[i], an$anchor_F2[i], an$TPR_base[i], an$anchor_TPR[i],
                an$TNR_base[i], an$anchor_TNR[i],
                if (isTRUE(an$anchor_pass[i])) "PASS" else "FAIL"))
  }
  cat(sprintf("  anchor: %d / %d cells PASS\n", sum(an$anchor_pass, na.rm = TRUE), nrow(an)))

  cat("\n==== PERTURBATION STABILITY ====\n")
  for (pt in unique(out$ptype)) {
    o <- out[out$ptype == pt, ]
    cat(sprintf("\n-- %s --\n", pt))
    cat(sprintf("  %-10s %-9s %6s %8s %8s %8s %7s %8s %8s %10s\n",
                "dataset", "method", "ARI", "Jaccard", "flip", "radCV", "ncls", "BA_sd", "F2_sd", "ncls_tab"))
    for (i in seq_len(nrow(o))) {
      cat(sprintf("  %-10s %-9s %6.3f %8.3f %8.4f %8.4f %7d %8.4f %8.4f %10s\n",
                  o$dataset[i], o$method[i], o$ari_vs_base_mean[i], o$jaccard_vs_base_mean[i],
                  o$flip_rate_mean[i], o$radius_mean_pointcv[i], o$n_clusters_base[i],
                  o$BA_sd[i], o$F2_sd[i], o$n_clusters_table[i]))
    }
  }
  invisible(out)
}

# ---------------------------------------------------------------------------

if (MODE_DET) {
  run_determinism()
} else if (MODE_SUM) {
  build_summary()
} else {
  run_grid()
}
