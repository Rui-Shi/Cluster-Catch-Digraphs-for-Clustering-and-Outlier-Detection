#!/usr/bin/env Rscript
# 33_final_baselines.R -- the five baselines at the FROZEN configuration.
#
#   Rscript 33_final_baselines.R <datasets>
#
# Configuration is fixed by BENCHMARK_EXPANSION_RULE.md and must not be
# adjusted after seeing results:
#
#   LOF      MinPts 11-30, threshold 1.5            (unchanged from published)
#   DBSCAN   MinPts = 4, eps = quantile(4-distance, 1 - n0/n), input X
#            -- the published protocol, with the Real_Data_DBSCAN.R:24 defect
#               corrected: that line passed `df`, the full frame INCLUDING the
#               ground-truth label column, so the label was an input feature.
#               Here the input is X. The oracle eps is retained deliberately,
#               as the documented protocol, and is disclosed as such.
#   MST      cont = 0.02, thresh = 1.2 for EVERY data set, no per-set override
#   ODIN     defaults: k = round(sqrt(n)), threshold = round(n^(1/3))
#   iForest  ntrees = 1000, sample_size = min(256, n), threshold 0.55, seed 1
#
# Writes results/final_baselines.csv, checkpointed on (dataset, method).

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))

args <- commandArgs(trailingOnly = TRUE)
DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else
              c("WBC", "WDBC", "lymphography", "pima", "shuffle",
                "PenDigits", "thyroid", "pageblocks")

OUT <- Sys.getenv("FINAL_BASE_OUT_CSV",
                  here::here("revision_experiments/results/tr1/final_baselines.csv"))

cat(sprintf("33_final_baselines.R: datasets = %s\n  output = %s\n",
            paste(DATASETS, collapse = ", "), OUT))

# --- the five scorers, each returning list(score, threshold, params) --------

run_LOF <- function(X, Y) {
  s <- LOF(as.matrix(X), L_MinPts = 11, U_MinPts = 30)
  list(score = s, thr = 1.5, params = "MinPts=11..30, max over k, threshold=1.5")
}

run_DBSCAN <- function(X, Y) {
  X <- as.matrix(X)
  cont <- sum(Y == 0) / length(Y)
  dM <- as.matrix(dist(X))
  kd <- apply(dM, 1, function(r) sort(r)[5])   # 4-distance, as DBSCAN.R computes it
  eps <- as.numeric(quantile(kd, 1 - cont))
  lab <- dbscan::dbscan(X, eps = eps, minPts = 4)$cluster
  list(score = ifelse(lab == 0, 1, 0), thr = 0.5,
       params = sprintf("minPts=4, eps=%.4f (4-dist at 1-contam=%.4f), input=X (df->X defect fixed)",
                        eps, 1 - cont))
}

run_MST <- function(X, Y) {
  lab <- MST_Outlier(as.matrix(X), cont = 0.02, thresh = 1.2)
  list(score = ifelse(lab == 0, 1, 0), thr = 0.5,
       params = "cont=0.02, thresh=1.2 (uniform across all data sets)")
}

run_ODIN <- function(X, Y) {
  X <- as.matrix(X)
  lab <- ODIN(X)
  list(score = ifelse(lab == 0, 1, 0), thr = 0.5,
       params = sprintf("k=round(sqrt(n))=%d, indeg_thr=round(n^(1/3))=%d (defaults)",
                        round(sqrt(nrow(X))), round(nrow(X)^(1/3))))
}

run_iForest <- function(X, Y) {
  X <- as.matrix(X)
  ss <- min(256, nrow(X))
  m <- isotree::isolation.forest(X, ntrees = 1000, sample_size = ss, seed = 1)
  list(score = as.numeric(predict(m, X)), thr = 0.55,
       params = sprintf("ntrees=1000, sample_size=%d, seed=1, threshold=0.55", ss))
}

SCORERS <- list(LOF = run_LOF, DBSCAN = run_DBSCAN, MST = run_MST,
                ODIN = run_ODIN, iForest = run_iForest)

for (ds in DATASETS) {
  dat <- load_real_dataset(ds)
  cat(sprintf("\n-- %s (n=%d, d=%d, n0=%d)\n", ds, dat$n, dat$d, sum(dat$Y == 0)))

  for (meth in names(SCORERS)) {
    if (isTRUE(has_result(OUT, c(dataset = ds, method = meth)))) {
      cat(sprintf("[skip] %s x %s\n", ds, meth)); next
    }
    t0 <- Sys.time()
    out <- tryCatch({
      r <- SCORERS[[meth]](dat$X, dat$Y)
      list(m = evaluate(dat$Y, r$score, r$thr), params = r$params,
           status = "ok", note = NA_character_)
    }, error = function(e) {
      list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
           params = NA_character_, status = "error", note = conditionMessage(e))
    })
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    append_result(OUT, list(
      dataset = ds, method = meth, params = out$params,
      n = dat$n, d = dat$d, n0 = sum(dat$Y == 0),
      TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
      BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
      t_total = wall, status = out$status, note = out$note,
      timestamp = format(Sys.time())))

    cat(sprintf("  %-8s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | %s | %.1fs\n",
                meth, out$m[["TPR"]], out$m[["TNR"]], out$m[["BA"]], out$m[["F2"]],
                out$status, wall))
  }
}
cat("\n33_final_baselines.R: done\n")
