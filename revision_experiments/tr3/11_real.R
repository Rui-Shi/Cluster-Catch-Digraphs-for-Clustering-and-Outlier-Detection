#!/usr/bin/env Rscript
# revision_experiments/tr3/11_real.R
#
# CCD rows of the real-data table (Real_Data_Result_Clustering in the
# manuscript), run the way
# simulations/clustering/<M>_CCDs/Simulations/Real_Data/Real_Data_<M>CCD.R
# ran them, through the same interface as 10_sim.R.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/11_real.R <method> [datasets] [variant]
#     <method>   UN-CCDs | RK-CCDs | KS-CCDs
#     [datasets] comma-separated, from iris, seeds, user_knowledge, wholesale,
#                breast_cancer, aggregation, asymmetric, R15, D31
#                (default iris,seeds,R15; user_knowledge needs readxl)
#     [variant]  published (the only one so far)
#
# What is reproduced from the script
#   - table, level and search direction per data set from
#     settings_realdata.csv (UN-CCDs: descend for seeds and breast_cancer);
#   - UN-CCDs and RK-CCDs run on the main process, one call per data set;
#   - KS-CCDs: foreach over the 200 grid values of m, si of each, then
#     m = grid[which.max(si)] and a final run at that m;
#   - ARI = adjustedRandIndex(label, Y), Sil = mean(silhouette(label,
#     dist(X))[, 3]) and k_hat = length(unique(label)). These differ from the
#     simulations, which average the method's own si and si.ind.
# The loader is Real_Data_Collection.R with its setwd() fixed in memory
# (tr3_load_real()).
#
# OUTPUT (revision_experiments/results/tr3/)
#   real.csv            one row per (method, data set): variant, method,
#                       dataset, n, d, ari, sil, k_hat, seconds, param (KS-CCDs:
#                       the chosen m), table
#   real_done.csv       finished rows; reruns skip them
#   real_ks_search.csv  KS-CCDs: si for every grid value

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript revision_experiments/tr3/11_real.R <method> [datasets] [variant]")
source(here::here("revision_experiments/tr3/tr3_lib.R"))

METHOD   <- args[1]
DATASETS <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("iris", "seeds", "R15")
VARIANT  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "published"
stopifnot(METHOD %in% TR3_METHODS)
if (VARIANT != "published") stop("11_real.R: variant '", VARIANT, "' has no protocol yet")

OUT_CSV    <- tr3_path("real.csv")
DONE_CSV   <- tr3_path("real_done.csv")
SEARCH_CSV <- tr3_path("real_ks_search.csv")

settings <- lapply(setNames(DATASETS, DATASETS), function(ds) tr3_real_setting(METHOD, ds))
tr3_source_method(METHOD, "realdata")
real <- tr3_load_real(DATASETS)

cores <- tr3_cores()
cat(sprintf("11_real.R: %s on %s, variant=%s\n", METHOD, paste(DATASETS, collapse = ", "), VARIANT))
if (METHOD == "KS-CCDs") {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
}

# The loop is at top level so that foreach exports the method to the workers
# (see 10_sim.R).
for (ds in DATASETS) {
  key <- list(variant = VARIANT, method = METHOD, dataset = ds)
  if (has_result(DONE_CSV, key)) {
    cat(sprintf("  %s: done, skipped\n", ds))
    next
  }
  s <- settings[[ds]]
  D <- real[[ds]]
  n_col <- ncol(D)
  X <- D[, -n_col]
  Y <- D[, n_col]
  t0 <- proc.time()[["elapsed"]]

  delta_best <- NULL
  if (METHOD == "KS-CCDs") {
    delta_range <- eval(parse(text = s$delta_grid))
    Sil_seq <- foreach(delta = delta_range, .packages = c("MASS", "cluster")) %dopar%
      ccd_cluster(METHOD, X, list(m = delta))$si
    Sil_seq <- unlist(Sil_seq)
    tr3_append(SEARCH_CSV, data.frame(key, delta_index = seq_along(delta_range),
                                      delta = fmt17(delta_range), si = fmt17(Sil_seq),
                                      stringsAsFactors = FALSE))
    delta_best <- delta_range[which.max(Sil_seq)]
    par <- list(m = delta_best)
  } else if (METHOD == "UN-CCDs") {
    par <- list(simul = tr3_load_table(s$table), direction = s$direction)
  } else {
    par <- list(simul = tr3_load_table(s$table), quant = as.numeric(s$quant))
  }
  res <- ccd_cluster(METHOD, X, par)
  label <- res$label

  n_clusters <- length(unique(label))
  ari <- adjustedRandIndex(label, Y)
  # silhouette() has no value for a single cluster, where the script would stop
  sil <- if (n_clusters > 1) mean(silhouette(label, dist(X))[, 3]) else NA_real_

  tr3_append(OUT_CSV, data.frame(key, n = nrow(X), d = ncol(X), ari = fmt17(ari), sil = fmt17(sil),
                                 k_hat = n_clusters,
                                 seconds = fmt17(proc.time()[["elapsed"]] - t0),
                                 param = if (is.null(delta_best)) "NA" else fmt17(delta_best),
                                 table = if (nzchar(s$table)) s$table else "NA",
                                 stringsAsFactors = FALSE))
  tr3_append(DONE_CSV, data.frame(key, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                  stringsAsFactors = FALSE))
  if (!is.null(delta_best)) print(delta_best)
  print(paste(ds, "x", paste0(METHOD, ": the ARI is"), ari, "the average silhouette index is", sil,
              "number of clusters detected", n_clusters))
}
