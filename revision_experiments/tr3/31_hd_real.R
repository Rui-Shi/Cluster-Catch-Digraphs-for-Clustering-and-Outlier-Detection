#!/usr/bin/env Rscript
# revision_experiments/tr3/31_hd_real.R
#
# Experiment 1, CCD part: UN-, RK- and KS-CCDs on six real data sets with
# 13 <= d <= 30 (wine, hepatitis, lymphography, wdbc, vehicle, segmentation),
# the counterpart of 11_real.R. The data come from tr3_load_hd() in
# 30_hd_load.R, the settings from settings_realdata_hd.csv (table, level,
# direction and grid by the rule in the header of 30_hd_load.R).
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/31_hd_real.R <method> [datasets]
#     <method>   UN-CCDs | RK-CCDs | KS-CCDs
#     [datasets] comma-separated, from wine, hepatitis, lymphography, wdbc,
#                vehicle, segmentation
#                (default: all six, in that order)
#   Environment variables
#     TR3_CORES        number of workers for KS-CCDs (default: all)
#     TR3_OUT_DIR      output folder (default: revision_experiments/results/tr3)
#     TR3_SMOKE_GRID   smoke tests only: KS-CCDs uses this many evenly spaced
#                      values of its grid instead of all 200. Rows are then
#                      written with variant "hd_smoke", never "hd". Unset for
#                      the real run.
#
# What is copied from 11_real.R
#   - UN-CCDs and RK-CCDs run on the main process, one call per data set, with
#     the table from tr3_load_table();
#   - KS-CCDs: foreach over the grid of m, si of each, then
#     m = grid[which.max(si)] and a final run at that m;
#   - ARI = adjustedRandIndex(label, Y), Sil = mean(silhouette(label,
#     dist(X))[, 3]) (NA for one cluster), k_hat = length(unique(label)).
# The method is loaded with tr3_source_method(method, "realdata").
#
# OUTPUT (OUT_DIR)
#   real_hd.csv            one row per (method, data set): variant, method,
#                          dataset, n, d, ari, sil, k_hat, seconds, param
#                          (KS-CCDs: the chosen m), table
#   real_hd_done.csv       finished rows; reruns skip them
#   real_hd_ks_search.csv  KS-CCDs: si for every grid value

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript revision_experiments/tr3/31_hd_real.R <method> [datasets]")
source(here::here("revision_experiments/tr3/tr3_lib.R"))
source(here::here("revision_experiments/tr3/30_hd_load.R"))

METHOD   <- args[1]
DATASETS <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else TR3_HD_DATASETS
stopifnot(METHOD %in% TR3_METHODS)

OUT_DIR <- sub("[/\\\\]+$", "", Sys.getenv("TR3_OUT_DIR", TR3_RES))
OUT_CSV    <- file.path(OUT_DIR, "real_hd.csv")
DONE_CSV   <- file.path(OUT_DIR, "real_hd_done.csv")
SEARCH_CSV <- file.path(OUT_DIR, "real_hd_ks_search.csv")

SMOKE_GRID <- suppressWarnings(as.integer(Sys.getenv("TR3_SMOKE_GRID", "")))
if (is.na(SMOKE_GRID)) SMOKE_GRID <- 0L
VARIANT <- if (SMOKE_GRID > 0) "hd_smoke" else "hd"
# A smoke grid must not write into the results folder.
if (SMOKE_GRID > 0 && identical(normalizePath(OUT_DIR, mustWork = FALSE), normalizePath(TR3_RES, mustWork = FALSE))) {
  stop("31_hd_real.R: TR3_SMOKE_GRID is set, so TR3_OUT_DIR must be a scratch folder, not results/tr3")
}

# Reader for settings_realdata_hd.csv (tr3_real_setting() reads the other file).
tr3_hd_setting <- function(method, dataset) {
  s <- utils::read.csv(file.path(TR3_DIR, "settings_realdata_hd.csv"),
                       stringsAsFactors = FALSE, colClasses = "character")
  r <- s[s$method == method & s$dataset == dataset, ]
  if (nrow(r) != 1) {
    stop(sprintf("tr3_hd_setting(): %d rows for method=%s dataset=%s", nrow(r), method, dataset))
  }
  as.list(r)
}

settings <- lapply(setNames(DATASETS, DATASETS), function(ds) tr3_hd_setting(METHOD, ds))
tr3_source_method(METHOD, "realdata")
real <- tr3_load_hd(DATASETS)

# The table named in the settings must be the one for the data set's final d,
# and UN-CCDs needs a direction (a blank one runs neither search).
for (ds in DATASETS) {
  s <- settings[[ds]]
  d <- ncol(real[[ds]]) - 1
  if (nzchar(s$table) && !grepl(sprintf("_%dd_", d), s$table, fixed = TRUE)) {
    stop(sprintf("31_hd_real.R: %s: table %s does not match d = %d", ds, s$table, d))
  }
  if (METHOD == "UN-CCDs" && !(s$direction %in% c("ascend", "descend"))) {
    stop(sprintf("31_hd_real.R: %s: UN-CCDs needs direction ascend or descend, the settings say '%s'", ds, s$direction))
  }
}

cores <- tr3_cores()
cat(sprintf("31_hd_real.R: %s on %s, variant=%s, out=%s\n", METHOD, paste(DATASETS, collapse = ", "),
            VARIANT, OUT_DIR))
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
    if (SMOKE_GRID > 0) {
      delta_range <- delta_range[unique(round(seq(1, length(delta_range), length.out = SMOKE_GRID)))]
    }
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
