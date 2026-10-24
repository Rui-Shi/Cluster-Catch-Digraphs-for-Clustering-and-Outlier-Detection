#!/usr/bin/env Rscript
# revision_experiments/tr3/13_real_baselines.R
#
# Baseline rows of the real-data table (DBSCAN, MST, K-means++, Spectral,
# Louvain), from the ORIGINAL code: each "##### <name> data #####" block of
# simulations/clustering/Algo_Compare_Clustering/<M>/Simulation/Real_Data/
# <script>.R is evaluated verbatim, in .GlobalEnv as the script ran it, after
# the script's own header lines (library calls and settings such as k (DBSCAN),
# thresh (MST) and sigma (Louvain); the DBSCAN header also defines the function
# knee_eps()). Same idea as 12_real_reference.R, for the five baselines.
# Nothing is written to the files 11_real.R / 12_real_reference.R own.
#
# What differs from running the script itself
#   - the script's `source(.../Real_Data_Collection.R)` line is replaced by
#     library() calls for the loader's packages (R.matlab, foreign, dplyr,
#     readxl, in the loader's order, at the line's position) and by
#     tr3_load_real(): the loader's setwd() names a folder that is gone;
#   - MST: the header's cores/makeCluster/registerDoParallel lines and the
#     closing stopCluster(cl) are dropped. The cluster is never used by the
#     blocks, and the machine is shared (at most 6 cores allowed). MST_Clustering.R
#     is sourced as the script did;
#   - DBSCAN: kNNdistplot()/abline() draw on a null pdf device instead of
#     writing Rplots.pdf;
#   - header lines are evaluated as complete expressions: a line that ends inside
#     an expression ("unexpected end of input") is joined with the next ones until
#     it parses (since 2026-10-23, because the DBSCAN header defines the multi-line
#     function knee_eps()); any other parse error stops the run;
#   - the driver's own variables live in local(), and only the header and the blocks
#     are evaluated in .GlobalEnv, so a loop variable of the driver cannot collide
#     with a variable of the script;
#   - package versions are the ones installed now (igraph 2.3.3, dbscan 1.2.5,
#     ClusterR 1.3.7, kernlab 0.9.33, mclust 6.1.3, cluster 2.1.8.2).
#
# Baseline scripts since 2026-10-23 (README.md, "Baseline fixes (2026-10-23)")
#   - every dataset block of every script starts with set.seed(1), so no cell
#     depends on the block order; K-means++ also passes seed = 1 to KMeans_rcpp();
#   - K-means++ on wholesale uses k = 3 (the script had k = 6);
#   - Spectral: specc(as.matrix(X), centers = k), the Gaussian (rbfdot) kernel on the
#     data with the width chosen by specc (kpar = "automatic"); the block stores it
#     in sigma_used. The script used to pass the distance matrix as.matrix(dist(X)),
#     which specc takes as a data matrix;
#   - DBSCAN: minPts = 5 (k = 4), Eps from the elbow rule of knee_eps(X, k), defined
#     in the script's header: sort the k-NN distances d_(1) <= ... <= d_(n), scale to
#     x_i = (i-1)/(n-1), y_i = (d_(i) - d_(1))/(d_(n) - d_(1)), Eps = d_(i*) * (1 + 1e-9)
#     with i* = argmax (x_i - y_i), the first index on ties. The factor keeps the
#     point that defines d_(i*) inside the radius. Eps used to be typed in;
#   - no set.seed() is added by this driver in the verbatim mode. In the sweep mode
#     the block's own `set.seed(1)` is replaced by set.seed(s); in the kmseed mode
#     `seed = 1)` becomes `seed = s)`. Both stop if the block lacks exactly one such
#     text. The sweep therefore varies only Spectral and Louvain: K-means++ passes
#     its own seed to KMeans_rcpp(), and DBSCAN and MST use no random numbers, so
#     sweep changes nothing for these three (only kmseed varies K-means++).
#
# USAGE (from the clone root, PowerShell; Rscript segfaults under Git Bash)
#   Rscript revision_experiments/tr3/13_real_baselines.R <method> [datasets] [mode]
#   method    DBSCAN | MST | K-means++ | Spectral | Louvain
#   datasets  comma list, default all of
#             iris,seeds,user_knowledge,wholesale,breast_cancer,aggregation,asymmetric,R15,D31
#   mode      "verbatim" (default): each block once, with its own set.seed(1)
#             "sweep:1:30": the block's set.seed(1) replaced by set.seed(s), s = 1..30
#             "kmseed:1:30" (K-means++ only): block text with the call's
#             `seed = 1` replaced by seed = s, s = 1..30 (set.seed(1) stays);
#             rows carry method "K-means++ (seed arg)"
#             "eps:<from>:<to>:<by>" (DBSCAN only, exploratory): eps replaced
#             by each value of the grid; ancillary, written to real_baselines_eps.csv
#
# OUTPUT  revision_experiments/results/tr3/real_baselines.csv        (verbatim)
#         revision_experiments/results/tr3/real_baselines_seeds.csv  (sweep, kmseed)
#         revision_experiments/results/tr3/real_baselines_eps.csv    (eps)
#         method, dataset, n, d, ari, sil, k_hat, seed, script   (eps file: + eps)
#         numbers with fmt17; seed is NA in verbatim mode.
#         Next to each output file (not in eps mode) a "<file stem>_params.csv" gets
#         method, dataset, param, value, seed, the parameter each block chose: DBSCAN
#         `eps`; Spectral `sigma` (specc's inverse kernel width, k = exp(-sigma ||x-y||^2))
#         and `bandwidth` = 1/sqrt(2 sigma) (the s of exp(-||x-y||^2 / (2 s^2))).
#         The log goes to results/tr3/_logs/real_baselines_<method>.log (by the caller).
#
# NOTES ON THE FILES
#   - Rows are appended, never replaced: a second verbatim run of a cell adds a second
#     row, so the row order within (method, dataset) is the run order. The verification
#     (real_verify_baselines.csv) takes the first row as "the rerun" and the second as the
#     run-to-run check.
#   - Parallel processes must not append to one file: their lines interleave on G:. Set
#     the environment variable TR3_BASELINES_OUT=<file name> to make a process write to
#     results/tr3/<file name> (a value that starts with a drive letter and a slash, or
#     with a slash, is used as the path itself), and concatenate the part files after the
#     processes ended. real_baselines_seeds.csv was built this way (part files part_*.csv,
#     since removed).
#   - Do not edit this file while a process runs it: Rscript reads the file as it goes, and
#     an edit makes the finished process end with a parse error (exit status 1) after the
#     last row was written.

source(here::here("revision_experiments/tr3/tr3_lib.R"))

local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("usage: Rscript revision_experiments/tr3/13_real_baselines.R <method> [datasets] [mode]")

  ALL_DS <- c("iris", "seeds", "user_knowledge", "wholesale", "breast_cancer", "aggregation",
              "asymmetric", "R15", "D31")
  BASE <- "simulations/clustering/Algo_Compare_Clustering"
  SCRIPTS <- c("DBSCAN"    = file.path(BASE, "DBSCAN/Simulation/Real_Data/Real_Data_DBSCAN.R"),
               "MST"       = file.path(BASE, "MST/Simulation/Real_Data/Real_Data_MST_fix_cutoff.R"),
               "K-means++" = file.path(BASE, "K-means++/Simulation/Real_Data/Real_Data_kmeanpp.R"),
               "Spectral"  = file.path(BASE, "Spectral_Clustering/Simulation/Real_Data/Real_Data_Spectrum_Clustering.R"),
               "Louvain"   = file.path(BASE, "Louvain/Simulation/Real_Data/Real_Data_Louvain.R"))

  METHOD   <- args[1]
  stopifnot(METHOD %in% names(SCRIPTS))
  DATASETS <- if (length(args) >= 2 && nzchar(args[2]) && args[2] != "all") strsplit(args[2], ",")[[1]] else ALL_DS
  MODE     <- if (length(args) >= 3) args[3] else "verbatim"
  stopifnot(all(DATASETS %in% ALL_DS))
  mp <- strsplit(MODE, ":")[[1]]
  stopifnot(mp[1] %in% c("verbatim", "sweep", "kmseed", "eps"))
  SEEDS <- if (mp[1] %in% c("sweep", "kmseed")) seq(as.integer(mp[2]), as.integer(mp[3])) else NA_integer_
  EPS   <- if (mp[1] == "eps") seq(as.numeric(mp[2]), as.numeric(mp[3]), by = as.numeric(mp[4])) else NA_real_
  if (mp[1] == "kmseed") stopifnot(METHOD == "K-means++")
  if (mp[1] == "eps")    stopifnot(METHOD == "DBSCAN")
  OUT_CSV <- tr3_path(switch(mp[1], verbatim = "real_baselines.csv", eps = "real_baselines_eps.csv",
                             "real_baselines_seeds.csv"))
  # Parallel processes must not append to one file (their lines interleave on G:). With
  # TR3_BASELINES_OUT=<file name> a process writes to results/tr3/<file name> instead (an
  # absolute path, "C:/..." or "/...", is used as it is); the part files are then
  # concatenated into the default file after all processes have ended.
  out_env <- Sys.getenv("TR3_BASELINES_OUT")
  if (nzchar(out_env)) OUT_CSV <- if (grepl("^([A-Za-z]:)?[/\\\\]", out_env)) out_env else tr3_path(out_env)

  SCRIPT <- SCRIPTS[[METHOD]]
  src <- readLines(here::here(SCRIPT), warn = FALSE)
  hdr <- grep("^##### \\S+( data)? #####\\s*$", src)
  blk <- sub("^##### (\\S+)( data)? #####\\s*$", "\\1", src[hdr])
  tail_line <- grep("^stopCluster\\(cl\\)", src)
  tail_line <- tail_line[tail_line > max(hdr)][1]
  ends <- c(hdr[-1] - 1, if (is.na(tail_line)) length(src) else tail_line - 1)
  absent <- setdiff(DATASETS, blk)
  if (length(absent)) stop("13_real_baselines.R: no block in ", SCRIPT, " for ", paste(absent, collapse = ", "))

  # ---- header: the script's lines before the first block --------------------
  header <- src[seq_len(hdr[1] - 1)]
  is_loader <- grepl("Real_Data_Collection\\.R", header)
  stopifnot(sum(is_loader) == 1)
  is_par <- grepl("^\\s*(cores\\s*=|cl\\s*=\\s*makeCluster|registerDoParallel)", header)
  cat(sprintf("13_real_baselines.R: %s, %s, mode %s\n  header lines: %d; loader line %d replaced; %d cluster line(s) dropped\n",
              METHOD, SCRIPT, MODE, length(header), which(is_loader), sum(is_par)))
  pdf(NULL)   # kNNdistplot() / abline() in the DBSCAN blocks
  pending <- character(0)   # header lines of an expression that spans several lines
  for (i in seq_along(header)) {
    if (is_loader[i]) {
      stopifnot(length(pending) == 0)
      for (p in c("R.matlab", "foreign", "dplyr", "readxl")) suppressPackageStartupMessages(library(p, character.only = TRUE))
      real <- tr3_load_real(DATASETS)
      for (ds in DATASETS) assign(ds, real[[ds]], envir = globalenv())
    } else if (!is_par[i]) {
      pending <- c(pending, header[i])
      ex <- tryCatch(parse(text = pending), error = function(e) e)
      if (inherits(ex, "error")) {
        # keep joining only while the expression is unfinished; any other error is real
        if (!grepl("unexpected end of input", conditionMessage(ex), fixed = TRUE)) {
          stop("13_real_baselines.R: header line ", i, " does not parse: ", conditionMessage(ex))
        }
      } else {
        eval(ex, envir = globalenv())
        pending <- character(0)
      }
    }
  }
  if (length(pending)) stop("13_real_baselines.R: the header ends inside an unfinished expression: ", pending[1])
  cat("  versions: R", as.character(getRversion()), "| igraph", as.character(packageVersion("igraph")),
      "| mclust", as.character(packageVersion("mclust")), "| cluster", as.character(packageVersion("cluster")),
      if (METHOD == "DBSCAN") paste("| dbscan", packageVersion("dbscan")) else "",
      if (METHOD == "K-means++") paste("| ClusterR", packageVersion("ClusterR")) else "",
      if (METHOD == "Spectral") paste("| kernlab", packageVersion("kernlab")) else "", "\n")

  G <- globalenv()
  for (ds in DATASETS) {
    i <- match(ds, blk)
    for (v in c("ari", "sil", "n_clusters", "eps", "sigma_used")) if (exists(v, envir = G, inherits = FALSE)) rm(list = v, envir = G)
    for (s in SEEDS) for (e in EPS) {
      text <- src[hdr[i]:ends[i]]
      if (mp[1] == "kmseed") {
        stopifnot(sum(grepl("seed = 1)", text, fixed = TRUE)) == 1)
        text <- sub("seed = 1)", sprintf("seed = %d)", s), text, fixed = TRUE)
      }
      if (mp[1] == "sweep") {
        stopifnot(sum(grepl("^set\\.seed\\(1\\)\\s*$", text)) == 1)
        text <- sub("^set\\.seed\\(1\\)\\s*$", sprintf("set.seed(%d)", s), text)
      }
      if (mp[1] == "eps") {
        stopifnot(sum(grepl("^eps = ", text)) == 1)
        text <- sub("^eps = .*$", sprintf("eps = %.17g", e), text)
      }
      t0 <- proc.time()[["elapsed"]]
      eval(parse(text = text), envir = G)
      ari <- get("ari", envir = G); sil <- get("sil", envir = G); n_clusters <- get("n_clusters", envir = G)
      row <- data.frame(method = if (mp[1] == "kmseed") "K-means++ (seed arg)" else METHOD, dataset = ds,
                        n = nrow(get("X", envir = G)), d = ncol(get("X", envir = G)),
                        ari = fmt17(ari), sil = fmt17(sil), k_hat = n_clusters,
                        seed = if (is.na(s)) "NA" else as.character(s),
                        script = SCRIPT, stringsAsFactors = FALSE)
      if (mp[1] == "eps") row$eps <- sprintf("%.17g", e)
      tr3_append(OUT_CSV, row)
      if (mp[1] != "eps") {   # the parameter the block chose, next to the main file
        par_rows <- if (METHOD == "DBSCAN") {
          data.frame(param = "eps", value = fmt17(get("eps", envir = G)))
        } else if (METHOD == "Spectral") {
          sg <- get("sigma_used", envir = G)
          data.frame(param = c("sigma", "bandwidth"), value = fmt17(c(sg, 1 / sqrt(2 * sg))))
        } else NULL
        if (!is.null(par_rows)) {
          tr3_append(sub("\\.csv$", "_params.csv", OUT_CSV),
                     data.frame(method = METHOD, dataset = ds, par_rows,
                                seed = if (is.na(s)) "NA" else as.character(s), stringsAsFactors = FALSE))
        }
      }
      cat(sprintf("  %s %s seed=%s%s: ari %.6f sil %.6f k_hat %d (%.1f s)\n", METHOD, ds,
                  if (is.na(s)) "NA" else s, if (mp[1] == "eps") paste0(" eps=", e) else "",
                  ari, sil, n_clusters, proc.time()[["elapsed"]] - t0))
    }
  }
})
