#!/usr/bin/env Rscript
# revision_experiments/tr3/37_sim_baselines.R
#
# Experiment 7 (NEW_EXPERIMENTS_PROTOCOL.md, "Experiment 7"): the five baselines DBSCAN, MST,
# K-means++, Spectral and Louvain on the simulated data of the UN-CCD paper (120 settings, 500
# replicates each) and on its noise data (d = 3, 10; four noise levels; 1000 replicates each).
# The protocol is binding; this header says what the driver does and every difference from the
# original real-data blocks. Layout as 10_sim.R (data on the main process, foreach at top level
# over the replicates of one cell, one cell appended when it finishes, then marked done).
#
# USAGE (from the clone root, in PowerShell; Rscript segfaults under Git Bash)
#   Rscript revision_experiments/tr3/37_sim_baselines.R <method|all> <part> [mode] [select]
#     method  DBSCAN | MST | K-means++ | Spectral | Louvain | all
#             all = one Rscript process per method, one after the other (as 32_hd_baselines.R):
#             never two methods attached in one session.
#     part    sim | noise | check
#     mode    run (default) | smoke
#             smoke = replicates 1..3 of every cell. It refuses to run unless TR3_OUT_DIR is set
#             and is not the folder results/tr3. Data are still generated in full (500 or 1000
#             replicates), so the random-number order is unchanged; the fingerprints of all of
#             them are compared, but only the first three are clustered and scored. At the end of a smoke
#             run the projected wall time of the full part is printed (formula below).
#     select  sim only, optional, e.g. "dist=uniform;d=2,10;n=50;k=3" (keys dist, d, n, k;
#             comma lists; default all 40 settings x 3 values of k = 120 cells per method).
#   Environment: TR3_OUT_DIR output folder (default revision_experiments/results/tr3);
#   TR3_CORES number of PSOCK workers (default detectCores()). The results do not depend on it:
#   every data set is clustered after set.seed(1) on its own.
#
# DATA (main process only; the workers never generate data)
#   sim    for each (dist in uniform, gaussian; d in 2, 3, 5, 10, 20; n in 50, 100, 200, 500):
#          s <- tr3_sim_setting("UN-CCDs", dist, d, n); data <- tr3_sim_data(s) (the lists with 2, 3
#          and 5 clusters, 500 replicates each) after tr3_source_method("UN-CCDs", "simulation"),
#          once at the start, as 10_sim.R. Cell = (method, dist, d, n, n_clusters); the true
#          labels and k come from tr3_sim_data (label_T, N). A setting is not generated when none
#          of its selected cells is pending.
#   noise  as 23_noise.R, published variant (r_noise 7, noise ball centred with sweep()): setting
#          = row "UN-CCDs gaussian <d> 200", 3 clusters of 66, 67, 67 points, the same centres,
#          Gaussian_cls() clusters, then rpoisball_unit(n_noise, d) * 7 centred on colMeans(MU),
#          set.seed(seed) once, then the four levels 0.03, 0.05, 0.10, 0.20 in this order,
#          1000 replicates each, n_noise = round(200 * level). Cell = (method, d, noise_level).
#   FINGERPRINTS fmt17(tr3_fp(x)) of every generated data set is compared, as text, with data_fp
#          of the stored rows. sim: rows of results/tr3/sim_replicates.csv with variant
#          "published" and the same dist, d, n, n_clusters, rep, from cells listed in sim_done.csv
#          (any method; keys on which methods disagree count as a mismatch). They exist for 66
#          settings (n = 50 at every d; d = 10, 20 at every n); where a setting has stored rows,
#          every replicate run must be among them and equal (fp_check "equal"); elsewhere fp_check
#          is "none". noise: rows of noise_replicates.csv with variant "published", r_noise "7",
#          center "sweep" and the same d, noise_level, rep; all must exist and be equal. Any
#          mismatch: that cell is not written (and not marked done), the driver goes on and exits
#          with status 1 at the end.
#
# METHODS (PSOCK workers; each worker attaches only the packages of its method)
#   DBSCAN     dbscan, mclust, cluster     eps = knee_eps(x, 4); dbscan(x, eps, minPts = 5)$cluster
#   MST        igraph, mclust, cluster     MST_Clustering(x, 2), MST_Clustering() from sourcing
#                                          .../MST/MST_Clustering.R on every worker
#   K-means++  ClusterR, mclust, cluster   KMeans_rcpp(x, clusters = k, initializer = "kmeans++",
#                                          verbose = F, seed = 1)$clusters
#   Spectral   kernlab, mclust, cluster    cl = specc(as.matrix(x), centers = k); cl@.Data;
#                                          sigma = kpar(kernelf(cl))$sigma
#   Louvain    igraph, mclust, cluster     S = exp(-as.matrix(dist(x))^2 / (2 * 2^2));
#                                          graph_from_adjacency_matrix(S, mode = "undirected",
#                                          weighted = T); membership(cluster_louvain(g))
#   set.seed(1) before each method on each data set. k = the true number of clusters (2, 3, 5;
#   noise: 3), used by K-means++ and Spectral only. knee_eps() is not retyped: its assignment is
#   parsed from the header of Real_Data_DBSCAN.R (k = 4 there, minPts = k + 1 = 5) and evaluated
#   on every worker. Uni-Gau_cls.R is never loaded on a worker (it masks igraph, cluster and
#   MASS functions).
#   Scoring per replicate: ind = all points (sim) or 1:200, the cluster points (noise).
#     ari  = adjustedRandIndex(label[ind], label_T) (DBSCAN label 0 = noise counts as one group)
#     sil  = mean(silhouette(label[ind], dist(x[ind, ]))[, 3]) if label[ind] has between 2 and
#            length(ind) - 1 distinct values, else NA
#     k_hat over ALL points: DBSCAN = number of distinct labels other than 0, others = number of
#            distinct labels; k_hat_block = the real-data block's own count: DBSCAN
#            length(unique(label)) - 1 (one too low when no point is noise), others = k_hat
#     n_label0 = DBSCAN: number of points with label 0 (all points); others NA
#     param = DBSCAN eps, Spectral sigma (fmt17), others NA
#     seconds = elapsed time of the method and the scoring on the worker, as a difference of
#            Sys.time() (proc.time() has a 10 to 20 ms step on Windows)
#   A method error on one replicate is caught: the row has NA metrics and the message in
#   column error; the errors of a cell are counted (n_errors). The replicate is not rerun.
#
# DIFFERENCES FROM THE ORIGINAL BLOCKS ("##### seeds data #####" of the five real-data scripts)
#   1. The block's computation is written as one function (sb_run) with the data, k and the
#      header's settings as arguments, calling the functions through their package namespaces
#      (dbscan::dbscan, ClusterR::KMeans_rcpp, kernlab::specc, igraph::..., mclust::...,
#      cluster::silhouette). The objects (distM, similarity matrix, graph, ...) are the block's;
#      the Louvain diagonal and the rest are not changed. The check part shows that this code
#      reproduces the stored real-data rows.
#   2. Packages: the scripts' headers also attach R.matlab, foreign, dplyr, readxl, ggplot2 and
#      factoextra (and the loader); the workers attach only the method's packages (list above).
#   3. ARI and Sil are computed on label[ind] (noise data: the 200 cluster points); on the real
#      data and the simulated data ind = all points, which is the block's computation.
#   4. Single-cluster silhouette guard: the block stops when silhouette() returns NA (1 or n
#      clusters); here sil is NA in that case (rule above), as 32_hd_baselines.R does.
#   5. k_hat for DBSCAN counts the labels other than 0; the block's count is k_hat_block.
#   6. Warnings of a method are muffled on the workers (they change no value).
#   7. The check part feeds the real data through the same worker set-up and the same
#      per-data-set function (sb_one) as the run, once as loaded (with dimnames) and once without.
#
# OUTPUT (OUT_DIR; numbers as fmt17 text; a cell is appended in one call, then its key to done)
#   simbase_replicates.csv   method, dist, d, n, n_clusters, mode, iteN, rep, ari, sil, k_hat,
#                            k_hat_block, n_label0, param, seconds, error, data_fp, fp_check
#   simbase_done.csv         method, dist, d, n, n_clusters, mode, iteN, n_rows, n_errors, cores,
#                            wall_secs, finished
#   noisebase_replicates.csv method, d, noise_level, r_noise, center, mode, iteN, rep, ari, sil,
#                            k_hat, k_hat_block, n_label0, param, seconds, error, data_fp, fp_check
#   noisebase_done.csv       method, d, noise_level, r_noise, center, mode, iteN, n_rows, n_errors,
#                            cores, wall_secs, finished
#   simbase_check.csv        check part, overwritten (see below)
#   iteN = the number of data sets generated for the setting (500 or 1000), also in smoke mode;
#   n_rows (done file) = the rows written. A rerun skips the cells whose key, mode included, is in
#   the done file. One log line per cell: method, cell, mean ARI (over defined values), SR (share
#   of replicates with k_hat == k, an error counts as a miss), errors, wall seconds.
#   Smoke projection per method and part, from the smoke rows and done rows in OUT_DIR:
#     projected wall = sum over cells of ( mean worker seconds per replicate x full replicates
#                      / cores + (cell wall seconds - longest replicate seconds of the cell) )
#   the second term being the measured per-cell overhead (fingerprints, dispatch, export, row
#   building, append); the cell wall time runs from the fingerprints to the append of the rows.
#   The main-process data generation of each setting is logged ("data generated in") but is not
#   part of the projection.
#
# CHECK (the gate before any full run; no data are simulated)
#   Runs the method on the 9 real sets of tr3_load_real(c("iris", "seeds", "user_knowledge",
#   "wholesale", "breast_cancer", "aggregation", "asymmetric", "R15", "D31")) and the 6 of
#   tr3_load_hd(TR3_HD_DATASETS) (features = all columns but the last, label = the last column,
#   k = number of distinct labels), each as loaded and without dimnames, through the run's worker
#   set-up and sb_one(). Compared as fmt17 text with the stored rows: ari, sil, k_hat (the stored
#   k_hat with k_hat_block) and DBSCAN eps / Spectral sigma with the params files. References:
#   the first row per (method, dataset) of real_baselines_fixed.csv and
#   real_baselines_fixed_params.csv (9 sets); the done rows of real_hd_baselines.csv and
#   real_hd_baselines_params.csv (6 sets). For DBSCAN the file also holds n_label0 and both
#   k counts for all 15 sets (reference "none", equal "info"). Rows go to OUT_DIR/simbase_check.csv
#   (overwritten; with `all`, each process writes simbase_check_<method>.csv and the parent joins
#   them). Exit status 0 only if every comparison is equal.

source(here::here("revision_experiments/tr3/tr3_lib.R"))

SB_METHODS <- c("DBSCAN", "MST", "K-means++", "Spectral", "Louvain")
SB_USAGE <- paste0("usage: Rscript revision_experiments/tr3/37_sim_baselines.R <method|all> <part> [mode] [select]\n",
                   "  method: ", paste(c(SB_METHODS, "all"), collapse = " | "), "\n",
                   "  part:   sim | noise | check\n",
                   "  mode:   run (default) | smoke\n",
                   "  select: sim only, e.g. \"dist=uniform;d=2,10;n=50;k=3\"")
SB_ARGS <- commandArgs(trailingOnly = TRUE)
if (length(SB_ARGS) < 2 || !(SB_ARGS[1] %in% c(SB_METHODS, "all")) || !(SB_ARGS[2] %in% c("sim", "noise", "check"))) {
  stop(SB_USAGE, call. = FALSE)
}
SB_METHOD <- SB_ARGS[1]
SB_PART   <- SB_ARGS[2]
SB_MODE   <- if (length(SB_ARGS) >= 3 && nzchar(SB_ARGS[3])) SB_ARGS[3] else "run"
SB_SELECT <- if (length(SB_ARGS) >= 4 && nzchar(SB_ARGS[4])) SB_ARGS[4] else ""
if (!(SB_MODE %in% c("run", "smoke"))) stop(SB_USAGE, call. = FALSE)
if (nzchar(SB_SELECT) && SB_PART != "sim") stop("37_sim_baselines.R: select applies to part sim only", call. = FALSE)
SB_OUT <- sub("[/\\\\]+$", "", Sys.getenv("TR3_OUT_DIR", TR3_RES))
if (SB_MODE == "smoke") {
  nz <- function(p) tolower(normalizePath(p, winslash = "/", mustWork = FALSE))
  if (!nzchar(Sys.getenv("TR3_OUT_DIR")) || identical(nz(SB_OUT), nz(TR3_RES))) {
    stop("37_sim_baselines.R: smoke mode needs TR3_OUT_DIR set to a folder other than results/tr3", call. = FALSE)
  }
}
SB_SIM_REP    <- file.path(SB_OUT, "simbase_replicates.csv")
SB_SIM_DONE   <- file.path(SB_OUT, "simbase_done.csv")
SB_NOISE_REP  <- file.path(SB_OUT, "noisebase_replicates.csv")
SB_NOISE_DONE <- file.path(SB_OUT, "noisebase_done.csv")
SB_CHECK_CSV  <- file.path(SB_OUT, "simbase_check.csv")
SB_CHECK_PART <- function(m) file.path(SB_OUT, sprintf("simbase_check_%s.csv", m))

# ---- all: one Rscript process per method, in sequence ----------------------------------------
if (SB_METHOD == "all") {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  me <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
  q <- function(x) paste0("\"", x, "\"")   # for the arguments only: system2() quotes the command itself
  Sys.setenv(TR3_SB_CHILD = "1")
  if (SB_PART == "check") for (f in SB_CHECK_PART(SB_METHODS)[file.exists(SB_CHECK_PART(SB_METHODS))]) file.remove(f)
  status <- 0L
  for (m in SB_METHODS) {
    cat(sprintf("37_sim_baselines.R all: starting %s %s %s\n", m, SB_PART, SB_MODE))
    st <- system2(rscript, c(q(me), q(m), q(SB_PART), q(SB_MODE), if (nzchar(SB_SELECT)) q(SB_SELECT)))
    cat(sprintf("37_sim_baselines.R all: %s exit status %d\n", m, st))
    if (st != 0L) status <- 1L
  }
  if (SB_PART == "check") {
    parts <- lapply(SB_CHECK_PART(SB_METHODS), tr3_read)
    none <- vapply(parts, is.null, logical(1))
    if (any(none)) {
      cat(sprintf("check summary: NO check file for %s (that process did not finish)\n",
                  paste(SB_METHODS[none], collapse = ", ")))
      status <- 1L
    }
    allc <- do.call(rbind, parts[!none])
    utils::write.csv(allc, SB_CHECK_CSV, row.names = FALSE)
    cmp <- allc[allc$equal %in% c("TRUE", "FALSE"), , drop = FALSE]
    for (m in unique(cmp$method)) {
      z <- cmp[cmp$method == m, , drop = FALSE]
      cat(sprintf("check summary %-9s: %d comparisons, %d equal, %d different\n", m, nrow(z),
                  sum(z$equal == "TRUE"), sum(z$equal != "TRUE")))
    }
    if (any(cmp$equal != "TRUE")) status <- 1L
    cat(sprintf("check summary: %d comparisons, %d different; written to %s\n", nrow(cmp),
                sum(cmp$equal != "TRUE"), SB_CHECK_CSV))
  }
  quit(save = "no", status = status)
}

# ---------------------------------------------------------------------------
# Worker code. Defined here in .GlobalEnv, exported to every worker by sb_make_cluster().
# ---------------------------------------------------------------------------

SB_PKGS <- list("DBSCAN"    = c("dbscan", "mclust", "cluster"),
                "MST"       = c("igraph", "mclust", "cluster"),
                "K-means++" = c("ClusterR", "mclust", "cluster"),
                "Spectral"  = c("kernlab", "mclust", "cluster"),
                "Louvain"   = c("igraph", "mclust", "cluster"))
SB_BASE <- "simulations/clustering/Algo_Compare_Clustering"
SB_DBSCAN_SCRIPT <- file.path(SB_BASE, "DBSCAN/Simulation/Real_Data/Real_Data_DBSCAN.R")
SB_MST_FILE      <- file.path(SB_BASE, "MST/MST_Clustering.R")

# Run once on every worker: attach the method's packages, in this order; define knee_eps() from
# its text (DBSCAN); source MST_Clustering.R (MST).
sb_worker_init <- function(pkgs, knee_text = NULL, mst_file = NULL) {
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  if (!is.null(knee_text)) eval(parse(text = knee_text), envir = globalenv())
  if (!is.null(mst_file)) source(mst_file, local = FALSE)
  invisible(TRUE)
}

sb_worker_info <- function(pkgs) {
  paste0("R ", getRversion(), " | ", paste(pkgs, vapply(pkgs, function(p) as.character(packageVersion(p)), ""),
                                           collapse = " "),
         " | attached: ", paste(sub("^package:", "", grep("^package:", search(), value = TRUE)), collapse = " "))
}

# The method of one block, on one data set x (matrix, one row per point); k = true number of
# clusters (K-means++ and Spectral). Returns the label of every point and the parameter.
sb_run <- function(method, x, k) {
  set.seed(1)
  switch(method,
    "DBSCAN" = {
      eps <- knee_eps(x, 4)
      label <- dbscan::dbscan(x, eps, minPts = 5)$cluster
      list(label = label, param = eps)
    },
    "MST" = list(label = MST_Clustering(x, 2), param = NA_real_),
    "K-means++" = list(label = ClusterR::KMeans_rcpp(x, clusters = k, initializer = "kmeans++", verbose = F,
                                                     seed = 1)$clusters, param = NA_real_),
    "Spectral" = {
      cl <- kernlab::specc(as.matrix(x), centers = k)
      list(label = cl@.Data, param = kernlab::kpar(kernlab::kernelf(cl))$sigma)
    },
    "Louvain" = {
      distM <- as.matrix(dist(x))
      S <- exp(-distM^2 / (2 * 2^2))
      g <- igraph::graph_from_adjacency_matrix(S, mode = "undirected", weighted = T)
      list(label = igraph::membership(igraph::cluster_louvain(g)), param = NA_real_)
    },
    stop("sb_run(): unknown method ", method))
}

sb_score <- function(method, label, param, x, label_T, ind) {
  lab <- label[ind]
  nd <- length(unique(lab))
  ari <- mclust::adjustedRandIndex(lab, label_T)
  sil <- if (nd >= 2 && nd <= length(ind) - 1) {
    mean(cluster::silhouette(lab, stats::dist(x[ind, , drop = FALSE]))[, 3])
  } else NA_real_
  ul <- unique(label)
  if (method == "DBSCAN") {
    k_hat <- sum(ul != 0); k_hat_block <- length(ul) - 1; n_label0 <- sum(label == 0)
  } else {
    k_hat <- length(ul); k_hat_block <- k_hat; n_label0 <- NA_real_
  }
  list(ari = as.numeric(ari), sil = as.numeric(sil), k_hat = as.numeric(k_hat),
       k_hat_block = as.numeric(k_hat_block), n_label0 = as.numeric(n_label0),
       param = as.numeric(param), error = "")
}

# One replicate (or one real data set): the per-data-set function of the run and of the check.
# nc = 0: ind = all points; nc > 0: ind = 1:nc.
sb_one <- function(x, method, k, label_T, nc = 0L) {
  t0 <- as.numeric(Sys.time())   # proc.time() elapsed has a 10 to 20 ms step on Windows
  ind <- if (nc > 0L) seq_len(nc) else seq_len(nrow(x))
  out <- tryCatch({
    r <- suppressWarnings(sb_run(method, x, k))
    suppressWarnings(sb_score(method, r$label, r$param, x, label_T, ind))
  }, error = function(e) {
    list(ari = NA_real_, sil = NA_real_, k_hat = NA_real_, k_hat_block = NA_real_, n_label0 = NA_real_,
         param = NA_real_, error = substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 300))
  })
  out$seconds <- as.numeric(Sys.time()) - t0
  out
}

# ---------------------------------------------------------------------------
# Main-process helpers
# ---------------------------------------------------------------------------

# knee_eps() text and the header's k, from the header of the original DBSCAN script.
SB_KNEE_TEXT <- NULL
if (SB_METHOD == "DBSCAN") {
  src <- readLines(here::here(SB_DBSCAN_SCRIPT), warn = FALSE)
  hdr <- grep("^##### \\S+( data)? #####\\s*$", src)
  ex <- parse(text = src[seq_len(hdr[1] - 1)], keep.source = TRUE)
  is_asg <- function(e, nm) {
    is.call(e) && (identical(e[[1]], as.name("=")) || identical(e[[1]], as.name("<-"))) &&
      is.name(e[[2]]) && identical(as.character(e[[2]]), nm)
  }
  i_k <- which(vapply(ex, is_asg, logical(1), nm = "k"))
  i_f <- which(vapply(ex, is_asg, logical(1), nm = "knee_eps"))
  stopifnot(length(i_k) == 1, length(i_f) == 1)
  if (!identical(eval(ex[[i_k]][[3]]), 4)) stop("37_sim_baselines.R: the header's k is not 4 (minPts = 5 assumed)")
  SB_KNEE_TEXT <- paste(as.character(attr(ex, "srcref")[[i_f]]), collapse = "\n")
  if (!is.function(eval(parse(text = SB_KNEE_TEXT)[[1]][[3]]))) stop("37_sim_baselines.R: knee_eps text is not a function")
  rm(src, hdr, ex, i_k, i_f)
}

sb_make_cluster <- function(method, cores) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterCall(cl, sb_worker_init, SB_PKGS[[method]], SB_KNEE_TEXT,
                        if (method == "MST") here::here(SB_MST_FILE) else NULL)
  parallel::clusterExport(cl, c("sb_run", "sb_score", "sb_one"), envir = globalenv())
  cat("  worker 1: ", parallel::clusterCall(cl[1], sb_worker_info, SB_PKGS[[method]])[[1]], "\n", sep = "")
  cl
}

# a result table of fmt17 text columns for the rows of one cell
sb_rows <- function(key, iteN, reps, r, fp, fp_check) {
  num <- function(f) vapply(r, function(z) as.numeric(z[[f]]), numeric(1))
  data.frame(key, iteN = iteN, rep = reps,
             ari = fmt17(num("ari")), sil = fmt17(num("sil")), k_hat = fmt17(num("k_hat")),
             k_hat_block = fmt17(num("k_hat_block")), n_label0 = fmt17(num("n_label0")),
             param = fmt17(num("param")), seconds = fmt17(num("seconds")),
             error = vapply(r, function(z) as.character(z$error), ""), data_fp = fp, fp_check = fp_check,
             stringsAsFactors = FALSE)
}

sb_done_row <- function(key, iteN, rows, cores, wall) {
  data.frame(key, iteN = iteN, n_rows = nrow(rows), n_errors = sum(nzchar(rows$error)), cores = cores,
             wall_secs = fmt17(wall), finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             stringsAsFactors = FALSE)
}

sb_log_cell <- function(label, rows, k, wall, fp_note) {
  ari <- suppressWarnings(as.numeric(rows$ari)); kh <- suppressWarnings(as.numeric(rows$k_hat))
  cat(sprintf("  %s: mean ARI %.4f, SR %.3f, errors %d, wall %.1f s%s\n", label,
              if (all(is.na(ari))) NA_real_ else mean(ari, na.rm = TRUE),
              sum(kh == k, na.rm = TRUE) / nrow(rows), sum(nzchar(rows$error)), wall, fp_note))
}

# Stored fingerprints. Each returns a named character vector (key -> data_fp) and, for sim, the
# keys of the settings that have stored rows.
sb_stored_sim <- function() {
  f <- file.path(TR3_RES, "sim_replicates.csv")
  sr <- utils::read.csv(f, colClasses = "character")
  done <- tr3_read(file.path(TR3_RES, "sim_done.csv"))
  cellk <- function(df) paste(df$variant, df$method, df$dist, df$d, df$n, df$n_clusters, sep = "|")
  sr <- sr[sr$variant == "published" & cellk(sr) %in% cellk(done), , drop = FALSE]
  key <- paste(sr$dist, sr$d, sr$n, sr$n_clusters, sr$rep, sep = "|")
  u <- unique(data.frame(key = key, fp = sr$data_fp, stringsAsFactors = FALSE))
  amb <- unique(u$key[duplicated(u$key)])
  u <- u[!duplicated(u$key), , drop = FALSE]
  fp <- stats::setNames(u$fp, u$key)
  fp[amb] <- "AMBIGUOUS"
  list(fp = fp, settings = unique(paste(sr$dist, sr$d, sr$n, sr$n_clusters, sep = "|")), n_ambiguous = length(amb))
}

sb_stored_noise <- function() {
  nr <- utils::read.csv(file.path(TR3_RES, "noise_replicates.csv"), colClasses = "character")
  done <- tr3_read(file.path(TR3_RES, "noise_done.csv"))
  cellk <- function(df) paste(df$variant, df$d, df$noise_level, df$r_noise, df$center, sep = "|")
  nr <- nr[nr$variant == "published" & nr$r_noise == "7" & nr$center == "sweep" & cellk(nr) %in% cellk(done), , drop = FALSE]
  key <- paste(nr$d, nr$noise_level, nr$rep, sep = "|")
  stopifnot(!anyDuplicated(key))
  stats::setNames(nr$data_fp, key)
}

# Smoke projection of the full part (formula in the header).
sb_projection <- function(rep_csv, done_csv, cell_cols, method, cells, cores, full_cells) {
  r <- tr3_read(rep_csv); d <- tr3_read(done_csv)
  if (is.null(r) || is.null(d)) { cat("  projection: no smoke rows\n"); return(invisible(NULL)) }
  r <- r[r$method == method & r$mode == "smoke", , drop = FALSE]
  d <- d[d$method == method & d$mode == "smoke", , drop = FALSE]
  ck <- function(df) do.call(paste, c(df[cell_cols], sep = "|"))
  r$cell <- ck(r); d$cell <- ck(d)
  d <- d[d$cell %in% cells & !duplicated(d$cell, fromLast = TRUE), , drop = FALSE]
  r <- r[r$cell %in% d$cell, , drop = FALSE]
  secs <- as.numeric(r$seconds)
  mean_s <- tapply(secs, r$cell, mean, na.rm = TRUE)[d$cell]
  max_s  <- tapply(secs, r$cell, max,  na.rm = TRUE)[d$cell]
  full   <- as.numeric(d$iteN)
  over   <- pmax(as.numeric(d$wall_secs) - max_s, 0)
  comp   <- mean_s * full / cores
  tot    <- sum(comp) + sum(over)
  cat(sprintf(paste0("  PROJECTION %s (%d of %d cells of the full part): sum over cells of (mean worker seconds per ",
                     "replicate x full replicates / cores + (cell wall - longest replicate)), cores = %d\n",
                     "    compute %.1f s + overhead %.2f s = %.1f s = %.3f h%s\n"),
              method, nrow(d), full_cells, cores, sum(comp), sum(over), tot, tot / 3600,
              if (nrow(d) < full_cells) " (selected cells only, not the full part)" else ""))
  invisible(tot)
}

SB_FAIL <- FALSE
SB_T0 <- Sys.time()
SB_CORES <- tr3_cores()
SB_LEVELS <- c(0.03, 0.05, 0.10, 0.20)

# ============================== sim ======================================================
if (SB_PART == "sim") {
  grid <- list(dist = c("uniform", "gaussian"), d = c("2", "3", "5", "10", "20"),
               n = c("50", "100", "200", "500"), k = c("2", "3", "5"))
  sel <- grid
  if (nzchar(SB_SELECT)) {
    for (p in strsplit(SB_SELECT, ";", fixed = TRUE)[[1]]) {
      kv <- strsplit(trimws(p), "=", fixed = TRUE)[[1]]
      if (length(kv) != 2 || !(kv[1] %in% names(sel))) stop("37_sim_baselines.R: bad select item '", p, "'")
      v <- trimws(strsplit(kv[2], ",", fixed = TRUE)[[1]])
      if (!all(v %in% grid[[kv[1]]])) stop("37_sim_baselines.R: select value not in the grid: ", kv[1], "=", kv[2])
      sel[[kv[1]]] <- v
    }
  }
  tr3_source_method("UN-CCDs", "simulation")
  stored <- sb_stored_sim()
  cat(sprintf("37_sim_baselines.R: %s sim, mode %s, out=%s, cores %d, select '%s'\n  stored fingerprints: %d keys in %d settings, %d ambiguous\n",
              SB_METHOD, SB_MODE, SB_OUT, SB_CORES, SB_SELECT, length(stored$fp), length(stored$settings), stored$n_ambiguous))
  cl <- sb_make_cluster(SB_METHOD, SB_CORES)
  doParallel::registerDoParallel(cl)
  sel_cells <- character(0)
  for (SB_DIST in sel$dist) for (SB_D in sel$d) for (SB_N in sel$n) {
    ks <- sel$k
    keys <- lapply(ks, function(k) list(method = SB_METHOD, dist = SB_DIST, d = SB_D, n = SB_N, n_clusters = k, mode = SB_MODE))
    sel_cells <- c(sel_cells, vapply(ks, function(k) paste(SB_DIST, SB_D, SB_N, k, sep = "|"), ""))
    pend <- ks[!vapply(keys, has_result, logical(1), csv_path = SB_SIM_DONE)]
    if (!length(pend)) { cat(sprintf("  %s d=%s n=%s: all selected cells done, skipped\n", SB_DIST, SB_D, SB_N)); next }
    s <- tr3_sim_setting("UN-CCDs", SB_DIST, SB_D, SB_N)
    tg <- Sys.time()
    data <- tr3_sim_data(s)
    cat(sprintf("  %s d=%s n=%s: data generated in %.1f s\n", SB_DIST, SB_D, SB_N, as.numeric(difftime(Sys.time(), tg, units = "secs"))))
    iteN <- as.integer(s$iteN)
    for (SB_KC in pend) {
      key <- keys[[match(SB_KC, ks)]]
      g <- data[[SB_KC]]
      reps <- if (SB_MODE == "smoke") seq_len(min(3L, iteN)) else seq_len(iteN)
      t0 <- Sys.time()   # the cell's wall time covers the fingerprints, the foreach and the append
      # the fingerprints of ALL generated data sets are compared, also in smoke mode
      all_reps <- seq_along(g$data)
      fp_all <- fmt17(vapply(g$data, tr3_fp, numeric(1)))
      sk <- paste(SB_DIST, SB_D, SB_N, SB_KC, sep = "|")
      if (sk %in% stored$settings) {
        st <- unname(stored$fp[paste(sk, all_reps, sep = "|")])
        bad <- is.na(st) | st != fp_all
        fp_check <- rep("equal", length(reps))
      } else {
        bad <- rep(FALSE, length(all_reps)); fp_check <- rep("none", length(reps))
      }
      data_list <- g$data[reps]
      fp <- fp_all[reps]
      label <- sprintf("%s %s d=%s n=%s k=%s", SB_METHOD, SB_DIST, SB_D, SB_N, SB_KC)
      if (any(bad)) {
        cat(sprintf("  %s: FINGERPRINT MISMATCH in %d of %d data sets (first: rep %s); cell not written\n",
                    label, sum(bad), length(all_reps), paste(head(all_reps[bad], 5), collapse = ",")), file = stderr())
        SB_FAIL <- TRUE
        next
      }
      SB_NCL <- as.numeric(SB_KC)
      SB_LABEL_T <- g$label_T
      r <- foreach(x = data_list) %dopar% sb_one(x, SB_METHOD, SB_NCL, SB_LABEL_T, 0L)
      rows <- sb_rows(key, iteN, reps, r, fp, fp_check)
      tr3_append(SB_SIM_REP, rows)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      tr3_append(SB_SIM_DONE, sb_done_row(key, iteN, rows, SB_CORES, wall))
      sb_log_cell(label, rows, SB_NCL, wall, sprintf(", fp %s", fp_check[1]))
    }
  }
  parallel::stopCluster(cl)
  if (SB_MODE == "smoke") {
    sb_projection(SB_SIM_REP, SB_SIM_DONE, c("dist", "d", "n", "n_clusters"), SB_METHOD, sel_cells, SB_CORES, 120L)
  }
}

# ============================== noise ====================================================
if (SB_PART == "noise") {
  tr3_source_method("UN-CCDs", "simulation")
  stored_fp <- sb_stored_noise()
  cat(sprintf("37_sim_baselines.R: %s noise, mode %s, out=%s, cores %d\n  stored fingerprints: %d\n",
              SB_METHOD, SB_MODE, SB_OUT, SB_CORES, length(stored_fp)))
  cl <- sb_make_cluster(SB_METHOD, SB_CORES)
  doParallel::registerDoParallel(cl)
  NC <- 200L
  sel_cells <- character(0)
  for (SB_D in c("3", "10")) {
    keys <- lapply(SB_LEVELS, function(l) list(method = SB_METHOD, d = SB_D, noise_level = as.character(l),
                                               r_noise = "7", center = "sweep", mode = SB_MODE))
    sel_cells <- c(sel_cells, vapply(SB_LEVELS, function(l) paste(SB_D, as.character(l), sep = "|"), ""))
    pend <- which(!vapply(keys, has_result, logical(1), csv_path = SB_NOISE_DONE))
    if (!length(pend)) { cat(sprintf("  d=%s: all noise levels done, skipped\n", SB_D)); next }
    # the setting and the data, as 23_noise.R (published variant); all four levels, in order, from one seed
    s <- tr3_sim_setting("UN-CCDs", "gaussian", SB_D, "200")
    dd <- as.numeric(SB_D)
    n_cls <- tr3_sizes(3, NC, as.numeric(s$n1_adj_3cls))
    MU <- rbind(rep(3, dd), c(9, rep(3, dd - 1)), c(3, 9, rep(3, dd - 2)))
    sig <- as.numeric(s$sig); rad <- c(as.numeric(s$r_min), as.numeric(s$r_max))
    label_T <- rep(0:2, n_cls)
    iteN <- 1000L
    noise_set <- function(i, n_noise) {
      x <- Gaussian_cls(n_cls = n_cls, d = dd, centers = MU, sigma = sig, noise = FALSE, rad_range = rad)
      z <- rpoisball_unit(n_noise, dd) * 7
      z <- sweep(z, 2, colMeans(MU), "+")
      rbind(x, z)
    }
    tg <- Sys.time()
    set.seed(as.numeric(s$seed))
    data <- lapply(SB_LEVELS, function(l) {
      nn <- round(NC * l)
      lapply(seq_len(iteN), noise_set, n_noise = nn)
    })
    names(data) <- as.character(SB_LEVELS)
    cat(sprintf("  d=%s: noise data generated in %.1f s\n", SB_D, as.numeric(difftime(Sys.time(), tg, units = "secs"))))
    for (i in pend) {
      key <- keys[[i]]
      lv <- key$noise_level
      reps <- if (SB_MODE == "smoke") seq_len(3L) else seq_len(iteN)
      t0 <- Sys.time()   # the cell's wall time covers the fingerprints, the foreach and the append
      all_reps <- seq_len(iteN)   # all generated data sets are compared, also in smoke mode
      fp_all <- fmt17(vapply(data[[lv]], tr3_fp, numeric(1)))
      st <- unname(stored_fp[paste(SB_D, lv, all_reps, sep = "|")])
      bad <- is.na(st) | st != fp_all
      data_list <- data[[lv]][reps]
      fp <- fp_all[reps]
      label <- sprintf("%s noise d=%s l=%s", SB_METHOD, SB_D, lv)
      if (any(bad)) {
        cat(sprintf("  %s: FINGERPRINT MISMATCH or missing in %d of %d data sets (first: rep %s); cell not written\n",
                    label, sum(bad), length(all_reps), paste(head(all_reps[bad], 5), collapse = ",")), file = stderr())
        SB_FAIL <- TRUE
        next
      }
      SB_NCL <- 3
      SB_LABEL_T <- label_T
      r <- foreach(x = data_list) %dopar% sb_one(x, SB_METHOD, SB_NCL, SB_LABEL_T, 200L)
      rows <- sb_rows(key, iteN, reps, r, fp, rep("equal", length(reps)))
      tr3_append(SB_NOISE_REP, rows)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      tr3_append(SB_NOISE_DONE, sb_done_row(key, iteN, rows, SB_CORES, wall))
      sb_log_cell(label, rows, SB_NCL, wall, ", fp equal")
    }
  }
  parallel::stopCluster(cl)
  if (SB_MODE == "smoke") {
    sb_projection(SB_NOISE_REP, SB_NOISE_DONE, c("d", "noise_level"), SB_METHOD, sel_cells, SB_CORES, 8L)
  }
}

# ============================== check ====================================================
if (SB_PART == "check") {
  library(foreach); library(doParallel)   # no Uni-Gau_cls.R and no simulated data here
  source(here::here("revision_experiments/tr3/30_hd_load.R"))
  REAL9 <- c("iris", "seeds", "user_knowledge", "wholesale", "breast_cancer", "aggregation", "asymmetric", "R15", "D31")
  real <- tr3_load_real(REAL9)
  hd <- tr3_load_hd(TR3_HD_DATASETS, verbose = FALSE)
  sets <- c(real[REAL9], hd[TR3_HD_DATASETS])
  group <- c(rep("real9", length(REAL9)), rep("hd6", length(TR3_HD_DATASETS)))
  names(group) <- names(sets)
  cat(sprintf("37_sim_baselines.R: %s check, out=%s, cores %d, %d data sets\n", SB_METHOD, SB_OUT, SB_CORES, length(sets)))

  # references: first row per (method, dataset) of the fixed file; done rows of the hd file
  rd <- function(rows) {
    rows <- rows[rows$method == SB_METHOD, , drop = FALSE]
    rows[!duplicated(rows$dataset), , drop = FALSE]
  }
  ref9 <- rd(tr3_read(file.path(TR3_RES, "real_baselines_fixed.csv")))
  par9 <- tr3_read(file.path(TR3_RES, "real_baselines_fixed_params.csv"))
  ref6 <- rd(tr3_read_done(file.path(TR3_RES, "real_hd_baselines.csv"), file.path(TR3_RES, "real_hd_baselines_done.csv"),
                           c("method", "dataset")))
  par6 <- tr3_read(file.path(TR3_RES, "real_hd_baselines_params.csv"))
  ref_name <- c(real9 = "real_baselines_fixed.csv", hd6 = "real_hd_baselines.csv")
  stopifnot(all(REAL9 %in% ref9$dataset), all(TR3_HD_DATASETS %in% ref6$dataset))

  tasks <- list()
  for (ds in names(sets)) {
    D0 <- as.matrix(sets[[ds]])
    for (form in c("asloaded", "plain")) {
      D <- D0
      if (form == "plain") dimnames(D) <- NULL
      yc <- ncol(D)
      tasks[[length(tasks) + 1]] <- list(ds = ds, form = form, x = D[, -yc, drop = FALSE], y = D[, yc],
                                         k = length(unique(D[, yc])))
    }
  }
  cl <- sb_make_cluster(SB_METHOD, SB_CORES)
  doParallel::registerDoParallel(cl)
  t0 <- Sys.time()
  res <- foreach(tk = tasks) %dopar% sb_one(tk$x, SB_METHOD, tk$k, tk$y, 0L)
  parallel::stopCluster(cl)
  cat(sprintf("  %d data sets x 2 forms in %.1f s\n", length(sets), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  txt <- function(x) ifelse(is.na(x), "NA", as.character(x))
  chk <- list()
  add <- function(tk, quantity, reference, stored, got, equal) {
    chk[[length(chk) + 1]] <<- data.frame(method = SB_METHOD, group = group[[tk$ds]], dataset = tk$ds, form = tk$form,
      n = nrow(tk$x), d = ncol(tk$x), K = tk$k, reference = reference, quantity = quantity,
      stored = stored, reproduced = got, equal = equal, stringsAsFactors = FALSE)
  }
  for (i in seq_along(tasks)) {
    tk <- tasks[[i]]; r <- res[[i]]
    refs <- if (group[[tk$ds]] == "real9") ref9 else ref6
    pars <- if (group[[tk$ds]] == "real9") par9 else par6
    rf <- ref_name[[group[[tk$ds]]]]
    row <- refs[refs$dataset == tk$ds, , drop = FALSE]
    if (nzchar(r$error)) {
      add(tk, "error", rf, "NA", r$error, "FALSE")
      next
    }
    cmp <- function(q, stored, got) add(tk, q, rf, txt(stored), got, as.character(identical(txt(stored), got)))
    cmp("ari", row$ari, fmt17(r$ari))
    cmp("sil", row$sil, fmt17(r$sil))
    cmp("k_hat", row$k_hat, fmt17(r$k_hat_block))
    pname <- c(DBSCAN = "eps", Spectral = "sigma")[SB_METHOD]
    if (!is.na(pname)) {
      pv <- pars$value[pars$method == SB_METHOD & pars$dataset == tk$ds & pars$param == pname]
      if (length(pv) < 1) stop("37_sim_baselines.R check: no stored ", pname, " for ", tk$ds)
      cmp(pname, pv[1], fmt17(r$param))
    }
    if (SB_METHOD == "DBSCAN" && tk$form == "asloaded") {
      for (q in c("n_label0", "k_hat", "k_hat_block")) {
        add(tk, if (q == "k_hat") "k_hat_excl0" else q, "none", "NA", fmt17(r[[q]]), "info")
      }
    }
  }
  chk <- do.call(rbind, chk)
  dir.create(SB_OUT, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(chk, SB_CHECK_PART(SB_METHOD), row.names = FALSE)
  if (!nzchar(Sys.getenv("TR3_SB_CHILD"))) utils::write.csv(chk, SB_CHECK_CSV, row.names = FALSE)
  z <- chk[chk$equal %in% c("TRUE", "FALSE"), , drop = FALSE]
  for (i in seq_len(nrow(z))) {
    if (z$equal[i] != "TRUE") {
      cat(sprintf("  DIFFERENT %-13s %-8s %-6s stored %-22s reproduced %s\n", z$dataset[i], z$form[i], z$quantity[i],
                  z$stored[i], z$reproduced[i]))
    }
  }
  if (SB_METHOD == "DBSCAN") {
    inf <- chk[chk$equal == "info", , drop = FALSE]
    w <- reshape(inf[, c("dataset", "quantity", "reproduced")], idvar = "dataset", timevar = "quantity", direction = "wide")
    cat("  DBSCAN noise points and cluster counts on the 15 sets (n_label0, k_hat without label 0, block's k_hat):\n")
    for (i in seq_len(nrow(w))) {
      cat(sprintf("    %-14s n_label0 %5s  k_hat %3s  k_hat_block %3s%s\n", w$dataset[i], w$reproduced.n_label0[i],
                  w$reproduced.k_hat_excl0[i], w$reproduced.k_hat_block[i],
                  if (w$reproduced.n_label0[i] == "0") "   <-- no noise point: the block's k_hat is one too low" else ""))
    }
  }
  ok <- nrow(z) > 0 && all(z$equal == "TRUE")
  cat(sprintf("37_sim_baselines.R check: %s: %d comparisons, %d equal, %d different: %s\n", SB_METHOD, nrow(z),
              sum(z$equal == "TRUE"), sum(z$equal != "TRUE"), if (ok) "PASS" else "FAIL"))
  quit(save = "no", status = if (ok) 0L else 1L)
}

cat(sprintf("37_sim_baselines.R: %s %s %s finished in %.1f s%s\n", SB_METHOD, SB_PART, SB_MODE,
            as.numeric(difftime(Sys.time(), SB_T0, units = "secs")),
            if (SB_FAIL) "; at least one cell failed its fingerprint check" else ""))
quit(save = "no", status = if (SB_FAIL) 1L else 0L)
