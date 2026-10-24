#!/usr/bin/env Rscript
# revision_experiments/tr3/tr3_lib.R
#
# Shared code for the TR3 drivers (the UN-CCD clustering paper, "Clustering
# with Uniformity- and Neighbor-Based Random Geometric Graphs"). Sourced by
# 10_sim.R, 11_real.R and 12_real_reference.R.
#
# This file defines helpers only. It does not source any method code: each
# driver process runs ONE method and sources it with tr3_source_method(), in
# the order its original script did. The original scripts never had two
# methods loaded at once, and Uni-Gau_cls.R (spatstat, hypercube -> plotly)
# masks functions of igraph, cluster and MASS, so the load order is part of
# what is being reproduced.
#
# Reads methods/, R/, data/ and simulations/; writes only to
# revision_experiments/results/tr3/.
#
# Provides
#   tr3_sim_setting(), tr3_real_setting()   rows of the declared settings tables
#   tr3_source_method(method, context)      load one method as its script did
#   tr3_load_table(file)                    quantile table, generation off
#   ccd_cluster(method, X, par)             the one interface to the three methods
#   tr3_sim_data(s)                         the scripts' data, in their RNG order
#   tr3_load_real(names)                    the real-data loader, patched in memory
#   tr3_append(), tr3_read(), fmt17()       per-cell CSV checkpointing

TR3_ROOT <- here::here()
TR3_DIR  <- file.path(TR3_ROOT, "revision_experiments", "tr3")
TR3_RES  <- file.path(TR3_ROOT, "revision_experiments", "results", "tr3")
TR3_METHODS <- c("UN-CCDs", "RK-CCDs", "KS-CCDs")

# A missing or short quantile table must stop the run. Generating one on the
# spot (R/ccds/quantile_table.R) reproduces published numbers only up to Monte
# Carlo error. The environment variable is set before any makeCluster(), so
# the PSOCK workers inherit it.
Sys.setenv(CCD_QUANTILE_GENERATE = "FALSE")
options(ccd.quantile.generate = FALSE)

# append_result() and has_result() from shared/harness.R. Sourcing the harness
# itself would load the outlier-detection methods into .GlobalEnv, so only
# these two definitions are evaluated.
local({
  exprs <- parse(file.path(TR3_ROOT, "revision_experiments", "shared", "harness.R"),
                 keep.source = FALSE)
  want <- c("append_result", "has_result")
  for (e in exprs) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% want) {
      eval(e, envir = globalenv())
    }
  }
  miss <- want[!vapply(want, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(miss)) stop("tr3_lib.R: not found in shared/harness.R: ", paste(miss, collapse = ", "))
})

# Numbers go to the CSVs as text with 17 significant digits, which round-trips
# a double exactly. write.table() would keep 15, and a mean recomputed from
# 15-digit values can differ from the logged mean in its last printed digit.
fmt17 <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))

tr3_cores <- function() as.integer(Sys.getenv("TR3_CORES", parallel::detectCores()))

# ---------------------------------------------------------------------------
# 1. Declared settings
# ---------------------------------------------------------------------------
#
# settings_simulation.csv and settings_realdata.csv hold the settings exactly
# as the original scripts ran them, not as the paper states them;
# 00_check_settings.R re-reads all 120 simulation drivers and the three
# real-data scripts and stops if any value differs. All columns are read as
# text and converted where they are used.

tr3_settings <- function(kind = c("simulation", "realdata")) {
  kind <- match.arg(kind)
  utils::read.csv(file.path(TR3_DIR, sprintf("settings_%s.csv", kind)),
                  stringsAsFactors = FALSE, colClasses = "character")
}

tr3_sim_setting <- function(method, dist, d, n) {
  s <- tr3_settings("simulation")
  r <- s[s$method == method & s$dist == dist & s$d == as.character(d) & s$n == as.character(n), ]
  if (nrow(r) != 1) {
    stop(sprintf("tr3_sim_setting(): %d rows for method=%s dist=%s d=%s n=%s",
                 nrow(r), method, dist, d, n))
  }
  as.list(r)
}

tr3_real_setting <- function(method, dataset) {
  s <- tr3_settings("realdata")
  r <- s[s$method == method & s$dataset == dataset, ]
  if (nrow(r) != 1) {
    stop(sprintf("tr3_real_setting(): %d rows for method=%s dataset=%s", nrow(r), method, dataset))
  }
  as.list(r)
}

# ---------------------------------------------------------------------------
# 2. Loading one method the way its script did
# ---------------------------------------------------------------------------

TR3_WRAPPER <- c("UN-CCDs" = "methods/clustering/UN-CCDs.R",
                 "RK-CCDs" = "methods/clustering/RK-CCDs.R",
                 "KS-CCDs" = "methods/clustering/KS-CCDs.R")

tr3_library <- function(pkgs) {
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
}

#' Source one method and the packages its script loaded, in the same order.
#'   simulation: <wrapper>, Uni-Gau_cls.R, quantile_table.R (not for KS), then
#'               mclust, parallel, doParallel, MASS, igraph
#'   realdata:   <wrapper>, the loader's packages (R.matlab, foreign, dplyr;
#'               readxl is attached only if the knowledge data is loaded), then
#'               mclust, parallel, doParallel, MASS, cluster, igraph
tr3_source_method <- function(method, context = c("simulation", "realdata")) {
  context <- match.arg(context)
  if (!(method %in% TR3_METHODS)) stop("tr3_source_method(): unknown method ", method)
  source(here::here(TR3_WRAPPER[[method]]))
  if (context == "simulation") {
    source(here::here("R/general_functions/Uni-Gau_cls.R"))
    if (method != "KS-CCDs") source(here::here("R/ccds/quantile_table.R"))
    tr3_library(c("mclust", "parallel", "doParallel", "MASS", "igraph"))
  } else {
    tr3_library(c("R.matlab", "foreign", "dplyr"))
    tr3_library(c("mclust", "parallel", "doParallel", "MASS", "cluster", "igraph"))
  }
  invisible(method)
}

#' Load a quantile table by the file name its script used. RK tables carry the
#' raw Monte Carlo draws (Kest.m, about 0.4 GB in memory), which
#' ccd.Kest.edge.quantile(), .ccd_qt_extent() and .ccd_qt_grow() never read;
#' they use only $r and $quan. Kest.m is dropped so that the table can be sent
#' to every worker.
tr3_load_table <- function(file) {
  variant <- if (startsWith(file, "RK-")) "RK" else if (startsWith(file, "NN-")) "NN" else
    stop("tr3_load_table(): not a quantile-table file name: ", file)
  path <- here::here(sprintf("R/%s-test_quantile/%s", variant, file))
  if (!file.exists(path)) stop("tr3_load_table(): missing table (generation is off): ", path)
  e <- new.env()
  load_quantile_table(path, envir = e)
  simul <- get("simul", envir = e)
  if (variant == "RK") simul$Kest.m <- NULL
  simul
}

# ---------------------------------------------------------------------------
# 3. The one interface
# ---------------------------------------------------------------------------
#
# ccd_cluster() calls the method's own wrapper with exactly the arguments the
# scripts used:
#   UN-CCDs  UNCCD_clustering(datax, simul, method = <direction>)
#   RK-CCDs  RKCCD_clustering(datax, simul, quant = <level>)
#   KS-CCDs  KSCCD_clustering(datax, m = <scaled density>)
# and returns the label of every point, the method's own average silhouette
# (si), its number of clusters (si.ind) and the elapsed seconds.

ccd_cluster <- function(method, X, par) {
  t0 <- proc.time()[["elapsed"]]
  res <- switch(method,
    "UN-CCDs" = UNCCD_clustering(datax = X, simul = par$simul, method = par$direction),
    "RK-CCDs" = RKCCD_clustering(datax = X, simul = par$simul, quant = par$quant),
    "KS-CCDs" = KSCCD_clustering(datax = X, m = par$m),
    stop("ccd_cluster(): unknown method ", method))
  list(label = res$label, si = res$si, si.ind = res$si.ind,
       seconds = proc.time()[["elapsed"]] - t0)
}

#' One simulated replicate, run on a worker: ARI against the true labels, the
#' method's si and si.ind (what the simulation scripts averaged), seconds.
tr3_one <- function(x, label_T, method, par) {
  r <- ccd_cluster(method, x, par)
  if (length(r$si) != 1 || length(r$si.ind) != 1) stop("tr3_one(): si or si.ind is not a scalar")
  c(ari = mclust::adjustedRandIndex(r$label, label_T), sil = r$si,
    k_hat = r$si.ind, seconds = r$seconds)
}

# ---------------------------------------------------------------------------
# 4. Simulated data, in the scripts' random-number order
# ---------------------------------------------------------------------------
#
# Every script calls set.seed(<seed>) once and then generates the 2-, 3- and
# 5-cluster lists of iteN data sets, in that order; nothing between them
# draws random numbers on the main process. tr3_sim_data() therefore always
# generates all three lists, even when a rerun skips finished cells: skipping
# one would shift every later data set without any error. With a different iteN
# the 3- and 5-cluster data differ too (as in the original 50-repetition RK-CCD
# runs at d = 3, n = 500; the scripts use 500 since 2026-10-18).

#' Cluster centres exactly as the scripts build them.
tr3_centers <- function(N, d, cls_dis) {
  mu1 <- rep(3, d)
  mu2 <- c(3 + cls_dis, rep(3, d - 1))
  if (N == 2) return(rbind(mu1, mu2))
  mu3 <- c(3 + cls_dis, 3 + cls_dis, rep(3, d - 2))
  if (N == 3) return(rbind(mu1, mu2, mu3))
  mu4 <- c(3, 3 + cls_dis - 0.5, rep(3, d - 2))
  mu5 <- c(3 + 2 * cls_dis - 0.5, rep(3, d - 1))
  if (N == 5) return(rbind(mu1, mu2, mu3, mu4, mu5))
  stop("tr3_centers(): N must be 2, 3 or 5")
}

#' Cluster sizes: round(n/N) each, except that the first of three clusters is
#' round(n/3) + n1_adj_3cls (+1 at n = 100, -1 otherwise, per script).
tr3_sizes <- function(N, n, n1_adj_3cls) {
  if (N == 3) return(c(round(n / N) + n1_adj_3cls, round(n / N), round(n / N)))
  rep(round(n / N), N)
}

tr3_sim_data <- function(s) {
  n <- as.numeric(s$n); d <- as.numeric(s$d); iteN <- as.numeric(s$iteN)
  cls_dis <- as.numeric(s$cls_dis); r_min <- as.numeric(s$r_min); r_max <- as.numeric(s$r_max)
  sig <- as.numeric(s$sig)
  set.seed(as.numeric(s$seed))
  out <- list()
  for (N in c(2, 3, 5)) {
    MU <- tr3_centers(N, d, cls_dis)
    n_cls <- tr3_sizes(N, n, as.numeric(s$n1_adj_3cls))
    data_list <- if (s$dist == "uniform") {
      lapply(1:iteN, uniform_cls, n_cls = n_cls, d = d, centers = MU, noise = F,
             rad_range = c(r_min, r_max))
    } else if (s$dist == "gaussian") {
      lapply(1:iteN, Gaussian_cls, n_cls = n_cls, d = d, centers = MU, sigma = sig, noise = F,
             rad_range = c(r_min, r_max))
    } else stop("tr3_sim_data(): unknown dist ", s$dist)
    out[[as.character(N)]] <- list(N = N, n_cls = n_cls,
                                   label_T = rep(seq_len(N) - 1, n_cls), data = data_list)
  }
  out
}

#' Fingerprint of one data set, written with each replicate so that pairing
#' across methods can be checked.
tr3_fp <- function(x) sum(x * seq_along(x))

# ---------------------------------------------------------------------------
# 5. Real data
# ---------------------------------------------------------------------------
#
# data/clustering/Real_Data_Collection.R is used unedited on disk. In memory:
#   - its setwd() points at data/clustering/ (the folder it names was moved
#     on 2026-05-10; CLUSTERING_REPRODUCTION_GUIDE.md section 6.1);
#   - library(readxl) is dropped, and the user_knowledge block, the only one
#     that needs it, runs only if readxl is installed and the set is asked for;
#   - only the requested blocks run. Each block reads its own file and
#     reassigns X and Y before using them, so running a subset gives the same
#     objects as running them all.
# The working directory is restored afterwards.

TR3_REAL_BLOCKS <- c("iris" = "iris", "seeds" = "seeds", "whole sale" = "wholesale",
                     "breast cancer" = "breast_cancer", "aggregation" = "aggregation",
                     "asymmetric" = "asymmetric", "R15" = "R15", "D31" = "D31",
                     "user_knowledge" = "user_knowledge")

tr3_load_real <- function(names) {
  src <- readLines(here::here("data/clustering/Real_Data_Collection.R"), warn = FALSE)
  hdr <- grep("^##### .* #####\\s*$", src)
  sec <- sub("^##### (.*) #####\\s*$", "\\1", src[hdr])
  if (!all(sec %in% names(TR3_REAL_BLOCKS))) {
    stop("tr3_load_real(): unexpected loader sections: ",
         paste(setdiff(sec, names(TR3_REAL_BLOCKS)), collapse = ", "))
  }
  unknown <- setdiff(names, TR3_REAL_BLOCKS)
  if (length(unknown)) stop("tr3_load_real(): unknown data sets: ", paste(unknown, collapse = ", "))
  pre <- src[seq_len(hdr[1] - 1)]
  pre <- pre[!grepl("^\\s*library\\(readxl\\)", pre)]
  sw <- grep("^\\s*setwd\\(", pre)
  if (length(sw) != 1) stop("tr3_load_real(): expected one setwd() line in the loader")
  pre[sw] <- 'setwd(here::here("data/clustering"))'

  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  env <- new.env(parent = globalenv())
  eval(parse(text = pre), envir = env)
  ends <- c(hdr[-1] - 1, length(src))
  out <- list()
  for (i in seq_along(hdr)) {
    obj <- TR3_REAL_BLOCKS[[sec[i]]]
    if (!(obj %in% names)) next
    if (obj == "user_knowledge") {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("tr3_load_real(): user_knowledge needs the readxl package, which is not installed")
      }
      tr3_library("readxl")
    }
    eval(parse(text = src[hdr[i]:ends[i]]), envir = env)
    out[[obj]] <- get(obj, envir = env)
  }
  out
}

# ---------------------------------------------------------------------------
# 6. Results
# ---------------------------------------------------------------------------
#
# G: has dropped off the bus under sustained writes, so every cell is written
# as soon as it finishes, never held until the end of a run. A cell is
# written in one append_result() call, and only then is its key added to the
# matching *_done.csv; has_result() on the done file decides what a rerun
# skips. Readers keep only cells whose key is in the done file, and the last
# copy of a row if a cell was written twice.

tr3_path <- function(name) file.path(TR3_RES, name)

#' If a drive drop cut the last line short, start the next append on a new
#' line, so that the cut line stays ragged (and is ignored) instead of being
#' joined to the first new row.
tr3_ensure_newline <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  sz <- file.size(path)
  if (is.na(sz) || sz == 0) return(invisible(NULL))
  con <- file(path, "rb")
  seek(con, sz - 1)
  last <- readBin(con, "raw", 1)
  close(con)
  if (!identical(last, as.raw(10))) cat("\n", file = path, append = TRUE)
  invisible(NULL)
}

tr3_append <- function(path, df) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tr3_ensure_newline(path)
  append_result(path, df)
}

#' Read a result CSV as text. Lines whose field count differs from the
#' header's (a write cut short) are dropped with a warning; they are filtered
#' before parsing, because read.csv() would wrap an over-long line into an
#' extra row.
tr3_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2) return(NULL)
  nf <- vapply(lines, function(l) {
    k <- tryCatch(utils::count.fields(textConnection(l), sep = ",", quote = "\""),
                  error = function(e) NA_integer_)
    if (length(k) == 1) k else NA_integer_
  }, integer(1), USE.NAMES = FALSE)
  bad <- which(is.na(nf) | nf != nf[1])
  if (length(bad)) {
    warning(sprintf("tr3_read(): %s: dropped %d line(s) with the wrong number of fields (%s)",
                    path, length(bad), paste(bad, collapse = ", ")), call. = FALSE)
    lines <- lines[-bad]
  }
  utils::read.csv(text = lines, stringsAsFactors = FALSE, colClasses = "character")
}

#' Rows of `path` that belong to cells listed in `done_path`, one copy per
#' (keys, `row_key`), the last one written.
tr3_read_done <- function(path, done_path, keys, row_key = NULL) {
  df <- tr3_read(path)
  done <- tr3_read(done_path)
  if (is.null(df) || is.null(done)) return(NULL)
  k_df <- do.call(paste, c(df[keys], sep = "\r"))
  k_done <- do.call(paste, c(done[keys], sep = "\r"))
  df <- df[k_df %in% k_done, , drop = FALSE]
  id <- do.call(paste, c(df[c(keys, row_key)], sep = "\r"))
  df[!duplicated(id, fromLast = TRUE), , drop = FALSE]
}

#' The scripts' printed summary line, so that a console log can be compared
#' with the slurm-*.out files by eye.
tr3_print_summary <- function(method, ari, sil, k_hat, N, delta = NULL) {
  msg <- paste(paste0(method, ": The mean ARI is"), mean(ari), ",",
               "The mean average silhouette index is", mean(sil), ",",
               "and, the percentatage of times identifying the correct number of clusters is",
               mean(k_hat == N) * 100, "%")
  if (!is.null(delta)) msg <- paste0(msg, ", best delta: ", delta)
  print(msg)
}
