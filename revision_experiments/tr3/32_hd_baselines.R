#!/usr/bin/env Rscript
# revision_experiments/tr3/32_hd_baselines.R
#
# Experiment 1, baseline part: DBSCAN, MST, K-means++, Spectral and Louvain on the six
# moderate-dimensional real data sets wine, hepatitis, lymphography, wdbc, vehicle and
# segmentation (13 <= d <= 30), the counterpart of 13_real_baselines.R for the data of
# 30_hd_load.R. The data come from tr3_load_hd(): numeric matrices, z-scored, the label
# (as.numeric(as.factor(class))) in the last column, no dimnames.
#
# DESIGN: the code that clusters is the original code, as text
#   For each method this driver takes the script of the paper's existing real-data baselines (the
#   SCRIPTS map below, the same files as 13_real_baselines.R, under
#   simulations/clustering/Algo_Compare_Clustering/), evaluates the script's header lines in
#   .GlobalEnv as 13_real_baselines.R does (library calls and settings: DBSCAN `k = 4` and
#   knee_eps(), MST `thresh = 2` and MST_Clustering(), Louvain `sigma = 2`), and then evaluates ONE
#   block of the script as text. The block "##### seeds data #####" (a z-scored data set) is the
#   template; every block of a script is the same text with another data name (TEXT CHECK below).
#   Only these substitutions are made on the template text:
#     (a) the word `seeds` becomes the name of the data set. The driver stops unless exactly 5 lines
#         contain the word (the block marker, `n_col = ncol(seeds)`, `X = seeds[, -n_col]`,
#         `Y = seeds[, n_col]` and the print line) and each of them has its expected text;
#     (b) K-means++ and Spectral only: the line `k=3` becomes `k=<K>`, K = number of distinct labels
#         of the data set (the number of true classes). The driver stops unless exactly one line
#         matches `^k=<digits>`. DBSCAN's `k = 4` (minPts - 1) is in the header and is not touched.
#         MST and Louvain have no k line (the driver stops if one appears).
#   The data matrix is assigned in .GlobalEnv under its name (`wine`, ...), where the loader had put
#   `seeds` for the original script. The block then runs unchanged: set.seed(1); DBSCAN minPts 5
#   with Eps = knee_eps(X, 4); MST_Clustering(X, 2); KMeans_rcpp(X, clusters = k, initializer =
#   "kmeans++", seed = 1); specc(as.matrix(X), centers = k) with the automatic kernel width;
#   Gaussian-kernel graph with sigma = 2 and cluster_louvain() at the default resolution;
#   ARI = adjustedRandIndex(label, Y), Sil = mean(silhouette(label, dist(X))[, 3]), k_hat =
#   n_clusters of the block (DBSCAN: number of clusters without the noise label, as the block has it).
#   After the block the driver checks that the block's X and Y equal the data set's feature columns
#   and label column (identical()), so a wrong substitution cannot go unnoticed.
#   `show` mode prints the exact text that will be evaluated.
#
# DIFFERENCES FROM THE ORIGINAL BLOCKS (all of them)
#   1. Block text: only (a) and (b) above, plus the marker line (a comment), which reads
#      `##### <data set> data #####`. The `seeds` block is the template for every data set; no other
#      block is used for the run.
#   2. Header, as in 13_real_baselines.R: the script's `source(.../Real_Data_Collection.R)` line is
#      replaced, at its position, by library() calls for the loader's packages (R.matlab, foreign,
#      dplyr, readxl, in the loader's order) and by the data load: tr3_load_hd() in run mode,
#      tr3_load_real() in check mode. The loader's other effects (its setwd(), the function scale_R,
#      X, Y and all nine data objects in .GlobalEnv) do not occur: only the data set of the current
#      cell is assigned, right before its block. MST: the header's cores / makeCluster /
#      registerDoParallel lines are dropped, and the script's closing `stopCluster(cl)` (after its
#      last block) is not run: no block uses the cluster, so nothing in this driver is parallel and
#      TR3_CORES has no effect. DBSCAN: kNNdistplot() and abline() draw on a null pdf device.
#      Header lines that end inside an expression (knee_eps()) are joined with the next ones until
#      they parse; any other parse error stops the run. Package versions are the installed ones
#      (printed at the start of a run).
#   3. Driver variables: everything of this driver but the two source() calls is inside local(), as
#      in 13_real_baselines.R, so no driver variable can coincide with a variable of a header line or
#      of a block (blocks and header lines are evaluated in .GlobalEnv). Checked on 2026-10-25 by
#      parsing the five scripts: none of the names they assign is an object that tr3_lib.R or
#      30_hd_load.R leave in .GlobalEnv. Before each cell the driver removes ari, sil, n_clusters, eps,
#      sigma_used, X, Y and n_col from .GlobalEnv (13_real_baselines.R removes the first five), so a
#      failed cell leaves nothing for the next; every block assigns n_col, X and Y before using them.
#   4. Data form: tr3_load_hd() returns a plain numeric matrix without dimnames; tr3_load_real()
#      returns matrices with dimnames (V1..., Y). Check mode runs the existing sets in both forms
#      ("asloaded" and "plain") to show that the form does not matter.
#   5. Single-cluster silhouette (a guard, not a change of the method): cluster::silhouette() returns
#      NA for a partition with one cluster or with n clusters, and the block's last line,
#      `mean(silhouette(...)[, 3])`, then stops with "incorrect number of dimensions" (checked in R
#      4.6.1 with cluster 2.1.8.2). If, and only if, the block stops with that error after `ari` and
#      `n_clusters` were computed and the block's label vector has 1 or n distinct values, the cell is
#      written with sil = "NA" and a NOTE is printed (31_hd_real.R writes NA for one cluster too).
#      Any other error stops that cell: it is not marked done, the driver goes on with the next cell
#      and exits with status 1 at the end.
#   6. Warnings of a block are printed with the cell (with their counts) and muffled; this changes
#      no value. R would have printed them at the end of the run.
#   7. `seed` column: "NA", the convention of 13_real_baselines.R in verbatim mode; the block's own
#      set.seed(1) (and `seed = 1` of KMeans_rcpp) is the seed. The `script` column is the original
#      script's path.
#   8. `all` runs each method in its own Rscript process, one after the other: the original scripts
#      never had two methods loaded at once, and attaching dbscan, kernlab, ClusterR and igraph in one
#      session would change which function masks which.
#
# USAGE (from the clone root, PowerShell; Rscript segfaults under Git Bash)
#   Rscript revision_experiments/tr3/32_hd_baselines.R <method|all> [datasets] [mode]
#     method    DBSCAN | MST | K-means++ | Spectral | Louvain | all
#     datasets  comma list (use "" to give a mode without datasets). Run and show: from wine, hepatitis,
#               lymphography, wdbc, vehicle, segmentation (default: all six, in that order). Check:
#               from iris, seeds, user_knowledge, wholesale, breast_cancer, aggregation, asymmetric, R15,
#               D31 (default iris, seeds, user_knowledge, breast_cancer)
#     mode      "run" (default) | "check" | "check:<reference file in results/tr3>" | "show"
#   Environment variables
#     TR3_OUT_DIR   output folder (default revision_experiments/results/tr3)
#     TR3_CORES     not used (nothing here is parallel; see difference 2)
#   Do not edit this file while a process runs it (Rscript reads the file as it goes, and `all`
#   starts a new process per method).
#
# OUTPUT (OUT_DIR), run mode; a cell is appended as soon as it is finished
#   real_hd_baselines.csv         method, dataset, n, d, ari, sil, k_hat, seed, script
#                                 (numbers with fmt17; the columns of real_baselines.csv)
#   real_hd_baselines_params.csv  method, dataset, param, value, seed: DBSCAN `eps`; Spectral `sigma`
#                                 (specc's inverse kernel width, k = exp(-sigma ||x-y||^2)) and
#                                 `bandwidth` = 1/sqrt(2 sigma), as 13_real_baselines.R writes them
#   real_hd_baselines_done.csv    method, dataset, finished: written last; a rerun skips these cells
#   The exit status is 1 if a cell failed (it is not marked done), else 0. Log lines: one per data set
#   (n, d, K, fingerprint tr3_fp() of the matrix) and one per cell (ARI, Sil, k_hat, seconds elapsed
#   and cpu).
#
# CHECK (proof that this code path is the original one; no result file is written but the check CSV)
#   Part 1, TEXT CHECK. For each of the nine existing data sets the template with (a) and (b) applied
#   (b with the original block's own k) is compared line by line with that data set's original block
#   (marker line and trailing blank lines left out). A difference in comment lines only is reported
#   and accepted; any other difference fails the check. Result on 2026-10-25: DBSCAN, K-means++,
#   Spectral and Louvain: all nine blocks identical. MST: four blocks identical; iris, aggregation,
#   asymmetric, R15 and D31 differ in one comment line ("# construct MST for clustering with the best
#   cutoff value" against "# construct MST for clustering"). For K-means++ and Spectral the block's k
#   equals the number of distinct labels in all nine data sets.
#   Part 2, numbers. The same code path (header, template, substitutions, K rule) runs on the existing
#   data sets from tr3_load_real(), each in both forms (difference 4), and ari, sil and k_hat (as
#   fmt17 strings; DBSCAN eps and Spectral sigma and bandwidth too, where a params file exists) are
#   compared with the FIRST row per (method, dataset) of each reference file in results/tr3:
#     real_baselines_fixed.csv (+ real_baselines_fixed_params.csv)   decides the exit status
#     real_baselines.csv                                             information only
#   Why two files: real_baselines.csv is the output of 13_real_baselines.R of 2026-10-19, BEFORE the
#   baseline fixes of 2026-10-23 (README.md, "Baseline fixes"); it is byte-identical to the copy in
#   _before_baseline_fix_20261023/. The scripts changed afterwards (DBSCAN: Eps from knee_eps() instead
#   of a typed-in value; Spectral: specc on the data instead of the distance matrix; Louvain and
#   Spectral: set.seed(1) per block; K-means++ on wholesale: k = 3), so the current code cannot equal
#   that file for those cells. real_baselines_fixed.csv holds the rows of the current scripts. The
#   output marks each difference against real_baselines.csv with whether the stored value differs
#   from real_baselines_fixed.csv in that cell (then the script changed there). `check:<file>` uses
#   that single file as the reference that decides.
#   Exit status 0 if the text check passes and every comparison with the deciding reference is equal.
#   Rows go to real_hd_baselines_check_<method>.csv in OUT_DIR (overwritten on each check). With
#   `all` the check files of the earlier run are deleted first, and the summary reads this run's.

source(here::here("revision_experiments/tr3/tr3_lib.R"))
source(here::here("revision_experiments/tr3/30_hd_load.R"))

# Everything below lives in local(): the blocks are evaluated in .GlobalEnv, and no variable of this driver
# can then coincide with a variable that a header line or a block assigns (as in 13_real_baselines.R).
local({
  ALL_METHODS <- c("DBSCAN", "MST", "K-means++", "Spectral", "Louvain")
  USAGE <- paste0("usage: Rscript revision_experiments/tr3/32_hd_baselines.R <method|all> [datasets] [mode]\n",
                  "  method:   ", paste(c(ALL_METHODS, "all"), collapse = " | "), "\n",
                  "  datasets: comma list; run and show: from ", paste(TR3_HD_DATASETS, collapse = ", "),
                  " (default all six);\n            check: from iris, seeds, user_knowledge, wholesale, breast_cancer, ",
                  "aggregation, asymmetric, R15, D31 (default iris, seeds, user_knowledge, breast_cancer)\n",
                  "  mode:     run (default) | check | check:<reference file in results/tr3> | show")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1 || !(args[1] %in% c(ALL_METHODS, "all"))) stop(USAGE, call. = FALSE)
  MODE_ARG <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "run"
  MODE_P   <- strsplit(MODE_ARG, ":", fixed = TRUE)[[1]]
  if (!(MODE_P[1] %in% c("run", "check", "show"))) stop(USAGE, call. = FALSE)
  CHECK    <- MODE_P[1] == "check"
  SHOW     <- MODE_P[1] == "show"
  # Reference files of check mode (results/tr3). The first one decides the exit status.
  CHECK_REFS <- if (CHECK && length(MODE_P) >= 2) MODE_P[2] else c("real_baselines_fixed.csv", "real_baselines.csv")
  CHECK_DS_ALL <- c("iris", "seeds", "user_knowledge", "wholesale", "breast_cancer", "aggregation",
                    "asymmetric", "R15", "D31")
  DS_DEFAULT <- if (CHECK) c("iris", "seeds", "user_knowledge", "breast_cancer") else TR3_HD_DATASETS
  DATASETS <- if (length(args) >= 2 && nzchar(args[2]) && args[2] != "all") strsplit(args[2], ",")[[1]] else DS_DEFAULT
  if (!all(DATASETS %in% if (CHECK) CHECK_DS_ALL else TR3_HD_DATASETS)) stop(USAGE, call. = FALSE)

  OUT_DIR <- sub("[/\\\\]+$", "", Sys.getenv("TR3_OUT_DIR", TR3_RES))
  OUT_CSV    <- file.path(OUT_DIR, "real_hd_baselines.csv")
  PARAMS_CSV <- file.path(OUT_DIR, "real_hd_baselines_params.csv")
  DONE_CSV   <- file.path(OUT_DIR, "real_hd_baselines_done.csv")

  # ---- all: one Rscript process per method, in sequence ----------------------------------
  if (args[1] == "all") {
    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    me <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
    q <- function(x) paste0("\"", x, "\"")   # for the arguments only: system2() quotes the command itself
    chk_files <- file.path(OUT_DIR, sprintf("real_hd_baselines_check_%s.csv", ALL_METHODS))
    # a check file read below must come from this run, not from an earlier one
    if (CHECK) for (f in chk_files[file.exists(chk_files)]) file.remove(f)
    status <- 0L
    for (m in ALL_METHODS) {
      cat(sprintf("32_hd_baselines.R all: starting %s\n", m))
      st <- system2(rscript, c(q(me), q(m), q(paste(DATASETS, collapse = ",")), q(MODE_ARG)))
      cat(sprintf("32_hd_baselines.R all: %s exit status %d\n", m, st))
      if (st != 0L) status <- 1L
    }
    if (CHECK) {
      parts <- lapply(chk_files, tr3_read)
      none <- vapply(parts, is.null, logical(1))
      if (any(none)) {
        cat(sprintf("check summary: NO check file for %s (that process did not finish)\n",
                    paste(ALL_METHODS[none], collapse = ", ")))
        status <- 1L
      }
      allc <- do.call(rbind, parts[!none])
      for (r in unique(allc$reference)) {
        z <- allc[allc$reference == r, , drop = FALSE]
        cat(sprintf("check summary, %d method(s) vs %s (%s): %d comparisons, %d equal, %d different%s\n",
                    length(unique(z$method)), r, if (r == CHECK_REFS[1]) "decides" else "information",
                    nrow(z), sum(z$equal == "TRUE"), sum(z$equal != "TRUE"),
                    if (r != CHECK_REFS[1] && any(z$equal != "TRUE"))
                      sprintf("; %d of them where the stored value differs from %s",
                              sum(z$equal != "TRUE" & z$stored_equals_first_ref != "TRUE"), CHECK_REFS[1]) else ""))
      }
    }
    quit(save = "no", status = status)
  }

  METHOD <- args[1]
  BASE <- "simulations/clustering/Algo_Compare_Clustering"
  SCRIPTS <- c("DBSCAN"    = file.path(BASE, "DBSCAN/Simulation/Real_Data/Real_Data_DBSCAN.R"),
               "MST"       = file.path(BASE, "MST/Simulation/Real_Data/Real_Data_MST_fix_cutoff.R"),
               "K-means++" = file.path(BASE, "K-means++/Simulation/Real_Data/Real_Data_kmeanpp.R"),
               "Spectral"  = file.path(BASE, "Spectral_Clustering/Simulation/Real_Data/Real_Data_Spectrum_Clustering.R"),
               "Louvain"   = file.path(BASE, "Louvain/Simulation/Real_Data/Real_Data_Louvain.R"))
  SCRIPT <- SCRIPTS[[METHOD]]
  USES_K <- METHOD %in% c("K-means++", "Spectral")
  TEMPLATE <- "seeds"

  # ---- header and template (as 13_real_baselines.R; its comments explain each step) --------
  src <- readLines(here::here(SCRIPT), warn = FALSE)
  hdr <- grep("^##### \\S+( data)? #####\\s*$", src)
  blk <- sub("^##### (\\S+)( data)? #####\\s*$", "\\1", src[hdr])
  if (!(TEMPLATE %in% blk)) stop("32_hd_baselines.R: no '", TEMPLATE, "' block in ", SCRIPT)
  tail_line <- grep("^stopCluster\\(cl\\)", src)
  tail_line <- tail_line[tail_line > max(hdr)][1]
  ends <- c(hdr[-1] - 1, if (is.na(tail_line)) length(src) else tail_line - 1)
  TEMPLATE_LINES <- src[hdr[match(TEMPLATE, blk)]:ends[match(TEMPLATE, blk)]]
  if (!grepl("^set\\.seed\\(1\\)\\s*$", TEMPLATE_LINES[2])) {
    stop("32_hd_baselines.R: the line after the block marker is not set.seed(1): ", TEMPLATE_LINES[2])
  }

  # k of the original block of data set `ds` (K-means++, Spectral), NA if none; information only
  orig_k <- function(ds) {
    i <- match(ds, blk)
    if (is.na(i) || !USES_K) return(NA_integer_)
    l <- grep("^k=\\d+", src[hdr[i]:ends[i]], value = TRUE)
    if (length(l) != 1) return(NA_integer_)
    as.integer(sub("^k=(\\d+).*$", "\\1", l))
  }

  # The block text for data set `ds` with K distinct labels: substitutions (a) and (b).
  make_block <- function(ds, K) {
    text <- TEMPLATE_LINES
    has <- grepl(sprintf("\\b%s\\b", TEMPLATE), text)
    if (sum(has) != 5) stop("32_hd_baselines.R: expected 5 template lines with the word '", TEMPLATE,
                            "', found ", sum(has))
    want <- c("^##### seeds data #####", "^n_col = ncol\\(seeds\\)$", "^X = seeds\\[, -n_col\\]",
              "^Y = seeds\\[, n_col\\]", "^print\\(paste\\(\"seeds x ")
    for (w in want) if (sum(grepl(w, text)) != 1) stop("32_hd_baselines.R: template line not found exactly once: ", w)
    text <- gsub(sprintf("\\b%s\\b", TEMPLATE), ds, text)
    kl <- grepl("^k\\s*=\\s*\\d+", text)
    if (USES_K) {
      if (sum(kl) != 1 || !grepl("^k=\\d+", text[kl])) stop("32_hd_baselines.R: expected exactly one 'k=<digits>' line in the template")
      text[kl] <- sub("^k=\\d+", sprintf("k=%d", K), text[kl])
    } else if (any(kl)) {
      stop("32_hd_baselines.R: this method's template has a k line, but the driver expects none")
    }
    text
  }

  # ---- show mode: print the block text that run mode evaluates, then stop ------------------
  if (SHOW) {
    shown <- tr3_load_hd(DATASETS, verbose = FALSE)
    for (ds in DATASETS) {
      K <- length(unique(shown[[ds]][, ncol(shown[[ds]])]))
      cat(sprintf("#### %s on %s (K = %d): the header of %s, then this block ####\n", METHOD, ds, K, SCRIPT))
      cat(make_block(ds, K), sep = "\n")
    }
    quit(save = "no", status = 0)
  }

  LOADED <- local({
    header <- src[seq_len(hdr[1] - 1)]
    is_loader <- grepl("Real_Data_Collection\\.R", header)
    stopifnot(sum(is_loader) == 1)
    is_par <- grepl("^\\s*(cores\\s*=|cl\\s*=\\s*makeCluster|registerDoParallel)", header)
    cat(sprintf("32_hd_baselines.R: %s, mode %s, out=%s\n  script %s\n  header lines: %d; loader line %d replaced; %d cluster line(s) dropped\n",
                METHOD, MODE_ARG, OUT_DIR, SCRIPT, length(header), which(is_loader), sum(is_par)))
    pdf(NULL)   # kNNdistplot() / abline() in the DBSCAN block
    pending <- character(0)
    loaded <- NULL
    for (i in seq_along(header)) {
      if (is_loader[i]) {
        stopifnot(length(pending) == 0)
        for (p in c("R.matlab", "foreign", "dplyr", "readxl")) suppressPackageStartupMessages(library(p, character.only = TRUE))
        # the data are loaded where the script sourced Real_Data_Collection.R (with the same packages attached):
        # run mode: the six new sets from tr3_load_hd(); check mode: all nine existing sets from tr3_load_real()
        loaded <- if (CHECK) tr3_load_real(CHECK_DS_ALL) else tr3_load_hd(DATASETS)
      } else if (!is_par[i]) {
        pending <- c(pending, header[i])
        ex <- tryCatch(parse(text = pending), error = function(e) e)
        if (inherits(ex, "error")) {
          if (!grepl("unexpected end of input", conditionMessage(ex), fixed = TRUE)) {
            stop("32_hd_baselines.R: header line ", i, " does not parse: ", conditionMessage(ex))
          }
        } else {
          eval(ex, envir = globalenv())
          pending <- character(0)
        }
      }
    }
    if (length(pending)) stop("32_hd_baselines.R: the header ends inside an unfinished expression: ", pending[1])
    cat("  versions: R", as.character(getRversion()), "| igraph", as.character(packageVersion("igraph")),
        "| mclust", as.character(packageVersion("mclust")), "| cluster", as.character(packageVersion("cluster")),
        if (METHOD == "DBSCAN") paste("| dbscan", packageVersion("dbscan")) else "",
        if (METHOD == "K-means++") paste("| ClusterR", packageVersion("ClusterR")) else "",
        if (METHOD == "Spectral") paste("| kernlab", packageVersion("kernlab")) else "", "\n")
    loaded
  })

  # ---- one cell: the block for `ds` on the data `D` (features, label in the last column) ----
  # Returns a list: row (the result row), params (data frame or NULL), seconds, cpu, note, warns, K, line.
  # The block is evaluated in .GlobalEnv, as the scripts ran.

  # The block's own variable that holds the partition (needed only by the single-cluster guard below).
  LABEL_VAR <- c("DBSCAN" = "db_label", "MST" = "MST_label", "K-means++" = "kmeanpp_label",
                 "Spectral" = "spectrum_label", "Louvain" = "louvain_label")

  run_cell <- function(ds, D) {
    G <- globalenv()
    Ycol <- ncol(D)
    K <- length(unique(D[, Ycol]))
    text <- make_block(ds, K)
    for (v in c("ari", "sil", "n_clusters", "eps", "sigma_used", "X", "Y", "n_col")) {
      if (exists(v, envir = G, inherits = FALSE)) rm(list = v, envir = G)
    }
    assign(ds, D, envir = G)
    t0 <- proc.time()
    note <- ""
    warns <- character(0)   # warnings of the block are logged with the cell (and muffled; they change nothing)
    err <- tryCatch(
      withCallingHandlers({ eval(parse(text = text), envir = G); NULL },
                          warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    if (!is.null(err)) {
      # The only error that is turned into a value: the block's last line (the silhouette), when the partition has
      # one cluster or n clusters. There cluster::silhouette() returns NA and `NA[, 3]` stops with "incorrect
      # number of dimensions" (checked: R 4.6.1, cluster 2.1.8.2). 31_hd_real.R writes NA for this case too.
      lab <- if (exists(LABEL_VAR[[METHOD]], envir = G, inherits = FALSE)) get(LABEL_VAR[[METHOD]], envir = G) else NULL
      n_lab <- length(unique(lab))
      in_sil <- exists("ari", envir = G, inherits = FALSE) && exists("n_clusters", envir = G, inherits = FALSE) &&
        !exists("sil", envir = G, inherits = FALSE) &&
        grepl("incorrect number of dimensions", conditionMessage(err), fixed = TRUE)
      if (!(in_sil && !is.null(lab) && (n_lab <= 1 || n_lab >= length(lab)))) stop(err)
      assign("sil", NA_real_, envir = G)
      note <- sprintf("NOTE: %s %s: silhouette undefined, the partition has %d distinct label(s) for %d points (n_clusters %d); sil written as NA",
                      METHOD, ds, n_lab, length(lab), get("n_clusters", envir = G))
    }
    tm <- proc.time() - t0
    secs <- tm[["elapsed"]]
    cpu <- tm[["user.self"]] + tm[["sys.self"]]
    # the block must have used exactly this data set
    if (!identical(unname(as.matrix(get("X", envir = G))), unname(as.matrix(D[, -Ycol, drop = FALSE]))) ||
        !identical(as.vector(get("Y", envir = G)), as.vector(D[, Ycol]))) {
      stop("32_hd_baselines.R: the block's X / Y differ from the data set ", ds)
    }
    ari <- get("ari", envir = G); sil <- get("sil", envir = G); n_clusters <- get("n_clusters", envir = G)
    X <- get("X", envir = G)
    row <- data.frame(method = METHOD, dataset = ds, n = nrow(X), d = ncol(X),
                      ari = fmt17(ari), sil = fmt17(sil), k_hat = n_clusters,
                      seed = "NA", script = SCRIPT, stringsAsFactors = FALSE)
    par_rows <- if (METHOD == "DBSCAN") {
      data.frame(param = "eps", value = fmt17(get("eps", envir = G)))
    } else if (METHOD == "Spectral") {
      sg <- get("sigma_used", envir = G)
      data.frame(param = c("sigma", "bandwidth"), value = fmt17(c(sg, 1 / sqrt(2 * sg))))
    } else NULL
    params <- if (is.null(par_rows)) NULL else
      data.frame(method = METHOD, dataset = ds, par_rows, seed = "NA", stringsAsFactors = FALSE)
    list(row = row, params = params, seconds = secs, cpu = cpu, note = note, warns = warns, K = K,
         line = sprintf("  %s %s (K %d): ari %.6f sil %.6f k_hat %d (%.1f s elapsed, %.1f s cpu)",
                        METHOD, ds, K, ari, sil, n_clusters, secs, cpu))
  }

  # the cell's note, its warnings (with counts) and its result line, to the log
  log_cell <- function(res, tag = "") {
    if (nzchar(res$note)) cat("  ", res$note, "\n", sep = "")
    if (length(res$warns)) {
      w <- table(res$warns)
      for (i in seq_along(w)) cat(sprintf("  WARNING (%dx): %s\n", w[[i]], names(w)[i]))
    }
    cat(res$line, tag, "\n", sep = "")
  }

  # ============================== run mode ==============================================
  if (!CHECK) {
    for (ds in DATASETS) {   # the fingerprint tr3_fp() allows comparing the matrix with the one 31_hd_real.R used
      D <- LOADED[[ds]]
      cat(sprintf("  data %-13s n %4d, d %2d, K %d, fingerprint %s\n", ds, nrow(D), ncol(D) - 1,
                  length(unique(D[, ncol(D)])), fmt17(tr3_fp(D))))
    }
    failed <- character(0)
    for (ds in DATASETS) {
      key <- list(method = METHOD, dataset = ds)
      if (has_result(DONE_CSV, key)) {
        cat(sprintf("  %s %s: done, skipped\n", METHOD, ds))
        next
      }
      res <- tryCatch(run_cell(ds, LOADED[[ds]]), error = function(e) e)
      if (inherits(res, "error")) {
        cat(sprintf("  %s %s: FAILED: %s\n", METHOD, ds, conditionMessage(res)), file = stderr())
        failed <- c(failed, ds)
        next
      }
      tr3_append(OUT_CSV, res$row)
      if (!is.null(res$params)) tr3_append(PARAMS_CSV, res$params)
      tr3_append(DONE_CSV, data.frame(key, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE))
      log_cell(res)
    }
    if (length(failed)) {
      cat(sprintf("32_hd_baselines.R: %s: %d cell(s) failed: %s\n", METHOD, length(failed), paste(failed, collapse = ", ")))
      quit(save = "no", status = 1)
    }
    quit(save = "no", status = 0)
  }

  # ============================== check mode ============================================

  # ---- 1. block text: what the two substitutions are worth -------------------------------------
  # For each of the nine existing data sets, the template block with the data name (and, for K-means++ and
  # Spectral, the block's own k) substituted must equal that data set's original block. The marker line is a
  # comment line and is left out; trailing blank lines do not count.
  trim_blank <- function(x) { while (length(x) && !nzchar(trimws(x[length(x)]))) x <- x[-length(x)]; x }
  text_ok <- TRUE
  cat("  block text, template with substitutions against the original block of each data set:\n")
  for (ds in CHECK_DS_ALL) {
    i <- match(ds, blk)
    if (is.na(i)) stop("32_hd_baselines.R check: no block for ", ds, " in ", SCRIPT)
    Dd <- LOADED[[ds]]
    K <- length(unique(Dd[, ncol(Dd)]))
    kb <- orig_k(ds)
    if (USES_K && is.na(kb)) stop("32_hd_baselines.R check: no 'k=<digits>' line in the original block of ", ds)
    gen <- trim_blank(make_block(ds, if (USES_K) kb else K)[-1])
    org <- trim_blank(src[hdr[i]:ends[i]][-1])
    dif <- if (length(gen) == length(org)) which(gen != org) else seq_len(max(length(gen), length(org)))
    comments_only <- length(dif) > 0 && length(gen) == length(org) &&
      all(grepl("^\\s*#", gen[dif]) & grepl("^\\s*#", org[dif]))
    verdict <- if (length(dif) == 0) "identical" else if (comments_only) sprintf("identical except %d comment line(s)", length(dif)) else "DIFFERENT"
    if (verdict == "DIFFERENT") text_ok <- FALSE
    cat(sprintf("    %-15s %-34s%s\n", ds, verdict,
                if (USES_K) sprintf("k in the block %2d, number of labels %2d%s", kb, K, if (kb == K) "" else "  <-- differ") else ""))
    if (comments_only) for (j in dif) cat(sprintf("      comment line, original: %s\n      comment line, template: %s\n", org[j], gen[j]))
  }

  # ---- 2. numbers: this code path on the existing data sets against the stored rows ------------
  read_ref <- function(f) {
    rows <- tr3_read(file.path(TR3_RES, f))
    if (is.null(rows)) stop("32_hd_baselines.R check: cannot read the reference ", file.path(TR3_RES, f))
    rows <- rows[rows$method == METHOD, , drop = FALSE]
    list(file = f, rows = rows[!duplicated(rows$dataset), , drop = FALSE],   # first row per (method, dataset)
         params = tr3_read(file.path(TR3_RES, sub("\\.csv$", "_params.csv", f))))
  }
  REFS <- lapply(CHECK_REFS, read_ref)
  for (ri in seq_along(REFS)) {
    cat(sprintf("  reference %s (%s): %d cell(s) of %s, first row of each; %s\n", REFS[[ri]]$file,
                if (ri == 1) "decides the exit status" else "information only", nrow(REFS[[ri]]$rows), METHOD,
                if (is.null(REFS[[ri]]$params)) "no params file" else
                  paste("params:", sub("\\.csv$", "_params.csv", REFS[[ri]]$file))))
  }
  check_rows <- list()
  add <- function(form, ds, ri, quantity, stored, got, first) {
    stored <- ifelse(is.na(stored), "NA", stored)   # tr3_read() reads the text "NA" as NA; fmt17() writes "NA"
    first  <- ifelse(is.na(first), "NA", first)
    check_rows[[length(check_rows) + 1]] <<- data.frame(method = METHOD, dataset = ds, form = form,
      reference = REFS[[ri]]$file, quantity = quantity, stored = stored, reproduced = got,
      equal = identical(stored, got), stored_equals_first_ref = identical(stored, first),
      stringsAsFactors = FALSE)
  }
  for (ds in DATASETS) {
    D0 <- LOADED[[ds]]
    cat(sprintf("  data set %s: class %s, dimnames %s; number of labels %d\n", ds, paste(class(D0), collapse = "/"),
                if (is.null(dimnames(D0))) "none" else "yes", length(unique(D0[, ncol(D0)]))))
    for (form in c("asloaded", "plain")) {
      D <- if (form == "plain") { m <- as.matrix(D0); dimnames(m) <- NULL; m } else D0
      res <- run_cell(ds, D)
      log_cell(res, sprintf(" [%s]", form))
      for (ri in seq_along(REFS)) {
        r  <- REFS[[ri]]$rows[REFS[[ri]]$rows$dataset == ds, , drop = FALSE]
        r1 <- REFS[[1]]$rows[REFS[[1]]$rows$dataset == ds, , drop = FALSE]
        if (nrow(r) != 1 || nrow(r1) != 1) stop("32_hd_baselines.R check: no reference row for ", METHOD, " ", ds)
        for (q in c("ari", "sil", "k_hat")) {
          add(form, ds, ri, q, as.character(r[[q]]), as.character(res$row[[q]]), as.character(r1[[q]]))
        }
        P <- REFS[[ri]]$params; P1 <- REFS[[1]]$params
        if (!is.null(res$params) && !is.null(P)) {
          for (p in res$params$param) {
            sp  <- P$value[P$method == METHOD & P$dataset == ds & P$param == p]
            sp1 <- if (is.null(P1)) NA_character_ else P1$value[P1$method == METHOD & P1$dataset == ds & P1$param == p][1]
            if (length(sp) >= 1) add(form, ds, ri, p, sp[1], res$params$value[res$params$param == p], sp1)
          }
        }
      }
    }
  }
  chk <- do.call(rbind, check_rows)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  out_chk <- file.path(OUT_DIR, sprintf("real_hd_baselines_check_%s.csv", METHOD))
  if (file.exists(out_chk)) file.remove(out_chk)
  for (i in seq_len(nrow(chk))) tr3_append(out_chk, chk[i, ])
  for (i in seq_len(nrow(chk))) {
    cat(sprintf("  %-9s %-14s %-8s %-9s %-24s stored %-22s reproduced %-22s %s\n", chk$method[i], chk$dataset[i],
                chk$form[i], chk$quantity[i], chk$reference[i], chk$stored[i], chk$reproduced[i],
                if (chk$equal[i]) "EQUAL" else if (chk$reference[i] != REFS[[1]]$file && !chk$stored_equals_first_ref[i])
                  sprintf("DIFFERENT (the stored value differs from %s: the script changed)", REFS[[1]]$file) else "DIFFERENT"))
  }
  for (ri in seq_along(REFS)) {
    z <- chk[chk$reference == REFS[[ri]]$file, , drop = FALSE]
    cat(sprintf("32_hd_baselines.R check: %s vs %s (%s): %d comparisons, %d equal, %d different%s\n", METHOD, REFS[[ri]]$file,
                if (ri == 1) "decides" else "information", nrow(z), sum(z$equal), sum(!z$equal),
                if (ri > 1 && any(!z$equal))
                  sprintf("; %d of them where the stored value differs from %s", sum(!z$equal & !z$stored_equals_first_ref), REFS[[1]]$file) else ""))
  }
  cat(sprintf("32_hd_baselines.R check: %s block text: %s\n", METHOD,
              if (text_ok) "all nine blocks reproduced (identical, or only comment lines differ)" else "DIFFERENT"))
  ok <- text_ok && all(chk$equal[chk$reference == REFS[[1]]$file])
  cat(sprintf("32_hd_baselines.R check: %s: %s\n", METHOD, if (ok) "PASS" else "FAIL"))
  quit(save = "no", status = if (ok) 0 else 1)
})
