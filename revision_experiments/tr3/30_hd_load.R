#!/usr/bin/env Rscript
# revision_experiments/tr3/30_hd_load.R
#
# Loader for the six moderate-dimensional real data sets (10 <= d <= 30) that
# the reviewers asked for (experiment 1). Sourced by 31_hd_real.R and by the
# baseline script for the same data, so that every method sees the same matrix.
#
# USAGE
#   source(here::here("revision_experiments/tr3/30_hd_load.R"))   # needs tr3_lib.R first
#   hd <- tr3_load_hd(c("wine", "vehicle"))        # named list of numeric matrices
#   hd <- tr3_load_hd(TR3_HD_DATASETS, verbose = FALSE)
# Run as a script it prints the summary of all six data sets:
#   Rscript revision_experiments/tr3/30_hd_load.R
# and, with the argument `settings`, (re)writes settings_realdata_hd.csv by the
# rule below:
#   Rscript revision_experiments/tr3/30_hd_load.R settings
#
# RETURNS
#   Each matrix has the features in the first d columns and the label
#   (as.numeric(as.factor(class))) in the LAST column, the shape that
#   tr3_load_real() returns.
#
# PREPROCESSING
#   The steps of data/clustering/Real_Data_Collection.R, in its order, the same
#   for every data set (its breast cancer block is the template):
#     1. rows that contain '?' are dropped (all fields tested, as in the block);
#     2. features to numeric;
#     3. constant feature columns are dropped (added step: a constant column
#        makes scale() return NaN and the NA filter would remove every row;
#        tested on the raw values, so no rounding enters);
#     4. each feature is z-scored with scale();
#     5. rows with NA are dropped (X and the label are subset together; the
#        original blocks subset X only, which differs only if a row is dropped);
#     6. label = as.numeric(as.factor(class));
#     7. cbind(X, Y) and duplicated rows removed with dplyr::distinct.
#
# RAW FILES
#   wine          revision_experiments/tr3_data/wine/wine.data            (class first)
#   vehicle       revision_experiments/tr3_data/vehicle/xaa.dat .. xai.dat (class last)
#   segmentation  revision_experiments/tr3_data/segmentation/segmentation.data
#                 and segmentation.test, together (class first, 5 lines to skip)
#   wdbc          revision_experiments/tr3_data/wdbc/wdbc.data  (column 1 = ID, dropped;
#                 column 2 = diagnosis)
#   hepatitis     data/clustering/hepatitis.data      (class first, '?' = missing)
#   lymphography  data/clustering/lymphography.data   (class first)
# The files in tr3_data/ are read as data only.
#
# SETTINGS RULE (authors' decision of 2026-10-25), by the final d
#   UN-CCDs  NN-test-simul_<d>d_99%.RData  if 10 <= d < 20 (alpha = 1%)
#            NN-test-simul_<d>d_999%.RData if d >= 20 (alpha = 0.1%); direction descend
#   RK-CCDs  RK-test-simul_<d>d_999%.RData, quant 0.999
#   KS-CCDs  the delta grid of settings_realdata.csv, verbatim
#   A final d below 10 stops the rule.

TR3_HD_DATASETS <- c("wine", "hepatitis", "lymphography", "wdbc", "vehicle", "segmentation")

# Expected sizes of the raw files (rows, raw fields per row), checked on reading.
TR3_HD_RAW <- list(
  wine         = c(rows = 178,  fields = 14),
  vehicle      = c(rows = 846,  fields = 19),
  segmentation = c(rows = 2310, fields = 20),
  wdbc         = c(rows = 569,  fields = 32),
  hepatitis    = c(rows = 155,  fields = 20),
  lymphography = c(rows = 148,  fields = 19))

# Raw data frame of character fields, class in `class_col`, features in
# `feat_cols`, as read from the file(s).
tr3_hd_raw <- function(name) {
  td <- function(...) here::here("revision_experiments", "tr3_data", ...)
  rd <- function(f, skip = 0) {
    utils::read.csv(f, header = FALSE, skip = skip, colClasses = "character",
                    stringsAsFactors = FALSE, comment.char = "", blank.lines.skip = TRUE)
  }
  switch(name,
    wine = list(df = rd(td("wine", "wine.data")), class_col = 1, feat_cols = 2:14),
    vehicle = {
      files <- td("vehicle", sprintf("xa%s.dat", letters[1:9]))
      rows <- unlist(lapply(files, function(f) trimws(readLines(f, warn = FALSE))))
      rows <- rows[nzchar(rows)]
      fld <- strsplit(rows, "[[:space:]]+")
      if (any(lengths(fld) != 19)) stop("tr3_hd_raw(): vehicle rows without 19 fields")
      df <- as.data.frame(do.call(rbind, fld), stringsAsFactors = FALSE)
      list(df = df, class_col = 19, feat_cols = 1:18)
    },
    segmentation = {
      a <- rd(td("segmentation", "segmentation.data"), skip = 5)
      b <- rd(td("segmentation", "segmentation.test"), skip = 5)
      list(df = rbind(a, b), class_col = 1, feat_cols = 2:20)
    },
    wdbc = list(df = rd(td("wdbc", "wdbc.data")), class_col = 2, feat_cols = 3:32),
    hepatitis = list(df = rd(here::here("data", "clustering", "hepatitis.data")),
                     class_col = 1, feat_cols = 2:20),
    lymphography = list(df = rd(here::here("data", "clustering", "lymphography.data")),
                        class_col = 1, feat_cols = 2:19),
    stop("tr3_load_hd(): unknown data set ", name))
}

# One data set through the preprocessing; returns the matrix and its summary.
tr3_hd_prepare <- function(name) {
  raw <- tr3_hd_raw(name)
  df <- raw$df
  exp <- TR3_HD_RAW[[name]]
  if (nrow(df) != exp[["rows"]] || ncol(df) != exp[["fields"]]) {
    stop(sprintf("tr3_hd_prepare(): %s: read %d x %d, expected %d x %d", name, nrow(df),
                 ncol(df), exp[["rows"]], exp[["fields"]]))
  }
  n_raw <- nrow(df)

  # 1. drop rows that contain '?'
  df <- df[rowSums(df == "?") == 0, , drop = FALSE]
  n_q <- nrow(df)

  # 2. features to numeric
  X <- as.matrix(as.data.frame(lapply(df[, raw$feat_cols, drop = FALSE], as.numeric)))
  if (anyNA(X)) stop("tr3_hd_prepare(): ", name, ": non-numeric feature field")
  cls <- df[[raw$class_col]]
  feat_names <- paste0("V", raw$feat_cols)
  colnames(X) <- feat_names

  # 3. drop constant columns
  const <- apply(X, 2, function(x) length(unique(x)) == 1)
  dropped <- feat_names[const]
  X <- X[, !const, drop = FALSE]

  # 4. z-score
  X <- apply(X, 2, scale)

  # 5. drop rows with NA
  keep <- rowSums(is.na(X)) == 0
  X <- X[keep, , drop = FALSE]
  cls <- cls[keep]

  # 6. label
  Y <- as.numeric(as.factor(cls))

  # 7. bind and remove duplicated rows
  D <- cbind(X, Y)
  D <- as.matrix(dplyr::distinct(as.data.frame(D)))
  dimnames(D) <- list(NULL, NULL)
  if (anyNA(D)) stop("tr3_hd_prepare(): ", name, ": NA in the final matrix")

  list(data = D,
       info = list(dataset = name, raw_rows = n_raw, rows_after_q = n_q,
                   dropped_constant = if (length(dropped)) paste(dropped, collapse = " ") else "none",
                   rows_after_na = nrow(X), n = nrow(D), d = ncol(D) - 1,
                   class_counts = paste(names(table(D[, ncol(D)])), as.vector(table(D[, ncol(D)])),
                                        sep = ":", collapse = " ")))
}

tr3_load_hd <- function(names = TR3_HD_DATASETS, verbose = TRUE) {
  unknown <- setdiff(names, TR3_HD_DATASETS)
  if (length(unknown)) stop("tr3_load_hd(): unknown data sets: ", paste(unknown, collapse = ", "))
  out <- list()
  for (nm in names) {
    p <- tr3_hd_prepare(nm)
    out[[nm]] <- p$data
    if (verbose) {
      i <- p$info
      cat(sprintf(paste0("%-13s raw rows %d | after '?' filter %d | dropped constant columns: %s | ",
                         "after NA filter %d | final n %d, final d %d | duplicates removed %d | classes %s\n"),
                  i$dataset, i$raw_rows, i$rows_after_q, i$dropped_constant, i$rows_after_na,
                  i$n, i$d, i$rows_after_na - i$n, i$class_counts))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# settings_realdata_hd.csv, by the rule in the header
# ---------------------------------------------------------------------------

TR3_HD_SCRIPT <- "revision_experiments/tr3/31_hd_real.R"

tr3_hd_grid <- function() {
  s <- utils::read.csv(file.path(TR3_DIR, "settings_realdata.csv"), stringsAsFactors = FALSE,
                       colClasses = "character")
  g <- unique(s$delta_grid[s$method == "KS-CCDs"])
  if (length(g) != 1 || !nzchar(g)) stop("tr3_hd_grid(): the KS-CCDs rows use more than one grid")
  g
}

tr3_hd_settings_rows <- function(dims) {
  # `dims`: named vector of final d per data set
  if (any(dims < 10)) {
    stop("tr3_hd_settings_rows(): final d below 10 needs a different rule: ",
         paste(names(dims)[dims < 10], dims[dims < 10], sep = "=", collapse = ", "))
  }
  grid <- tr3_hd_grid()
  rows <- list()
  for (ds in names(dims)) {
    d <- dims[[ds]]
    un <- sprintf("NN-test-simul_%dd_%s.RData", d, if (d < 20) "99%" else "999%")
    rk <- sprintf("RK-test-simul_%dd_999%%.RData", d)
    rows[[length(rows) + 1]] <- c("UN-CCDs", ds, un, "", "descend", "")
    rows[[length(rows) + 1]] <- c("RK-CCDs", ds, rk, "0.999", "", "")
    rows[[length(rows) + 1]] <- c("KS-CCDs", ds, "", "", "", grid)
  }
  m <- do.call(rbind, rows)
  m <- cbind(m, TR3_HD_SCRIPT)
  colnames(m) <- c("method", "dataset", "table", "quant", "direction", "delta_grid", "script")
  m
}

# Same quoting as settings_realdata.csv: only fields that contain a comma.
tr3_hd_write_settings <- function(path, hd) {
  dims <- vapply(hd, function(D) ncol(D) - 1, numeric(1))
  m <- tr3_hd_settings_rows(dims)
  q <- function(x) ifelse(grepl(",", x), paste0("\"", x, "\""), x)
  lines <- c(paste(colnames(m), collapse = ","),
             apply(m, 1, function(r) paste(q(r), collapse = ",")))
  writeLines(lines, path, useBytes = TRUE)
  invisible(m)
}

if (sys.nframe() == 0 && identical(basename(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "30_hd_load.R")) {
  source(here::here("revision_experiments/tr3/tr3_lib.R"))
  suppressPackageStartupMessages(library(dplyr))
  args <- commandArgs(trailingOnly = TRUE)
  hd <- tr3_load_hd(TR3_HD_DATASETS)
  if (length(args) >= 1 && args[1] == "settings") {
    out <- file.path(TR3_DIR, "settings_realdata_hd.csv")
    tr3_hd_write_settings(out, hd)
    cat("wrote", out, "\n")
  }
}
