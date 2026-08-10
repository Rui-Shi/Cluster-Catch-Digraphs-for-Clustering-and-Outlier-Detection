#!/usr/bin/env Rscript
# 51_wp2a_smin_bigfour.R -- finishes WP2(a)'s S_min sensitivity sweep.
#
# 44_wp2a_smin_sweep16.R's runner died after completing 227 of the intended
# 320 cells (16 datasets x 2 methods x 10-point grid). What is done (12
# datasets fully, PenDigits partially) sits, untouched, in
# results/tr1/wp2a_smin_sweep16.csv, which this script only ever READS -- it
# is never appended to or rewritten (see CLAUDE.md hard rule and the task
# brief: 227 completed cells must not be disturbed). This script writes a
# SEPARATE, NEW file, results/tr1/wp2a_smin_bigfour.csv, using the identical
# column schema (read directly off wp2a_smin_sweep16.csv's own header, not
# retyped by hand, so the two are guaranteed to concatenate cleanly later).
#
# Only SU-MCCD and SUN-MCCD take min.cls (a PROPORTION of n -- see
# R/ccds/RK_CCD_New.R / R/ccds/UN_CCD.R; passing a raw count instead of a
# proportion collapses lenD to 0, see wp0_mccd_methods.R's own "UNITS FIX"
# comment). U-MCCD / UN-MCCD are out of scope, exactly as in 44.
#
# WHAT IS MISSING (per the task brief, cross-checked against the actual CSV
# on 2026-08-09 -- see the Rscript one-liners in the session log; PenDigits'
# specific missing cells are NOT hardcoded here, they are computed at
# runtime by reading wp2a_smin_sweep16.csv, per instruction):
#
#   pageblocks, thyroid, waveform, wilt -- NOTHING done. These four are run
#     on a DELIBERATELY REDUCED grid, S_min in {0.01, 0.0625, 0.15, 0.30}
#     (4 points, not the full 10), because the full 10-point grid is not
#     affordable at these n (4795-4819 rows; PenDigits alone, at n=3200,
#     already ran 190-355s/cell in the base sweep -- see wp2a_smin_sweep16.csv
#     elapsed_sec for PenDigits). Every row for these four datasets carries a
#     `note` recording the reduction, so partial coverage here is never later
#     mistaken for the full 10-point grid the other 12 datasets have.
#
#   PenDigits -- partially done (confirmed by reading the CSV at runtime, not
#     assumed): SU-MCCD has 4/10 (0.005, 0.02, 0.0625, 0.30 present),
#     SUN-MCCD has 3/10 (0.005, 0.02, 0.0625 present). This script fills in
#     exactly the missing points of the FULL 10-point grid
#     {0.005,0.01,0.02,0.03,0.05,0.0625,0.10,0.15,0.20,0.30} for PenDigits --
#     no reduction, no note needed on success.
#
# CHECKPOINTING. A cell (dataset, method, s_min) is skipped if it already
# appears in EITHER wp2a_smin_sweep16.csv (read-only -- this is how
# PenDigits' already-done cells are recognized and not redone) OR this
# script's own wp2a_smin_bigfour.csv (read-write -- this is the restart
# path). The float-safe `cell_done()` check is copied from
# 44_wp2a_smin_sweep16.R verbatim (harness.R's has_result() does a string
# compare and is unsafe for floats whose printed form drops trailing zeros,
# e.g. "0.10" -> "0.1").
#
# Usage:
#   Rscript 51_wp2a_smin_bigfour.R --smoke
#     Runs ONE cheap, already-known cell (hepatitis x SU-MCCD x s_min=0.0625,
#     ~2s in the original sweep) through the exact same run_cell() path,
#     against the real wp2a_smin_bigfour.csv, so its BA/F2 can be diffed
#     against wp2a_smin_sweep16.csv's existing row for the same key. Run
#     twice in a row to see the second invocation skip via resume. THIS ROW
#     MUST BE DELETED (or the whole file removed, since it would be the only
#     row) before the real run below, since hepatitis is not part of this
#     script's intended scope.
#
#   Rscript 51_wp2a_smin_bigfour.R [datasets] [methods] [smins]
#     The real run. All three positional args are optional:
#       datasets  comma list, or "ALL" (default) = the 5 in-scope datasets
#                 (pageblocks,thyroid,waveform,wilt,PenDigits), or "BIGFOUR"
#                 for just the reduced-grid four, or "PENDIGITS" for just
#                 PenDigits.
#       methods   comma list, default "SU-MCCD,SUN-MCCD"
#       smins     comma list of numeric S_min values, or "AUTO" (default):
#                 for pageblocks/thyroid/waveform/wilt, AUTO resolves to the
#                 reduced grid {0.01,0.0625,0.15,0.30}; for PenDigits, AUTO
#                 resolves to whatever full-grid points
#                 wp2a_smin_sweep16.csv does not yet have for that
#                 (dataset,method) pair, read fresh at runtime. An explicit
#                 list overrides AUTO uniformly for every dataset in the run
#                 (used by --smoke and for manual re-runs of a single point).
#
# Per-cell timeout: 2400s (task spec), recorded as status="timeout" with the
# measured elapsed time -- never silently dropped. CSV rows are
# written/appended via harness.R's append_result() (write.csv/write.table,
# no held-open connection), so a killed process loses at most the cell in
# flight.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

SWEEP16_CSV <- here::here("revision_experiments/results/tr1/wp2a_smin_sweep16.csv")  # READ-ONLY, never written
OUT_CSV     <- here::here("revision_experiments/results/tr1/wp2a_smin_bigfour.csv")  # NEW file, ours

stopifnot(file.exists(SWEEP16_CSV))
# Schema is read off the existing sweep, not retyped, so the two files are
# guaranteed to concatenate cleanly (same column names, same order).
ROW_COLS <- names(read.csv(SWEEP16_CSV, nrows = 0, stringsAsFactors = FALSE))

METHODS_ALL <- c("SU-MCCD", "SUN-MCCD")
stopifnot(all(METHODS_ALL %in% names(METHOD_REGISTRY)))

FULL_GRID    <- c(0.005, 0.01, 0.02, 0.03, 0.05, 0.0625, 0.10, 0.15, 0.20, 0.30)
REDUCED_GRID <- c(0.01, 0.0625, 0.15, 0.30)

BIGFOUR <- c("pageblocks", "thyroid", "waveform", "wilt")
DATASETS_DEFAULT <- c(BIGFOUR, "PenDigits")

REDUCTION_NOTE <- paste(
  "reduced grid: S_min in {0.01, 0.0625, 0.15, 0.30} only (4 of the full 10",
  "points) -- the full 10-point grid is not affordable at this n; see",
  "51_wp2a_smin_bigfour.R header for the measured-cost rationale."
)

CELL_TIMEOUT_SEC <- 2400   # per task spec

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
collapse_num <- function(v) if (length(v) == 0) "" else paste(v, collapse = ";")

#' Float-safe existence check (copied verbatim from 44_wp2a_smin_sweep16.R).
cell_done <- function(csv_path, dataset, method, s_min) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  any(df$dataset == dataset & df$method == method & abs(df$s_min - s_min) < 1e-9)
}

#' A cell counts as done if it is in EITHER the read-only base sweep or our
#' own (resumable) output file. `ignore_base` skips the base-sweep half of
#' the check -- used ONLY by --smoke below, to force a live recomputation of
#' a cell that (by design) already exists in wp2a_smin_sweep16.csv, so its
#' output can be diffed against that existing row. The real run never sets
#' this (default FALSE): it must never redo a cell the base sweep already has.
cell_done_anywhere <- function(dataset, method, s_min, ignore_base = FALSE) {
  (!ignore_base && cell_done(SWEEP16_CSV, dataset, method, s_min)) ||
    cell_done(OUT_CSV, dataset, method, s_min)
}

#' PenDigits-only: which FULL_GRID points does wp2a_smin_sweep16.csv NOT yet
#' have for this method? Read fresh every call -- never hardcoded.
pendigits_missing_smins <- function(method) {
  df <- read.csv(SWEEP16_CSV, stringsAsFactors = FALSE)
  have <- df$s_min[df$dataset == "PenDigits" & df$method == method]
  FULL_GRID[!sapply(FULL_GRID, function(g) any(abs(have - g) < 1e-9))]
}

run_cell <- function(dataset, method, s_min, extra_note = NA_character_, ignore_base = FALSE) {
  if (isTRUE(cell_done_anywhere(dataset, method, s_min, ignore_base = ignore_base))) {
    cat(sprintf("[skip, already recorded] %-12s x %-9s s_min=%.4f\n", dataset, method, s_min))
    return(invisible("skip"))
  }

  dat <- load_real_dataset(dataset)
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n

  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT_SEC, transient = TRUE)
    res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = s_min)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

    if (length(res$score) != n || anyNA(res$score)) {
      stop(sprintf("score sanity check failed: length=%d (n=%d), any NA=%s",
                    length(res$score), n, anyNA(res$score)))
    }
    thr <- REAL_DATA_THRESHOLDS[[method]]
    m <- evaluate(Y, res$score, thr)

    cl <- res$cluster
    cl <- cl[!is.na(cl)]
    tab <- table(cl)
    csizes <- as.integer(tab)

    list(m = m, score = res$score, n_clusters = length(tab),
         median_cluster_size = if (length(csizes)) median(csizes) else NA_real_,
         cluster_sizes = csizes,
         status = "ok", note = NA_character_)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    is_timeout <- grepl("elapsed time limit|reached CPU time limit|reached elapsed", msg)
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
         score = NULL, n_clusters = NA_integer_, median_cluster_size = NA_real_,
         cluster_sizes = integer(0),
         status = if (is_timeout) "timeout" else "error",
         note = substr(msg, 1, 300))
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  flagged <- if (is.null(out$score)) integer(0) else which(out$score == 1)

  note_final <- if (!is.na(out$note) && !is.na(extra_note)) paste0(extra_note, " | ", out$note)
                else if (!is.na(out$note)) out$note
                else extra_note

  row <- as.data.frame(list(
    dataset = dataset, method = method, s_min = s_min, n = n, d = d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
    n_flagged = if (is.null(out$score)) NA_integer_ else length(flagged),
    n_clusters = out$n_clusters, median_cluster_size = out$median_cluster_size,
    cluster_sizes = collapse_num(out$cluster_sizes),
    flagged_idx = collapse_idx(flagged),
    elapsed_sec = wall, status = out$status, note = note_final,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ), stringsAsFactors = FALSE)[ROW_COLS]

  append_result(OUT_CSV, row)

  fields <- sprintf("dataset=%s method=%s param=s_min:%.4f status=%s sec=%.1f",
                     dataset, method, s_min, out$status, wall)
  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE %s BA=%.3f F2=%.3f nflag=%s ncls=%s\n",
                fields,
                if (is.na(out$m[["BA"]])) NA_real_ else out$m[["BA"]],
                if (is.na(out$m[["F2"]])) NA_real_ else out$m[["F2"]],
                if (is.null(out$score)) "NA" else length(flagged),
                ifelse(is.na(out$n_clusters), "NA", out$n_clusters)))
  } else {
    note_1line <- gsub("[\r\n]+", " ", substr(as.character(note_final), 1, 150))
    cat(sprintf("CELL_FAIL %s note=%s\n", fields, note_1line))
  }
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# --smoke mode
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (any(args == "--smoke")) {
  # ignore_base=TRUE: hepatitis/SU-MCCD/0.0625 already exists in the
  # read-only wp2a_smin_sweep16.csv by design (it is one of the 12 fully
  # completed datasets) -- the smoke test needs a LIVE recomputation of that
  # exact cell so its output can be diffed against the existing row, not a
  # skip. The second call below does NOT set ignore_base and is skipped
  # purely because the FIRST call just wrote it to our own OUT_CSV -- that is
  # the resume-logic proof, isolated from the (also-true, but different)
  # fact that the base sweep already has this cell.
  cat("==== 51_wp2a_smin_bigfour --smoke: hepatitis x SU-MCCD x s_min=0.0625 (forced live run) ====\n")
  run_cell("hepatitis", "SU-MCCD", 0.0625, ignore_base = TRUE)
  cat("\n==== 51_wp2a_smin_bigfour --smoke: re-invoking the SAME cell to demonstrate resume/skip (own-file check) ====\n")
  run_cell("hepatitis", "SU-MCCD", 0.0625)
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Real run
# ---------------------------------------------------------------------------
resolve_datasets <- function(s) {
  if (is.null(s) || !nzchar(s) || toupper(s) == "ALL") return(DATASETS_DEFAULT)
  if (toupper(s) == "BIGFOUR")   return(BIGFOUR)
  if (toupper(s) == "PENDIGITS") return("PenDigits")
  strsplit(s, ",")[[1]]
}
resolve_methods <- function(s) if (is.null(s) || !nzchar(s) || toupper(s) == "ALL") METHODS_ALL else strsplit(s, ",")[[1]]

#' Resolve the S_min list for ONE dataset, honouring an explicit override if
#' given, else AUTO per-dataset logic (reduced grid for the big four,
#' missing-from-sweep16 for PenDigits, full grid as a generic fallback).
resolve_smins_for_dataset <- function(dataset, method, s) {
  if (!is.null(s) && nzchar(s) && toupper(s) != "AUTO") return(as.numeric(strsplit(s, ",")[[1]]))
  if (dataset %in% BIGFOUR) return(REDUCED_GRID)
  if (dataset == "PenDigits") return(pendigits_missing_smins(method))
  FULL_GRID
}

DATASETS <- resolve_datasets(if (length(args) >= 1) args[1] else NULL)
METHODS  <- resolve_methods(if (length(args) >= 2) args[2] else NULL)
SMINS_ARG <- if (length(args) >= 3) args[3] else NULL

cat(sprintf("51_wp2a_smin_bigfour: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("51_wp2a_smin_bigfour: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("51_wp2a_smin_bigfour: smins arg = %s ; out = %s\n",
            if (is.null(SMINS_ARG)) "AUTO" else SMINS_ARG, OUT_CSV))

for (ds in DATASETS) {
  is_reduced <- ds %in% BIGFOUR
  for (method in METHODS) {
    smins <- resolve_smins_for_dataset(ds, method, SMINS_ARG)
    cat(sprintf("\n== %s x %s : s_min grid = %s ==\n", ds, method, paste(smins, collapse = ", ")))
    for (s_min in smins) {
      run_cell(ds, method, s_min, extra_note = if (is_reduced) REDUCTION_NOTE else NA_character_)
    }
  }
}

cat("\n51_wp2a_smin_bigfour: run done.\n")
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  cat(sprintf("%s now has %d rows.\n", OUT_CSV, nrow(final)))
  cat(sprintf("status counts: %s\n", paste(names(table(final$status)), table(final$status), sep = "=", collapse = ", ")))
}
cat("ALL_CELLS_COMPLETE\n")
