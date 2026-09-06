#!/usr/bin/env Rscript
# revision_experiments/tr1/84_wp5_highd.R
#
# WP5 (AE.4, R1.8, R3.8, R5.6): the four proposed MCCD detectors plus the
# five original baselines, on the four real data sets above d = 21 (letter
# d=32, mnist d=100, musk d=166, arrhythmia d=274).
#
# Everything this script does is declared in advance in
# revision_experiments/tr1/WP5_PROTOCOL.md. Read that first; this file is the
# implementation, not the specification. Data comes from
# revision_experiments/tr1/84a_wp5_fetch_convert.py's output
# (results/tr1/wp5/data/<dataset>.csv, feature columns V1..Vd + a trailing
# `label` column, 1 = regular / 0 = outlier -- this repo's convention, not
# PyOD's), not from load_real_dataset() (RealData_Collection.R never saw
# these four data sets).
#
# METHOD SETTINGS -- unchanged from every other WP that has used them:
#   * S_min = 0.05, passed as min.cls to SU-MCCD/SUN-MCCD only (U-MCCD/
#     UN-MCCD take no min.cls argument at all).
#   * alpha resolved by shared/harness.R's own three resolvers
#     (rk_quant_label_paper, nn_quant_label_paper_UN,
#     nn_quant_label_paper_SUN), called with no override -- at every d in
#     this WP (32, 100, 166, 274) all three resolve to "999" (see
#     WP5_PROTOCOL.md S4.1).
#   * Baselines (LOF, DBSCAN, MST, ODIN, iForest) at METHOD_REGISTRY's own
#     defaults.
#
# THE N/A RULE (WP5_PROTOCOL.md S5): U-MCCD and SU-MCCD are RK-based. No
# production RK quantile table exists at d=166 (musk) or d=274 (arrhythmia)
# at all -- this is checked by file existence BEFORE calling the detector,
# not by attempt-then-catch, because a missing table there is a known
# structural fact (100% zero quantiles, confirmed by probe tables), not an
# error condition. Those cells are written with status="n/a" and an explicit
# reason string. At d=32 (letter) and d=100 (mnist), RK tables exist and the
# methods run, but the row is flagged with the measured zero-quantile
# fraction (51.276% / 91.4%) so the degeneracy trend is visible across all
# four running/n/a cells together -- see rk_degeneracy.csv, written once per
# invocation regardless of scope.
#
# CHECKPOINTING: per-cell append via the harness's own append_result()/
# has_result(), keyed on (dataset, method) -- resumable, because the volume
# this repo lives on drops off the bus under sustained write load (CLAUDE.md).
#
# USAGE
#   Rscript "revision_experiments/tr1/84_wp5_highd.R" [--datasets=a,b,...] [--methods=a,b,...] [--smoke]
#
#   --datasets=...   comma-separated subset of {letter, mnist, musk,
#                    arrhythmia}. Default: all four.
#   --methods=...    comma-separated subset of {U-MCCD, SU-MCCD, UN-MCCD,
#                    SUN-MCCD, LOF, DBSCAN, MST, ODIN, iForest}. Default: all
#                    nine.
#   --smoke          overrides both of the above to a single cheap cell
#                    (arrhythmia x LOF) and writes to results/tr1/wp5/smoke/
#                    instead of the production output -- proves the driver
#                    wires correctly without attempting the full grid, per
#                    the WP5 task brief ("do NOT run the full configuration").
#
# This script does NOT run the full configuration by default just because it
# is invoked with no arguments -- there is no scheduling/budget wrapper here
# the way 81_wp4_baselines.py has one, because the WP0 gate already measured
# this paper's four detectors to be cheap (tens of seconds at n ~1000-1600,
# WP5_INVENTORY.md S5); a human runs this WP with an explicit --datasets/
# --methods scope per invocation to stay inside the 10-minute foreground cap,
# the same discipline 13_wp0_gate.R uses for its own datasets/methods args.

suppressPackageStartupMessages({
  library(here)
})

source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))

REPO_ROOT <- here::here()
WP5_DIR   <- file.path(REPO_ROOT, "revision_experiments/results/tr1/wp5")
DATA_DIR  <- file.path(WP5_DIR, "data")

S_MIN <- 0.05
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
RK_METHODS      <- c("U-MCCD", "SU-MCCD")
NN_METHODS      <- c("UN-MCCD", "SUN-MCCD")
PROPOSED_METHODS  <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
BASELINE_METHODS  <- c("LOF", "DBSCAN", "MST", "ODIN", "iForest")
ALL_METHODS       <- c(PROPOSED_METHODS, BASELINE_METHODS)
ALL_DATASETS      <- c("letter", "mnist", "musk", "arrhythmia")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", flag, "="), "", hit[[length(hit)]])
}
flag_set <- function(flag) any(args == paste0("--", flag))

SMOKE <- flag_set("smoke")
if (SMOKE) {
  DATASETS <- "arrhythmia"
  METHODS  <- "LOF"
  OUT_CSV  <- file.path(WP5_DIR, "smoke", "wp5_highd_smoke.csv")
  cat("84_wp5_highd.R: --smoke set -- forcing datasets=arrhythmia, methods=LOF, ",
      "output redirected to results/tr1/wp5/smoke/\n", sep = "")
} else {
  ds_arg <- opt("datasets", paste(ALL_DATASETS, collapse = ","))
  me_arg <- opt("methods", paste(ALL_METHODS, collapse = ","))
  DATASETS <- strsplit(ds_arg, ",")[[1]]
  METHODS  <- strsplit(me_arg, ",")[[1]]
  OUT_CSV  <- file.path(WP5_DIR, "wp5_highd_results.csv")
}

stopifnot(all(DATASETS %in% ALL_DATASETS))
stopifnot(all(METHODS %in% ALL_METHODS))

cat(sprintf("84_wp5_highd.R: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("84_wp5_highd.R: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("84_wp5_highd.R: output   = %s\n", OUT_CSV))

# ---------------------------------------------------------------------------
# data loader -- CSV convention written by 84a_wp5_fetch_convert.py
# ---------------------------------------------------------------------------
load_wp5_dataset <- function(name) {
  path <- file.path(DATA_DIR, paste0(name, ".csv"))
  if (!file.exists(path)) {
    stop(sprintf("load_wp5_dataset(): %s not found. Run 84a_wp5_fetch_convert.py first.", path))
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(names(df)[ncol(df)] == "label")
  d <- ncol(df) - 1
  X <- as.matrix(df[, seq_len(d), drop = FALSE])
  Y <- df[["label"]]
  stopifnot(all(Y %in% c(0, 1)))
  list(X = X, Y = Y, d = d, n = nrow(df))
}

# ---------------------------------------------------------------------------
# RK-degeneracy measurement (WP5_PROTOCOL.md S5) -- written once per
# invocation, unconditional of scope, since it is cheap (a handful of
# .RData loads) and every n/a / caveated row cites it.
# ---------------------------------------------------------------------------
WP5_D_GRID <- c(32, 100, 166, 274)

rk_degeneracy_row <- function(d) {
  q <- rk_quant_label_paper(d)
  prod_path <- file.path(RK_QUANT_TABLE_DIR, sprintf("RK-test-simul_%dd_%s%%.RData", d, q))
  probe_path <- file.path(REPO_ROOT, "revision_experiments/results/tr2/probes",
                           sprintf("RK-test-simul_%dd_%s%%.RData", d, q))
  path <- NA_character_; source_tag <- NA_character_
  frac <- NA_real_; extent <- NA_integer_
  if (file.exists(prod_path)) {
    path <- prod_path; source_tag <- "production table"
  } else if (file.exists(probe_path)) {
    path <- probe_path; source_tag <- "probe table (niter=20)"
  }
  if (!is.na(path)) {
    e <- new.env(); load(path, envir = e)
    simul <- get("simul", envir = e)
    m <- simul$quan[[1]]
    frac <- mean(m == 0)
    extent <- nrow(m)
  } else {
    source_tag <- "no table available at this dimension"
  }
  data.frame(d = d, quant_label = q, zero_frac = frac, extent = extent,
             source = source_tag, file = if (is.na(path)) NA_character_ else path,
             stringsAsFactors = FALSE)
}

rk_degen <- do.call(rbind, lapply(WP5_D_GRID, rk_degeneracy_row))
RK_DEGEN_CSV <- file.path(WP5_DIR, "rk_degeneracy.csv")
dir.create(dirname(RK_DEGEN_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(rk_degen, RK_DEGEN_CSV, row.names = FALSE)
cat("84_wp5_highd.R: wrote ", RK_DEGEN_CSV, "\n", sep = "")
print(rk_degen, row.names = FALSE)

rk_table_exists <- function(d) {
  q <- rk_quant_label_paper(d)
  file.exists(file.path(RK_QUANT_TABLE_DIR, sprintf("RK-test-simul_%dd_%s%%.RData", d, q)))
}

# ---------------------------------------------------------------------------
# per-cell driver
# ---------------------------------------------------------------------------
ROW_COLS <- c("dataset", "method", "n", "d", "min_cls", "TPR", "TNR", "BA", "F2",
              "unassigned_rows", "quant_used", "quant_label", "quant_source",
              "rk_zero_frac", "t_construct", "t_total", "status", "reason",
              "timestamp")

na_row <- function(dataset, method, n, d, min_cls, status, reason) {
  data.frame(
    dataset = dataset, method = method, n = n, d = d, min_cls = min_cls,
    TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
    unassigned_rows = NA_integer_, quant_used = NA_real_,
    quant_label = NA_character_, quant_source = NA_character_,
    rk_zero_frac = NA_real_, t_construct = NA_real_, t_total = NA_real_,
    status = status, reason = reason,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}
TIMEOUT_SEC <- as.numeric(Sys.getenv("WP5_HIGHD_TIMEOUT_SEC", "540"))

run_cell <- function(dataset, method) {
  keys <- c(dataset = dataset, method = method)
  if (has_result(OUT_CSV, keys)) {
    cat(sprintf("[skip, already recorded] %s x %s\n", dataset, method))
    return(invisible(NULL))
  }

  dat <- load_wp5_dataset(dataset)
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n

  # n/a rule (WP5_PROTOCOL.md S5): checked BEFORE any call, not caught after.
  if (method %in% RK_METHODS && !rk_table_exists(d)) {
    row <- na_row(dataset, method, n, d,
                   min_cls = if (method %in% MIN_CLS_METHODS) S_MIN else NA_real_,
                   status = "n/a",
                   reason = sprintf(
                     "RK envelope 100%% zero quantiles; no production RK quantile table at d=%d (see rk_degeneracy.csv)",
                     d))
    stopifnot(identical(names(row), ROW_COLS))
    append_result(OUT_CSV, row)
    cat(sprintf("  %-11s x %-9s d=%-3d | n/a: %s\n", dataset, method, d, row$reason))
    return(invisible(row))
  }

  min_cls_arg <- if (method %in% MIN_CLS_METHODS) S_MIN else NA_real_
  extra <- if (method %in% MIN_CLS_METHODS) list(min.cls = S_MIN) else list()

  cell_result <- tryCatch({
    with_timeout(function() {
      res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra))
      if (length(res$score) != n || anyNA(res$score)) {
        stop(sprintf("score sanity check failed: length=%d (n=%d), any NA=%s",
                      length(res$score), n, anyNA(res$score)))
      }
      thr <- REAL_DATA_THRESHOLDS[[method]]
      m <- evaluate(Y, res$score, thr)

      rk_zero_frac <- if (method %in% RK_METHODS) {
        rk_degen$zero_frac[rk_degen$d == d]
      } else NA_real_

      row <- data.frame(
        dataset = dataset, method = method, n = n, d = d, min_cls = min_cls_arg,
        TPR = unname(m["TPR"]), TNR = unname(m["TNR"]), BA = unname(m["BA"]), F2 = unname(m["F2"]),
        unassigned_rows = if (!is.null(res$unassigned_rows)) res$unassigned_rows else NA_integer_,
        quant_used = if (!is.null(res$quant_used)) res$quant_used else NA_real_,
        quant_label = if (!is.null(res$quant_label)) res$quant_label else NA_character_,
        quant_source = if (!is.null(res$quant_source)) res$quant_source else NA_character_,
        rk_zero_frac = rk_zero_frac,
        t_construct = if (!is.null(res$t_construct)) res$t_construct else NA_real_,
        t_total = res$t_total,
        status = "ok", reason = "",
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      cat(sprintf("  %-11s x %-9s n=%-5d d=%-3d TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | quant=%s | t=%.2fs\n",
                  dataset, method, n, d, m["TPR"], m["TNR"], m["BA"], m["F2"],
                  row$quant_label, row$t_total))
      row
    }, TIMEOUT_SEC)
  }, error = function(e) {
    is_timeout <- grepl("elapsed time limit|reached CPU time limit", conditionMessage(e))
    status <- if (is_timeout) sprintf("timeout(>%.0fs)", TIMEOUT_SEC) else "error"
    cat(sprintf("  %-11s x %-9s %s: %s\n", dataset, method,
                if (is_timeout) "TIMEOUT" else "ERROR", conditionMessage(e)))
    na_row(dataset, method, n, d, min_cls_arg, status = status,
           reason = substr(conditionMessage(e), 1, 300))
  })

  stopifnot(identical(names(cell_result), ROW_COLS))
  append_result(OUT_CSV, cell_result)
  invisible(cell_result)
}

for (dataset in DATASETS) {
  for (method in METHODS) {
    run_cell(dataset, method)
  }
}

cat("\n84_wp5_highd.R: done.\n")
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  cat(sprintf("%s now has %d rows (status counts: %s).\n", OUT_CSV, nrow(final),
              paste(sprintf("%s=%d", names(table(final$status)), table(final$status)), collapse = ", ")))
}
