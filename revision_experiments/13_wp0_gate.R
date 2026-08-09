#!/usr/bin/env Rscript
# revision_experiments/13_wp0_gate.R
#
# WP0 reproduction gate: runs the 4 newly-wired MCCD detectors (U-MCCD,
# SU-MCCD, UN-MCCD, SUN-MCCD; see wp0_mccd_methods.R) on real datasets and
# compares TPR/TNR/BA/F2 against revision_experiments/published_realdata_
# truth.csv to 3 decimals. Does NOT tune anything to force a match -- a
# mismatch is written to the output CSV as-is (match_3dp = FALSE) and
# printed; it is a finding, not a bug to paper over here.
#
# PARAMETER FIX (this revision): the prior run of this gate used
# wp0_mccd_methods.R wrappers that fell back to get_simul()'s generic
# quantile buckets and to min.cls=0. Those wrappers now default to the
# manuscript's own alpha schedule (main text L880; Supplementary Material
# L655-677) for `quant`, and expose `min.cls` for the (genuinely ambiguous,
# main text L605) S_min rule -- see wp0_mccd_methods.R's "Paper-faithful
# parameter resolver" section. This script sweeps every distinct min.cls
# reading from mccd_min_cls_readings() for SU-MCCD/SUN-MCCD and records
# which reading(s) produced each value, rather than picking one.
#
# SCHEMA CHANGE: adds `min_cls`, `min_cls_readings`, `quant_used`,
# `quant_label`, `param_source` columns to record which paper-faithful
# parameter reading produced each row. Because this changes the row schema,
# results are written to a NEW file, results/wp0_gate_v2.csv -- NOT
# appended to the old results/wp0_gate.csv (which used the wrong
# quantile-bucket rule and min.cls=0 throughout; left untouched here, for
# the record of what the pre-fix wiring produced).
#
# USAGE
#   Rscript "revision_experiments/13_wp0_gate.R" [datasets] [methods] [min_cls_reading]
#
#   [datasets]        comma-separated dataset names as accepted by
#                      load_real_dataset(), e.g. "hepatitis,glass,ecoli".
#                      Default: "hepatitis,glass,ecoli" (the WP0 smoke set).
#   [methods]         comma-separated method names as they appear in
#                      METHOD_REGISTRY. Default: "U-MCCD,SU-MCCD,UN-MCCD,SUN-MCCD".
#   [min_cls_reading] "all" (default) sweeps every distinct min.cls reading
#                      from mccd_min_cls_readings() for SU-MCCD/SUN-MCCD; or
#                      one of "a","b","c","d" to run only the single
#                      reading with that label (used for the vowels
#                      runtime-anchor run, where re-running the full sweep
#                      would multiply an already-expensive cell for no
#                      reproduction benefit -- see task brief).
#
# Resumable: checkpoints to results/wp0_gate_v2.csv via the harness's own
# append_result()/has_result(), keyed on (dataset, method, min_cls) -- so
# distinct min.cls readings for the same (dataset, method) do not collide
# or get skipped as "already recorded".
#
# Per-cell isolation: each (dataset, method, min_cls) cell is wrapped in
# tryCatch AND a wall-clock time cap (default 600s/cell, override via
# WP0_GATE_TIMEOUT_SEC env var -- e.g. 1800 for the vowels 30-minute cap)
# so one hang or error does not take down the whole run. A failed or
# timed-out cell is recorded with status columns left NA and is reported,
# not silently dropped.

source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

REPO_ROOT   <- here::here()
# OUT_CSV: overridable via WP0_GATE_OUT_CSV so a concurrent gate run writing
# to the original wp0_gate_v2.csv is never touched by this one (two writers
# on one CSV can interleave partial appended lines). Default unchanged for
# any other caller that doesn't set the env var.
OUT_CSV     <- Sys.getenv("WP0_GATE_OUT_CSV", file.path(REPO_ROOT, "revision_experiments/results/wp0_gate_v2.csv"))
TRUTH_CSV   <- file.path(REPO_ROOT, "revision_experiments/published_realdata_truth.csv")
TIMEOUT_SEC <- as.numeric(Sys.getenv("WP0_GATE_TIMEOUT_SEC", "600"))

args <- commandArgs(trailingOnly = TRUE)
DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else c("hepatitis", "glass", "ecoli")
METHODS  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
MIN_CLS_READING <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "all"
# "half_contam"/"full_contam": UNITS FIX readings -- min.cls passed as a
# PROPORTION of n (0.5*(n0/n) and n0/n respectively) via
# mccd_min_cls_proportion_readings(), instead of the pre-fix "all"/"a".."d"
# sweep's raw counts (see wp0_mccd_methods.R's "UNITS FIX" section).
stopifnot(MIN_CLS_READING %in% c("all", "a", "b", "c", "d", "half_contam", "full_contam"))

stopifnot(all(METHODS %in% names(METHOD_REGISTRY)))

cat(sprintf("13_wp0_gate.R: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("13_wp0_gate.R: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("13_wp0_gate.R: min_cls_reading = %s\n", MIN_CLS_READING))
cat(sprintf("13_wp0_gate.R: per-cell timeout = %.0fs, output = %s\n", TIMEOUT_SEC, OUT_CSV))

if (!file.exists(TRUTH_CSV)) stop("Published truth CSV not found: ", TRUTH_CSV)
truth <- read.csv(TRUTH_CSV, stringsAsFactors = FALSE)

get_truth <- function(dataset, method) {
  sub <- truth[truth$dataset == dataset & truth$method == method, ]
  if (nrow(sub) == 0) return(NULL)
  v <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
  for (i in seq_len(nrow(sub))) {
    if (sub$metric[i] %in% names(v)) v[sub$metric[i]] <- sub$value[i]
  }
  v
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}

# ---------------------------------------------------------------------------
# Table 5 (tab:Real_Data, main text lines 1063-1080) published contamination
# fractions -- needed ONLY for min.cls reading (d). Known to disagree with
# the loader's actual n0 for glass (Table 5: 10 outliers; loader and the
# published results tables agree on 9) and ecoli (Table 5: 9; loader and
# results tables agree on 8) -- see task brief. `n` used in reading (d) is
# always the LOADER's n, not Table 5's (which also differs slightly for
# glass: 214 vs 213 after de-duplication).
#
# The four datasets added below were checked directly against
# load_real_dataset() before adding these entries (see WP0 report):
# vertebral (n=240, n0=30), stamps (n=340, n0=31) and waveform (n=3443,
# n0=100) match Table 5's n AND outlier count exactly. wilt does NOT: Table 5
# states n=4735, n0=257 (5.4%), but load_real_dataset("wilt") returns n=4819,
# n0=257 -- the outlier count matches but the loader's arff file
# (Wilt_withoutdupl_05.arff) has 84 more regular-class rows than Table 5
# implies (4562 vs 4478), and distinct() removes zero duplicate rows from it
# (verified directly), so this is not a de-duplication difference of the
# kind seen for glass/ecoli. Reading (d) below uses Table 5's published
# fraction 0.054 regardless, per this function's own contract; the
# discrepancy is reported, not resolved, here.
TABLE5_CONTAMINATION_PCT <- c(
  hepatitis = 0.095,
  glass     = 0.045,
  ecoli     = 0.026,
  vowels    = 0.034,
  vertebral = 0.125,
  stamps    = 0.091,
  waveform  = 0.029,
  wilt      = 0.054
)

MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")

ROW_COLS <- c("dataset", "method", "n", "d", "min_cls", "min_cls_readings", "min_cls_units",
              "TPR", "TNR", "BA", "F2",
              "published_TPR", "published_TNR", "published_BA", "published_F2",
              "max_abs_diff", "match_3dp", "unassigned_rows",
              "quant_used", "quant_label", "param_source",
              "t_construct", "t_total", "status", "timestamp")

na_row <- function(dataset, method, n, d, min_cls, min_cls_readings, min_cls_units, status) {
  data.frame(
    dataset = dataset, method = method, n = n, d = d,
    min_cls = min_cls, min_cls_readings = min_cls_readings, min_cls_units = min_cls_units,
    TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
    published_TPR = NA_real_, published_TNR = NA_real_, published_BA = NA_real_, published_F2 = NA_real_,
    max_abs_diff = NA_real_, match_3dp = NA, unassigned_rows = NA_integer_,
    quant_used = NA_real_, quant_label = NA_character_, param_source = NA_character_,
    t_construct = NA_real_, t_total = NA_real_, status = status,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

run_cell <- function(dataset, method, min_cls = NA_integer_, min_cls_readings = NA_character_, min_cls_units = NA_character_) {
  keys <- c(dataset = dataset, method = method, min_cls = as.character(min_cls))
  if (has_result(OUT_CSV, keys)) {
    cat(sprintf("[skip, already recorded] %s x %s (min_cls=%s)\n", dataset, method, min_cls))
    return(invisible(NULL))
  }

  cell_result <- tryCatch({
    with_timeout(function() {
      dat <- load_real_dataset(dataset)
      X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n

      t0 <- Sys.time()
      if (method %in% MIN_CLS_METHODS && !is.na(min_cls)) {
        res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = min_cls)
      } else {
        res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y)
      }
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (length(res$score) != n || anyNA(res$score)) {
        stop(sprintf("score sanity check failed: length=%d (n=%d), any NA=%s",
                      length(res$score), n, anyNA(res$score)))
      }

      thr <- REAL_DATA_THRESHOLDS[[method]]
      m <- evaluate(Y, res$score, thr)

      param_source <- sprintf("quant=%s (%s)",
                               if (!is.null(res$quant_label)) res$quant_label else NA_character_,
                               if (!is.null(res$quant_source)) res$quant_source else NA_character_)

      # NOTE: min_cls is kept as a double, not coerced to integer. The
      # pre-fix code below this comment used as.integer(min_cls), which is
      # harmless for the old count-valued sweep (a/b/c/d, always whole
      # numbers) but would truncate a proportion reading like 0.0625 to 0 --
      # exactly the kind of unit confusion this re-run exists to eliminate.
      min_cls_out <- if (is.na(min_cls)) NA_real_ else as.numeric(min_cls)

      base_fields <- list(
        dataset = dataset, method = method, n = n, d = d,
        min_cls = min_cls_out, min_cls_readings = min_cls_readings, min_cls_units = min_cls_units,
        TPR = unname(m["TPR"]), TNR = unname(m["TNR"]), BA = unname(m["BA"]), F2 = unname(m["F2"]),
        unassigned_rows = if (!is.null(res$unassigned_rows)) res$unassigned_rows else NA_integer_,
        quant_used = if (!is.null(res$quant_used)) res$quant_used else NA_real_,
        quant_label = if (!is.null(res$quant_label)) res$quant_label else NA_character_,
        param_source = param_source,
        t_construct = res$t_construct, t_total = res$t_total,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )

      pub <- get_truth(dataset, method)
      if (is.null(pub)) {
        row <- as.data.frame(c(base_fields, list(
          published_TPR = NA_real_, published_TNR = NA_real_, published_BA = NA_real_, published_F2 = NA_real_,
          max_abs_diff = NA_real_, match_3dp = NA, status = "no_published_truth"
        )), stringsAsFactors = FALSE)[ROW_COLS]
        cat(sprintf("  %-10s x %-9s min_cls=%-4s n=%-4d d=%-3d TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | NO PUBLISHED TRUTH | quant=%s | wall=%.2fs\n",
                    dataset, method, min_cls, n, d, m["TPR"], m["TNR"], m["BA"], m["F2"], base_fields$quant_label, wall))
        return(row)
      }

      diffs <- abs(c(m["TPR"], m["TNR"], m["BA"], m["F2"]) - pub[c("TPR", "TNR", "BA", "F2")])
      max_diff <- max(diffs)
      match_3dp <- all(round(c(m["TPR"], m["TNR"], m["BA"], m["F2"]), 3) == round(pub[c("TPR", "TNR", "BA", "F2")], 3))

      row <- as.data.frame(c(base_fields, list(
        published_TPR = unname(pub["TPR"]), published_TNR = unname(pub["TNR"]),
        published_BA = unname(pub["BA"]), published_F2 = unname(pub["F2"]),
        max_abs_diff = unname(max_diff), match_3dp = match_3dp, status = "ok"
      )), stringsAsFactors = FALSE)[ROW_COLS]

      cat(sprintf("  %-10s x %-9s min_cls=%-4s n=%-4d d=%-3d TPR=%.3f/%.3f TNR=%.3f/%.3f BA=%.3f/%.3f F2=%.3f/%.3f | %s | quant=%s | wall=%.2fs\n",
                  dataset, method, min_cls, n, d, m["TPR"], pub["TPR"], m["TNR"], pub["TNR"],
                  m["BA"], pub["BA"], m["F2"], pub["F2"],
                  if (match_3dp) "MATCH" else "MISMATCH", base_fields$quant_label, wall))
      row
    }, TIMEOUT_SEC)
  }, error = function(e) {
    is_timeout <- grepl("elapsed time limit|reached CPU time limit", conditionMessage(e))
    status_prefix <- if (is_timeout) sprintf("timeout(>%.0fs): ", TIMEOUT_SEC) else "error: "
    cat(sprintf("  %-10s x %-9s min_cls=%s %s: %s\n", dataset, method, min_cls,
                if (is_timeout) "TIMEOUT" else "ERROR", conditionMessage(e)))
    na_row(dataset, method, NA_integer_, NA_integer_,
           min_cls = if (is.na(min_cls)) NA_real_ else as.numeric(min_cls),
           min_cls_readings = min_cls_readings,
           min_cls_units = min_cls_units,
           status = paste0(status_prefix, substr(conditionMessage(e), 1, 200)))
  })

  stopifnot(identical(names(cell_result), ROW_COLS))
  append_result(OUT_CSV, cell_result)
  invisible(cell_result)
}

for (dataset in DATASETS) {
  dat_meta <- load_real_dataset(dataset)  # cheap: RealData_Collection.R is only source()'d once per session
  n0 <- sum(dat_meta$Y == 0)
  n  <- dat_meta$n
  contam_pct <- TABLE5_CONTAMINATION_PCT[[dataset]]

  for (method in METHODS) {
    if (method %in% MIN_CLS_METHODS) {
      if (MIN_CLS_READING %in% c("half_contam", "full_contam")) {
        # UNITS FIX readings: min.cls as a PROPORTION of n (see
        # wp0_mccd_methods.R's mccd_min_cls_proportion_readings()). Does not
        # touch the pre-fix mccd_min_cls_readings() sweep below at all.
        prop_readings <- mccd_min_cls_proportion_readings(n0, n)
        run_cell(dataset, method,
                  min_cls = prop_readings[[MIN_CLS_READING]],
                  min_cls_readings = MIN_CLS_READING,
                  min_cls_units = "proportion")
      } else {
        if (is.null(contam_pct)) {
          stop(sprintf("No Table 5 contamination fraction registered for dataset '%s' -- required for min.cls readings on method '%s'. Add it to TABLE5_CONTAMINATION_PCT.",
                        dataset, method))
        }
        readings <- mccd_min_cls_readings(n0, n, contam_pct)  # value(as char) -> reading labels

        if (MIN_CLS_READING == "all") {
          for (val_chr in names(readings)) {
            run_cell(dataset, method,
                      min_cls = as.integer(val_chr),
                      min_cls_readings = paste(readings[[val_chr]], collapse = ","),
                      min_cls_units = "count")
          }
        } else {
          match_val <- NULL
          for (val_chr in names(readings)) {
            if (MIN_CLS_READING %in% readings[[val_chr]]) { match_val <- as.integer(val_chr); break }
          }
          if (is.null(match_val)) {
            stop(sprintf("min_cls_reading '%s' not found among computed readings for dataset '%s' (readings: %s)",
                          MIN_CLS_READING, dataset, paste(names(readings), collapse = ", ")))
          }
          run_cell(dataset, method, min_cls = match_val, min_cls_readings = MIN_CLS_READING, min_cls_units = "count")
        }
      }
    } else {
      run_cell(dataset, method, min_cls = NA_integer_, min_cls_readings = NA_character_, min_cls_units = NA_character_)
    }
  }
}

cat("\n13_wp0_gate.R: done.\n")
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  cat(sprintf("%s now has %d rows.\n", OUT_CSV, nrow(final)))
  n_ok <- sum(final$status == "ok", na.rm = TRUE)
  n_match <- sum(final$match_3dp, na.rm = TRUE)
  cat(sprintf("status=ok: %d/%d ; match_3dp=TRUE: %d/%d\n", n_ok, nrow(final), n_match, nrow(final)))
}
