#!/usr/bin/env Rscript
# 40_wp2a_tolerance_sweep.R -- WP2(a): sensitivity of the four MCCD detectors
# to the ONE hard-coded knob inside the KS-CCD MCG validation step.
#
#   Rscript 40_wp2a_tolerance_sweep.R [datasets] [methods] [tolerances]
#   Rscript 40_wp2a_tolerance_sweep.R --summary
#
# THE KNOB.  R/ccds/mKNN_CCD_functions.R:50-88 defines connected.ksccd.m(t),
# which recovers, per macro-cluster, the largest density parameter delta for
# which that cluster's KS-CCD is still connected. It does so by a
# bracketing + geometric-refinement search whose convergence tolerance is the
# literal 0.01 in the outer loop guard
#
#     while (m/step <= 10 | step > 0.01)          # <- the tolerance
#
# and it is not exposed as a function argument anywhere. All four detectors
# call the function as sapply(catch.info, connected.ksccd.m) with no way to
# reach it.
#
# THE CLAIM UNDER TEST.  The recovered delta moves with the tolerance, but the
# outlier LABELS do not. This script is written so that a negative result is
# recorded as plainly as a positive one; nothing here is tuned toward label
# stability.
#
# HOW THE SWEEP IS DONE.  connected.ksccd.m.tol() below is a byte-faithful copy
# of the production body with the literal 0.01 replaced by a `tol` argument.
# Nothing in the shared codebase is modified. The production function is then
# shadowed by a global binding installed AFTER every source() has run:
# methods/outlier_detection/*.R each source mKNN_CCD_functions.R at their
# line 2, and wp0_mccd_methods.R sources those, so an override installed any
# earlier would be silently clobbered and every tolerance would secretly run
# at 0.01 -- a failure mode indistinguishable from a perfect stability result.
# The override is therefore (re)installed immediately before every single cell,
# and 41_* / the --summary mode assert that recovered deltas actually differ
# across tolerances before any stability claim is allowed to stand.
#
# The shadow also acts as the delta tap: the detectors do not return delta.max
# / delta_max, so the only way to observe the object the tolerance directly
# controls is to record it as it is produced.
#
# Checkpointed on (dataset, method, tol): a restart skips completed cells.
#
# CAVEAT ON delta_values / max_abs_delta_diff.  The CSV stores delta rounded to
# 6 dp, so max_abs_delta_diff = 0 means "equal to 6 dp", NOT "bit-identical".
# One cell in the 8-dataset grid, WBC x SU-MCCD, is exactly that case: delta is
# 9.0 to 6 dp at every tolerance, yet the label set moves. `--probe WBC SU-MCCD`
# resolves it -- delta is 8.9999999999999964 at tol in {1, 0.1} and
# 8.9999999999999911 at tol in {0.01, 0.001}, a 5.3e-15 difference, and WBC's
# integer-valued features put pairwise distances exactly on the boundary of
# ksccd.connected's `<=` radius test, so 3 points change side. The same probe
# repeats the reference tolerance and confirms the pipeline is deterministic,
# so this is floating-point sensitivity, not run-to-run noise.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- only now is it safe to define and install the override -----------------

args <- commandArgs(trailingOnly = TRUE)
SUMMARY_ONLY <- any(args == "--summary")
args <- args[args != "--summary"]
PROBE_AT <- which(args == "--probe")
PROBE_ARGS <- if (length(PROBE_AT)) args[PROBE_AT + 1:2] else NULL
if (length(PROBE_AT)) args <- args[-(PROBE_AT + 0:2)]

DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else
              c("hepatitis", "glass", "vertebral", "stamps")
METHODS  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else
              c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
TOLS     <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(strsplit(args[3], ",")[[1]]) else
              c(1, 0.1, 0.01, 0.001)

REF_TOL      <- 0.01      # production reference
S_MIN        <- 0.0625    # proportion of n, not a count (see 28_regen_final.R)
CELL_TIMEOUT <- 300       # seconds; a cell over this is recorded as timed out

OUT_CSV     <- here::here("revision_experiments/results/tr1/wp2a_tolerance_sweep.csv")
SUMMARY_CSV <- here::here("revision_experiments/results/tr1/wp2a_tolerance_summary.csv")

MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")

# ---------------------------------------------------------------------------
# Tolerance-parameterised copy of connected.ksccd.m (R/ccds/mKNN_CCD_functions.R
# lines 50-88). The ONLY edit is `step > tol` in place of `step > 0.01`.
# ---------------------------------------------------------------------------
connected.ksccd.m.tol <- function(t, tol) {
  source(here::here("R/ccds/ccd_ks_NEW.R"))
  source(here::here("R/ccds/ccdfunctions.R"))
  if (is.null(nrow(t))) {
    m <- 1
  } else {
    m.intial <- 1
    m <- m.intial
    member <- ksccd.connected(t, m, sequential = FALSE, alpha = 0.05)$member
    while (length(unique(member)) == 1) {
      m <- m * 10
      member <- ksccd.connected(t, m, sequential = FALSE, alpha = 0.05)$member
    }

    step <- m

    while (m / step <= 10 | step > tol) {      # <-- swept knob

      step <- step / 10
      while (length(unique(member)) > 1) {
        m <- m - step
        member <- ksccd.connected(t, m, sequential = FALSE, alpha = 0.05)$member
      }

      step <- step / 10
      while (length(unique(member)) == 1) {
        m <- m + step
        member <- ksccd.connected(t, m, sequential = FALSE, alpha = 0.05)$member
      }
      m <- m - step
      if (m / step >= 1e+12) break
    }
  }
  return(m)
}

# Delta tap + current tolerance. An environment, not a bare global, so the
# recursive source() calls inside connected.ksccd.m.tol cannot stomp it.
WP2A <- new.env(parent = emptyenv())
WP2A$tol    <- REF_TOL
WP2A$deltas <- numeric(0)

install_override <- function(tol) {
  WP2A$tol    <- tol
  WP2A$deltas <- numeric(0)
  assign("connected.ksccd.m", function(t) {
    v <- connected.ksccd.m.tol(t, WP2A$tol)
    WP2A$deltas <- c(WP2A$deltas, v)
    v
  }, envir = globalenv())
}

collapse_num <- function(v) if (length(v) == 0) "" else paste(round(v, 6), collapse = ";")
collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
parse_num    <- function(s) if (is.na(s) || !nzchar(s)) numeric(0) else as.numeric(strsplit(s, ";")[[1]])
parse_idx    <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else as.integer(strsplit(s, ";")[[1]])

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
run_grid <- function() {
  cat(sprintf("40_wp2a: datasets  = %s\n", paste(DATASETS, collapse = ", ")))
  cat(sprintf("40_wp2a: methods   = %s\n", paste(METHODS, collapse = ", ")))
  cat(sprintf("40_wp2a: tols      = %s   (reference %g)\n", paste(TOLS, collapse = ", "), REF_TOL))
  cat(sprintf("40_wp2a: S_min     = %g   out = %s\n", S_MIN, OUT_CSV))

  for (ds in DATASETS) {
    dat <- load_real_dataset(ds)
    for (meth in METHODS) {
      for (tol in TOLS) {
        keys <- c(dataset = ds, method = meth, tol = as.character(tol))
        if (isTRUE(has_result(OUT_CSV, keys))) {
          cat(sprintf("[skip] %s x %s x tol=%g\n", ds, meth, tol)); next
        }

        install_override(tol)
        t0 <- Sys.time()
        out <- tryCatch({
          setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
          res <- if (meth %in% MIN_CLS_METHODS) {
            METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
          } else {
            METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
          }
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
               score = res$score, status = "ok", note = NA_character_)
        }, error = function(e) {
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          msg <- conditionMessage(e)
          list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
               score = NULL,
               status = if (grepl("elapsed time limit", msg)) "timeout" else "error",
               note = msg)
        })
        wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        deltas <- WP2A$deltas
        flagged <- if (is.null(out$score)) integer(0) else which(out$score > 0.5)

        append_result(OUT_CSV, list(
          dataset = ds, method = meth, tol = tol, n = dat$n, d = dat$d,
          TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
          BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
          n_flagged   = if (is.null(out$score)) NA_integer_ else length(flagged),
          n_clusters  = length(deltas),
          delta_values = collapse_num(deltas),
          flagged_idx  = collapse_idx(flagged),
          elapsed_sec  = wall,
          s_min        = if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_,
          status = out$status, note = out$note, timestamp = format(Sys.time())
        ))

        cat(sprintf("  %-10s x %-9s tol=%-6g  BA=%s F2=%s  nflag=%-4s ncls=%-2d  delta=[%s]  %s  %.1fs\n",
                    ds, meth, tol,
                    if (is.na(out$m[["BA"]])) "  NA " else sprintf("%.3f", out$m[["BA"]]),
                    if (is.na(out$m[["F2"]])) "  NA " else sprintf("%.3f", out$m[["F2"]]),
                    if (is.null(out$score)) "NA" else length(flagged),
                    length(deltas), collapse_num(deltas), out$status, wall))
        flush.console()
      }
    }
  }
  cat("40_wp2a: grid done\n")
}

# ---------------------------------------------------------------------------
# Summary vs. the tol = REF_TOL reference
# ---------------------------------------------------------------------------
build_summary <- function() {
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE, colClasses = c(
    delta_values = "character", flagged_idx = "character"))
  df$delta_values[is.na(df$delta_values)] <- ""
  df$flagged_idx[is.na(df$flagged_idx)]   <- ""

  rows <- list()
  for (ds in unique(df$dataset)) {
    for (meth in unique(df$method[df$dataset == ds])) {
      sub <- df[df$dataset == ds & df$method == meth, ]
      ref <- sub[abs(sub$tol - REF_TOL) < 1e-12, ]
      if (nrow(ref) == 0) next
      ref <- ref[nrow(ref), ]
      ref_idx <- parse_idx(ref$flagged_idx)
      ref_del <- parse_num(ref$delta_values)

      for (i in seq_len(nrow(sub))) {
        r <- sub[i, ]
        cur_idx <- parse_idx(r$flagged_idx)
        cur_del <- parse_num(r$delta_values)

        inter <- length(intersect(cur_idx, ref_idx))
        uni   <- length(union(cur_idx, ref_idx))
        jac   <- if (uni == 0) 1 else inter / uni
        flips <- length(setdiff(cur_idx, ref_idx)) + length(setdiff(ref_idx, cur_idx))
        ncdiff <- length(cur_del) != length(ref_del)
        mad <- if (ncdiff || length(cur_del) == 0) NA_real_ else max(abs(cur_del - ref_del))

        ok <- r$status == "ok" && ref$status == "ok"
        rows[[length(rows) + 1]] <- data.frame(
          dataset = ds, method = meth, tol = r$tol, ref_tol = REF_TOL,
          n = r$n, d = r$d,
          jaccard = if (ok) jac else NA_real_,
          n_label_flips = if (ok) flips else NA_integer_,
          n_flagged = r$n_flagged, ref_n_flagged = ref$n_flagged,
          max_abs_delta_diff = if (ok) mad else NA_real_,
          identical_labels = if (ok) identical(cur_idx, ref_idx) else NA,
          n_clusters = length(cur_del), ref_n_clusters = length(ref_del),
          n_clusters_differs = ncdiff,
          delta_identical = identical(r$delta_values, ref$delta_values),
          BA = r$BA, ref_BA = ref$BA, F2 = r$F2, ref_F2 = ref$F2,
          status = r$status, ref_status = ref$status, note = r$note,
          elapsed_sec = r$elapsed_sec,
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$dataset, out$method, out$tol), ]
  dir.create(dirname(SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("40_wp2a: wrote %s (%d rows)\n", SUMMARY_CSV, nrow(out)))

  # --- OVERRIDE GUARD -------------------------------------------------------
  # If deltas are bit-identical across every tolerance, the override did NOT
  # take effect and any "labels are stable" reading is meaningless.
  cat("\n==== OVERRIDE GUARD: did delta actually move with tol? ====\n")
  moved <- 0L; checked <- 0L
  for (ds in unique(df$dataset)) {
    for (meth in unique(df$method[df$dataset == ds])) {
      sub <- df[df$dataset == ds & df$method == meth & df$status == "ok", ]
      lo <- sub[abs(sub$tol - 1) < 1e-12, ]
      hi <- sub[abs(sub$tol - 0.001) < 1e-12, ]
      if (nrow(lo) == 0 || nrow(hi) == 0) next
      checked <- checked + 1L
      same <- identical(lo$delta_values[1], hi$delta_values[1])
      if (!same) moved <- moved + 1L
      cat(sprintf("  %-10s %-9s tol=1    delta=[%s]\n", ds, meth, lo$delta_values[1]))
      cat(sprintf("  %-10s %-9s tol=.001 delta=[%s]   -> %s\n", ds, meth,
                  hi$delta_values[1], if (same) "IDENTICAL" else "DIFFERS"))
    }
  }
  cat(sprintf("  %d of %d (dataset,method) pairs show delta moving between tol=1 and tol=0.001\n",
              moved, checked))
  if (checked > 0 && moved == 0) {
    cat("  *** GUARD FAILED: delta never moved. The override did not take effect. ***\n")
  }

  # --- Label stability ------------------------------------------------------
  cat("\n==== LABEL STABILITY vs tol = 0.01 ====\n")
  nonref <- out[abs(out$tol - REF_TOL) > 1e-12, ]
  ok <- nonref[nonref$status == "ok" & nonref$ref_status == "ok", ]
  cat(sprintf("  cells compared (excluding the reference itself): %d\n", nrow(ok)))
  cat(sprintf("  identical labels: %d / %d\n", sum(ok$identical_labels), nrow(ok)))
  mv <- ok[!ok$identical_labels, ]
  if (nrow(mv) == 0) {
    cat("  no cell moved its labels.\n")
  } else {
    cat("  cells where labels MOVED:\n")
    for (i in seq_len(nrow(mv))) {
      cat(sprintf("    %-10s %-9s tol=%-6g jaccard=%.4f flips=%d  nflag %d vs ref %d\n",
                  mv$dataset[i], mv$method[i], mv$tol[i], mv$jaccard[i],
                  mv$n_label_flips[i], mv$n_flagged[i], mv$ref_n_flagged[i]))
    }
  }
  dd <- ok$max_abs_delta_diff[!is.na(ok$max_abs_delta_diff)]
  if (length(dd)) {
    k <- which.max(ok$max_abs_delta_diff)
    cat(sprintf("  largest max_abs_delta_diff: %g  (%s / %s / tol=%g)\n",
                max(dd), ok$dataset[k], ok$method[k], ok$tol[k]))
  }
  bad <- out[out$status != "ok", ]
  cat(sprintf("\n==== NON-OK CELLS: %d ====\n", nrow(bad)))
  if (nrow(bad)) for (i in seq_len(nrow(bad))) {
    cat(sprintf("    %-10s %-9s tol=%-6g %s :: %s\n", bad$dataset[i], bad$method[i],
                bad$tol[i], bad$status[i], bad$note[i]))
  }
  cat(sprintf("\n==== CLUSTER-COUNT CHANGES: %d ====\n", sum(out$n_clusters_differs)))
  cc <- out[out$n_clusters_differs, ]
  if (nrow(cc)) for (i in seq_len(nrow(cc))) {
    cat(sprintf("    %-10s %-9s tol=%-6g ncls=%d vs ref %d\n", cc$dataset[i], cc$method[i],
                cc$tol[i], cc$n_clusters[i], cc$ref_n_clusters[i]))
  }
  invisible(out)
}

# ---------------------------------------------------------------------------
# Probe: full-precision deltas + a same-tolerance repeat, for one cell.
#
# Needed because the grid turned up a cell (WBC x SU-MCCD) whose recorded
# delta is identical to 6 dp across all four tolerances while the label set
# still moves. Two candidate causes -- (i) the pipeline is not deterministic,
# so the difference has nothing to do with the tolerance; (ii) delta differs
# only in floating-point noise below the 6 dp the CSV records, and lands on an
# exact distance tie in ksccd.connected's `<=` radius test. The repeat at a
# fixed tolerance separates them: identical repeats rule out (i).
#
#   Rscript 40_wp2a_tolerance_sweep.R --probe <dataset> <method>
# ---------------------------------------------------------------------------
run_probe <- function(ds, meth) {
  dat <- load_real_dataset(ds)
  cat(sprintf("PROBE %s x %s  (n=%d d=%d)\n", ds, meth, dat$n, dat$d))
  prev <- NULL
  for (tol in c(TOLS, REF_TOL)) {   # last entry repeats the reference
    install_override(tol)
    res <- if (meth %in% MIN_CLS_METHODS) {
      METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
    } else {
      METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
    }
    fl <- which(res$score > 0.5)
    cat(sprintf("  tol=%-7g nflag=%-4d delta(17 sig)= %s\n", tol, length(fl),
                paste(format(WP2A$deltas, digits = 17), collapse = "; ")))
    if (!is.null(prev) && abs(tol - REF_TOL) < 1e-12) {
      cat(sprintf("  repeat-at-reference vs first reference run: labels %s\n",
                  if (identical(prev, fl)) "IDENTICAL (deterministic)" else "DIFFER (NON-DETERMINISTIC)"))
    }
    if (abs(tol - REF_TOL) < 1e-12 && is.null(prev)) prev <- fl
  }
}

if (!is.null(PROBE_ARGS)) {
  run_probe(PROBE_ARGS[1], PROBE_ARGS[2])
} else if (SUMMARY_ONLY) {
  build_summary()
} else {
  run_grid()
  build_summary()
}
