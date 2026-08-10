#!/usr/bin/env Rscript
# 43_wp2a_direction_sweep.R -- WP2(a) radius-search direction sweep (ascend vs
# descend), the last row of the parameter-sensitivity table. Reviewer R3.10a
# also lists this as an ablation toggle, so one run serves both points.
#
# ---------------------------------------------------------------------------
# STEP 1 FINDINGS (read from source before any code below was written):
#
# THE TOGGLE. UNMCCD_outlier(datax, simul, method="ascend", niter)
# (methods/outlier_detection/UN-MCCD.R L7) and SUNMCCD_outlier(datax, simul,
# min.cls, method="ascend", low.num) (methods/outlier_detection/SUN-MCCD.R
# L9) both pass `method` straight into nnccd_clustering_quantile(), which
# passes it straight into nnccd.radi(dx, quantile, method, low.num, quant,
# simul, niter, scores) (R/ccds/UN_CCD.R L233-323). harness.R's own wrappers
# already expose it: wp0_mccd_methods.R L316 (unmccd_method) and L333
# (sunmccd_method) both declare `method = "ascend"` as a pass-through
# parameter of the METHOD_REGISTRY entry. No override of harness.R or
# wp0_mccd_methods.R is needed to flip direction -- it is a normal call
# argument: METHOD_REGISTRY[["UN-MCCD"]](X=, d=, Y=, method="descend").
#
# WHAT ASCEND/DESCEND ACTUALLY DO (R/ccds/UN_CCD.R L233-323). For each point
# i, "ascend" (L245-264) starts from a ball holding only its `low.num`
# nearest neighbours and GROWS the ball outward one point at a time
# (increasing radius, j = low.num..n), testing at each step whether the
# within-ball NN-distance statistic (mean and median, center point excluded)
# has fallen below the CSR lower envelope; it stops and reports the radius
# from the PREVIOUS step (j-1) the first time that happens -- i.e. it takes
# the SMALLEST radius at which the neighbourhood first looks significantly
# too tight relative to a Poisson null. "descend" (L265-279, L303-317) starts
# from the ball containing ALL n-1 other points (the largest possible
# radius) and SHRINKS it inward one point at a time (j = 1..n-low.num,
# points dropped in decreasing-distance order), stopping and reporting the
# CURRENT radius the first time the statistic is already back ABOVE the
# envelope -- i.e. it takes the LARGEST radius that still looks
# CSR-consistent. Both searches are hunting for the same
# CSR/non-CSR crossing along the distance-sorted neighbour sequence for
# point i, approached from opposite ends. If that crossing is monotonic
# (statistic drops below the envelope exactly once as the ball grows and
# never recovers), ascend and descend land on the same radius. If it is not
# monotonic -- the statistic dips below the envelope and climbs back above
# it more than once, which can happen with the small, noisy `low.num`-sized
# balls this test uses -- the two directions can lock onto genuinely
# different radii for the same point. That non-monotonicity is the entire
# mechanism by which this sweep can find a difference.
#
# RK-BASED DETECTORS HAVE NO ANALOGOUS TOGGLE. RUMCCD_outlier
# (methods/outlier_detection/RU-MCCDs.R) and SUMCCD_outlier
# (methods/outlier_detection/SU-MCCDs.R) take no `method` argument at all
# (verified against their full signatures). Their underlying radius routine,
# RKCCD_correct_quant -> rccd.clustering_correct_quantile (R/ccds/RK_CCD_New.R
# L19-46), does have a `method` argument, but it selects between
# "non-dynamic"/"dynamic" edge correction for the K-function estimate inside
# ccd.Kest.edge.quantile() (L158-224) -- an unrelated axis (edge-correction
# scheme, not search direction) -- and RUMCCD_outlier/SUMCCD_outlier both
# hard-code it to "non-dynamic" via RKCCD_correct_quant's default
# (RK_CCD_New.R L225-247) with no path from the outlier-detector call to
# override it. wp0_mccd_methods.R's own wrappers confirm this from the
# harness side: umccd_method (L282) and sumccd_method (L299) declare no
# `method` parameter at all, unlike unmccd_method/sunmccd_method. So the
# ascend/descend toggle -- and R3.10a's ablation -- reaches only half of the
# 2x2 design family (the NND-based half: UN-MCCD, SUN-MCCD), not
# U-MCCD/SU-MCCD. This asymmetry is itself worth stating in the response
# letter.
#
# VERIFICATION THAT THE ARGUMENT GENUINELY REACHES nnccd.radi. nnccd.radi is
# shadowed by a thin tap installed AFTER wp0_mccd_methods.R has fully sourced
# (so the tap is not clobbered by a later source() -- neither UN-MCCD.R nor
# SUN-MCCD.R nor UN_CCD.R is re-sourced per detector call, only once at
# startup via wp0_mccd_methods.R's `if (!exists(...))` guards). The tap calls
# straight through to the REAL nnccd.radi (captured into ORIG_NNCCD_RADI
# before the shadow is installed) and records the `method` value it actually
# received plus the resulting R vector. Every CSV row therefore carries
# nnradi_method_received (must equal the requested direction) and a radii
# signature (radii_sig, radii_mean); the analysis step below refuses to
# report "no effect" unless it can show at least one cell where ascend and
# descend produced different radii_sig -- same discipline
# 40_wp2a_tolerance_sweep.R used for its own override guard.
#
# Usage:
#   Rscript 43_wp2a_direction_sweep.R [datasets] [methods] [directions] [time_budget_sec]
#   Rscript 43_wp2a_direction_sweep.R --summary
#
# Checkpointed on (dataset, method, direction): a restart skips completed
# cells. TIME_BUDGET (default 480s) stops the grid early and prints what is
# left, so a single invocation fits inside the 600s PowerShell tool timeout;
# re-invoke (same command) to continue -- has_result() skips finished cells.
# ---------------------------------------------------------------------------

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- capture the REAL nnccd.radi before installing any tap ------------------
ORIG_NNCCD_RADI <- get("nnccd.radi", envir = globalenv())

NNRADI_TAP <- new.env(parent = emptyenv())
NNRADI_TAP$calls <- list()

install_nnradi_tap <- function() {
  NNRADI_TAP$calls <- list()
  assign("nnccd.radi", function(dx, quantile = "lower", method = "ascend", low.num, quant,
                                 simul = NULL, niter, scores = F) {
    res <- ORIG_NNCCD_RADI(dx = dx, quantile = quantile, method = method, low.num = low.num,
                            quant = quant, simul = simul, niter = niter, scores = scores)
    NNRADI_TAP$calls[[length(NNRADI_TAP$calls) + 1]] <- list(method_received = method, R = res$R)
    res
  }, envir = globalenv())
}
install_nnradi_tap()

# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
SUMMARY_ONLY <- any(args == "--summary")
MARK_NOT_RUN <- any(args == "--mark-not-run")
args <- args[!(args %in% c("--summary", "--mark-not-run"))]

ALL_DATASETS <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli", "stamps", "WDBC",
                   "pima", "shuffle", "vowels", "PenDigits", "waveform", "thyroid", "pageblocks", "wilt")

DATASETS   <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else ALL_DATASETS
METHODS    <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("UN-MCCD", "SUN-MCCD")
DIRECTIONS <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else c("ascend", "descend")
TIME_BUDGET <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 480

S_MIN <- 0.0625     # manuscript's min.cls proportion (see CLAUDE.md, 40_wp2a)
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
CELL_TIMEOUT <- 480  # per-cell hard cap (s); large real sets can be slow

OUT_CSV      <- here::here("revision_experiments/results/tr1/wp2a_direction_sweep.csv")
FINDINGS_MD  <- here::here("revision_experiments/results/tr1/wp2a_direction_findings.md")

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
radii_sig    <- function(R) if (is.null(R) || length(R) == 0) NA_character_ else paste(round(sort(R), 8), collapse = ";")
parse_idx    <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else as.integer(strsplit(s, ";")[[1]])

# ---------------------------------------------------------------------------
# --mark-not-run: record cells that were deliberately not attempted, with the
# measured per-cell cost that justified stopping (honest partial coverage,
# not silent truncation). Hardcoded n/d from the task's orientation table.
# ---------------------------------------------------------------------------
mark_not_run <- function() {
  info <- list(vowels = c(1452, 12), PenDigits = c(3200, 16), waveform = c(3443, 21),
               thyroid = c(3656, 6), pageblocks = c(4795, 10), wilt = c(4819, 5), shuffle = c(1013, 9))
  note_ascend <- paste("not_run: ascend cost scales steeply with n (measured ok: PenDigits n=3200",
    "194-355s, waveform n=3443 262-446s, thyroid n=3656 265-351s; pageblocks n=4795 UN-MCCD ascend",
    "hit the 480s per-cell cap and was killed -- see pageblocks/UN-MCCD/ascend status=timeout row).",
    "wilt (n=4819, larger than pageblocks) and pageblocks/SUN-MCCD were not attempted on that basis.")
  note_descend <- paste("not_run: descend cost scales far worse than ascend (measured: shuffle n=1013",
    "UN-MCCD descend hit the 480s per-cell cap and was killed -- see shuffle/UN-MCCD/descend",
    "status=timeout row; the 8 small sets n<=367 already showed descend running 10-96s vs <1s for",
    "ascend at the same n). No descend cell at n>=1013 was attempted beyond that probe.")
  rows <- list(
    list("pageblocks", "SUN-MCCD", "ascend", note_ascend), list("wilt", "UN-MCCD", "ascend", note_ascend),
    list("wilt", "SUN-MCCD", "ascend", note_ascend), list("shuffle", "SUN-MCCD", "descend", note_descend),
    list("vowels", "UN-MCCD", "descend", note_descend), list("vowels", "SUN-MCCD", "descend", note_descend),
    list("PenDigits", "UN-MCCD", "descend", note_descend), list("PenDigits", "SUN-MCCD", "descend", note_descend),
    list("waveform", "UN-MCCD", "descend", note_descend), list("waveform", "SUN-MCCD", "descend", note_descend),
    list("thyroid", "UN-MCCD", "descend", note_descend), list("thyroid", "SUN-MCCD", "descend", note_descend),
    list("pageblocks", "UN-MCCD", "descend", note_descend), list("pageblocks", "SUN-MCCD", "descend", note_descend),
    list("wilt", "UN-MCCD", "descend", note_descend), list("wilt", "SUN-MCCD", "descend", note_descend))
  for (r in rows) {
    ds <- r[[1]]; meth <- r[[2]]; dir <- r[[3]]; note <- r[[4]]
    keys <- c(dataset = ds, method = meth, direction = dir)
    if (isTRUE(has_result(OUT_CSV, keys))) { cat(sprintf("[already present, skip] %s\n", paste(keys, collapse = " x "))); next }
    nd <- info[[ds]]
    append_result(OUT_CSV, list(
      dataset = ds, method = meth, direction = dir, n = nd[1], d = nd[2],
      TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
      n_flagged = NA_integer_, n_clusters = NA_integer_, flagged_idx = "",
      elapsed_sec = NA_real_, status = "not_run", note = note,
      nnradi_calls = NA_integer_, nnradi_method_received = NA_character_,
      radii_sig = NA_character_, radii_mean = NA_real_, timestamp = format(Sys.time())))
    cat(sprintf("[appended not_run] %s\n", paste(keys, collapse = " x ")))
  }
}

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------
run_grid <- function() {
  cat(sprintf("43_wp2a: datasets   = %s\n", paste(DATASETS, collapse = ", ")))
  cat(sprintf("43_wp2a: methods    = %s\n", paste(METHODS, collapse = ", ")))
  cat(sprintf("43_wp2a: directions = %s\n", paste(DIRECTIONS, collapse = ", ")))
  cat(sprintf("43_wp2a: time_budget = %gs   out = %s\n", TIME_BUDGET, OUT_CSV))
  t_start <- Sys.time()

  for (ds in DATASETS) {
    dat <- tryCatch(load_real_dataset(ds), error = function(e) NULL)
    if (is.null(dat)) { cat(sprintf("[dataset load FAILED] %s\n", ds)); next }

    for (meth in METHODS) {
      for (dir in DIRECTIONS) {
        keys <- c(dataset = ds, method = meth, direction = dir)
        if (isTRUE(has_result(OUT_CSV, keys))) {
          cat(sprintf("[skip] %s x %s x %s\n", ds, meth, dir)); next
        }

        elapsed_so_far <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        if (elapsed_so_far > TIME_BUDGET) {
          cat(sprintf("[time budget %gs reached after %.1fs -- stopping; re-invoke to continue] next up: %s x %s x %s\n",
                      TIME_BUDGET, elapsed_so_far, ds, meth, dir))
          return(invisible("budget_stop"))
        }

        install_nnradi_tap()
        t0 <- Sys.time()
        out <- tryCatch({
          setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
          res <- if (meth %in% MIN_CLS_METHODS) {
            METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, method = dir, min.cls = S_MIN)
          } else {
            METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, method = dir)
          }
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
               score = res$score, cluster = res$cluster, status = "ok", note = NA_character_)
        }, error = function(e) {
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          msg <- conditionMessage(e)
          list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
               score = NULL, cluster = NULL,
               status = if (grepl("elapsed time limit|reached elapsed time limit", msg)) "timeout" else "error",
               note = msg)
        })
        wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        flagged    <- if (is.null(out$score)) integer(0) else which(out$score > 0.5)
        n_clusters <- if (is.null(out$cluster)) NA_integer_ else length(unique(stats::na.omit(out$cluster)))

        calls <- NNRADI_TAP$calls
        n_calls <- length(calls)
        method_received <- if (n_calls >= 1) calls[[n_calls]]$method_received else NA_character_
        Rvec <- if (n_calls >= 1) calls[[n_calls]]$R else NULL

        append_result(OUT_CSV, list(
          dataset = ds, method = meth, direction = dir, n = dat$n, d = dat$d,
          TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]), BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
          n_flagged  = if (is.null(out$score)) NA_integer_ else length(flagged),
          n_clusters = n_clusters,
          flagged_idx = collapse_idx(flagged),
          elapsed_sec = wall,
          status = out$status, note = out$note,
          nnradi_calls = n_calls, nnradi_method_received = method_received,
          radii_sig  = radii_sig(Rvec),
          radii_mean = if (is.null(Rvec)) NA_real_ else mean(Rvec),
          timestamp  = format(Sys.time())
        ))

        cat(sprintf("  %-14s x %-9s x %-8s  BA=%s F2=%s nflag=%-4s ncls=%-3s recv=%-8s %-7s %.1fs\n",
                    ds, meth, dir,
                    if (is.na(out$m[["BA"]])) "  NA " else sprintf("%.3f", out$m[["BA"]]),
                    if (is.na(out$m[["F2"]])) "  NA " else sprintf("%.3f", out$m[["F2"]]),
                    if (is.null(out$score)) "NA" else length(flagged),
                    n_clusters, method_received, out$status, wall))
        flush.console()
      }
    }
  }
  cat("43_wp2a: grid done (all requested cells covered)\n")
  invisible("done")
}

# ---------------------------------------------------------------------------
# Analysis: per dataset x method, compare ascend vs descend
# ---------------------------------------------------------------------------
build_summary <- function() {
  if (!file.exists(OUT_CSV)) { cat("43_wp2a: no results yet\n"); return(invisible(NULL)) }
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE, colClasses = c(flagged_idx = "character",
                                                                    radii_sig = "character"))
  df$flagged_idx[is.na(df$flagged_idx)] <- ""

  rows <- list()
  for (ds in unique(df$dataset)) {
    for (meth in unique(df$method[df$dataset == ds])) {
      sub <- df[df$dataset == ds & df$method == meth, ]
      a <- sub[sub$direction == "ascend", ]
      b <- sub[sub$direction == "descend", ]
      if (nrow(a) == 0 || nrow(b) == 0) next
      a <- a[nrow(a), ]; b <- b[nrow(b), ]

      ok <- identical(a$status, "ok") && identical(b$status, "ok")
      a_idx <- parse_idx(a$flagged_idx); b_idx <- parse_idx(b$flagged_idx)
      inter <- length(intersect(a_idx, b_idx)); uni <- length(union(a_idx, b_idx))
      jac <- if (uni == 0) 1 else inter / uni
      identical_labels <- identical(a_idx, b_idx)
      radii_differ <- !identical(a$radii_sig, b$radii_sig)
      ncl_differ <- !identical(a$n_clusters, b$n_clusters)

      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, method = meth, n = a$n, d = a$d,
        status_ascend = a$status, status_descend = b$status,
        recv_ascend = a$nnradi_method_received, recv_descend = b$nnradi_method_received,
        radii_differ = if (ok) radii_differ else NA,
        n_clusters_ascend = a$n_clusters, n_clusters_descend = b$n_clusters,
        n_clusters_differ = if (ok) ncl_differ else NA,
        n_flagged_ascend = a$n_flagged, n_flagged_descend = b$n_flagged,
        jaccard = if (ok) jac else NA_real_,
        identical_labels = if (ok) identical_labels else NA,
        TPR_ascend = a$TPR, TPR_descend = b$TPR, dTPR = if (ok) b$TPR - a$TPR else NA_real_,
        TNR_ascend = a$TNR, TNR_descend = b$TNR, dTNR = if (ok) b$TNR - a$TNR else NA_real_,
        BA_ascend = a$BA, BA_descend = b$BA, dBA = if (ok) b$BA - a$BA else NA_real_,
        F2_ascend = a$F2, F2_descend = b$F2, dF2 = if (ok) b$F2 - a$F2 else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows) == 0) { cat("43_wp2a: no complete (ascend,descend) pairs yet\n"); return(invisible(NULL)) }
  out <- do.call(rbind, rows)
  out <- out[order(out$dataset, out$method), ]

  # ---- console report ----
  cat("\n==== ARGUMENT-REACHES-nnccd.radi GUARD ====\n")
  # Only judge cells where the tap actually recorded a value (status ok or
  # timeout both still populate recv_* from the last completed/attempted
  # call attempt; not_run rows never called the detector at all and have
  # NA here by construction -- NA must not be misread as a mismatch).
  checkable <- !is.na(out$recv_ascend) & !is.na(out$recv_descend)
  bad_recv <- out[checkable & ((out$recv_ascend != "ascend") | (out$recv_descend != "descend")), ]
  cat(sprintf("  cells with a recorded nnccd.radi call on both sides: %d / %d pairs (rest are not_run, no call made)\n",
              sum(checkable), nrow(out)))
  if (nrow(bad_recv) == 0) {
    cat("  every checkable cell's nnccd.radi tap received the direction that was requested (nnradi_method_received matches).\n")
  } else {
    cat("  *** MISMATCH: requested direction did not reach nnccd.radi in some cells ***\n")
    print(bad_recv[, c("dataset", "method", "recv_ascend", "recv_descend")])
  }
  n_radii_differ <- sum(out$radii_differ, na.rm = TRUE)
  cat(sprintf("  cells where ascend/descend produced DIFFERENT radii_sig: %d / %d ok pairs\n",
              n_radii_differ, sum(!is.na(out$radii_differ))))
  if (n_radii_differ == 0) {
    cat("  *** GUARD FAILED: radii never differ -- 'no effect' would be indistinguishable from a plumbing failure ***\n")
  }

  cat("\n==== LABEL / METRIC COMPARISON: ascend vs descend ====\n")
  ok_rows <- out[!is.na(out$identical_labels), ]
  cat(sprintf("  cells compared: %d\n", nrow(ok_rows)))
  cat(sprintf("  identical flagged sets: %d / %d\n", sum(ok_rows$identical_labels), nrow(ok_rows)))
  changed <- ok_rows[!ok_rows$identical_labels, ]
  cat(sprintf("  cells whose labels CHANGED: %d\n", nrow(changed)))
  if (nrow(changed) > 0) {
    for (i in seq_len(nrow(changed))) {
      r <- changed[i, ]
      cat(sprintf("    %-14s %-9s jaccard=%.4f dBA=%+.4f dF2=%+.4f nflag %d(asc) vs %d(desc)\n",
                  r$dataset, r$method, r$jaccard, r$dBA, r$dF2, r$n_flagged_ascend, r$n_flagged_descend))
    }
  }
  dBA_valid <- ok_rows$dBA[!is.na(ok_rows$dBA)]
  if (length(dBA_valid)) {
    k <- which.max(abs(ok_rows$dBA))
    cat(sprintf("  largest |dBA| = %.4f at %s / %s (ascend BA=%.4f, descend BA=%.4f)\n",
                abs(ok_rows$dBA[k]), ok_rows$dataset[k], ok_rows$method[k],
                ok_rows$BA_ascend[k], ok_rows$BA_descend[k]))
    cat(sprintf("  mean dBA = %+.5f, mean dF2 = %+.5f (descend - ascend, over %d ok pairs)\n",
                mean(ok_rows$dBA, na.rm = TRUE), mean(ok_rows$dF2, na.rm = TRUE), nrow(ok_rows)))
  }

  cat(sprintf("\n==== CLUSTER-COUNT CHANGES: %d / %d ok pairs ====\n",
              sum(ok_rows$n_clusters_differ, na.rm = TRUE), sum(!is.na(ok_rows$n_clusters_differ))))
  cc <- ok_rows[!is.na(ok_rows$n_clusters_differ) & ok_rows$n_clusters_differ, ]
  if (nrow(cc)) {
    for (i in seq_len(nrow(cc))) {
      cat(sprintf("    %-14s %-9s n_clusters ascend=%d descend=%d\n",
                  cc$dataset[i], cc$method[i], cc$n_clusters_ascend[i], cc$n_clusters_descend[i]))
    }
  } else cat("  none: cluster count never moved between directions.\n")

  bad <- out[out$status_ascend != "ok" | out$status_descend != "ok", ]
  cat(sprintf("\n==== NON-OK PAIRS: %d ====\n", nrow(bad)))
  if (nrow(bad)) print(bad[, c("dataset", "method", "status_ascend", "status_descend")])

  not_run <- setdiff(DATASETS, unique(df$dataset))
  if (length(not_run)) cat(sprintf("\nNOT YET RUN (time budget): %s\n", paste(not_run, collapse = ", ")))

  # ---- findings.md ----
  hdr <- paste(readLines(here::here("revision_experiments/43_wp2a_direction_sweep.R"), n = 87)[3:87], collapse = "\n")
  md <- c(
    "# WP2(a) radius-search direction sweep (ascend vs descend) -- findings",
    "",
    "Also serves reviewer point R3.10a (ablation toggle).",
    "",
    "## Step 1: what the toggle does and its reach",
    "",
    "```",
    hdr,
    "```",
    "",
    "## Step 2/3: measured comparison (ascend vs descend), per dataset x method",
    "",
    sprintf("Cells run: %d rows in `wp2a_direction_sweep.csv`. Complete (ascend,descend) pairs analysed: %d.",
            nrow(df), nrow(out)),
    sprintf("Datasets not yet run (time budget exhausted): %s.",
            if (length(not_run)) paste(not_run, collapse = ", ") else "none -- all requested datasets covered."),
    "",
    "### Argument-reaches-nnccd.radi guard",
    "",
    sprintf("- Every cell's `nnradi_method_received` matched the requested direction: %s.",
            if (nrow(bad_recv) == 0) sprintf("YES (0 mismatches across %d checkable pairs; %d not_run pairs made no call and are excluded)", sum(checkable), nrow(out) - sum(checkable))
            else sprintf("NO -- %d mismatches out of %d checkable pairs, see CSV", nrow(bad_recv), sum(checkable))),
    sprintf("- Cells where ascend and descend produced a different `radii_sig` (the actual radius vector nnccd.radi returned): %d / %d ok pairs.",
            n_radii_differ, sum(!is.na(out$radii_differ))),
    if (n_radii_differ == 0)
      "- **GUARD FAILED**: radii are byte-identical in every cell. Any 'no effect on labels' claim below would be indistinguishable from a plumbing failure and should NOT be reported as a scientific null result without further investigation."
    else
      "- Guard passed: the two directions are demonstrably reaching different code paths inside nnccd.radi and producing different radii, so any label/metric equality found below is a genuine null result, not a wiring failure.",
    "",
    "### Label and metric comparison",
    "",
    sprintf("- Cells compared: %d.", nrow(ok_rows)),
    sprintf("- Identical flagged sets: %d / %d.", sum(ok_rows$identical_labels), nrow(ok_rows)),
    sprintf("- Cells whose labels changed: %d.", nrow(changed)),
    if (length(dBA_valid)) sprintf("- Largest |dBA| = %.4f at %s / %s (ascend BA=%.4f, descend BA=%.4f).",
                                    abs(ok_rows$dBA[k]), ok_rows$dataset[k], ok_rows$method[k],
                                    ok_rows$BA_ascend[k], ok_rows$BA_descend[k]) else "",
    if (length(dBA_valid)) sprintf("- Mean dBA (descend-ascend) = %+.5f, mean dF2 = %+.5f over %d ok pairs.",
                                    mean(ok_rows$dBA, na.rm = TRUE), mean(ok_rows$dF2, na.rm = TRUE), nrow(ok_rows)) else "",
    sprintf("- Cluster-count changes: %d / %d ok pairs.",
            sum(ok_rows$n_clusters_differ, na.rm = TRUE), sum(!is.na(ok_rows$n_clusters_differ))),
    "",
    "### Per-cell table",
    "",
    "| dataset | method | n | d | radii_differ | ncls_asc | ncls_desc | nflag_asc | nflag_desc | jaccard | identical | dBA | dF2 |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    apply(out, 1, function(r) sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                                       r["dataset"], r["method"], r["n"], r["d"], r["radii_differ"],
                                       r["n_clusters_ascend"], r["n_clusters_descend"],
                                       r["n_flagged_ascend"], r["n_flagged_descend"],
                                       ifelse(is.na(r["jaccard"]), "NA", sprintf("%.4f", as.numeric(r["jaccard"]))),
                                       r["identical_labels"],
                                       ifelse(is.na(r["dBA"]), "NA", sprintf("%+.4f", as.numeric(r["dBA"]))),
                                       ifelse(is.na(r["dF2"]), "NA", sprintf("%+.4f", as.numeric(r["dF2"]))))),
    "",
    "### Headline verdict",
    "",
    if (n_radii_differ == 0) "Cannot state a verdict: the direction argument never demonstrably changed the radius search (see guard above)." else
    if (nrow(changed) == 0) "Ascend and descend produce IDENTICAL flagged sets (and hence identical TPR/TNR/BA/F2) in every cell tested, even though the radii themselves differ internally -- the discrete outlier decision is not sensitive to search direction on these real data sets, despite the toggle demonstrably reaching the radius search." else
    sprintf("Ascend and descend disagree in %d of %d cells tested. %s", nrow(changed), nrow(ok_rows),
            if (!is.na(mean(ok_rows$dBA, na.rm = TRUE)) && abs(mean(ok_rows$dBA, na.rm = TRUE)) > 0.01)
              sprintf("Mean dBA (descend-ascend) = %+.4f suggests a systematic direction, not a wash.", mean(ok_rows$dBA, na.rm = TRUE))
            else "The differences are small and not systematically in one direction -- effectively a wash."),
    ""
  )
  dir.create(dirname(FINDINGS_MD), recursive = TRUE, showWarnings = FALSE)
  writeLines(md, FINDINGS_MD)
  cat(sprintf("\n43_wp2a: wrote %s\n", FINDINGS_MD))
  invisible(out)
}

if (MARK_NOT_RUN) {
  mark_not_run()
  build_summary()
} else if (SUMMARY_ONLY) {
  build_summary()
} else {
  run_grid()
  build_summary()
}
