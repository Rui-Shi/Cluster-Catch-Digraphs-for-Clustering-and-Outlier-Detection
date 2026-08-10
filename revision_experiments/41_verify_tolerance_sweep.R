#!/usr/bin/env Rscript
# 41_verify_tolerance_sweep.R -- independent verification of 40_wp2a's result.
#
# 40_wp2a_tolerance_sweep.R reports that the outlier labels are NOT invariant to
# the KS-CCD cascade tolerance. That is a negative result, so it gets checked
# rather than taken on trust. Three checks, in increasing order of strength.
#
# (A) RECOMPUTE the stability statistics from the raw per-cell CSV's
#     flagged_idx column, ignoring wp2a_tolerance_summary.csv entirely. If 40's
#     summary builder has a bug, A and the summary disagree.
#
# (B) ANCHOR the tol=0.01 rows against results/tr1/final_comparison.csv -- the
#     regenerated real-data table the manuscript now reports. 40 runs every
#     cell through a SHADOWED connected.ksccd.m; the production pipeline runs
#     the pristine one. At tol = 0.01 the two must agree, because 0.01 is the
#     literal in the production body. This is the check 40 never did: its own
#     guard only proves delta MOVES between tolerances, which a faithful copy
#     and a subtly mistranscribed copy would both pass. A uniform shift caused
#     by a bad copy shows up here and nowhere else.
#
# (C) LIVE CONTROL: call the four detectors with NO override installed at all,
#     so the pristine connected.ksccd.m runs, and compare the flagged sets to
#     40's tol=0.01 rows. (B) is free and covers all 8 datasets; (C) costs a
#     few minutes and proves it against the actual function object rather than
#     against a stored table.
#
# Read-only with respect to every result file. Writes nothing.

suppressMessages(library(here))

SWEEP <- here::here("revision_experiments/results/tr1/wp2a_tolerance_sweep.csv")
SUMM  <- here::here("revision_experiments/results/tr1/wp2a_tolerance_summary.csv")
FINAL <- here::here("revision_experiments/results/tr1/final_comparison.csv")
REF   <- 0.01

args      <- commandArgs(trailingOnly = TRUE)
RUN_LIVE  <- any(args == "--live")     # (C) is opt-in; it re-runs detectors

df <- read.csv(SWEEP, stringsAsFactors = FALSE,
               colClasses = c(delta_values = "character", flagged_idx = "character"))
df$flagged_idx[is.na(df$flagged_idx)]   <- ""
df$delta_values[is.na(df$delta_values)] <- ""

pidx <- function(s) if (!nzchar(s)) integer(0) else sort(as.integer(strsplit(s, ";")[[1]]))

cat(sprintf("rows: %d | datasets: %d | methods: %d | tols: %s\n",
            nrow(df), length(unique(df$dataset)), length(unique(df$method)),
            paste(sort(unique(df$tol)), collapse = ", ")))
cat(sprintf("status values: %s\n\n", paste(unique(df$status), collapse = ", ")))

# ---- (A) recompute stability from flagged_idx -----------------------------
cat("=== (A) recomputed from raw flagged_idx ===\n")
pairs <- unique(df[, c("dataset", "method")])
rec <- list()
for (i in seq_len(nrow(pairs))) {
  ds <- pairs$dataset[i]; me <- pairs$method[i]
  sub <- df[df$dataset == ds & df$method == me, ]
  r <- sub[abs(sub$tol - REF) < 1e-12, ]
  if (!nrow(r)) next
  ridx <- pidx(r$flagged_idx[1])
  for (j in seq_len(nrow(sub))) {
    if (abs(sub$tol[j] - REF) < 1e-12) next
    cidx <- pidx(sub$flagged_idx[j])
    u <- length(union(cidx, ridx))
    rec[[length(rec) + 1]] <- data.frame(
      dataset = ds, method = me, tol = sub$tol[j],
      same    = identical(cidx, ridx),
      jaccard = if (u == 0) 1 else length(intersect(cidx, ridx)) / u,
      flips   = length(setdiff(cidx, ridx)) + length(setdiff(ridx, cidx)),
      dBA     = abs(sub$BA[j] - r$BA[1]),
      dF2     = abs(sub$F2[j] - r$F2[1]),
      stringsAsFactors = FALSE)
  }
}
rec <- do.call(rbind, rec)
cat(sprintf("  non-reference cells      : %d\n", nrow(rec)))
cat(sprintf("  identical to tol=0.01    : %d  (moved: %d)\n",
            sum(rec$same), sum(!rec$same)))
for (t in sort(unique(rec$tol))) {
  s <- rec[rec$tol == t, ]
  cat(sprintf("    tol=%-6g identical %2d/%2d   max|dBA|=%.4f  max|dF2|=%.4f\n",
              t, sum(s$same), nrow(s), max(s$dBA, na.rm = TRUE), max(s$dF2, na.rm = TRUE)))
}
cat(sprintf("  mean|dBA| over movers    : %.4f\n", mean(rec$dBA[!rec$same], na.rm = TRUE)))
k <- which.max(rec$dBA)
cat(sprintf("  largest |dBA| overall    : %.4f  (%s / %s / tol=%g)\n",
            rec$dBA[k], rec$dataset[k], rec$method[k], rec$tol[k]))

# does tol=1 give the same labels as tol=0.1, as 40 claims?
n_same <- 0L; n_chk <- 0L
for (i in seq_len(nrow(pairs))) {
  sub <- df[df$dataset == pairs$dataset[i] & df$method == pairs$method[i], ]
  a <- sub[abs(sub$tol - 1)   < 1e-12, ]; b <- sub[abs(sub$tol - 0.1) < 1e-12, ]
  if (!nrow(a) || !nrow(b)) next
  n_chk <- n_chk + 1L
  if (identical(pidx(a$flagged_idx[1]), pidx(b$flagged_idx[1]))) n_same <- n_same + 1L
}
cat(sprintf("  tol=1 vs tol=0.1 identical labels: %d / %d\n", n_same, n_chk))

# cross-check against 40's own summary file
if (file.exists(SUMM)) {
  sm <- read.csv(SUMM, stringsAsFactors = FALSE)
  sm <- sm[abs(sm$tol - REF) > 1e-12, ]
  cat(sprintf("  40's summary says identical: %d / %d  -> %s\n",
              sum(sm$identical_labels, na.rm = TRUE), nrow(sm),
              if (sum(sm$identical_labels, na.rm = TRUE) == sum(rec$same) &&
                  nrow(sm) == nrow(rec)) "AGREES" else "DISAGREES"))
}

# ---- (B) anchor tol=0.01 to the regenerated production table --------------
cat("\n=== (B) tol=0.01 vs final_comparison.csv (production pipeline) ===\n")
if (!file.exists(FINAL)) {
  cat("  final_comparison.csv not found -- skipped\n")
} else {
  fc <- read.csv(FINAL, stringsAsFactors = FALSE)
  ref <- df[abs(df$tol - REF) < 1e-12, ]
  bad <- 0L; chk <- 0L
  for (i in seq_len(nrow(ref))) {
    r  <- ref[i, ]
    fr <- fc[fc$dataset == r$dataset & fc$method == r$method, ]
    if (!nrow(fr)) { cat(sprintf("  %-10s %-9s  NOT IN final_comparison\n",
                                 r$dataset, r$method)); next }
    chk <- chk + 1L
    d <- c(TPR = abs(r$TPR - fr$TPR[1]), TNR = abs(r$TNR - fr$TNR[1]),
           BA  = abs(r$BA  - fr$BA[1]),  F2  = abs(r$F2  - fr$F2[1]))
    if (any(d > 5.1e-4, na.rm = TRUE)) {
      bad <- bad + 1L
      cat(sprintf("  MISMATCH %-10s %-9s  sweep %.3f/%.3f/%.3f/%.3f  final %.3f/%.3f/%.3f/%.3f\n",
                  r$dataset, r$method, r$TPR, r$TNR, r$BA, r$F2,
                  fr$TPR[1], fr$TNR[1], fr$BA[1], fr$F2[1]))
    }
  }
  cat(sprintf("  %d reference cells checked, %d mismatches -> %s\n", chk, bad,
              if (bad == 0) "OVERRIDE IS FAITHFUL AT tol=0.01" else "OVERRIDE IS NOT FAITHFUL"))
}

# ---- (C) live control against the pristine function -----------------------
cat("\n=== (C) live control: pristine connected.ksccd.m, no override ===\n")
if (!RUN_LIVE) {
  cat("  skipped (pass --live to run)\n")
} else {
  source(here::here("revision_experiments/harness.R"))
  source(here::here("revision_experiments/wp0_mccd_methods.R"))
  # deliberately install NOTHING: connected.ksccd.m is the shipped one.
  stopifnot("connected.ksccd.m must be the pristine 1-arg version" =
              length(formals(connected.ksccd.m)) == 1)
  S_MIN <- 0.0625
  bad <- 0L; chk <- 0L
  for (ds in c("hepatitis", "glass", "vertebral")) {
    dat <- load_real_dataset(ds)
    for (me in c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")) {
      res <- if (me %in% c("SU-MCCD", "SUN-MCCD"))
        METHOD_REGISTRY[[me]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
      else
        METHOD_REGISTRY[[me]](X = dat$X, d = dat$d, Y = dat$Y)
      live <- sort(which(res$score > 0.5))
      r <- df[df$dataset == ds & df$method == me & abs(df$tol - REF) < 1e-12, ]
      swept <- pidx(r$flagged_idx[1])
      chk <- chk + 1L
      ok <- identical(as.integer(live), as.integer(swept))
      if (!ok) bad <- bad + 1L
      cat(sprintf("  %-10s %-9s live nflag=%-4d sweep nflag=%-4d  %s\n",
                  ds, me, length(live), length(swept), if (ok) "MATCH" else "*** DIFFER ***"))
    }
  }
  cat(sprintf("  %d cells, %d differ -> %s\n", chk, bad,
              if (bad == 0) "SHADOWED PATH == PRODUCTION PATH" else "PATHS DIVERGE"))
}

cat("\ndone\n")
