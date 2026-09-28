#!/usr/bin/env Rscript
# 52_direction_summary.R -- analysis of 43_wp2a_direction_sweep.R.
#
# 43's runner stalled before producing its own summary, so this does that work
# from the per-cell CSV. Independent of anything 43 computed beyond the raw
# rows.
#
# The question: does the NN radius-search direction (ascend vs descend) change
# anything? It is the last row of the WP2(a) sensitivity table and also a WP9
# ablation toggle, so one measurement serves R3.2/R5.2 and R3.10a.
#
# Only UN-MCCD and SUN-MCCD are in scope -- the RK detectors expose no
# direction argument. That asymmetry is itself reportable: the toggle covers
# only half the 2x2 design family.
#
# Guard against the null-result trap: "identical everywhere" is only a finding
# if the argument actually reached the radius search. 43 instrumented
# nnccd.radi to record the method string it received (nnradi_method_received)
# and a signature of the radii it produced (radii_sig). Both are checked here
# before any "no effect" is reported.
#
# Read-only. Writes one summary CSV.

suppressMessages(library(here))
SWEEP <- here::here("revision_experiments/results/tr1/wp2a_direction_sweep.csv")
OUT   <- here::here("revision_experiments/results/tr1/wp2a_direction_summary.csv")

df <- read.csv(SWEEP, stringsAsFactors = FALSE, colClasses = c(flagged_idx = "character"))
df$flagged_idx[is.na(df$flagged_idx)] <- ""

cat(sprintf("rows=%d  ok=%d  timeout=%d  not_run=%d\n", nrow(df),
            sum(df$status == "ok"), sum(df$status == "timeout"),
            sum(df$status == "not_run")))

ok <- df[df$status == "ok", ]

# ---- plumbing guard --------------------------------------------------------
cat("\n=== plumbing: did the toggle reach nnccd.radi? ===\n")
if ("nnradi_method_received" %in% names(ok)) {
  bad <- ok[!is.na(ok$nnradi_method_received) & ok$nnradi_method_received != ok$direction, ]
  cat(sprintf("  cells where received != requested: %d\n", nrow(bad)))
  cat(sprintf("  cells with a recorded value: %d / %d\n",
              sum(!is.na(ok$nnradi_method_received)), nrow(ok)))
  if (nrow(bad) > 0) print(bad[, c("dataset","method","direction","nnradi_method_received")])
} else {
  cat("  column absent -- cannot verify; treat any null result as unproven\n")
}

pidx <- function(s) if (!nzchar(s)) integer(0) else sort(as.integer(strsplit(s, ";")[[1]]))
jac  <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

pairs <- unique(ok[, c("dataset", "method")])
rows <- list()
for (i in seq_len(nrow(pairs))) {
  ds <- pairs$dataset[i]; me <- pairs$method[i]
  a <- ok[ok$dataset == ds & ok$method == me & ok$direction == "ascend", ]
  b <- ok[ok$dataset == ds & ok$method == me & ok$direction == "descend", ]
  if (nrow(a) != 1 || nrow(b) != 1) next          # incomplete pair, skip
  ia <- pidx(a$flagged_idx[1]); ib <- pidx(b$flagged_idx[1])
  rows[[length(rows) + 1]] <- data.frame(
    dataset = ds, method = me, n = a$n[1], d = a$d[1],
    labels_identical = identical(ia, ib),
    jaccard          = jac(ia, ib),
    n_flag_asc = length(ia), n_flag_desc = length(ib),
    k_asc = a$n_clusters[1], k_desc = b$n_clusters[1],
    k_changed = a$n_clusters[1] != b$n_clusters[1],
    radii_sig_same = if ("radii_sig" %in% names(ok)) identical(a$radii_sig[1], b$radii_sig[1]) else NA,
    dTPR = b$TPR[1] - a$TPR[1], dTNR = b$TNR[1] - a$TNR[1],
    dBA  = b$BA[1]  - a$BA[1],  dF2  = b$F2[1]  - a$F2[1],
    BA_asc = a$BA[1], BA_desc = b$BA[1],
    sec_asc = a$elapsed_sec[1], sec_desc = b$elapsed_sec[1],
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
res <- res[order(-abs(res$dBA)), ]
write.csv(res, OUT, row.names = FALSE)

cat(sprintf("\n=== complete ascend/descend pairs: %d ===\n", nrow(res)))
cat(sprintf("  labels identical      : %d / %d\n", sum(res$labels_identical), nrow(res)))
cat(sprintf("  cluster count changed : %d / %d\n", sum(res$k_changed), nrow(res)))
if (!all(is.na(res$radii_sig_same)))
  cat(sprintf("  radii signature same  : %d / %d  (if identical everywhere, the toggle is inert downstream)\n",
              sum(res$radii_sig_same, na.rm = TRUE), sum(!is.na(res$radii_sig_same))))
cat(sprintf("  mean |dBA| over movers: %.4f\n",
            if (any(!res$labels_identical)) mean(abs(res$dBA[!res$labels_identical])) else 0))
cat(sprintf("  largest |dBA|         : %.4f (%s / %s)\n",
            abs(res$dBA[1]), res$dataset[1], res$method[1]))

# is one direction systematically better?
w <- sum(res$dBA > 1e-9); l <- sum(res$dBA < -1e-9); t <- sum(abs(res$dBA) <= 1e-9)
cat(sprintf("\n  descend better on BA: %d   ascend better: %d   tied: %d\n", w, l, t))
if (w + l > 0) {
  p <- stats::binom.test(w, w + l, 0.5)$p.value
  cat(sprintf("  sign test on the %d non-tied pairs: p = %.3f -> %s\n", w + l, p,
              if (p < 0.05) "one direction is systematically better" else "no systematic winner"))
}
cat(sprintf("  mean dBA (descend - ascend) = %+.4f\n", mean(res$dBA)))

cat("\n  per pair:\n")
print(res[, c("dataset","method","d","labels_identical","jaccard","k_asc","k_desc",
              "BA_asc","BA_desc","dBA")], row.names = FALSE, digits = 4)

miss <- df[df$status != "ok", c("dataset","method","direction","status")]
if (nrow(miss)) {
  cat("\n  cells not completed:\n")
  print(miss, row.names = FALSE)
}
cat(sprintf("\nwrote %s\n", OUT))
