#!/usr/bin/env Rscript
# 73_niter_sensitivity.R -- how much does the 0.1% table move if niter is cut?
#
# d=16 and d=21 at niter=10000 are roughly 4 and 6 days of wall clock, against
# a 2026-08-24 deadline. Halving niter halves both. The question is what that
# costs statistically, and it should be measured rather than argued: at
# niter=10000 the 0.1% quantile is the 10th order statistic, at 5000 the 5th,
# at 2500 the 2nd or 3rd -- and the original defect this whole work package
# repairs came from a table whose "0.1%" was effectively a sample minimum.
#
# The d=12 draws are already on disk at the full 10000 iterations. Treat the
# 10000-iteration table as truth, draw random subsets of the iterations at each
# candidate niter, reduce each subset at 0.999, and report how far the result
# lands from truth. The comparison is scale-free (median absolute relative
# difference over subsample sizes) so it transfers to other dimensions, and it
# is on the same footing as 72's T1 acceptance threshold of 2%.
#
# Read-only apart from its own CSV.

suppressMessages(library(here))
source(here::here("revision_experiments", "01i_nn_multiquant_table.R"))

CHK   <- here::here("R/NN-test_quantile_regen999/12d/chunks")
OUT   <- here::here("revision_experiments/results/tr1/wp2c_niter_sensitivity.csv")
CAND  <- c(2500L, 5000L, 7500L)
NREP  <- 5L
SEED  <- 20260811L

fs <- sort(list.files(CHK, pattern = "^draws_chunk[0-9]+\\.RData$", full.names = TRUE))
cat(sprintf("pooling %d chunk files from %s\n", length(fs), CHK))
ave <- do.call(rbind, lapply(fs, function(p) {
  e <- new.env(); load(p, envir = e); get("NN.dist.ave.mat", envir = e)
}))
cat(sprintf("pooled draws: %d x %d\n", nrow(ave), ncol(ave)))

#' Lower-tail quantile column by column -- same reduction as the generator.
red999 <- function(M) apply(M, 2, function(col) unname(quantile(col, 0.001)))

truth <- red999(ave)
K <- 2:ncol(ave)                      # index 1 is the structural zero

set.seed(SEED)
rows <- list()
for (m in CAND) {
  d <- vapply(seq_len(NREP), function(r) {
    idx <- sample.int(nrow(ave), m)
    est <- red999(ave[idx, , drop = FALSE])
    median(abs(est[K] - truth[K]) / truth[K])
  }, numeric(1))
  rows[[length(rows) + 1]] <- data.frame(
    niter = m, order_stat = m * 0.001, reps = NREP,
    med_rel_diff_mean = mean(d), med_rel_diff_max = max(d),
    wall_frac = m / 10000, stringsAsFactors = FALSE)
  cat(sprintf("  niter=%5d  (0.1%% = order statistic %4.1f)  median |rel diff| vs 10000: mean %.4f  max %.4f\n",
              m, m * 0.001, mean(d), max(d)))
}
res <- do.call(rbind, rows)

cat("\nreference: 72's T1 accepted the new 1% tables at median |rel diff| < 0.02\n")
cat("           vs the shipped tables, so anything well under that is noise\n")
cat("           of the same order as the generator's own reproducibility.\n\n")
print(res, row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------------
# The subsampling figures above are BIASED OPTIMISTIC and must not be quoted on
# their own. Each subset is drawn from the same pool as the "truth" it is
# compared against, so a 7500-iteration subset shares 75% of its draws with the
# reference -- and a 0.1% quantile is determined by a handful of the smallest
# draws, which are largely the SAME handful in both. The two estimates are
# therefore strongly coupled, and the difference between them understates what
# an independent rerun at that niter would produce.
#
# Disjoint splits remove the coupling: cut the 10000 iterations into
# non-overlapping blocks of size m and compare blocks to EACH OTHER. Two
# independent m-iteration runs differ by this much, which is the quantity that
# actually matters when deciding whether to generate d=16 and d=21 at a reduced
# niter.
# ---------------------------------------------------------------------------
cat("\n=== disjoint splits: two INDEPENDENT runs of size m, compared to each other ===\n")
drows <- list()
for (m in CAND) {
  nblk <- floor(nrow(ave) / m)
  if (nblk < 2) { cat(sprintf("  niter=%5d  only %d disjoint block(s) available -- skipped\n", m, nblk)); next }
  set.seed(SEED + m)
  perm <- sample.int(nrow(ave))
  tabs <- lapply(seq_len(nblk), function(b)
    red999(ave[perm[((b - 1) * m + 1):(b * m)], , drop = FALSE]))
  prs <- combn(nblk, 2)
  dd <- apply(prs, 2, function(p)
    median(abs(tabs[[p[1]]][K] - tabs[[p[2]]][K]) / tabs[[p[2]]][K]))
  drows[[length(drows) + 1]] <- data.frame(
    niter = m, blocks = nblk, pairs = ncol(prs),
    med_rel_diff_mean = mean(dd), med_rel_diff_max = max(dd),
    stringsAsFactors = FALSE)
  cat(sprintf("  niter=%5d  %d disjoint blocks, %d pairs: median |rel diff| mean %.4f  max %.4f\n",
              m, nblk, ncol(prs), mean(dd), max(dd)))
}
dres <- do.call(rbind, drows)
cat("\nUse THESE numbers, not the subsampling ones, to judge a niter cut.\n")

write.csv(cbind(res, method = "subsample_biased"), OUT, row.names = FALSE)
if (!is.null(dres))
  write.table(cbind(dres, method = "disjoint_honest"), OUT, sep = ",",
              row.names = FALSE, col.names = FALSE, append = TRUE)
cat(sprintf("\nwrote %s\n", OUT))
