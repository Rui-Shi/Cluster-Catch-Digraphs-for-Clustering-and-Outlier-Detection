#!/usr/bin/env Rscript
# 13_glass_os_regen.R -- regenerate the four OS-method rows for glass only.
#
# data/outlier_detection/RealData_Collection.R sorts glass by column 9 (a
# feature; the label is column 10), so after that sort the 9 true glass
# outliers sit at rows 182-190 of 213 rather than at the end. The count_*
# helpers in R/general_functions/count.R score POSITIONALLY (they assume the
# last n0 rows are the outliers), so any glass number obtained by counting
# positionally is measured against the wrong ground truth.
#
# harness.R's evaluate() fixes this by jointly reordering (Y, score) before
# counting. This script uses evaluate() (not the raw count_* path) for the
# four OS methods -- RKCCD-OOS, RKCCD-IOS, UNCCD-OOS, UNCCD-IOS -- on glass,
# and writes the corrected TPR/TNR/BA/F2 to results/tr2/glass_os_regen.csv.
#
# RealData_Collection.R itself is NOT modified (it is shared with another
# project).

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))

set.seed(1)  # defensive; the four OS scorers (greedy2 domination) are
             # deterministic, but this is cheap insurance.

OUT_CSV <- here::here("revision_experiments/results/tr2/glass_os_regen.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT_CSV)) file.remove(OUT_CSV)  # fresh run, no stale rows

METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")

dat <- load_real_dataset("glass")
cat(sprintf("glass: n=%d, d=%d, n0=%d\n", dat$n, dat$d, sum(dat$Y == 0)))

for (meth in METHODS) {
  t0 <- Sys.time()
  res <- METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
  m <- evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]])
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  row <- list(dataset = "glass", method = meth,
              TPR = unname(m[["TPR"]]), TNR = unname(m[["TNR"]]),
              BA  = unname(m[["BA"]]),  F2  = unname(m[["F2"]]))
  append_result(OUT_CSV, row)

  cat(sprintf("  %-9s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f  (%.1fs)\n",
              meth, m[["TPR"]], m[["TNR"]], m[["BA"]], m[["F2"]], wall))
}
cat("13_glass_os_regen.R: done\n")
