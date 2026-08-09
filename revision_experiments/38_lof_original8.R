#!/usr/bin/env Rscript
# 38_lof_original8.R -- compute LOF on the original eight data sets.
#
# Why, given that LOF was deliberately not re-run:
#
# LOF's CONFIGURATION was never in doubt -- LOF.R defaults to MinPts 11..30 and
# Real_Data_LOF_1.5.R calls it with exactly those bounds and Thresh = 1.5,
# matching main text L1007. That is why the published rows were retained.
#
# What was in doubt, and is now settled, is the COUNTING. The retained rows
# came from published_realdata_truth.csv by transcription
# (34_final_summary.R:59), making LOF x original-8 the only block of the 16x9
# table not produced by this pipeline. Two independent problems have since
# surfaced in that block:
#
#   1. 37_arith_audit.R finds LOF/ecoli and LOF/vertebral report a BA that
#      does not equal (TPR + TNR)/2 -- 0.711 vs 0.715 and 0.488 vs 0.485.
#      Every other cell of the 144 satisfies the identity to rounding.
#   2. glass and ecoli are the two data sets RealData_Collection.R mis-sorts,
#      and the published counting code slices positionally. Whatever those two
#      LOF rows measure, it is not performance against the true labels.
#
# So this script changes nothing about how LOF is run -- same function, same
# bounds, same threshold -- and only changes how the output is scored: through
# harness.R's evaluate(), which reorders (Y, score) jointly before counting.
#
# Writes results/tr1/lof_original8.csv. Touches no existing file.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))

ORIG <- c("hepatitis", "glass", "vertebral", "ecoli",
          "stamps", "vowels", "waveform", "wilt")

OUT <- here::here("revision_experiments/results/tr1/lof_original8.csv")

TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
pub <- function(ds, mt) {
  v <- TRUTH$value[TRUTH$dataset == ds & TRUTH$method == "LOF" & TRUTH$metric == mt]
  if (length(v) == 0) NA_real_ else v[1]
}

cat(sprintf("%-11s %5s %4s %4s | %-27s | %-27s | %s\n",
            "dataset", "n", "d", "n0",
            "recomputed TPR/TNR/BA/F2", "published TPR/TNR/BA/F2", "max|d|"))
cat(strrep("-", 108), "\n")

for (ds in ORIG) {
  dat <- load_real_dataset(ds)
  n0  <- sum(dat$Y == 0)

  t0 <- Sys.time()
  score <- LOF(as.matrix(dat$X), L_MinPts = 11, U_MinPts = 30)
  m <- evaluate(dat$Y, score, 1.5)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  p <- sapply(c("TPR", "TNR", "BA", "F2"), function(x) pub(ds, x))
  mx <- suppressWarnings(max(abs(round(unname(m), 3) - unname(p)), na.rm = TRUE))

  cat(sprintf("%-11s %5d %4d %4d | %6.3f %6.3f %6.3f %6.3f | %6.3f %6.3f %6.3f %6.3f | %.3f\n",
              ds, dat$n, dat$d, n0,
              m[["TPR"]], m[["TNR"]], m[["BA"]], m[["F2"]],
              p[["TPR"]], p[["TNR"]], p[["BA"]], p[["F2"]], mx))

  append_result(OUT, list(
    dataset = ds, method = "LOF",
    params = "MinPts=11..30, max over k, threshold=1.5",
    n = dat$n, d = dat$d, n0 = n0,
    TPR = unname(m[["TPR"]]), TNR = unname(m[["TNR"]]),
    BA = unname(m[["BA"]]), F2 = unname(m[["F2"]]),
    published_TPR = unname(p[["TPR"]]), published_TNR = unname(p[["TNR"]]),
    published_BA = unname(p[["BA"]]), published_F2 = unname(p[["F2"]]),
    max_abs_diff = mx, match_3dp = isTRUE(mx < 0.0005),
    t_total = wall, status = "ok", note = NA_character_,
    timestamp = format(Sys.time())))
}

cat("\nWritten:", OUT, "\n")
