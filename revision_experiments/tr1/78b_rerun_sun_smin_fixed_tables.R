#!/usr/bin/env Rscript
# 78b_rerun_sun_smin_fixed_tables.R -- rerun the SUN-MCCD S_min sensitivity cells
# with the repaired 0.1% NND quantile tables (71_regen_nnd_alpha001.R,
# 75_install_regen999.R, 77_purge_wrong_tables.R).
#
# Why: an outside check (2026-09-22) found that the SUN-MCCD columns of
#   results/tr1/wp2c_smin_zero_real.csv (60_*.R; supplement S_min=0 table)
#   and results/tr1/wp2a_smin_bigfour.csv (51_*.R; supplement plateau table)
# predate the alpha repair for hepatitis, vowels and waveform (and PenDigits
# for all but S_min = 0.05), so they disagree with the published per-data-set
# table. SU-MCCD is unaffected (RK tables were never mislabelled).
#
# Part A: all 11 small sets x S_min in {0, 0.03, 0.05, 0.0625}  (same grid as 60)
# Part B: waveform, PenDigits x S_min in {0.01, 0.05, 0.0625, 0.15, 0.30} (grid of 51)
#
# Sanity anchor: every S_min = 0.05 cell must reproduce the published SUN-MCCD
# row (e.g. hepatitis BA 0.705, vowels 0.658, waveform 0.641, PenDigits 0.5395).
#
# Appends one row per cell to a NEW file (G: drops off the bus under load; never
# buffer). Existing result files are read-only.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

OUT <- here::here("revision_experiments/results/tr1/wp2c_sun_smin_rerun_fixed.csv")
args <- commandArgs(trailingOnly = TRUE)
PART <- if (length(args) >= 1) args[1] else "A"

if (PART == "A") {
  SETS <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
            "stamps", "WDBC", "pima", "shuffle", "vowels")
  SMIN <- c(0, 0.03, 0.05, 0.0625)
} else {
  SETS <- c("waveform", "PenDigits")
  SMIN <- c(0.01, 0.05, 0.0625, 0.15, 0.30)
}

done <- if (file.exists(OUT)) read.csv(OUT, stringsAsFactors = FALSE) else NULL
for (ds in SETS) {
  dat <- load_real_dataset(ds)
  for (sm in SMIN) {
    if (!is.null(done) && any(done$dataset == ds & abs(done$s_min - sm) < 1e-12 &
                              done$status == "ok")) next
    t0 <- Sys.time()
    out <- tryCatch({
      res <- METHOD_REGISTRY[["SUN-MCCD"]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = sm)
      m <- evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[["SUN-MCCD"]])
      list(m = m, res = res, status = "ok", note = NA_character_)
    }, error = function(e) list(m = setNames(rep(NA_real_, 4), c("TPR","TNR","BA","F2")),
                                res = NULL, status = "error", note = conditionMessage(e)))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ncls <- if (is.null(out$res)) NA_integer_ else
      length(unique(out$res$cluster[!is.na(out$res$cluster)]))
    row <- data.frame(
      part = PART, dataset = ds, method = "SUN-MCCD", s_min = sm,
      n = nrow(dat$X), d = dat$d, k_threshold = round(sm * nrow(dat$X)),
      TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
      BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
      n_clusters = ncls,
      quant_label = if (is.null(out$res)) NA_character_ else out$res$quant_label,
      elapsed_sec = el, status = out$status, note = out$note,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE)
    write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
                col.names = !file.exists(OUT), qmethod = "double")
    cat(sprintf("CELL_DONE %s s_min=%.4f BA=%.4f F2=%.4f ncls=%s q=%s sec=%.1f\n",
                ds, sm, ifelse(is.na(row$BA), -1, row$BA), ifelse(is.na(row$F2), -1, row$F2),
                ncls, row$quant_label, el))
    flush.console()
  }
}
cat("ALL_CELLS_COMPLETE\n")
