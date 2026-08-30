#!/usr/bin/env Rscript
# 15_os_repro_audit.R -- do the four OS rows of the manuscript's real-data
# tables reproduce under the current harness?
#
# Motivation (2026-08-16): the WP3 sensitivity sweep is to swap Arrhythmia for
# a real dataset that RK-CCD can score. Adding WDBC tripped the new
# reproduction gate in 10_wp3_real.R: none of its four published OS values came
# back. Re-checking Thyroid showed the two RK-CCD values do not reproduce
# either, while both UN-CCD values do. WBC reproduces on all four.
#
# This script settles how far that goes, so the replacement dataset can be
# chosen on evidence rather than on the two data points we happen to have.
#
# For each dataset: score with METHOD_REGISTRY, evaluate() at the published
# cutoff of 2, and compare against the manuscript's Tables
# Real_Data_Result_OS1/OS2 (values passed in via tr1/published_realdata_truth.csv
# if present, else the hard-coded PUBLISHED list below, transcribed from
# CCDwScores.tex on 2026-08-16).
#
# wilt (n=4819) is excluded by default: RK-CCD on it needs the dedicated
# launchers (revision_experiments/tr1/regen_wilt_rk_launcher.ps1). Pass it explicitly to include it.
#
# CLI:  Rscript revision_experiments/tr2/15_os_repro_audit.R [dataset ...]

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments", "shared", "harness.R")))

OUT_CSV <- here::here("revision_experiments/results/tr2/os_repro_audit.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
TOL <- 0.005 + 1e-9  # manuscript rounds to 3 decimals

# Manuscript values: dataset -> method -> c(TPR, TNR, BA, F2).
PUBLISHED <- list(
  hepatitis = list("RKCCD-OOS" = c(0.000, 0.866, 0.433, 0.000),
                   "RKCCD-IOS" = c(0.142, 0.925, 0.534, 0.147),
                   "UNCCD-OOS" = c(0.286, 0.866, 0.576, 0.256),
                   "UNCCD-IOS" = c(0.571, 0.925, 0.748, 0.541)),
  lymphography = list("RKCCD-OOS" = c(0.500, 0.725, 0.613, 0.227),
                      "RKCCD-IOS" = c(0.833, 0.789, 0.811, 0.424),
                      "UNCCD-OOS" = c(0.833, 0.697, 0.765, 0.347),
                      "UNCCD-IOS" = c(0.833, 0.873, 0.853, 0.532)),
  glass = list("RKCCD-OOS" = c(0.333, 0.789, 0.561, 0.183),
               "RKCCD-IOS" = c(0.444, 0.770, 0.607, 0.230),
               "UNCCD-OOS" = c(0.444, 0.838, 0.641, 0.274),
               "UNCCD-IOS" = c(0.222, 0.931, 0.577, 0.192)),
  WBC = list("RKCCD-OOS" = c(0.200, 0.850, 0.525, 0.135),
             "RKCCD-IOS" = c(1.000, 0.779, 0.890, 0.515),
             "UNCCD-OOS" = c(0.200, 0.897, 0.548, 0.156),
             "UNCCD-IOS" = c(1.000, 0.807, 0.904, 0.549)),
  vertebral = list("RKCCD-OOS" = c(0.133, 0.905, 0.519, 0.139),
                   "RKCCD-IOS" = c(0.067, 0.843, 0.455, 0.064),
                   "UNCCD-OOS" = c(0.100, 0.905, 0.502, 0.104),
                   "UNCCD-IOS" = c(0.100, 0.810, 0.455, 0.092)),
  stamps = list("RKCCD-OOS" = c(0.258, 0.864, 0.561, 0.230),
                "RKCCD-IOS" = c(0.226, 0.838, 0.532, 0.193),
                "UNCCD-OOS" = c(0.194, 0.887, 0.540, 0.181),
                "UNCCD-IOS" = c(0.258, 0.838, 0.548, 0.220)),
  WDBC = list("RKCCD-OOS" = c(0.700, 0.807, 0.753, 0.302),
              "RKCCD-IOS" = c(0.700, 0.913, 0.807, 0.449),
              "UNCCD-OOS" = c(0.700, 0.874, 0.787, 0.380),
              "UNCCD-IOS" = c(0.300, 0.913, 0.607, 0.203)),
  vowels = list("RKCCD-OOS" = c(0.326, 0.900, 0.613, 0.221),
                "RKCCD-IOS" = c(0.783, 0.898, 0.840, 0.496),
                "UNCCD-OOS" = c(0.413, 0.927, 0.670, 0.310),
                "UNCCD-IOS" = c(0.848, 0.903, 0.875, 0.542)),
  thyroid = list("RKCCD-OOS" = c(0.280, 0.885, 0.582, 0.161),
                 "RKCCD-IOS" = c(0.828, 0.842, 0.835, 0.381),
                 "UNCCD-OOS" = c(0.247, 0.909, 0.578, 0.160),
                 "UNCCD-IOS" = c(0.989, 0.827, 0.908, 0.425)),
  wilt = list("RKCCD-OOS" = c(0.105, 0.912, 0.509, 0.092),
              "RKCCD-IOS" = c(0.304, 0.810, 0.557, 0.197),
              "UNCCD-OOS" = c(0.206, 0.911, 0.559, 0.178),
              "UNCCD-IOS" = c(0.288, 0.832, 0.560, 0.198))
)

# R object name in RealData_Collection.R (differs from the table label for
# lymphography, which the manuscript prints as "lymph").
ROBJ <- c(hepatitis = "hepatitis", lymphography = "lymphography",
          glass = "glass", WBC = "WBC", vertebral = "vertebral",
          stamps = "stamps", WDBC = "WDBC", vowels = "vowels",
          thyroid = "thyroid", wilt = "wilt")

DEFAULT_SETS <- c("hepatitis", "lymphography", "glass", "vertebral",
                  "stamps", "WDBC", "vowels")

args <- commandArgs(trailingOnly = TRUE)
SETS <- if (length(args) > 0) args else DEFAULT_SETS
stopifnot(all(SETS %in% names(PUBLISHED)))

rows <- list()
for (ds in SETS) {
  dat <- load_real_dataset(ROBJ[[ds]])
  n0 <- sum(dat$Y == 0)
  cat(sprintf("\n== %s: n=%d, d=%d, n0=%d  (RK quant %s, NN quant %s)\n",
              ds, dat$n, dat$d, n0, rk_quant_for_d(dat$d), nn_quant_for_d(dat$d)))

  for (meth in METHODS) {
    t0 <- Sys.time()
    got <- tryCatch({
      res <- METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
      evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]])
    }, error = function(e) {
      cat(sprintf("  %-10s ERROR: %s\n", meth, conditionMessage(e)))
      NULL
    })
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(got)) next

    g <- c(got[["TPR"]], got[["TNR"]], got[["BA"]], got[["F2"]])
    p <- PUBLISHED[[ds]][[meth]]
    ok <- all(abs(g - p) <= TOL)
    cat(sprintf("  %-10s got %.3f/%.3f/%.3f/%.3f | pub %.3f/%.3f/%.3f/%.3f | %-8s (%.1fs)\n",
                meth, g[1], g[2], g[3], g[4], p[1], p[2], p[3], p[4],
                if (ok) "MATCH" else "MISMATCH", wall))

    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, n = dat$n, d = dat$d, n0 = n0, method = meth,
      TPR = g[1], TNR = g[2], BA = g[3], F2 = g[4],
      pub_TPR = p[1], pub_TNR = p[2], pub_BA = p[3], pub_F2 = p[4],
      match = ok, secs = round(wall, 1), stringsAsFactors = FALSE)
  }
}

out <- do.call(rbind, rows)
write.csv(out, OUT_CSV, row.names = FALSE)

cat("\n---- summary: matching OS rows per dataset (of 4) ----\n")
for (ds in unique(out$dataset)) {
  s <- out[out$dataset == ds, ]
  cat(sprintf("  %-13s d=%-3d %d/%d match%s\n", ds, s$d[1], sum(s$match), nrow(s),
              if (all(s$match)) "   <- fully reproduces" else ""))
}
cat(sprintf("\nWrote %s\n15_os_repro_audit.R done.\n", OUT_CSV))
