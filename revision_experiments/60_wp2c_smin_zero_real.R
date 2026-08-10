#!/usr/bin/env Rscript
# 60_wp2c_smin_zero_real.R -- the missing cell in the WP2(c) comparison.
#
# WP2(c) evaluated `min.cls = 0` on SIMULATED data, where it never collapses a
# cluster and recovers the true cluster count 95.1% of the time -- the joint
# best of every rule tested. It was never run on REAL data, so the three-way
# choice between 0.05, 0.03 and 0 rests on an incomplete table.
#
# That matters more than a missing row usually would, because `min.cls = 0` is
# not an exotic setting: it is the SHIPPED DEFAULT of both detectors
# (SUMCCD_outlier and SUNMCCD_outlier declare `min.cls = 0` in their
# signatures). Adopting it would remove S_min from the method's parameter list
# altogether, which is the strongest available reply to AE.2's
# "over-parameterised" charge. If it also holds up on real data, it dominates
# both fixed constants on the argument even where it ties on the numbers.
#
# Restricted to the 11 data sets under n = 1500. The four largest are being
# swept concurrently by job 51 and cost 350-1600 s per cell; they are recorded
# as not_run rather than silently omitted.
#
# Read-only with respect to every existing result file.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

OUT <- here::here("revision_experiments/results/tr1/wp2c_smin_zero_real.csv")

SETS <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
          "stamps", "WDBC", "pima", "shuffle", "vowels")
BIG  <- c("PenDigits", "waveform", "thyroid", "pageblocks", "wilt")
METH <- c("SU-MCCD", "SUN-MCCD")
SMIN <- c(0, 0.03, 0.05, 0.0625)   # 0 is the gap; the rest anchor against
                                   # values WP2(c) already measured, so a
                                   # disagreement on those exposes a setup
                                   # difference rather than a real effect.

rows <- list()
for (ds in SETS) {
  dat <- load_real_dataset(ds)
  for (me in METH) {
    for (sm in SMIN) {
      t0 <- Sys.time()
      out <- tryCatch({
        res <- METHOD_REGISTRY[[me]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = sm)
        m <- evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[me]])
        list(m = m, res = res, status = "ok", note = NA_character_)
      }, error = function(e) list(m = setNames(rep(NA_real_, 4), c("TPR","TNR","BA","F2")),
                                  res = NULL, status = "error",
                                  note = conditionMessage(e)))
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      ncls <- if (is.null(out$res)) NA_integer_ else
        length(unique(out$res$cluster[!is.na(out$res$cluster)]))
      nfl <- if (is.null(out$res)) NA_integer_ else sum(out$res$score > 0.5)
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, method = me, s_min = sm, n = nrow(dat$X), d = dat$d,
        k_threshold = round(sm * nrow(dat$X)),
        TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
        BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
        n_flagged = nfl, n_clusters = ncls, elapsed_sec = el,
        status = out$status, note = out$note, stringsAsFactors = FALSE)
      cat(sprintf("CELL_DONE %s %s s_min=%.4f k=%d BA=%.3f F2=%.3f ncls=%s sec=%.1f\n",
                  ds, me, sm, round(sm * nrow(dat$X)),
                  ifelse(is.na(out$m[["BA"]]), -1, out$m[["BA"]]),
                  ifelse(is.na(out$m[["F2"]]), -1, out$m[["F2"]]),
                  ifelse(is.na(ncls), "NA", ncls), el))
      flush.console()
    }
  }
}
for (ds in BIG) for (me in METH) rows[[length(rows) + 1]] <- data.frame(
  dataset = ds, method = me, s_min = NA_real_, n = NA_integer_, d = NA_integer_,
  k_threshold = NA_integer_, TPR = NA_real_, TNR = NA_real_, BA = NA_real_,
  F2 = NA_real_, n_flagged = NA_integer_, n_clusters = NA_integer_,
  elapsed_sec = NA_real_, status = "not_run",
  note = "n > 1500; 350-1600 s per cell and the volume is busy with job 51",
  stringsAsFactors = FALSE)

df <- do.call(rbind, rows)
write.csv(df, OUT, row.names = FALSE)

ok <- df[df$status == "ok", ]
cat("\n=== mean over the 11 data sets, by method and s_min ===\n")
ag <- aggregate(cbind(BA, F2) ~ method + s_min, data = ok, FUN = mean)
ag$frac_ncls_1 <- aggregate(cbind(z = n_clusters == 1) ~ method + s_min,
                            data = ok, FUN = mean)$z
print(ag[order(ag$method, ag$s_min), ], row.names = FALSE, digits = 4)

cat("\n=== per data set: s_min = 0 against 0.0625 (the manuscript value) ===\n")
for (me in METH) {
  cat(sprintf("\n  %s\n", me))
  a <- ok[ok$method == me & ok$s_min == 0, ]
  b <- ok[ok$method == me & ok$s_min == 0.0625, ]
  m <- merge(a[c("dataset","BA","F2","n_clusters")],
             b[c("dataset","BA","F2","n_clusters")],
             by = "dataset", suffixes = c("_0", "_0625"))
  m$dBA <- m$BA_0 - m$BA_0625
  m$identical <- abs(m$dBA) < 1e-9 & abs(m$F2_0 - m$F2_0625) < 1e-9
  print(m[order(m$dBA), c("dataset","BA_0","BA_0625","dBA",
                          "n_clusters_0","n_clusters_0625","identical")],
        row.names = FALSE, digits = 4)
  cat(sprintf("    identical on %d/%d;  mean dBA = %+.4f;  clusters>1 at s_min=0: %d/%d\n",
              sum(m$identical), nrow(m), mean(m$dBA),
              sum(m$n_clusters_0 > 1), nrow(m)))
}
cat(sprintf("\nwrote %s\nALL_CELLS_COMPLETE\n", OUT))
