#!/usr/bin/env Rscript
# 76_aggregate_after_alpha_fix.R -- what the corrected alpha does to the
# aggregate table, and to the headline claim resting on it.
#
# Six real-data cells were recomputed against genuine 0.1% NND quantile tables
# (WP2c): SUN-MCCD on vowels, PenDigits, lymphography, hepatitis and waveform,
# and UN-MCCD on waveform. Every other cell is untouched -- UN-MCCD resolves to
# the "99" token below d=20, and those tables were never in question.
#
# The claim under test is the manuscript's headline, stated in the abstract and
# the conclusion: across the sixteen data sets SUN-MCCD ranks FIRST on mean BA,
# median BA, mean F2 and median F2, ahead of all five baselines. Five of the
# sixteen SUN-MCCD rows just moved. Whether that claim survives is not
# something to assume in either direction, so it is recomputed here from the
# same table the manuscript prints.
#
# Read-only apart from its own CSV.

suppressMessages(library(here))
FC   <- here::here("revision_experiments/results/tr1/final_comparison.csv")
NEW  <- here::here("revision_experiments/results/tr1/wp2c_rerun_alpha001.csv")
OUT  <- here::here("revision_experiments/results/tr1/wp2c_aggregate_after_fix.csv")

fc  <- read.csv(FC,  stringsAsFactors = FALSE)
nw  <- read.csv(NEW, stringsAsFactors = FALSE)
nw  <- nw[nw$status == "ok", ]

cat("=== cell-level change ===\n")
chg <- list()
for (i in seq_len(nrow(nw))) {
  ds <- nw$dataset[i]; me <- nw$method[i]
  j <- which(tolower(fc$dataset) == tolower(ds) & fc$method == me)
  if (length(j) != 1) { cat(sprintf("  %-13s %-9s NOT FOUND in final_comparison\n", ds, me)); next }
  chg[[length(chg) + 1]] <- data.frame(
    dataset = ds, method = me,
    old_BA = fc$BA[j], new_BA = round(nw$BA[i], 3), d_BA = round(nw$BA[i] - fc$BA[j], 3),
    old_F2 = fc$F2[j], new_F2 = round(nw$F2[i], 3), d_F2 = round(nw$F2[i] - fc$F2[j], 3),
    stringsAsFactors = FALSE)
}
cd <- do.call(rbind, chg)
print(cd, row.names = FALSE)
cat(sprintf("\n  materially changed (|delta| >= 0.005 on either metric): %d of %d\n",
            sum(abs(cd$d_BA) >= 0.005 | abs(cd$d_F2) >= 0.005), nrow(cd)))

# ---- apply the corrections -------------------------------------------------
fx <- fc
for (i in seq_len(nrow(nw))) {
  j <- which(tolower(fx$dataset) == tolower(nw$dataset[i]) & fx$method == nw$method[i])
  if (length(j) == 1) {
    fx$BA[j]  <- round(nw$BA[i],  3); fx$F2[j]  <- round(nw$F2[i],  3)
    fx$TPR[j] <- round(nw$TPR[i], 3); fx$TNR[j] <- round(nw$TNR[i], 3)
  }
}

agg <- function(df) {
  s <- do.call(rbind, lapply(split(df, df$method), function(g) data.frame(
    method = g$method[1], n = nrow(g),
    mean_BA = round(mean(g$BA), 3), median_BA = round(median(g$BA), 3),
    mean_F2 = round(mean(g$F2), 3), median_F2 = round(median(g$F2), 3),
    stringsAsFactors = FALSE)))
  s[order(-s$mean_BA), ]
}
before <- agg(fc); after <- agg(fx)

cat("\n=== aggregate table BEFORE (as printed in the manuscript) ===\n")
print(before, row.names = FALSE)
cat("\n=== aggregate table AFTER the alpha correction ===\n")
print(after, row.names = FALSE)

cat("\n=== the headline claim: does SUN-MCCD still rank first on all four? ===\n")
verdict <- TRUE
for (m in c("mean_BA", "median_BA", "mean_F2", "median_F2")) {
  v <- after[[m]]; names(v) <- after$method
  best <- names(v)[which.max(v)]
  sun  <- v[["SUN-MCCD"]]
  rank <- sum(v > sun) + 1
  tied <- sum(v == sun) - 1
  ok <- rank == 1
  verdict <- verdict && ok
  cat(sprintf("  %-10s SUN-MCCD = %.3f, rank %d of %d%s, best = %s (%.3f)  -> %s\n",
              m, sun, rank, length(v),
              if (tied > 0) sprintf(" (tied with %d)", tied) else "",
              best, max(v), if (ok) "HOLDS" else "*** FAILS ***"))
}
cat(sprintf("\nVERDICT: the headline claim %s\n",
            if (verdict) "SURVIVES the correction" else "*** DOES NOT SURVIVE -- the manuscript must change ***"))

# margin to the nearest competitor, since "first" can be first by nothing
cat("\n=== margin over the best non-proposed baseline ===\n")
BASE <- c("LOF", "DBSCAN", "MST", "ODIN", "iForest")
for (m in c("mean_BA", "median_BA", "mean_F2", "median_F2")) {
  v <- after[[m]]; names(v) <- after$method
  cat(sprintf("  %-10s SUN-MCCD %.3f vs best baseline %.3f (%s)  margin %+.3f\n",
              m, v[["SUN-MCCD"]], max(v[BASE]), names(which.max(v[BASE])),
              v[["SUN-MCCD"]] - max(v[BASE])))
}

write.csv(cbind(after, stage = "after_alpha_fix"), OUT, row.names = FALSE)
cat(sprintf("\nwrote %s\n", OUT))
