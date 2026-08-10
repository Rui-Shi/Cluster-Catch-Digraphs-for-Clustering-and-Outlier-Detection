#!/usr/bin/env Rscript
# 61_aggregate_impact_smin005.R -- exactly which printed table values change
# when the real-data operating point moves from S_min = 0.0625 to 0.05.
#
# Across all 16 data sets and both shape-adaptive methods, 31 of 32 cells are
# byte-identical at the two settings. The single exception is
# vertebral / SU-MCCD. This computes the knock-on effect on the aggregate table
# to three decimals, the precision the manuscript prints, so the edit is made
# from arithmetic rather than from a guess about rounding.
#
# Note the change is NOT a pure improvement: BA rises 0.388 -> 0.493 but F2
# falls 0.140 -> 0.129. Both directions must be reported.
#
# Read-only.

suppressMessages(library(here))
fc <- read.csv(here::here("revision_experiments/results/tr1/final_comparison.csv"),
               stringsAsFactors = FALSE)

NEW_BA <- 0.492857142857143
NEW_F2 <- 0.129032258064516

for (meth in c("SU-MCCD", "SUN-MCCD")) {
  s <- fc[fc$method == meth, ]
  i <- which(s$dataset == "vertebral")
  oldBA <- s$BA; oldF2 <- s$F2
  newBA <- oldBA; newF2 <- oldF2
  changed <- (meth == "SU-MCCD")
  if (changed) { newBA[i] <- NEW_BA; newF2[i] <- NEW_F2 }
  cat(sprintf("\n=== %s (n = %d data sets)%s ===\n", meth, nrow(s),
              if (changed) "  [vertebral swapped]" else "  [unaffected]"))
  if (changed)
    cat(sprintf("  vertebral: BA %.4f -> %.4f ; F2 %.4f -> %.4f\n",
                oldBA[i], NEW_BA, oldF2[i], NEW_F2))
  f <- function(lab, o, n) cat(sprintf(
    "  %-9s %.5f -> %.5f    printed %.3f -> %.3f    %s\n", lab, o, n,
    round(o, 3), round(n, 3),
    if (round(o, 3) == round(n, 3)) "unchanged at 3dp" else "*** TABLE EDIT ***"))
  f("mean BA",   mean(oldBA),   mean(newBA))
  f("median BA", median(oldBA), median(newBA))
  f("mean F2",   mean(oldF2),   mean(newF2))
  f("median F2", median(oldF2), median(newF2))
}

# Does the aggregate table's ordering (by mean F2) change?
cat("\n=== ordering by mean F2, before and after ===\n")
agg <- aggregate(cbind(BA, F2) ~ method, data = fc, FUN = mean)
agg2 <- agg
su <- fc[fc$method == "SU-MCCD", ]
i <- which(su$dataset == "vertebral")
f2new <- su$F2; f2new[i] <- NEW_F2
banew <- su$BA; banew[i] <- NEW_BA
agg2$F2[agg2$method == "SU-MCCD"] <- mean(f2new)
agg2$BA[agg2$method == "SU-MCCD"] <- mean(banew)
o1 <- agg$method[order(-agg$F2)]; o2 <- agg2$method[order(-agg2$F2)]
cat(sprintf("  before: %s\n", paste(o1, collapse = " > ")))
cat(sprintf("  after : %s\n", paste(o2, collapse = " > ")))
cat(sprintf("  -> ordering %s\n", if (identical(o1, o2)) "UNCHANGED" else "CHANGES"))
cat(sprintf("\n  SUN-MCCD still first on mean F2: %s ; on mean BA: %s\n",
            o2[1] == "SUN-MCCD",
            agg2$method[order(-agg2$BA)][1] == "SUN-MCCD"))
cat("\ndone\n")
