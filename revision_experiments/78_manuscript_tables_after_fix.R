#!/usr/bin/env Rscript
# 78_manuscript_tables_after_fix.R -- regenerate the two real-data tables the
# main text prints, with the alpha correction applied.
#
# Neither table can be patched cell by cell. tab:Real_Data_Result_Summary
# carries three DERIVED columns -- best overall method, best proposed method,
# and SUN-MCCD's rank -- so changing one F2 can flip a winner or move a rank on
# a row whose own numbers did not change. tab:Real_Data_Aggregate is ordered by
# mean F2, so a shifted value can reorder rows. Both are therefore rebuilt from
# the corrected per-cell table, and the script prints what it would change so
# the edit is made from arithmetic, not from reading down a column.
#
# Rounding follows the manuscript: values are printed to three decimals, and
# ranks and winners are computed on the ROUNDED values, because that is what a
# reader comparing printed numbers would do. Ties take the best shared rank,
# as the existing caption states.
#
# Read-only apart from its own CSV.

suppressMessages(library(here))
FC  <- here::here("revision_experiments/results/tr1/final_comparison.csv")
NEW <- here::here("revision_experiments/results/tr1/wp2c_rerun_alpha001.csv")
OUT <- here::here("revision_experiments/results/tr1/wp2c_manuscript_tables.csv")

fc <- read.csv(FC, stringsAsFactors = FALSE)
nw <- read.csv(NEW, stringsAsFactors = FALSE); nw <- nw[nw$status == "ok", ]

# ROUNDING. The manuscript rounds half UP; R's round() rounds half to even.
# The difference is not cosmetic here: SUN-MCCD's median F2 is exactly 0.3285,
# printed as 0.329 in tab:Real_Data_Aggregate and returned as 0.328 by round().
# Using R's default would have looked like a value the correction moved, when
# nothing about it changed.
r3 <- function(x) as.numeric(sprintf("%.3f", x))

# STALE CELL. final_comparison.csv still holds vertebral/SU-MCCD at the old
# S_min = 0.0625 operating point (BA 0.388, F2 0.140). The manuscript already
# carries the S_min = 0.05 values from 61_aggregate_impact_smin005.R. Left
# unpatched, the SU-MCCD aggregate row recomputes to numbers that disagree with
# the printed table for a reason unrelated to alpha.
k <- which(tolower(fc$dataset) == "vertebral" & fc$method == "SU-MCCD")
stopifnot(length(k) == 1)
fc$BA[k] <- 0.493; fc$F2[k] <- 0.129

PROPOSED <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
ORDER <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
           "stamps", "WDBC", "pima", "Shuttle", "vowels", "PenDigits",
           "waveform", "thyroid", "pageblocks", "wilt")

apply_fix <- function(df) {
  for (i in seq_len(nrow(nw))) {
    j <- which(tolower(df$dataset) == tolower(nw$dataset[i]) & df$method == nw$method[i])
    if (length(j) == 1) {
      df$F2[j] <- round(nw$F2[i], 3); df$BA[j] <- round(nw$BA[i], 3)
      df$TPR[j] <- round(nw$TPR[i], 3); df$TNR[j] <- round(nw$TNR[i], 3)
    }
  }
  df
}
fx <- apply_fix(fc)

# ---- tab:Real_Data_Result_Summary ------------------------------------------
summ <- function(df) {
  do.call(rbind, lapply(ORDER, function(ds) {
    g <- df[tolower(df$dataset) == tolower(ds), ]
    if (!nrow(g)) return(NULL)
    f <- round(g$F2, 3); names(f) <- g$method
    bo <- names(f)[f == max(f)]
    p  <- f[names(f) %in% PROPOSED]
    bp <- names(p)[p == max(p)]
    data.frame(dataset = ds,
               best_overall = paste(bo, collapse = "/"), f2_overall = max(f),
               best_proposed = paste(bp, collapse = "/"), f2_proposed = max(p),
               sun_rank = sum(f > f[["SUN-MCCD"]]) + 1,
               stringsAsFactors = FALSE)
  }))
}
s0 <- summ(fc); s1 <- summ(fx)

cat("=== tab:Real_Data_Result_Summary -- rows that change ===\n")
chg <- which(s0$best_overall != s1$best_overall | s0$f2_overall != s1$f2_overall |
             s0$best_proposed != s1$best_proposed | s0$f2_proposed != s1$f2_proposed |
             s0$sun_rank != s1$sun_rank)
if (!length(chg)) cat("  none\n") else for (i in chg) {
  cat(sprintf("  %-13s BEFORE %-18s %.3f | %-18s %.3f | rank %d\n", s0$dataset[i],
              s0$best_overall[i], s0$f2_overall[i], s0$best_proposed[i], s0$f2_proposed[i], s0$sun_rank[i]))
  cat(sprintf("  %-13s AFTER  %-18s %.3f | %-18s %.3f | rank %d\n", "",
              s1$best_overall[i], s1$f2_overall[i], s1$best_proposed[i], s1$f2_proposed[i], s1$sun_rank[i]))
}

cat("\n=== corrected LaTeX rows (paste into tab:Real_Data_Result_Summary) ===\n")
for (i in chg)
  cat(sprintf("  %-12s & %-16s & %.3f & %-16s & %.3f & %d \\\\ \\hline\n",
              s1$dataset[i], s1$best_overall[i], s1$f2_overall[i],
              s1$best_proposed[i], s1$f2_proposed[i], s1$sun_rank[i]))

# ---- tab:Real_Data_Aggregate -----------------------------------------------
agg <- function(df) {
  a <- do.call(rbind, lapply(split(df, df$method), function(g) data.frame(
    method = g$method[1],
    mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
    mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
    stringsAsFactors = FALSE)))
  a[order(-a$mean_F2), ]
}
a0 <- agg(fc); a1 <- agg(fx)

# ---- does the BEFORE state reproduce what the manuscript actually prints? ---
# If it does not, the AFTER numbers cannot be trusted to replace them. Values
# transcribed from tab:Real_Data_Aggregate, lines 1191-1199.
PRINTED <- data.frame(
  method    = c("SUN-MCCD","SU-MCCD","U-MCCD","iForest","LOF","MST","DBSCAN","ODIN","UN-MCCD"),
  mean_F2   = c(0.321, 0.302, 0.291, 0.260, 0.254, 0.237, 0.234, 0.229, 0.229),
  median_F2 = c(0.329, 0.301, 0.280, 0.115, 0.195, 0.237, 0.133, 0.224, 0.209),
  mean_BA   = c(0.722, 0.696, 0.709, 0.632, 0.654, 0.606, 0.599, 0.619, 0.658),
  median_BA = c(0.726, 0.692, 0.708, 0.526, 0.577, 0.599, 0.546, 0.585, 0.634),
  stringsAsFactors = FALSE)
cat("\n=== does BEFORE reproduce the printed table? ===\n")
m <- merge(PRINTED, a0, by = "method", suffixes = c("_printed", "_computed"))
bad <- 0L
for (col in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  d <- abs(m[[paste0(col, "_printed")]] - m[[paste0(col, "_computed")]])
  if (any(d > 0.0005)) {
    bad <- bad + sum(d > 0.0005)
    for (i in which(d > 0.0005))
      cat(sprintf("  MISMATCH %-9s %-10s printed %.3f, computed %.3f\n",
                  m$method[i], col, m[[paste0(col, "_printed")]][i], m[[paste0(col, "_computed")]][i]))
  }
}
cat(if (bad == 0) "  all 36 printed values reproduced exactly\n"
    else sprintf("  *** %d values do not reproduce -- do not edit the manuscript from this ***\n", bad))
cat("\n=== tab:Real_Data_Aggregate BEFORE ===\n"); print(a0, row.names = FALSE)
cat("\n=== tab:Real_Data_Aggregate AFTER ===\n");  print(a1, row.names = FALSE)
cat("\n  row order changed: ", !identical(a0$method, a1$method), "\n", sep = "")

cat("\n=== corrected LaTeX rows (paste into tab:Real_Data_Aggregate) ===\n")
for (i in seq_len(nrow(a1))) {
  m <- a1$method[i]
  bold <- m == "SUN-MCCD"
  fmt <- function(x) if (bold) sprintf("\\textbf{%.3f}", x) else sprintf("%.3f", x)
  cat(sprintf("  %s & %s & %s & %s & %s \\\\ \\hline\n",
              if (bold) sprintf("\\textbf{%s}", m) else sprintf("%-8s", m),
              fmt(a1$mean_F2[i]), fmt(a1$median_F2[i]), fmt(a1$mean_BA[i]), fmt(a1$median_BA[i])))
}

# ---- the counting claims in the surrounding prose --------------------------
cat("\n=== claims in the prose ===\n")
cnt <- function(s) c(prop_wins = sum(sapply(strsplit(s$best_overall, "/"), function(x) any(x %in% PROPOSED))),
                     sun_best_proposed = sum(sapply(strsplit(s$best_proposed, "/"), function(x) "SUN-MCCD" %in% x)))
c0 <- cnt(s0); c1 <- cnt(s1)
cat(sprintf("  a proposed method attains the highest F2 on : %d of 16  ->  %d of 16\n", c0[1], c1[1]))
cat(sprintf("  SUN-MCCD strongest of the four proposed on  : %d of 16  ->  %d of 16\n", c0[2], c1[2]))
for (m in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  v <- a1[[m]]; names(v) <- a1$method
  cat(sprintf("  SUN-MCCD first on %-10s: %s (%.3f vs next %.3f)\n", m,
              which.max(v) == which(a1$method == "SUN-MCCD"),
              v[["SUN-MCCD"]], max(v[names(v) != "SUN-MCCD"])))
}

write.csv(merge(s1, a1, all = TRUE), OUT, row.names = FALSE)
cat(sprintf("\nwrote %s\n", OUT))
