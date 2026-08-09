#!/usr/bin/env Rscript
# 35_print_tables.R -- emit the regenerated real-data numbers in the exact
# row/column order the manuscript tables need, so the LaTeX can be written
# straight from this output without hand-transcription.
#
# Read-only.

suppressMessages(library(here))
cmp <- read.csv(here::here("revision_experiments/results/final_comparison.csv"),
                stringsAsFactors = FALSE)

# Manuscript order: original eight first (as currently printed), then the new
# eight in ascending n.
DS <- c("hepatitis", "glass", "vertebral", "ecoli", "stamps", "vowels", "waveform", "wilt",
        "lymphography", "WBC", "WDBC", "pima", "shuffle", "PenDigits", "thyroid", "pageblocks")
M  <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD", "LOF", "DBSCAN", "MST", "ODIN", "iForest")

g <- function(ds, m, k) {
  r <- cmp[cmp$dataset == ds & cmp$method == m, ]
  if (!nrow(r)) return(NA_real_)
  as.numeric(r[[k]][1])
}
f3 <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)

for (pair in list(c("TPR", "TNR"), c("BA", "F2"))) {
  cat(sprintf("\n########## %s / %s -- supplement table rows ##########\n", pair[1], pair[2]))
  for (ds in DS) {
    cat(sprintf("  %-13s", ds))
    for (m in M) cat(sprintf(" & %s & %s", f3(g(ds, m, pair[1])), f3(g(ds, m, pair[2]))))
    cat(" \\\\ \\hline\n")
  }
}

# ---- Table 6: best overall, best proposed, rank of best proposed ----------
cat("\n########## Table 6 summary rows ##########\n")
PROP <- M[1:4]
for (ds in DS) {
  v  <- setNames(sapply(M,    function(m) g(ds, m, "F2")), M)
  vp <- setNames(sapply(PROP, function(m) g(ds, m, "F2")), PROP)
  bo <- names(which.max(v)); bp <- names(which.max(vp))
  # Column 5 of Table 6 is SUN-MCCD's rank among ALL NINE methods (ties take
  # the best shared rank), not the best-proposed method's rank.
  rk <- rank(-v, ties.method = "min")[["SUN-MCCD"]]
  # name every method tied for best overall, so "UN-MCCD/SUN-MCCD" style
  # entries survive as in the published table
  bo_all <- paste(names(v)[!is.na(v) & v >= max(v, na.rm = TRUE) - 1e-9], collapse = "/")
  bp_all <- paste(names(vp)[!is.na(vp) & vp >= max(vp, na.rm = TRUE) - 1e-9], collapse = "/")
  cat(sprintf("  %-13s & %s & %s & %s & %s & %d \\\\ \\hline\n",
              ds, bo_all, f3(max(v, na.rm = TRUE)), bp_all, f3(max(vp, na.rm = TRUE)), rk))
}

# ---- means, for the strengthened abstract claim --------------------------
cat("\n########## mean F2 / BA over all 16, ranked ##########\n")
mf <- sapply(M, function(m) mean(sapply(DS, function(d) g(d, m, "F2")), na.rm = TRUE))
mb <- sapply(M, function(m) mean(sapply(DS, function(d) g(d, m, "BA")), na.rm = TRUE))
o  <- order(-mf)
for (i in o) cat(sprintf("  %-9s  mean F2 %.3f   mean BA %.3f\n", M[i], mf[i], mb[i]))
cat(sprintf("\n  SUN-MCCD F2 rank: %d of %d;  BA rank: %d of %d\n",
            rank(-mf)[["SUN-MCCD"]], length(M), rank(-mb)[["SUN-MCCD"]], length(M)))

# median too -- a mean can be driven by one data set, and we should know
cat("\n########## median F2 / BA over all 16, ranked ##########\n")
df_ <- sapply(M, function(m) median(sapply(DS, function(d) g(d, m, "F2")), na.rm = TRUE))
db_ <- sapply(M, function(m) median(sapply(DS, function(d) g(d, m, "BA")), na.rm = TRUE))
for (i in order(-df_)) cat(sprintf("  %-9s  median F2 %.3f   median BA %.3f\n", M[i], df_[i], db_[i]))
cat(sprintf("\n  SUN-MCCD median F2 rank: %d of %d;  median BA rank: %d of %d\n",
            rank(-df_)[["SUN-MCCD"]], length(M), rank(-db_)[["SUN-MCCD"]], length(M)))
