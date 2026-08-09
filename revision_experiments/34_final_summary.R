#!/usr/bin/env Rscript
# 34_final_summary.R -- the complete regenerated Section 6 comparison,
# 16 data sets x 9 methods, at the configuration frozen in
# BENCHMARK_EXPANSION_RULE.md.
#
# Sources:
#   proposed   regen_final_*.csv                     (S_min = 0.0625)
#   DBSCAN     original 8: regen_baselines.csv variant "oracle-X"
#              new 8:      final_baselines.csv
#              (= published protocol, df->X label-leak defect fixed)
#   MST        original 8: regen_baselines.csv variant "uniform-1.2"
#              new 8:      final_baselines.csv       (cont=0.02, thresh=1.2)
#   ODIN       "default" / final_baselines.csv
#   iForest    "seed1"   / final_baselines.csv
#   LOF        original 8: PUBLISHED (not re-run, per user decision)
#              new 8:      final_baselines.csv
#
# Writes results/final_comparison.csv and prints the tables.

suppressMessages(library(here))
RES <- here::here("revision_experiments/results")

ORIG <- c("hepatitis", "glass", "vertebral", "ecoli", "stamps", "vowels", "waveform", "wilt")
NEW  <- c("lymphography", "WBC", "WDBC", "pima", "shuffle", "PenDigits", "thyroid", "pageblocks")
ALL_DS   <- c(ORIG, NEW)
PROPOSED <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
ALL_M    <- c(PROPOSED, "LOF", "DBSCAN", "MST", "ODIN", "iForest")

rd <- function(f) {
  p <- file.path(RES, f)
  if (!file.exists(p)) { message("  (missing: ", f, ")"); return(NULL) }
  read.csv(p, stringsAsFactors = FALSE)
}

prop <- do.call(rbind, Filter(Negate(is.null), lapply(
  c("regen_final_small.csv", "regen_final_vowels.csv", "regen_final_waveform.csv",
    "regen_final_wilt.csv", "regen_final_new_small.csv", "regen_final_pendigits.csv",
    "regen_final_thyroid.csv", "regen_final_pageblocks.csv"),
  function(f) { d <- rd(f); if (is.null(d)) NULL else d[, c("dataset","method","TPR","TNR","BA","F2")] })))

bo <- rd("regen_baselines.csv")            # original 8, variant-keyed
bn <- rd("final_baselines.csv")            # new 8, frozen config

TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
pub <- function(ds, m, mt) {
  v <- TRUTH$value[TRUTH$dataset == ds & TRUTH$method == m & TRUTH$metric == mt]
  if (length(v) == 0) NA_real_ else v[1]
}

VARIANT <- c(DBSCAN = "oracle-X", MST = "uniform-1.2", ODIN = "default", iForest = "seed1")

val <- function(ds, m, mt) {
  if (m %in% PROPOSED) {
    r <- prop[prop$dataset == ds & prop$method == m, ]
    return(if (nrow(r)) round(as.numeric(r[[mt]][1]), 3) else NA_real_)
  }
  if (ds %in% ORIG) {
    if (m == "LOF") return(pub(ds, m, mt))                    # not re-run
    if (is.null(bo)) return(NA_real_)
    r <- bo[bo$dataset == ds & bo$method == m & bo$variant == VARIANT[[m]], ]
  } else {
    if (is.null(bn)) return(NA_real_)
    r <- bn[bn$dataset == ds & bn$method == m, ]
  }
  if (nrow(r)) round(as.numeric(r[[mt]][1]), 3) else NA_real_
}

# ---- full tables --------------------------------------------------------
out <- list()
for (mt in c("F2", "BA")) {
  cat(sprintf("\n===== %s, 16 data sets x 9 methods =====\n", mt))
  cat(sprintf("%-13s", "dataset")); for (m in ALL_M) cat(sprintf("%9s", m)); cat("\n")
  cat(strrep("-", 13 + 9 * length(ALL_M)), "\n")
  for (ds in ALL_DS) {
    if (ds == NEW[1]) cat(sprintf("%-13s%s\n", "-- new --", strrep(" ", 9 * length(ALL_M))))
    cat(sprintf("%-13s", ds))
    for (m in ALL_M) {
      v <- val(ds, m, mt)
      cat(if (is.na(v)) sprintf("%9s", "-") else sprintf("%9.3f", v))
      if (mt == "F2") out[[length(out) + 1]] <-
        data.frame(dataset = ds, method = m, F2 = v,
                   BA = val(ds, m, "BA"), TPR = val(ds, m, "TPR"),
                   TNR = val(ds, m, "TNR"),
                   group = if (ds %in% ORIG) "original" else "new",
                   stringsAsFactors = FALSE)
    }
    cat("\n")
  }
}

cmp <- do.call(rbind, out)
write.csv(cmp, file.path(RES, "final_comparison.csv"), row.names = FALSE)

# ---- who wins -----------------------------------------------------------
cat("\n===== best F2 per data set =====\n")
cat(sprintf("%-13s %-24s %-24s %s\n", "dataset", "best overall", "best proposed", "proposed?"))
wins <- 0; sun <- 0; nprop <- 0
for (ds in ALL_DS) {
  v  <- sapply(ALL_M,    function(m) val(ds, m, "F2"))
  vp <- sapply(PROPOSED, function(m) val(ds, m, "F2"))
  if (all(is.na(vp))) { cat(sprintf("%-13s (proposed not yet run)\n", ds)); next }
  nprop <- nprop + 1
  bo_ <- names(which.max(v)); bp_ <- names(which.max(vp))
  won <- bo_ %in% PROPOSED; wins <- wins + won
  if (!is.na(vp[["SUN-MCCD"]]) && vp[["SUN-MCCD"]] >= max(vp, na.rm = TRUE) - 1e-9) sun <- sun + 1
  cat(sprintf("%-13s %-24s %-24s %s\n", ds,
              sprintf("%s (%.3f)", bo_, max(v, na.rm = TRUE)),
              sprintf("%s (%.3f)", bp_, max(vp, na.rm = TRUE)),
              if (won) "YES" else ""))
}
cat(sprintf("\n  proposed methods win %d of %d data sets\n", wins, nprop))
cat(sprintf("  SUN-MCCD best-or-tied among the four on %d of %d\n", sun, nprop))

# ---- mean performance ---------------------------------------------------
cat("\n===== mean F2 / BA across data sets =====\n")
cat(sprintf("%-9s %18s %18s %18s\n", "method", "all 16", "original 8", "new 8"))
for (m in ALL_M) {
  f_all <- mean(sapply(ALL_DS, function(d) val(d, m, "F2")), na.rm = TRUE)
  b_all <- mean(sapply(ALL_DS, function(d) val(d, m, "BA")), na.rm = TRUE)
  f_o   <- mean(sapply(ORIG,   function(d) val(d, m, "F2")), na.rm = TRUE)
  b_o   <- mean(sapply(ORIG,   function(d) val(d, m, "BA")), na.rm = TRUE)
  f_n   <- mean(sapply(NEW,    function(d) val(d, m, "F2")), na.rm = TRUE)
  b_n   <- mean(sapply(NEW,    function(d) val(d, m, "BA")), na.rm = TRUE)
  cat(sprintf("%-9s   F2 %.3f BA %.3f   F2 %.3f BA %.3f   F2 %.3f BA %.3f\n",
              m, f_all, b_all, f_o, b_o, f_n, b_n))
}

# ---- dimension trend for the proposed four ------------------------------
cat("\n===== proposed-method F2 by dimension (AE.4) =====\n")
dims <- c(hepatitis=19, glass=9, vertebral=6, ecoli=7, stamps=9, vowels=12,
          waveform=21, wilt=5, lymphography=18, WBC=9, WDBC=30, pima=8,
          shuffle=9, PenDigits=16, thyroid=6, pageblocks=10)
ord <- names(sort(dims))
cat(sprintf("%-13s %4s", "dataset", "d")); for (m in PROPOSED) cat(sprintf("%10s", m)); cat("      best baseline\n")
for (ds in ord) {
  vp <- sapply(PROPOSED, function(m) val(ds, m, "F2"))
  if (all(is.na(vp))) next
  vb <- sapply(setdiff(ALL_M, PROPOSED), function(m) val(ds, m, "F2"))
  cat(sprintf("%-13s %4d", ds, dims[[ds]]))
  for (m in PROPOSED) { v <- vp[[m]]; cat(if (is.na(v)) sprintf("%10s","-") else sprintf("%10.3f", v)) }
  cat(sprintf("      %.3f\n", max(vb, na.rm = TRUE)))
}

cat("\nWritten: results/final_comparison.csv\n")
