#!/usr/bin/env Rscript
# 31_headline_impact.R -- does regeneration help or hurt the four CCD methods?
#
# Assembles the regenerated Section 6 comparison from:
#   proposed   regen_final_{small,vowels,waveform,wilt}.csv  (S_min = 0.0625)
#   DBSCAN     regen_baselines.csv, variant "knee-4"       (label-free knee)
#   MST        regen_baselines.csv, variant "uniform-1.2"  (one setting, all sets)
#   ODIN       regen_baselines.csv, variant "default"
#   iForest    regen_baselines.csv, variant "seed1"
#   LOF        published (not re-run -- already label-free and correctly configured)
#
# and answers three questions:
#   1. per method, does F2/BA go up or down on average?
#   2. how much of that is waveform, the one dataset we cannot explain?
#   3. does the paper's competitive standing change -- who wins each dataset,
#      and on how many is SUN-MCCD the best of the four?
#
# Read-only. Writes nothing.

suppressMessages(library(here))
RES <- here::here("revision_experiments/results/tr1")

DATASETS <- c("hepatitis", "glass", "vertebral", "ecoli",
              "stamps", "vowels", "waveform", "wilt")
PROPOSED <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")

prop <- do.call(rbind, lapply(
  file.path(RES, c("regen_final_small.csv", "regen_final_vowels.csv",
                   "regen_final_waveform.csv", "regen_final_wilt.csv")),
  function(f) read.csv(f, stringsAsFactors = FALSE)))

base <- read.csv(file.path(RES, "regen_baselines.csv"), stringsAsFactors = FALSE)
pick <- rbind(
  base[base$method == "DBSCAN"  & base$variant == "knee-4", ],
  base[base$method == "MST"     & base$variant == "uniform-1.2", ],
  base[base$method == "ODIN"    & base$variant == "default", ],
  base[base$method == "iForest" & base$variant == "seed1", ])

TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
pub <- function(ds, meth, metric) {
  v <- TRUTH$value[TRUTH$dataset == ds & TRUTH$method == meth & TRUTH$metric == metric]
  if (length(v) == 0) NA_real_ else v[1]
}

get_new <- function(ds, meth, metric) {
  if (meth == "LOF") return(pub(ds, meth, metric))          # not re-run
  src <- if (meth %in% PROPOSED) prop else pick
  r <- src[src$dataset == ds & src$method == meth, ]
  if (nrow(r) == 0) return(NA_real_)
  round(as.numeric(r[[metric]][1]), 3)
}

# ---- 1 & 2: per-method change ------------------------------------------
for (metric in c("F2", "BA")) {
  cat(sprintf("\n===== %s: regenerated minus published =====\n", metric))
  cat(sprintf("%-10s %8s %8s %8s %8s\n", "dataset", PROPOSED[1], PROPOSED[2], PROPOSED[3], PROPOSED[4]))
  D <- matrix(NA_real_, length(DATASETS), 4, dimnames = list(DATASETS, PROPOSED))
  for (ds in DATASETS) {
    for (m in PROPOSED) D[ds, m] <- get_new(ds, m, metric) - pub(ds, m, metric)
    cat(sprintf("%-10s %+8.3f %+8.3f %+8.3f %+8.3f\n", ds, D[ds, 1], D[ds, 2], D[ds, 3], D[ds, 4]))
  }
  cat(strrep("-", 47), "\n")
  cat(sprintf("%-10s %+8.3f %+8.3f %+8.3f %+8.3f\n", "mean all", mean(D[, 1]), mean(D[, 2]), mean(D[, 3]), mean(D[, 4])))
  k <- rownames(D) != "waveform"
  cat(sprintf("%-10s %+8.3f %+8.3f %+8.3f %+8.3f\n", "ex-wavefm",
              mean(D[k, 1]), mean(D[k, 2]), mean(D[k, 3]), mean(D[k, 4])))
  k2 <- !(rownames(D) %in% c("waveform", "ecoli"))
  cat(sprintf("%-10s %+8.3f %+8.3f %+8.3f %+8.3f\n", "ex-wf+eco",
              mean(D[k2, 1]), mean(D[k2, 2]), mean(D[k2, 3]), mean(D[k2, 4])))
}

# ---- 3: competitive standing -------------------------------------------
ALL <- c(PROPOSED, "LOF", "DBSCAN", "MST", "ODIN", "iForest")

standing <- function(getter, label) {
  cat(sprintf("\n===== %s: best F2 per data set =====\n", label))
  cat(sprintf("%-10s %-22s %-22s %s\n", "dataset", "best overall", "best proposed", "proposed wins?"))
  wins <- 0; sun_best <- 0
  for (ds in DATASETS) {
    v  <- sapply(ALL,      function(m) getter(ds, m, "F2"))
    vp <- sapply(PROPOSED, function(m) getter(ds, m, "F2"))
    bo <- names(which.max(v)); bp <- names(which.max(vp))
    won <- bo %in% PROPOSED
    wins <- wins + won
    # SUN-MCCD best-or-tied among the four?
    if (!is.na(vp[["SUN-MCCD"]]) && vp[["SUN-MCCD"]] >= max(vp, na.rm = TRUE) - 1e-9) sun_best <- sun_best + 1
    cat(sprintf("%-10s %-22s %-22s %s\n", ds,
                sprintf("%s (%.3f)", bo, max(v, na.rm = TRUE)),
                sprintf("%s (%.3f)", bp, max(vp, na.rm = TRUE)),
                if (won) "yes" else "no"))
  }
  cat(sprintf("\n  proposed methods win %d of 8\n", wins))
  cat(sprintf("  SUN-MCCD is best-or-tied among the four on %d of 8\n", sun_best))
}

standing(function(ds, m, mt) pub(ds, m, mt),   "PUBLISHED")
standing(function(ds, m, mt) get_new(ds, m, mt), "REGENERATED")

cat("\ndone\n")
