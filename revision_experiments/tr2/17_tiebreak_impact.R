#!/usr/bin/env Rscript
# 17_tiebreak_impact.R -- how much does the tie-breaking rule move the
# published real-data OS numbers?
#
# WHY (2026-08-22): the shipped tie-breaking loop computes
#     IOS_(k-1) + (IOS_(k+m) - IOS_(k-1)) * rho_i / sum(rho),
# whereas manuscript Eq. (10) is
#     IOS_(k+m) - (IOS_(k+m) - IOS_(k-1)) * rho_i / sum(rho),
# the mirror image about the midpoint of the bracket. The within-block order is
# therefore reversed relative to what the manuscript says it is ("points with
# smaller vicinity density receive larger IOS values within the tie block").
# break_ties() in Outlyingness_Score.R now offers both, plus "none".
#
# This script bounds the cost of the fix before anyone decides whether the
# reported tables must be regenerated. The CCD construction is paid ONCE per
# (dataset, method) -- the tie-break is a pure post-processing step on the
# score vector -- so all modes are compared on identical digraphs and identical
# raw scores. Any difference reported here is attributable to the tie rule and
# nothing else.
#
# Modes compared:
#   legacy       as published (bit-exact reproduction of the old loop)
#   eq10         manuscript Eq. (10), same scope as legacy (whole dataset)
#   none         ties retained, i.e. Ceyhan's suggested wording
#   eq10_cluster Eq. (10) applied WITHIN each cluster, which is the scope the
#                manuscript actually states. IOS only; the manuscript defines
#                no tie rule for OOS, so OOS rows record NA here.
#
# Reported per (dataset, method, mode): TPR/TNR/BA/F2 at the published cutoff
# of 2, plus a tie-aware rank AUC (tie-breaking changes the ranking, which the
# fixed-cutoff metrics can hide), plus tie counts.
#
# CLI (PowerShell; Rscript segfaults under Bash):
#   Rscript "revision_experiments/tr2/17_tiebreak_impact.R" [dataset ...]
# With no arguments it runs the seven fast datasets; thyroid and wilt are slow
# (n = 3772 / 4819 under RK-CCD) and are opt-in by name.

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments", "shared", "harness.R")))

OUT_CSV <- here::here("revision_experiments/results/tr2/tiebreak_impact.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
MODES   <- c("legacy", "eq10", "none", "eq10_cluster")
CUTOFF  <- 2   # REAL_DATA_THRESHOLDS for all four OS methods

# Manuscript Tables 12/13, transcribed 2026-08-16 (same list as
# 15_os_repro_audit.R). dataset -> method -> c(TPR, TNR, BA, F2).
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

# R object name in RealData_Collection.R (the manuscript prints lymphography
# as "lymph").
ROBJ <- c(hepatitis = "hepatitis", lymphography = "lymphography",
          glass = "glass", WBC = "WBC", vertebral = "vertebral",
          stamps = "stamps", WDBC = "WDBC", vowels = "vowels",
          thyroid = "thyroid", wilt = "wilt")

FAST_SETS <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral",
               "stamps", "WDBC", "vowels")

args <- commandArgs(trailingOnly = TRUE)
SETS <- if (length(args) > 0) args else FAST_SETS
stopifnot(all(SETS %in% names(PUBLISHED)))

# --- helpers ----------------------------------------------------------------

#' Tie-aware rank AUC (Mann-Whitney with midranks), the "standard treatment of
#' tied scores" the manuscript's ties paragraph refers to. Y: 1 = regular,
#' 0 = outlier; higher score = more outlying.
rank_auc <- function(Y, score) {
  if (any(is.na(score))) return(NA_real_)
  pos <- score[Y == 0]; neg <- score[Y == 1]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

#' Number of observations that share their value with at least one other
#' observation (i.e. sit inside a tie block).
n_tied <- function(x) {
  tb <- table(x)
  sum(tb[tb > 1])
}

#' Eq. (10) applied within each cluster, the scope the manuscript states.
break_ties_within_clusters <- function(scores, vd, label) {
  out <- scores
  for (m in unique(label)) {
    idx <- which(label == m)
    if (length(idx) < 2) next
    out[idx] <- break_ties(scores[idx], vd[idx], mode = "eq10")
  }
  out
}

#' One (dataset, method) cell: pay the construction once, return the parts.
score_parts <- function(meth, X, d) {
  tab <- if (startsWith(meth, "RKCCD")) get_simul("RK", d) else get_simul("NN", d)
  if (meth == "RKCCD-OOS") {
    RKCCD_OOS(datax = X, simul = tab$simul, d = d, quant = tab$quant, parts = TRUE)
  } else if (meth == "RKCCD-IOS") {
    RKCCD_IOS(datax = X, simul = tab$simul, d = d, quant = tab$quant,
              min.cls = 0, parts = TRUE)
  } else if (meth == "UNCCD-OOS") {
    NNCCD_OOS(datax = X, simul = tab$simul, method = "ascend", d = d, parts = TRUE)
  } else {
    NNCCD_IOS(datax = X, simul = tab$simul, method = "ascend", d = d,
              min.cls = 0, parts = TRUE)
  }
}

# --- main -------------------------------------------------------------------

rows <- list()
for (ds in SETS) {
  dat <- load_real_dataset(ROBJ[[ds]])
  cat(sprintf("\n== %s: n=%d, d=%d, outliers=%d\n", ds, dat$n, dat$d, sum(dat$Y == 0)))

  for (meth in METHODS) {
    t0 <- Sys.time()
    got <- tryCatch(score_parts(meth, dat$X, dat$d),
                    error = function(e) { cat(sprintf("  %-10s ERROR: %s\n", meth,
                                                      conditionMessage(e))); NULL })
    if (is.null(got)) next
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    raw <- got$raw
    variants <- list(
      legacy       = break_ties(raw, got$vd, mode = "legacy"),
      eq10         = break_ties(raw, got$vd, mode = "eq10"),
      none         = raw,
      eq10_cluster = if (is.null(got$label)) NULL
                     else break_ties_within_clusters(raw, got$vd, got$label)
    )

    pub <- PUBLISHED[[ds]][[meth]]
    base <- NULL
    cat(sprintf("  %-10s n_tied(raw)=%d of %d  (%.1fs)\n", meth, n_tied(raw), dat$n, wall))

    for (md in MODES) {
      sc <- variants[[md]]
      if (is.null(sc)) next
      v <- evaluate(dat$Y, sc, CUTOFF)
      a <- rank_auc(dat$Y, sc)
      if (md == "legacy") base <- v
      cat(sprintf("      %-13s TPR %.3f TNR %.3f BA %.3f F2 %.3f AUC %.3f | dBA %+0.3f dF2 %+0.3f | tied %d\n",
                  md, v[["TPR"]], v[["TNR"]], v[["BA"]], v[["F2"]], a,
                  v[["BA"]] - base[["BA"]], v[["F2"]] - base[["F2"]], n_tied(sc)))

      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, n = dat$n, d = dat$d, method = meth, mode = md,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]],
        AUC = a,
        dBA_vs_legacy = v[["BA"]] - base[["BA"]],
        dF2_vs_legacy = v[["F2"]] - base[["F2"]],
        n_tied_raw = n_tied(raw), n_tied_after = n_tied(sc),
        n_changed_vs_raw = sum(sc != raw),
        pub_TPR = pub[1], pub_TNR = pub[2], pub_BA = pub[3], pub_F2 = pub[4],
        matches_pub = all(abs(c(v[["TPR"]], v[["TNR"]], v[["BA"]], v[["F2"]]) - pub) <= 0.005 + 1e-9),
        secs = round(wall, 1), stringsAsFactors = FALSE)
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, OUT_CSV, row.names = FALSE)

cat("\n==== how far each mode moves the published-path numbers ====\n")
for (md in setdiff(MODES, "legacy")) {
  s <- out[out$mode == md, ]
  if (!nrow(s)) next
  cat(sprintf("  %-13s max |dBA| %.3f, max |dF2| %.3f, cells changed %d of %d\n",
              md, max(abs(s$dBA_vs_legacy)), max(abs(s$dF2_vs_legacy)),
              sum(s$dBA_vs_legacy != 0 | s$dF2_vs_legacy != 0), nrow(s)))
}
cat("\n---- cells where the tie rule changes BA or F2 ----\n")
ch <- out[out$mode != "legacy" & (out$dBA_vs_legacy != 0 | out$dF2_vs_legacy != 0), ]
if (nrow(ch) == 0) {
  cat("  none\n")
} else {
  for (i in seq_len(nrow(ch))) {
    cat(sprintf("  %-13s %-10s %-13s BA %.3f -> %.3f, F2 %.3f -> %.3f\n",
                ch$dataset[i], ch$method[i], ch$mode[i],
                ch$BA[i] - ch$dBA_vs_legacy[i], ch$BA[i],
                ch$F2[i] - ch$dF2_vs_legacy[i], ch$F2[i]))
  }
}
cat(sprintf("\nWrote %s\n17_tiebreak_impact.R done.\n", OUT_CSV))
