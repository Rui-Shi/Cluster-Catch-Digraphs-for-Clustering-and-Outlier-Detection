#!/usr/bin/env Rscript
# 82_wp4_metrics.R -- turn the WP4 raw scores into TPR/TNR/BA/F2, and rebuild
# the two main-text real-data tables with the enlarged method set.
#
# WP4 (tr1/81_wp4_baselines.py) deliberately computes no metrics: it emits raw
# scores and, where a method has one, a native threshold-free label vector.
# Everything downstream of that lives here, so that the competitors and the
# nine methods already in the study are scored by the SAME evaluate() from
# shared/harness.R.
#
# The protocol this follows is tr1/WP4_PROTOCOL.md, declared before any
# competitor was run. Nothing in it is reinterpreted here.
#
# ---------------------------------------------------------------------------
# LABEL POLARITY -- two opposite conventions, never silently mixed
#
#   this repo (data/<ds>.csv label column, evaluate(), count_scores2)
#       1 = REGULAR, 0 = OUTLIER
#   pyod / sklearn / hdbscan (every is_outlier column WP4 wrote)
#       1 = OUTLIER, 0 = REGULAR
#
# ---------------------------------------------------------------------------
# THRESHOLDING -- pyod's rule, and the > vs >= mismatch with evaluate()
#
# pyod sets
#     threshold_ = np.percentile(decision_scores_, 100 * (1 - contamination))
#     labels_    = (decision_scores_ > threshold_).astype('int')
# np.percentile's default interpolation is linear, which is exactly R's
# quantile(type = 7). So the threshold is reproduced by
#     quantile(score, 1 - contamination, type = 7, names = FALSE).
#
# The comparison is STRICTLY GREATER. evaluate() -> count_scores2() uses
#     label_pred[which(score >= threshold)] = 0   # 0 == outlier
# i.e. GREATER OR EQUAL. The two differ on every point whose score equals the
# threshold, and that is not a corner case here: mutual-kNN scores are
# integers, GLOSH is exactly 0 for most core points, and SNN densities are
# integers, so hundreds of points can sit on the threshold at once.
#
# Resolution: derive the labels here under pyod's own strict rule, then hand
# evaluate() the BINARY indicator with threshold 0.5. `indicator >= 0.5` is
# true exactly when the pyod label is 1, so the confusion counts are pyod's
# while the TPR/TNR/BA/F2 arithmetic is the study's. TIES AT THE THRESHOLD ARE
# THEREFORE CALLED REGULAR, which is pyod's behaviour and is why the flagged
# fraction can fall well below the nominal contamination; n_flagged is written
# to the long file for every cell so this is visible rather than inferred.
#
# The native threshold-free variants (HDBSCAN noise, OPTICS noise,
# mutual-kNN m_i = 0) take the same binary-indicator route with no threshold
# at all.
#
# ---------------------------------------------------------------------------
# WHAT IS REBUILT, AND THE GATE IN FRONT OF IT
#
# The two printed tables cannot be patched cell by cell: the summary table has
# three DERIVED columns (best overall, best proposed, SUN-MCCD's rank) and the
# aggregate table is ORDERED by mean F2, so adding eight methods can move rows
# whose own numbers did not change. Both are rebuilt from per-cell values,
# following 78_manuscript_tables_after_fix.R exactly -- same sprintf("%.3f")
# rounding (the manuscript rounds half up, R's round() rounds half to even,
# and SUN-MCCD's median F2 is exactly 0.3285), same stale-cell patch, same
# rank/winner computation on ROUNDED values, same tie handling.
#
# Before any of that is believed, the nine-method table as published must
# reproduce the 36 values printed in tab:Real_Data_Aggregate. That gate runs
# first and the script stops if it fails.
#
# This script writes only its own CSVs and its own findings file.

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments/shared/harness.R")))
options(width = 200)   # the comparison tables are wide; do not let print() wrap

WP4  <- here::here("revision_experiments/results/tr1/wp4")
DATA <- file.path(WP4, "data")
SCO  <- file.path(WP4, "scores")
FC   <- here::here("revision_experiments/results/tr1/final_comparison.csv")
NEW  <- here::here("revision_experiments/results/tr1/wp2c_rerun_alpha001.csv")

OUT_LONG <- file.path(WP4, "wp4_metrics_long.csv")
OUT_MAIN <- file.path(WP4, "wp4_metrics_main.csv")
OUT_MD   <- file.path(WP4, "WP4_FINDINGS.md")

ORDER <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
           "stamps", "WDBC", "pima", "Shuttle", "vowels", "PenDigits",
           "waveform", "thyroid", "pageblocks", "wilt")
PROPOSED <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
KS <- c(5, 10, 15, 20, 30)
SEEDS <- 1:5

# manuscript rounding: half UP, not R's half-to-even
r3 <- function(x) as.numeric(sprintf("%.3f", x))

# capture everything printed, so the findings file and the console agree
CON <- textConnection("REPORT_LINES", "w", local = FALSE)
say <- function(...) { txt <- paste0(...); cat(txt)
                       writeLines(sub("\n$", "", txt), CON) }
sayf <- function(fmt, ...) say(sprintf(fmt, ...))
saydf <- function(df) {
  o <- capture.output(print(df, row.names = FALSE))
  cat(paste0(o, collapse = "\n"), "\n", sep = ""); writeLines(o, CON)
}

# ---------------------------------------------------------------------------
# 0. inputs
# ---------------------------------------------------------------------------
man <- read.csv(file.path(DATA, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(all(man$table_match), setequal(man$dataset, ORDER))

labs <- setNames(lapply(ORDER, function(ds) {
  d <- read.csv(file.path(DATA, paste0(ds, ".csv")), stringsAsFactors = FALSE)
  y <- d[[ncol(d)]]
  stopifnot(names(d)[ncol(d)] == "label", all(y %in% c(0, 1)))
  m <- man[man$dataset == ds, ]
  stopifnot(length(y) == m$n, sum(y == 0) == m$n_outliers)
  y
}), ORDER)

cont_true <- setNames(man$contamination[match(ORDER, man$dataset)], ORDER)

read_score <- function(ds, tag) {
  f <- file.path(SCO, sprintf("%s_%s.csv", ds, tag))
  if (!file.exists(f)) stop("missing score file: ", f)
  s <- read.csv(f, stringsAsFactors = FALSE)$score
  if (length(s) != length(labs[[ds]]))
    stop(sprintf("length mismatch %s_%s: %d scores, %d labels", ds, tag,
                 length(s), length(labs[[ds]])))
  if (any(!is.finite(s))) stop("non-finite score in ", f)
  s
}
read_native <- function(ds, tag) {
  f <- file.path(SCO, sprintf("%s_%s_labels.csv", ds, tag))
  if (!file.exists(f)) stop("missing label file: ", f)
  v <- read.csv(f, stringsAsFactors = FALSE)$is_outlier
  stopifnot(length(v) == length(labs[[ds]]), all(v %in% c(0, 1)))
  v
}

# pyod thresholding, then evaluate() on the binary indicator (see header)
metrics_from_score <- function(ds, score, contamination) {
  thr <- quantile(score, 1 - contamination, type = 7, names = FALSE)
  is_out <- as.numeric(score > thr)              # STRICT >, pyod's rule
  v <- evaluate(labs[[ds]], is_out, 0.5)
  c(v, n_flagged = sum(is_out), n_tied_at_thr = sum(score == thr),
    threshold = thr)
}
metrics_from_labels <- function(ds, is_out) {
  v <- evaluate(labs[[ds]], as.numeric(is_out), 0.5)
  c(v, n_flagged = sum(is_out), n_tied_at_thr = NA_real_, threshold = NA_real_)
}

# ---------------------------------------------------------------------------
# 1-3. every cell, every regime, every seed / k
# ---------------------------------------------------------------------------
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(...,
                                                    stringsAsFactors = FALSE)

emit <- function(ds, method, tag, k = NA, seed = NA) {
  s <- read_score(ds, tag)
  for (rg in c("T1", "T2")) {
    cn <- if (rg == "T1") 0.10 else cont_true[[ds]]
    v <- metrics_from_score(ds, s, cn)
    add(dataset = ds, method = method, variant = "score", regime = rg,
        contamination = cn, k = k, seed = seed,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]],
        n_flagged = v[["n_flagged"]], n_tied_at_thr = v[["n_tied_at_thr"]],
        threshold = v[["threshold"]])
  }
}

for (ds in ORDER) {
  emit(ds, "ECOD",  "ECOD")
  emit(ds, "COPOD", "COPOD")
  for (s in SEEDS) emit(ds, "DIF",   sprintf("DIF_seed%d", s),   seed = s)
  for (s in SEEDS) emit(ds, "LUNAR", sprintf("LUNAR_seed%d", s), seed = s)
  emit(ds, "GLOSH",  "HDBSCAN")            # hdbscan outlier_scores_
  emit(ds, "OPTICS", "OPTICS")             # reachability_, inf filled 1.01x
  for (k in KS) emit(ds, "MutualKNN", sprintf("MutualKNN_k%d", k), k = k)
  for (k in KS) emit(ds, "SNN",       sprintf("SNN_k%d", k),       k = k)

  # threshold-free native variants
  nat <- list(c("HDBSCAN-noise", "HDBSCAN"), c("OPTICS-noise", "OPTICS"))
  for (p in nat) {
    v <- metrics_from_labels(ds, read_native(ds, p[2]))
    add(dataset = ds, method = p[1], variant = "native", regime = "native",
        contamination = NA, k = NA, seed = NA,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]],
        n_flagged = v[["n_flagged"]], n_tied_at_thr = NA, threshold = NA)
  }
  for (k in KS) {
    v <- metrics_from_labels(ds, read_native(ds, sprintf("MutualKNN_k%d", k)))
    add(dataset = ds, method = "MutualKNN-m0", variant = "native",
        regime = "native", contamination = NA, k = k, seed = NA,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]],
        n_flagged = v[["n_flagged"]], n_tied_at_thr = NA, threshold = NA)
  }
}
long <- do.call(rbind, rows)
write.csv(long, OUT_LONG, row.names = FALSE)

# --- seeded methods: mean and sd over the 5 seeds; the MEAN is reported ------
seed_summary <- do.call(rbind, lapply(c("DIF", "LUNAR"), function(m) {
  do.call(rbind, lapply(c("T1", "T2"), function(rg) {
    g <- long[long$method == m & long$regime == rg, ]
    do.call(rbind, lapply(ORDER, function(ds) {
      h <- g[g$dataset == ds, ]; stopifnot(nrow(h) == 5)
      data.frame(dataset = ds, method = m, regime = rg,
                 TPR = mean(h$TPR), TNR = mean(h$TNR),
                 BA = mean(h$BA), F2 = mean(h$F2),
                 sd_TPR = sd(h$TPR), sd_TNR = sd(h$TNR),
                 sd_BA = sd(h$BA), sd_F2 = sd(h$F2),
                 k = NA_real_, stringsAsFactors = FALSE)
    }))
  }))
}))

# --- k-swept methods: oracle-best k by F2 (declared), and fixed k = 10 -------
pick_k <- function(m, rg, mode) {
  do.call(rbind, lapply(ORDER, function(ds) {
    h <- long[long$method == m & long$regime == rg & long$dataset == ds, ]
    stopifnot(nrow(h) == length(KS))
    i <- if (mode == "oracle") which.max(h$F2) else which(h$k == 10)
    data.frame(dataset = ds, method = m, regime = rg,
               TPR = h$TPR[i], TNR = h$TNR[i], BA = h$BA[i], F2 = h$F2[i],
               sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
               k = h$k[i], stringsAsFactors = FALSE)
  }))
}
# which.max takes the FIRST maximum, i.e. the smallest k among F2 ties.
oracle_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) pick_k(m, rg, "oracle")))))
fixed_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) pick_k(m, rg, "fixed")))))
fixed_k$method <- paste0(fixed_k$method, "-k10")

plain <- do.call(rbind, lapply(c("ECOD", "COPOD", "GLOSH", "OPTICS"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) {
    g <- long[long$method == m & long$regime == rg, ]
    g <- g[match(ORDER, g$dataset), ]
    data.frame(dataset = g$dataset, method = m, regime = rg,
               TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
               sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
               k = NA_real_, stringsAsFactors = FALSE)
  }))))

native_rows <- do.call(rbind, lapply(c("HDBSCAN-noise", "OPTICS-noise"), function(m) {
  g <- long[long$method == m, ]; g <- g[match(ORDER, g$dataset), ]
  data.frame(dataset = g$dataset, method = m, regime = "native",
             TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
             sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
             k = NA_real_, stringsAsFactors = FALSE)
}))
# mutual-kNN m_i = 0 also sweeps k; report it oracle-k on F2 for symmetry
native_rows <- rbind(native_rows, do.call(rbind, lapply(ORDER, function(ds) {
  h <- long[long$method == "MutualKNN-m0" & long$dataset == ds, ]
  i <- which.max(h$F2)
  data.frame(dataset = ds, method = "MutualKNN-m0", regime = "native",
             TPR = h$TPR[i], TNR = h$TNR[i], BA = h$BA[i], F2 = h$F2[i],
             sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
             k = h$k[i], stringsAsFactors = FALSE)
})))

main <- rbind(plain, seed_summary, oracle_k, fixed_k, native_rows)
main <- main[order(match(main$dataset, ORDER), main$method, main$regime), ]
write.csv(main, OUT_MAIN, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. the corrected nine-method table (78's logic, verbatim in effect)
# ---------------------------------------------------------------------------
fc <- read.csv(FC, stringsAsFactors = FALSE)
nw <- read.csv(NEW, stringsAsFactors = FALSE); nw <- nw[nw$status == "ok", ]

# NAME. final_comparison.csv calls the Shuttle set by the loader's R object
# name, "shuffle". 78_manuscript_tables_after_fix.R matched data set names
# against the manuscript's ORDER vector and returned NULL when nothing matched,
# so its rebuilt tab:Real_Data_Result_Summary silently carried 15 rows, not 16
# -- Shuttle was dropped from it. The aggregate table was unaffected (it splits
# by method, not by data set), which is why the 36-value gate still passed.
# Renamed here so the summary table is complete.
stopifnot(sum(fc$dataset == "shuffle") == 9)
fc$dataset[fc$dataset == "shuffle"] <- "Shuttle"

# STALE CELL. final_comparison.csv still holds vertebral/SU-MCCD at the old
# S_min = 0.0625 operating point (BA 0.388, F2 0.140). The manuscript already
# prints the S_min = 0.05 values.
kk <- which(tolower(fc$dataset) == "vertebral" & fc$method == "SU-MCCD")
stopifnot(length(kk) == 1)
fc$BA[kk] <- 0.493; fc$F2[kk] <- 0.129

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
fx <- apply_fix(fc)   # the corrected nine-method table (hepatitis/SUN included)

agg <- function(df) {
  a <- do.call(rbind, lapply(split(df, df$method), function(g) data.frame(
    method = g$method[1], n = nrow(g),
    mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
    mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
    stringsAsFactors = FALSE)))
  a[order(-a$mean_F2), ]
}

# --- GATE: does the PUBLISHED nine-method table reproduce what is printed? ---
PRINTED <- data.frame(
  method    = c("SUN-MCCD","SU-MCCD","U-MCCD","iForest","LOF","MST","DBSCAN","ODIN","UN-MCCD"),
  mean_F2   = c(0.321, 0.302, 0.291, 0.260, 0.254, 0.237, 0.234, 0.229, 0.229),
  median_F2 = c(0.329, 0.301, 0.280, 0.115, 0.195, 0.237, 0.133, 0.224, 0.209),
  mean_BA   = c(0.722, 0.696, 0.709, 0.632, 0.654, 0.606, 0.599, 0.619, 0.658),
  median_BA = c(0.726, 0.692, 0.708, 0.526, 0.577, 0.599, 0.546, 0.585, 0.634),
  stringsAsFactors = FALSE)
a_pub <- agg(fc)
say("=== GATE: does the published nine-method table reproduce tab:Real_Data_Aggregate? ===\n")
mm <- merge(PRINTED, a_pub, by = "method", suffixes = c("_printed", "_computed"))
bad <- 0L
for (col in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  d <- abs(mm[[paste0(col, "_printed")]] - mm[[paste0(col, "_computed")]])
  for (i in which(d > 0.0005)) {
    bad <- bad + 1L
    sayf("  MISMATCH %-9s %-10s printed %.3f, computed %.3f\n", mm$method[i],
         col, mm[[paste0(col, "_printed")]][i], mm[[paste0(col, "_computed")]][i])
  }
}
if (bad == 0) say("  all 36 printed values reproduced exactly\n") else
  stop(sprintf("%d printed values do not reproduce -- nothing downstream is trusted", bad))

# ---------------------------------------------------------------------------
# 5. BEFORE / AFTER
# ---------------------------------------------------------------------------
new_cells <- function(rg) {
  keep <- c("ECOD", "COPOD", "DIF", "LUNAR", "GLOSH", "OPTICS",
            "MutualKNN", "SNN")
  g <- main[main$regime == rg & main$method %in% keep, ]
  data.frame(dataset = g$dataset, method = g$method,
             F2 = r3(g$F2), BA = r3(g$BA), TPR = r3(g$TPR), TNR = r3(g$TNR),
             group = "wp4", stringsAsFactors = FALSE)
}
combined <- function(rg) rbind(fx[, c("dataset","method","F2","BA","TPR","TNR","group")],
                               new_cells(rg))

summ <- function(df) {
  do.call(rbind, lapply(ORDER, function(ds) {
    g <- df[tolower(df$dataset) == tolower(ds), ]
    stopifnot(nrow(g) > 0, "SUN-MCCD" %in% g$method)
    f <- round(g$F2, 3); names(f) <- g$method
    bo <- names(f)[f == max(f)]
    p  <- f[names(f) %in% PROPOSED]
    bp <- names(p)[p == max(p)]
    data.frame(dataset = ds,
               best_overall = paste(bo, collapse = "/"), f2_overall = max(f),
               best_proposed = paste(bp, collapse = "/"), f2_proposed = max(p),
               sun_rank = sum(f > f[["SUN-MCCD"]]) + 1, n_methods = length(f),
               stringsAsFactors = FALSE)
  }))
}

s_before <- summ(fx); a_before <- agg(fx)
c1 <- combined("T1"); s_after <- summ(c1); a_after <- agg(c1)
c2 <- combined("T2"); s_after2 <- summ(c2); a_after2 <- agg(c2)

say("\n=== tab:Real_Data_Aggregate -- BEFORE (9 methods, alpha-corrected) ===\n")
saydf(a_before[, -2])
say("\n=== tab:Real_Data_Aggregate -- AFTER (9 + 8 WP4 baselines), T1 = contamination 0.1 ===\n")
saydf(a_after[, -2])
say("\n=== tab:Real_Data_Aggregate -- T2 = ORACLE regime (contamination = true rate) ===\n")
say("    reported as a sensitivity only; the true outlier rate is not available to an\n")
say("    unsupervised detector, and the nine incumbent rows are NOT recomputed at T2.\n")
saydf(a_after2[, -2])

verdict <- function(a, tag) {
  say(sprintf("\n--- SUN-MCCD's standing, %s ---\n", tag))
  for (m in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
    v <- a[[m]]; names(v) <- a$method
    pos <- sum(v > v[["SUN-MCCD"]]) + 1
    best_other <- max(v[names(v) != "SUN-MCCD"])
    who <- paste(names(v)[names(v) != "SUN-MCCD" & v == best_other], collapse = "/")
    sayf("  %-10s rank %2d of %2d | SUN-MCCD %.3f vs best other %.3f (%s) | margin %+0.3f\n",
         m, pos, length(v), v[["SUN-MCCD"]], best_other, who, v[["SUN-MCCD"]] - best_other)
  }
}
verdict(a_before, "BEFORE, 9 methods"); verdict(a_after, "AFTER, T1")
verdict(a_after2, "AFTER, T2 oracle")

cnt <- function(s) c(
  prop_wins = sum(sapply(strsplit(s$best_overall, "/"), function(x) any(x %in% PROPOSED))),
  sun_wins  = sum(sapply(strsplit(s$best_overall, "/"), function(x) "SUN-MCCD" %in% x)),
  sun_sole  = sum(s$best_overall == "SUN-MCCD"))
cb <- cnt(s_before); ca <- cnt(s_after); ca2 <- cnt(s_after2)
say("\n=== outright wins on F2, per data set ===\n")
sayf("  a PROPOSED method attains the highest F2 : %2d of 16 (before) -> %2d of 16 (T1) -> %2d of 16 (T2)\n", cb[1], ca[1], ca2[1])
sayf("  SUN-MCCD attains the highest F2 (incl. ties): %2d of 16 (before) -> %2d of 16 (T1) -> %2d of 16 (T2)\n", cb[2], ca[2], ca2[2])
sayf("  SUN-MCCD sole winner                        : %2d of 16 (before) -> %2d of 16 (T1) -> %2d of 16 (T2)\n", cb[3], ca[3], ca2[3])

say("\n=== tab:Real_Data_Result_Summary -- BEFORE (9) vs AFTER (17), T1 ===\n")
chg <- which(s_before$best_overall != s_after$best_overall |
             s_before$f2_overall  != s_after$f2_overall |
             s_before$sun_rank    != s_after$sun_rank)
for (i in seq_len(nrow(s_before))) {
  mark <- if (i %in% chg) "*" else " "
  sayf("%s %-13s BEFORE %-20s %.3f rank %d/%d | AFTER %-20s %.3f rank %2d/%2d\n",
       mark, s_before$dataset[i], s_before$best_overall[i], s_before$f2_overall[i],
       s_before$sun_rank[i], s_before$n_methods[i],
       s_after$best_overall[i], s_after$f2_overall[i], s_after$sun_rank[i], s_after$n_methods[i])
}
sayf("  rows whose best-overall method or SUN-MCCD rank changes: %d of 16\n", length(chg))

say("\n=== rebuilt tab:Real_Data_Result_Summary (all 16 rows, AFTER / T1) ===\n")
saydf(s_after[, setdiff(names(s_after), "n_methods")])

# ---------------------------------------------------------------------------
# 6. R3.3 -- mutual reciprocity vs mutual CATCH
# ---------------------------------------------------------------------------
say("\n=== R3.3: oracle-k mutual-kNN vs the two NND-based proposed methods (T1) ===\n")
mk  <- main[main$method == "MutualKNN" & main$regime == "T1", ]
mk10 <- main[main$method == "MutualKNN-k10" & main$regime == "T1", ]
mk0 <- main[main$method == "MutualKNN-m0", ]
g <- function(ds, m, col) fx[[col]][tolower(fx$dataset) == tolower(ds) & fx$method == m]
mkt <- do.call(rbind, lapply(ORDER, function(ds) data.frame(
  dataset = ds,
  k_star  = mk$k[mk$dataset == ds],
  mkNN_F2 = r3(mk$F2[mk$dataset == ds]), mkNN_BA = r3(mk$BA[mk$dataset == ds]),
  mkNN_k10_F2 = r3(mk10$F2[mk10$dataset == ds]),
  mkNN_m0_F2  = r3(mk0$F2[mk0$dataset == ds]),
  SUN_F2 = g(ds, "SUN-MCCD", "F2"), SUN_BA = g(ds, "SUN-MCCD", "BA"),
  UN_F2  = g(ds, "UN-MCCD",  "F2"), UN_BA  = g(ds, "UN-MCCD",  "BA"),
  stringsAsFactors = FALSE)))
saydf(mkt)
sayf("  SUN-MCCD F2 > oracle-k mutual-kNN F2 on %d of 16; mutual-kNN ahead on %d; tied %d\n",
     sum(mkt$SUN_F2 > mkt$mkNN_F2), sum(mkt$mkNN_F2 > mkt$SUN_F2), sum(mkt$SUN_F2 == mkt$mkNN_F2))
sayf("  SUN-MCCD BA > oracle-k mutual-kNN BA on %d of 16; mutual-kNN ahead on %d; tied %d\n",
     sum(mkt$SUN_BA > mkt$mkNN_BA), sum(mkt$mkNN_BA > mkt$SUN_BA), sum(mkt$SUN_BA == mkt$mkNN_BA))
sayf("  UN-MCCD  F2 > oracle-k mutual-kNN F2 on %d of 16; mutual-kNN ahead on %d; tied %d\n",
     sum(mkt$UN_F2 > mkt$mkNN_F2), sum(mkt$mkNN_F2 > mkt$UN_F2), sum(mkt$UN_F2 == mkt$mkNN_F2))
sayf("  mean F2: mutual-kNN oracle-k %.3f, fixed k=10 %.3f, m_i=0 native %.3f, SUN-MCCD %.3f, UN-MCCD %.3f\n",
     r3(mean(mkt$mkNN_F2)), r3(mean(mkt$mkNN_k10_F2)), r3(mean(mkt$mkNN_m0_F2)),
     r3(mean(mkt$SUN_F2)), r3(mean(mkt$UN_F2)))
sayf("  mean BA: mutual-kNN oracle-k %.3f, SUN-MCCD %.3f, UN-MCCD %.3f\n",
     r3(mean(mkt$mkNN_BA)), r3(mean(mkt$SUN_BA)), r3(mean(mkt$UN_BA)))
sayf("  oracle-k gain over fixed k = 10, mean F2: %+0.3f (the oracle's own advantage)\n",
     r3(mean(mkt$mkNN_F2)) - r3(mean(mkt$mkNN_k10_F2)))

say("\n=== SNN: oracle-k vs fixed k = 10 (T1) ===\n")
sn <- main[main$method == "SNN" & main$regime == "T1", ]
sn10 <- main[main$method == "SNN-k10" & main$regime == "T1", ]
snt <- data.frame(dataset = ORDER, k_star = sn$k[match(ORDER, sn$dataset)],
                  SNN_F2 = r3(sn$F2[match(ORDER, sn$dataset)]),
                  SNN_k10_F2 = r3(sn10$F2[match(ORDER, sn10$dataset)]),
                  SNN_BA = r3(sn$BA[match(ORDER, sn$dataset)]),
                  stringsAsFactors = FALSE)
saydf(snt)
sayf("  mean F2: SNN oracle-k %.3f, SNN fixed k=10 %.3f\n",
     r3(mean(snt$SNN_F2)), r3(mean(snt$SNN_k10_F2)))

# ---------------------------------------------------------------------------
# 7. R1.3 -- density clustering with varying-density support
# ---------------------------------------------------------------------------
say("\n=== R1.3: HDBSCAN/GLOSH and OPTICS vs the two shape-adaptive proposed methods (T1) ===\n")
fl <- read.csv(file.path(WP4, "fit_log.csv"), stringsAsFactors = FALSE)
nc <- function(ds) {
  n <- fl$note[fl$dataset == ds & fl$method == "HDBSCAN"][1]
  as.integer(sub(".*n_clusters=([0-9]+).*", "\\1", n))
}
gl <- main[main$method == "GLOSH" & main$regime == "T1", ]
hn <- main[main$method == "HDBSCAN-noise", ]
op <- main[main$method == "OPTICS" & main$regime == "T1", ]
hb <- do.call(rbind, lapply(ORDER, function(ds) data.frame(
  dataset = ds, hdb_clusters = nc(ds),
  GLOSH_F2 = r3(gl$F2[gl$dataset == ds]), GLOSH_BA = r3(gl$BA[gl$dataset == ds]),
  HDBnoise_F2 = r3(hn$F2[hn$dataset == ds]), HDBnoise_BA = r3(hn$BA[hn$dataset == ds]),
  OPTICS_F2 = r3(op$F2[op$dataset == ds]),
  SU_F2 = g(ds, "SU-MCCD", "F2"), SU_BA = g(ds, "SU-MCCD", "BA"),
  SUN_F2 = g(ds, "SUN-MCCD", "F2"), SUN_BA = g(ds, "SUN-MCCD", "BA"),
  flag = ifelse(nc(ds) == 0, "HDBSCAN found NO cluster: every point noise", ""),
  stringsAsFactors = FALSE)))
saydf(hb)
sayf("  SUN-MCCD F2 > GLOSH F2 on %d of 16; GLOSH ahead on %d; tied %d\n",
     sum(hb$SUN_F2 > hb$GLOSH_F2), sum(hb$GLOSH_F2 > hb$SUN_F2), sum(hb$SUN_F2 == hb$GLOSH_F2))
sayf("  SU-MCCD  F2 > GLOSH F2 on %d of 16; GLOSH ahead on %d; tied %d\n",
     sum(hb$SU_F2 > hb$GLOSH_F2), sum(hb$GLOSH_F2 > hb$SU_F2), sum(hb$SU_F2 == hb$GLOSH_F2))
sayf("  mean F2: GLOSH %.3f, HDBSCAN-noise %.3f, OPTICS %.3f, SU-MCCD %.3f, SUN-MCCD %.3f\n",
     r3(mean(hb$GLOSH_F2)), r3(mean(hb$HDBnoise_F2)), r3(mean(hb$OPTICS_F2)),
     r3(mean(hb$SU_F2)), r3(mean(hb$SUN_F2)))
sayf("  mean BA: GLOSH %.3f, HDBSCAN-noise %.3f, SU-MCCD %.3f, SUN-MCCD %.3f\n",
     r3(mean(hb$GLOSH_BA)), r3(mean(hb$HDBnoise_BA)), r3(mean(hb$SU_BA)), r3(mean(hb$SUN_BA)))
sayf("  HDBSCAN found no cluster on: %s (its native labels flag 100%% of the data there)\n",
     paste(hb$dataset[hb$hdb_clusters == 0], collapse = ", "))

say("\n=== threshold-free native variants, aggregate (no contamination argument) ===\n")
nat_agg <- do.call(rbind, lapply(c("HDBSCAN-noise", "OPTICS-noise", "MutualKNN-m0"),
  function(m) { g <- main[main$method == m & main$regime == "native", ]
    data.frame(method = m, mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
               mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
               stringsAsFactors = FALSE) }))
saydf(nat_agg[order(-nat_agg$mean_F2), ])
say("  these are NOT in the AFTER aggregate table above: they take no contamination\n")
say("  argument, so they are not on the T1 footing the other 17 rows share.\n")

# ---------------------------------------------------------------------------
# 2 (report). seed stability
# ---------------------------------------------------------------------------
say("\n=== seed stability, 5 seeds, sd across seeds per data set (T1) ===\n")
for (m in c("DIF", "LUNAR")) {
  h <- seed_summary[seed_summary$method == m & seed_summary$regime == "T1", ]
  sayf("  %-6s sd(F2) min %.4f max %.4f mean %.4f (max on %s) | sd(BA) min %.4f max %.4f mean %.4f\n",
       m, min(h$sd_F2), max(h$sd_F2), mean(h$sd_F2), h$dataset[which.max(h$sd_F2)],
       min(h$sd_BA), max(h$sd_BA), mean(h$sd_BA))
}
say("  per-data-set mean and sd over the 5 seeds (T1):\n")
ss <- function(m, col) { h <- seed_summary[seed_summary$method == m & seed_summary$regime == "T1", ]
                         h[[col]][match(ORDER, h$dataset)] }
saydf(data.frame(dataset = ORDER,
                 DIF_meanF2 = r3(ss("DIF", "F2")),     DIF_sdF2 = round(ss("DIF", "sd_F2"), 4),
                 DIF_meanBA = r3(ss("DIF", "BA")),     DIF_sdBA = round(ss("DIF", "sd_BA"), 4),
                 LUNAR_meanF2 = r3(ss("LUNAR", "F2")), LUNAR_sdF2 = round(ss("LUNAR", "sd_F2"), 4),
                 LUNAR_meanBA = r3(ss("LUNAR", "BA")), LUNAR_sdBA = round(ss("LUNAR", "sd_BA"), 4),
                 stringsAsFactors = FALSE))

# ---------------------------------------------------------------------------
# 8. findings file
# ---------------------------------------------------------------------------
headline_ok <- function(a) {
  sapply(c("mean_F2","median_F2","mean_BA","median_BA"), function(m) {
    v <- a[[m]]; names(v) <- a$method; v[["SUN-MCCD"]] == max(v) })
}
h1 <- headline_ok(a_after)
say("\n=== headline claim ===\n")
say("  claim as published: \"SUN-MCCD ranks first on mean and median BA and mean and\n")
say("  median F2 across the 16 real data sets, ahead of all five baselines.\"\n")
sayf("  against the enlarged set (9 + 8 = %d methods, T1): first on %d of the 4 aggregates\n",
     nrow(a_after), sum(h1))
for (m in names(h1)) {
  v <- a_after[[m]]; names(v) <- a_after$method
  bo <- max(v[names(v) != "SUN-MCCD"])
  who <- paste(names(v)[names(v) != "SUN-MCCD" & v == bo], collapse = "/")
  sayf("    %-10s %s (SUN-MCCD %.3f, best other %.3f = %s, %+0.3f)\n", m,
       ifelse(h1[[m]], "HOLDS ", "LOST  "), v[["SUN-MCCD"]], bo, who, v[["SUN-MCCD"]] - bo)
}

say("\n=== what the numbers do and do not license ===\n")
say("  1. The published claim does not survive the enlarged comparison at T1. SUN-MCCD\n")
say("     is fourth on mean F2 and second on the other three aggregates, and every\n")
say("     position it loses, it loses to OPTICS run as a SCORE (reachability at\n")
say("     contamination 0.1), not to a deep detector.\n")
say("  2. OPTICS is not uniformly better: its own threshold-free variant, the one an\n")
say("     analyst gets without being told the contamination rate, averages F2 0.186 and\n")
say("     BA 0.476 -- below every proposed method. The strong row above is the\n")
say("     reachability score cut at a contamination the detector was handed.\n")
say("  3. The enlarged set moves SUN-MCCD from a first place to a top-quartile place;\n")
say("     it does not overturn the paper's internal comparison. SUN-MCCD still leads\n")
say("     the four proposed detectors on all four aggregates, and UN-MCCD -> SUN-MCCD,\n")
say("     the shape-adaptive step the paper argues for, is unchanged.\n")
say("  4. R3.3 has an answer that does not depend on how k was chosen: an oracle-k\n")
say("     mutual-kNN detector, given the true labels to pick k, still averages below\n")
say("     SUN-MCCD on both F2 (0.237 vs 0.317) and BA (0.628 vs 0.718), and loses on\n")
say("     11 of 16 sets on each. Reciprocity alone does not reproduce the gain. The\n")
say("     same comparison against UN-MCCD goes the OTHER way -- oracle-k mutual-kNN is\n")
say("     ahead on 9 of 16 and on mean F2 (0.237 vs 0.229) -- so what beats plain\n")
say("     reciprocity is the shape-adaptive coverage, not the NND test on its own.\n")
say("  5. The nine incumbent rows are NOT recomputed at T2. The T2 block therefore\n")
say("     compares oracle-thresholded competitors against fixed-threshold incumbents\n")
say("     and is a sensitivity on the competitors only, not a like-for-like table.\n")

close(CON)
md <- c("# WP4 findings -- modern baselines on the sixteen real data sets",
        "",
        "Generated by `tr1/82_wp4_metrics.R` from the raw scores of",
        "`tr1/81_wp4_baselines.py`, under the pre-declared `tr1/WP4_PROTOCOL.md`.",
        "Metrics come from the same `evaluate()` in `shared/harness.R` that scored the",
        "nine methods already in the study. Thresholds follow pyod exactly:",
        "`threshold = quantile(score, 1 - contamination, type = 7)` and",
        "`is_outlier = score > threshold` (strictly greater, so ties at the threshold",
        "are called regular). T1 = contamination 0.1, the main-table setting; T2 =",
        "the data set's true outlier rate, an oracle reported only as a sensitivity.",
        "k for mutual-kNN and SNN is chosen per data set by F2 on the true labels --",
        "an oracle concession to the competitor, declared in advance.",
        "", "```", REPORT_LINES, "```", "")
writeLines(md, OUT_MD)

sayf2 <- function(...) cat(sprintf(...))
sayf2("\nwrote %s\n", OUT_LONG)
sayf2("wrote %s\n", OUT_MAIN)
sayf2("wrote %s\n", OUT_MD)
