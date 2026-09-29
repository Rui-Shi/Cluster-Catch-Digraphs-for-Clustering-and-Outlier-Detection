#!/usr/bin/env Rscript
# 83_wp4_latex_rows.R -- print the LaTeX table rows for the WP4 write-up.
#
# Every number that goes into the manuscript or the supplement for WP4 is
# emitted here and pasted, so that no value is typed by hand. The inputs are
# the CSVs written by tr1/82_wp4_metrics.R (read-only) plus the two files that
# hold the nine incumbent methods, reconstructed exactly as 82 reconstructs
# them (Shuttle/shuffle rename, the stale vertebral/SU-MCCD cell, and the
# alpha-corrected rerun).
#
# Rounding is sprintf("%.3f") throughout, matching 82 and the manuscript
# (round half up), not R's round() (round half to even).
#
# This script writes nothing; it prints.

suppressMessages(library(here))
options(width = 250)

WP4  <- here::here("revision_experiments/results/tr1/wp4")
FC   <- here::here("revision_experiments/results/tr1/final_comparison.csv")
NEW  <- here::here("revision_experiments/results/tr1/wp2c_rerun_alpha001.csv")

ORDER <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
           "stamps", "WDBC", "pima", "Shuttle", "vowels", "PenDigits",
           "waveform", "thyroid", "pageblocks", "wilt")
PROPOSED <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
NEWM <- c("ECOD", "COPOD", "DIF", "LUNAR", "GLOSH", "OPTICS", "MutualKNN", "SNN")
INCUMBENT <- c(PROPOSED, "LOF", "DBSCAN", "MST", "ODIN", "iForest")
# block 1 = takes no contamination fraction as input; block 2 = is given one
BLOCK1 <- c(PROPOSED, "LOF", "iForest", "MST", "ODIN")

r3 <- function(x) as.numeric(sprintf("%.3f", x))
f3 <- function(x) sprintf("%.3f", x)

# WHY THE PER-CELL METRICS ARE RECOMPUTED RATHER THAN READ
#
# Both CSVs written by 82_wp4_metrics.R store 15 significant digits, which is
# not a lossless round trip for a double. That is invisible everywhere except
# at an exact half: DIF on PenDigits has mean F2 exactly 0.0725 and LUNAR on
# PenDigits exactly 0.1425, and a change in the last bit moves sprintf("%.3f")
# by 0.001. Reading either CSV back therefore fails to reproduce 82's printed
# third decimal on those two cells.
#
# The fix is to rebuild the doubles rather than parse them. TPR and TNR are
# ratios of small integers, so the integer counts are recovered exactly by
# rounding, and BA and F2 are then recomputed by the arithmetic of
# count_scores2() in R/general_functions/count.R, verbatim:
#     precision = n0*TPR / (n0*TPR + (1 - TNR)*(n - n0))
#     F2        = 5*precision*recall / (4*precision + recall)
# This reproduces 82 bit for bit. Both CSVs are still read, as gates.
long <- read.csv(file.path(WP4, "wp4_metrics_long.csv"), stringsAsFactors = FALSE)
KS <- c(5, 10, 15, 20, 30)
mcols <- c("TPR", "TNR", "BA", "F2")

man <- read.csv(file.path(WP4, "data/manifest.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(man$dataset, ORDER), all(man$table_match))
nn <- setNames(man$n, man$dataset); n0 <- setNames(man$n_outliers, man$dataset)
KSQ <- round(sqrt(nn))   # label-free k per data set (ODIN's rule), as in 82

exact <- function(df) {
  n <- nn[df$dataset]; k <- n0[df$dataset]
  tp <- round(df$TPR * k); tn <- round(df$TNR * (n - k))
  stopifnot(max(abs(df$TPR * k - tp)) < 1e-6, max(abs(df$TNR * (n - k) - tn)) < 1e-6)
  TPR <- tp / k; TNR <- tn / (n - k)
  prec <- k * TPR / (k * TPR + (1 - TNR) * (n - k))
  F2 <- 5 * prec * TPR / (4 * prec + TPR); F2[is.na(F2)] <- 0
  df$TPR <- TPR; df$TNR <- TNR; df$BA <- (TPR + TNR) / 2; df$F2 <- F2
  df
}
raw <- long
long <- exact(long)
for (cl in mcols) stopifnot(max(abs(long[[cl]] - raw[[cl]])) < 1e-9)
cat("%% gate: recomputed per-cell metrics agree with wp4_metrics_long.csv to < 1e-9\n")

seed_summary <- do.call(rbind, lapply(c("DIF", "LUNAR"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg)
    do.call(rbind, lapply(ORDER, function(ds) {
      h <- long[long$method == m & long$regime == rg & long$dataset == ds, ]
      stopifnot(nrow(h) == 5)
      data.frame(dataset = ds, method = m, regime = rg,
                 TPR = mean(h$TPR), TNR = mean(h$TNR),
                 BA = mean(h$BA), F2 = mean(h$F2),
                 sd_TPR = sd(h$TPR), sd_TNR = sd(h$TNR),
                 sd_BA = sd(h$BA), sd_F2 = sd(h$F2), k = NA_real_,
                 stringsAsFactors = FALSE)
    }))))))

# --- k-swept methods (2026-09-26): the PRIMARY rows use the label-free rule
# k = round(sqrt(n)), the rule ODIN uses in this study; the label-chosen
# (oracle-best-F2) k over the declared grid KS is kept as an upper bound under
# the method name "<m>-oracle"; the fixed k = 10 rows stay as "<m>-k10".
pick_k <- function(m, rg, mode) {
  do.call(rbind, lapply(ORDER, function(ds) {
    h <- long[long$method == m & long$regime == rg & long$dataset == ds, ]
    stopifnot(all(KS %in% h$k), KSQ[[ds]] %in% h$k)
    g <- which(h$k %in% KS)
    i <- switch(mode,
                oracle = g[which.max(h$F2[g])],   # first maximum = smallest k among F2 ties
                fixed  = which(h$k == 10),
                sqrt   = which(h$k == KSQ[[ds]]))
    stopifnot(length(i) == 1)
    data.frame(dataset = ds, method = m, regime = rg,
               TPR = h$TPR[i], TNR = h$TNR[i], BA = h$BA[i], F2 = h$F2[i],
               sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
               k = h$k[i], stringsAsFactors = FALSE)
  }))
}
sqrt_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) pick_k(m, rg, "sqrt")))))
oracle_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) pick_k(m, rg, "oracle")))))
oracle_k$method <- paste0(oracle_k$method, "-oracle")
fixed_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) pick_k(m, rg, "fixed")))))
fixed_k$method <- paste0(fixed_k$method, "-k10")
plain <- do.call(rbind, lapply(c("ECOD", "COPOD", "GLOSH", "OPTICS"), function(m)
  do.call(rbind, lapply(c("T1", "T2"), function(rg) {
    g <- long[long$method == m & long$regime == rg, ]; g <- g[match(ORDER, g$dataset), ]
    data.frame(dataset = g$dataset, method = m, regime = rg,
               TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
               sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA, k = NA_real_,
               stringsAsFactors = FALSE)}))))
native_rows <- do.call(rbind, lapply(c("HDBSCAN-noise", "OPTICS-noise"), function(m) {
  g <- long[long$method == m, ]; g <- g[match(ORDER, g$dataset), ]
  data.frame(dataset = g$dataset, method = m, regime = "native",
             TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
             sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA, k = NA_real_,
             stringsAsFactors = FALSE)}))
# mutual-kNN m_i = 0: primary at the sqrt rule, label-chosen k kept as "-oracle"
m0_row <- function(ds, mode) {
  h <- long[long$method == "MutualKNN-m0" & long$dataset == ds, ]
  g <- which(h$k %in% KS)
  i <- if (mode == "oracle") g[which.max(h$F2[g])] else which(h$k == KSQ[[ds]])
  stopifnot(length(i) == 1)
  data.frame(dataset = ds, method = if (mode == "oracle") "MutualKNN-m0-oracle" else "MutualKNN-m0",
             regime = "native",
             TPR = h$TPR[i], TNR = h$TNR[i], BA = h$BA[i], F2 = h$F2[i],
             sd_TPR = NA, sd_TNR = NA, sd_BA = NA, sd_F2 = NA,
             k = h$k[i], stringsAsFactors = FALSE)
}
native_rows <- rbind(native_rows,
                     do.call(rbind, lapply(ORDER, m0_row, mode = "sqrt")),
                     do.call(rbind, lapply(ORDER, m0_row, mode = "oracle")))
main <- rbind(plain, seed_summary, sqrt_k, oracle_k, fixed_k, native_rows)

# gate: the rebuild must agree with 82's own summary CSV to numerical noise
chk <- read.csv(file.path(WP4, "wp4_metrics_main.csv"), stringsAsFactors = FALSE)
key <- function(d) paste(d$dataset, d$method, d$regime)
stopifnot(setequal(key(main), key(chk)))
j <- match(key(main), key(chk))
for (cl in mcols) stopifnot(max(abs(main[[cl]] - chk[[cl]][j])) < 1e-9)
cat("%% gate: rebuilt metrics agree with wp4_metrics_main.csv to < 1e-9\n")

# --- the nine incumbent methods, rebuilt as in 82_wp4_metrics.R -------------
fc <- read.csv(FC, stringsAsFactors = FALSE)
nw <- read.csv(NEW, stringsAsFactors = FALSE); nw <- nw[nw$status == "ok", ]
stopifnot(sum(fc$dataset == "shuffle") == 9)
fc$dataset[fc$dataset == "shuffle"] <- "Shuttle"
kk <- which(tolower(fc$dataset) == "vertebral" & fc$method == "SU-MCCD")
stopifnot(length(kk) == 1)
fc$BA[kk] <- 0.493; fc$F2[kk] <- 0.129
for (i in seq_len(nrow(nw))) {
  j <- which(tolower(fc$dataset) == tolower(nw$dataset[i]) & fc$method == nw$method[i])
  if (length(j) == 1) {
    fc$F2[j] <- round(nw$F2[i], 3); fc$BA[j] <- round(nw$BA[i], 3)
    fc$TPR[j] <- round(nw$TPR[i], 3); fc$TNR[j] <- round(nw$TNR[i], 3)
  }
}
inc <- function(ds, m, col) {
  v <- fc[[col]][tolower(fc$dataset) == tolower(ds) & fc$method == m]
  stopifnot(length(v) == 1); r3(v)
}

# --- the eight new methods, T1 ---------------------------------------------
wp4v <- function(ds, m, col, regime = "T1") {
  v <- main[[col]][main$dataset == ds & main$method == m & main$regime == regime]
  stopifnot(length(v) == 1); r3(v)
}

# combined per-cell F2/BA table over all 17 methods, T1
cell <- do.call(rbind, lapply(ORDER, function(ds) {
  do.call(rbind, lapply(c(INCUMBENT, NEWM), function(m) data.frame(
    dataset = ds, method = m,
    F2 = if (m %in% INCUMBENT) inc(ds, m, "F2") else wp4v(ds, m, "F2"),
    BA = if (m %in% INCUMBENT) inc(ds, m, "BA") else wp4v(ds, m, "BA"),
    stringsAsFactors = FALSE)))
}))

hdr <- function(s) cat("\n%% ---------- ", s, " ----------\n", sep = "")

# ===========================================================================
# 1. main text: tab:Real_Data_Aggregate, two blocks
# ===========================================================================
aggr <- function(ms) {
  a <- do.call(rbind, lapply(ms, function(m) {
    g <- cell[cell$method == m, ]
    data.frame(method = m, mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
               mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
               stringsAsFactors = FALSE)
  }))
  a[order(-a$mean_F2), ]
}
BLOCK2 <- setdiff(c(INCUMBENT, NEWM), BLOCK1)
a1 <- aggr(BLOCK1); a2 <- aggr(BLOCK2); a17 <- aggr(c(INCUMBENT, NEWM))

hdr("main text tab:Real_Data_Aggregate -- BLOCK 1 (no contamination input)")
for (i in seq_len(nrow(a1))) {
  b <- a1$method[i] == "SUN-MCCD"
  w <- function(x) if (b) sprintf("\\textbf{%s}", f3(x)) else f3(x)
  nm <- if (b) "\\textbf{SUN-MCCD}" else a1$method[i]
  cat(sprintf("  %-18s & %s & %s & %s & %s \\\\ \\hline\n", nm,
              w(a1$mean_F2[i]), w(a1$median_F2[i]), w(a1$mean_BA[i]), w(a1$median_BA[i])))
}
hdr("main text tab:Real_Data_Aggregate -- BLOCK 2 (given a contamination fraction)")
for (i in seq_len(nrow(a2)))
  cat(sprintf("  %-18s & %s & %s & %s & %s \\\\ \\hline\n", a2$method[i],
              f3(a2$mean_F2[i]), f3(a2$median_F2[i]), f3(a2$mean_BA[i]), f3(a2$median_BA[i])))

cat("\n%% block-1 check: is SUN-MCCD first on all four?\n")
for (col in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  v <- a1[[col]]; names(v) <- a1$method
  cat(sprintf("%%   %-10s SUN-MCCD %.3f, best other %.3f (%s) -> %s\n", col,
              v[["SUN-MCCD"]], max(v[names(v) != "SUN-MCCD"]),
              paste(names(v)[names(v) != "SUN-MCCD" &
                             v == max(v[names(v) != "SUN-MCCD"])], collapse = "/"),
              ifelse(v[["SUN-MCCD"]] == max(v), "FIRST", "NOT FIRST")))
}
cat("%% overall rank of SUN-MCCD among the 17:\n")
for (col in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  v <- a17[[col]]; names(v) <- a17$method
  bo <- max(v[names(v) != "SUN-MCCD"])
  cat(sprintf("%%   %-10s rank %2d of 17 | %.3f vs %.3f (%s) | %+0.3f\n", col,
              sum(v > v[["SUN-MCCD"]]) + 1, v[["SUN-MCCD"]], bo,
              paste(names(v)[names(v) != "SUN-MCCD" & v == bo], collapse = "/"),
              v[["SUN-MCCD"]] - bo))
}

# ===========================================================================
# 2. main text: tab:Real_Data_Result_Summary, rebuilt over 17 methods
# ===========================================================================
hdr("main text tab:Real_Data_Result_Summary (17 methods, T1)")
for (ds in ORDER) {
  g <- cell[cell$dataset == ds, ]
  f <- g$F2; names(f) <- g$method
  bo <- names(f)[f == max(f)]
  p <- f[names(f) %in% PROPOSED]; bp <- names(p)[p == max(p)]
  cat(sprintf("  %-12s & %-12s & %s & %-8s & %s & %d \\\\ \\hline\n", ds,
              paste(bo, collapse = "/"), f3(max(f)),
              paste(bp, collapse = "/"), f3(max(p)),
              sum(f > f[["SUN-MCCD"]]) + 1))
}
prop_win <- sum(sapply(ORDER, function(ds) {
  g <- cell[cell$dataset == ds, ]; any(g$method[g$F2 == max(g$F2)] %in% PROPOSED) }))
sun_sole <- sum(sapply(ORDER, function(ds) {
  g <- cell[cell$dataset == ds, ]; w <- g$method[g$F2 == max(g$F2)]
  length(w) == 1 && w == "SUN-MCCD" }))
bp_sun <- sum(sapply(ORDER, function(ds) {
  g <- cell[cell$dataset == ds & cell$method %in% PROPOSED, ]
  "SUN-MCCD" %in% g$method[g$F2 == max(g$F2)] }))
cat(sprintf("%%   a proposed method attains the highest F2 on %d of 16\n", prop_win))
cat(sprintf("%%   SUN-MCCD is the sole best overall on %d of 16\n", sun_sole))
cat(sprintf("%%   SUN-MCCD is the best proposed method on %d of 16\n", bp_sun))
lowd <- ORDER[c(3,4,5,6,7,9,10,14,15,16)]   # d <= 10
hid  <- setdiff(ORDER, lowd)                 # d >= 12
for (nmn in names(list(lowd = lowd, hid = hid))) {
  s <- list(lowd = lowd, hid = hid)[[nmn]]
  n <- sum(sapply(s, function(ds) { g <- cell[cell$dataset == ds, ]
    any(g$method[g$F2 == max(g$F2)] %in% PROPOSED) }))
  cat(sprintf("%%   proposed best on %d of %d data sets (%s)\n", n, length(s), nmn))
}

# ===========================================================================
# 3. main-text prose numbers: the best competitor on named data sets
# ===========================================================================
hdr("prose: best non-proposed method per data set (T1)")
for (ds in ORDER) {
  g <- cell[cell$dataset == ds & !(cell$method %in% PROPOSED), ]
  b <- g[g$F2 == max(g$F2), ]
  cat(sprintf("%%   %-12s best competitor %-14s F2 %s | SUN-MCCD %s | best proposed %s\n",
              ds, paste(b$method, collapse = "/"), f3(max(g$F2)),
              f3(cell$F2[cell$dataset == ds & cell$method == "SUN-MCCD"]),
              f3(max(cell$F2[cell$dataset == ds & cell$method %in% PROPOSED]))))
}

# ===========================================================================
# 4. supplement: undivided 17-method aggregate, T1
# ===========================================================================
hdr("supplement: undivided 17-method aggregate, T1, ordered by mean F2")
for (i in seq_len(nrow(a17)))
  cat(sprintf("  %-10s & %s & %s & %s & %s \\\\ \\hline\n", a17$method[i],
              f3(a17$mean_F2[i]), f3(a17$median_F2[i]), f3(a17$mean_BA[i]), f3(a17$median_BA[i])))

# ===========================================================================
# 5. supplement: T2 oracle sensitivity aggregate
# ===========================================================================
cell2 <- do.call(rbind, lapply(ORDER, function(ds) {
  do.call(rbind, lapply(c(INCUMBENT, NEWM), function(m) data.frame(
    dataset = ds, method = m,
    F2 = if (m %in% INCUMBENT) inc(ds, m, "F2") else wp4v(ds, m, "F2", "T2"),
    BA = if (m %in% INCUMBENT) inc(ds, m, "BA") else wp4v(ds, m, "BA", "T2"),
    stringsAsFactors = FALSE)))
}))
a17b <- do.call(rbind, lapply(c(INCUMBENT, NEWM), function(m) {
  g <- cell2[cell2$method == m, ]
  data.frame(method = m, mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
             mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
             stringsAsFactors = FALSE) }))
a17b <- a17b[order(-a17b$mean_F2), ]
hdr("supplement: T2 oracle-regime aggregate (competitors re-thresholded, incumbents NOT)")
for (i in seq_len(nrow(a17b))) {
  star <- if (a17b$method[i] %in% NEWM) "$^{*}$" else ""
  cat(sprintf("  %-10s%-6s & %s & %s & %s & %s \\\\ \\hline\n", a17b$method[i], star,
              f3(a17b$mean_F2[i]), f3(a17b$median_F2[i]), f3(a17b$mean_BA[i]), f3(a17b$median_BA[i])))
}

# ===========================================================================
# 6. supplement: per-cell TPR/TNR and BA/F2 for the eight new methods
# ===========================================================================
hdr("supplement: TPR and TNR, eight new methods, T1")
for (ds in ORDER) {
  v <- unlist(lapply(NEWM, function(m) c(f3(wp4v(ds, m, "TPR")), f3(wp4v(ds, m, "TNR")))))
  cat(sprintf("  %-12s & %s \\\\ \\hline\n", ds, paste(v, collapse = " & ")))
}
hdr("supplement: BA and F2, eight new methods, T1")
for (ds in ORDER) {
  v <- unlist(lapply(NEWM, function(m) c(f3(wp4v(ds, m, "BA")), f3(wp4v(ds, m, "F2")))))
  cat(sprintf("  %-12s & %s \\\\ \\hline\n", ds, paste(v, collapse = " & ")))
}

# ===========================================================================
# 7. supplement: mutual-kNN comparison
# ===========================================================================
hdr("supplement: oracle-k mutual-kNN vs SUN-MCCD and UN-MCCD, T1")
mk   <- main[main$method == "MutualKNN"     & main$regime == "T1", ]
mk10 <- main[main$method == "MutualKNN-k10" & main$regime == "T1", ]
mk0  <- main[main$method == "MutualKNN-m0", ]
mkt <- do.call(rbind, lapply(ORDER, function(ds) data.frame(
  dataset = ds, kstar = mk$k[mk$dataset == ds],
  mF2 = r3(mk$F2[mk$dataset == ds]),   mBA = r3(mk$BA[mk$dataset == ds]),
  k10 = r3(mk10$F2[mk10$dataset == ds]), m0 = r3(mk0$F2[mk0$dataset == ds]),
  sF2 = inc(ds, "SUN-MCCD", "F2"), sBA = inc(ds, "SUN-MCCD", "BA"),
  uF2 = inc(ds, "UN-MCCD",  "F2"), uBA = inc(ds, "UN-MCCD",  "BA"),
  stringsAsFactors = FALSE)))
for (i in seq_len(nrow(mkt)))
  cat(sprintf("  %-12s & %2d & %s & %s & %s & %s & %s & %s & %s & %s \\\\ \\hline\n",
      mkt$dataset[i], mkt$kstar[i], f3(mkt$mF2[i]), f3(mkt$mBA[i]), f3(mkt$k10[i]),
      f3(mkt$m0[i]), f3(mkt$sF2[i]), f3(mkt$sBA[i]), f3(mkt$uF2[i]), f3(mkt$uBA[i])))
cat(sprintf("%%   SUN-MCCD F2 ahead on %d of 16, BA ahead on %d of 16\n",
            sum(mkt$sF2 > mkt$mF2), sum(mkt$sBA > mkt$mBA)))
cat(sprintf("%%   UN-MCCD F2 ahead on %d of 16; mutual-kNN ahead on %d; tied %d\n",
            sum(mkt$uF2 > mkt$mF2), sum(mkt$mF2 > mkt$uF2), sum(mkt$uF2 == mkt$mF2)))
cat(sprintf("%%   mean F2: oracle-k %.3f, k=10 %.3f, m_i=0 %.3f, SUN %.3f, UN %.3f\n",
            r3(mean(mkt$mF2)), r3(mean(mkt$k10)), r3(mean(mkt$m0)),
            r3(mean(mkt$sF2)), r3(mean(mkt$uF2))))
cat(sprintf("%%   mean BA: oracle-k %.3f, SUN %.3f, UN %.3f\n",
            r3(mean(mkt$mBA)), r3(mean(mkt$sBA)), r3(mean(mkt$uBA))))
sn   <- main[main$method == "SNN"     & main$regime == "T1", ]
sn10 <- main[main$method == "SNN-k10" & main$regime == "T1", ]
cat(sprintf("%%   SNN mean F2: oracle-k %.3f, fixed k=10 %.3f\n",
            r3(mean(r3(sn$F2[match(ORDER, sn$dataset)]))),
            r3(mean(r3(sn10$F2[match(ORDER, sn10$dataset)])))))   # 82 rounds per cell first

# ===========================================================================
# 8. supplement: threshold-free native variants
# ===========================================================================
hdr("supplement: threshold-free native variants (no contamination argument)")
nat <- do.call(rbind, lapply(c("HDBSCAN-noise", "OPTICS-noise", "MutualKNN-m0"),
  function(m) { g <- main[main$method == m & main$regime == "native", ]
    data.frame(method = m, mean_F2 = r3(mean(g$F2)), median_F2 = r3(median(g$F2)),
               mean_BA = r3(mean(g$BA)), median_BA = r3(median(g$BA)),
               stringsAsFactors = FALSE) }))
nat <- nat[order(-nat$mean_F2), ]
for (i in seq_len(nrow(nat)))
  cat(sprintf("  %-16s & %s & %s & %s & %s \\\\ \\hline\n", nat$method[i],
              f3(nat$mean_F2[i]), f3(nat$median_F2[i]), f3(nat$mean_BA[i]), f3(nat$median_BA[i])))

# ===========================================================================
# 9. supplement: DIF and LUNAR seed stability
# ===========================================================================
hdr("supplement: DIF and LUNAR, mean and s.d. over five seeds, T1")
sd4 <- function(x) sprintf("%.4f", x)
for (ds in ORDER) {
  d <- main[main$dataset == ds & main$method == "DIF"   & main$regime == "T1", ]
  l <- main[main$dataset == ds & main$method == "LUNAR" & main$regime == "T1", ]
  cat(sprintf("  %-12s & %s & %s & %s & %s & %s & %s & %s & %s \\\\ \\hline\n", ds,
      f3(r3(d$F2)), sd4(d$sd_F2), f3(r3(d$BA)), sd4(d$sd_BA),
      f3(r3(l$F2)), sd4(l$sd_F2), f3(r3(l$BA)), sd4(l$sd_BA)))
}
for (m in c("DIF", "LUNAR")) {
  h <- main[main$method == m & main$regime == "T1", ]
  cat(sprintf("%%   %-5s sd(F2) %.4f--%.4f mean %.4f (max on %s); sd(BA) %.4f--%.4f mean %.4f\n",
      m, min(h$sd_F2), max(h$sd_F2), mean(h$sd_F2), h$dataset[which.max(h$sd_F2)],
      min(h$sd_BA), max(h$sd_BA), mean(h$sd_BA)))
}

# ===========================================================================
# 10. R1.3: HDBSCAN/GLOSH against the shape-adaptive pair
# ===========================================================================
hdr("R1.3 counts: GLOSH vs SU-MCCD and SUN-MCCD, T1")
gl <- sapply(ORDER, function(ds) wp4v(ds, "GLOSH", "F2"))
su <- sapply(ORDER, function(ds) inc(ds, "SU-MCCD",  "F2"))
sn2 <- sapply(ORDER, function(ds) inc(ds, "SUN-MCCD", "F2"))
cat(sprintf("%%   SU-MCCD F2 ahead of GLOSH on %d of 16; GLOSH ahead on %d; tied %d\n",
            sum(su > gl), sum(gl > su), sum(su == gl)))
cat(sprintf("%%   SUN-MCCD F2 ahead of GLOSH on %d of 16; GLOSH ahead on %d; tied %d\n",
            sum(sn2 > gl), sum(gl > sn2), sum(sn2 == gl)))
cat(sprintf("%%   mean F2: GLOSH %.3f, OPTICS %.3f, SU-MCCD %.3f, SUN-MCCD %.3f\n",
            r3(mean(gl)), r3(mean(sapply(ORDER, function(ds) wp4v(ds, "OPTICS", "F2")))),
            r3(mean(su)), r3(mean(sn2))))

# ===========================================================================
# 2026-09-26: label-free k for mutual-kNN and SNN (round(sqrt(n)), ODIN's rule)
# is the primary result above; the label-chosen k is printed here as an upper
# bound, unranked, for the main-text aggregate table and the supplement k table.
# ===========================================================================
hdr("main text tab:Real_Data_Aggregate -- upper-bound rows (k chosen with the true labels; not ranked)")
for (m in c("MutualKNN-oracle", "SNN-oracle")) {
  g <- main[main$method == m & main$regime == "T1", ]; g <- g[match(ORDER, g$dataset), ]
  f2 <- sapply(ORDER, function(ds) wp4v(ds, m, "F2")); ba <- sapply(ORDER, function(ds) wp4v(ds, m, "BA"))
  cat(sprintf("  %-18s & %s & %s & %s & %s \\\\ \\hline\n", m,
              f3(mean(f2)), f3(median(f2)), f3(mean(ba)), f3(median(ba))))
}
hdr("supplement: mutual-kNN at k = round(sqrt(n)) and at the label-chosen k, vs SUN-MCCD and UN-MCCD, T1")
for (ds in ORDER) {
  kq <- main$k[main$dataset == ds & main$method == "MutualKNN" & main$regime == "T1"]
  ko <- main$k[main$dataset == ds & main$method == "MutualKNN-oracle" & main$regime == "T1"]
  cat(sprintf("  %-12s & %2d & %s & %s & %2d & %s & %s & %s & %s & %s & %s \\\\ \\hline\n", ds,
              kq, f3(wp4v(ds, "MutualKNN", "F2")), f3(wp4v(ds, "MutualKNN", "BA")),
              ko, f3(wp4v(ds, "MutualKNN-oracle", "F2")), f3(wp4v(ds, "MutualKNN-oracle", "BA")),
              f3(inc(ds, "SUN-MCCD", "F2")), f3(inc(ds, "SUN-MCCD", "BA")),
              f3(inc(ds, "UN-MCCD", "F2")), f3(inc(ds, "UN-MCCD", "BA"))))
}
cnt <- function(a, b) sprintf("%d ahead, %d behind, %d tied", sum(a > b), sum(b > a), sum(a == b))
for (m in c("MutualKNN", "MutualKNN-oracle")) for (p in c("SUN-MCCD", "UN-MCCD")) for (col in c("F2", "BA")) {
  a <- sapply(ORDER, function(ds) inc(ds, p, col)); b <- sapply(ORDER, function(ds) wp4v(ds, m, col))
  cat(sprintf("%%   %-8s vs %-16s on %s: %s\n", p, m, col, cnt(a, b)))
}
for (m in c("MutualKNN", "MutualKNN-oracle", "MutualKNN-k10", "SNN", "SNN-oracle", "SNN-k10")) {
  f2 <- sapply(ORDER, function(ds) wp4v(ds, m, "F2")); ba <- sapply(ORDER, function(ds) wp4v(ds, m, "BA"))
  cat(sprintf("%%   %-17s mean F2 %s  mean BA %s\n", m, f3(mean(f2)), f3(mean(ba))))
}
hdr("supplement: mutual-kNN m_i = 0 (native), sqrt rule and label-chosen k")
for (m in c("MutualKNN-m0", "MutualKNN-m0-oracle")) {
  g <- main[main$method == m & main$regime == "native", ]
  cat(sprintf("  %-20s & %s & %s & %s & %s \\\\ \\hline\n", m, f3(mean(r3(g$F2))), f3(median(r3(g$F2))),
              f3(mean(r3(g$BA))), f3(median(r3(g$BA)))))
}
