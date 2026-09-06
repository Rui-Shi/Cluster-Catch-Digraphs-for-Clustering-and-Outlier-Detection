#!/usr/bin/env Rscript
# revision_experiments/tr1/84b_wp5_metrics.R
#
# WP5 (AE.4, R1.8, R3.8, R5.6): turns the WP4-competitor raw scores (from
# 81_wp4_baselines.py, called with --data-dir/--out-dir pointed at the WP5
# folders -- see WP5_PROTOCOL.md S7) into TPR/TNR/BA/F2, and merges them with
# 84_wp5_highd.R's own MCCD + baseline output into one per-(dataset, method)
# table for the four real data sets above d = 21.
#
# THRESHOLDING LOGIC -- reused verbatim from tr1/82_wp4_metrics.R:138-149,
# not reimplemented independently:
#
#   thr <- quantile(score, 1 - contamination, type = 7, names = FALSE)
#   is_out <- as.numeric(score > thr)              # STRICT >, pyod's rule
#   v <- evaluate(labs[[ds]], is_out, 0.5)
#
# i.e. pyod sets threshold_ = np.percentile(decision_scores_, 100*(1-cont)),
# which is R's quantile(..., type = 7); labels_ = (decision_scores_ >
# threshold_), strictly greater; and the resulting binary indicator is handed
# to evaluate() at threshold 0.5 so the confusion counts are pyod's own while
# the TPR/TNR/BA/F2 arithmetic is this study's evaluate(). T1 (contamination
# = 0.1) is WP4's main-table setting and the only regime WP5 computes (no T2
# oracle sensitivity here -- WP5_PROTOCOL.md S4.3).
#
# The native threshold-free variants (HDBSCAN noise, OPTICS noise, mutual-kNN
# m_i = 0) and the oracle-k selection on F2 for mutual-kNN/SNN follow
# 82_wp4_metrics.R's own logic (its lines 190-261) at the same settings.
#
# INPUTS
#   results/tr1/wp5/data/<dataset>.csv, manifest.csv    (84a's output)
#   results/tr1/wp5/scores/<dataset>_<tag>.csv           (81's output, WP5 folder)
#   results/tr1/wp5/wp5_highd_results.csv                (84's output: 4 MCCD + 5 baselines)
#
# OUTPUTS
#   results/tr1/wp5/wp5_metrics_main.csv   one row per (dataset, method), all 17 methods
#   results/tr1/wp5/WP5_FINDINGS.md        narrative summary
#
# USAGE
#   Rscript "revision_experiments/tr1/84b_wp5_metrics.R"
#
# This script computes no thresholds of its own invention and reinterprets
# nothing from WP4_PROTOCOL.md or WP5_PROTOCOL.md -- it is downstream
# arithmetic only.

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments/shared/harness.R")))
options(width = 200)

WP5  <- here::here("revision_experiments/results/tr1/wp5")
DATA <- file.path(WP5, "data")
SCO  <- file.path(WP5, "scores")
HIGHD_CSV <- file.path(WP5, "wp5_highd_results.csv")

OUT_MAIN <- file.path(WP5, "wp5_metrics_main.csv")
OUT_MD   <- file.path(WP5, "WP5_FINDINGS.md")

ORDER <- c("letter", "mnist", "musk", "arrhythmia")
KS <- c(5, 10, 15, 20, 30)
SEEDS <- 1:5
r3 <- function(x) as.numeric(sprintf("%.3f", x))

CON <- textConnection("REPORT_LINES", "w", local = FALSE)
say <- function(...) { txt <- paste0(...); cat(txt); writeLines(sub("\n$", "", txt), CON) }
sayf <- function(fmt, ...) say(sprintf(fmt, ...))
saydf <- function(df) {
  o <- capture.output(print(df, row.names = FALSE))
  cat(paste0(o, collapse = "\n"), "\n", sep = ""); writeLines(o, CON)
}

# ---------------------------------------------------------------------------
# 0. inputs
# ---------------------------------------------------------------------------
man <- read.csv(file.path(DATA, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(man$dataset, ORDER))

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
    stop(sprintf("length mismatch %s_%s: %d scores, %d labels", ds, tag, length(s), length(labs[[ds]])))
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

# pyod thresholding at T1 (contamination = 0.1), then evaluate() on the
# binary indicator -- 82_wp4_metrics.R:138-149, verbatim.
metrics_from_score <- function(ds, score, contamination) {
  thr <- quantile(score, 1 - contamination, type = 7, names = FALSE)
  is_out <- as.numeric(score > thr)
  v <- evaluate(labs[[ds]], is_out, 0.5)
  c(v, n_flagged = sum(is_out))
}
metrics_from_labels <- function(ds, is_out) {
  v <- evaluate(labs[[ds]], as.numeric(is_out), 0.5)
  c(v, n_flagged = sum(is_out))
}

# ---------------------------------------------------------------------------
# 1. WP4-competitor cells, T1 only (WP5_PROTOCOL.md S4.3: no T2 here)
# ---------------------------------------------------------------------------
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

emit <- function(ds, method, tag, k = NA, seed = NA) {
  s <- read_score(ds, tag)
  v <- metrics_from_score(ds, s, 0.10)
  add(dataset = ds, method = method, k = k, seed = seed,
      TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]],
      n_flagged = v[["n_flagged"]])
}

for (ds in ORDER) {
  emit(ds, "ECOD", "ECOD")
  emit(ds, "COPOD", "COPOD")
  for (s in SEEDS) emit(ds, "DIF", sprintf("DIF_seed%d", s), seed = s)
  for (s in SEEDS) emit(ds, "LUNAR", sprintf("LUNAR_seed%d", s), seed = s)
  emit(ds, "GLOSH", "HDBSCAN")
  emit(ds, "OPTICS", "OPTICS")
  for (k in KS) emit(ds, "MutualKNN", sprintf("MutualKNN_k%d", k), k = k)
  for (k in KS) emit(ds, "SNN", sprintf("SNN_k%d", k), k = k)

  nat <- list(c("HDBSCAN-noise", "HDBSCAN"), c("OPTICS-noise", "OPTICS"))
  for (p in nat) {
    v <- metrics_from_labels(ds, read_native(ds, p[2]))
    add(dataset = ds, method = p[1], k = NA, seed = NA,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]], n_flagged = v[["n_flagged"]])
  }
  for (k in KS) {
    v <- metrics_from_labels(ds, read_native(ds, sprintf("MutualKNN_k%d", k)))
    add(dataset = ds, method = "MutualKNN-m0", k = k, seed = NA,
        TPR = v[["TPR"]], TNR = v[["TNR"]], BA = v[["BA"]], F2 = v[["F2"]], n_flagged = v[["n_flagged"]])
  }
}
long <- do.call(rbind, rows)

# Recommended (verifier pass, 2026-09-05): keep the per-seed/per-k long frame
# on disk, not just in memory -- it is what the seed-stability and oracle-k
# read-offs below are computed from, and it should be inspectable independent
# of this script's console/markdown summary.
OUT_LONG <- file.path(WP5, "wp5_metrics_long.csv")
write.csv(long, OUT_LONG, row.names = FALSE)

# seeded methods: mean over 5 seeds (82_wp4_metrics.R's own convention)
seed_summary <- do.call(rbind, lapply(c("DIF", "LUNAR"), function(m) {
  g <- long[long$method == m, ]
  do.call(rbind, lapply(ORDER, function(ds) {
    h <- g[g$dataset == ds, ]; stopifnot(nrow(h) == 5)
    data.frame(dataset = ds, method = m,
               TPR = mean(h$TPR), TNR = mean(h$TNR), BA = mean(h$BA), F2 = mean(h$F2),
               sd_F2 = sd(h$F2), k = NA_real_, stringsAsFactors = FALSE)
  }))
}))

# k-swept methods: oracle-best k by F2 (82_wp4_metrics.R's declared oracle)
pick_k <- function(m) {
  do.call(rbind, lapply(ORDER, function(ds) {
    h <- long[long$method == m & long$dataset == ds, ]
    stopifnot(nrow(h) == length(KS))
    i <- which.max(h$F2)
    data.frame(dataset = ds, method = m, TPR = h$TPR[i], TNR = h$TNR[i], BA = h$BA[i], F2 = h$F2[i],
               sd_F2 = NA, k = h$k[i], stringsAsFactors = FALSE)
  }))
}
oracle_k <- do.call(rbind, lapply(c("MutualKNN", "SNN"), pick_k))

plain <- do.call(rbind, lapply(c("ECOD", "COPOD", "GLOSH", "OPTICS"), function(m) {
  g <- long[long$method == m, ]; g <- g[match(ORDER, g$dataset), ]
  data.frame(dataset = g$dataset, method = m, TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
             sd_F2 = NA, k = NA_real_, stringsAsFactors = FALSE)
}))

native_rows <- do.call(rbind, lapply(c("HDBSCAN-noise", "OPTICS-noise"), function(m) {
  g <- long[long$method == m, ]; g <- g[match(ORDER, g$dataset), ]
  data.frame(dataset = g$dataset, method = m, TPR = g$TPR, TNR = g$TNR, BA = g$BA, F2 = g$F2,
             sd_F2 = NA, k = NA_real_, stringsAsFactors = FALSE)
}))
native_rows <- rbind(native_rows, do.call(rbind, lapply(ORDER, function(ds) {
  h <- long[long$method == "MutualKNN-m0" & long$dataset == ds, ]
  i <- which.max(h$F2)
  data.frame(dataset = ds, method = "MutualKNN-m0", TPR = h$TPR[i], TNR = h$TNR[i],
             BA = h$BA[i], F2 = h$F2[i], sd_F2 = NA, k = h$k[i], stringsAsFactors = FALSE)
})))

wp4_competitors <- rbind(plain, seed_summary, oracle_k, native_rows)
wp4_competitors <- wp4_competitors[order(match(wp4_competitors$dataset, ORDER), wp4_competitors$method), ]

# ---------------------------------------------------------------------------
# 2. merge with 84's own MCCD + baseline output
# ---------------------------------------------------------------------------
if (!file.exists(HIGHD_CSV)) {
  stop("Missing ", HIGHD_CSV, " -- run 84_wp5_highd.R first (this script only scores WP4 competitors and merges).")
}
highd <- read.csv(HIGHD_CSV, stringsAsFactors = FALSE)

highd_slim <- data.frame(
  dataset = highd$dataset, method = highd$method,
  TPR = r3(highd$TPR), TNR = r3(highd$TNR), BA = r3(highd$BA), F2 = r3(highd$F2),
  status = highd$status, reason = highd$reason,
  stringsAsFactors = FALSE
)
wp4_slim <- data.frame(
  dataset = wp4_competitors$dataset, method = wp4_competitors$method,
  TPR = r3(wp4_competitors$TPR), TNR = r3(wp4_competitors$TNR),
  BA = r3(wp4_competitors$BA), F2 = r3(wp4_competitors$F2),
  status = "ok", reason = "", stringsAsFactors = FALSE
)

main <- rbind(highd_slim, wp4_slim)
main <- main[order(match(main$dataset, ORDER), main$method), ]

# Recommended (verifier pass, 2026-09-05): the merge above combines two
# independently-built tables (84's own MCCD+baseline output and the WP4
# competitor summary); assert the merge did not silently duplicate a
# (dataset, method) cell before writing it out as the study's table.
dup <- duplicated(main[, c("dataset", "method")])
if (any(dup)) {
  stop(sprintf("84b_wp5_metrics.R: %d duplicated (dataset, method) row(s) after merge: %s",
               sum(dup), paste(sprintf("%s/%s", main$dataset[dup], main$method[dup]), collapse = "; ")))
}

write.csv(main, OUT_MAIN, row.names = FALSE)

say("=== WP5: real data above d = 21 -- merged per-(dataset, method) table ===\n")
say("  4 proposed (U-MCCD/SU-MCCD n/a at musk/arrhythmia, per rk_degeneracy.csv) +\n")
say("  5 original baselines + 8 WP4 competitors (T1, contamination = 0.1)\n\n")
for (ds in ORDER) {
  say(sprintf("--- %s ---\n", ds))
  saydf(main[main$dataset == ds, c("method", "TPR", "TNR", "BA", "F2", "status", "reason")])
}

na_rows <- main[main$status == "n/a", ]
if (nrow(na_rows) > 0) {
  say("\n=== n/a rows (RK envelope degenerate at this dimension) ===\n")
  saydf(na_rows[, c("dataset", "method", "reason")])
}

# ---------------------------------------------------------------------------
# 3 (recommended, verifier pass 2026-09-05). R3.3 and R1.3 read-offs on the
# four WP5 data sets, same comparisons as 82_wp4_metrics.R's sections 6-7,
# extended above d = 21 (WP5_PROTOCOL.md S9's last bullet).
# ---------------------------------------------------------------------------
gv <- function(ds, m, col) {
  v <- main[[col]][main$dataset == ds & main$method == m]
  if (length(v) == 0) NA else v[1]
}

say("\n=== R3.3 read-off: oracle-k mutual-kNN vs SUN-MCCD / UN-MCCD (T1) ===\n")
mk_oracle <- oracle_k[oracle_k$method == "MutualKNN", ]
r33 <- do.call(rbind, lapply(ORDER, function(ds) data.frame(
  dataset = ds,
  k_star = mk_oracle$k[mk_oracle$dataset == ds],
  mkNN_F2 = r3(mk_oracle$F2[mk_oracle$dataset == ds]),
  mkNN_BA = r3(mk_oracle$BA[mk_oracle$dataset == ds]),
  SUN_F2 = gv(ds, "SUN-MCCD", "F2"), SUN_BA = gv(ds, "SUN-MCCD", "BA"),
  UN_F2 = gv(ds, "UN-MCCD", "F2"), UN_BA = gv(ds, "UN-MCCD", "BA"),
  stringsAsFactors = FALSE)))
saydf(r33)

say("\n=== R1.3 read-off: GLOSH (HDBSCAN) vs SU-MCCD / SUN-MCCD (T1) ===\n")
gl <- plain[plain$method == "GLOSH", ]
r13 <- do.call(rbind, lapply(ORDER, function(ds) data.frame(
  dataset = ds,
  GLOSH_F2 = r3(gl$F2[gl$dataset == ds]), GLOSH_BA = r3(gl$BA[gl$dataset == ds]),
  SU_F2 = gv(ds, "SU-MCCD", "F2"), SU_BA = gv(ds, "SU-MCCD", "BA"),
  SUN_F2 = gv(ds, "SUN-MCCD", "F2"), SUN_BA = gv(ds, "SUN-MCCD", "BA"),
  stringsAsFactors = FALSE)))
saydf(r13)

close(CON)
md <- c("# WP5 findings -- real data above d = 21",
        "",
        "Generated by `tr1/84b_wp5_metrics.R` from `tr1/84_wp5_highd.R`'s output",
        "(4 proposed MCCD detectors + 5 original baselines) and",
        "`tr1/81_wp4_baselines.py`'s output (8 WP4 competitors, run against the",
        "WP5 data folder via `--data-dir`/`--out-dir`, WP5_PROTOCOL.md S7).",
        "Thresholding for the WP4 competitors reuses `82_wp4_metrics.R`'s own T1",
        "logic verbatim (lines 138-149): `threshold = quantile(score, 0.9, type = 7)`,",
        "`is_outlier = score > threshold`, scored by the same `evaluate()`.",
        "", "```", REPORT_LINES, "```", "")
writeLines(md, OUT_MD)

cat(sprintf("\nwrote %s\nwrote %s\n", OUT_MAIN, OUT_MD))
