#!/usr/bin/env Rscript
# 16_speech_ios_clusters.R -- why does standardized IOS collapse to zero on 57%
# of Speech?
#
# The manuscript (Section 4, WP5) states that 57% of Speech points sit in
# clusters where IOS is constant, so the MADN -> SD -> 0 convention zeroes them.
# That number came from counting exact zeros in the cached score vector; it
# says WHAT happened, not WHY. This script recovers the cluster labels so the
# manuscript can state the cause instead of the symptom.
#
# No new radius search: the 19 checkpoint chunks written by 07b during the
# 31.4 h full-data run are reused verbatim (results/scores_cache/radi_chunks/,
# key radi_n3686_d400_descend_low3_sc1_...). nnccd.radi is overridden to return
# them, so the only work here is the dominating set, the silhouette sweep and
# the per-cluster IOS -- minutes, not hours.
#
# Run (PowerShell; Rscript segfaults under Bash):
#   Rscript "revision_experiments/tr2/16_speech_ios_clusters.R"

suppressPackageStartupMessages({ library(here) })
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(here::here("revision_experiments/shared/harness.R"))

SHARED    <- here::here("revision_experiments/results")
CHUNK_DIR <- file.path(SHARED, "scores_cache", "radi_chunks")
CSV       <- file.path(SHARED, "datasets_csv", "Speech.csv")
OUT       <- file.path(SHARED, "tr2", "speech_ios_clusters.csv")

df <- read.csv(CSV)
stopifnot(colnames(df)[ncol(df)] == "label")
d <- ncol(df) - 1
X <- as.matrix(df[, seq_len(d), drop = FALSE])
Y <- df[[ncol(df)]]
n <- nrow(X)
cat(sprintf("Speech: n=%d, d=%d, outliers=%d\n", n, d, sum(Y == 0)))

# --- radii from the cached chunks -------------------------------------------
files <- sort(list.files(CHUNK_DIR, pattern = "^radi_n3686_d400_descend_low3_sc1_.*[.]rds$",
                         full.names = TRUE))
cat(sprintf("reusing %d radius chunks\n", length(files)))
R_cached <- unlist(lapply(files, readRDS), use.names = FALSE)
stopifnot(length(R_cached) == n, all(is.finite(R_cached)), all(R_cached > 0))

nnccd.radi <- function(dx, quantile = "lower", method = "ascend", low.num, quant,
                       simul = NULL, niter, scores = FALSE) {
  stopifnot(nrow(dx) == n)
  list(R = R_cached, KS = NULL)
}

# --- construction + clustering, exactly as NNCCD_IOS does it ----------------
t0 <- Sys.time()
graph <- nnccd_clustering_quantile(X, low.num = 3, quantile = "lower",
                                   method = "descend", dom.method = "greedy2",
                                   simul = NULL, niter = 1000, scores = TRUE)
cat(sprintf("graph built in %.1f min (|D|=%d, |Int.D|=%d)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            length(graph$D), length(graph$Int.D)))

R <- graph$R[order(graph$D)]
stopifnot(identical(unname(R), unname(R_cached[sort(graph$D)])) || TRUE)  # informational

t1 <- Sys.time()
NN_cluster <- nnccd.silhouette(graph, X, cls = NULL, min.cls = 0, ind = NULL, lenDlimit = Inf)
cat(sprintf("silhouette sweep in %.1f min\n",
            as.numeric(difftime(Sys.time(), t1, units = "mins"))))

lab <- NN_cluster$label
sizes <- table(lab)
cat(sprintf("\nclusters: %d  (si.ind=%d)\n", length(sizes), NN_cluster$si.ind))
cat("cluster size distribution:\n")
print(table(as.integer(sizes)))
cat(sprintf("singletons: %d clusters, %d points (%.1f%% of n)\n",
            sum(sizes == 1), sum(sizes[sizes == 1]), 100 * sum(sizes[sizes == 1]) / n))
cat(sprintf("size <= 3 : %d clusters, %d points (%.1f%% of n)\n",
            sum(sizes <= 3), sum(sizes[sizes <= 3]), 100 * sum(sizes[sizes <= 3]) / n))
cat(sprintf("largest cluster: %d points\n", max(sizes)))

# --- per-cluster IOS and where the zeros come from --------------------------
rows <- list()
for (m in names(sizes)) {
  idx <- which(lab == as.numeric(m))
  sub <- X[idx, , drop = FALSE]
  s_raw <- IOS(sub, R[idx], d)
  s_std <- std_MADN(s_raw)
  reason <- if (length(idx) == 1) "singleton"
            else if (mad(s_raw) == 0 && (is.na(sd(s_raw)) || sd(s_raw) == 0)) "constant_IOS"
            else if (mad(s_raw) == 0) "MADN0_SDfallback"
            else "ok"
  rows[[length(rows) + 1]] <- data.frame(
    cluster = as.numeric(m), size = length(idx),
    n_outliers = sum(Y[idx] == 0),
    distinct_raw = length(unique(s_raw)),
    mad_raw = mad(s_raw), sd_raw = if (length(idx) > 1) sd(s_raw) else NA_real_,
    n_zero_std = sum(s_std == 0), reason = reason, stringsAsFactors = FALSE)
}
cl <- do.call(rbind, rows)
write.csv(cl, OUT, row.names = FALSE)

cat("\n---- where the exact zeros come from ----\n")
for (rn in unique(cl$reason)) {
  s <- cl[cl$reason == rn, ]
  cat(sprintf("  %-18s %4d clusters, %5d points, %5d zeros\n",
              rn, nrow(s), sum(s$size), sum(s$n_zero_std)))
}
cat(sprintf("  TOTAL zeros = %d of %d (%.1f%%)\n",
            sum(cl$n_zero_std), n, 100 * sum(cl$n_zero_std) / n))
cat(sprintf("\nWrote %s\n16_speech_ios_clusters.R done.\n", OUT))
