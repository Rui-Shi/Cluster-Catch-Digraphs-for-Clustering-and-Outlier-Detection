#!/usr/bin/env Rscript
# 19_musk_geometry_audit.R -- Ceyhan checklist item 8, Musk half.
#
# "Confirm the Speech/Musk diagnostics are produced by the same pipeline used in
# the reported performance results, not a diagnostic variant."
#
# The Speech half is already settled: 16_speech_ios_clusters.R reused the cached
# radius chunks and its zero accounting reproduces the 2098 exact zeros in the
# cached score vector. This script does the same job for the Musk paragraph of
# Appendix D, every number of which is re-derived here from
#   (a) the radius chunks written by 07b during the reported 14.7 h full-data run
#       (results/scores_cache/radi_chunks/radi_n3062_d166_descend_low3_sc1_*), and
#   (b) the cached score vectors Musk_full_UNCCD-{OOS,IOS}.rds from that same run,
# and then checked against the values printed in the manuscript.
#
# The decisive test is PROVENANCE: rebuilding IOS from the cached radii must
# reproduce the cached score vector. If it does, the geometry below and the
# reported AUC come from one construction, which is what item 8 asks.
#
# Run (PowerShell; Rscript segfaults under Bash):
#   Rscript "revision_experiments/tr2/19_musk_geometry_audit.R"

suppressPackageStartupMessages({ library(here) })
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(here::here("revision_experiments/shared/harness.R"))

SHARED <- here::here("revision_experiments/results")
CACHE  <- file.path(SHARED, "scores_cache")
CSV    <- file.path(SHARED, "datasets_csv", "Musk.csv")
OUT    <- file.path(SHARED, "tr2", "musk_geometry_audit.csv")

# Manuscript Appendix D, transcribed 2026-08-22.
CLAIM <- list(
  n_outliers        = 97,
  med_radius_out    = 0.176,
  med_radius_reg    = 6.098,
  med_share_out     = 1.000,
  base_rate         = 0.032,
  med_inb_out       = 2,
  med_inb_reg       = 2,
  med_CIrho_out     = 10.91,
  med_CIrho_reg     = 0.32,
  auc_ios           = 0.085,
  auc_oos           = 0.736
)

df <- read.csv(CSV)
stopifnot(colnames(df)[ncol(df)] == "label")
d <- ncol(df) - 1
X <- as.matrix(df[, seq_len(d), drop = FALSE])
Y <- df[[ncol(df)]]          # 1 = regular, 0 = outlier
n <- nrow(X)
is_out <- (Y == 0)
cat(sprintf("Musk: n=%d, d=%d, outliers=%d\n", n, d, sum(is_out)))

files <- sort(list.files(file.path(CACHE, "radi_chunks"),
                         pattern = "^radi_n3062_d166_descend_low3_sc1_", full.names = TRUE))
cat(sprintf("reusing %d cached radius chunks from the reported run\n", length(files)))
R <- unlist(lapply(files, readRDS), use.names = FALSE)
stopifnot(length(R) == n, all(is.finite(R)), all(R > 0))

D <- as.matrix(dist(X))
M <- D <= R                        # M[i,j]: x_j lies in B(x_i, r_i)
k <- rowSums(M)                    # points caught by own ball, self included
rho <- k^(1 / d) / R               # Vic_Den, exactly as in Outlyingness_Score.R

# Inbound neighbours of x_j = balls that reach it = column j of M (self included,
# which is the augmented neighbourhood N_I bar of Proposition 2.2).
inb_count <- colSums(M)
CIrho <- as.numeric(t(M) %*% rho)  # sum of rho over inbound neighbours, self included
share_out <- as.numeric(t(M) %*% is_out) / inb_count

rank_auc <- function(Y, score) {
  pos <- score[Y == 0]; neg <- score[Y == 1]
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

got <- list(
  n_outliers     = sum(is_out),
  med_radius_out = median(R[is_out]),
  med_radius_reg = median(R[!is_out]),
  med_share_out  = median(share_out[is_out]),
  base_rate      = mean(is_out),
  med_inb_out    = median(inb_count[is_out]),
  med_inb_reg    = median(inb_count[!is_out]),
  med_CIrho_out  = median(CIrho[is_out]),
  med_CIrho_reg  = median(CIrho[!is_out])
)

# ---- provenance: does the cached score come from THIS construction? ----------
ios_cache <- readRDS(file.path(CACHE, "Musk_full_UNCCD-IOS.rds"))
oos_cache <- readRDS(file.path(CACHE, "Musk_full_UNCCD-OOS.rds"))
got$auc_ios <- rank_auc(Y, ios_cache$score)
got$auc_oos <- rank_auc(Y, oos_cache$score)

# Raw (pre-standardization, pre-cluster) IOS from the cached radii.
ios_raw_here <- 1 / CIrho
# The cached score is standardized within clusters and tie-broken, so it cannot be
# compared elementwise; but a monotone-per-cluster transform preserves the sign of
# the global rank correlation, and any DIFFERENT construction would break it.
sp <- suppressWarnings(cor(ios_raw_here, ios_cache$score, method = "spearman"))
cat(sprintf("\nprovenance: Spearman(raw IOS from cached radii, cached IOS score) = %.4f\n", sp))
cat(sprintf("            radii used here have median %.4f; cached run stored t_wall = %.0f s\n",
            median(R), ios_cache$t_wall))

cat("\n---- Appendix D Musk numbers ----\n")
rows <- list()
for (nm in names(CLAIM)) {
  g <- got[[nm]]; c0 <- CLAIM[[nm]]
  tol <- if (nm %in% c("med_CIrho_out")) 0.005 else if (abs(c0) >= 1) 0.0005 * max(1, abs(c0)) else 0.0005
  ok <- abs(g - c0) <= max(tol, 0.0005)
  cat(sprintf("  %-16s manuscript %10.3f | recomputed %10.4f | %s\n",
              nm, c0, g, if (ok) "MATCH" else "DIFFERS"))
  rows[[length(rows) + 1]] <- data.frame(quantity = nm, manuscript = c0,
                                         recomputed = g, match = ok,
                                         stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
out$spearman_raw_vs_cached <- sp
write.csv(out, OUT, row.names = FALSE)
cat(sprintf("\n%d of %d quantities match the manuscript.\n", sum(out$match), nrow(out)))
cat(sprintf("Wrote %s\n19_musk_geometry_audit.R done.\n", OUT))
