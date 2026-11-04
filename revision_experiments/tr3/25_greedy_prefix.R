#!/usr/bin/env Rscript
# revision_experiments/tr3/25_greedy_prefix.R  (2026-10-25)
# Real data: does stopping the first greedy selection at the first covering prefix (the
# outdegree greedy algorithm the paper used to describe) change any CCD result, compared with
# dom.method = "greedy2" as shipped (every vertex ranked by outdegree)? R/ccds/ is not edited:
# the prefix version replaces dominate.mat.greedy2 in the global environment for variant "prefix".
# USAGE (clone root):
#   Rscript revision_experiments/tr3/25_greedy_prefix.R <UN-CCDs|RK-CCDs> <orig|prefix> <out.csv> [datasets]
# The 2026-10-25 outputs are in results/tr3/greedy_prefix/ (un_orig.csv, un_prefix.csv, rk_prefix.csv).
args <- commandArgs(trailingOnly = TRUE)
METHOD <- args[1]; VARIANT <- args[2]; OUT <- args[3]
DATASETS <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else
  c("iris", "seeds", "user_knowledge", "wholesale", "asymmetric", "R15", "D31")
source(here::here("revision_experiments/tr3/tr3_lib.R"))
settings <- lapply(setNames(DATASETS, DATASETS), function(ds) tr3_real_setting(METHOD, ds))
tr3_source_method(METHOD, "realdata")
real <- tr3_load_real(DATASETS)

shipped_greedy2 <- dominate.mat.greedy2
prefix_greedy2 <- function(A) {
  S <- shipped_greedy2(A)                      # every vertex, by outdegree
  if (!is.matrix(A)) A <- matrix(A, nrow = 1)
  covered <- rep(FALSE, nrow(A))
  for (m in seq_along(S)) {
    covered <- covered | (A[S[m], ] == TRUE)
    if (all(covered)) break
  }
  S[seq_len(m)]                                # shortest covering prefix
}
if (VARIANT == "prefix") assign("dominate.mat.greedy2", prefix_greedy2, envir = globalenv())
stopifnot(VARIANT %in% c("orig", "prefix"))

for (ds in DATASETS) {
  s <- settings[[ds]]
  D <- real[[ds]]; X <- D[, -ncol(D)]; Y <- D[, ncol(D)]
  t0 <- proc.time()[["elapsed"]]
  res <- if (METHOD == "UN-CCDs") {
    UNCCD_clustering(datax = X, simul = tr3_load_table(s$table), method = s$direction)
  } else {
    RKCCD_clustering(datax = X, simul = tr3_load_table(s$table), quant = as.numeric(s$quant))
  }
  lab <- res$label
  k <- length(unique(lab))
  row <- data.frame(method = METHOD, variant = VARIANT, dataset = ds,
                    n_D = if (!is.null(res$D)) length(res$D) else NA,
                    n_IntD = if (!is.null(res$Int.D)) length(res$Int.D) else NA,
                    ari = sprintf("%.15g", mclust::adjustedRandIndex(lab, Y)),
                    sil = sprintf("%.15g", if (k > 1) mean(cluster::silhouette(lab, dist(X))[, 3]) else NA),
                    k_hat = k, labels = paste(lab, collapse = " "),
                    seconds = round(proc.time()[["elapsed"]] - t0, 2))
  write.table(row, OUT, sep = ",", append = file.exists(OUT), col.names = !file.exists(OUT),
              row.names = FALSE, qmethod = "double")
  cat(METHOD, VARIANT, ds, "k =", k, "ARI =", row$ari, "|D| =", row$n_D, "|Int.D| =", row$n_IntD, "\n")
}
