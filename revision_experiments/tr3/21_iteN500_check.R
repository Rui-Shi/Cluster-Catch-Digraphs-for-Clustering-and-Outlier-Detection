#!/usr/bin/env Rscript
# revision_experiments/tr3/21_iteN500_check.R
#
# Checks of the iteN500 rerun (ITEN500_PROTOCOL.md): RK-CCDs, d = 3, n = 500,
# uniform and Gaussian clusters, 500 repetitions instead of the scripts' 50.
# Stops with an error if any check fails; nothing from the rerun is to be used
# before all of them pass.
#
# Checks, per distribution
#   1. Complete: 500 replicates for each of 2, 3 and 5 clusters, no NA.
#   2. Reproduction: the 2-cluster data sets 1-50 are the published 50
#      (same seed, generated first), so the means over replicates 1-50 must
#      equal the 2-cluster line of the published log to all 15 printed digits:
#        uniform   simulations/clustering/RK_CCDs/Simulations/Uniform/3d/slurm-1442239.out
#        Gaussian  simulations/clustering/RK_CCDs/Simulations/Gaussian/3d/slurm-1443038.out
#   3. Pairing: every data fingerprint equals that of the same replicate in
#      the UN-CCD and KS-CCD cell (data regenerated here with their settings).
#
# USAGE   Rscript revision_experiments/tr3/21_iteN500_check.R
# OUTPUT  revision_experiments/results/tr3/iteN500_check.csv (rewritten):
#         the new 500-repetition means next to the values now in Tables S3/S7

source(here::here("revision_experiments/tr3/tr3_lib.R"))
tr3_source_method("UN-CCDs", "simulation")   # for the data generators only

CELLS <- data.frame(
  dist = c("uniform", "gaussian"),
  log = c("simulations/clustering/RK_CCDs/Simulations/Uniform/3d/slurm-1442239.out",
          "simulations/clustering/RK_CCDs/Simulations/Gaussian/3d/slurm-1443038.out"),
  table = c("S3", "S7"), stringsAsFactors = FALSE)
# Supplementary_Material.tex, tables tab:uni_RK-CCDs_Apdx and tab:Gau_RK-CCDs_Apdx,
# n = 500, d = 3 (the current 50-repetition values)
PAPER <- utils::read.csv(colClasses = "character", text = "
dist,n_clusters,ari,sil,sr
uniform,2,1.000,0.669,1.000
uniform,3,1.000,0.654,1.000
uniform,5,0.999,0.600,1.000
gaussian,2,0.990,0.633,1.000
gaussian,3,0.987,0.624,1.000
gaussian,5,0.976,0.581,1.000")

same15 <- function(got, expected) isTRUE(as.numeric(sprintf("%.15g", got)) == as.numeric(sprintf("%.15g", as.numeric(expected))))
keys <- c("variant", "method", "dist", "d", "n", "n_clusters")
reps <- tr3_read_done(tr3_path("sim_replicates.csv"), tr3_path("sim_done.csv"), keys, "rep")
if (is.null(reps)) stop("21_iteN500_check.R: no results yet")

out <- list()
for (i in seq_len(nrow(CELLS))) {
  C <- CELLS[i, ]
  r <- reps[reps$variant == "iteN500" & reps$method == "RK-CCDs" & reps$dist == C$dist &
              reps$d == "3" & reps$n == "500", ]
  # 1. complete
  for (N in c(2, 3, 5)) {
    rn <- r[r$n_clusters == as.character(N), ]
    if (!identical(sort(as.integer(rn$rep)), 1:500) || any(is.na(as.numeric(rn$ari))) ||
        any(is.na(as.numeric(rn$sil))) || any(is.na(as.numeric(rn$k_hat)))) {
      stop(sprintf("21_iteN500_check.R: %s, %d clusters is not complete (%d rows)", C$dist, N, nrow(rn)))
    }
  }
  # 2. reproduction of the published 2-cluster line from replicates 1-50
  txt <- readLines(here::here(C$log), warn = FALSE)
  m <- regmatches(txt, regexec(paste0("The mean ARI is ([^ ]+) , The mean average silhouette index is ([^ ]+) , ",
                                      "and, the percentatage of times identifying the correct number of clusters is ([^ ]+) %"), txt))
  m <- Filter(length, m)
  stopifnot(length(m) == 3)
  r2 <- r[r$n_clusters == "2" & as.integer(r$rep) <= 50, ]
  got <- c(mean(as.numeric(r2$ari)), mean(as.numeric(r2$sil)), mean(as.numeric(r2$k_hat) == 2) * 100)
  ok <- mapply(same15, got, m[[1]][2:4])
  cat(sprintf("%s: replicates 1-50, 2 clusters: ARI %s / Sil %s / SR %s; log %s / %s / %s -> %s\n", C$dist,
              sprintf("%.15g", got[1]), sprintf("%.15g", got[2]), sprintf("%.15g", got[3]),
              m[[1]][2], m[[1]][3], m[[1]][4], if (all(ok)) "match" else "MISMATCH"))
  if (!all(ok)) stop("21_iteN500_check.R: check 2 failed for ", C$dist)
  # 3. pairing with the UN-CCD / KS-CCD data
  for (other in c("UN-CCDs", "KS-CCDs")) {
    s <- tr3_sim_setting(other, C$dist, "3", "500")
    stopifnot(s$iteN == "500")
    dat <- tr3_sim_data(s)
    for (N in c(2, 3, 5)) {
      rn <- r[r$n_clusters == as.character(N), ]
      rn <- rn[order(as.integer(rn$rep)), ]
      fp <- vapply(dat[[as.character(N)]]$data, tr3_fp, numeric(1))
      if (!identical(as.numeric(rn$data_fp), fp)) {
        stop(sprintf("21_iteN500_check.R: check 3 failed: %s, %d clusters, data differ from %s", C$dist, N, other))
      }
    }
  }
  cat(sprintf("%s: all 1500 data sets equal the UN-CCD and KS-CCD data sets\n", C$dist))
  # summary
  for (N in c(2, 3, 5)) {
    rn <- r[r$n_clusters == as.character(N), ]
    p <- PAPER[PAPER$dist == C$dist & PAPER$n_clusters == as.character(N), ]
    ari <- mean(as.numeric(rn$ari)); sil <- mean(as.numeric(rn$sil)); sr <- mean(as.numeric(rn$k_hat) == N)
    out[[length(out) + 1]] <- data.frame(
      table = C$table, dist = C$dist, n_clusters = N,
      ari_500 = fmt17(ari), sil_500 = fmt17(sil), sr_500 = fmt17(sr),
      ari_500_3dp = sprintf("%.3f", ari), sil_500_3dp = sprintf("%.3f", sil), sr_500_3dp = sprintf("%.3f", sr),
      ari_paper = p$ari, sil_paper = p$sil, sr_paper = p$sr, stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, out)
utils::write.csv(res, tr3_path("iteN500_check.csv"), row.names = FALSE)
print(res[, c("table", "dist", "n_clusters", "ari_500_3dp", "ari_paper", "sil_500_3dp", "sil_paper",
              "sr_500_3dp", "sr_paper")], row.names = FALSE)
cat("21_iteN500_check.R: all checks pass; written", tr3_path("iteN500_check.csv"), "\n")
