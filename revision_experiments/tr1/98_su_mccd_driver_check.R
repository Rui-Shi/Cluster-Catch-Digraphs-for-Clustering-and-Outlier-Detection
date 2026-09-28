#!/usr/bin/env Rscript
# 98_su_mccd_driver_check.R -- rerun one SU-MCCD two-cluster simulation cell
# (Uniform or Gaussian) with an explicit dimension and RK quantile table.
#
# Why: two groups of published SU-MCCD drivers do not run what their file
# names say.
#   simulations/outlier_detection/SU-MCCDs/Simulation/{Uniform,Gaussian}/3d/3d_2cls_n1000_cont5%.R
#     set d = 2 and load the 2d table; their logs (slurm-1027222, -1014603)
#     equal the d = 2, n = 1000 logs (slurm-1027221, -1014601) exactly.
#   simulations/outlier_detection/SU-MCCDs/Simulation/{Uniform,Gaussian}/50d/50d_2cls_n50_cont5%.R
#     load RK-test-simul_2d_99%.RData, which holds only the 0.99 level, and
#     ask for 0.999; simul$quan[["0.999"]] is NULL, NULL[j, ] is NULL, the
#     test never rejects, and every covering radius is set to 0.
# Everything else is copied from those drivers: generator, cluster geometry,
# contamination 5%, set.seed(123) before generating all data sets,
# S_min = max(0.04, cont/2) for d <= 5 and 0 otherwise, level 0.99 for d < 10
# and 0.999 otherwise, low.num = 2, LDOF k = max(d, 5). n2 is
# round(n*(1-cont)*0.5) in the n = 50 and n = 1000 drivers and one less in
# the others (--n2adj). The drivers source methods/outlier_detection/SU-MCCDs.R,
# which defines SUMCCD_outlier; they call it by its earlier name, MFCCD_outlier
# (the file was M-FCCDs.R on the cluster; commit 389ec3d), so it is aliased
# here. --impl=ldof uses the LDOF variant in
# simulations/outlier_detection/SU-MCCDs/M-FCCDs_LDOF.R instead.
#
# Output (appended per chunk): results/tr1/su_driver_fix/98_<setting>_d<d>_n<n>_<table>.csv
# Usage: Rscript 98_su_mccd_driver_check.R --setting=Uniform --d=50 --n=50
#          --table=RK-test-simul_50d_999%.RData [--n2adj=0] [--reps=1000] [--cores=22]

suppressPackageStartupMessages({ library(here); library(parallel); library(MASS) })
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  h <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", name, "="), "", h) else default
}
SETTING <- match.arg(arg("setting", "Uniform"), c("Uniform", "Gaussian"))
d       <- as.integer(arg("d", NA)); n <- as.integer(arg("n", NA))
TABLE   <- arg("table", NA)
N2ADJ   <- as.integer(arg("n2adj", 0))
REPS    <- as.integer(arg("reps", 1000))
CORES   <- as.integer(arg("cores", max(1L, detectCores() - 2L)))
CHUNK   <- as.integer(arg("chunk", CORES))
stopifnot(!is.na(d), !is.na(n), !is.na(TABLE))

IMPL <- match.arg(arg("impl", "su"), c("su", "ldof"))
METHOD_FILE <- here::here(c(su = "methods/outlier_detection/SU-MCCDs.R",
                            ldof = "simulations/outlier_detection/SU-MCCDs/M-FCCDs_LDOF.R")[[IMPL]])
source(METHOD_FILE)  # also defines rpoisball.unit
if (IMPL == "su") MFCCD_outlier <- SUMCCD_outlier
source(here::here("R/general_functions/count.R"))

cont <- 0.05
min.cls <- if (d <= 5) max(0.04, cont / 2) else 0
quant <- if (d < 10) 0.99 else 0.999
iteN <- 1000; cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1)); mu <- colMeans(rbind(mu1, mu2))
n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) + N2ADJ; n0 <- round(n * cont)
if (SETTING == "Gaussian") { noise_level <- 0.01; sigma <- 1 / sqrt(qchisq(1 - noise_level, d)) }

set.seed(123)
data.list <- lapply(1:iteN, function(x) {
  if (SETTING == "Uniform") {
    data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
    data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  } else {
    data1 <- mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
    data2 <- mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  }
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1
    }
  }
  rownames(outlier) <- NULL
  rbind(data1, data2, outlier)
})
# count() uses the global n, n0 and d, with n left at its nominal value as in the drivers

TABPATH <- here::here("R/RK-test_quantile", TABLE)
OUTDIR <- here::here("revision_experiments/results/tr1/su_driver_fix")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, sprintf("98_%s_d%d_n%d_%s_%s.csv", SETTING, d, n, sub("\\.RData$", "", TABLE), IMPL))
done <- if (file.exists(OUT)) read.csv(OUT)$rep else integer(0)
todo <- setdiff(seq_len(min(REPS, iteN)), done)
cat(sprintf("SU-MCCD %s d=%d n=%d (n1=%d n2=%d n0=%d) table=%s quant=%s min.cls=%s: %d replicates to run, %d cores\n",
            SETTING, d, n, n1, n2, n0, TABLE, quant, min.cls, length(todo), CORES))

cl <- makeCluster(CORES)
invisible(clusterCall(cl, function(tab, mf, impl) {
  suppressPackageStartupMessages({ library(MASS); library(cluster); library(igraph) })
  source(mf)
  if (impl == "su") assign("MFCCD_outlier", SUMCCD_outlier, envir = globalenv())
  source(here::here("R/ccds/quantile_table.R"))   # a missing table is generated on the spot
  load_quantile_table(tab, envir = globalenv())
  NULL
}, TABPATH, METHOD_FILE, IMPL))
clusterExport(cl, c("min.cls", "quant"))

t0 <- Sys.time()
for (chunk in split(todo, ceiling(seq_along(todo) / CHUNK))) {
  res <- parLapply(cl, data.list[chunk], function(x) MFCCD_outlier(x, simul, min.cls = min.cls, quant = quant))
  outlier.result <- vector("list", iteN); outlier.result[chunk] <- res
  rows <- do.call(rbind, lapply(chunk, function(i) {
    v <- count(i)
    data.frame(rep = i, TPR = v[["success_rate"]], TNR = v[["true_positive"]], BA = v[["BA"]], F2 = v[["F2_score"]])
  }))
  write.table(rows, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT), col.names = !file.exists(OUT))
  cat(sprintf("  reps %d-%d done (%.1f min elapsed)\n", min(chunk), max(chunk),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
stopCluster(cl)

x <- read.csv(OUT); x <- x[!duplicated(x$rep), ]
cat(sprintf("RESULT SU-MCCD %s d=%d n=%d %s reps=%d: TPR %.6f TNR %.6f BA %.6f F2 %.6f\n",
            SETTING, d, n, TABLE, nrow(x), mean(x$TPR), mean(x$TNR), mean(x$BA), mean(x$F2)))
