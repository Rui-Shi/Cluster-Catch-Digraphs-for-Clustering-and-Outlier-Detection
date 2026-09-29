#!/usr/bin/env Rscript
# 97_sun_ns_d10_rerun.R -- rerun SUN-MCCD in the three Neyman-Scott settings
# (Matern, Thomas, mixed) at d = 10.
#
# Why: the published drivers
#   simulations/outlier_detection/SUN-MCCDs/Simulation/Complex_Clusters/{Matern,Thomas,Mix}/10d_clx_cls.R
# load NN-test-simul_10d_99%.RData (alpha = 1%), while the stated schedule
# (main text, tab:alpha_sim) and every other SUN-MCCD d = 10 driver use the
# 0.1% table. Pass "99" reproduces the published run (check against
# slurm-1011335/1011342/1011349.out); pass "999" is the corrected run.
#
# Everything else is copied from the drivers: generator Uni.Gau_cls, the
# settings below, set.seed(1234) before generating all 1000 data sets,
# S_min = 0, descending radius search, low.num = 3, LDOF k = max(d, 5).
# The drivers source methods/outlier_detection/SUN-MCCD.R, which does not
# define the MFNNCCD_outlier they call; the definition they ran is in
# simulations/outlier_detection/SUN-MCCDs/M-FNNCCDs_LDOF.R, sourced here.
# The radius search is the incremental one (R/ccds/UN_CCD.R), validated
# bit-identical to the recomputing search (tr1/92, tr1/93).
#
# Output (appended per chunk): results/tr1/ns_sun_d10/97_<process>_<level>.csv
# Usage: Rscript 97_sun_ns_d10_rerun.R --process=Matern --level=999 [--cores=22] [--chunk=44]

suppressPackageStartupMessages({ library(here); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default) {
  h <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", name, "="), "", h) else default
}
PROCESS <- match.arg(arg("process", "Matern"), c("Matern", "Thomas", "Mix"))
LEVEL   <- match.arg(arg("level", "999"), c("99", "999"))
CORES   <- as.integer(arg("cores", max(1L, detectCores() - 2L)))
CHUNK   <- as.integer(arg("chunk", 2L * CORES))

# ---- settings copied from the drivers --------------------------------------
# --d selects a sibling driver (d = 2, 3, 5, 20) to check that this script
# reproduces a published cell whose level was never in doubt; the default is
# d = 10, the cell being corrected. Per-d table, search direction, S_min and
# cluster intensity index are those of {Matern,Thomas,Mix}/<d>d_clx_cls.R.
d <- as.integer(arg("d", "10")); iteN <- 1000; cont <- 0.10
SET <- list(`2`  = list(tab = "2d_85%",  method = "ascend",  min.cls = max(0.04, cont / 2), idx = 1),
            `3`  = list(tab = "3d_90%",  method = "ascend",  min.cls = max(0.04, cont / 2), idx = 2),
            `5`  = list(tab = "5d_95%",  method = "ascend",  min.cls = 0, idx = 3),
            `10` = list(tab = sprintf("10d_%s%%", LEVEL), method = "descend", min.cls = 0, idx = 4),
            `20` = list(tab = "20d_999%", method = "descend", min.cls = 0, idx = 5))[[as.character(d)]]
min.cls <- SET$min.cls; method <- SET$method
if (d != 10) LEVEL <- sub("^.*_", "", sub("%$", "", SET$tab))
ratio_file <- c(Matern = "ratio1.R", Thomas = "ratio2.R", Mix = "ratio3.R")[[PROCESS]]
kappa1 <- c(Matern = 6, Thomas = 0, Mix = 3)[[PROCESS]]
kappa2 <- c(Matern = 0, Thomas = 6, Mix = 3)[[PROCESS]]

# the method file comes first, as in the drivers: it also defines rpoisball.unit,
# which the generator calls
IMPL <- match.arg(arg("impl", "sun"), c("sun", "ldof"))
METHOD_FILE <- here::here(c(sun = "methods/outlier_detection/SUN-MCCD.R",
                            ldof = "simulations/outlier_detection/SUN-MCCDs/M-FNNCCDs_LDOF.R")[[IMPL]])
source(METHOD_FILE)
if (IMPL == "sun") MFNNCCD_outlier <- SUNMCCD_outlier
source(here::here("R/general_functions/count.R"))
source(here::here("R/general_functions/Uni-Gau_cls.R"))
source(here::here("R/general_functions", ratio_file))
mu1 <- ratio[SET$idx]; expand1 <- 0; r <- 0.1
scale <- 0.005; mu2 <- ratio[SET$idx]; expand2 <- 0; slen <- 1; kappa_O <- 20

set.seed(1234)
data.listNum <- lapply(1:iteN, function(x) {
  data_simu <- Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 <- data_simu$Matérn_children
  cls2 <- data_simu$Thomas_children
  cls3 <- data_simu$noise
  outlier <- data_simu$Outlier
  c(data = list(rbind(cls1, cls2, cls3, outlier)), num = list(data_simu$num))
})
data.list <- lapply(data.listNum, function(z) z$data)
data.num  <- lapply(data.listNum, function(z) z$num)

TABLE <- here::here(sprintf("R/NN-test_quantile/NN-test-simul_%s.RData", SET$tab))
OUTDIR <- here::here("revision_experiments/results/tr1/ns_sun_d10")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, sprintf("97_%s_d%d_%s_%s.csv", PROCESS, d, LEVEL, IMPL))
done <- if (file.exists(OUT)) read.csv(OUT)$rep else integer(0)
todo <- setdiff(seq_len(iteN), done)
cat(sprintf("%s level=%s%% table=%s: %d of %d replicates to run, %d cores\n",
            PROCESS, LEVEL, basename(TABLE), length(todo), iteN, CORES))

cl <- makeCluster(CORES)
invisible(clusterCall(cl, function(tab, mf, impl) {
  suppressPackageStartupMessages({ library(MASS); library(cluster); library(igraph) })
  source(mf)
  if (impl == "sun") assign("MFNNCCD_outlier", SUNMCCD_outlier, envir = globalenv())
  source(here::here("R/ccds/quantile_table.R"))   # a missing table is generated on the spot
  load_quantile_table(tab, envir = globalenv())
  NULL
}, TABLE, METHOD_FILE, IMPL))
clusterExport(cl, c("min.cls", "method"))

t0 <- Sys.time()
for (chunk in split(todo, ceiling(seq_along(todo) / CHUNK))) {
  res <- parLapply(cl, data.list[chunk], function(x) MFNNCCD_outlier(x, simul, min.cls = min.cls, method = method))
  outlier.result <- vector("list", iteN); outlier.result[chunk] <- res
  rows <- do.call(rbind, lapply(chunk, function(i) {
    v <- count1(i)
    data.frame(rep = i, n = sum(data.num[[i]]), n0 = data.num[[i]][2],
               TPR = v[["success_rate"]], FPR = v[["false_positive"]], F2 = v[["F2_score"]])
  }))
  write.table(rows, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT), col.names = !file.exists(OUT))
  cat(sprintf("  reps %d-%d done (%.1f min elapsed)\n", min(chunk), max(chunk),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
stopCluster(cl)

x <- read.csv(OUT)
x <- x[!duplicated(x$rep), ]
stopifnot(nrow(x) == iteN)
m <- c(TPR = mean(x$TPR), TNR = 1 - mean(x$FPR), F2 = mean(x$F2))
m <- c(m, BA = (m[["TPR"]] + m[["TNR"]]) / 2)
cat(sprintf("RESULT %s d=%d level=%s%%: TPR %.6f TNR %.6f BA %.6f F2 %.6f\n",
            PROCESS, d, LEVEL, m[["TPR"]], m[["TNR"]], m[["BA"]], m[["F2"]]))
