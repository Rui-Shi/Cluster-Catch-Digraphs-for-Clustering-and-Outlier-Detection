#!/usr/bin/env Rscript
# revision_experiments/tr3/23_noise.R
#
# The noise experiment of the UN-CCD clustering paper: main text Section
# "Clustering with noise" (tab:Noise_Level_Clustering), Supplement S1.3
# (tab:Noise_Level_Clustering_Apdx). The original driver is not in the
# repository, so this one was written from the supplement's settings
# (README.md, "Noise experiment (2026-10-24)" says which setting comes from
# where and how the result was validated). Same layout as 10_sim.R.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/23_noise.R <d> [variant] [r_noise] [center] [iteN]
#     <d>        3 or 10
#     [variant]  published (default): the committed code.
#                oldindex: d = 10 only, for validation. nnccd.radi.descend.one is
#                replaced in the global environment (main process and workers) by a
#                copy that reads env.*.rev[j + 2] as the code did until 2026-10-23;
#                nothing in R/ccds/ is edited.
#                allpoints: experiment 8 (NEW_EXPERIMENTS_PROTOCOL.md). The published
#                computation plus the all-points choice of the number of clusters scored
#                as the baselines of 37_sim_baselines.R (ari_all, sil200_all, nlab_all).
#                Every published value is first compared with noise_replicates.csv (gate).
#                Writes noise_allpoints_replicates.csv and noise_allpoints_done.csv.
#     [r_noise]  radius of the ball of noise points (default 7)
#     [center]   sweep (default): the ball is centred on the mean of the cluster
#                centres. recycle: the noise branch of Gaussian_cls() as written, which
#                adds that vector by recycling down an n_noise x d matrix and so
#                misplaces the noise unless all its coordinates are equal.
#     [iteN]     replicates per noise level (default 1000)
#   TR3_CORES overrides the number of workers (default detectCores()). The results
#   do not depend on it: the data are generated on the main process and the
#   clustering is deterministic.
#
# Settings (supplement S1.3, items i to vii)
#   d = 3, 10; noise levels l = 3%, 5%, 10%, 20% (n_noise = round(200 * l));
#   200 cluster points + n_noise noise points; 3 equal clusters; centres
#   (3,..,3), (9,3,..,3), (3,9,3,..,3); covariance Delta * I_d with the radius
#   factor ~ U[0.8, 1.2], as Gaussian_cls() draws it; 1000 repetitions.
# Not stated in the supplement, taken from the matching non-noise driver
# (UN_CCDs/Simulations/Gaussian/<d>d/<d>d_235clusters_n200.R, row of
# settings_simulation.csv): seed 2024, sigma 1, cluster sizes 66, 67, 67
# (first of three is round(n/3) - 1), quantile table NN-test-simul_3d_90% with the
# ascending search at d = 3 and NN-test-simul_10d_99% with the descending search
# at d = 10, low.num 3, niter 1000, greedy2 dominating set, silhouette-based
# choice of the number of clusters.
# Not stated anywhere and chosen here: r_noise, the centre of the noise ball,
# and what ARI and Sil are computed on (see below). These choices were made by
# looking at the published d = 3 row (the demo figure suggests a radius of about
# 6 in 2 dimensions), so they are inferences, not a reproduction. The check
# against the published table FAILED (README.md, "Noise experiment
# (2026-10-24)"): the published numbers are not reproduced at d = 10 with the old
# index. On 2026-10-24 the authors decided to adopt this script's rows
# (defaults: r_noise 7, sweep, 1000; d = 3 published, d = 10 published) for the paper.
#
# What is computed per replicate
#   The clustering runs on all 200 + n_noise points, as UNCCD_clustering() does
#   (nnccd_clustering_quantile() then nnccd.silhouette()), except that the
#   silhouette used to pick the number of clusters, and the reported Sil, are
#   computed on the 200 cluster points only (nnccd.silhouette(ind = 1:200)), and
#   ARI compares the labels of those 200 points with the truth. Noise points
#   are clustered but are not scored. With the silhouette over all points the
#   reported Sil falls from 0.6 to 0.5 as l goes from 3% to 20%, which the
#   published table does not show. k_hat is the chosen number of clusters, SR the
#   share of replicates with k_hat = 3. Columns sil_all and k_all give the
#   all-points selection for comparison.
#
# OUTPUT (revision_experiments/results/tr3/)
#   noise_replicates.csv  one row per replicate: variant, d, noise_level, r_noise,
#                         center, iteN, rep, ari, sil, k_hat, sil_all, k_all,
#                         seconds, data_fp
#   noise_done.csv        one row per finished cell; reruns skip these cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: Rscript revision_experiments/tr3/23_noise.R <d> [variant] [r_noise] [center] [iteN]")
}
source(here::here("revision_experiments/tr3/tr3_lib.R"))

D       <- args[1]
VARIANT <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "published"
R_NOISE <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 7
CENTER  <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "sweep"
ITEN    <- if (length(args) >= 5 && nzchar(args[5])) as.integer(args[5]) else 1000L
LEVELS  <- c(0.03, 0.05, 0.10, 0.20)
stopifnot(D %in% c("3", "10"), VARIANT %in% c("published", "oldindex", "allpoints"), CENTER %in% c("sweep", "recycle"),
          is.finite(R_NOISE), R_NOISE > 0, ITEN >= 1)
if (VARIANT == "oldindex" && D != "10") stop("23_noise.R: oldindex applies to the descending search, d = 10")

DONE_CSV <- tr3_path("noise_done.csv")
REP_CSV  <- tr3_path("noise_replicates.csv")
# TR3_SMOKE = m > 0: generate all data as usual, run only the first m replicates of each level,
# print (and gate), write nothing
SMOKE <- as.integer(Sys.getenv("TR3_SMOKE", "0"))
if (VARIANT == "allpoints") {
  stopifnot(ITEN == 1000L)  # the published data: 1000 replicates per level, drawn in order
  PUB_CSV  <- REP_CSV
  DONE_CSV <- tr3_path("noise_allpoints_done.csv")
  REP_CSV  <- tr3_path("noise_allpoints_replicates.csv")
  PUB <- utils::read.csv(PUB_CSV, colClasses = "character")
  PUB <- PUB[PUB$variant == "published" & PUB$d == D & PUB$r_noise == as.character(R_NOISE) &
               PUB$center == CENTER, , drop = FALSE]
  stopifnot(nrow(PUB) == 4L * ITEN)
}

# the non-noise Gaussian n = 200 row gives table, direction, seed, sigma, radius range, sizes
s <- tr3_sim_setting("UN-CCDs", "gaussian", D, "200")
tr3_source_method("UN-CCDs", "simulation")
par <- list(simul = tr3_load_table(s$table), direction = s$direction)
stopifnot(s$direction == if (D == "3") "ascend" else "descend")

if (VARIANT == "oldindex") {
  src <- deparse(nnccd.radi.descend.one)
  stopifnot(sum(grepl("env.ave.rev[j + 1]", src, fixed = TRUE)) == 1,
            sum(grepl("env.med.rev[j + 1]", src, fixed = TRUE)) == 1)
  src <- gsub("env.ave.rev[j + 1]", "env.ave.rev[j + 2]", src, fixed = TRUE)
  src <- gsub("env.med.rev[j + 1]", "env.med.rev[j + 2]", src, fixed = TRUE)
  old_descend <- eval(parse(text = src), envir = globalenv())
  environment(old_descend) <- globalenv()
  assign("nnccd.radi.descend.one", old_descend, envir = globalenv())
}

# which index the function that runs uses: 1 = j + 2 (old), 0 = j + 1 (committed)
descend_is_old <- function() {
  as.integer(any(grepl("env.ave.rev[j + 2]", deparse(nnccd.radi.descend.one), fixed = TRUE)))
}

NC <- 200L
n_cls <- tr3_sizes(3, NC, as.numeric(s$n1_adj_3cls))
dd <- as.numeric(D)
MU <- rbind(rep(3, dd), c(9, rep(3, dd - 1)), c(3, 9, rep(3, dd - 2)))
sig <- as.numeric(s$sig); rad <- c(as.numeric(s$r_min), as.numeric(s$r_max))
label_T <- rep(0:2, n_cls)

# One data set. The clusters come from Gaussian_cls() unchanged. The noise is
# what its noise branch draws (rpoisball_unit(n_noise, d) * r_noise, drawn after
# the clusters), centred on colMeans(MU) with sweep(), or as written (recycle).
noise_set <- function(i, n_noise) {
  x <- Gaussian_cls(n_cls = n_cls, d = dd, centers = MU, sigma = sig, noise = FALSE, rad_range = rad)
  z <- rpoisball_unit(n_noise, dd) * R_NOISE
  z <- if (CENTER == "sweep") sweep(z, 2, colMeans(MU), "+") else z + colMeans(MU)
  rbind(x, z)
}

# One replicate, on a worker.
tr3_noise_one <- function(x, label_T, nc, par) {
  t0 <- proc.time()[["elapsed"]]
  g <- nnccd_clustering_quantile(x, low.num = 3, quantile = "lower", method = par$direction,
                                 dom.method = "greedy2", simul = par$simul, niter = 1000)
  ri <- nnccd.silhouette(g, x, cls = NULL, ind = seq_len(nc), lenDlimit = Inf)
  ra <- nnccd.silhouette(g, x, cls = NULL, ind = NULL, lenDlimit = Inf)
  # experiment 8: the all-points choice, scored as 37_sim_baselines.R scores the baselines
  la <- ra$label[seq_len(nc)]
  nd <- length(unique(la))
  s200 <- if (nd >= 2 && nd <= nc - 1) mean(cluster::silhouette(la, stats::dist(x[seq_len(nc), , drop = FALSE]))[, 3]) else NA_real_
  c(ari = mclust::adjustedRandIndex(ri$label[seq_len(nc)], label_T), sil = ri$si, k_hat = ri$si.ind,
    sil_all = ra$si, k_all = ra$si.ind, seconds = proc.time()[["elapsed"]] - t0,
    old = descend_is_old(),
    ari_all = mclust::adjustedRandIndex(la, label_T), sil200_all = s200, nlab_all = length(unique(ra$label)))
}

t_start <- Sys.time()
cores <- tr3_cores()
cat(sprintf("23_noise.R: d=%s variant=%s r_noise=%g center=%s iteN=%d cores=%d table=%s direction=%s\n",
            D, VARIANT, R_NOISE, CENTER, ITEN, cores, s$table, s$direction))

# All four noise levels are always generated, in this order, from one seed, so
# that skipping a finished cell cannot shift the data of the later ones.
set.seed(as.numeric(s$seed))
data <- lapply(LEVELS, function(l) {
  nn <- round(NC * l)
  lapply(seq_len(ITEN), noise_set, n_noise = nn)
})
names(data) <- as.character(LEVELS)

cl <- makeCluster(cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

# foreach is called at top level, as in 10_sim.R, so that it exports the global
# functions (including the replaced nnccd.radi.descend.one) to the workers.
for (l in LEVELS) {
  key <- list(variant = VARIANT, d = D, noise_level = as.character(l), r_noise = as.character(R_NOISE),
              center = CENTER)
  if (has_result(DONE_CSV, key)) {
    cat(sprintf("  l = %g: done, skipped\n", l))
    next
  }
  data_list <- data[[as.character(l)]]
  if (SMOKE > 0) data_list <- data_list[seq_len(SMOKE)]
  t0 <- Sys.time()
  r <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% tr3_noise_one(x, label_T, NC, par)
  r <- do.call(rbind, r)
  if (VARIANT == "oldindex" && !all(r[, "old"] == 1)) stop("23_noise.R: the old index did not reach every worker")
  if (VARIANT %in% c("published", "allpoints") && D == "10" && !all(r[, "old"] == 0)) stop("23_noise.R: an old index is loaded")

  rows <- data.frame(key, iteN = ITEN, rep = seq_len(nrow(r)),
                     ari = fmt17(r[, "ari"]), sil = fmt17(r[, "sil"]), k_hat = fmt17(r[, "k_hat"]),
                     sil_all = fmt17(r[, "sil_all"]), k_all = fmt17(r[, "k_all"]),
                     seconds = fmt17(r[, "seconds"]),
                     data_fp = fmt17(vapply(data_list, tr3_fp, numeric(1))),
                     stringsAsFactors = FALSE)
  if (VARIANT == "allpoints") {
    # gate: every published value of these replicates equals the stored published row
    pub <- PUB[PUB$noise_level == as.character(l), , drop = FALSE]
    pub <- pub[match(as.character(rows$rep), pub$rep), , drop = FALSE]
    cols <- c("data_fp", "ari", "sil", "k_hat", "sil_all", "k_all")
    bad <- vapply(cols, function(cc) sum(is.na(pub[[cc]]) | pub[[cc]] != rows[[cc]]), numeric(1))
    cat(sprintf("  l = %g: gate on %d replicates, differences %s\n", l, nrow(rows),
                paste(names(bad), bad, sep = "=", collapse = " ")))
    if (any(bad > 0)) stop("23_noise.R allpoints: the published values are not reproduced (gate)")
    rows$ari_all    <- fmt17(r[, "ari_all"])
    rows$sil200_all <- fmt17(r[, "sil200_all"])
    rows$nlab_all   <- fmt17(r[, "nlab_all"])
    rows$gate       <- "equal"
    rows <- rows[, c(names(key), "iteN", "rep", "ari", "sil", "k_hat", "sil_all", "k_all",
                     "ari_all", "sil200_all", "nlab_all", "seconds", "data_fp", "gate")]
    cat(sprintf("  l = %g: all points: mean ARI %.4f, SR %.3f, nlab_all != k_all in %d\n", l,
                mean(r[, "ari_all"]), mean(r[, "nlab_all"] == 3), sum(r[, "nlab_all"] != r[, "k_all"])))
  }
  if (SMOKE > 0) {
    cat(sprintf("  l = %g: smoke, %d replicates, %.1f s wall, nothing written\n", l, nrow(rows),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    next
  }
  tr3_append(REP_CSV, rows)
  tr3_append(DONE_CSV, data.frame(key, iteN = ITEN, n_rows = nrow(rows), cores = cores,
                                  wall_secs = fmt17(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                                  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                  stringsAsFactors = FALSE))
  print(paste0("UN-CCDs d=", D, " l=", l, ": The mean ARI is ", mean(r[, "ari"]),
               ", the mean average silhouette index is ", mean(r[, "sil"]),
               ", and the share of times identifying the correct number of clusters is ",
               mean(r[, "k_hat"] == 3)))
}

print(Sys.time() - t_start)
