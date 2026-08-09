#!/usr/bin/env Rscript
# revision_experiments/09_wp3_synthetic.R
#
# WP3 threshold-sensitivity study, synthetic part (Task T7 Phase A).
# Reviewers (R1 #3, R2 3rd comment) find the score-cutoff calibration "too
# empirical"; this script quantifies how BA and F2 respond to perturbing the
# calibrated cutoff, for 3 synthetic settings x 4 CCD-based OS methods.
#
# Settings (data-generating code copied verbatim from the first-cycle
# scripts; only change: seed varies per replicate, seed = base + rep, and one
# dataset is generated per call instead of a length-iteN lapply):
#   gaussian: simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/
#             Gaussian/10d/10d_2cls_n500_cont5%.R   (d=10, 2 cls, n=500, 5%)
#   uniform:  .../RKCCD_OOS_IOS/Simulation/Uniform/10d/10d_2cls_n500_cont5%.R
#   matern:   .../RKCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/
#             10d_clx_cls.R (random-cluster process, variable n per rep)
# The UNCCD twins of these scripts generate IDENTICAL data (same code, same
# base seeds: 123 for gaussian/uniform, 1234 for matern), so one dataset per
# replicate serves all four methods.
#
# Method parameters at d=10, matched to the original scripts:
#   RKCCD: quant = 0.999 (table RK-test-simul_10d_999%.RData)
#   UNCCD: quant = 0.99  (table NN-test-simul_10d_99%.RData), method="descend"
#   IOS min.cls: 0 for gaussian/uniform (script default), 0.04 for matern
#
# Calibrated cutoffs: read at runtime from each method-family's own
# Threshold.R; d=10 corresponds to index 4 (vectors ordered for
# d = 2, 3, 5, 10, 20, 50, 100 -- verified against the if/else chain that
# maps d==10 -> Threshold_*[4] in every one of the six simulation scripts).
#
# Cutoff grid: 17 multipliers seq(0.5, 2, length.out = 17) x the calibrated
# value, per method, PLUS multiplier 1.0 (the seq grid has step 0.09375 and
# does not contain 1.0 exactly; the calibrated point itself is needed for
# the reviewer-facing figures and the reproduction cross-check) -> 18
# multipliers. Scores are computed ONCE per (setting, method, rep) and all
# 18 cutoffs are applied to the same score vector.
#
# NOTE: the 10-rep probe of 2026-07-10 ran with the original 17-point grid
# (68 rows per rep). After this grid change those probe reps no longer match
# EXPECTED_ROWS_PER_REP (72), so pending_reps() automatically drops and
# re-runs them at the next (production) launch. This is intentional.
#
# Checkpointing: rows appended per (setting, rep) chunk to
# results/wp3_synthetic_raw.csv; a (setting, rep) with a complete 4x17 block
# is skipped on rerun; partial blocks (crash mid-write) are dropped and
# re-run. Aggregation at the end writes
# results/wp3_sensitivity_synthetic.csv.
#
# CLI:  Rscript "revision_experiments/09_wp3_synthetic.R" <n_reps> <cores>
#       defaults: n_reps = 10, cores = 20

args <- commandArgs(trailingOnly = TRUE)
N_REPS <- if (length(args) >= 1) as.integer(args[1]) else 10L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 20L
stopifnot(!is.na(N_REPS), N_REPS >= 1, !is.na(CORES), CORES >= 1)

suppressMessages(source(here::here("revision_experiments/harness.R")))
# Uni.Gau_cls (Matern/Thomas random-cluster generator) + the `ratio` vector
# it is parameterized with, exactly as sourced by the original clx_cls scripts.
source(here::here("R/general_functions/Uni-Gau_cls.R"))
source(here::here("R/general_functions/ratio1.R"))

suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
})

REPO_ROOT <- here::here()
RAW_CSV <- file.path(REPO_ROOT, "revision_experiments/results/tr2/wp3_synthetic_raw.csv")
AGG_CSV <- file.path(REPO_ROOT, "revision_experiments/results/tr2/wp3_sensitivity_synthetic.csv")

D <- 10L
D_INDEX <- 4L  # d=10 -> index 4 (vectors ordered d = 2,3,5,10,20,50,100)
MULTIPLIERS <- sort(unique(c(seq(0.5, 2, length.out = 17), 1)))  # 18 (incl. calibrated 1.0)
OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
EXPECTED_ROWS_PER_REP <- length(OS_METHODS) * length(MULTIPLIERS)  # 72

# ---------------------------------------------------------------------------
# Calibrated cutoffs, read from each family/setting's own Threshold.R
# ---------------------------------------------------------------------------
read_thresholds <- function(rel_path) {
  e <- new.env()
  sys.source(here::here(rel_path), envir = e)
  list(OOS = get("Threshold_OOS", envir = e),
       IOS = get("Threshold_IOS", envir = e))
}

th_rk_gau <- read_thresholds("simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/Threshold.R")
th_rk_uni <- read_thresholds("simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Uniform/Threshold.R")
th_rk_mat <- read_thresholds("simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/Threshold.R")
th_nn_gau <- read_thresholds("simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/Gaussian/Threshold.R")
th_nn_uni <- read_thresholds("simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/Uniform/Threshold.R")
th_nn_mat <- read_thresholds("simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/Threshold.R")

CUTOFFS <- list(
  gaussian = c("RKCCD-OOS" = th_rk_gau$OOS[D_INDEX], "RKCCD-IOS" = th_rk_gau$IOS[D_INDEX],
               "UNCCD-OOS" = th_nn_gau$OOS[D_INDEX], "UNCCD-IOS" = th_nn_gau$IOS[D_INDEX]),
  uniform  = c("RKCCD-OOS" = th_rk_uni$OOS[D_INDEX], "RKCCD-IOS" = th_rk_uni$IOS[D_INDEX],
               "UNCCD-OOS" = th_nn_uni$OOS[D_INDEX], "UNCCD-IOS" = th_nn_uni$IOS[D_INDEX]),
  matern   = c("RKCCD-OOS" = th_rk_mat$OOS[D_INDEX], "RKCCD-IOS" = th_rk_mat$IOS[D_INDEX],
               "UNCCD-OOS" = th_nn_mat$OOS[D_INDEX], "UNCCD-IOS" = th_nn_mat$IOS[D_INDEX])
)

# Hard cross-check against the values read manually from the Threshold.R
# files (evidence recorded in the T7 report). Fails loudly if the files ever
# change underneath us.
stopifnot(
  identical(unname(CUTOFFS$gaussian), c(3.5, 6.5, 3.5, 6.5)),
  identical(unname(CUTOFFS$uniform),  c(4.0, 5.0, 3.0, 3.5)),
  identical(unname(CUTOFFS$matern),   c(4.0, 5.0, 3.0, 3.5))
)

# Quantile tables (master-side sanity check; workers load their own copies).
RK_TAB <- get_simul("RK", D)
NN_TAB <- get_simul("NN", D)
stopifnot(RK_TAB$quant_label == "999", NN_TAB$quant_label == "99")

# ---------------------------------------------------------------------------
# Data generators (verbatim from the first-cycle scripts, one dataset per
# call, seed argument replaces the original global set.seed)
# ---------------------------------------------------------------------------

# From RKCCD_OOS_IOS/Simulation/Gaussian/10d/10d_2cls_n500_cont5%.R
# (identical generator in the UNCCD twin).
gen_gaussian <- function(seed) {
  n = 500; d = 10; cont = 0.05
  cls_dis = 3   # the distances between each cluster center
  otl_dis = 2   # the minimal distances of outliers to cluster centers
  r_min = 0.7; r_max = 1.3

  mu1 = rep(3, d)
  mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = rbind(mu1, mu2)
  mu = apply(mu, 2, mean)

  n1 = round(n * (1 - cont) * 0.5)
  n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)

  noise_level = 0.01
  sigma = 1 / sqrt(qchisq(1 - noise_level, d))

  set.seed(seed)
  data1 = mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 = mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  i = 0
  outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    r1 = sqrt(sum((temp - mu1)^2))
    r2 = sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) {
      outlier = rbind(outlier, temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

# From RKCCD_OOS_IOS/Simulation/Uniform/10d/10d_2cls_n500_cont5%.R
# (identical generator in the UNCCD twin).
gen_uniform <- function(seed) {
  n = 500; d = 10; cont = 0.05
  cls_dis = 3
  otl_dis = 2
  r_min = 0.7; r_max = 1.3

  mu1 = rep(3, d)
  mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = rbind(mu1, mu2)
  mu = apply(mu, 2, mean)

  n1 = round(n * (1 - cont) * 0.5)
  n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)

  set.seed(seed)
  data1 = rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = T)
  data2 = rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = T)
  i = 0
  outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    r1 = sqrt(sum((temp - mu1)^2))
    r2 = sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) {
      outlier = rbind(outlier, temp)
      i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

# From RKCCD_OOS_IOS/Simulation/Complex_Clusters/Matern/10d_clx_cls.R
# (identical generator + parameters in the UNCCD twin; both seed the
# generation phase with 1234). Uni.Gau_cls returns, in order:
# [[1]] Matérn_parents, [[2]] Matérn_children, [[3]] Thomas_parents,
# [[4]] Thomas_children, [[5]] noise, [[6]] Outlier, [[7]] num.
# Positional access avoids the non-ASCII "Matérn" name.
gen_matern <- function(seed) {
  d = 10
  # simulation settings (verbatim)
  kappa1 = 6
  mu1 = ratio[4]
  expand1 = 0
  r = 0.1
  kappa2 = 0
  scale = 0.005
  mu2 = ratio[4]
  expand2 = 0
  slen = 1
  kappa_O = 20

  set.seed(seed)
  data_simu = Uni.Gau_cls(d, kappa1, r, mu1, expand1, kappa2, scale, mu2, expand2, slen, kappa_O)
  cls1 = data_simu[[2]]      # Matérn_children
  cls2 = data_simu[[4]]      # Thomas_children (empty here: kappa2 = 0)
  cls3 = data_simu$noise
  outlier = data_simu$Outlier
  num = data_simu$num        # c(n_regular, n0)
  list(X = rbind(cls1, cls2, cls3, outlier), n = sum(num), n0 = num[2])
}

SETTINGS <- list(
  gaussian = list(gen = "gen_gaussian", base_seed = 123L,  min_cls = 0),
  uniform  = list(gen = "gen_uniform",  base_seed = 123L,  min_cls = 0),
  matern   = list(gen = "gen_matern",   base_seed = 1234L, min_cls = 0.04)
)

# ---------------------------------------------------------------------------
# Per-replicate worker: generate data, score each method ONCE, apply all
# 17 cutoffs to the same scores.
# ---------------------------------------------------------------------------
score_one_rep <- function(setting, rep_id) {
  cfg <- SETTINGS[[setting]]
  seed <- cfg$base_seed + rep_id
  dat <- do.call(cfg$gen, list(seed = seed))
  X <- dat$X; n <- dat$n; n0 <- dat$n0
  stopifnot(nrow(X) == n)
  Y <- c(rep(1, n - n0), rep(0, n0))  # outliers are the trailing n0 rows

  # UNCCD at d=10 uses method="descend" in all three original 10d scripts.
  score_calls <- list(
    "RKCCD-OOS" = function() RKCCD_OOS(datax = X, simul = RK_TAB$simul, d = D, quant = RK_TAB$quant),
    "RKCCD-IOS" = function() RKCCD_IOS(datax = X, simul = RK_TAB$simul, d = D, quant = RK_TAB$quant, min.cls = cfg$min_cls),
    "UNCCD-OOS" = function() NNCCD_OOS(datax = X, simul = NN_TAB$simul, method = "descend", d = D),
    "UNCCD-IOS" = function() NNCCD_IOS(datax = X, simul = NN_TAB$simul, method = "descend", d = D, min.cls = cfg$min_cls)
  )

  rows <- vector("list", length(OS_METHODS))
  for (k in seq_along(OS_METHODS)) {
    m <- OS_METHODS[k]
    t0 <- Sys.time()
    score <- score_calls[[m]]()
    t_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (length(score) != n || any(is.na(score))) {
      stop(sprintf("[%s | %s | rep %d] bad score vector (len %d vs n %d, NA %s)",
                   setting, m, rep_id, length(score), n, any(is.na(score))))
    }
    calib <- CUTOFFS[[setting]][[m]]
    block <- lapply(MULTIPLIERS, function(mult) {
      v <- evaluate(Y, score, mult * calib)
      data.frame(setting = setting, method = m, rep = rep_id, seed = seed,
                 n = n, n0 = n0, multiplier = mult, cutoff = mult * calib,
                 TPR = unname(v["TPR"]), TNR = unname(v["TNR"]),
                 BA = unname(v["BA"]), F2 = unname(v["F2"]),
                 t_score = round(t_score, 3),
                 stringsAsFactors = FALSE)
    })
    rows[[k]] <- do.call(rbind, block)
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Checkpointing helpers
# ---------------------------------------------------------------------------
pending_reps <- function(setting, n_reps) {
  if (!file.exists(RAW_CSV)) return(seq_len(n_reps))
  df <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
  sub <- df[df$setting == setting, , drop = FALSE]
  if (nrow(sub) == 0) return(seq_len(n_reps))
  cnt <- table(sub$rep)
  done <- as.integer(names(cnt)[cnt == EXPECTED_ROWS_PER_REP])
  partial <- as.integer(names(cnt)[cnt != EXPECTED_ROWS_PER_REP])
  if (length(partial) > 0) {
    cat(sprintf("[%s] dropping %d partially-written rep(s) from %s: %s\n",
                setting, length(partial), basename(RAW_CSV),
                paste(partial, collapse = ",")))
    keep <- !(df$setting == setting & df$rep %in% partial)
    write.csv(df[keep, , drop = FALSE], RAW_CSV, row.names = FALSE)
  }
  setdiff(seq_len(n_reps), done)
}

append_rows <- function(csv_path, df) {
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(csv_path)) {
    write.csv(df, csv_path, row.names = FALSE)
  } else {
    write.table(df, csv_path, sep = ",", col.names = FALSE, row.names = FALSE,
                append = TRUE, qmethod = "double")
  }
}

# ---------------------------------------------------------------------------
# Main: parallel sweep over replicates
# ---------------------------------------------------------------------------
cat(sprintf("09_wp3_synthetic.R: n_reps = %d, cores = %d\n", N_REPS, CORES))
cat("Calibrated cutoffs (d=10):\n")
for (s in names(CUTOFFS)) {
  cat(sprintf("  %-8s %s\n", s,
              paste(sprintf("%s=%.1f", names(CUTOFFS[[s]]), CUTOFFS[[s]]), collapse = "  ")))
}

t_start <- Sys.time()

cl <- makeCluster(CORES)
registerDoParallel(cl)
clusterExport(cl, "REPO_ROOT")
invisible(clusterEvalQ(cl, {
  setwd(REPO_ROOT)
  suppressMessages(source(file.path(REPO_ROOT, "revision_experiments/harness.R")))
  source(file.path(REPO_ROOT, "R/general_functions/Uni-Gau_cls.R"))
  source(file.path(REPO_ROOT, "R/general_functions/ratio1.R"))
  RK_TAB <- get_simul("RK", 10)
  NN_TAB <- get_simul("NN", 10)
  TRUE
}))
clusterExport(cl, c("D", "MULTIPLIERS", "OS_METHODS", "CUTOFFS", "SETTINGS",
                    "gen_gaussian", "gen_uniform", "gen_matern", "score_one_rep"))

run_chunk <- function(st, chunk) {
  foreach(r = chunk, .combine = rbind,
          .packages = c("MASS", "cluster", "igraph")) %dopar% score_one_rep(st, r)
}

n_run <- 0L
for (setting in names(SETTINGS)) {
  pend <- pending_reps(setting, N_REPS)
  if (length(pend) == 0) {
    cat(sprintf("[%s] all %d reps already recorded; skipping.\n", setting, N_REPS))
    next
  }
  cat(sprintf("[%s] %d/%d reps pending.\n", setting, length(pend), N_REPS))
  chunks <- split(pend, ceiling(seq_along(pend) / CORES))
  for (ci in seq_along(chunks)) {
    tc0 <- Sys.time()
    res <- run_chunk(setting, chunks[[ci]])
    append_rows(RAW_CSV, res)
    n_run <- n_run + length(chunks[[ci]])
    cat(sprintf("[%s] chunk %d/%d (%d reps) done in %.1f s (elapsed %.1f min)\n",
                setting, ci, length(chunks), length(chunks[[ci]]),
                as.numeric(difftime(Sys.time(), tc0, units = "secs")),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

stopCluster(cl)

wall_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("Sweep finished: %d newly-run reps in %.1f min (%.1f reps/hour on %d cores)\n",
            n_run, wall_min,
            if (wall_min > 0) n_run / wall_min * 60 else NA, CORES))

# ---------------------------------------------------------------------------
# Aggregation -> wp3_sensitivity_synthetic.csv
# ---------------------------------------------------------------------------
raw <- read.csv(RAW_CSV, stringsAsFactors = FALSE)

agg <- do.call(rbind, lapply(split(raw, list(raw$setting, raw$method, raw$multiplier), drop = TRUE), function(g) {
  data.frame(
    setting = g$setting[1], method = g$method[1],
    multiplier = g$multiplier[1], cutoff = g$cutoff[1],
    TPR_mean = mean(g$TPR), TNR_mean = mean(g$TNR),
    BA_mean = mean(g$BA), F2_mean = mean(g$F2),
    TPR_sd = sd(g$TPR), TNR_sd = sd(g$TNR),
    BA_sd = sd(g$BA), F2_sd = sd(g$F2),
    n_reps = nrow(g),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$setting, agg$method, agg$multiplier), ]
rownames(agg) <- NULL
write.csv(agg, AGG_CSV, row.names = FALSE)
cat(sprintf("Aggregated %d raw rows into %d (setting, method, multiplier) cells -> %s\n",
            nrow(raw), nrow(agg), AGG_CSV))

# Quick sanity summary at the calibrated cutoff (multiplier = 1)
cat("\nBA_mean / F2_mean at multiplier = 1.0 (calibrated cutoff):\n")
at1 <- agg[abs(agg$multiplier - 1) < 1e-9, ]
for (i in seq_len(nrow(at1))) {
  cat(sprintf("  %-8s %-10s BA=%.3f F2=%.3f (n_reps=%d)\n",
              at1$setting[i], at1$method[i], at1$BA_mean[i], at1$F2_mean[i], at1$n_reps[i]))
}

stopifnot(!anyNA(raw[, c("TPR", "TNR", "BA", "F2")]),
          all(raw$TPR >= 0 & raw$TPR <= 1), all(raw$TNR >= 0 & raw$TNR <= 1),
          all(raw$BA >= 0 & raw$BA <= 1), all(raw$F2 >= 0 & raw$F2 <= 1))
cat("\nAll metrics in [0,1], no NA. 09_wp3_synthetic.R done.\n")
