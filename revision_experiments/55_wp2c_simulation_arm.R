#!/usr/bin/env Rscript
# revision_experiments/55_wp2c_simulation_arm.R
#
# WP2(c), Step 3: the SIMULATION arm -- the actual deliverable.
#
# Main text L605: SU-MCCD is evaluated "under the same simulation settings ...
# with S_min set to half the contamination level". That is the oracle rule
# R3.2/R5.2 object to: it reads the true contamination, i.e. the labels. This
# script re-runs the paper's synthetic settings under that oracle rule and
# under each label-free candidate, and measures the paired difference.
#
# --------------------------------------------------------------------------
# WHAT IS RE-RUN, AND HOW IT MAPS TO THE MANUSCRIPT
# --------------------------------------------------------------------------
# Generators: gen_uniform() and gen_gaussian() are lifted from
# revision_experiments/09_wp3_synthetic.R lines 124-195, which in turn copies
# them verbatim from the first-cycle simulation drivers
#   simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Uniform/10d/10d_2cls_n500_cont5%.R
#   simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/10d/10d_2cls_n500_cont5%.R
# The ONLY change here is that n, d and the contamination level are arguments
# instead of literals -- every other constant (2 clusters, inter-centre
# distance cls_dis = 3, outlier stand-off otl_dis = 2, radius jitter
# r_min/r_max = 0.7/1.3, Gaussian noise_level = 0.01, outlier radius 5) is
# byte-identical to the drivers. The manuscript's Section 5 sweeps
# d in {2,3,5,10,20,50,100} and n in {50,...,1000}; this arm fixes n = 500
# (the drivers' own n) and takes d in {5, 10} -- one below and one at the
# alpha-schedule breakpoint (RK 1% -> 0.1%, main text L880) -- and four
# contamination levels, because the oracle rule's entire content is that it
# MOVES with contamination and a single level cannot test it.
#
#   generator  in {uniform, gaussian}      2
#   d          in {5, 10}                  2
#   contam     in {0.02, 0.05, 0.10, 0.20} 4    -> 16 settings
#   n = 500, 2 clusters, 100 replicates, seeds base + rep (base 123 for both
#   generators, as in the drivers and in 09_wp3_synthetic.R).
#
# Methods: SU-MCCD and SUN-MCCD only. U-MCCD and UN-MCCD take no min.cls
# argument, so no S_min rule can change them.
#
# Metrics via harness.R's evaluate(Y, score, 0.5); Y == 1 regular, Y == 0
# outlier; larger score = more outlying. Never bypassed.
#
# --------------------------------------------------------------------------
# RULES COMPARED
# --------------------------------------------------------------------------
#   oracle_half_contamination  S_min = 0.5 * contamination        LABEL-USING
#   min_cls_zero               S_min = 0 (the detectors' shipped default)
#   fixed_0.0625               the manuscript's real-data value
#   fixed_0.0500               the other end of the [0.05, 0.0625] bracket
#                              that reproduces the published SUN-MCCD numbers
#   fixed_0.0300               the constant with the most plateau hits on real
#                              data (54_wp2c_smin_rule.R candidates)
#   adaptive_one_step_rho0.5   S_min = 0.5 * median(cluster size at
#                              min.cls = 0) / n -- circularity resolved by
#                              taking the partition once at the neutral seed
#                              min.cls = 0 and NOT iterating (iterating drives
#                              it to a degenerate fixed point; measured on
#                              real data in wp2c_adaptive_real.csv)
#
# COMPUTE SAVING, and why it is exact rather than an approximation: min.cls
# reaches the detectors only through k = round(min.cls * n)
# (R/ccds/RK_CCD_New.R:113, R/ccds/UN_CCD.R:102 -- see the header of
# 54_wp2c_smin_rule.R). Two rules that select the same k on a given replicate
# give bit-identical output, so the detector is run ONCE PER DISTINCT k and
# the result is shared by every rule that selected it.
#
# --------------------------------------------------------------------------
# OUTPUTS (results/tr1/)
#   wp2c_simulation_arm.csv      one row per (setting, rep, method, rule)
#   wp2c_simulation_summary.csv  one row per (setting, method, rule), with the
#                                paired label-free-minus-oracle deltas
#
# USAGE
#   Rscript 55_wp2c_simulation_arm.R probe                 cost probe, 2 reps
#   Rscript 55_wp2c_simulation_arm.R run [settings] [reps] [budget] [cores]
#       settings: comma list of setting_ids, or ALL
#       reps:     target replicate count (default 100)
#       budget:   stop starting new chunks after this many seconds (default 480)
#       cores:    default 6 (three other jobs share this 24-core box)
#   Rscript 55_wp2c_simulation_arm.R summarize
#   Rscript 55_wp2c_simulation_arm.R status

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments", "harness.R")))
suppressMessages(source(here::here("revision_experiments", "wp0_mccd_methods.R")))
suppressPackageStartupMessages({ library(parallel); library(doParallel); library(foreach) })

REPO   <- here::here()
RESDIR <- file.path(REPO, "revision_experiments/results/tr1")
RAW    <- file.path(RESDIR, "wp2c_simulation_arm.csv")
SUMM   <- file.path(RESDIR, "wp2c_simulation_summary.csv")

BASE_SEED <- 123L
METHODS   <- c("SU-MCCD", "SUN-MCCD")
RULES     <- c("oracle_half_contamination", "min_cls_zero",
               "fixed_0.0625", "fixed_0.0500", "fixed_0.0300",
               "adaptive_one_step_rho0.5")
EXPECTED_ROWS_PER_REP <- length(METHODS) * length(RULES)   # 12

# Measured per-detector cost (this box, one core, n/d as shown, min.cls=0.05):
#   n=200 d=5  SU 9.01 s  SUN 1.17 s      n=500 d=5  SU 24.92 s SUN 2.89 s
#   n=200 d=10 SU 4.05 s  SUN 0.88 s      n=500 d=10 SU  7.23 s SUN 4.56 s
#   n=200 d=20 SU 2.75 s  SUN 0.34 s
# Six rules resolve to ~5-6 distinct k, so a replicate costs ~6x the pair
# cost. n = 200 is one of the manuscript's own sample sizes
# (n in {50,100,200,500,1000}, Section 5) and is what makes a 100-replicate,
# 4-contamination-level, 2-dimension, 2-generator design affordable at 6
# cores; the drivers' exact n = 500, d = 10, 5% configuration is retained as
# a separate fidelity anchor so the reduction can be checked rather than
# assumed.
SETTINGS <- rbind(
  do.call(rbind, lapply(c("uniform", "gaussian"), function(gen)
    do.call(rbind, lapply(c(5L, 10L, 20L), function(d)
      do.call(rbind, lapply(c(0.02, 0.05, 0.10, 0.20), function(ct)
        data.frame(generator = gen, n_nominal = 200L, d = d, contam = ct,
                   tier = if (d == 5L) 3L else if (d == 10L) 1L else 2L,
                   stringsAsFactors = FALSE))))))),
  data.frame(generator = "uniform",  n_nominal = 500L, d = 10L, contam = 0.05, tier = 1L,
             stringsAsFactors = FALSE),
  data.frame(generator = "gaussian", n_nominal = 500L, d = 10L, contam = 0.05, tier = 2L,
             stringsAsFactors = FALSE))
SETTINGS$setting_id <- sprintf("%s_d%d_n%d_c%02d", SETTINGS$generator, SETTINGS$d,
                               SETTINGS$n_nominal, round(SETTINGS$contam * 100))

# ---------------------------------------------------------------------------
# Generators (09_wp3_synthetic.R:124-195, n/d/cont promoted to arguments)
# ---------------------------------------------------------------------------
gen_uniform <- function(seed, n, d, cont) {
  cls_dis = 3; otl_dis = 2; r_min = 0.7; r_max = 1.3
  mu1 = rep(3, d); mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = apply(rbind(mu1, mu2), 2, mean)
  n1 = round(n * (1 - cont) * 0.5); n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)
  set.seed(seed)
  data1 = rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 = rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i = 0; outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier = rbind(outlier, temp); i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

gen_gaussian <- function(seed, n, d, cont) {
  cls_dis = 3; otl_dis = 2; r_min = 0.7; r_max = 1.3
  mu1 = rep(3, d); mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = apply(rbind(mu1, mu2), 2, mean)
  n1 = round(n * (1 - cont) * 0.5); n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)
  noise_level = 0.01
  sigma = 1 / sqrt(qchisq(1 - noise_level, d))
  set.seed(seed)
  data1 = mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 = mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  i = 0; outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier = rbind(outlier, temp); i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

# ---------------------------------------------------------------------------
# One replicate: all methods, all rules, one detector call per distinct k
# ---------------------------------------------------------------------------
run_one_rep <- function(setting_id, generator, n_nominal, d, contam, rep_id) {
  seed <- BASE_SEED + rep_id
  dat  <- if (generator == "uniform") gen_uniform(seed, n_nominal, d, contam)
          else                        gen_gaussian(seed, n_nominal, d, contam)
  X <- dat$X; n <- dat$n; n0 <- dat$n0
  Y <- c(rep(1, n - n0), rep(0, n0))   # outliers are the trailing n0 rows

  out <- list()
  for (meth in METHODS) {
    cache <- list()   # key = as.character(k) -> result list
    call_at <- function(s_min) {
      k <- round(s_min * n); key <- as.character(k)
      if (!is.null(cache[[key]])) return(cache[[key]])
      t0 <- Sys.time()
      r <- tryCatch({
        res <- METHOD_REGISTRY[[meth]](X = X, d = d, Y = Y, min.cls = s_min)
        stopifnot(length(res$score) == n, !anyNA(res$score))
        m  <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[[meth]])
        cl <- res$cluster; cl <- cl[!is.na(cl)]
        tab <- table(cl); cs <- as.integer(tab)
        list(TPR = unname(m[["TPR"]]), TNR = unname(m[["TNR"]]),
             BA = unname(m[["BA"]]), F2 = unname(m[["F2"]]),
             n_flagged = sum(res$score == 1), n_clusters = length(tab),
             mcs = if (length(cs)) median(cs) else NA_real_, status = "ok")
      }, error = function(e)
        list(TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
             n_flagged = NA_integer_, n_clusters = NA_integer_, mcs = NA_real_,
             status = paste0("error: ", substr(conditionMessage(e), 1, 200))))
      r$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      r$k <- k
      cache[[key]] <<- r
      r
    }

    # seed run first: min.cls = 0 both is a rule in its own right and the
    # input the one-step adaptive rule needs.
    r0 <- call_at(0)
    s_adapt <- if (is.na(r0$mcs)) NA_real_ else 0.5 * r0$mcs / n

    rule_smin <- c(
      oracle_half_contamination = 0.5 * contam,
      min_cls_zero              = 0,
      fixed_0.0625              = 0.0625,
      fixed_0.0500              = 0.0500,
      fixed_0.0300              = 0.0300,
      adaptive_one_step_rho0.5  = s_adapt)

    for (rl in RULES) {
      s <- unname(rule_smin[[rl]])
      r <- if (is.na(s)) list(TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
                              n_flagged = NA_integer_, n_clusters = NA_integer_, mcs = NA_real_,
                              status = "skipped: seed run gave no cluster sizes",
                              elapsed_sec = NA_real_, k = NA_integer_)
           else call_at(s)
      out[[length(out) + 1L]] <- data.frame(
        setting_id = setting_id, generator = generator, d = d, contam = contam,
        n_nominal = n_nominal, n = n, n0 = n0, rep = rep_id, seed = seed,
        method = meth, rule = rl,
        label_free = rl != "oracle_half_contamination",
        s_min = s, k = r$k,
        TPR = r$TPR, TNR = r$TNR, BA = r$BA, F2 = r$F2,
        n_flagged = r$n_flagged, n_clusters = r$n_clusters,
        median_cluster_size = r$mcs,
        seed_median_cluster_size = r0$mcs, seed_n_clusters = r0$n_clusters,
        elapsed_sec = r$elapsed_sec, status = r$status,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

append_rows <- function(path, df) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) write.csv(df, path, row.names = FALSE)
  else write.table(df, path, sep = ",", col.names = FALSE, row.names = FALSE,
                   append = TRUE, qmethod = "double")
}

pending_reps <- function(setting_id, n_reps) {
  if (!file.exists(RAW)) return(seq_len(n_reps))
  df <- tryCatch(read.csv(RAW, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(seq_len(n_reps))
  sub <- df[df$setting_id == setting_id, , drop = FALSE]
  if (nrow(sub) == 0) return(seq_len(n_reps))
  cnt  <- table(sub$rep)
  done <- as.integer(names(cnt)[cnt == EXPECTED_ROWS_PER_REP])
  setdiff(seq_len(n_reps), done)
}

# ---------------------------------------------------------------------------
# modes
# ---------------------------------------------------------------------------
do_probe <- function() {
  for (i in seq_len(nrow(SETTINGS))) {
    s <- SETTINGS[i, ]
    t0 <- Sys.time()
    r <- run_one_rep(s$setting_id, s$generator, s$n_nominal, s$d, s$contam, 9001L)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("%-20s n=%d n0=%d  distinct_k=%s  1 rep = %.1f s  -> 100 reps on 6 cores ~ %.1f min\n",
                s$setting_id, r$n[1], r$n0[1],
                paste(sort(unique(r$k[!is.na(r$k)])), collapse = ","),
                el, el * 100 / 6 / 60))
    flush.console()
  }
}

resolve_settings <- function(sel) {
  if (is.null(sel)) return(SETTINGS[order(SETTINGS$tier), ])
  if (length(sel) == 1 && grepl("^TIER[0-9]$", toupper(sel)))
    return(SETTINGS[SETTINGS$tier == as.integer(sub("TIER", "", toupper(sel))), ])
  SETTINGS[match(sel, SETTINGS$setting_id), ]   # keep the caller's order
}

do_run <- function(sel, n_reps, budget, cores) {
  st <- resolve_settings(sel)
  st <- st[!is.na(st$setting_id), ]
  stopifnot(nrow(st) > 0)
  cl <- makeCluster(cores); registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, "REPO", envir = environment())
  invisible(clusterEvalQ(cl, {
    setwd(REPO)
    suppressMessages(source(file.path(REPO, "revision_experiments/harness.R")))
    suppressMessages(source(file.path(REPO, "revision_experiments/wp0_mccd_methods.R")))
    TRUE
  }))
  clusterExport(cl, c("BASE_SEED", "METHODS", "RULES",
                      "gen_uniform", "gen_gaussian", "run_one_rep"),
                envir = environment())

  t_start <- Sys.time()
  for (i in seq_len(nrow(st))) {
    s <- st[i, ]
    pend <- pending_reps(s$setting_id, n_reps)
    if (!length(pend)) { cat(sprintf("[%s] complete (%d reps)\n", s$setting_id, n_reps)); next }
    cat(sprintf("[%s] %d/%d reps pending\n", s$setting_id, length(pend), n_reps)); flush.console()
    chunks <- split(pend, ceiling(seq_along(pend) / cores))
    for (ci in seq_along(chunks)) {
      if (as.numeric(difftime(Sys.time(), t_start, units = "secs")) > budget) {
        cat("[budget exhausted -- rerun to continue]\n"); return(invisible(NULL))
      }
      tc <- Sys.time()
      res <- foreach(rr = chunks[[ci]], .combine = rbind,
                     .packages = c("MASS", "cluster", "igraph"),
                     .errorhandling = "remove") %dopar%
        run_one_rep(s$setting_id, s$generator, s$n_nominal, s$d, s$contam, rr)
      if (!is.null(res) && nrow(res)) append_rows(RAW, res)   # on disk per chunk
      cat(sprintf("  [%s] chunk %d/%d (%d reps) %.1f s | total %.1f min\n",
                  s$setting_id, ci, length(chunks), length(chunks[[ci]]),
                  as.numeric(difftime(Sys.time(), tc, units = "secs")),
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
      flush.console()
    }
  }
}

do_status <- function(n_reps = 100L) {
  if (!file.exists(RAW)) { cat("no raw file yet\n"); return(invisible(NULL)) }
  df <- read.csv(RAW, stringsAsFactors = FALSE)
  tb <- sapply(SETTINGS$setting_id, function(sid) {
    sub <- df[df$setting_id == sid, ]
    if (!nrow(sub)) return(0L)
    cnt <- table(sub$rep); sum(cnt == EXPECTED_ROWS_PER_REP)
  })
  print(data.frame(setting_id = SETTINGS$setting_id, complete_reps = as.integer(tb),
                   row.names = NULL))
  cat(sprintf("total rows %d; settings complete at %d reps: %d/%d\n",
              nrow(df), n_reps, sum(tb >= n_reps), nrow(SETTINGS)))
  bad <- df[df$status != "ok", ]
  cat(sprintf("non-ok rows: %d\n", nrow(bad)))
  if (nrow(bad)) print(table(bad$setting_id, bad$status))
}

do_summarize <- function() {
  df <- read.csv(RAW, stringsAsFactors = FALSE)
  df <- df[df$status == "ok", ]
  rows <- list()
  for (sid in unique(df$setting_id)) for (meth in METHODS) {
    blk <- df[df$setting_id == sid & df$method == meth, ]
    if (!nrow(blk)) next
    orc <- blk[blk$rule == "oracle_half_contamination", ]
    for (rl in RULES) {
      sub <- blk[blk$rule == rl, ]
      if (!nrow(sub)) next
      # pair on rep
      mm <- merge(sub, orc[, c("rep", "BA", "F2", "TPR", "TNR", "n_clusters")],
                  by = "rep", suffixes = c("", "_orc"))
      rows[[length(rows) + 1L]] <- data.frame(
        setting_id = sid, generator = sub$generator[1], d = sub$d[1],
        contam = sub$contam[1], method = meth, rule = rl,
        label_free = sub$label_free[1],
        n_reps = nrow(sub), n_paired = nrow(mm),
        s_min_mean = mean(sub$s_min), k_mean = mean(sub$k),
        BA_mean = mean(sub$BA),  BA_sd = sd(sub$BA),
        F2_mean = mean(sub$F2),  F2_sd = sd(sub$F2),
        TPR_mean = mean(sub$TPR), TPR_sd = sd(sub$TPR),
        TNR_mean = mean(sub$TNR), TNR_sd = sd(sub$TNR),
        n_clusters_mean = mean(sub$n_clusters),
        n_clusters_median = median(sub$n_clusters),
        n_clusters_min = min(sub$n_clusters), n_clusters_max = max(sub$n_clusters),
        frac_single_cluster = mean(sub$n_clusters == 1),
        # paired label-free minus oracle
        dBA_mean = mean(mm$BA - mm$BA_orc), dBA_sd = sd(mm$BA - mm$BA_orc),
        dBA_se   = sd(mm$BA - mm$BA_orc) / sqrt(nrow(mm)),
        dF2_mean = mean(mm$F2 - mm$F2_orc), dF2_sd = sd(mm$F2 - mm$F2_orc),
        dF2_se   = sd(mm$F2 - mm$F2_orc) / sqrt(nrow(mm)),
        dTPR_mean = mean(mm$TPR - mm$TPR_orc), dTPR_sd = sd(mm$TPR - mm$TPR_orc),
        dTNR_mean = mean(mm$TNR - mm$TNR_orc), dTNR_sd = sd(mm$TNR - mm$TNR_orc),
        frac_identical_k = mean(sub$k == orc$k[match(sub$rep, orc$rep)], na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  S <- S[order(S$generator, S$d, S$contam, S$method, S$rule), ]
  write.csv(S, SUMM, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMM, nrow(S)))

  cat("\n=== paired label-free minus oracle, BA (mean +/- sd over reps) ===\n")
  lf <- S[S$label_free, ]
  pr <- lf[, c("setting_id", "method", "rule", "n_paired", "BA_mean", "dBA_mean", "dBA_sd", "dF2_mean", "dF2_sd", "frac_single_cluster")]
  pr[, 5:10] <- round(pr[, 5:10], 4)
  print(pr, row.names = FALSE)

  cat("\n=== pooled over all settings, per rule ===\n")
  key <- paste(lf$rule)
  po <- data.frame(
    rule = tapply(lf$rule, key, function(x) x[1]),
    n_cells = tapply(lf$dBA_mean, key, length),
    mean_dBA = round(tapply(lf$dBA_mean, key, mean), 4),
    worst_dBA = round(tapply(lf$dBA_mean, key, min), 4),
    best_dBA = round(tapply(lf$dBA_mean, key, max), 4),
    mean_dF2 = round(tapply(lf$dF2_mean, key, mean), 4),
    worst_dF2 = round(tapply(lf$dF2_mean, key, min), 4),
    mean_BA = round(tapply(lf$BA_mean, key, mean), 4),
    mean_frac_single_cluster = round(tapply(lf$frac_single_cluster, key, mean), 3),
    stringsAsFactors = FALSE)
  print(po[order(-po$mean_dBA), ], row.names = FALSE)

  orc <- S[!S$label_free, ]
  cat(sprintf("\noracle: mean BA over %d cells = %.4f, mean F2 = %.4f, mean frac single cluster = %.3f\n",
              nrow(orc), mean(orc$BA_mean), mean(orc$F2_mean), mean(orc$frac_single_cluster)))
  invisible(S)
}

# ---------------------------------------------------------------------------
# report: the four things the response letter needs, printed from the CSVs.
#
# NOTE ON `k`. k is round(s_min * n) -- the integer cluster-size floor the
# detectors actually apply (R/ccds/RK_CCD_New.R:113, R/ccds/UN_CCD.R:102). It
# is NOT a cluster count. 16 of 10608 rows fail a naive re-derivation
# k == round(s_min*n) after the CSV round trip; all 16 are
# adaptive_one_step_rho0.5 rows whose product lands exactly on a .5 boundary,
# where R's half-to-even rounding plus a 1-ULP shift in the reread s_min flips
# the result. The recorded k is the one the detector used.
#
# NOTE ON THE TRUE CLUSTER COUNT. There is no ground-truth-k column in this
# file, and none is needed: both generators build exactly two clusters
# (data1, data2) plus outliers in every setting, so true k = 2 throughout.
# n_clusters is the ESTIMATED macro-cluster count. k-hat accuracy is therefore
# computable from this file as mean(n_clusters == 2) -- reported below.
do_report <- function() {
  df <- read.csv(RAW, stringsAsFactors = FALSE)
  df <- df[df$status == "ok", ]
  TRUE_K <- 2L

  cat("\n############ 1. PAIRED label-free MINUS oracle, per setting x method ############\n")
  cat("dBA / dF2 are per-replicate differences (rule minus oracle_half_contamination),\n")
  cat("averaged over the paired replicates; sd is across replicates.\n\n")
  for (sid in sort(unique(df$setting_id))) {
    blk <- df[df$setting_id == sid, ]
    cat(sprintf("--- %s   (n=%d, n0=%d, contamination=%.2f) ---\n",
                sid, blk$n[1], blk$n0[1], blk$contam[1]))
    for (meth in METHODS) {
      b <- blk[blk$method == meth, ]
      orc <- b[b$rule == "oracle_half_contamination", ]
      cat(sprintf("  %-9s oracle: S_min=%.4f k=%d  BA=%.4f(sd %.4f) F2=%.4f(sd %.4f) TPR=%.3f TNR=%.3f  ncls med=%d, ==1 in %.0f%%, ==2 in %.0f%%\n",
                  meth, orc$s_min[1], orc$k[1],
                  mean(orc$BA), sd(orc$BA), mean(orc$F2), sd(orc$F2),
                  mean(orc$TPR), mean(orc$TNR),
                  median(orc$n_clusters), 100 * mean(orc$n_clusters == 1),
                  100 * mean(orc$n_clusters == TRUE_K)))
      for (rl in setdiff(RULES, "oracle_half_contamination")) {
        s <- b[b$rule == rl, ]
        if (!nrow(s)) next
        mm <- merge(s, orc[, c("rep", "BA", "F2")], by = "rep", suffixes = c("", "_o"))
        dB <- mm$BA - mm$BA_o; dF <- mm$F2 - mm$F2_o
        cat(sprintf("    %-26s S_min=%-7.4f dBA=%+.4f (sd %.4f, se %.4f) dF2=%+.4f (sd %.4f) | BA=%.4f F2=%.4f | ncls==1 %.0f%%, ==2 %.0f%% | reps=%d\n",
                    rl, mean(s$s_min), mean(dB), sd(dB), sd(dB) / sqrt(length(dB)),
                    mean(dF), sd(dF), mean(s$BA), mean(s$F2),
                    100 * mean(s$n_clusters == 1), 100 * mean(s$n_clusters == TRUE_K),
                    nrow(mm)))
      }
    }
    cat("\n")
  }

  cat("\n############ 2. POOLED over all settings x methods, per rule ############\n")
  lf <- df[df$label_free, ]
  cells <- unique(df[, c("setting_id", "method")])
  per_cell <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    b <- df[df$setting_id == cells$setting_id[i] & df$method == cells$method[i], ]
    orc <- b[b$rule == "oracle_half_contamination", ]
    do.call(rbind, lapply(setdiff(RULES, "oracle_half_contamination"), function(rl) {
      s <- b[b$rule == rl, ]; if (!nrow(s)) return(NULL)
      mm <- merge(s, orc[, c("rep", "BA", "F2")], by = "rep", suffixes = c("", "_o"))
      data.frame(setting_id = cells$setting_id[i], method = cells$method[i], rule = rl,
                 dBA = mean(mm$BA - mm$BA_o), dF2 = mean(mm$F2 - mm$F2_o),
                 BA = mean(s$BA), F2 = mean(s$F2),
                 frac1 = mean(s$n_clusters == 1), frac2 = mean(s$n_clusters == TRUE_K),
                 stringsAsFactors = FALSE)
    }))
  }))
  po <- do.call(rbind, lapply(split(per_cell, per_cell$rule), function(g)
    data.frame(rule = g$rule[1], n_cells = nrow(g),
               mean_dBA = mean(g$dBA), worst_dBA = min(g$dBA), best_dBA = max(g$dBA),
               cells_worse = sum(g$dBA < -1e-9), cells_better = sum(g$dBA > 1e-9),
               mean_dF2 = mean(g$dF2), worst_dF2 = min(g$dF2),
               mean_BA = mean(g$BA), mean_F2 = mean(g$F2),
               frac_single_cluster = mean(g$frac1), khat_correct = mean(g$frac2),
               stringsAsFactors = FALSE)))
  po <- po[order(-po$mean_dBA), ]
  print(format(po, digits = 4), row.names = FALSE)
  orc_all <- df[!df$label_free, ]
  cat(sprintf("\nORACLE reference: BA=%.4f F2=%.4f | n_clusters==1 in %.1f%% of runs, ==2 (correct) in %.1f%%\n",
              mean(orc_all$BA), mean(orc_all$F2),
              100 * mean(orc_all$n_clusters == 1), 100 * mean(orc_all$n_clusters == TRUE_K)))
  cat("worst cell per rule:\n")
  for (rl in unique(per_cell$rule)) {
    g <- per_cell[per_cell$rule == rl, ]; w <- g[which.min(g$dBA), ]
    cat(sprintf("  %-26s %s / %-9s dBA=%+.4f dF2=%+.4f (BA %.4f)\n",
                rl, w$setting_id, w$method, w$dBA, w$dF2, w$BA))
  }

  cat("\n############ 3. CLUSTER-COUNT COLLAPSE, by rule (true k = 2 everywhere) ############\n")
  cc <- do.call(rbind, lapply(split(df, df$rule), function(g)
    data.frame(rule = g$rule[1], n_runs = nrow(g),
               mean_S_min = mean(g$s_min), mean_k = mean(g$k, na.rm = TRUE),
               frac_ncls_1 = mean(g$n_clusters == 1),
               frac_ncls_2 = mean(g$n_clusters == TRUE_K),
               frac_ncls_ge3 = mean(g$n_clusters >= 3),
               median_ncls = median(g$n_clusters), stringsAsFactors = FALSE)))
  print(format(cc[order(cc$frac_ncls_1), ], digits = 4), row.names = FALSE)
  cat("\nby rule x method:\n")
  cm <- do.call(rbind, lapply(split(df, list(df$rule, df$method), drop = TRUE), function(g)
    data.frame(rule = g$rule[1], method = g$method[1], n = nrow(g),
               frac_ncls_1 = mean(g$n_clusters == 1),
               frac_ncls_2 = mean(g$n_clusters == TRUE_K), stringsAsFactors = FALSE)))
  print(format(cm[order(cm$method, cm$frac_ncls_1), ], digits = 3), row.names = FALSE)
  cat("\nfixed_0.0625 only, by setting x method (the real-data value):\n")
  f6 <- df[df$rule == "fixed_0.0625", ]
  fm <- do.call(rbind, lapply(split(f6, list(f6$setting_id, f6$method), drop = TRUE), function(g)
    data.frame(setting_id = g$setting_id[1], method = g$method[1], k = g$k[1],
               frac_ncls_1 = mean(g$n_clusters == 1),
               frac_ncls_2 = mean(g$n_clusters == TRUE_K), stringsAsFactors = FALSE)))
  print(format(fm[order(fm$method, fm$setting_id), ], digits = 3), row.names = FALSE)
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "status"
sel  <- if (length(args) >= 2 && nzchar(args[2]) && toupper(args[2]) != "ALL")
          strsplit(args[2], ",")[[1]] else NULL
nrp  <- if (length(args) >= 3) as.integer(args[3]) else 100L
bud  <- if (length(args) >= 4) as.numeric(args[4]) else 480
crs  <- if (length(args) >= 5) as.integer(args[5]) else 6L

cat(sprintf("55_wp2c_simulation_arm.R: mode=%s settings=%s reps=%d budget=%.0fs cores=%d\n",
            mode, if (is.null(sel)) "ALL" else paste(sel, collapse = ","), nrp, bud, crs))
switch(mode,
  probe     = do_probe(),
  run       = do_run(sel, nrp, bud, crs),
  status    = do_status(nrp),
  summarize = invisible(do_summarize()),
  report    = do_report(),
  stop(sprintf("unknown mode '%s'", mode)))
cat("55_wp2c_simulation_arm.R: done.\n")
