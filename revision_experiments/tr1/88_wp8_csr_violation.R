#!/usr/bin/env Rscript
# revision_experiments/tr1/88_wp8_csr_violation.R
#
# WP8 experiment 4 -- false positives under violated CSR, no outliers (R3.4;
# also grounds the R3.5 disclaimer). See WP8_PROTOCOL.md for the pre-declared
# design.
#
# R3.4: the radius-selection procedure and within-cluster validation rely on
# complete spatial randomness (a homogeneous Poisson process); the target
# scenarios (Gaussian clusters, non-uniform densities, density gradients)
# create a potential mismatch. This experiment generates NO-OUTLIER data
# under four single-population regimes -- a genuine CSR control, a density
# gradient (Beta(2,5) per coordinate), an isotropic Gaussian (`gauss_equal`,
# 2026-09-05 addition, see below) and an anisotropic ellipsoid
# (`gauss_unequal`) -- and reports the raw flag rate.
#
# n0 = 0 THROUGHOUT. evaluate()/count_scores2 indexes
# label_pred[(n-n0+1):n] for the TPR numerator; at n0=0 that is
# (n+1):n, a reversed two-element sequence reading past the vector's end.
# evaluate() is therefore NOT called here -- see WP8_PROTOCOL.md experiment 4.
# The reported quantity is flag_rate = mean(score >= threshold) over all n
# points.
#
# Settings: generator in {uniform_control, beta_gradient, gauss_equal,
# gauss_unequal} x d in {3, 10}, n = 200, n0 = 0. 8 settings, 100 reps, 9
# default methods + 4 MST-threshold-sweep variants per rep. Single row per
# cell -- has_result() on the main CSV gates skip/restart directly.
#
# Usage:
#   Rscript 88_wp8_csr_violation.R --smoke
#     One rep of one setting (uniform_control, d=3), all default + MST-sweep
#     methods, written to results/tr1/wp8/smoke/88_wp8_csr_violation.csv.
#   Rscript 88_wp8_csr_violation.R --summarize
#     Mean/SE of flag_rate per (generator, d, method), plus the paired
#     generator-minus-uniform_control delta at the same seed with its own
#     SE, and a best-MST-threshold-per-(generator,d) note. Writes
#     results/tr1/wp8/88_csr_violation_summary.csv.
#   Rscript 88_wp8_csr_violation.R [settings] [reps] [budget] [methods]

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
                      "LOF", "DBSCAN", "MST", "ODIN", "iForest")
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

# MST threshold sweep (2026-09-05, WP8_VERIFICATION.md item 88.3): the
# paper's real-data script tunes `thresh` per data set (1.05-1.6); this
# experiment sweeps the same range plus 2.0 (cheap, ~0.07 s/fit) and reports
# the best per generator ALONGSIDE the fixed 1.2 value used by "MST" in
# DEFAULT_METHODS (not duplicated -- 1.2 is covered by "MST" itself).
MST_SWEEP_VALS    <- c(1.05, 1.4, 1.6, 2.0)
MST_SWEEP_METHODS <- paste0("MST@", sprintf("%.2f", MST_SWEEP_VALS))
ALL_METHODS <- c(DEFAULT_METHODS, MST_SWEEP_METHODS)

ALL_THRESHOLDS <- REAL_DATA_THRESHOLDS
ALL_THRESHOLDS[MST_SWEEP_METHODS] <- 0.5  # same binarized-score threshold as "MST"

BASE_SEED <- 8801L
N_NOMINAL <- 200L

RESDIR    <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV  <- file.path(RESDIR, "88_csr_violation.csv")
CLUST_CSV <- file.path(RESDIR, "88_csr_violation_clusters.csv")
SUMM_CSV  <- file.path(RESDIR, "88_csr_violation_summary.csv")
SMOKE_DIR <- file.path(RESDIR, "smoke")
SMOKE_MAIN  <- file.path(SMOKE_DIR, "88_wp8_csr_violation.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "88_wp8_csr_violation_clusters.csv")

# ---------------------------------------------------------------------------
# Four single-population, no-outlier generators (WP8_PROTOCOL.md
# experiment 4, `gauss_equal` added 2026-09-05). Each returns n = N_NOMINAL
# points, no cluster/outlier structure -- true_cluster is always 1.
# ---------------------------------------------------------------------------
gen_uniform_control <- function(seed, n, d) {
  set.seed(seed)
  matrix(runif(n * d, 0, 1), nrow = n, ncol = d)
}

gen_beta_gradient <- function(seed, n, d) {
  set.seed(seed)
  matrix(rbeta(n * d, 2, 5), nrow = n, ncol = d)
}

# gauss_equal (2026-09-05, WP8_VERIFICATION.md item 88.1): an ISOTROPIC
# Gaussian (sigma = 1 in every coordinate), added so anisotropy is isolated
# by the gauss_unequal-minus-gauss_equal contrast; gauss_unequal alone
# confounded anisotropy with the Gaussian's own radial density gradient
# (any Gaussian has non-uniform density even with equal variances).
gen_gauss_equal <- function(seed, n, d) {
  set.seed(seed)
  mvrnorm(n, mu = rep(0, d), Sigma = diag(d))
}

gen_gauss_unequal <- function(seed, n, d) {
  set.seed(seed)
  sigma <- seq(0.3, 1.5, length.out = d)
  mvrnorm(n, mu = rep(0, d), Sigma = diag(sigma^2))
}

GENERATORS <- list(uniform_control = gen_uniform_control,
                    beta_gradient   = gen_beta_gradient,
                    gauss_equal     = gen_gauss_equal,
                    gauss_unequal   = gen_gauss_unequal)

build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(generator = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$generator, s$d)
  s
}
FULL_DIMS <- c(3L, 10L)
CANON <- build_settings(FULL_DIMS)

resolve_settings <- function(sel) {
  if (is.null(sel) || (length(sel) == 1 && toupper(sel) == "ALL")) return(CANON)
  out <- CANON[match(sel, CANON$setting_id), ]
  if (anyNA(out$setting_id)) stop(sprintf("unknown setting_id(s): %s",
                                          paste(sel[is.na(match(sel, CANON$setting_id))], collapse = ", ")))
  out
}

# ---------------------------------------------------------------------------
# One (setting, d, rep, method) cell. `method` may be an MST-sweep label
# ("MST@1.05" etc.); the underlying registry call is always "MST" with
# `thresh` parsed out of the label.
# ---------------------------------------------------------------------------
run_cell <- function(setting_id, generator, d, rep_id, seed, method, main_csv, clust_csv) {
  keys <- c(setting_id = setting_id, d = d, rep = rep_id, method = method)
  if (isTRUE(has_result(main_csv, keys))) {
    cat(sprintf("[skip] %s d=%d rep=%d method=%s\n", setting_id, d, rep_id, method))
    return(invisible("skip"))
  }

  X <- GENERATORS[[generator]](seed, N_NOMINAL, d)
  n <- nrow(X)
  Y <- rep(1L, n)  # no outliers

  is_mst_sweep <- grepl("^MST@", method)
  reg_name <- if (is_mst_sweep) "MST" else method
  extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN)
           else if (is_mst_sweep) list(thresh = as.numeric(sub("^MST@", "", method)))
           else list()

  out <- tryCatch({
    res <- do.call(METHOD_REGISTRY[[reg_name]], c(list(X = X, d = d, Y = Y), extra))
    stopifnot(length(res$score) == n, !anyNA(res$score))
    flagged <- res$score >= ALL_THRESHOLDS[[method]]
    # note = "-" (never "" or NA) on success -- see 86_wp8_outlier_types.R
    # for the full has_result()/type.convert() explanation. DBSCAN gets a
    # descriptive note instead (2026-09-05, WP8_VERIFICATION.md item 88.2):
    # it is fed Y = rep(1,n) so its oracle contamination is 0 and its flag
    # rate is 0 BY CONSTRUCTION, not because it detected anything -- kept in
    # the method list, footnoted wherever its numbers are tabulated.
    note <- if (identical(method, "DBSCAN"))
      "oracle contamination = 0; flag rate 0 by construction" else "-"
    list(flag_rate = mean(flagged), n_flagged = sum(flagged), res = res, status = "ok", note = note)
  }, error = function(e) list(flag_rate = NA_real_, n_flagged = NA_integer_, res = NULL,
                              status = "error", note = substr(conditionMessage(e), 1, 200)))

  row <- list(setting_id = setting_id, generator = generator, d = d, rep = rep_id, seed = seed,
              method = method, flag_rate = out$flag_rate, n_flagged = out$n_flagged,
              status = out$status, note = out$note)
  append_result(main_csv, row)

  if (identical(out$status, "ok") && method %in% MCCD_METHODS && !is.null(out$res$cluster)) {
    # Block write: all n cluster rows for this cell in ONE append_result call.
    idx <- seq_len(n)
    dc <- out$res$cluster[idx]
    append_result(clust_csv, list(
      setting_id = rep(setting_id, n), d = rep(d, n), rep = rep(rep_id, n), seed = rep(seed, n),
      method = rep(method, n), row_index = idx, true_cluster = rep(1L, n),
      detected_cluster = ifelse(is.na(dc), NA_integer_, dc)))
  }

  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE %s d=%d rep=%d method=%s flag_rate=%.3f\n",
                setting_id, d, rep_id, method, out$flag_rate))
  } else {
    cat(sprintf("CELL_FAIL %s d=%d rep=%d method=%s note=%s\n",
                setting_id, d, rep_id, method, out$note))
  }
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# --summarize: mean/SE of flag_rate per (generator, d, method); the paired
# generator-minus-uniform_control delta at the same (d, rep) with its own
# SE; a "best MST threshold per (generator,d)" note, where "best" =
# minimizes |mean paired delta vs uniform_control| among the swept
# thresholds (closest match to the CSR control's own flag rate = least
# distorted by the CSR violation), reported alongside the fixed-1.2 "MST"
# row for comparison, not in place of it.
# ---------------------------------------------------------------------------
do_summarize <- function() {
  if (!file.exists(MAIN_CSV)) { cat("no results yet\n"); return(invisible(NULL)) }
  df_all <- read.csv(MAIN_CSV, stringsAsFactors = FALSE)
  key <- paste(df_all$setting_id, df_all$d, df_all$rep, df_all$method, sep = "\r")
  df_all <- df_all[!duplicated(key, fromLast = TRUE), ]
  df <- df_all[df_all$status == "ok", ]

  rows <- list()
  for (dd in unique(df$d)) for (meth in unique(df$method[df$d == dd])) {
    ctrl <- df[df$d == dd & df$method == meth & df$generator == "uniform_control", ]
    for (gen in unique(df$generator[df$d == dd])) {
      sub <- df[df$d == dd & df$method == meth & df$generator == gen, ]
      if (!nrow(sub)) next
      mm <- merge(sub, ctrl[, c("rep", "flag_rate")], by = "rep", suffixes = c("", "_ctrl"))
      delta <- mm$flag_rate - mm$flag_rate_ctrl
      rows[[length(rows) + 1L]] <- data.frame(
        generator = gen, d = dd, method = meth, n_reps = nrow(sub),
        flag_rate_mean = mean(sub$flag_rate), flag_rate_se = sd(sub$flag_rate) / sqrt(nrow(sub)),
        n_paired = nrow(mm),
        delta_vs_control_mean = if (nrow(mm)) mean(delta) else NA_real_,
        delta_vs_control_se = if (nrow(mm) > 1) sd(delta) / sqrt(nrow(mm)) else NA_real_,
        n_error = sum(df_all$d == dd & df_all$method == meth & df_all$generator == gen & df_all$status != "ok"),
        stringsAsFactors = FALSE)
    }
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  S <- S[order(S$generator, S$d, S$method), ]
  write.csv(S, SUMM_CSV, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMM_CSV, nrow(S)))

  cat("\n=== MST: fixed 1.2 vs best-of-sweep per (generator, d), non-control only ===\n")
  mst_rows <- S[grepl("^MST", S$method) & S$generator != "uniform_control", ]
  for (gen in unique(mst_rows$generator)) for (dd in unique(mst_rows$d[mst_rows$generator == gen])) {
    blk <- mst_rows[mst_rows$generator == gen & mst_rows$d == dd, ]
    fixed <- blk[blk$method == "MST", ]
    sweep <- blk[blk$method != "MST", ]
    if (!nrow(sweep)) next
    best <- sweep[which.min(abs(sweep$delta_vs_control_mean)), ]
    cat(sprintf("  %-15s d=%-2d  fixed MST(1.2): flag_rate=%.3f delta=%+.3f  |  best %-9s flag_rate=%.3f delta=%+.3f\n",
                gen, dd,
                ifelse(nrow(fixed), fixed$flag_rate_mean, NA), ifelse(nrow(fixed), fixed$delta_vs_control_mean, NA),
                best$method, best$flag_rate_mean, best$delta_vs_control_mean))
  }
  invisible(S)
}

# ---------------------------------------------------------------------------
# --run: full (subset) grid, checkpointed by budget
# ---------------------------------------------------------------------------
do_run <- function(sel, n_reps, budget, methods) {
  st <- resolve_settings(sel)
  cat(sprintf("88_wp8_csr_violation: %d settings x %d reps x %d methods (budget %.0fs)\n",
              nrow(st), n_reps, length(methods), budget))
  t_start <- Sys.time()
  for (row in seq_len(nrow(st))) {
    s <- st[row, ]
    s_i <- match(s$setting_id, CANON$setting_id)
    for (rep_id in seq_len(n_reps)) {
      if (as.numeric(difftime(Sys.time(), t_start, units = "secs")) > budget) {
        cat(sprintf("[budget %gs exhausted -- rerun to continue] next up: %s rep %d\n",
                    budget, s$setting_id, rep_id))
        return(invisible("budget_stop"))
      }
      seed <- BASE_SEED + 100000L * s_i + rep_id
      for (m in methods) run_cell(s$setting_id, s$generator, s$d, rep_id, seed, m, MAIN_CSV, CLUST_CSV)
    }
  }
  cat("88_wp8_csr_violation: run complete.\n")
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
MODE_SUMM  <- "--summarize" %in% args
args <- args[!args %in% c("--smoke", "--summarize")]

if (!isTRUE(getOption("wp8.no_main", FALSE))) {
  if (MODE_SMOKE) {
    dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
    s <- CANON[CANON$generator == "uniform_control" & CANON$d == 3L, ][1, ]
    cat(sprintf("==== 88_wp8_csr_violation --smoke: %s, 1 rep, %d methods (incl. %d MST-sweep) ====\n",
                s$setting_id, length(ALL_METHODS), length(MST_SWEEP_METHODS)))
    s_i <- match(s$setting_id, CANON$setting_id)
    seed <- BASE_SEED + 100000L * s_i + 9001L
    for (m in ALL_METHODS) run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
    cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
    for (m in ALL_METHODS) run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
    cat("ALL_CELLS_COMPLETE\n")
    quit(save = "no", status = 0)
  } else if (MODE_SUMM) {
    do_summarize()
    quit(save = "no", status = 0)
  } else {
    sel     <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
    n_reps  <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 100L
    budget  <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 480
    methods <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else ALL_METHODS
    do_run(sel, n_reps, budget, methods)
  }
}
