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
# under three single-population regimes -- a genuine CSR control, a density
# gradient (Beta(2,5) per coordinate), and an anisotropic ellipsoid (unequal
# per-axis Gaussian variance) -- and reports the raw flag rate.
#
# n0 = 0 THROUGHOUT. evaluate()/count_scores2 indexes
# label_pred[(n-n0+1):n] for the TPR numerator; at n0=0 that is (n+1):n, a
# reversed two-element sequence reading past the vector's end. evaluate() is
# therefore NOT called here -- see WP8_PROTOCOL.md experiment 4. The reported
# quantity is flag_rate = mean(score >= threshold) over all n points.
#
# Settings: generator in {uniform_control, beta_gradient, gauss_unequal} x
# d in {3, 10}, n = 200, n0 = 0. 6 settings, 100 reps, 9 methods (4 MCCD + 5
# baselines) per rep. Single row per cell -- has_result() on the main CSV
# gates skip/restart directly.
#
# Usage:
#   Rscript 88_wp8_csr_violation.R --smoke
#     One rep of one setting (uniform_control, d=3), all 9 default methods,
#     written to results/tr1/wp8/smoke/88_wp8_csr_violation.csv.
#   Rscript 88_wp8_csr_violation.R [reps] [dims] [methods]

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

ALL_THRESHOLDS <- REAL_DATA_THRESHOLDS
DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
                      "LOF", "DBSCAN", "MST", "ODIN", "iForest")
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

BASE_SEED <- 8801L
N_NOMINAL <- 200L

RESDIR    <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV  <- file.path(RESDIR, "88_csr_violation.csv")
CLUST_CSV <- file.path(RESDIR, "88_csr_violation_clusters.csv")
SMOKE_DIR <- file.path(RESDIR, "smoke")
SMOKE_MAIN  <- file.path(SMOKE_DIR, "88_wp8_csr_violation.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "88_wp8_csr_violation_clusters.csv")

# ---------------------------------------------------------------------------
# Three single-population, no-outlier generators (WP8_PROTOCOL.md
# experiment 4). Each returns n = N_NOMINAL points, no cluster/outlier
# structure -- true_cluster is always 1.
# ---------------------------------------------------------------------------
gen_uniform_control <- function(seed, n, d) {
  set.seed(seed)
  matrix(runif(n * d, 0, 1), nrow = n, ncol = d)
}

gen_beta_gradient <- function(seed, n, d) {
  set.seed(seed)
  matrix(rbeta(n * d, 2, 5), nrow = n, ncol = d)
}

gen_gauss_unequal <- function(seed, n, d) {
  set.seed(seed)
  sigma <- seq(0.3, 1.5, length.out = d)
  mvrnorm(n, mu = rep(0, d), Sigma = diag(sigma^2))
}

GENERATORS <- list(uniform_control = gen_uniform_control,
                    beta_gradient   = gen_beta_gradient,
                    gauss_unequal   = gen_gauss_unequal)

build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(generator = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$generator, s$d)
  s
}

# ---------------------------------------------------------------------------
# One (setting, d, rep, method) cell
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

  extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN) else list()
  out <- tryCatch({
    res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra))
    stopifnot(length(res$score) == n, !anyNA(res$score))
    flagged <- res$score >= ALL_THRESHOLDS[[method]]
    # note = "-" (never "" or NA) on success: has_result() treats ANY NA
    # in a non-key column as an incomplete row, and an all-"" column
    # round-trips through read.csv's type.convert() as logical NA (see
    # 86_wp8_outlier_types.R for the full explanation), so a non-empty
    # placeholder is required, not just a non-NA one.
    list(flag_rate = mean(flagged), n_flagged = sum(flagged), res = res, status = "ok",
         note = "-")
  }, error = function(e) list(flag_rate = NA_real_, n_flagged = NA_integer_, res = NULL,
                              status = "error", note = substr(conditionMessage(e), 1, 200)))

  row <- list(setting_id = setting_id, generator = generator, d = d, rep = rep_id, seed = seed,
              method = method, flag_rate = out$flag_rate, n_flagged = out$n_flagged,
              status = out$status, note = out$note)
  append_result(main_csv, row)

  if (identical(out$status, "ok") && method %in% MCCD_METHODS && !is.null(out$res$cluster)) {
    for (i in seq_len(n)) {
      append_result(clust_csv, list(
        setting_id = setting_id, d = d, rep = rep_id, seed = seed, method = method,
        row_index = i, true_cluster = 1L,
        detected_cluster = if (is.na(out$res$cluster[i])) NA_integer_ else out$res$cluster[i]))
    }
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
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
args <- args[args != "--smoke"]

if (MODE_SMOKE) {
  dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
  st <- build_settings(3L)
  s <- st[st$generator == "uniform_control", ][1, ]
  cat(sprintf("==== 88_wp8_csr_violation --smoke: %s, 1 rep, %d methods ====\n",
              s$setting_id, length(DEFAULT_METHODS)))
  seed <- BASE_SEED + 100000L * 1L + 9001L
  for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
  cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
  for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$generator, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

N_REPS  <- if (length(args) >= 1 && nzchar(args[1])) as.integer(args[1]) else 100L
DIMS    <- if (length(args) >= 2 && nzchar(args[2])) as.integer(strsplit(args[2], ",")[[1]]) else c(3L, 10L)
METHODS <- if (length(args) >= 3 && nzchar(args[3])) strsplit(args[3], ",")[[1]] else DEFAULT_METHODS

SETTINGS <- build_settings(DIMS)
cat(sprintf("88_wp8_csr_violation: %d settings x %d reps x %d methods\n",
            nrow(SETTINGS), N_REPS, length(METHODS)))

for (s_i in seq_len(nrow(SETTINGS))) {
  s <- SETTINGS[s_i, ]
  for (rep_id in seq_len(N_REPS)) {
    seed <- BASE_SEED + 100000L * s_i + rep_id
    for (m in METHODS) run_cell(s$setting_id, s$generator, s$d, rep_id, seed, m, MAIN_CSV, CLUST_CSV)
  }
}
cat("88_wp8_csr_violation: done.\n")
