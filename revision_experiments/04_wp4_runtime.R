#!/usr/bin/env Rscript
# revision_experiments/04_wp4_runtime.R
#
# WP4: wall-clock runtime / scalability experiments (task T4; reviewers R1 #2,
# R2 5th comment). Times all 9 R methods from the harness registry on
# synthetic Gaussian 2-cluster data (5% contamination) over two grids:
#
#   Grid A ("n"): n in {250, 500, 1000, 2000, 4000} at d = 20
#   Grid B ("d"): d in {10, 50, 100, 500, 1000}     at n = 1000
#
# USAGE
#   Rscript "revision_experiments/04_wp4_runtime.R" <reps> <mode>
#
#   <reps>  number of repetitions per (cell, method), e.g. 10
#   <mode>  one of:
#     all            export rep-1 datasets for every cell, then run the full
#                    R grid (both grids, cheapest cell first), then aggregate.
#     r_only         same grid run + aggregation, but skips the dataset-export
#                    step (use for reruns when results/wp4_data/ already
#                    exists, e.g. after the d=500/1000 quantile tables arrive).
#     cell=GRID:VAL  run a single cell. GRID is A|n (value = n at d=20) or
#                    B|d (value = d at n=1000). Examples: cell=A:2000,
#                    cell=n:4000, cell=B:500, cell=d:1000. Exports that
#                    cell's rep-1 dataset if missing, runs it, aggregates.
#     micro          hidden validation mode: three tiny cells
#                    (n=100,d=10), (n=200,d=10), (n=100,d=40). The d=40 cell
#                    exercises the SKIPPED_NO_TABLE path (no RK/NN quantile
#                    table exists at d=40). Micro results go to SEPARATE files
#                    (wp4_runtime_micro_raw.csv / wp4_runtime_micro.csv) so
#                    the production raw file stays clean.
#
# DESIGN NOTES (read before launching the full grid)
#
# Data generation. Same generator idiom as the first-cycle scripts, taken
# from simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/20d/
# 20d_2cls_n500_cont5%.R, parameterized by (n, d): two Gaussian clusters at
# mu1 = (3,...,3) and mu2 = mu1 + (3,0,...,0), per-cluster sigma scaled so
# ~1% of cluster mass falls outside the unit-ball-inspired radius, cluster
# radii jittered by runif(1, 0.7, 1.3), and 5% outliers drawn uniformly from
# a radius-5 ball around the mid-point, rejected until > 2 away from both
# cluster centers. One deliberate generalization: the reference script
# hard-codes n2 = round(n*(1-cont)*0.5) - 1 (an n=500-specific adjustment);
# here n2 = n - n1 - n0 so the three block sizes sum to exactly n for EVERY
# n. Row order is (cluster1, cluster2, outliers), labels 1 = regular,
# 0 = outlier (repo convention).
#
# Seeding rule (documented, deterministic): seed = 1000*rep + cell_index,
# where cell_index is the fixed enumeration below (Grid A cells 1-5 in
# ascending n, Grid B cells 6-10 in ascending d; micro cells 101-103). Each
# R rep therefore times a FRESH dataset drawn from the same distribution.
#
# Dataset export (FIRST ACTION of a run): each in-scope cell's rep-1 dataset
# (seed = 1000*1 + cell_index) is written to
#   results/wp4_data/<grid>_<cell_value>_rep1.csv   (V1..Vd + label)
# BEFORE any timing starts, so 05_wp4_runtime_pyod.py times the IDENTICAL
# data. Asymmetry (documented): only rep 1 per cell is exported; PyOD runs
# its 10 timing reps on that one file (its reps vary model seeds and system
# noise only), while the R reps vary the dataset draw as well. Existing
# files are not rewritten (the generator is deterministic, so a rewrite
# would be byte-identical anyway).
#
# Timing protocol: sequential execution in a single R process; no doParallel
# anywhere in the timed path (the harness registry is serial); OMP/MKL/
# OPENBLAS thread env vars pinned to 1; the registry's iForest entry is
# replaced by a local single-threaded variant (isotree defaults to all
# cores, which would break single-worker comparability -- same
# hyperparameters, nthreads = 1). Dataset generation and quantile-table
# loading are NOT timed. gc() runs before each timed rep.
#
# Reported columns per rep:
#   t_construct  CCD digraph-construction time alone (NA for baselines, and
#                NA when construction timing is skipped, see below)
#   t_total      the full single-pass pipeline time = THE runtime number
#                (for CCD methods the scorer internally re-runs
#                construction, so t_total alone is the honest end-to-end
#                cost; summing t_construct + t_total would double-count)
#   seconds      outer wall clock actually spent on the rep (= t_construct
#                run + t_total run for CCDs when both are timed)
#
# Construction-timing skip (>5 min rule): the registry's CCD entries run
# construction twice (once for t_construct, once inside the scorer). If any
# prior rep of the same (cell, method) took > 300 s wall clock, remaining
# reps call the scorer directly (single pass): t_construct = NA, status
# OK_NOCONSTRUCT. This avoids doubling already-long runs.
#
# Timeout discipline: each rep runs under setTimeLimit(elapsed = 1800,
# transient); if a rep exceeds 30 min it is aborted, a FLAGGED_TIMEOUT row
# is recorded (seconds = elapsed at abort), and the remaining reps of that
# (cell, method) are dropped -- completed reps stay, so the cell reduces to
# the reps already done (min 1 recorded row; a cell is NEVER silently
# dropped). On resume, a FLAGGED_TIMEOUT row blocks re-running that
# (cell, method) (delete its rows from the raw CSV to force a retry).
# Note: setTimeLimit fires at R-level interrupt points; long uninterruptible
# C calls (e.g. one dist() call) can overshoot slightly.
#
# Missing quantile tables: tables for d in {500, 1000} do not exist yet.
# CCD cells at those d get ONE marker row (rep = 0, status =
# SKIPPED_NO_TABLE) and are skipped. SKIPPED_NO_TABLE rows do NOT block
# resume: when the tables arrive, rerunning (mode r_only) fills in the real
# reps; the marker row is ignored by the timing aggregation.
#
# Checkpointing: every rep appends immediately to the raw CSV
# (results/wp4_runtime_raw.csv; micro mode uses wp4_runtime_micro_raw.csv)
# keyed by (grid, cell_value, method, rep). Fully resumable: on restart,
# reps already recorded with a non-SKIPPED status are skipped.
#
# Aggregation (end of every invocation): mean/sd per (cell, method) over
# status OK/OK_NOCONSTRUCT rows, plus a semicolon-joined status summary:
#   grid "n"     -> results/wp4_runtime_n.csv
#   grid "d"     -> results/wp4_runtime_d.csv
#   grid "micro" -> results/wp4_runtime_micro.csv
#
# UN-CCD radius-search direction follows the original simulation scripts:
# method = "ascend" for d <= 5, "descend" for d >= 10 (grep over
# simulations/outlyingness_scores/UNCCD_OOS_IOS/Simulation/Gaussian/*).
# MST uses cont = 0.05 (the datasets' true contamination); DBSCAN uses the
# registry's oracle-contamination convention (needs Y).

suppressPackageStartupMessages({
  library(here)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript \"revision_experiments/04_wp4_runtime.R\" <reps> <mode: all|r_only|cell=GRID:VALUE|micro>",
       call. = FALSE)
}
REPS <- as.integer(args[[1]])
MODE <- args[[2]]
stopifnot("reps must be a positive integer" = is.finite(REPS) && REPS >= 1)

TIMEOUT_SEC        <- 30 * 60   # hard per-rep cap
CONSTRUCT_SKIP_SEC <- 5 * 60    # prior-rep threshold for skipping t_construct

# Pin math-library threading BEFORE any heavy computation (comparability;
# R's own interpreter is single-threaded already).
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

source(here::here("revision_experiments/harness.R"))

RESULTS_DIR <- here::here("revision_experiments/results")
DATA_DIR    <- file.path(RESULTS_DIR, "wp4_data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Cell enumeration (FIXED -- the seeding rule depends on cell_index)
# ---------------------------------------------------------------------------
CELLS_MAIN <- data.frame(
  cell_index = 1:10,
  grid       = c(rep("n", 5), rep("d", 5)),
  cell_value = c("250", "500", "1000", "2000", "4000",
                 "10", "50", "100", "500", "1000"),
  n          = c(250, 500, 1000, 2000, 4000, rep(1000, 5)),
  d          = c(rep(20, 5), 10, 50, 100, 500, 1000),
  stringsAsFactors = FALSE
)

CELLS_MICRO <- data.frame(
  cell_index = c(101L, 102L, 103L),
  grid       = "micro",
  cell_value = c("n100_d10", "n200_d10", "n100_d40"),
  n          = c(100, 200, 100),
  d          = c(10, 10, 40),
  stringsAsFactors = FALSE
)

# Within-cell method order: cheap baselines first so partial results
# accumulate quickly; CCD methods (the potentially long ones) last.
METHOD_ORDER <- c("LOF", "DBSCAN", "MST", "ODIN", "iForest",
                  "RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
CCD_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")

seed_for <- function(cell_index, rep) 1000L * rep + cell_index

unccd_dir_for_d <- function(d) if (d <= 5) "ascend" else "descend"

# ---------------------------------------------------------------------------
# Mode -> cells + raw-CSV path
# ---------------------------------------------------------------------------
parse_cell_spec <- function(spec) {
  m <- regmatches(spec, regexec("^cell=([AaBbNnDd]):([0-9]+)$", spec))[[1]]
  if (length(m) != 3) {
    stop("Bad cell spec '", spec, "'. Expected cell=GRID:VALUE, e.g. cell=A:2000 or cell=d:500.", call. = FALSE)
  }
  g <- tolower(m[2]); v <- m[3]
  grid <- if (g %in% c("a", "n")) "n" else "d"
  hit <- CELLS_MAIN[CELLS_MAIN$grid == grid & CELLS_MAIN$cell_value == v, ]
  if (nrow(hit) != 1) {
    stop("cell=", spec, " does not match a defined cell. Grid A values: 250/500/1000/2000/4000; Grid B values: 10/50/100/500/1000.", call. = FALSE)
  }
  hit
}

if (MODE == "all") {
  CELLS <- CELLS_MAIN; DO_EXPORT <- TRUE
  RAW_CSV <- file.path(RESULTS_DIR, "wp4_runtime_raw.csv")
} else if (MODE == "r_only") {
  CELLS <- CELLS_MAIN; DO_EXPORT <- FALSE
  RAW_CSV <- file.path(RESULTS_DIR, "wp4_runtime_raw.csv")
} else if (grepl("^cell=", MODE)) {
  CELLS <- parse_cell_spec(MODE); DO_EXPORT <- TRUE
  RAW_CSV <- file.path(RESULTS_DIR, "wp4_runtime_raw.csv")
} else if (MODE == "micro") {
  CELLS <- CELLS_MICRO; DO_EXPORT <- TRUE
  RAW_CSV <- file.path(RESULTS_DIR, "wp4_runtime_micro_raw.csv")
} else {
  stop("Unknown mode '", MODE, "'. Use all | r_only | cell=GRID:VALUE | micro.", call. = FALSE)
}

cat(sprintf("=== 04_wp4_runtime.R === reps=%d mode=%s cells=%d raw=%s\n",
            REPS, MODE, nrow(CELLS), RAW_CSV))

# ---------------------------------------------------------------------------
# Data generator (parameterized copy of the reference block; see header)
# ---------------------------------------------------------------------------
gen_gaussian_2cls_runtime <- function(n, d, cont = 0.05, seed,
                                       cls_dis = 3, otl_dis = 2,
                                       r_min = 0.7, r_max = 1.3,
                                       noise_level = 0.01) {
  n1 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)
  n2 <- n - n1 - n0   # exact total (generalizes the reference's n=500 "-1")

  mu1 <- rep(3, d)
  mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- colMeans(rbind(mu1, mu2))

  sigma <- 1 / sqrt(qchisq(1 - noise_level, d))

  set.seed(seed)
  data1 <- MASS::mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 <- MASS::mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)

  i <- 0
  outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu   # rpoisball.unit: R/ccds/Kest.R (sourced via harness)
    r1 <- sqrt(sum((temp - mu1)^2))
    r2 <- sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis && r2 > otl_dis) {
      outlier <- rbind(outlier, temp)
      i <- i + 1
    }
  }
  rownames(outlier) <- NULL

  X <- rbind(data1, data2, outlier)
  colnames(X) <- paste0("V", seq_len(d))
  Y <- c(rep(1, n1 + n2), rep(0, n0))   # 1 = regular, 0 = outlier
  list(X = X, Y = Y)
}

# ---------------------------------------------------------------------------
# Rep-1 dataset export (FIRST ACTION; consumed by 05_wp4_runtime_pyod.py)
# ---------------------------------------------------------------------------
export_path_for <- function(cell) {
  file.path(DATA_DIR, sprintf("%s_%s_rep1.csv", cell$grid, cell$cell_value))
}

export_cells_rep1 <- function(cells) {
  for (k in seq_len(nrow(cells))) {
    cell <- cells[k, ]
    path <- export_path_for(cell)
    if (file.exists(path)) {
      cat(sprintf("[export] exists, kept: %s\n", basename(path)))
      next
    }
    dat <- gen_gaussian_2cls_runtime(cell$n, cell$d, seed = seed_for(cell$cell_index, 1L))
    write.csv(data.frame(dat$X, label = dat$Y), path, row.names = FALSE)
    cat(sprintf("[export] wrote %s (n=%d, d=%d)\n", basename(path), cell$n, cell$d))
  }
}

# ---------------------------------------------------------------------------
# Quantile-table availability (mirrors harness get_simul() path resolution)
# ---------------------------------------------------------------------------
ccd_variant_for_method <- function(m) {
  if (startsWith(m, "RKCCD")) "RK" else if (startsWith(m, "UNCCD")) "NN" else NA_character_
}

has_quant_table <- function(variant, d) {
  if (variant == "RK") {
    file.exists(file.path(RK_QUANT_TABLE_DIR,
                          sprintf("RK-test-simul_%dd_%s%%.RData", d, rk_quant_for_d(d))))
  } else {
    file.exists(file.path(NN_QUANT_TABLE_DIR,
                          sprintf("NN-test-simul_%dd_%s%%.RData", d, nn_quant_for_d(d))))
  }
}

# ---------------------------------------------------------------------------
# Method dispatch
# ---------------------------------------------------------------------------

# Local single-threaded iForest (registry's entry lets isotree grab all
# cores; identical hyperparameters otherwise). Timed like the registry.
iforest_1thread <- function(X, d, Y = NULL, seed = 1) {
  X <- as.matrix(X)
  t0 <- Sys.time()
  set.seed(seed)
  model <- isotree::isolation.forest(X, ntrees = 1000,
                                     sample_size = min(256, nrow(X)),
                                     nthreads = 1)
  score <- as.numeric(predict(model, X))
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

# Single-pass CCD scorer calls (used when construction timing is skipped by
# the >5 min rule). Same arguments the registry entries pass internally.
ccd_single_pass <- function(m, X, d) {
  X <- as.matrix(X)
  variant <- ccd_variant_for_method(m)
  tab <- get_simul(variant, d)   # table load NOT timed (matches registry)
  um <- unccd_dir_for_d(d)
  t0 <- Sys.time()
  score <- switch(m,
    "RKCCD-OOS" = RKCCD_OOS(datax = X, simul = tab$simul, d = d, quant = tab$quant),
    "RKCCD-IOS" = RKCCD_IOS(datax = X, simul = tab$simul, d = d, quant = tab$quant, min.cls = 0),
    "UNCCD-OOS" = NNCCD_OOS(datax = X, simul = tab$simul, method = um, d = d),
    "UNCCD-IOS" = NNCCD_IOS(datax = X, simul = tab$simul, method = um, d = d, min.cls = 0),
    stop("ccd_single_pass(): not a CCD method: ", m)
  )
  t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  list(score = score, t_construct = NA_real_, t_total = t_total)
}

run_method_rep <- function(m, X, Y, d, seed, skip_construct) {
  if (m == "iForest") return(iforest_1thread(X, d, Y, seed = seed))
  if (m %in% CCD_METHODS && skip_construct) return(ccd_single_pass(m, X, d))
  extra <- list()
  if (m %in% c("UNCCD-OOS", "UNCCD-IOS")) extra$method <- unccd_dir_for_d(d)
  if (m == "MST") extra$cont <- 0.05
  do.call(METHOD_REGISTRY[[m]], c(list(X = X, d = d, Y = Y), extra))
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}

# ---------------------------------------------------------------------------
# Checkpoint history (status-aware; harness has_result() is key-only)
# ---------------------------------------------------------------------------
ROW_COLS <- c("grid", "cell_value", "n", "d", "method", "rep", "seed",
              "t_construct", "t_total", "seconds", "status", "timestamp")

HIST <- new.env()
HIST$df <- if (file.exists(RAW_CSV)) {
  read.csv(RAW_CSV, stringsAsFactors = FALSE)
} else {
  NULL
}

hist_subset <- function(grid, cell_value, method) {
  df <- HIST$df
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df[df$grid == grid & as.character(df$cell_value) == as.character(cell_value) &
       df$method == method, , drop = FALSE]
}

rep_recorded <- function(grid, cell_value, method, rep) {
  sub <- hist_subset(grid, cell_value, method)
  if (is.null(sub) || nrow(sub) == 0) return(FALSE)
  # SKIPPED_NO_TABLE never blocks (rerun picks the cell up once tables exist)
  any(sub$rep == rep & sub$status != "SKIPPED_NO_TABLE")
}

method_flagged_timeout <- function(grid, cell_value, method) {
  sub <- hist_subset(grid, cell_value, method)
  if (is.null(sub) || nrow(sub) == 0) return(FALSE)
  any(sub$status == "FLAGGED_TIMEOUT")
}

skip_marker_exists <- function(grid, cell_value, method) {
  sub <- hist_subset(grid, cell_value, method)
  if (is.null(sub) || nrow(sub) == 0) return(FALSE)
  any(sub$status == "SKIPPED_NO_TABLE")
}

max_prior_seconds <- function(grid, cell_value, method) {
  sub <- hist_subset(grid, cell_value, method)
  if (is.null(sub) || nrow(sub) == 0) return(0)
  s <- suppressWarnings(as.numeric(sub$seconds))
  s <- s[is.finite(s)]
  if (length(s) == 0) 0 else max(s)
}

append_and_track <- function(row) {
  stopifnot(identical(names(row), ROW_COLS))
  append_result(RAW_CSV, row)
  new_df <- as.data.frame(row, stringsAsFactors = FALSE)
  HIST$df <- if (is.null(HIST$df)) new_df else rbind(HIST$df, new_df)
}

make_row <- function(cell, method, rep, seed, t_construct, t_total, seconds, status) {
  list(grid = cell$grid, cell_value = cell$cell_value,
       n = cell$n, d = cell$d, method = method,
       rep = as.integer(rep), seed = if (is.na(seed)) NA_integer_ else as.integer(seed),
       t_construct = if (is.na(t_construct)) NA_real_ else round(t_construct, 4),
       t_total = if (is.na(t_total)) NA_real_ else round(t_total, 4),
       seconds = if (is.na(seconds)) NA_real_ else round(seconds, 4),
       status = status,
       timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
GRID_AGG_FILE <- c(n = "wp4_runtime_n.csv", d = "wp4_runtime_d.csv",
                   micro = "wp4_runtime_micro.csv")

aggregate_wp4 <- function(raw_csv) {
  if (!file.exists(raw_csv)) {
    cat("[aggregate] no raw file yet; nothing to aggregate\n")
    return(invisible(NULL))
  }
  df <- read.csv(raw_csv, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(invisible(NULL))
  for (col in c("t_construct", "t_total", "seconds")) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  groups <- split(seq_len(nrow(df)),
                  interaction(df$grid, df$cell_value, df$method, drop = TRUE))
  agg <- do.call(rbind, lapply(groups, function(idx) {
    sub <- df[idx, ]
    ok <- sub[sub$status %in% c("OK", "OK_NOCONSTRUCT"), ]
    msd <- function(x) {
      x <- x[is.finite(x)]
      c(mean = if (length(x)) mean(x) else NA_real_,
        sd   = if (length(x) > 1) sd(x) else NA_real_)
    }
    tt <- msd(ok$t_total); tc <- msd(ok$t_construct); ss <- msd(ok$seconds)
    data.frame(grid = sub$grid[1], cell_value = sub$cell_value[1],
               n = sub$n[1], d = sub$d[1], method = sub$method[1],
               n_reps_ok = nrow(ok),
               mean_t_total = round(tt["mean"], 4), sd_t_total = round(tt["sd"], 4),
               mean_t_construct = round(tc["mean"], 4), sd_t_construct = round(tc["sd"], 4),
               mean_seconds = round(ss["mean"], 4), sd_seconds = round(ss["sd"], 4),
               statuses = paste(sort(unique(sub$status)), collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  rownames(agg) <- NULL
  for (g in unique(agg$grid)) {
    sub <- agg[agg$grid == g, ]
    ord_val <- suppressWarnings(as.numeric(sub$cell_value))
    if (any(is.na(ord_val))) ord_val <- rank(sub$cell_value)  # micro labels
    sub <- sub[order(ord_val, match(sub$method, METHOD_ORDER)), ]
    out_name <- GRID_AGG_FILE[[g]]
    if (is.null(out_name)) out_name <- sprintf("wp4_runtime_%s.csv", g)
    out_path <- file.path(RESULTS_DIR, out_name)
    write.csv(sub, out_path, row.names = FALSE)
    cat(sprintf("[aggregate] wrote %s (%d rows)\n", out_path, nrow(sub)))
  }
  invisible(agg)
}

# ---------------------------------------------------------------------------
# Main loop (cells in given order = cheapest first by construction)
# ---------------------------------------------------------------------------
if (DO_EXPORT) {
  cat("---- exporting rep-1 datasets (first action) ----\n")
  export_cells_rep1(CELLS)
} else {
  missing_exports <- !file.exists(vapply(seq_len(nrow(CELLS)),
                                          function(k) export_path_for(CELLS[k, ]), ""))
  if (any(missing_exports)) {
    cat("[warn] r_only mode but some rep-1 exports are missing; run mode 'all' (or 'cell=') to create them for the PyOD script.\n")
  }
}

run_start <- Sys.time()

for (k in seq_len(nrow(CELLS))) {
  cell <- CELLS[k, ]
  cat(sprintf("\n==== cell %d/%d: grid=%s value=%s (n=%d, d=%d) ====\n",
              k, nrow(CELLS), cell$grid, cell$cell_value, cell$n, cell$d))

  for (m in METHOD_ORDER) {
    variant <- ccd_variant_for_method(m)

    # SKIPPED_NO_TABLE path (CCD methods at d without a quantile table)
    if (!is.na(variant) && !has_quant_table(variant, cell$d)) {
      if (!skip_marker_exists(cell$grid, cell$cell_value, m)) {
        append_and_track(make_row(cell, m, rep = 0L, seed = NA,
                                  t_construct = NA, t_total = NA, seconds = NA,
                                  status = "SKIPPED_NO_TABLE"))
        cat(sprintf("  %-10s SKIPPED_NO_TABLE (no %s table at d=%d; rerun after table generation)\n",
                    m, variant, cell$d))
      } else {
        cat(sprintf("  %-10s SKIPPED_NO_TABLE (marker already recorded)\n", m))
      }
      next
    }

    if (method_flagged_timeout(cell$grid, cell$cell_value, m)) {
      cat(sprintf("  %-10s previously FLAGGED_TIMEOUT -- skipping (delete its rows in %s to retry)\n",
                  m, basename(RAW_CSV)))
      next
    }

    for (r in seq_len(REPS)) {
      if (rep_recorded(cell$grid, cell$cell_value, m, r)) {
        cat(sprintf("  %-10s rep %d/%d [checkpoint skip]\n", m, r, REPS))
        next
      }
      seed <- seed_for(cell$cell_index, r)
      dat <- gen_gaussian_2cls_runtime(cell$n, cell$d, seed = seed)  # untimed
      skip_construct <- (m %in% CCD_METHODS) &&
        (max_prior_seconds(cell$grid, cell$cell_value, m) > CONSTRUCT_SKIP_SEC)

      invisible(gc(FALSE))
      t0 <- Sys.time()
      res <- tryCatch(
        with_timeout(function() run_method_rep(m, dat$X, dat$Y, cell$d, seed, skip_construct),
                     TIMEOUT_SEC),
        error = function(e) e
      )
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (inherits(res, "error")) {
        msg <- conditionMessage(res)
        is_timeout <- grepl("elapsed time limit|CPU time limit", msg) ||
          elapsed >= TIMEOUT_SEC - 1
        status <- if (is_timeout) "FLAGGED_TIMEOUT" else
          paste0("ERROR: ", substr(gsub("[\r\n,]+", " ", msg), 1, 160))
        append_and_track(make_row(cell, m, r, seed, NA, NA, elapsed, status))
        cat(sprintf("  %-10s rep %d/%d %s after %.1fs -- remaining reps of this (cell, method) dropped\n",
                    m, r, REPS, status, elapsed))
        break   # never silently dropped: the row above records the event
      }

      if (length(res$score) != cell$n) {
        append_and_track(make_row(cell, m, r, seed, NA, NA, elapsed,
                                  sprintf("ERROR: score length %d != n %d", length(res$score), cell$n)))
        cat(sprintf("  %-10s rep %d/%d bad score length -- dropped remaining reps\n", m, r, REPS))
        break
      }

      status <- if (skip_construct) "OK_NOCONSTRUCT" else "OK"
      append_and_track(make_row(cell, m, r, seed,
                                res$t_construct, res$t_total, elapsed, status))
      cat(sprintf("  %-10s rep %d/%d t_construct=%s t_total=%.3fs wall=%.3fs [%s]\n",
                  m, r, REPS,
                  ifelse(is.na(res$t_construct), "NA", sprintf("%.3fs", res$t_construct)),
                  res$t_total, elapsed, status))
    }
  }
}

cat(sprintf("\n---- grid done in %.1f min; aggregating ----\n",
            as.numeric(difftime(Sys.time(), run_start, units = "mins"))))
aggregate_wp4(RAW_CSV)
cat("=== 04_wp4_runtime.R DONE ===\n")
