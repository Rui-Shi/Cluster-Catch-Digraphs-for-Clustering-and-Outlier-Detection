#!/usr/bin/env Rscript
# revision_experiments/tr1/91_wp6_runtime.R
#
# WP6: runtime and memory grid (AE.5, R3.10 runtime/memory/scalability half,
# R5.5, R5.6 scalability angle). See tr1/WP6_PROTOCOL.md, declared before
# this file was written -- that document is the specification; this is the
# implementation and must not diverge from it silently.
#
# Times the paper's 4 proposed detectors (U-MCCD, SU-MCCD, UN-MCCD,
# SUN-MCCD, via tr1/wp0_mccd_methods.R) plus the 5 "old" baselines (LOF,
# DBSCAN, MST, ODIN, iForest) on synthetic uniform two-cluster data (5%
# contamination), over two grids:
#   Sweep A ("n"): n in {100, 250, 500, 1000, 2000} at d = 10
#   Sweep B ("d"): d in {5, 10, 50, 100}            at n = 500
# (n=500, d=10) is shared by both sweeps and computed once. Also exports one
# rep-1 dataset per cell to results/tr1/wp6/data/ for
# 91b_wp6_runtime_py.py's 8 WP4-competitor timing grid.
#
# USAGE
#   Rscript 91_wp6_runtime.R --smoke
#     One cell (n=100, d=10), 1 rep, all 9 R-side methods, written to
#     results/tr1/wp6/smoke/. Safe to run twice -- the second call should
#     print "[checkpoint skip]" for every (method, rep).
#   Rscript 91_wp6_runtime.R --idle-check
#     Lists other Rscript/Rterm/R/python/pythonw processes (by PID, excluding
#     this process) and exits -- 0 if none found, 1 otherwise. Runs nothing.
#   Rscript 91_wp6_runtime.R --summarize
#     Reads the raw CSV, writes results/tr1/wp6/91_wp6_runtime_n.csv,
#     _d.csv (median/CV per cell/method) and 91_wp6_slope.csv (log-log
#     slope +/- SE per MCCD method, sweep A only; see WP6_PROTOCOL.md #7).
#   Rscript 91_wp6_runtime.R [--grid=n|d|all] [--reps=10] [--force]
#     Runs the real grid (default: both sweeps, 10 reps). Refuses to start
#     unless the idle-check passes, or --force is given.
#
# NOT RUN FOR REAL by the session that wrote this file (see
# WP6_PROTOCOL.md #10): needs an idle machine, and four background R jobs
# were running at write time. Only --smoke was exercised.

suppressPackageStartupMessages(library(here))

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
raw_args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(name) any(raw_args == name)
opt_val  <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), raw_args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

MODE <- if (has_flag("--smoke")) {
  "smoke"
} else if (has_flag("--summarize")) {
  "summarize"
} else if (has_flag("--idle-check")) {
  "idle_check"
} else if (has_flag("--check-exports")) {
  "check_exports"     # D fix dry check (WP6_VERIFICATION.md): no idle-check, no compute
} else {
  "run"
}
GRID_SEL <- opt_val("grid", "all")           # n | d | all
REPS     <- as.integer(opt_val("reps", "10"))
FORCE    <- has_flag("--force")

stopifnot("REPS must be a positive integer" = is.finite(REPS) && REPS >= 1)
stopifnot("--grid must be one of n, d, all" = GRID_SEL %in% c("n", "d", "all"))

TIMEOUT_SEC <- 30 * 60   # per-rep cap; WP6_PROTOCOL.md #8 (defensive only)

# ---------------------------------------------------------------------------
# Single-threading (before harness.R is sourced -- matches
# 04_wp4_runtime.R's ordering). See WP6_PROTOCOL.md #4 for the verification
# this is a defensive no-op given this R install's BLAS.
# ---------------------------------------------------------------------------
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
THREADS <- 1L

source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))

# ---------------------------------------------------------------------------
# A1 (WP6_VERIFICATION.md, blocker) -- memoize get_simul() in THIS DRIVER
# ONLY; harness.R is untouched. umccd_method/sumccd_method/unmccd_method/
# sunmccd_method (wp0_mccd_methods.R:276,293,...) each call get_simul()
# themselves before their own internal t0 -- but get_simul() (harness.R
# :204-243) load()s the table's .RData file from disk on EVERY call, with no
# cache, and this driver's own outer proc.time() (the timed region, below)
# wraps the whole wrapper call, so an uncached load inside the wrapper still
# lands inside the span this driver reports as time_s. Measured before this
# fix: U-MCCD outer (this driver's) elapsed 3.96s vs the wrapper's own
# internal t_total 1.28s at n=250 -- 2.68s of load time counted as "runtime";
# the original smoke row showed time_s 2.95 with table_load_s (this driver's
# separate, already-untimed measurement) 2.66. A constant ~2.7s offset like
# that turns a true log-log slope of ~1.0-1.6 into a measured 0.36. One
# cache entry only (~763 MB per RK table); a new (variant, d, quant) key
# evicts the old entry first so memory never compounds across cells. The
# short-table abort (check_simul_extent) still runs on every call, cached or
# not, per get_simul()'s own contract.
local({
  .orig <- get_simul
  .cache <- new.env(parent = emptyenv())
  get_simul <<- function(variant = c("RK", "NN"), d, quant = NULL, n = NULL) {
    variant <- match.arg(variant)
    key <- paste(variant, d, quant, sep = "|")
    if (!exists(key, envir = .cache, inherits = FALSE)) {
      rm(list = ls(.cache, all.names = TRUE), envir = .cache)  # one entry: ~763 MB each
      assign(key, .orig(variant, d, quant, n = NULL), envir = .cache)
    }
    tab <- get(key, envir = .cache, inherits = FALSE)
    check_simul_extent(tab$simul, variant, d, n, tab$file)   # keep the loud short-table abort
    tab
  }
})

cat(sprintf("[config] mode=%s grid=%s reps=%d force=%s threads=%d\n",
            MODE, GRID_SEL, REPS, FORCE, THREADS))

# ---------------------------------------------------------------------------
# Idle-machine check (WP6_PROTOCOL.md #6)
# ---------------------------------------------------------------------------
idle_check <- function() {
  self_pid <- Sys.getpid()
  # On Windows, `Rscript.exe` is a launcher that spawns the actual R process, so
  # Sys.getpid() is the CHILD pid and the launcher itself shows up in tasklist as
  # an "other" Rscript.exe. Found 2026-09-06: the check could never pass. Walk
  # the parent chain (PowerShell CIM) and exclude every ancestor pid too.
  ancestors <- integer(0)
  pid_cur <- self_pid
  for (k in 1:4) {
    pp <- tryCatch(suppressWarnings(as.integer(system2(
      "powershell", args = c("-NoProfile", "-Command",
        sprintf("(Get-CimInstance Win32_Process -Filter \"ProcessId=%d\").ParentProcessId", pid_cur)),
      stdout = TRUE, stderr = FALSE))), error = function(e) NA_integer_)
    pp <- pp[!is.na(pp)]
    if (length(pp) == 0 || pp[1] <= 0) break
    ancestors <- c(ancestors, pp[1]); pid_cur <- pp[1]
  }
  exclude <- c(self_pid, ancestors)
  out <- tryCatch(system2("tasklist", args = c("/FO", "CSV", "/NH"),
                          stdout = TRUE, stderr = TRUE),
                   error = function(e) character(0))
  if (length(out) == 0) {
    cat("[idle-check] tasklist unavailable -- cannot verify; proceeding as if idle.\n")
    return(invisible(TRUE))
  }
  watch <- c("Rscript.exe", "Rterm.exe", "R.exe", "python.exe", "pythonw.exe")
  hits <- list()
  for (ln in out) {
    fields <- tryCatch(unlist(utils::read.csv(text = ln, header = FALSE,
                                              stringsAsFactors = FALSE)),
                       error = function(e) NULL)
    if (is.null(fields) || length(fields) < 2) next
    img <- trimws(fields[[1]]); pid <- suppressWarnings(as.integer(fields[[2]]))
    if (is.na(pid) || pid %in% exclude) next
    if (img %in% watch) hits[[length(hits) + 1L]] <- c(image = img, pid = pid)
  }
  if (length(hits) == 0) {
    cat("[idle-check] no other Rscript/Rterm/R/python/pythonw processes found. Machine looks idle.\n")
    return(invisible(TRUE))
  }
  cat("[idle-check] OTHER PROCESSES FOUND (machine is NOT idle):\n")
  for (h in hits) cat(sprintf("    %s  PID=%s\n", h[["image"]], h[["pid"]]))
  invisible(FALSE)
}

if (MODE == "idle_check") {
  ok <- idle_check()
  quit(status = if (isTRUE(ok)) 0L else 1L, save = "no")
}

if (MODE == "run" && !FORCE) {
  ok <- idle_check()
  if (!isTRUE(ok)) {
    stop("Refusing to start the real WP6 grid: other Rscript/python processes are present ",
        "(see list above). Wait for them to finish, or pass --force to override.",
        call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Generator (WP6_PROTOCOL.md #2 -- reproduced from 55_wp2c_simulation_arm.R's
# gen_uniform via 85_wp8_boundary_fp.R's gen_uniform_boundary, generalized to
# an exact n1+n2+n0 = n split at every n, same fix 04_wp4_runtime.R applied).
# rpoisball.unit comes from R/ccds/Kest.R, already in scope via harness.R's
# sourcing of methods/outlyingness_scores/{RKCCD,UNCCD}_OOS_IOS.R.
# ---------------------------------------------------------------------------
gen_uniform_wp6 <- function(seed, n, d, cont = 0.05) {
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- (mu1 + mu2) / 2
  n1 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)
  n2 <- n - n1 - n0
  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) +
             matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) +
             matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i <- 0L; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis && sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1L
    }
  }
  if (!is.null(outlier)) rownames(outlier) <- NULL
  X <- rbind(data1, data2, outlier)
  colnames(X) <- paste0("V", seq_len(d))
  Y <- c(rep(1L, n1 + n2), rep(0L, n0))
  list(X = X, Y = Y)
}

# ---------------------------------------------------------------------------
# Cell table (WP6_PROTOCOL.md #3). cell_index is FIXED -- the seeding rule
# depends on it; never reorder or renumber this table.
# ---------------------------------------------------------------------------
# D fix (WP6_VERIFICATION.md, blocker): cell_value used to be n's value for
# EVERY cell (including the three d-only cells 6-8), so export_path_for()
# wrote "d_500_rep1.csv" for all three -- one file, silently overwritten
# twice -- and the Python d-sweep collapsed to a single point (d=5). The
# d-only cells now carry their own d value as cell_value.
CELLS <- data.frame(
  cell_index = 1:8,
  n          = c(100, 250, 500, 1000, 2000, 500, 500, 500),
  d          = c(10, 10, 10, 10, 10, 5, 50, 100),
  cell_value = c("100", "250", "500", "1000", "2000", "5", "50", "100"),
  in_n       = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  in_d       = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

seed_for <- function(cell_index, rep) 60000L + cell_index * 100L + rep

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------
# iForest: local single-threaded override (harness.R's iforest_method does
# not pin nthreads; isotree defaults to all cores). Same hyperparameters as
# the registry entry (ntrees=1000, sample_size=min(256,n)), + nthreads=1.
# This is an addition inside THIS script, not an edit to harness.R.
iforest_1thread_wp6 <- function(X, d, Y = NULL, seed = 1, ...) {
  X <- as.matrix(X)
  set.seed(seed)
  sample_size <- min(256, nrow(X))
  model <- isotree::isolation.forest(X, ntrees = 1000, sample_size = sample_size,
                                     nthreads = 1)
  score <- as.numeric(predict(model, X))
  list(score = score)
}

# R-side method dispatch. MCCD methods get S_min/method/quant per
# WP6_PROTOCOL.md #3; DBSCAN needs Y for its oracle-contamination
# convention (matches harness.R's dbscan_method); iForest uses the local
# single-threaded override above, not METHOD_REGISTRY[["iForest"]].
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
BASELINE_METHODS <- c("LOF", "DBSCAN", "MST", "ODIN", "iForest")
R_METHOD_ORDER <- c(BASELINE_METHODS, MCCD_METHODS)   # cheap-first
# --methods=UN-MCCD,SUN-MCCD restricts the run to a subset (2026-09-23: rerun of
# the NND-based pair after the incremental nearest-neighbour update in
# R/ccds/UN_CCD.R; see tr1/92b_validate_incremental_radi.R).
METHODS_SEL <- opt_val("methods", "")
if (nzchar(METHODS_SEL)) {
  sel <- strsplit(METHODS_SEL, ",", fixed = TRUE)[[1]]
  stopifnot("--methods: unknown method" = all(sel %in% R_METHOD_ORDER))
  R_METHOD_ORDER <- R_METHOD_ORDER[R_METHOD_ORDER %in% sel]
}
S_MIN <- 0.05
# --nn-direction=descend times the descending radius search of UN-MCCD and
# SUN-MCCD (the direction the d >= 10 simulation scripts use); default ascend.
NN_DIRECTION <- opt_val("nn-direction", "ascend")
stopifnot("--nn-direction must be ascend or descend" = NN_DIRECTION %in% c("ascend", "descend"))

call_method <- function(method, X, d, Y, seed) {
  if (method == "U-MCCD")   return(list(res = umccd_method(X = X, d = d, Y = Y)))
  if (method == "SU-MCCD")  return(list(res = sumccd_method(X = X, d = d, Y = Y, min.cls = S_MIN)))
  if (method == "UN-MCCD")  return(list(res = unmccd_method(X = X, d = d, Y = Y, method = NN_DIRECTION)))
  if (method == "SUN-MCCD") return(list(res = sunmccd_method(X = X, d = d, Y = Y, method = NN_DIRECTION, min.cls = S_MIN)))
  if (method == "LOF")      return(list(res = METHOD_REGISTRY[["LOF"]](X = X, d = d, Y = Y)))
  if (method == "DBSCAN")   return(list(res = METHOD_REGISTRY[["DBSCAN"]](X = X, d = d, Y = Y)))
  if (method == "MST")      return(list(res = METHOD_REGISTRY[["MST"]](X = X, d = d, Y = Y, cont = 0.05)))
  if (method == "ODIN")     return(list(res = METHOD_REGISTRY[["ODIN"]](X = X, d = d, Y = Y)))
  if (method == "iForest")  return(list(res = iforest_1thread_wp6(X = X, d = d, Y = Y, seed = seed)))
  stop("call_method(): unknown method: ", method)
}

# Table-load timing for the 4 MCCD methods only (WP6_PROTOCOL.md #4): loaded
# again by this driver, using the same resolvers wp0_mccd_methods.R's
# wrappers use internally, purely to time it. Harmless double load -- both
# occur before either's own timed region.
table_load_seconds <- function(method, d, n) {
  if (!(method %in% MCCD_METHODS)) return(NA_real_)
  t0 <- proc.time()
  if (method %in% c("U-MCCD", "SU-MCCD")) {
    invisible(get_simul("RK", d, quant = rk_quant_label_paper(d), n = n))
  } else if (method == "UN-MCCD") {
    invisible(get_simul("NN", d, quant = nn_quant_label_paper_UN(d), n = n))
  } else if (method == "SUN-MCCD") {
    invisible(get_simul("NN", d, quant = nn_quant_label_paper_SUN(d), n = n))
  }
  as.numeric((proc.time() - t0)[["elapsed"]])
}

# ---------------------------------------------------------------------------
# Memory: gc() high-water mark (WP6_PROTOCOL.md #5)
# ---------------------------------------------------------------------------
mem_peak_mb_of <- function(gc_after) sum(gc_after[, 6])

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
# --resdir=wp6_incr writes to results/tr1/wp6_incr instead (the original wp6
# results stay untouched).
RESDIR      <- here::here("revision_experiments/results/tr1", opt_val("resdir", "wp6"))
DATA_DIR    <- file.path(RESDIR, "data")
SMOKE_DIR   <- file.path(RESDIR, "smoke")
SMOKE_DATA  <- file.path(SMOKE_DIR, "data")

if (MODE == "smoke") {
  RAW_CSV  <- file.path(SMOKE_DIR, "91_wp6_runtime_raw.csv")
  DONE_CSV <- file.path(SMOKE_DIR, "91_wp6_runtime_done.csv")
  EXPORT_DIR <- SMOKE_DATA
} else {
  RAW_CSV  <- file.path(RESDIR, "91_wp6_runtime_raw.csv")
  DONE_CSV <- file.path(RESDIR, "91_wp6_runtime_done.csv")
  EXPORT_DIR <- DATA_DIR
}
dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# B1 (WP6_VERIFICATION.md, blocker): mem_delta_mb added beside mem_peak_mb --
# see the gc_before/gc_after computation below for why mem_peak_mb alone is
# an absolute peak (baseline included), not a delta.
ROW_COLS <- c("method", "n", "d", "rep", "seed", "time_s", "mem_peak_mb",
              "mem_delta_mb", "table_load_s", "threads", "status")
DONE_COLS <- c("method", "n", "d", "rep")

export_path_for <- function(cell) {
  grid_tag <- if (isTRUE(cell$in_n) && isTRUE(cell$in_d)) "n" else if (isTRUE(cell$in_n)) "n" else "d"
  file.path(EXPORT_DIR, sprintf("%s_%s_rep1.csv", grid_tag, cell$cell_value))
}

export_rep1 <- function(cell) {
  paths <- export_path_for(cell)
  # D fix (WP6_VERIFICATION.md): cell 3 (n=500, d=10) is the ONE cell shared
  # by both sweeps (in_n && in_d both TRUE), but export_path_for() always
  # tags a shared cell "n", so only n_500_rep1.csv was ever written and the
  # Python d-sweep never saw a d=10 point. Write the SAME rep-1 dataset a
  # second time under the d-sweep's own naming convention too.
  if (isTRUE(cell$in_n) && isTRUE(cell$in_d)) {
    paths <- c(paths, file.path(EXPORT_DIR, sprintf("d_%d_rep1.csv", cell$d)))
  }
  todo <- paths[!file.exists(paths)]
  if (length(todo) == 0) return(invisible(paths))
  dat <- gen_uniform_wp6(seed = seed_for(cell$cell_index, 1L), n = cell$n, d = cell$d)
  df <- data.frame(dat$X, label = dat$Y)
  for (p in todo) {
    write.csv(df, p, row.names = FALSE)
    cat(sprintf("[export] wrote %s (n=%d, d=%d)\n", basename(p), cell$n, cell$d))
  }
  invisible(paths)
}

# D fix, dry check (WP6_VERIFICATION.md): --check-exports resolves every
# CELLS row (including cell 3's two names) to a file BASENAME, without
# generating data or writing anything, and asserts the result is exactly 9
# distinct names -- catches a naming collision like the original bug without
# needing a full grid run. Exits 0/1; runs no idle-check (nothing computed).
if (MODE == "check_exports") {
  all_paths <- character(0)
  for (k in seq_len(nrow(CELLS))) {
    cell <- CELLS[k, ]
    all_paths <- c(all_paths, export_path_for(cell))
    if (isTRUE(cell$in_n) && isTRUE(cell$in_d)) {
      all_paths <- c(all_paths, file.path(EXPORT_DIR, sprintf("d_%d_rep1.csv", cell$d)))
    }
  }
  names_only <- basename(all_paths)
  cat(sprintf("[check-exports] %d CELLS rows resolve to %d export names (%d unique):\n",
              nrow(CELLS), length(names_only), length(unique(names_only))))
  for (nm in sort(unique(names_only))) cat("   ", nm, "\n")
  if (length(unique(names_only)) != 9L) {
    stop(sprintf("[check-exports] FAIL: expected 9 distinct export names, got %d.",
                 length(unique(names_only))))
  }
  cat("[check-exports] OK: 9 distinct export names confirmed.\n")
  quit(status = 0L, save = "no")
}

# ---------------------------------------------------------------------------
# --summarize
# ---------------------------------------------------------------------------
if (MODE == "summarize") {
  if (!file.exists(RAW_CSV)) stop("No raw CSV at ", RAW_CSV, " -- run the grid first.")
  df <- read.csv(RAW_CSV, stringsAsFactors = FALSE)
  for (col in c("n", "d", "rep", "time_s", "mem_peak_mb", "mem_delta_mb", "table_load_s")) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  ok <- df[df$status == "OK", ]

  # Summariser fix (WP6_VERIFICATION.md): mean_time_s added beside
  # median_time_s -- CV = sd/mean needs mean as its own centre, not median,
  # to be interpretable as a coefficient of variation. mem_delta_mb columns
  # added alongside the legacy mem_peak_mb ones (B1).
  EMPTY_AGG <- data.frame(n = integer(0), d = integer(0), method = character(0),
                          n_reps = integer(0),
                          median_time_s = numeric(0), mean_time_s = numeric(0),
                          cv_time = numeric(0),
                          median_mem_peak_mb = numeric(0), cv_mem = numeric(0),
                          median_mem_delta_mb = numeric(0), cv_mem_delta = numeric(0),
                          median_table_load_s = numeric(0),
                          stringsAsFactors = FALSE)

  summarize_grid <- function(sub, out_path) {
    if (nrow(sub) == 0) {
      write.csv(EMPTY_AGG, out_path, row.names = FALSE)
      cat(sprintf("[summarize] wrote %s (0 rows -- no matching cells in the raw CSV yet)\n", out_path))
      return(EMPTY_AGG)
    }
    groups <- split(seq_len(nrow(sub)), interaction(sub$n, sub$d, sub$method, drop = TRUE))
    agg <- do.call(rbind, lapply(groups, function(idx) {
      g <- sub[idx, ]
      t <- g$time_s[is.finite(g$time_s)]
      m <- g$mem_peak_mb[is.finite(g$mem_peak_mb)]
      md <- g$mem_delta_mb[is.finite(g$mem_delta_mb)]
      tl <- g$table_load_s[is.finite(g$table_load_s)]
      data.frame(n = g$n[1], d = g$d[1], method = g$method[1],
                n_reps = length(t),
                median_time_s = if (length(t)) stats::median(t) else NA_real_,
                mean_time_s = if (length(t)) mean(t) else NA_real_,
                cv_time = if (length(t) > 1) stats::sd(t) / mean(t) else NA_real_,
                median_mem_peak_mb = if (length(m)) stats::median(m) else NA_real_,
                cv_mem = if (length(m) > 1) stats::sd(m) / mean(m) else NA_real_,
                median_mem_delta_mb = if (length(md)) stats::median(md) else NA_real_,
                cv_mem_delta = if (length(md) > 1) stats::sd(md) / mean(md) else NA_real_,
                median_table_load_s = if (length(tl)) stats::median(tl) else NA_real_,
                stringsAsFactors = FALSE)
    }))
    rownames(agg) <- NULL
    agg <- agg[order(agg$n, agg$d, agg$method), ]
    write.csv(agg, out_path, row.names = FALSE)
    cat(sprintf("[summarize] wrote %s (%d rows)\n", out_path, nrow(agg)))
    agg
  }

  agg_n <- summarize_grid(ok[ok$d == 10, ], file.path(RESDIR, "91_wp6_runtime_n.csv"))
  agg_d <- summarize_grid(ok[ok$n == 500, ], file.path(RESDIR, "91_wp6_runtime_d.csv"))

  # ---------------------------------------------------------------------
  # Log-log slope fits (WP6_PROTOCOL.md #7, extended by B2): every rep in
  # sweep A, not medians. One row per (method, metric): the plain time
  # slope for all 4 MCCD methods; a folded time slope (log(time_s/log(n))
  # ~ log(n), compared to 3 directly instead of "3 with a log factor") for
  # the two RK-based methods only, since folding is only meaningful when
  # the stated bound HAS a log n factor to fold out; and a memory slope
  # (log(mem_delta_mb) ~ log(n), compared to the manuscript's O(n^2) space
  # claim) for all 4.
  # ---------------------------------------------------------------------
  MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
  RK_BASED <- c("U-MCCD", "SU-MCCD")   # stated bound carries a log n factor

  fit_loglog <- function(sub, yvar, method, metric, stated_bound) {
    if (nrow(sub) < 3 || length(unique(sub$n)) < 2) {
      return(data.frame(method = method, metric = metric, n_points = nrow(sub),
                        slope = NA_real_, se = NA_real_, r_squared = NA_real_,
                        stated_bound = stated_bound, stringsAsFactors = FALSE))
    }
    fit <- stats::lm(stats::as.formula(sprintf("%s ~ log(n)", yvar)), data = sub)
    cf <- summary(fit)$coefficients
    data.frame(method = method, metric = metric, n_points = nrow(sub),
              slope = cf["log(n)", "Estimate"], se = cf["log(n)", "Std. Error"],
              r_squared = summary(fit)$r.squared, stated_bound = stated_bound,
              stringsAsFactors = FALSE)
  }

  slope_rows <- list()
  for (m in MCCD_METHODS) {
    time_sub <- ok[ok$d == 10 & ok$method == m & is.finite(ok$time_s) & ok$time_s > 0, ]
    time_sub$log_time <- if (nrow(time_sub)) log(time_sub$time_s) else numeric(0)
    stated <- if (m %in% RK_BASED) "n^3*log(n)" else "n^3"
    slope_rows[[length(slope_rows) + 1L]] <-
      fit_loglog(time_sub, "log_time", m, "time", stated)

    if (m %in% RK_BASED) {
      folded <- time_sub[time_sub$n > 1, ]   # log(n)=0 at n=1; not in this grid, defensive only
      folded$log_time_folded <- if (nrow(folded)) log(folded$time_s / log(folded$n)) else numeric(0)
      slope_rows[[length(slope_rows) + 1L]] <-
        fit_loglog(folded, "log_time_folded", m, "time_folded", "n^3")
    }

    mem_sub <- ok[ok$d == 10 & ok$method == m & is.finite(ok$mem_delta_mb) & ok$mem_delta_mb > 0, ]
    mem_sub$log_mem <- if (nrow(mem_sub)) log(mem_sub$mem_delta_mb) else numeric(0)
    slope_rows[[length(slope_rows) + 1L]] <-
      fit_loglog(mem_sub, "log_mem", m, "mem", "n^2")
  }
  slope_df <- do.call(rbind, slope_rows)
  slope_out <- file.path(RESDIR, "91_wp6_slope.csv")
  write.csv(slope_df, slope_out, row.names = FALSE)
  cat(sprintf("[summarize] wrote %s\n", slope_out))
  print(slope_df)

  # ---------------------------------------------------------------------
  # Python raw CSV summary (WP6_VERIFICATION.md summariser item): the 8
  # WP4 competitors, with the PRE-DECLARED collapse rules from
  # WP6_PROTOCOL.md #3 applied before aggregation, not after -- the raw
  # CSV keeps every k / seed so the supplement can show the full range;
  # this summary is the headline table only.
  #   - MutualKNN, SNN: report k=10 only (median over reps); full k in
  #     {5,10,15,20,30} stays in the raw CSV / supplement.
  #   - DIF, LUNAR: median over ALL seed x rep fits pooled together (5
  #     seeds x --reps reps), not per-seed.
  #   - ECOD, COPOD, HDBSCAN, OPTICS: one variant (deterministic fit);
  #     median over reps.
  # ---------------------------------------------------------------------
  py_raw_path <- file.path(RESDIR, "91b_wp6_runtime_py_raw.csv")
  if (file.exists(py_raw_path)) {
    pdf <- read.csv(py_raw_path, stringsAsFactors = FALSE)
    for (col in c("n", "d", "rep", "time_s", "mem_peak_mb")) {
      pdf[[col]] <- suppressWarnings(as.numeric(pdf[[col]]))
    }
    pok <- pdf[pdf$status == "OK", ]
    K10_METHODS <- c("MutualKNN", "SNN")
    POOL_METHODS <- c("DIF", "LUNAR")
    pok_collapsed <- pok[
      (pok$method %in% K10_METHODS & pok$variant == "k10") |
      (pok$method %in% POOL_METHODS) |
      (!(pok$method %in% c(K10_METHODS, POOL_METHODS))),
    ]
    groups <- split(seq_len(nrow(pok_collapsed)),
                    interaction(pok_collapsed$n, pok_collapsed$d, pok_collapsed$method, drop = TRUE))
    py_agg <- do.call(rbind, lapply(groups, function(idx) {
      g <- pok_collapsed[idx, ]
      t <- g$time_s[is.finite(g$time_s)]
      m <- g$mem_peak_mb[is.finite(g$mem_peak_mb)]
      data.frame(n = g$n[1], d = g$d[1], method = g$method[1],
                collapse_rule = if (g$method[1] %in% K10_METHODS) "k=10 only, median over reps"
                                else if (g$method[1] %in% POOL_METHODS) "median over all seed x rep fits"
                                else "single variant, median over reps",
                n_points = length(t),
                median_time_s = if (length(t)) stats::median(t) else NA_real_,
                mean_time_s = if (length(t)) mean(t) else NA_real_,
                cv_time = if (length(t) > 1) stats::sd(t) / mean(t) else NA_real_,
                median_mem_peak_mb = if (length(m)) stats::median(m) else NA_real_,
                cv_mem = if (length(m) > 1) stats::sd(m) / mean(m) else NA_real_,
                stringsAsFactors = FALSE)
    }))
    rownames(py_agg) <- NULL
    py_agg <- py_agg[order(py_agg$n, py_agg$d, py_agg$method), ]
    py_out <- file.path(RESDIR, "91b_wp6_runtime_py_summary.csv")
    write.csv(py_agg, py_out, row.names = FALSE)
    cat(sprintf("[summarize] wrote %s (%d rows, collapse rules applied)\n", py_out, nrow(py_agg)))
  } else {
    cat(sprintf("[summarize] no Python raw CSV at %s yet -- skipping the Python summary.\n", py_raw_path))
  }

  quit(status = 0L, save = "no")
}

# ---------------------------------------------------------------------------
# Real / smoke run
# ---------------------------------------------------------------------------
CELLS_RUN <- if (MODE == "smoke") {
  CELLS[CELLS$cell_index == 1, ]      # n=100, d=10
} else if (GRID_SEL == "n") {
  CELLS[CELLS$in_n, ]
} else if (GRID_SEL == "d") {
  CELLS[CELLS$in_d, ]
} else {
  CELLS
}
REPS_RUN <- if (MODE == "smoke") 1L else REPS

HIST <- if (file.exists(DONE_CSV)) read.csv(DONE_CSV, stringsAsFactors = FALSE) else NULL

is_done <- function(method, n, d, rep) {
  if (is.null(HIST) || nrow(HIST) == 0) return(FALSE)
  isTRUE(any(HIST$method == method & HIST$n == n & HIST$d == d & HIST$rep == rep))
}

mark_done <- function(method, n, d, rep) {
  row <- list(method = method, n = n, d = d, rep = as.integer(rep))
  append_result(DONE_CSV, row)
  new_df <- as.data.frame(row, stringsAsFactors = FALSE)
  HIST <<- if (is.null(HIST)) new_df else rbind(HIST, new_df)
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}

run_start <- Sys.time()
cat(sprintf("---- exporting rep-1 datasets (%d cells) ----\n", nrow(CELLS_RUN)))
for (k in seq_len(nrow(CELLS_RUN))) export_rep1(CELLS_RUN[k, ])

for (k in seq_len(nrow(CELLS_RUN))) {
  cell <- CELLS_RUN[k, ]
  cat(sprintf("\n==== cell %d/%d: n=%d d=%d ====\n", k, nrow(CELLS_RUN), cell$n, cell$d))

  for (method in R_METHOD_ORDER) {
    # A2 (WP6_VERIFICATION.md, blocker): table_load_s is now measured ONCE
    # per (cell, method), not once per rep. With A1's memoizer in place, this
    # one call both times the load (cold, since the previous method's key
    # was just evicted from the cache) AND warms the cache that every rep
    # below will then hit -- so all reps of this (cell, method) share one
    # table_load_s figure instead of each re-triggering its own load.
    tls <- table_load_seconds(method, cell$d, cell$n)
    for (r in seq_len(REPS_RUN)) {
      if (is_done(method, cell$n, cell$d, r)) {
        cat(sprintf("  %-9s rep %d/%d [checkpoint skip]\n", method, r, REPS_RUN))
        next
      }
      seed <- seed_for(cell$cell_index, r)
      dat <- gen_uniform_wp6(seed = seed, n = cell$n, d = cell$d)   # untimed

      # B1 (WP6_VERIFICATION.md, blocker): mem_peak_mb (gc_after's high-water
      # column, col 6) is an ABSOLUTE peak, not a delta -- it includes
      # whatever was already resident before this call, e.g. every RK table
      # carries an unused 763 MB Kest.m matrix that the wrapper's own
      # get_simul() call reads in (kept warm across reps by A1's cache, but
      # allocated before THIS gc(reset=TRUE)). mem_delta_mb below subtracts
      # the pre-call baseline (col 2, "used (Mb)" as of the reset) from the
      # post-call peak, isolating what this one call allocated.
      gc_before <- gc(reset = TRUE)
      t0 <- proc.time()
      out <- tryCatch(
        with_timeout(function() call_method(method, dat$X, cell$d, dat$Y, seed), TIMEOUT_SEC),
        error = function(e) e
      )
      elapsed <- as.numeric((proc.time() - t0)[["elapsed"]])
      gc_after <- gc(reset = FALSE)
      mem_mb <- mem_peak_mb_of(gc_after)
      mem_delta_mb <- sum(gc_after[, 6]) - sum(gc_before[, 2])

      if (inherits(out, "error")) {
        msg <- conditionMessage(out)
        is_timeout <- grepl("elapsed time limit|CPU time limit", msg) || elapsed >= TIMEOUT_SEC - 1
        status <- if (is_timeout) "FLAGGED_TIMEOUT" else
          paste0("ERROR: ", substr(gsub("[\r\n,]+", " ", msg), 1, 160))
        row <- list(method = method, n = cell$n, d = cell$d, rep = as.integer(r), seed = as.integer(seed),
                    time_s = round(elapsed, 4), mem_peak_mb = round(mem_mb, 3),
                    mem_delta_mb = round(mem_delta_mb, 3),
                    table_load_s = if (is.na(tls)) NA_real_ else round(tls, 4),
                    threads = THREADS, status = status)
        append_result(RAW_CSV, row)
        cat(sprintf("  %-9s rep %d/%d %s after %.2fs -- remaining reps of this (cell, method) dropped\n",
                    method, r, REPS_RUN, status, elapsed))
        break
      }

      score_len <- length(out$res$score)
      if (score_len != cell$n) {
        row <- list(method = method, n = cell$n, d = cell$d, rep = as.integer(r), seed = as.integer(seed),
                    time_s = round(elapsed, 4), mem_peak_mb = round(mem_mb, 3),
                    mem_delta_mb = round(mem_delta_mb, 3),
                    table_load_s = if (is.na(tls)) NA_real_ else round(tls, 4),
                    threads = THREADS,
                    status = sprintf("ERROR: score length %d != n %d", score_len, cell$n))
        append_result(RAW_CSV, row)
        cat(sprintf("  %-9s rep %d/%d bad score length -- dropped remaining reps\n", method, r, REPS_RUN))
        break
      }

      row <- list(method = method, n = cell$n, d = cell$d, rep = as.integer(r), seed = as.integer(seed),
                  time_s = round(elapsed, 4), mem_peak_mb = round(mem_mb, 3),
                  mem_delta_mb = round(mem_delta_mb, 3),
                  table_load_s = if (is.na(tls)) NA_real_ else round(tls, 4),
                  threads = THREADS, status = "OK")
      append_result(RAW_CSV, row)
      mark_done(method, cell$n, cell$d, r)
      cat(sprintf("  %-9s rep %d/%d time_s=%.4f mem_peak_mb=%.2f mem_delta_mb=%.2f table_load_s=%s [OK]\n",
                  method, r, REPS_RUN, elapsed, mem_mb, mem_delta_mb,
                  if (is.na(tls)) "NA" else sprintf("%.4f", tls)))
    }
  }
}

cat(sprintf("\n---- done in %.2f min ----\n", as.numeric(difftime(Sys.time(), run_start, units = "mins"))))
if (MODE == "smoke") {
  cat("Smoke run complete. Re-run with --smoke to confirm the checkpoint skip fires.\n")
}
