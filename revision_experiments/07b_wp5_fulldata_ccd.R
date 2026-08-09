# 07b_wp5_fulldata_ccd.R -- FULL-data Musk/Speech UNCCD via parallel radius
# search + spliced exact-n quantile tables.
#
# WHY: 07_wp5_subsample_ccd.R ran n=1000 subsamples because the tables ended
# at n=1000. The spliced tables (01g, FINDINGS "Design upgrade 2026-07-17")
# extend them to n=3062 (Musk) / n=3686 (Speech). Serial construction at
# those n was forecast at 37-65+ h/method; this runner removes that wall by
# (a) parallelizing nnccd.radi's per-point radius loop (deterministic,
# independent across points -- identical arithmetic, distributed), and
# (b) memoizing nnccd_clustering_quantile so OOS and IOS share ONE
# construction per process (the wrappers otherwise rebuild the digraph).
#
# CORRECTNESS: the parallel override reproduces harness.R's clamped
# nnccd.radi VERBATIM per point (ascend + descend + scores-mode zero-radius
# fallback). `validate` mode proves it: radii and OOS scores must be
# IDENTICAL (bitwise) to the serial path on Arrhythmia (452x274).
#
# MODES (PowerShell; Rscript under Bash segfaults):
#   Rscript "revision_experiments/07b_wp5_fulldata_ccd.R" validate [cores]
#   Rscript "revision_experiments/07b_wp5_fulldata_ccd.R" run <musk|speech> [cores]
#
# run appends per-(dataset,method) rows to results/wp5_fulldata_raw.csv,
# caches scores to results/scores_cache/<Dataset>_full_<method>.rds, prints
# "DONE run <dataset>" / "ERROR ..." terminal lines. Threshold + metrics
# identical to 06/07 (REAL_DATA_THRESHOLDS, evaluate()).

suppressPackageStartupMessages({ library(parallel); library(here) })

args  <- commandArgs(trailingOnly = TRUE)
MODE  <- if (length(args) >= 1) args[[1]] else stop("mode required: validate | run")
CORES <- {
  ci <- suppressWarnings(as.integer(args[[length(args)]]))
  if (!is.na(ci)) ci else 20L
}

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

source(here::here("revision_experiments/harness.R"))

SHARED_DIR       <- here::here("revision_experiments/results")       # shared infra: not project-specific
RESULTS_DIR      <- file.path(SHARED_DIR, "tr2")
DATASETS_CSV_DIR <- file.path(SHARED_DIR, "datasets_csv")
CACHE_DIR        <- file.path(SHARED_DIR, "scores_cache")
RAW_CSV          <- file.path(RESULTS_DIR, "wp5_fulldata_raw.csv")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

FULL_SPECS <- list(
  musk   = list(csv = "Musk",   d = 166, nmax = 3062,
                spliced = "NN-test-simul_166d_999%_n3062_spliced.RData"),
  speech = list(csv = "Speech", d = 400, nmax = 3686,
                spliced = "NN-test-simul_400d_999%_n3686_spliced.RData")
)

# ---------------------------------------------------------------------------
# (a) Parallel nnccd.radi override. Arithmetic per point is VERBATIM from
# harness.R's clamped serial version (which itself only added the pmin index
# clamp to the UN_CCD.R original). Serial version kept for validate mode.
# ---------------------------------------------------------------------------
nnccd.radi.serial <- nnccd.radi

PAR_RADI_CORES <- CORES   # read by the override below

nnccd.radi <- function(dx, quantile = "lower", method = "ascend", low.num, quant,
                       simul = NULL, niter, scores = F) {
  ddx <- as.matrix(dist(dx))
  n <- nrow(dx)
  d <- ncol(dx)

  if (quantile != "lower") stop("parallel nnccd.radi: only quantile='lower' supported (as used by all callers)")

  if (!is.null(simul)) {
    avg_len <- length(simul$average)
    med_len <- length(simul$median)
    NN.envelop <- list(
      average = simul$average[pmin(1:n, avg_len)],
      median  = simul$median[pmin(1:n, med_len)]
    )
  } else {
    NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)
  }

  radi_one <- function(i) {
    R_i <- 0
    if (method == "ascend") {
      o.d <- order(ddx[i, ])
      for (j in low.num:n) {
        r <- ddx[i, o.d[j]]
        NN.dist.obs <- NNDest.dist.f(ddx[o.d[2:j], o.d[2:j]], r)
        lower.bound.ave <- NN.envelop$average[j - 1]
        lower.bound.med <- NN.envelop$median[j - 1]
        if (NN.dist.obs$averge < lower.bound.ave | NN.dist.obs$median < lower.bound.med) {
          if (j == low.num) R_i <- 0
          else R_i <- ddx[i, o.d[j - 1]]
          break
        }
      }
    }
    if (method == "descend") {
      o.d <- order(ddx[i, ], decreasing = T)
      for (j in 1:(n - low.num)) {
        r <- ddx[i, o.d[j]]
        NN.dist.obs <- NNDest.dist.f(ddx[o.d[j:(n - 1)], o.d[j:(n - 1)]], r)
        lower.bound.ave <- rev(NN.envelop$average)[j + 2]
        lower.bound.med <- rev(NN.envelop$median)[j + 2]
        if (NN.dist.obs$averge > lower.bound.ave & NN.dist.obs$median > lower.bound.med) {
          R_i <- r
          break
        }
      }
    }
    if (scores && R_i == 0) R_i <- sort(ddx[i, ])[2]   # zero-radius fallback (scores mode only)
    R_i
  }

  # Chunked with per-chunk disk checkpoints: the 2026-07-19 Windows-Update
  # reboot killed a 13-h all-or-nothing search. Chunks make an interruption
  # cost at most one chunk and give live progress in the log. Values are
  # bit-exact vs the unchunked call: radi_one per point is unchanged and
  # parSapplyLB preserves input order within each chunk.
  fp  <- gsub("[^0-9a-zA-Z]", "", sprintf("%.10e", sum(ddx)))  # data fingerprint
  key <- sprintf("radi_n%d_d%d_%s_low%d_sc%d_%s", n, d, method, low.num,
                 as.integer(isTRUE(scores)), fp)
  chunk_dir <- file.path(CACHE_DIR, "radi_chunks")
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  chunks <- split(1:n, ceiling((1:n) / 200))

  cl <- makeCluster(PAR_RADI_CORES)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("ddx", "n", "low.num", "method", "scores", "NN.envelop", "NNDest.dist.f"),
                envir = environment())
  t0 <- Sys.time()
  R <- numeric(n)
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    cf  <- file.path(chunk_dir, sprintf("%s_c%03d.rds", key, ci))
    if (file.exists(cf)) {
      R[idx] <- readRDS(cf)
      cat(sprintf("  [radi] chunk %d/%d reused from checkpoint\n", ci, length(chunks)))
    } else {
      R[idx] <- unname(parSapplyLB(cl, idx, radi_one))  # unname: parSapplyLB attaches
                                                        # names; serial is unnamed and
                                                        # identical() checks attributes
      saveRDS(R[idx], cf)
      cat(sprintf("  [radi] chunk %d/%d done, elapsed %.1f min\n", ci, length(chunks),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
  cat(sprintf("  [parallel radi] n=%d cores=%d wall=%.1f s\n", n, PAR_RADI_CORES,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  return(list(R = R, KS = NULL))
}

# ---------------------------------------------------------------------------
# (b) Memoize nnccd_clustering_quantile: OOS and IOS call it with identical
# arguments; first call computes, second returns the cached graph.
# ---------------------------------------------------------------------------
nnccd_clustering_quantile_orig <- nnccd_clustering_quantile
.graph_cache <- new.env()
nnccd_clustering_quantile <- function(datax, ...) {
  key <- sprintf("%d_%d", nrow(datax), ncol(datax))
  if (!is.null(.graph_cache[[key]])) {
    cat("  [memo] reusing cached construction\n")
    return(.graph_cache[[key]])
  }
  g <- nnccd_clustering_quantile_orig(datax, ...)
  .graph_cache[[key]] <- g
  g
}

# ---------------------------------------------------------------------------
# Shared helpers (identical conventions to 07_wp5_subsample_ccd.R)
# ---------------------------------------------------------------------------
load_csv <- function(name) {
  path <- file.path(DATASETS_CSV_DIR, paste0(name, ".csv"))
  stopifnot("dataset CSV not found" = file.exists(path))
  df <- read.csv(path)
  stopifnot(colnames(df)[ncol(df)] == "label")
  d <- ncol(df) - 1
  list(X = as.matrix(df[, seq_len(d), drop = FALSE]), Y = df[[ncol(df)]],
       d = d, n = nrow(df), n_outliers = sum(df[[ncol(df)]] == 0))
}

score_problem <- function(score, n, method) {
  if (!is.numeric(score)) return("score not numeric")
  if (length(score) != n) return(sprintf("score length %d != n %d", length(score), n))
  if (any(is.na(score))) return(sprintf("%d NA/NaN scores", sum(is.na(score))))
  bad <- is.infinite(score) & !(identical(method, "UNCCD-OOS") & score > 0)
  if (any(bad)) return(sprintf("%d unexpected non-finite scores", sum(bad)))
  NULL
}

append_row <- function(row) append_result(RAW_CSV, row)

# ---------------------------------------------------------------------------
if (MODE == "validate") {
  cat(sprintf("=== 07b validate === Arrhythmia serial-vs-parallel, cores=%d\n", CORES))
  ds <- load_csv("Arrhythmia")
  tab <- get_simul("NN", ds$d)
  dir <- if (ds$d <= 5) "ascend" else "descend"

  cat("serial radii...\n"); t0 <- Sys.time()
  Rs <- nnccd.radi.serial(ds$X, quantile = "lower", method = dir, low.num = 3,
                          quant = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE)
  cat(sprintf("  serial wall %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  cat("parallel radii...\n")
  Rp <- nnccd.radi(ds$X, quantile = "lower", method = dir, low.num = 3,
                   quant = tab$quant, simul = tab$simul, niter = 1000, scores = TRUE)

  ok_r <- identical(Rs$R, unname(Rp$R)) || identical(unname(Rs$R), unname(Rp$R))
  cat(sprintf("[%s] radii identical (n=%d)\n", if (ok_r) "PASS" else "FAIL", ds$n))

  cat("full OOS score via parallel construction...\n")
  sc_p <- NNCCD_OOS(datax = ds$X, simul = tab$simul, method = dir, d = ds$d)
  cat("full OOS score via serial construction...\n")
  .graph_cache <- new.env()                      # drop memo so serial path recomputes
  nnccd.radi <- nnccd.radi.serial                # temporarily restore serial
  sc_s <- NNCCD_OOS(datax = ds$X, simul = tab$simul, method = dir, d = ds$d)
  ok_s <- identical(sc_s, sc_p)
  cat(sprintf("[%s] OOS scores identical\n", if (ok_s) "PASS" else "FAIL"))

  if (ok_r && ok_s) cat("VALIDATE PASS\n") else { cat("VALIDATE FAIL\n"); quit(status = 1) }

} else if (MODE == "run") {
  spec_name <- args[[2]]
  spec <- FULL_SPECS[[spec_name]]
  stopifnot(!is.null(spec))
  cat(sprintf("=== 07b run === %s full data, cores=%d\nStart: %s\n",
              spec_name, CORES, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  e <- new.env(); load(file.path(here::here("R/NN-test_quantile"), spec$spliced), envir = e)
  simul_spliced <- e$simul
  stopifnot(length(simul_spliced$average) == spec$nmax)

  ds <- load_csv(spec$csv)
  stopifnot(ds$n == spec$nmax, ds$d == spec$d)
  cat(sprintf("  loaded: n=%d, d=%d, n_outliers=%d (%.2f%%)\n",
              ds$n, ds$d, ds$n_outliers, 100 * ds$n_outliers / ds$n))
  dir <- "descend"

  for (m in c("UNCCD-OOS", "UNCCD-IOS")) {
    res <- tryCatch({
      t0 <- Sys.time()
      sc <- if (m == "UNCCD-OOS") NNCCD_OOS(datax = ds$X, simul = simul_spliced, method = dir, d = ds$d)
            else                  NNCCD_IOS(datax = ds$X, simul = simul_spliced, method = dir, d = ds$d, min.cls = 0)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      prob <- score_problem(sc, ds$n, m)
      if (!is.null(prob)) stop(prob, call. = FALSE)
      saveRDS(list(score = sc, dataset = paste0(spec$csv, "_full"), method = m,
                   n = ds$n, d = ds$d, t_wall = wall),
              file.path(CACHE_DIR, sprintf("%s_full_%s.rds", spec$csv, m)))
      thr <- REAL_DATA_THRESHOLDS[[m]]
      v <- evaluate(ds$Y, sc, thr)
      append_row(list(dataset = paste0(spec$csv, "_full"), method = m, n = ds$n, d = ds$d,
                      n_outliers = ds$n_outliers, cores = CORES,
                      wall_seconds = round(wall, 2),
                      TPR = round(v[["TPR"]], 4), TNR = round(v[["TNR"]], 4),
                      BA = round(v[["BA"]], 4), F2 = round(v[["F2"]], 4),
                      status = "OK", timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      cat(sprintf("  OK %s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f wall=%.1fs\n",
                  m, v["TPR"], v["TNR"], v["BA"], v["F2"], wall))
      NULL
    }, error = function(err) err)
    if (!is.null(res)) {
      msg <- gsub("[\r\n]+", " ", conditionMessage(res))
      cat(sprintf("ERROR run %s %s: %s\n", spec_name, m, msg))
      append_row(list(dataset = paste0(spec$csv, "_full"), method = m, n = ds$n, d = ds$d,
                      n_outliers = ds$n_outliers, cores = CORES, wall_seconds = NA,
                      TPR = NA, TNR = NA, BA = NA, F2 = NA,
                      status = paste0("ERROR: ", substr(msg, 1, 200)),
                      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    }
  }
  cat(sprintf("DONE run %s\n", spec_name))

} else {
  stop("Unknown mode. Use validate [cores] | run <musk|speech> [cores]")
}
