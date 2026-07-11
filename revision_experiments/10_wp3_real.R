#!/usr/bin/env Rscript
# revision_experiments/10_wp3_real.R
#
# WP3 threshold-sensitivity study, real-data part (Task T7 Phase A).
# For each dataset x each of the 4 CCD-based OS methods: compute the score
# vector ONCE via the harness METHOD_REGISTRY (the exact configuration whose
# WBC row reproduced the manuscript table in 03_smoke.R: UNCCD
# method="ascend", min.cls=0, RK/NN quantile tables resolved by get_simul),
# cache it to results/scores_cache/<dataset>_<method>.rds, then sweep an
# absolute cutoff grid and write results/wp3_sensitivity_real.csv.
#
# Datasets: WBC (n=223, d=9), Thyroid (n=3656, d=6) now; Musk (d=166) is
# auto-included once results/scores_cache/musk_*.rds exist (produced by a
# later task) and skipped with a log line until then.
#
# Cutoff grid: seq(1, 4, length.out = 17) -- spanning 0.5x to 2x of the
# calibrated real-data threshold of 2 (manuscript Section 4). NOTE:
# seq(1, 4, length.out = 17) has step 0.1875 and does NOT contain 2.0
# (neighbors 1.9375 and 2.125), so 2.0 is added explicitly -> 18 cutoffs.
#
# The scoring pipeline is deterministic given the dataset and the
# precomputed quantile tables (no RNG in RKCCD/UNCCD score construction),
# so a single scoring pass per dataset x method suffices.
#
# Rerun-safe: expensive scores come from the RDS cache; the sensitivity CSV
# is rewritten in full on every run (cheap).
#
# Run:  Rscript "revision_experiments/10_wp3_real.R"

suppressMessages(source(here::here("revision_experiments/harness.R")))

REPO_ROOT  <- here::here()
CACHE_DIR  <- file.path(REPO_ROOT, "revision_experiments/results/scores_cache")
OUT_CSV    <- file.path(REPO_ROOT, "revision_experiments/results/wp3_sensitivity_real.csv")
FULL_GRID  <- sort(unique(c(seq(1, 4, length.out = 17), 2)))  # 18 cutoffs incl. calibrated 2.0
CALIBRATED <- 2  # manuscript: "The four OSs use a threshold of 2" (Section 4)
OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# dataset key -> loader. WBC and Thyroid load through the harness
# (RealData_Collection.R objects "WBC" and "thyroid" -- same objects the
# smoke test validated). Musk is score-cache-only: this script never
# computes Musk scores itself (d=166 tables + scoring are a later task);
# it only sweeps cutoffs over cached vectors when they appear.
DATASETS <- list(
  wbc     = list(label = "WBC",     robj = "WBC",     compute = TRUE),
  thyroid = list(label = "Thyroid", robj = "thyroid", compute = TRUE),
  musk    = list(label = "Musk",    robj = NA,        compute = FALSE)
)

cache_path <- function(ds_key, method) file.path(CACHE_DIR, sprintf("%s_%s.rds", ds_key, method))

get_scores <- function(ds_key, ds_cfg) {
  # Returns list(per-method list(score, t_construct, t_total), Y, n, d)
  # or NULL if the dataset must be skipped.
  paths <- vapply(OS_METHODS, function(m) cache_path(ds_key, m), character(1))
  cached <- file.exists(paths)

  if (!ds_cfg$compute && !all(cached)) {
    cat(sprintf("[%s] scores_cache/%s_*.rds not (fully) present (%d/4 cached) -- skipping; rerun after the scoring task lands.\n",
                ds_cfg$label, ds_key, sum(cached)))
    return(NULL)
  }

  need_data <- ds_cfg$compute && !all(cached)
  X <- NULL; Y <- NULL; d <- NA_integer_; n <- NA_integer_
  if (ds_cfg$compute) {
    dat <- load_real_dataset(ds_cfg$robj)
    X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
    cat(sprintf("[%s] loaded: n=%d, d=%d, outliers=%d (%.1f%%)\n",
                ds_cfg$label, n, d, sum(Y == 0), 100 * mean(Y == 0)))
  }

  out <- list()
  for (m in OS_METHODS) {
    p <- cache_path(ds_key, m)
    if (file.exists(p)) {
      out[[m]] <- readRDS(p)
      cat(sprintf("[%s] %-10s loaded from cache (%s)\n", ds_cfg$label, m, basename(p)))
    } else {
      t0 <- Sys.time()
      res <- METHOD_REGISTRY[[m]](X = X, d = d, Y = Y)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rec <- list(score = res$score, t_construct = res$t_construct,
                  t_total = res$t_total, t_wall = wall,
                  dataset = ds_cfg$label, method = m, n = n, d = d)
      saveRDS(rec, p)
      out[[m]] <- rec
      cat(sprintf("[%s] %-10s scored: t_construct=%.1fs t_total=%.1fs wall=%.1fs -> cached %s\n",
                  ds_cfg$label, m, res$t_construct, res$t_total, wall, basename(p)))
    }
  }

  # Y for the sweep: from the loaded dataset when available, otherwise from
  # the cache records (cached for Musk by the later scoring task; all our
  # own records store n but not Y -- for cache-only datasets we fall back to
  # the CSV export, which uses the same row ordering).
  if (is.null(Y)) {
    csv_path <- file.path(REPO_ROOT, "revision_experiments/results/datasets_csv",
                          paste0(ds_cfg$label, ".csv"))
    stopifnot(file.exists(csv_path))
    df <- read.csv(csv_path)
    Y <- df$label
    n <- nrow(df); d <- ncol(df) - 1
  }
  list(scores = out, Y = Y, n = n, d = d)
}

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
all_rows <- list()
timing_rows <- list()

for (ds_key in names(DATASETS)) {
  ds_cfg <- DATASETS[[ds_key]]
  got <- get_scores(ds_key, ds_cfg)
  if (is.null(got)) next

  for (m in OS_METHODS) {
    rec <- got$scores[[m]]
    stopifnot(length(rec$score) == got$n, !any(is.na(rec$score)))
    for (cutoff in FULL_GRID) {
      v <- evaluate(got$Y, rec$score, cutoff)
      all_rows[[length(all_rows) + 1]] <- data.frame(
        dataset = ds_cfg$label, method = m, n = got$n, d = got$d,
        cutoff = cutoff, calibrated = CALIBRATED,
        TPR = unname(v["TPR"]), TNR = unname(v["TNR"]),
        BA = unname(v["BA"]), F2 = unname(v["F2"]),
        stringsAsFactors = FALSE
      )
    }
    timing_rows[[length(timing_rows) + 1]] <- data.frame(
      dataset = ds_cfg$label, method = m,
      t_construct = round(rec$t_construct, 2), t_total = round(rec$t_total, 2),
      stringsAsFactors = FALSE
    )
  }
}

res <- do.call(rbind, all_rows)
write.csv(res, OUT_CSV, row.names = FALSE)
cat(sprintf("\nWrote %d rows (%d dataset(s) x 4 methods x %d cutoffs) -> %s\n",
            nrow(res), length(unique(res$dataset)), length(FULL_GRID), OUT_CSV))

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
stopifnot(!anyNA(res[, c("TPR", "TNR", "BA", "F2")]),
          all(res$TPR >= 0 & res$TPR <= 1), all(res$TNR >= 0 & res$TNR <= 1),
          all(res$BA >= 0 & res$BA <= 1), all(res$F2 >= 0 & res$F2 <= 1))
cat("All metrics in [0,1], no NA.\n")

cat("\nScoring timings (s):\n")
print(do.call(rbind, timing_rows), row.names = FALSE)

# Reproduction check at the calibrated cutoff (2.0) against the manuscript's
# published WBC values (Tables Real_Data_Result_OS1/OS2; identical numbers
# reproduced by 03_smoke.R Part C).
PUBLISHED_WBC <- list(
  "RKCCD-OOS" = c(TPR = 0.200, TNR = 0.850, BA = 0.525, F2 = 0.135),
  "RKCCD-IOS" = c(TPR = 1.000, TNR = 0.779, BA = 0.890, F2 = 0.515),
  "UNCCD-OOS" = c(TPR = 0.200, TNR = 0.897, BA = 0.548, F2 = 0.156),
  "UNCCD-IOS" = c(TPR = 1.000, TNR = 0.807, BA = 0.904, F2 = 0.549)
)
cat("\nWBC @ cutoff = 2.0 vs manuscript:\n")
wbc2 <- res[res$dataset == "WBC" & abs(res$cutoff - 2) < 1e-9, ]
gate_ok <- TRUE
for (m in OS_METHODS) {
  r <- wbc2[wbc2$method == m, ]
  pub <- PUBLISHED_WBC[[m]]
  diffs <- abs(c(r$TPR, r$TNR, r$BA, r$F2) - pub)
  ok <- all(diffs <= 0.005 + 1e-9)  # manuscript rounds to 3 decimals
  if (!ok) gate_ok <- FALSE
  cat(sprintf("  %-10s got TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | pub %.3f/%.3f/%.3f/%.3f | %s\n",
              m, r$TPR, r$TNR, r$BA, r$F2, pub["TPR"], pub["TNR"], pub["BA"], pub["F2"],
              if (ok) "MATCH" else "MISMATCH"))
}
cat(sprintf("WBC reproduction gate at cutoff=2: %s\n", if (gate_ok) "PASS" else "FAIL"))
if (!gate_ok) stop("WBC values at cutoff=2 do not match the manuscript -- investigate before using wp3_sensitivity_real.csv")

cat("\n10_wp3_real.R done.\n")
