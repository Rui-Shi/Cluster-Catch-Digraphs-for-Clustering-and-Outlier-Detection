#!/usr/bin/env Rscript
# revision_experiments/21_regen_smin_grid.R
#
# Agent A / Task 2 (AE.2 parameter-sensitivity): sweeps S_min (the `min.cls`
# argument, a PROPORTION of n) for the two shape-adaptive proposed detectors,
# SU-MCCD and SUN-MCCD, over the grid {0.01, 0.02, 0.03, 0.05, 0.0625, 0.10,
# 0.15, 0.20}, on the five small real datasets (hepatitis, glass, vertebral,
# ecoli, stamps). U-MCCD and UN-MCCD are NOT swept: they use uniform coverage
# and their wrapper functions (umccd_method/unmccd_method in
# wp0_mccd_methods.R) do not even accept a min.cls argument.
#
# Modeled on 13_wp0_gate.R's run_cell() (tryCatch per cell, append_result(),
# has_result() for restart) but this is a NEW file -- 13_wp0_gate.R,
# harness.R and wp0_mccd_methods.R are read-only per REGEN_SPEC.md rule 0.
#
# Output schema is REGEN_SPEC.md section 3's SHARED schema (dataset, method,
# variant, params, n, d, TPR, TNR, BA, F2, published_*, max_abs_diff,
# match_3dp, t_total, status, note, timestamp) -- NOT 13_wp0_gate.R's own
# (min_cls/min_cls_readings/quant_used/...) schema. `variant` carries the
# swept min.cls value as a string, e.g. "smin=0.0625".
#
# Paper-faithful settings (quant label, low.num) are left at each method
# wrapper's own default (resolved from d via rk_quant_label_paper() /
# nn_quant_label_paper_SUN()) -- only min.cls is swept here.
#
# USAGE
#   Rscript revision_experiments/21_regen_smin_grid.R
#
# Resumable: checkpoints to results/regen_smin_grid.csv via has_result(),
# keyed on (dataset, method, variant).

source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

REPO_ROOT   <- here::here()
OUT_CSV     <- file.path(REPO_ROOT, "revision_experiments/results/regen_smin_grid.csv")
TRUTH_CSV   <- file.path(REPO_ROOT, "revision_experiments/published_realdata_truth.csv")
TIMEOUT_SEC <- 600

DATASETS <- c("hepatitis", "glass", "vertebral", "ecoli", "stamps")
METHODS  <- c("SU-MCCD", "SUN-MCCD")
# Named so the label preserves the grid's own printed form (e.g. "0.10", not
# "0.1") -- purely cosmetic, does not affect the numeric value passed in.
SMIN_GRID <- c("0.01" = 0.01, "0.02" = 0.02, "0.03" = 0.03, "0.05" = 0.05,
                "0.0625" = 0.0625, "0.10" = 0.10, "0.15" = 0.15, "0.20" = 0.20)

stopifnot(all(METHODS %in% names(METHOD_REGISTRY)))

cat(sprintf("21_regen_smin_grid.R: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("21_regen_smin_grid.R: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("21_regen_smin_grid.R: smin grid = %s\n", paste(names(SMIN_GRID), collapse = ", ")))
cat(sprintf("21_regen_smin_grid.R: output = %s\n", OUT_CSV))

if (!file.exists(TRUTH_CSV)) stop("Published truth CSV not found: ", TRUTH_CSV)
truth <- read.csv(TRUTH_CSV, stringsAsFactors = FALSE)

get_truth <- function(dataset, method) {
  sub <- truth[truth$dataset == dataset & truth$method == method, ]
  if (nrow(sub) == 0) return(NULL)
  v <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
  for (i in seq_len(nrow(sub))) {
    if (sub$metric[i] %in% names(v)) v[sub$metric[i]] <- sub$value[i]
  }
  v
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}

ROW_COLS <- c("dataset", "method", "variant", "params", "n", "d",
              "TPR", "TNR", "BA", "F2",
              "published_TPR", "published_TNR", "published_BA", "published_F2",
              "max_abs_diff", "match_3dp", "t_total", "status", "note", "timestamp")

na_row <- function(dataset, method, variant, params, status, note) {
  data.frame(
    dataset = dataset, method = method, variant = variant, params = params,
    n = NA_integer_, d = NA_integer_,
    TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
    published_TPR = NA_real_, published_TNR = NA_real_, published_BA = NA_real_, published_F2 = NA_real_,
    max_abs_diff = NA_real_, match_3dp = NA, t_total = NA_real_,
    status = status, note = note,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

run_cell <- function(dataset, method, smin_label, smin) {
  variant <- sprintf("smin=%s", smin_label)
  keys <- c(dataset = dataset, method = method, variant = variant)
  if (has_result(OUT_CSV, keys)) {
    cat(sprintf("[skip, already recorded] %s x %s (%s)\n", dataset, method, variant))
    return(invisible(NULL))
  }

  cell_result <- tryCatch({
    with_timeout(function() {
      dat <- load_real_dataset(dataset)
      X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n

      t0 <- Sys.time()
      res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = smin)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      if (length(res$score) != n || anyNA(res$score)) {
        stop(sprintf("score sanity check failed: length=%d (n=%d), any NA=%s",
                      length(res$score), n, anyNA(res$score)))
      }

      thr <- REAL_DATA_THRESHOLDS[[method]]
      m <- evaluate(Y, res$score, thr)

      params <- sprintf("min.cls=%s (proportion of n; round(min.cls*n)=%d), quant=%s (%s), low.num=%s",
                         smin_label, round(smin * n),
                         if (!is.null(res$quant_label)) res$quant_label else NA_character_,
                         if (!is.null(res$quant_source)) res$quant_source else NA_character_,
                         if (method == "SU-MCCD") "2 (default)" else "3 (default)")

      base <- list(dataset = dataset, method = method, variant = variant, params = params,
                   n = n, d = d,
                   TPR = unname(m["TPR"]), TNR = unname(m["TNR"]), BA = unname(m["BA"]), F2 = unname(m["F2"]),
                   t_total = wall,
                   timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

      pub <- get_truth(dataset, method)
      if (is.null(pub)) {
        row <- as.data.frame(c(base, list(
          published_TPR = NA_real_, published_TNR = NA_real_, published_BA = NA_real_, published_F2 = NA_real_,
          max_abs_diff = NA_real_, match_3dp = NA, status = "no_published_truth", note = NA_character_
        )), stringsAsFactors = FALSE)[ROW_COLS]
        cat(sprintf("  %-10s x %-9s %-12s n=%-4d d=%-3d TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | NO PUBLISHED TRUTH | wall=%.2fs\n",
                    dataset, method, variant, n, d, m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
        return(row)
      }

      diffs <- abs(c(m["TPR"], m["TNR"], m["BA"], m["F2"]) - pub[c("TPR", "TNR", "BA", "F2")])
      max_diff <- max(diffs)
      match_3dp <- all(round(c(m["TPR"], m["TNR"], m["BA"], m["F2"]), 3) == round(pub[c("TPR", "TNR", "BA", "F2")], 3))

      row <- as.data.frame(c(base, list(
        published_TPR = unname(pub["TPR"]), published_TNR = unname(pub["TNR"]),
        published_BA = unname(pub["BA"]), published_F2 = unname(pub["F2"]),
        max_abs_diff = unname(max_diff), match_3dp = match_3dp, status = "ok", note = NA_character_
      )), stringsAsFactors = FALSE)[ROW_COLS]

      cat(sprintf("  %-10s x %-9s %-12s n=%-4d d=%-3d TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | wall=%.2fs\n",
                  dataset, method, variant, n, d, m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
      row
    }, TIMEOUT_SEC)
  }, error = function(e) {
    is_timeout <- grepl("elapsed time limit|reached CPU time limit", conditionMessage(e))
    status_val <- if (is_timeout) sprintf("timeout(>%.0fs)", TIMEOUT_SEC) else "error"
    cat(sprintf("  %-10s x %-9s %-12s %s: %s\n", dataset, method, variant,
                if (is_timeout) "TIMEOUT" else "ERROR", conditionMessage(e)))
    na_row(dataset, method, variant,
           params = sprintf("min.cls=%s (proportion of n)", smin_label),
           status = status_val, note = substr(conditionMessage(e), 1, 200))
  })

  stopifnot(identical(names(cell_result), ROW_COLS))
  append_result(OUT_CSV, cell_result)
  invisible(cell_result)
}

for (dataset in DATASETS) {
  for (method in METHODS) {
    for (smin_label in names(SMIN_GRID)) {
      run_cell(dataset, method, smin_label, SMIN_GRID[[smin_label]])
    }
  }
}

cat("\n21_regen_smin_grid.R: done.\n")
if (file.exists(OUT_CSV)) {
  final <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  cat(sprintf("%s now has %d rows.\n", OUT_CSV, nrow(final)))
  n_ok <- sum(final$status == "ok", na.rm = TRUE)
  cat(sprintf("status=ok: %d/%d\n", n_ok, nrow(final)))
}
