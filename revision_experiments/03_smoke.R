#!/usr/bin/env Rscript
# revision_experiments/03_smoke.R
#
# Smoke test + reproduction gate for the 9-method registry defined in
# revision_experiments/harness.R (Task T2).
#
#   Part A: all 9 methods on WBC (real data)              -> smoke_wbc.csv
#   Part B: all 9 methods on synthetic Gaussian d=20 data  -> smoke_synthetic.csv
#   Part C: REPRODUCTION GATE -- compare Part-A WBC numbers against the
#           published manuscript table.
#
# Run with:  Rscript "revision_experiments/03_smoke.R"

source(here::here("revision_experiments/harness.R"))

cat("\n================ 03_smoke.R : START ================\n")

RESULTS_DIR   <- here::here("revision_experiments/results/tr2")
WBC_CSV       <- file.path(RESULTS_DIR, "smoke_wbc.csv")
SYNTH_CSV     <- file.path(RESULTS_DIR, "smoke_synthetic.csv")
METHOD_NAMES  <- names(METHOD_REGISTRY)

# ---------------------------------------------------------------------------
# Sanity-check helpers
# ---------------------------------------------------------------------------

assert_scores_sane <- function(score, n, method_name, allow_pos_inf = FALSE) {
  if (!is.numeric(score)) stop(sprintf("[%s] score is not numeric", method_name))
  if (length(score) != n) stop(sprintf("[%s] score length %d != n %d", method_name, length(score), n))
  if (any(is.na(score))) stop(sprintf("[%s] score contains NA/NaN", method_name))  # is.na() also catches NaN
  bad_inf <- is.infinite(score) & !(allow_pos_inf & score > 0)
  if (any(bad_inf)) {
    stop(sprintf("[%s] score contains an unexpected non-finite value (-Inf always disallowed; +Inf only allowed for OOS methods)", method_name))
  }
  invisible(TRUE)
}

assert_metrics_sane <- function(v, method_name) {
  if (any(is.na(v))) stop(sprintf("[%s] metrics contain NA: %s", method_name, paste(v, collapse = ",")))
  if (any(v < -1e-9 | v > 1 + 1e-9)) {
    stop(sprintf("[%s] metrics out of [0,1]: %s", method_name, paste(round(v, 4), collapse = ",")))
  }
  invisible(TRUE)
}

is_oos_method <- function(method_name) method_name %in% c("RKCCD-OOS", "UNCCD-OOS")

# ---------------------------------------------------------------------------
# Runner: apply every registry method to (X, Y, d) and checkpoint to csv_path
# ---------------------------------------------------------------------------
run_all_methods <- function(dataset_label, X, Y, d, thresholds, csv_path, seed_label,
                             method_args = list()) {
  n <- nrow(X)
  out <- list()
  for (m in METHOD_NAMES) {
    keys <- c(dataset = dataset_label, method = m, seed = seed_label)
    if (has_result(csv_path, keys)) {
      cat(sprintf("  [skip, already recorded] %-10s x %s\n", m, dataset_label))
      # Reconstruct metrics (but not the raw score, which isn't persisted) from
      # the checkpoint so downstream consumers (e.g. the Part-C reproduction
      # gate) still see a result after a restart, not just a fresh run.
      prev <- read.csv(csv_path, stringsAsFactors = FALSE)
      match_row <- prev$dataset == keys["dataset"] & prev$method == keys["method"] &
        as.character(prev$seed) == keys["seed"]
      r <- prev[which(match_row)[1], ]
      out[[m]] <- list(
        score = NULL,
        metrics = c(TPR = r$TPR, TNR = r$TNR, BA = r$BA, F2 = r$F2),
        t_construct = r$t_construct, t_total = r$t_total
      )
      next
    }
    extra_args <- method_args[[m]]
    if (is.null(extra_args)) extra_args <- list()

    t_wall0 <- Sys.time()
    res <- tryCatch(
      do.call(METHOD_REGISTRY[[m]], c(list(X = X, d = d, Y = Y), extra_args)),
      error = function(e) {
        cat(sprintf("  [ERROR] %s on %s: %s\n", m, dataset_label, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) next

    assert_scores_sane(res$score, n, m, allow_pos_inf = is_oos_method(m))

    thresh <- thresholds[[m]]
    metrics <- evaluate(Y, res$score, thresh)
    assert_metrics_sane(metrics, m)

    out[[m]] <- list(score = res$score, metrics = metrics,
                      t_construct = res$t_construct, t_total = res$t_total)

    row <- c(
      dataset = dataset_label, method = m, seed = seed_label,
      n = n, d = d, threshold = thresh,
      TPR = unname(metrics["TPR"]), TNR = unname(metrics["TNR"]),
      BA = unname(metrics["BA"]), F2 = unname(metrics["F2"]),
      t_construct = ifelse(is.na(res$t_construct), NA, round(res$t_construct, 4)),
      t_total = round(res$t_total, 4)
    )
    append_result(csv_path, row)

    cat(sprintf("  %-10s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f  t_construct=%s t_total=%.3fs\n",
                m, metrics["TPR"], metrics["TNR"], metrics["BA"], metrics["F2"],
                ifelse(is.na(res$t_construct), "NA", sprintf("%.3fs", res$t_construct)),
                res$t_total))
  }
  out
}

# ---------------------------------------------------------------------------
# PART A: WBC (real data)
# ---------------------------------------------------------------------------
cat("\n---- Part A: WBC (real data), all 9 registry methods ----\n")

wbc <- load_real_dataset("WBC")
cat(sprintf("WBC: n=%d, d=%d, #outliers=%d (%.1f%%)\n",
            wbc$n, wbc$d, sum(wbc$Y == 0), 100 * mean(wbc$Y == 0)))
# Sanity check against the manuscript's dataset table (Table \ref{tab:Real_Data_OS},
# CCDwScores.tex ~line 1111): WBC n=223, d=9, 10 outliers (4.5%).
stopifnot(wbc$n == 223, wbc$d == 9, sum(wbc$Y == 0) == 10)

partA <- run_all_methods(
  dataset_label = "WBC", X = wbc$X, Y = wbc$Y, d = wbc$d,
  thresholds = REAL_DATA_THRESHOLDS, csv_path = WBC_CSV, seed_label = "1"
)

# ---------------------------------------------------------------------------
# PART A2: Glass ordering-guard sanity check (LOF only)
# ---------------------------------------------------------------------------
#
# RealData_Collection.R's glass block sorts by glass[,9] (a feature, Fe)
# instead of glass[,10] (the label), so -- unlike every other dataset --
# glass's outliers are NOT positionally last. evaluate() now guards against
# this by jointly reordering (Y, score) regulars-first before count_scores2.
# This section runs one cheap method (LOF, Thresh = 1.5 per Real_Data_LOF.R's
# glass block) through the guarded path to pin a known-good Glass reference.
cat("\n---- Part A2: Glass ordering-guard sanity (LOF) ----\n")

glass_ds <- load_real_dataset("glass")
outlier_pos <- which(glass_ds$Y == 0)
cat(sprintf("glass: n=%d, d=%d, #outliers=%d; outlier row positions: %s\n",
            glass_ds$n, glass_ds$d, length(outlier_pos),
            paste(range(outlier_pos), collapse = "-")))
if (max(outlier_pos) == glass_ds$n && min(outlier_pos) == glass_ds$n - length(outlier_pos) + 1) {
  cat("glass outliers ARE positionally last (guard is a no-op here)\n")
} else {
  cat("glass outliers are NOT positionally last -- evaluate()'s reorder guard is load-bearing\n")
}

GLASS_CSV <- file.path(RESULTS_DIR, "smoke_glass.csv")
glass_keys <- c(dataset = "glass", method = "LOF", seed = "1")
if (has_result(GLASS_CSV, glass_keys)) {
  cat("  [skip, already recorded] LOF x glass\n")
  prev <- read.csv(GLASS_CSV, stringsAsFactors = FALSE)
  r <- prev[prev$dataset == "glass" & prev$method == "LOF" & as.character(prev$seed) == "1", ][1, ]
  glass_lof_metrics <- c(TPR = r$TPR, TNR = r$TNR, BA = r$BA, F2 = r$F2)
} else {
  res_glass <- lof_method(glass_ds$X, glass_ds$d, glass_ds$Y)
  assert_scores_sane(res_glass$score, glass_ds$n, "LOF(glass)")
  glass_lof_metrics <- evaluate(glass_ds$Y, res_glass$score, 1.5)
  assert_metrics_sane(glass_lof_metrics, "LOF(glass)")
  append_result(GLASS_CSV, c(
    dataset = "glass", method = "LOF", seed = "1",
    n = glass_ds$n, d = glass_ds$d, threshold = 1.5,
    TPR = unname(glass_lof_metrics["TPR"]), TNR = unname(glass_lof_metrics["TNR"]),
    BA = unname(glass_lof_metrics["BA"]), F2 = unname(glass_lof_metrics["F2"]),
    t_construct = NA, t_total = round(res_glass$t_total, 4)
  ))
}
cat(sprintf("  LOF x glass (guarded): TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f\n",
            glass_lof_metrics["TPR"], glass_lof_metrics["TNR"],
            glass_lof_metrics["BA"], glass_lof_metrics["F2"]))

# ---------------------------------------------------------------------------
# PART B: synthetic Gaussian, d=20, n=250, 2 clusters + 5% contamination
# ---------------------------------------------------------------------------
cat("\n---- Part B: synthetic Gaussian d=20, n=250, 2cls, 5% contamination ----\n")

# Generator idiom copied from
# simulations/outlyingness_scores/RKCCD_OOS_IOS/Simulation/Gaussian/20d/20d_2cls_n200_cont5%.R
# (and the UNCCD twin), adapted from n=200 to n=250.
gen_gaussian_2cls <- function(n, d, cont, cls_dis = 3, otl_dis = 2,
                               r_min = 0.7, r_max = 1.3, noise_level = 0.01, seed = 123) {
  n1 <- round(n * (1 - cont) * 0.5)
  n2 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)

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
    temp <- rpoisball.unit(1, d) * 5 + mu
    r1 <- sqrt(sum((temp - mu1)^2))
    r2 <- sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) {
      outlier <- rbind(outlier, temp)
      i <- i + 1
    }
  }
  rownames(outlier) <- NULL

  X <- rbind(data1, data2, outlier)
  Y <- c(rep(1, n1 + n2), rep(0, n0))
  list(X = as.data.frame(X), Y = Y, n1 = n1, n2 = n2, n0 = n0)
}

synth <- gen_gaussian_2cls(n = 250, d = 20, cont = 0.05, seed = 123)
cat(sprintf("synthetic: n=%d (n1=%d, n2=%d, n0=%d), d=20\n",
            nrow(synth$X), synth$n1, synth$n2, synth$n0))

# Threshold_OOS/Threshold_IOS index 5 == d=20 (vectors are ordered for
# d = 2,3,5,10,20,50,100), per
# simulations/outlyingness_scores/{RKCCD,UNCCD}_OOS_IOS/Simulation/Gaussian/Threshold.R
RK_THRESH_OOS_20D <- 3.5; RK_THRESH_IOS_20D <- 2.5   # RKCCD_OOS_IOS/Simulation/Gaussian/Threshold.R
NN_THRESH_OOS_20D <- 3.0; NN_THRESH_IOS_20D <- 6.5   # UNCCD_OOS_IOS/Simulation/Gaussian/Threshold.R

synthetic_thresholds <- list(
  "RKCCD-OOS" = RK_THRESH_OOS_20D, "RKCCD-IOS" = RK_THRESH_IOS_20D,
  "UNCCD-OOS" = NN_THRESH_OOS_20D, "UNCCD-IOS" = NN_THRESH_IOS_20D,
  "LOF" = 1.5, "DBSCAN" = 0.5, "MST" = 0.5, "ODIN" = 0.5, "iForest" = 0.55
)

# method="descend" for UN-CCD at d=20 matches
# UNCCD_OOS_IOS/Simulation/Gaussian/20d/20d_2cls_n200_cont5%.R exactly.
synthetic_method_args <- list(
  "UNCCD-OOS" = list(method = "descend"),
  "UNCCD-IOS" = list(method = "descend"),
  "MST" = list(cont = 0.05)  # match the dataset's actual 5% contamination level
)

partB <- run_all_methods(
  dataset_label = "synthetic_gaussian_20d_n250_cont5", X = synth$X, Y = synth$Y, d = 20,
  thresholds = synthetic_thresholds, csv_path = SYNTH_CSV, seed_label = "123",
  method_args = synthetic_method_args
)

# ---------------------------------------------------------------------------
# PART C: REPRODUCTION GATE (WBC vs. published manuscript numbers)
# ---------------------------------------------------------------------------
cat("\n---- Part C: REPRODUCTION GATE (WBC vs. published manuscript table) ----\n")

# Published WBC column, extracted read-only from
#   G:\Submissions\TR2\Pattern Recognition (Elsevier) - Resubmit\CCDwScores.tex
# (sibling of this CLONE; the Box path named in the task instructions does not
# exist in this environment -- this is the equivalent manuscript file, same
# relative project layout).
# Table \ref{Real_Data_Result_OS1} (TPR, TNR), lines 1127-1147, WBC column.
# Table \ref{Real_Data_Result_OS2} (BA, F2), lines 1149-1169, WBC column.
PUBLISHED_WBC <- list(
  "RKCCD-OOS" = c(TPR = 0.200, TNR = 0.850, BA = 0.525, F2 = 0.135),
  "RKCCD-IOS" = c(TPR = 1.000, TNR = 0.779, BA = 0.890, F2 = 0.515),
  "UNCCD-OOS" = c(TPR = 0.200, TNR = 0.897, BA = 0.548, F2 = 0.156),
  "UNCCD-IOS" = c(TPR = 1.000, TNR = 0.807, BA = 0.904, F2 = 0.549),
  "LOF"       = c(TPR = 1.000, TNR = 0.793, BA = 0.897, F2 = 0.532),
  "DBSCAN"    = c(TPR = 0.600, TNR = 1.000, BA = 0.800, F2 = 0.652),
  "MST"       = c(TPR = 0.700, TNR = 0.756, BA = 0.728, F2 = 0.354),
  "ODIN"      = c(TPR = 0.500, TNR = 0.869, BA = 0.684, F2 = 0.342),
  "iForest"   = c(TPR = 0.800, TNR = 0.939, BA = 0.869, F2 = 0.656)
)

TOL_PASS <- 0.10
TOL_IDEAL <- 0.05

cat(sprintf("\n%-10s %-6s %10s %10s %10s %8s\n", "method", "metric", "published", "reproduced", "abs_diff", "status"))
cat(strrep("-", 60), "\n")

all_pass <- TRUE
for (m in METHOD_NAMES) {
  pub <- PUBLISHED_WBC[[m]]
  rep_metrics <- if (!is.null(partA[[m]])) partA[[m]]$metrics else NULL
  if (is.null(rep_metrics)) {
    cat(sprintf("%-10s %-6s %10s %10s %10s %8s\n", m, "*", "-", "MISSING", "-", "FAIL"))
    all_pass <- FALSE
    next
  }
  for (metric in c("TPR", "TNR", "BA", "F2")) {
    diff <- abs(pub[metric] - rep_metrics[metric])
    status <- if (diff <= TOL_PASS) "PASS" else "FAIL"
    if (status == "FAIL") all_pass <- FALSE
    flag <- if (diff <= TOL_IDEAL) "" else if (diff <= TOL_PASS) " (drift)" else ""
    cat(sprintf("%-10s %-6s %10.3f %10.3f %10.3f %8s%s\n",
                m, metric, pub[metric], rep_metrics[metric], diff, status, flag))
  }
}

cat(strrep("-", 60), "\n")
cat(sprintf("\nOVERALL REPRODUCTION GATE: %s\n", if (all_pass) "PASS" else "FAIL"))

cat("\n================ 03_smoke.R : DONE ================\n")
