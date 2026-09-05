# 08b_wp6_metrics.R
#
# WP6 (task T6) metrics pass for the PR-D-26-05767 revision-experiments
# pipeline. Reads the per-(dataset, method, seed) raw PyOD decision scores
# produced by 08_wp6_pyod_baselines.py and computes TPR/TNR/BA/F2 via this repo's own
# R/general_functions/count.R::count_scores2(), for two label thresholds:
#   (a) contamination = 0.1        (PyOD's own default)
#   (b) contamination = true rate  (n_outliers / n, from the dataset itself)
#
# IMPORTANT -- count_scores2's implicit contract:
#   count_scores2(Y, score, threshold) does NOT use Y positionally to decide
#   which rows are "supposed to be" regular vs outlier. It only uses Y to
#   compute n0 <- sum(Y == 0), then ASSUMES the first (n - n0) rows of the
#   vectors it is given are the regular observations and the last n0 rows
#   are the outliers:
#     TNR <- length(which(label_pred[1:(n - n0)] == 1)) / (n - n0)
#     TPR <- length(which(label_pred[(n - n0 + 1):n] == 0)) / n0
#   Every one of this repo's own datasets_csv exports satisfies that
#   block-order convention (regulars first, outliers last) EXCEPT Glass.csv,
#   where the 9 outlier rows sit at positions 182-190 of 213 (not the final
#   9 rows) -- verified directly against the CSV before writing this script.
#   Calling count_scores2 on Glass's native row order would therefore
#   silently score against the wrong rows. To use count_scores2 UNMODIFIED
#   (per task instructions, it must not be touched) while still getting
#   correct numbers on every dataset, we jointly reorder (Y, score) by Y so
#   that all regular rows precede all outlier rows before every
#   count_scores2 call, for every dataset -- a no-op everywhere except
#   Glass, and asserted below to actually produce a block-ordered result.
#
# Threshold derivation ("careful with ties"): the score quantile threshold
# for a target flagged count k is taken as the k-th largest score (order
# statistic), i.e. threshold <- sort(score, decreasing = TRUE)[k]. Because
# count_scores2 flags score >= threshold as outlier, ties AT the threshold
# value can inflate the flagged count above k; we detect and report (via a
# console warning) any case where the actual flagged count exceeds 1.5x the
# target k, but do not otherwise adjust the threshold (a tie-breaking rule
# beyond this is not specified and the inflation is diagnostic, not fatal).
#
# Aggregation: LUNAR and AutoEncoder have 5 seeds each; ECOD has exactly one
# (deterministic) fit, saved under seed slot 0. sd is computed with base R's
# sd(), which returns NA for a vector of length 1 (ECOD) -- left as NA rather
# than forced to 0, since NA correctly signals "no variability was
# estimated" as opposed to 0 which would misleadingly claim exact
# reproducibility was measured across repeated seeds.
#
# Run from anywhere:
#   Rscript "revision_experiments/tr2/08b_wp6_metrics.R"
#
# Output:
#   revision_experiments/results/tr2/wp6_pyod_metrics.csv
#     dataset, method, labeling, TPR_mean, TPR_sd, TNR_mean, TNR_sd,
#     BA_mean, BA_sd, F2_mean, F2_sd, n_seeds

suppressPackageStartupMessages({
  library(dplyr)
})

# ---------------------------------------------------------------------------
# 0. Path resolution (robust to invocation directory; mirrors 00_env_check.R)
# ---------------------------------------------------------------------------
this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) == 1) normalizePath(file_arg) else NA_character_
}, error = function(e) NA_character_)

script_dir <- if (!is.na(this_file)) dirname(this_file) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

count_r_path <- file.path(repo_root, "R", "general_functions", "count.R")
stopifnot("count.R not found" = file.exists(count_r_path))
source(count_r_path)
stopifnot("count_scores2 not defined after sourcing count.R" = exists("count_scores2"))

data_dir <- file.path(repo_root, "revision_experiments", "results", "datasets_csv")
scores_dir <- file.path(repo_root, "revision_experiments", "results", "tr2", "wp6_scores")
out_path <- file.path(repo_root, "revision_experiments", "results", "tr2", "wp6_pyod_metrics.csv")

stopifnot("datasets_csv directory not found" = dir.exists(data_dir))
stopifnot("wp6_scores directory not found (run 08_wp6_pyod_baselines.py first)" = dir.exists(scores_dir))

# ---------------------------------------------------------------------------
# 1. Manifest + dataset enumeration
# ---------------------------------------------------------------------------
manifest <- read.csv(file.path(data_dir, "manifest.csv"), stringsAsFactors = FALSE)

dataset_files <- list.files(data_dir, pattern = "\\.csv$", full.names = FALSE)
dataset_names <- sort(setdiff(sub("\\.csv$", "", dataset_files), "manifest"))
cat(sprintf("Found %d datasets: %s\n", length(dataset_names), paste(dataset_names, collapse = ", ")))

methods <- c("ECOD", "LUNAR", "AutoEncoder")

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------

# k-th largest order statistic as the "quantile that yields the intended
# flagged count k" threshold. score >= threshold flags a point as outlier
# under count_scores2's convention.
kth_largest_threshold <- function(score, k) {
  n <- length(score)
  k <- max(1L, min(n, as.integer(round(k))))
  sorted_desc <- sort(score, decreasing = TRUE)
  threshold <- sorted_desc[k]
  n_flagged <- sum(score >= threshold)
  list(threshold = threshold, k = k, n_flagged = n_flagged)
}

# List available seeds for a (dataset, method) from the score files on disk
# (raw score files only, i.e. NOT the *_labels.csv companions).
available_seeds <- function(dataset, method) {
  pattern <- sprintf("^%s_%s_seed([0-9]+)\\.csv$", dataset, method)
  files <- list.files(scores_dir, pattern = pattern, full.names = FALSE)
  if (length(files) == 0) return(integer(0))
  seeds <- as.integer(sub(pattern, "\\1", files))
  sort(seeds)
}

# ---------------------------------------------------------------------------
# 3. Main loop: per (dataset, method, seed, labeling) -> TPR/TNR/BA/F2
# ---------------------------------------------------------------------------
per_seed_rows <- list()
missing_combos <- character(0)

for (dataset in dataset_names) {
  ds_path <- file.path(data_dir, paste0(dataset, ".csv"))
  ds <- read.csv(ds_path)
  stopifnot("expected final column named 'label'" = colnames(ds)[ncol(ds)] == "label")
  Y <- as.integer(ds$label)  # 1 = regular, 0 = outlier
  n <- length(Y)
  n0_data <- sum(Y == 0)

  manifest_row <- manifest[manifest$dataset == dataset, ]
  if (nrow(manifest_row) == 1 && manifest_row$n_outliers != n0_data) {
    cat(sprintf(
      "WARNING: %s manifest n_outliers=%d != data label==0 count=%d; using data-derived value\n",
      dataset, manifest_row$n_outliers, n0_data
    ))
  }

  # Joint reorder: regulars (Y==1) first, outliers (Y==0) last. Required
  # unconditionally for count_scores2's positional contract (see header);
  # a no-op for every dataset except Glass.
  ord <- order(Y == 0)  # FALSE (Y==1, regular) sorts before TRUE (Y==0, outlier)
  Y_ord <- Y[ord]
  stopifnot(
    "post-reorder block structure violated (regulars-first)" =
      all(Y_ord[seq_len(n - n0_data)] == 1),
    "post-reorder block structure violated (outliers-last)" =
      if (n0_data > 0) all(Y_ord[(n - n0_data + 1):n] == 0) else TRUE
  )

  k_contam01 <- n * 0.1
  k_truerate <- n0_data

  for (method in methods) {
    seeds <- available_seeds(dataset, method)
    if (length(seeds) == 0) {
      missing_combos <- c(missing_combos, sprintf("%s / %s (no seeds found)", dataset, method))
      next
    }

    for (seed in seeds) {
      score_path <- file.path(scores_dir, sprintf("%s_%s_seed%d.csv", dataset, method, seed))
      score_df <- read.csv(score_path)
      stopifnot("score CSV must have a single 'score' column" = "score" %in% colnames(score_df))
      score <- score_df$score
      stopifnot("score length mismatch vs dataset" = length(score) == n)
      score_ord <- score[ord]

      th01 <- kth_largest_threshold(score_ord, k_contam01)
      thtrue <- kth_largest_threshold(score_ord, k_truerate)

      if (th01$n_flagged > 1.5 * th01$k) {
        cat(sprintf(
          "WARNING: tie inflation, %s/%s/seed%d contam0.1: target k=%d, actual flagged=%d\n",
          dataset, method, seed, th01$k, th01$n_flagged
        ))
      }
      if (thtrue$n_flagged > 1.5 * thtrue$k) {
        cat(sprintf(
          "WARNING: tie inflation, %s/%s/seed%d true_rate: target k=%d, actual flagged=%d\n",
          dataset, method, seed, thtrue$k, thtrue$n_flagged
        ))
      }

      res01 <- count_scores2(Y_ord, score_ord, th01$threshold)   # c(TPR, TNR, BA, F2)
      restrue <- count_scores2(Y_ord, score_ord, thtrue$threshold)

      per_seed_rows[[length(per_seed_rows) + 1]] <- data.frame(
        dataset = dataset, method = method, labeling = "contam_0.1", seed = seed,
        TPR = res01[1], TNR = res01[2], BA = res01[3], F2 = res01[4]
      )
      per_seed_rows[[length(per_seed_rows) + 1]] <- data.frame(
        dataset = dataset, method = method, labeling = "true_rate", seed = seed,
        TPR = restrue[1], TNR = restrue[2], BA = restrue[3], F2 = restrue[4]
      )
    }
  }
}

if (length(missing_combos) > 0) {
  cat("\nMissing (dataset, method) combos (no score files found):\n")
  cat(paste(" -", missing_combos, collapse = "\n"), "\n")
}

stopifnot("no per-seed metric rows computed -- check wp6_scores/ contents" = length(per_seed_rows) > 0)
per_seed <- do.call(rbind, per_seed_rows)

# ---------------------------------------------------------------------------
# 4. Aggregate across seeds (mean, sd; sd = NA for n_seeds == 1, see header)
# ---------------------------------------------------------------------------
agg <- per_seed %>%
  group_by(dataset, method, labeling) %>%
  summarise(
    TPR_mean = mean(TPR), TPR_sd = sd(TPR),
    TNR_mean = mean(TNR), TNR_sd = sd(TNR),
    BA_mean = mean(BA), BA_sd = sd(BA),
    F2_mean = mean(F2), F2_sd = sd(F2),
    n_seeds = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(dataset, method, labeling)

write.csv(agg, out_path, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(agg), out_path))

# ---------------------------------------------------------------------------
# 5. Sanity checks
# ---------------------------------------------------------------------------
cat("\n=== Sanity checks ===\n")

metric_cols <- c("TPR_mean", "TNR_mean", "BA_mean", "F2_mean")
out_of_range <- agg[apply(agg[, metric_cols], 1, function(r) any(r < 0 | r > 1 | is.nan(r))), ]
if (nrow(out_of_range) > 0) {
  cat("WARNING: metrics outside [0,1] or NaN found:\n")
  print(out_of_range)
} else {
  cat("All TPR/TNR/BA/F2 means are within [0,1] and non-NaN.\n")
}

na_check <- agg[is.na(agg$TPR_mean) | is.na(agg$TNR_mean) | is.na(agg$BA_mean) | is.na(agg$F2_mean), ]
if (nrow(na_check) > 0) {
  cat("WARNING: NA found in *_mean columns:\n")
  print(na_check)
} else {
  cat("No NA in *_mean columns.\n")
}

musk_row <- agg[agg$dataset == "Musk" & agg$labeling == "true_rate", ]
cat("\nMusk / true_rate BA_mean by method (expect at least one method near ceiling):\n")
print(musk_row[, c("method", "BA_mean", "F2_mean")])

speech_row <- agg[agg$dataset == "Speech" & agg$labeling == "true_rate", ]
cat("\nSpeech / true_rate BA_mean by method (expect uniformly weak performance):\n")
print(speech_row[, c("method", "BA_mean", "F2_mean")])

cat("\nDone.\n")
