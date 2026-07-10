# 02_load_data.R
#
# Task T1 (high-dimensional dataset acquisition) for the PR-D-26-05767
# revision-experiments pipeline.
#
# Loads the 10 existing real-data outlier-detection benchmarks (via
# RealData_Collection.R) plus 4 new high-dimensional benchmarks (Musk,
# Arrhythmia, Speech, InternetAds), applies the SAME label-polarity and
# row-ordering conventions used by RealData_Collection.R, and exports all
# 14 datasets to CSV under revision_experiments/results/datasets_csv/,
# along with a manifest.csv summarizing n, d, n_outliers, source_url, and
# normalization for each.
#
# Convention (matches RealData_Collection.R): final column = label,
# 1 = regular observation, 0 = outlier. Row order: regular rows first,
# outliers last (decreasing sort on the label column).
#
# Run from the CLONE root:
#   Rscript "revision_experiments/02_load_data.R"
#
# Idempotent: re-running overwrites the same CSV outputs; no state persists
# between runs other than those files.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(foreign)
  library(R.matlab)
})

repo_root <- here::here()
cat("Repo root (here::here()):", repo_root, "\n")

# ---------------------------------------------------------------------------
# 1. Load the 10 existing datasets via the canonical loader.
#    NOTE: RealData_Collection.R calls setwd(here::here("data/outlier_detection"))
#    internally and does not restore the working directory afterward. We
#    explicitly reset it below. here::here() itself is unaffected by setwd()
#    because the `here` package resolves the project root once, from
#    project-root markers (e.g. .git), independent of the current working
#    directory at call time.
# ---------------------------------------------------------------------------
cat("\n=== Sourcing RealData_Collection.R (10 existing datasets) ===\n")
source(here::here("data", "outlier_detection", "RealData_Collection.R"))

# Restore working directory to repo root for the rest of this script.
setwd(repo_root)
cat("Working directory reset to:", getwd(), "\n")

# The 10 datasets used in the manuscript (CLAUDE.md Section 6):
# Hepatitis, Lymphography, Glass, WBC, Vertebral, Stamps, WDBC, Vowels,
# Thyroid, Wilt. RealData_Collection.R additionally builds ecoli,
# pageblocks, PenDigits, pima, shuffle (shuttle), waveform as a side effect
# of sourcing -- these are NOT part of the manuscript's 10-dataset real-data
# benchmark set and are intentionally not exported here.
existing_datasets <- list(
  Hepatitis     = hepatitis,
  Lymphography  = lymphography,
  Glass         = glass,
  WBC           = WBC,
  Vertebral     = vertebral,
  Stamps        = stamps,
  WDBC          = WDBC,
  Vowels        = vowels,
  Thyroid       = thyroid,
  Wilt          = wilt
)

existing_sources <- list(
  Hepatitis    = "data/outlier_detection/Hepatitis_withoutdupl_10_v01.arff (DAMI/ODDS collection)",
  Lymphography = "data/outlier_detection/Lymphography_withoutdupl_norm_idf.arff (DAMI/ODDS collection)",
  Glass        = "data/outlier_detection/glass.mat (ODDS)",
  WBC          = "data/outlier_detection/WBC_withoutdupl_norm_v02.arff (DAMI/ODDS collection)",
  Vertebral    = "data/outlier_detection/vertebral.mat (ODDS)",
  Stamps       = "data/outlier_detection/Stamps_withoutdupl_norm_09.arff (DAMI/ODDS collection)",
  WDBC         = "data/outlier_detection/WDBC_withoutdupl_norm_v02.arff (DAMI/ODDS collection)",
  Vowels       = "data/outlier_detection/vowels.mat (ODDS)",
  Thyroid      = "data/outlier_detection/thyroid.mat (ODDS)",
  Wilt         = "data/outlier_detection/Wilt_withoutdupl_05.arff (DAMI/ODDS collection)"
)

existing_normalization <- list(
  Hepatitis    = "z-score (scale) per column [RealData_Collection.R]",
  Lymphography = "none (pre-normalized source arff) [RealData_Collection.R]",
  Glass        = "z-score (scale) per column [RealData_Collection.R]",
  WBC          = "none (pre-normalized source arff) [RealData_Collection.R]",
  Vertebral    = "robust median/MADN (scale_R) per column [RealData_Collection.R]",
  Stamps       = "none (pre-normalized source arff) [RealData_Collection.R]",
  WDBC         = "none (pre-normalized source arff) [RealData_Collection.R]",
  Vowels       = "robust median/MADN (scale_R) per column [RealData_Collection.R]",
  Thyroid      = "robust median/MADN (scale_R) per column [RealData_Collection.R]",
  Wilt         = "robust median/MADN (scale_R) per column [RealData_Collection.R]"
)

# ---------------------------------------------------------------------------
# 2. Load the 4 new high-dimensional datasets from the raw CSVs produced by
#    revision_experiments/convert_highdim.py, then apply THE SAME
#    normalization philosophy, dedup, and row-ordering conventions as
#    RealData_Collection.R -- with one necessary generalization, documented
#    below.
#
#    LABEL POLARITY: convert_highdim.py writes `label` using the SOURCE
#    (ODDS/ADBench) convention: 1 = outlier, 0 = regular. This repo's
#    convention (and RealData_Collection.R's own output) is the OPPOSITE:
#    1 = regular, 0 = outlier. We flip polarity below with the identical
#    idiom RealData_Collection.R uses for its own raw .mat sources:
#    ifelse(y == 1, 0, 1).
#
#    NORMALIZATION: RealData_Collection.R's scale_R() = (x - median(x)) /
#    mad(x) is applied blindly, column-wise, to its raw .mat sources
#    (vertebral/vowels/thyroid: 6-12 features each). It has no fallback for
#    mad(x) == 0 because none of the existing low-dimensional datasets ever
#    hit that edge case. Our new high-dimensional datasets DO hit it:
#      - InternetAds: ALL 1555 feature columns are binary {0,1}; every
#        single column has mad(x) == 0 (per-task expectation).
#      - Musk: 3 of 166 continuous columns have mad(x) == 0 despite nonzero
#        variance (skewed/near-constant distributions).
#      - Arrhythmia: 149 of 274 columns have mad(x) == 0 (17 of those are
#        fully constant; the rest are skewed clinical indicators where the
#        mode covers >50% of rows).
#      - Speech: no zero-MADN columns; scale_R applies cleanly.
#    To avoid injecting NaNs (forbidden by the T1 gate) while never
#    silently dropping a column (d must stay exact), we use a per-column
#    safe fallback in the SAME spirit as this project's other documented
#    MADN-fallback convention (CLAUDE.md: IOS standardization falls back to
#    SD when MADN==0, then to a constant when SD==0 too):
#      mad(x) > 0  -> (x - median(x)) / mad(x)          [scale_R, unchanged]
#      mad(x) == 0, sd(x) > 0 -> (x - median(x)) / sd(x) [SD fallback]
#      mad(x) == 0, sd(x) == 0 (constant column) -> x - median(x) (= 0)
#    For datasets where EVERY feature column is binary {0,1} (InternetAds),
#    we instead keep the entire feature matrix unscaled, per the task's
#    explicit instruction, since scaling a purely binary matrix column-wise
#    with 100% MADN==0 columns is not a fallback-of-a-few-columns situation
#    but a systematic mismatch between the normalization and the data type.
# ---------------------------------------------------------------------------
cat("\n=== Loading 4 new high-dimensional datasets from raw CSVs ===\n")

raw_csv_dir <- here::here("data", "outlier_detection")

stopifnot("scale_R() not found after sourcing RealData_Collection.R" = exists("scale_R"))

scale_R_safe <- function(x) {
  M <- median(x)
  madn <- mad(x)
  if (madn > 0) {
    return((x - M) / madn)
  }
  sdx <- sd(x)
  if (sdx > 0) {
    return((x - M) / sdx)
  }
  return(x - M)  # constant column -> all zeros, no NaN
}

load_highdim <- function(name, csv_name) {
  path <- file.path(raw_csv_dir, csv_name)
  stopifnot("raw CSV not found" = file.exists(path))
  df <- read.csv(path)

  label_col <- ncol(df)
  stopifnot("expected final column named 'label'" = colnames(df)[label_col] == "label")

  X_raw <- as.matrix(df[, -label_col, drop = FALSE])
  y_raw <- df[, label_col]  # source convention: 1 = outlier
  d <- ncol(X_raw)

  col_is_binary <- apply(X_raw, 2, function(col) all(col %in% c(0, 1)))
  all_binary <- all(col_is_binary)

  if (all_binary) {
    col_madn <- apply(X_raw, 2, mad)
    n_zero_madn <- sum(col_madn == 0)
    cat(sprintf(
      "  [%s] all %d feature columns are binary {0,1}; %d/%d columns have MADN==0 (scale_R would produce NaN there).\n",
      name, d, n_zero_madn, d
    ))
    cat(sprintf("  [%s] DECISION: keeping features UNSCALED (raw 0/1) to avoid systematic NaN injection.\n", name))
    X <- X_raw
    normalization_used <- "none (all features binary 0/1; scale_R would yield NaN on every column -- kept raw per task instruction)"
  } else {
    col_madn <- apply(X_raw, 2, mad)
    n_zero_madn <- sum(col_madn == 0)
    if (n_zero_madn > 0) {
      cat(sprintf(
        "  [%s] %d/%d columns have MADN==0; applying SD fallback (or 0 for constant columns) on those columns only.\n",
        name, n_zero_madn, d
      ))
    }
    X <- apply(X_raw, 2, scale_R_safe)
    if (n_zero_madn > 0) {
      normalization_used <- sprintf(
        "robust median/MADN (scale_R) per column, with SD/constant fallback on %d/%d zero-MADN columns",
        n_zero_madn, d
      )
    } else {
      normalization_used <- "robust median/MADN (scale_R) per column [RealData_Collection.R convention]"
    }
  }

  # Flip label polarity to repo convention: 1 = regular, 0 = outlier.
  y <- ifelse(y_raw == 1, 0, 1)

  mat <- cbind(X, y)
  colnames(mat)[ncol(mat)] <- "label"

  # Never silently drop columns: assert d is unchanged.
  stopifnot("feature count changed unexpectedly" = ncol(mat) - 1 == d)

  # Remove duplicated rows, exactly as RealData_Collection.R does for every dataset.
  n_before <- nrow(mat)
  mat <- as.matrix(distinct(as.data.frame(mat)))
  n_after <- nrow(mat)
  if (n_after != n_before) {
    cat(sprintf(
      "  [%s] WARNING: distinct() removed %d duplicate row(s) (%d -> %d).\n",
      name, n_before - n_after, n_before, n_after
    ))
  }

  # Move outliers to the end (decreasing sort on label -> regular(1) first, outlier(0) last).
  label_idx <- ncol(mat)
  mat <- mat[order(mat[, label_idx], decreasing = TRUE), ]

  has_na <- anyNA(mat)
  cat(sprintf(
    "  [%s] n=%d d=%d n_outliers=%d anyNA=%s normalization=%s\n",
    name, nrow(mat), ncol(mat) - 1, sum(mat[, label_idx] == 0), has_na, normalization_used
  ))
  if (has_na) stop(sprintf("[%s] NA present after normalization", name))

  list(mat = mat, normalization = normalization_used)
}

musk_res        <- load_highdim("Musk",        "musk_raw.csv")
arrhythmia_res  <- load_highdim("Arrhythmia",  "arrhythmia_raw.csv")
speech_res      <- load_highdim("Speech",      "speech_raw.csv")
internetads_res <- load_highdim("InternetAds", "internetads_raw.csv")

new_datasets <- list(
  Musk        = musk_res$mat,
  Arrhythmia  = arrhythmia_res$mat,
  Speech      = speech_res$mat,
  InternetAds = internetads_res$mat
)

new_sources <- list(
  Musk        = "https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/25_musk.npz (ADBench mirror of ODDS musk.mat; repo commit 3dac8221081e190f157d78e93bfa8867f90d0965, fetched 2026-07-10)",
  Arrhythmia  = "https://raw.githubusercontent.com/BELLoney/Outlier-detection/master/ODDS/original/arrhythmiaori.mat (community ODDS mirror; repo commit a90424a2322f863209534a23e125ce175c6a75d5, fetched 2026-07-10; odds.cs.stonybrook.edu unreachable via TLS handshake failure at acquisition time, and Arrhythmia is absent from ADBench's Classical set)",
  Speech      = "https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/36_speech.npz (ADBench mirror of ODDS speech.mat; repo commit 3dac8221081e190f157d78e93bfa8867f90d0965, fetched 2026-07-10)",
  InternetAds = "https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical/17_InternetAds.npz (ADBench mirror of ODDS InternetAds.mat; repo commit 3dac8221081e190f157d78e93bfa8867f90d0965, fetched 2026-07-10)"
)

new_normalization <- list(
  Musk        = musk_res$normalization,
  Arrhythmia  = arrhythmia_res$normalization,
  Speech      = speech_res$normalization,
  InternetAds = internetads_res$normalization
)

# ---------------------------------------------------------------------------
# 3. Export all 14 datasets to CSV + manifest.
# ---------------------------------------------------------------------------
out_dir <- here::here("revision_experiments", "results", "datasets_csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_datasets <- c(existing_datasets, new_datasets)
all_sources <- c(existing_sources, new_sources)
all_normalization <- c(existing_normalization, new_normalization)

cat("\n=== Exporting all 14 datasets to CSV ===\n")

manifest_rows <- list()
for (nm in names(all_datasets)) {
  mat <- all_datasets[[nm]]
  d <- ncol(mat) - 1
  n <- nrow(mat)
  label_col <- ncol(mat)
  n_outliers <- sum(mat[, label_col] == 0)

  df_out <- as.data.frame(mat)
  colnames(df_out) <- c(paste0("V", seq_len(d)), "label")

  out_path <- file.path(out_dir, paste0(nm, ".csv"))
  write.csv(df_out, out_path, row.names = FALSE)

  na_count <- sum(is.na(df_out[, seq_len(d), drop = FALSE]))

  cat(sprintf(
    "  %-12s n=%-5d d=%-5d n_outliers=%-4d NA_in_features=%-4d -> %s\n",
    nm, n, d, n_outliers, na_count, out_path
  ))

  manifest_rows[[nm]] <- data.frame(
    dataset = nm,
    n = n,
    d = d,
    n_outliers = n_outliers,
    source_url = all_sources[[nm]],
    normalization = all_normalization[[nm]],
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest_rows)
rownames(manifest) <- NULL

manifest_path <- file.path(out_dir, "manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
cat("\nManifest written to:", manifest_path, "\n")

# ---------------------------------------------------------------------------
# 4. Verification
# ---------------------------------------------------------------------------
cat("\n=== Manifest ===\n")
print(manifest[, c("dataset", "n", "d", "n_outliers")], row.names = FALSE)

cat("\n=== Verification: expected gate values for the 4 new datasets ===\n")
expected <- list(
  Musk        = c(n = 3062, d = 166,  n_outliers = 97),
  Arrhythmia  = c(n = 452,  d = 274,  n_outliers = 66),
  Speech      = c(n = 3686, d = 400,  n_outliers = 61),
  InternetAds = c(n = 1966, d = 1555, n_outliers = 368)
)

all_gate_ok <- TRUE
for (nm in names(expected)) {
  row <- manifest[manifest$dataset == nm, ]
  exp <- expected[[nm]]
  ok <- (row$n == exp["n"]) && (row$d == exp["d"]) && (row$n_outliers == exp["n_outliers"])
  cat(sprintf(
    "  %-12s expected n=%d d=%d n_outliers=%d | got n=%d d=%d n_outliers=%d | MATCH=%s\n",
    nm, exp["n"], exp["d"], exp["n_outliers"], row$n, row$d, row$n_outliers, ok
  ))
  if (!ok) all_gate_ok <- FALSE
}
stopifnot("One or more new datasets FAILED the expected (n, d, n_outliers) gate" = all_gate_ok)

cat("\n=== Verification: no NA in any exported feature matrix; row/col counts ===\n")
all_na_ok <- TRUE
for (nm in names(all_datasets)) {
  path <- file.path(out_dir, paste0(nm, ".csv"))
  df_check <- read.csv(path)
  feat_cols <- setdiff(colnames(df_check), "label")
  na_present <- anyNA(df_check[, feat_cols, drop = FALSE])
  n_rows <- nrow(df_check)
  n_cols <- ncol(df_check)
  expected_d <- ncol(all_datasets[[nm]]) - 1
  expected_n <- nrow(all_datasets[[nm]])
  dims_ok <- ((n_cols - 1) == expected_d) && (n_rows == expected_n)
  cat(sprintf(
    "  %-12s rows=%-5d cols=%-5d (features=%-5d) anyNA=%-5s dims_ok=%s\n",
    nm, n_rows, n_cols, n_cols - 1, na_present, dims_ok
  ))
  if (na_present || !dims_ok) all_na_ok <- FALSE
}
stopifnot("NA present in an exported feature matrix, or dimension mismatch" = all_na_ok)

cat("\n=== Spot-check: WBC round-trips correctly ===\n")
wbc_csv <- read.csv(file.path(out_dir, "WBC.csv"))
cat(sprintf(
  "  WBC.csv: n=%d d=%d n_outliers=%d anyNA=%s\n",
  nrow(wbc_csv), ncol(wbc_csv) - 1, sum(wbc_csv$label == 0), anyNA(wbc_csv)
))
stopifnot(
  "WBC row count mismatch vs in-memory object" = nrow(wbc_csv) == nrow(existing_datasets$WBC),
  "WBC column count mismatch vs in-memory object" = ncol(wbc_csv) == ncol(existing_datasets$WBC)
)

cat("\nALL CHECKS PASSED.\n")
cat("02_load_data.R completed successfully.\n")
