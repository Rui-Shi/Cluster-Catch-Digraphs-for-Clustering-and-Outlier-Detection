# 80_export_datasets_wp4.R
#
# WP4 (modern baselines, AE.3 / R1.3 / R1.5 / R3.3 / R6.5) -- data export.
#
# Exports the sixteen Section 6 real-data sets to CSV so that the Python
# competitor driver (tr1/81_wp4_baselines.py) sees EXACTLY the feature
# matrices the four proposed detectors saw. Preprocessing is entirely the R
# loader's: nothing is scaled, imputed, deduplicated or reordered here.
#
# Convention of the exported `label` column: 1 = regular, 0 = outlier
# (RealData_Collection.R / harness.R convention). PyOD's convention is the
# OPPOSITE; the Python driver never writes into this column and names its own
# native-label columns `is_outlier` (1 = outlier). See tr1/WP4_PROTOCOL.md.
#
# Gate: (n, d, n_outliers) must match Table tab:Real_Data of
# CCD_OutlierDetection_Neurocomputing.tex for all sixteen sets. Any single
# disagreement aborts the script before anything is written to the manifest.
#
# Run from the CLONE root:
#   Rscript "revision_experiments/tr1/80_export_datasets_wp4.R"
#
# Idempotent: re-running overwrites the same CSVs.

suppressPackageStartupMessages({
  library(here)
})

source(here::here("revision_experiments", "shared", "harness.R"))  # read-only

out_dir <- here::here("revision_experiments", "results", "tr1", "wp4", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Export name -> R object name created by RealData_Collection.R.
# NOTE: "Shuttle" is the object `shuffle` (the loader's own spelling).
obj_of <- c(
  hepatitis    = "hepatitis",
  lymphography = "lymphography",
  glass        = "glass",
  WBC          = "WBC",
  vertebral    = "vertebral",
  ecoli        = "ecoli",
  stamps       = "stamps",
  WDBC         = "WDBC",
  pima         = "pima",
  Shuttle      = "shuffle",
  vowels       = "vowels",
  PenDigits    = "PenDigits",
  waveform     = "waveform",
  thyroid      = "thyroid",
  pageblocks   = "pageblocks",
  wilt         = "wilt"
)

# Scaling regime applied by RealData_Collection.R, transcribed from the loader.
scaling_of <- c(
  hepatitis    = "classical z-score (scale)",
  lymphography = "none (source arff pre-normalized)",
  glass        = "classical z-score (scale)",
  WBC          = "none (source arff pre-normalized)",
  vertebral    = "robust median/MADN (scale_R)",
  ecoli        = "none (raw csv features)",
  stamps       = "none (source arff pre-normalized)",
  WDBC         = "none (source arff pre-normalized)",
  pima         = "none (source arff pre-normalized)",
  Shuttle      = "classical z-score (scale)",
  vowels       = "robust median/MADN (scale_R)",
  PenDigits    = "robust median/MADN (scale_R)",
  waveform     = "robust median/MADN (scale_R)",
  thyroid      = "robust median/MADN (scale_R)",
  pageblocks   = "robust median/MADN (scale_R)",
  wilt         = "robust median/MADN (scale_R)"
)

# Table tab:Real_Data, CCD_OutlierDetection_Neurocomputing.tex lines 1135-1150.
truth <- data.frame(
  dataset    = c("hepatitis", "lymphography", "glass", "WBC", "vertebral",
                 "ecoli", "stamps", "WDBC", "pima", "Shuttle", "vowels",
                 "PenDigits", "waveform", "thyroid", "pageblocks", "wilt"),
  n          = c(74, 148, 213, 223, 240, 336, 340, 367, 555, 1013, 1452,
                 3200, 3443, 3656, 4795, 4819),
  d          = c(19, 18, 9, 9, 6, 7, 9, 30, 8, 9, 12, 16, 21, 6, 10, 5),
  n_outliers = c(7, 6, 9, 10, 30, 8, 31, 10, 55, 13, 46, 20, 100, 93, 510, 257),
  stringsAsFactors = FALSE
)

cat("=== WP4 data export: 16 real-data sets ===\n")
cat("Output dir:", out_dir, "\n\n")

rows <- list()
all_ok <- TRUE

for (nm in truth$dataset) {
  dat <- load_real_dataset(obj_of[[nm]])
  X <- dat$X
  Y <- dat$Y
  n <- nrow(X)
  d <- ncol(X)
  n0 <- sum(Y == 0)

  exp_row <- truth[truth$dataset == nm, ]
  ok <- (n == exp_row$n) && (d == exp_row$d) && (n0 == exp_row$n_outliers)
  if (!ok) all_ok <- FALSE

  na_feat <- sum(is.na(X))
  if (na_feat > 0) all_ok <- FALSE

  df_out <- as.data.frame(X)
  colnames(df_out) <- paste0("V", seq_len(d))
  df_out$label <- as.integer(Y)

  write.csv(df_out, file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)

  cat(sprintf(
    "  %-13s n=%-5d d=%-3d n0=%-4d | table n=%-5d d=%-3d n0=%-4d | MATCH=%-5s NA=%d\n",
    nm, n, d, n0, exp_row$n, exp_row$d, exp_row$n_outliers, ok, na_feat
  ))

  rows[[nm]] <- data.frame(
    dataset          = nm,
    r_object         = obj_of[[nm]],
    n                = n,
    d                = d,
    n_outliers       = n0,
    contamination    = n0 / n,
    scaling          = scaling_of[[nm]],
    table_n          = exp_row$n,
    table_d          = exp_row$d,
    table_n_outliers = exp_row$n_outliers,
    table_match      = ok,
    na_in_features   = na_feat,
    label_convention = "1 = regular, 0 = outlier",
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, rows)
rownames(manifest) <- NULL

if (!all_ok) {
  print(manifest[, c("dataset", "n", "d", "n_outliers", "table_n", "table_d",
                     "table_n_outliers", "table_match", "na_in_features")])
  stop("STOP: at least one data set disagrees with Table tab:Real_Data, or ",
       "carries NA features. No manifest written.")
}

write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

cat("\nAll 16 data sets match Table tab:Real_Data exactly; no NA features.\n")
cat("Manifest written to:", file.path(out_dir, "manifest.csv"), "\n")
