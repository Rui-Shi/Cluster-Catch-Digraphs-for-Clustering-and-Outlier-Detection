# 26_loader_sort_audit.R -- how badly does the loader's "move outliers to the
# end" step actually fail, and on which datasets?
#
# RealData_Collection.R ends each dataset block with
#     <ds> = <ds>[order(<ds>[,K], decreasing=T),]
# intending K to be the label column (= d+1). Two blocks use the wrong K:
#     glass  line 34  sorts by column 9 (Fe, a feature); label is column 10
#     ecoli  line 112 sorts by column 6 (alm1, a feature); label is column 8
#
# This matters because the PUBLISHED numbers were produced with count_scores2 /
# count_DBSCAN / count_MST2 / count_ODIN, which never look up Y per row -- they
# slice positionally, treating the first n-n0 rows as regular and the last n0
# as outliers. If the sort did not actually put the outliers last, the
# published row for that dataset is scored against the wrong ground truth.
#
# harness.R's evaluate() reorders (Y, score) jointly before counting, so OUR
# numbers are unaffected. This script quantifies the damage to the published
# ones.
#
# Read-only. Writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))

DATASETS <- c("hepatitis", "glass", "vertebral", "ecoli",
              "stamps", "vowels", "waveform", "wilt")

cat(sprintf("\n%-11s %6s %4s %5s %10s %10s %14s\n",
            "dataset", "n", "d", "n0", "sorted?", "last-n0 TP", "positional TPR cap"))
cat(strrep("-", 74), "\n")

for (ds in DATASETS) {
  dat <- load_real_dataset(ds)
  Y <- dat$Y
  n <- length(Y); n0 <- sum(Y == 0)

  # Is the data actually arranged regulars-first, outliers-last?
  ideal <- c(rep(1, n - n0), rep(0, n0))
  sorted_ok <- identical(as.numeric(Y), as.numeric(ideal))

  # Of the last n0 rows -- the ones count_scores2 CALLS outliers -- how many
  # really are outliers? This is the ceiling on the TPR any positional count
  # can report, no matter how good the detector is.
  tail_true <- sum(Y[(n - n0 + 1):n] == 0)

  cat(sprintf("%-11s %6d %4d %5d %10s %10s %14.3f\n",
              ds, n, dat$d, n0,
              if (sorted_ok) "OK" else "MIS-SORTED",
              sprintf("%d/%d", tail_true, n0),
              tail_true / n0))
}

cat("\n'positional TPR cap' is the largest TPR the published counting code\n",
    "could report for that dataset even with a perfect detector.\n", sep = "")

# Where the outliers actually sit, for the mis-sorted ones.
cat("\n--- outlier row positions on the mis-sorted datasets ---\n")
for (ds in DATASETS) {
  dat <- load_real_dataset(ds)
  Y <- dat$Y; n <- length(Y); n0 <- sum(Y == 0)
  if (identical(as.numeric(Y), as.numeric(c(rep(1, n - n0), rep(0, n0))))) next
  pos <- which(Y == 0)
  cat(sprintf("  %-10s n=%d, n0=%d, outliers at rows: %s\n",
              ds, n, n0,
              paste(pos, collapse = ", ")))
}
cat("\ndone\n")
