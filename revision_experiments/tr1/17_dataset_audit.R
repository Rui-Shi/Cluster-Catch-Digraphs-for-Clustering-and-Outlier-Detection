# 17_dataset_audit.R -- ground truth for Table 5 (tab:Real_Data).
# Prints n, d, n0 and contamination for each of the eight real datasets exactly
# as load_real_dataset() returns them. Writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))

cat(sprintf("\n%-11s %6s %4s %6s %8s\n", "dataset", "n", "d", "n0", "pct"))
for (ds in c("hepatitis", "glass", "vertebral", "ecoli", "stamps", "vowels", "waveform", "wilt")) {
  d <- load_real_dataset(ds)
  n0 <- sum(d$Y == 0)
  cat(sprintf("%-11s %6d %4d %6d %7.1f%%\n", ds, d$n, d$d, n0, 100 * n0 / d$n))
}
cat("\ndone\n")
