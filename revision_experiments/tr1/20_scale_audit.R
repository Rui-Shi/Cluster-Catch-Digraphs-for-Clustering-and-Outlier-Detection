# 20_scale_audit.R -- what scale are the features actually on, per dataset?
#
# RealData_Collection.R applies scale_R (median/MADN) to vertebral, vowels,
# waveform and wilt; plain scale (z-score) to glass and hepatitis; and nothing
# to ecoli and stamps. Before writing a corrected preprocessing sentence we need
# to know whether the two unscaled sets are already on comparable scales (which
# would justify leaving them) or simply unscaled.
#
# Writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))

cat(sprintf("\n%-11s %8s %9s %9s %9s %9s\n",
            "dataset", "d", "min", "max", "max|sd|", "sd ratio"))
for (ds in c("hepatitis", "glass", "vertebral", "ecoli",
             "stamps", "vowels", "waveform", "wilt")) {
  X  <- as.matrix(load_real_dataset(ds)$X)
  sds <- apply(X, 2, sd)
  cat(sprintf("%-11s %8d %9.3f %9.3f %9.3f %9.1f\n",
              ds, ncol(X), min(X), max(X), max(sds),
              max(sds) / max(min(sds), 1e-12)))
}
cat("\nsd ratio = largest column sd / smallest column sd.",
    "\nNear 1 means columns are already on a common scale.\n")
