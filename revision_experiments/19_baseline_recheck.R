# 19_baseline_recheck.R -- resolve three internally inconsistent baseline rows
# found by 18_table_consistency.R.
#
#   DBSCAN / wilt      TPR 0.000, TNR 0.959, BA 0.673, F2 0.381
#                      -> impossible: TPR 0 forces BA = TNR/2 = 0.480 and F2 = 0
#   LOF / ecoli        BA off by 0.004 from (TPR+TNR)/2
#   LOF / vertebral    BA off by 0.002
#
# These are baselines in METHOD_REGISTRY, so they can simply be re-run.
# Writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))

TRUTH <- read.csv(here::here("revision_experiments", "published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

CASES <- list(c("wilt", "DBSCAN"), c("ecoli", "LOF"), c("vertebral", "LOF"))

for (cs in CASES) {
  ds <- cs[1]; meth <- cs[2]
  out <- tryCatch({
    dat <- load_real_dataset(ds)
    res <- METHOD_REGISTRY[[meth]](dat$X, dat$d, dat$Y)
    thr <- REAL_DATA_THRESHOLDS[[meth]]
    ev  <- evaluate(dat$Y, res$score, thr)
    pub <- published(ds, meth)
    cat(sprintf("\n%s / %s   (n=%d, d=%d, n0=%d, threshold=%s)\n",
                ds, meth, dat$n, dat$d, sum(dat$Y == 0), format(thr)))
    cat(sprintf("  published  TPR %.3f  TNR %.3f  BA %.3f  F2 %.3f   -> BA check %.4f\n",
                pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]],
                (pub[["TPR"]] + pub[["TNR"]]) / 2))
    cat(sprintf("  recomputed TPR %.3f  TNR %.3f  BA %.3f  F2 %.3f   -> BA check %.4f\n",
                ev[["TPR"]], ev[["TNR"]], ev[["BA"]], ev[["F2"]],
                (ev[["TPR"]] + ev[["TNR"]]) / 2))
    TRUE
  }, error = function(e) {
    cat(sprintf("\n%s / %s  ERROR: %s\n", ds, meth, conditionMessage(e))); FALSE
  })
}

cat("\ndone\n")
