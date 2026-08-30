# 27_positional_vs_correct.R -- the decisive test for glass and ecoli.
#
# 26_loader_sort_audit.R established that the loader mis-sorts glass (sorts by
# column 9, a feature; label is column 10) and ecoli (sorts by column 6; label
# is column 8), so on those two datasets the last n0 rows contain NO true
# outliers at all. The published numbers were produced by count_scores2, which
# scores POSITIONALLY -- it calls the last n0 rows the outliers without ever
# looking at Y.
#
# So there are two possible numbers for every cell:
#   correct     evaluate(Y, score, thr)         -- reorders (Y, score) first
#   positional  count_scores2(Y, score, thr)    -- takes the loader's row order
#                                                  at face value, as published
#
# If the published row equals the POSITIONAL number computed from our scores,
# then our detectors reproduce the published run exactly and the published
# table is simply scored against the wrong ground truth. That fully closes
# glass and ecoli, and it means the regenerated rows are corrections, not
# disagreements.
#
# Read-only. Writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

TRUTH <- read.csv(here::here("revision_experiments/tr1/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  if (nrow(s) == 0) return(NULL)
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

# The published counting path, verbatim: no reordering.
positional <- function(Y, score, threshold) {
  v <- count_scores2(Y, score, threshold)
  names(v) <- c("TPR", "TNR", "BA", "F2")
  v
}

S_MIN <- 0.0625   # recovered in 25_smin_match.R
METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")

for (ds in c("glass", "ecoli")) {
  dat <- load_real_dataset(ds)
  cat(sprintf("\n===== %s (n=%d, d=%d, n0=%d) =====\n",
              ds, dat$n, dat$d, sum(dat$Y == 0)))
  cat(sprintf("%-9s %-11s %7s %7s %7s %7s   %s\n",
              "method", "counting", "TPR", "TNR", "BA", "F2", "vs published"))
  cat(strrep("-", 72), "\n")

  for (meth in METHODS) {
    res <- tryCatch({
      if (meth %in% c("SU-MCCD", "SUN-MCCD")) {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
      } else {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
      }
    }, error = function(e) { cat(sprintf("  %-9s ERROR: %s\n", meth, conditionMessage(e))); NULL })
    if (is.null(res)) next

    thr <- REAL_DATA_THRESHOLDS[[meth]]
    pub <- published(ds, meth)

    cor_m <- evaluate(dat$Y, res$score, thr)
    pos_m <- positional(dat$Y, res$score, thr)

    for (nm in c("correct", "positional")) {
      m <- if (nm == "correct") cor_m else pos_m
      dif <- if (is.null(pub)) NA_real_ else max(abs(round(m, 3) - pub))
      tag <- if (is.na(dif)) "" else if (dif < 5e-4) "<<< MATCHES PUBLISHED" else sprintf("max diff %.3f", dif)
      cat(sprintf("%-9s %-11s %7.3f %7.3f %7.3f %7.3f   %s\n",
                  if (nm == "correct") meth else "", nm,
                  m[["TPR"]], m[["TNR"]], m[["BA"]], m[["F2"]], tag))
    }
    if (!is.null(pub)) {
      cat(sprintf("%-9s %-11s %7.3f %7.3f %7.3f %7.3f\n",
                  "", "PUBLISHED", pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]]))
    }
    cat("\n")
  }
}
cat("done\n")
