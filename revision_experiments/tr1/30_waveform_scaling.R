#!/usr/bin/env Rscript
# 30_waveform_scaling.R -- last standing hypothesis for waveform.
#
# Established so far:
#   - all four detectors under-flag on waveform (TPR 0.100-0.360 below
#     published) while TNR sits slightly ABOVE published
#   - the loader sorts waveform correctly, so the glass/ecoli positional
#     counting bug does not apply here
#   - alpha is irrelevant: UN-MCCD and SUN-MCCD give BIT-IDENTICAL results at
#     the 99% and 999% NN tables (29_waveform_alpha.R). At d=21 the
#     spatial-randomness test is saturated -- the threshold no longer moves
#     the decision.
#
# That leaves preprocessing. RealData_Collection.R MADN-scales waveform
# (scale_R). If the published run used the classical z-score, or the raw
# features, the covering radii and hence the whole clustering would differ.
#
# Runs U-MCCD (the uniform-coverage reference) and SUN-MCCD (the headline) on
# waveform under raw / z-score / MADN, re-reading the source .arff so the raw
# variant is genuinely unscaled.
#
# Writes results/waveform_scaling.csv.

suppressMessages({
  library(here); library(foreign); library(dplyr)
})
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

OUT   <- here::here("revision_experiments/results/tr1/waveform_scaling.csv")
S_MIN <- 0.0625

TRUTH <- read.csv(here::here("revision_experiments/tr1/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
published <- function(meth) {
  s <- TRUTH[TRUTH$dataset == "waveform" & TRUTH$method == meth, ]
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

# scale_R, verbatim from RealData_Collection.R:10 (median / normalized MAD).
scale_R <- function(x) (x - median(x)) / (mad(x) )

owd <- setwd(here::here("data/outlier_detection"))
raw <- foreign::read.arff("Waveform_withoutdupl_v01.arff")
setwd(owd)
raw$id <- NULL
raw$outlier <- ifelse(raw$outlier == "yes", 0, 1)

build <- function(mode) {
  Xr <- subset(raw, select = -outlier)
  Xs <- switch(mode,
               raw    = as.matrix(Xr),
               zscore = apply(Xr, 2, scale),
               madn   = apply(Xr, 2, scale_R))
  df <- as.matrix(distinct(as.data.frame(cbind(Xs, raw$outlier))))
  df <- df[order(df[, ncol(df)], decreasing = TRUE), ]
  list(X = df[, 1:(ncol(df) - 1), drop = FALSE], Y = df[, ncol(df)],
       d = ncol(df) - 1, n = nrow(df))
}

for (mode in c("raw", "zscore", "madn")) {
  dat <- build(mode)
  cat(sprintf("\n-- %s: n=%d d=%d n0=%d, feature range [%.3f, %.3f]\n",
              mode, dat$n, dat$d, sum(dat$Y == 0), min(dat$X), max(dat$X)))

  for (meth in c("U-MCCD", "SUN-MCCD")) {
    keys <- c(dataset = "waveform", method = meth, scaling = mode)
    if (isTRUE(has_result(OUT, keys))) { cat(sprintf("[skip] %s %s\n", meth, mode)); next }

    t0 <- Sys.time()
    out <- tryCatch({
      res <- if (meth == "SUN-MCCD") {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
      } else {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
      }
      list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
           status = "ok", note = NA_character_)
    }, error = function(e) {
      list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
           status = "error", note = conditionMessage(e))
    })
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    pub <- published(meth)
    dif <- if (any(is.na(out$m))) NA_real_ else max(abs(round(out$m, 3) - pub))

    append_result(OUT, list(
      dataset = "waveform", method = meth, scaling = mode, s_min = S_MIN,
      n = dat$n, d = dat$d,
      TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
      BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
      published_TPR = unname(pub[["TPR"]]), published_TNR = unname(pub[["TNR"]]),
      published_BA = unname(pub[["BA"]]), published_F2 = unname(pub[["F2"]]),
      max_abs_diff = dif, match_3dp = !is.na(dif) && dif < 5e-4,
      t_total = wall, status = out$status, note = out$note,
      timestamp = format(Sys.time())
    ))

    cat(sprintf("  %-9s %-7s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | pub %.3f/%.3f/%.3f/%.3f | %s | %.1fs\n",
                meth, mode, out$m[["TPR"]], out$m[["TNR"]], out$m[["BA"]], out$m[["F2"]],
                pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]],
                if (!is.na(dif) && dif < 5e-4) "MATCH" else sprintf("diff %.3f", dif), wall))
  }
}
cat("\n30_waveform_scaling.R: done\n")
