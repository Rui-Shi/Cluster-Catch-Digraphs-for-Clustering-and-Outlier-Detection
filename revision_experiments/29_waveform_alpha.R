#!/usr/bin/env Rscript
# 29_waveform_alpha.R -- is waveform's shortfall an alpha-level artefact?
#
# All four detectors under-flag on waveform: our TPR is well below the
# published one on every method (U-MCCD 0.770 vs 0.870, SU-MCCD 0.470 vs
# 0.830, UN-MCCD 0.540 vs 0.860). A uniform under-flagging across methods
# points at the spatial-randomness test being too permissive, not at any one
# detector.
#
# waveform is d=21, and the paper-faithful resolvers force the 999% table
# there (nn_quant_label_paper_UN: d>=20 -> 999; nn_quant_label_paper_SUN:
# d>=10 -> 999; rk_quant_label_paper: d>=10 -> 999). 999% is the most
# permissive level in the family: the test rarely rejects, clusters stay
# coarse, and fewer points fall outside their cluster's mutual-catch-graph
# component -- exactly the observed signature.
#
# A table inventory shows the option was not symmetric:
#   RK at d=21: only 999% exists
#   NN at d=21: BOTH 99% and 999% exist
# so the two NND-based detectors could have been run at 99%, and if the
# published run was, that alone would explain their shortfall.
#
# This script runs UN-MCCD and SUN-MCCD on waveform at quant="99" and compares
# against both the published row and our 999% result.
#
# Writes results/waveform_alpha.csv.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

OUT <- here::here("revision_experiments/results/tr1/waveform_alpha.csv")

TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

dat   <- load_real_dataset("waveform")
S_MIN <- 0.0625

CELLS <- list(
  list(meth = "UN-MCCD",  quant = "99"),
  list(meth = "SUN-MCCD", quant = "99")
)

for (cl in CELLS) {
  meth <- cl$meth; q <- cl$quant
  keys <- c(dataset = "waveform", method = meth, quant_label = q)
  if (isTRUE(has_result(OUT, keys))) { cat(sprintf("[skip] %s q=%s\n", meth, q)); next }

  t0 <- Sys.time()
  out <- tryCatch({
    res <- if (meth == "SUN-MCCD") {
      METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN, quant = q)
    } else {
      METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, quant = q)
    }
    list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
         status = "ok", note = NA_character_)
  }, error = function(e) {
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
         status = "error", note = conditionMessage(e))
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  pub <- published("waveform", meth)
  dif <- if (any(is.na(out$m))) NA_real_ else max(abs(round(out$m, 3) - pub))

  append_result(OUT, list(
    dataset = "waveform", method = meth, quant_label = q, s_min = S_MIN,
    n = dat$n, d = dat$d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
    published_TPR = unname(pub[["TPR"]]), published_TNR = unname(pub[["TNR"]]),
    published_BA = unname(pub[["BA"]]), published_F2 = unname(pub[["F2"]]),
    max_abs_diff = dif, match_3dp = !is.na(dif) && dif < 5e-4,
    t_total = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time())
  ))

  cat(sprintf("  waveform x %-9s q=%-4s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | pub %.3f/%.3f/%.3f/%.3f | %s | %.1fs\n",
              meth, q, out$m[["TPR"]], out$m[["TNR"]], out$m[["BA"]], out$m[["F2"]],
              pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]],
              if (!is.na(dif) && dif < 5e-4) "MATCH" else sprintf("diff %.3f", dif), wall))
}
cat("29_waveform_alpha.R: done\n")
