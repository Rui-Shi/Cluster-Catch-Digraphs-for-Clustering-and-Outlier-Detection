#!/usr/bin/env Rscript
# 28_regen_final.R -- the final regeneration pass at a single fixed S_min.
#
#   Rscript 28_regen_final.R <datasets> <methods> [s_min]
#
# 25_smin_match.R recovered S_min = 0.0625 as the unique value in [0.01, 0.20]
# reproducing every reproducible published cell; vertebral pins it exactly
# (round(S_min*240) must equal 15, so S_min = 15/240 = 0.0625). A single fixed
# constant is also what WP2c needs: it uses no label information, so it answers
# the fairness objection three reviewers raised against deriving S_min from the
# contamination rate.
#
# This script therefore runs the four proposed detectors at ONE documented
# S_min for every dataset -- the configuration that goes in the paper.
#
# Checkpointed on (dataset, method, s_min): a restart skips completed cells,
# so this is safe to re-run after an interruption.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

args     <- commandArgs(trailingOnly = TRUE)
DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else
              c("hepatitis", "glass", "vertebral", "ecoli", "stamps", "vowels", "waveform", "wilt")
METHODS  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else
              c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN    <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 0.0625

OUT_CSV <- Sys.getenv("REGEN_FINAL_OUT_CSV",
                      here::here("revision_experiments/results/regen_final.csv"))

stopifnot(all(METHODS %in% names(METHOD_REGISTRY)))
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")

TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)
published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  if (nrow(s) == 0) return(setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")))
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

cat(sprintf("28_regen_final.R: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("28_regen_final.R: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("28_regen_final.R: S_min    = %g   output = %s\n", S_MIN, OUT_CSV))

for (ds in DATASETS) {
  dat <- load_real_dataset(ds)
  for (meth in METHODS) {
    # min.cls is a no-op for the uniform-coverage pair; record it as NA so the
    # row is not duplicated across S_min values.
    smin_cell <- if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_
    keys <- c(dataset = ds, method = meth, s_min = as.character(smin_cell))
    if (isTRUE(has_result(OUT_CSV, keys))) {
      cat(sprintf("[skip] %s x %s (S_min=%s)\n", ds, meth, smin_cell)); next
    }

    t0  <- Sys.time()
    out <- tryCatch({
      res <- if (meth %in% MIN_CLS_METHODS) {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y, min.cls = S_MIN)
      } else {
        METHOD_REGISTRY[[meth]](X = dat$X, d = dat$d, Y = dat$Y)
      }
      list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
           q = if (!is.null(res$quant_label)) res$quant_label else NA_character_,
           status = "ok", note = NA_character_)
    }, error = function(e) {
      list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
           q = NA_character_, status = "error", note = conditionMessage(e))
    })
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    pub <- published(ds, meth)
    dif <- if (any(is.na(out$m)) || any(is.na(pub))) NA_real_ else max(abs(round(out$m, 3) - pub))

    append_result(OUT_CSV, list(
      dataset = ds, method = meth, s_min = smin_cell,
      n = dat$n, d = dat$d,
      TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
      BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
      published_TPR = unname(pub[["TPR"]]), published_TNR = unname(pub[["TNR"]]),
      published_BA  = unname(pub[["BA"]]),  published_F2  = unname(pub[["F2"]]),
      max_abs_diff = dif, match_3dp = !is.na(dif) && dif < 5e-4,
      quant_label = out$q, t_total = wall,
      status = out$status, note = out$note, timestamp = format(Sys.time())
    ))

    cat(sprintf("  %-10s x %-9s S_min=%-7s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | pub %.3f/%.3f/%.3f/%.3f | %s | q=%s | %.1fs\n",
                ds, meth, format(smin_cell),
                out$m[["TPR"]], out$m[["TNR"]], out$m[["BA"]], out$m[["F2"]],
                pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]],
                if (!is.na(dif) && dif < 5e-4) "MATCH" else sprintf("diff %.3f", dif),
                out$q, wall))
  }
}
cat("28_regen_final.R: done\n")
