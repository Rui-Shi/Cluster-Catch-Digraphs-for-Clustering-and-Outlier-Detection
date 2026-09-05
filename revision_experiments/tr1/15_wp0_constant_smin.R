# 15_wp0_constant_smin.R -- extend the constant-S_min arm.
#
# 14_wp0_probe.R showed a constant S_min = 0.05 reproduces SUN-MCCD on glass,
# vertebral and stamps, where neither 0.5*(n0/n) nor n0/n does. This extends the
# arm to the datasets that probe did not cover, and sweeps a small grid so the
# plateau can be bracketed rather than asserted from a single value.
#
# Fast datasets only (hepatitis, ecoli, vowels). waveform and wilt are deferred
# so this does not contend with the gate run currently holding cores.
#
# Writes nothing. Prints only.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

TRUTH <- read.csv(here::here("revision_experiments", "published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)

published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

GRID <- c(0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.100)

for (ds in c("hepatitis", "ecoli", "vowels")) {
  dat <- load_real_dataset(ds)
  n0  <- sum(dat$Y == 0)
  cat(sprintf("\n=== %s  n=%d d=%d n0=%d   half=%.4f full=%.4f\n",
              ds, dat$n, dat$d, n0, 0.5 * n0 / dat$n, n0 / dat$n))
  for (meth in c("SU-MCCD", "SUN-MCCD")) {
    pub <- published(ds, meth)
    cat(sprintf("  %-9s published TPR %.3f TNR %.3f BA %.3f F2 %.3f\n",
                meth, pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]]))
    for (mc in GRID) {
      r <- tryCatch({
        res <- METHOD_REGISTRY[[meth]](dat$X, dat$d, dat$Y, min.cls = mc)
        ev  <- evaluate(dat$Y, res$score, 0.5)
        dif <- max(abs(round(ev, 3) - pub), na.rm = TRUE)
        sprintf("    S_min=%.3f  TPR %.3f  TNR %.3f  BA %.3f  F2 %.3f  maxdiff %.4f %s",
                mc, ev[["TPR"]], ev[["TNR"]], ev[["BA"]], ev[["F2"]], dif,
                if (dif < 0.0005) "<<< MATCH" else if (dif < 0.003) "<<< rounding-only" else "")
      }, error = function(e) sprintf("    S_min=%.3f  ERROR: %s", mc, conditionMessage(e)))
      cat(r, "\n")
    }
  }
}

cat("\ndone\n")
