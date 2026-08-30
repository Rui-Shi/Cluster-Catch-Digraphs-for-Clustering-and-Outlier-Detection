# 16_wp0_su_highsmin.R -- one open question from the constant-S_min arm.
#
# SUN-MCCD reproduces under a constant S_min in [0.05, 0.0625]. SU-MCCD is
# consistent with the same constant but needs more of it (hepatitis matches at
# >=0.05, vertebral at >=0.0625). Glass and stamps SU-MCCD were only ever tested
# at or below 0.05, where they fail. If they reproduce above 0.0625, both
# shape-adaptive methods shared one constant; if not, they were run differently.
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

GRID <- c(0.0625, 0.070, 0.080, 0.100, 0.150, 0.200, 0.300)

for (ds in c("glass", "stamps", "vertebral")) {
  dat <- load_real_dataset(ds)
  pub <- published(ds, "SU-MCCD")
  cat(sprintf("\n=== %s SU-MCCD   published TPR %.3f TNR %.3f BA %.3f F2 %.3f\n",
              ds, pub[["TPR"]], pub[["TNR"]], pub[["BA"]], pub[["F2"]]))
  for (mc in GRID) {
    line <- tryCatch({
      res <- METHOD_REGISTRY[["SU-MCCD"]](dat$X, dat$d, dat$Y, min.cls = mc)
      ev  <- evaluate(dat$Y, res$score, 0.5)
      dif <- max(abs(round(ev, 3) - pub), na.rm = TRUE)
      sprintf("  S_min=%.4f  TPR %.3f  TNR %.3f  BA %.3f  F2 %.3f  maxdiff %.4f %s",
              mc, ev[["TPR"]], ev[["TNR"]], ev[["BA"]], ev[["F2"]], dif,
              if (dif < 0.0005) "<<< MATCH" else if (dif < 0.003) "<<< rounding-only" else "")
    }, error = function(e) sprintf("  S_min=%.4f  ERROR: %s", mc, conditionMessage(e)))
    cat(line, "\n")
  }
}

cat("\ndone\n")
