# 14_wp0_probe.R -- two targeted WP0 hypotheses, small datasets only.
#
# H1 (alpha interpolation): the manuscript fixes SUN-MCCD's alpha only at
#     d in {2,3,5,10,20,50,100}. Our resolver forces the 999% table at d>=10,
#     which is an interpolation the paper never states. NN tables exist at BOTH
#     99% and 999% for d=12 and d=21, and every waveform cell under-flags
#     (TPR low, TNR high) -- the signature of an alpha that is too stringent.
#     Test: rerun SUN-MCCD at 99% and see whether it moves toward published.
#
# H2 (constant S_min): no single multiple of the contamination rate reproduces
#     the published rows -- glass SUN-MCCD needs n0/n, vertebral SUN-MCCD needs
#     0.5*n0/n. Values that do reproduce cluster in 0.042-0.063, consistent with
#     a fixed default near 0.05. Test a constant against the contamination-derived
#     readings on the datasets that discriminate.
#
# Read-only with respect to every existing file: prints to stdout, writes nothing.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

TRUTH <- read.csv(here::here("revision_experiments", "published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)

published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

probe <- function(ds, meth, label, ...) {
  out <- tryCatch({
    dat <- load_real_dataset(ds)
    t0  <- Sys.time()
    res <- METHOD_REGISTRY[[meth]](dat$X, dat$d, dat$Y, ...)
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ev  <- evaluate(dat$Y, res$score, 0.5)
    pub <- published(ds, meth)
    dif <- max(abs(round(ev, 3) - pub), na.rm = TRUE)
    cat(sprintf("%-10s %-9s %-26s q=%-4s TPR %.3f/%.3f  TNR %.3f/%.3f  BA %.3f/%.3f  F2 %.3f/%.3f  maxdiff %.4f %s  [%.0fs]\n",
                ds, meth, label,
                if (!is.null(res$quant_label)) res$quant_label else "?",
                ev[["TPR"]], pub[["TPR"]], ev[["TNR"]], pub[["TNR"]],
                ev[["BA"]],  pub[["BA"]],  ev[["F2"]],  pub[["F2"]],
                dif, if (dif < 0.0005) "MATCH" else "", el))
    invisible(dif)
  }, error = function(e) {
    cat(sprintf("%-10s %-9s %-26s ERROR: %s\n", ds, meth, label, conditionMessage(e)))
    invisible(NA)
  })
  out
}

half <- function(ds) { d <- load_real_dataset(ds); 0.5 * sum(d$Y == 0) / d$n }
full <- function(ds) { d <- load_real_dataset(ds); sum(d$Y == 0) / d$n }

cat("\n=== H1: alpha interpolation, SUN-MCCD at 99% vs the resolver's 999% ===\n")
cat("--- vowels (d=12, both NN tables exist); resolver picks 999\n")
for (mc in c(half("vowels"), full("vowels"), 0.05)) {
  probe("vowels", "SUN-MCCD", sprintf("min.cls=%.4f q=999", mc), min.cls = mc)
  probe("vowels", "SUN-MCCD", sprintf("min.cls=%.4f q=99",  mc), min.cls = mc, quant = "99")
}

cat("\n=== H2: constant S_min = 0.05 vs contamination-derived ===\n")
for (ds in c("glass", "vertebral", "stamps")) {
  for (m in c("SU-MCCD", "SUN-MCCD")) {
    probe(ds, m, sprintf("half=%.4f", half(ds)), min.cls = half(ds))
    probe(ds, m, sprintf("full=%.4f", full(ds)), min.cls = full(ds))
    probe(ds, m, "const=0.0500",                 min.cls = 0.05)
  }
}

cat("\ndone\n")
