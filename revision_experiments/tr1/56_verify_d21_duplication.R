#!/usr/bin/env Rscript
# 56_verify_d21_duplication.R -- did the GENERATOR produce the duplicate NN
# quantile tables, or did a file copy?
#
# 21 pairs of NN quantile tables on disk are byte-identical (d = 6-9, 11-19,
# 21-28). The manuscript's withdrawn saturation claim rested on one of them.
# Knowing WHICH mechanism produced them matters for what else is suspect:
#
#   - If the generator ignores its `quant` argument, then every NN table ever
#     produced by it is at an unknown level, and the fault is systematic.
#   - If the generator honours `quant` and the duplicates came from copying
#     files, the fault is clerical, confined to the affected filenames, and
#     every table generated in a single-level run is sound.
#
# The two are distinguished by a design point that is easy to get wrong:
# running the generator twice at two levels produces two DIFFERENT random draw
# sets, so the outputs would differ even if `quant` were ignored entirely.
# That comparison proves nothing. The decisive test holds the draws FIXED and
# varies only the reduction.
#
# TEST A -- on-disk baseline. Record that the shipped d=21 pair is identical,
#   so the comparison below has a reference point.
#
# TEST B -- SAME DRAWS, two levels (decisive). One Monte Carlo run reduced at
#   0.99 and at 0.999 by 01i_nn_multiquant_table.R, whose reduction was
#   validated bitwise against the untouched original function. If these two
#   come out identical, the reduction conflates the levels and the fault is in
#   the generator. If they differ, the generator is sound.
#
# TEST C -- INDEPENDENT RUNS, mimicking production. Two separate invocations of
#   the ORIGINAL 01_gen_quantile_table.R at 0.99 and 0.999, exactly as the
#   production tables were made. These must differ (different draws AND
#   different levels). If they differ here while the shipped pair is identical,
#   the shipped duplication cannot have come from this code path.
#
# Interpreting B and C together:
#   B differs, C differs  -> generator sound; shipped duplicates are a FILE COPY
#   B identical           -> generator conflates levels; systematic fault
#   B differs, C identical-> would be bizarre; investigate before concluding
#
# Sizing, split by what each test needs:
#   TEST B is generated at the shipped n=5000. That is more than the identity
#   question strictly requires, but it makes the output a genuine drop-in
#   replacement, which TEST E below depends on: waveform has n=3443, and a
#   table built at n=1000 would trip nnccd.radi's silent clamp
#   (pmin(1:n, len), harness.R:157-158), reusing the last table entry for two
#   thirds of the points with no warning.
#   TEST C is generated at n=1000. Its only job is to show that two
#   independent runs differ, which does not depend on n, and at n=5000 it
#   would cost another two full generations.
#
# Read-only with respect to every shipped table. Writes only under
# R/NN-test_quantile_d21_regen/.

suppressMessages(library(here))

SHIPPED <- here::here("R/NN-test_quantile")
REGEN   <- here::here("R/NN-test_quantile_d21_regen")
OUT     <- here::here("revision_experiments/results/tr1/wp2a_d21_duplication.csv")

load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

cmp <- function(label, pa, pb) {
  if (!file.exists(pa) || !file.exists(pb)) {
    cat(sprintf("  %-34s MISSING (%s / %s)\n", label,
                file.exists(pa), file.exists(pb)))
    return(data.frame(test = label, ok = NA, identical_avg = NA, identical_med = NA,
                      max_abs_diff_avg = NA_real_, mean_abs_diff_avg = NA_real_,
                      n_avg = NA_integer_, stringsAsFactors = FALSE))
  }
  a <- load1(pa); b <- load1(pb)
  m <- min(length(a$average), length(b$average))
  ia <- identical(a$average, b$average); im <- identical(a$median, b$median)
  mx <- max(abs(a$average[1:m] - b$average[1:m]), na.rm = TRUE)
  mn <- mean(abs(a$average[1:m] - b$average[1:m]), na.rm = TRUE)
  cat(sprintf("  %-34s identical(avg)=%-5s identical(med)=%-5s  max|d|=%.6g  mean|d|=%.6g  n=%d\n",
              label, ia, im, mx, mn, m))
  data.frame(test = label, ok = TRUE, identical_avg = ia, identical_med = im,
             max_abs_diff_avg = mx, mean_abs_diff_avg = mn, n_avg = m,
             stringsAsFactors = FALSE)
}

rows <- list()

cat("=== TEST A: shipped d=21 pair (baseline) ===\n")
rows[[1]] <- cmp("A shipped 99 vs 999",
                 file.path(SHIPPED, "NN-test-simul_21d_99%.RData"),
                 file.path(SHIPPED, "NN-test-simul_21d_999%.RData"))

cat("\n=== TEST B: same draws, two reductions (DECISIVE) ===\n")
rows[[2]] <- cmp("B same-draws 99 vs 999",
                 file.path(REGEN, "samedraws", "NN-test-simul_21d_99%.RData"),
                 file.path(REGEN, "samedraws", "NN-test-simul_21d_999%.RData"))

cat("\n=== TEST C: two independent generator runs (production-style) ===\n")
rows[[3]] <- cmp("C independent 99 vs 999",
                 file.path(REGEN, "indep99",  "NN-test-simul_21d_99%.RData"),
                 file.path(REGEN, "indep999", "NN-test-simul_21d_999%.RData"))

# A same-level control: if two independent runs at the SAME level differ only
# by Monte Carlo noise, that calibrates how big the level effect in C is.
cat("\n=== TEST D: control, independent run vs same-draws run at 0.99 ===\n")
rows[[4]] <- cmp("D indep99 vs samedraws99",
                 file.path(REGEN, "indep99",   "NN-test-simul_21d_99%.RData"),
                 file.path(REGEN, "samedraws", "NN-test-simul_21d_99%.RData"))

res <- do.call(rbind, rows)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write.csv(res, OUT, row.names = FALSE)

cat("\n=== VERDICT ===\n")
B <- res[res$test == "B same-draws 99 vs 999", ]
C <- res[res$test == "C independent 99 vs 999", ]
if (isTRUE(is.na(B$ok)) || isTRUE(is.na(C$ok))) {
  cat("  Tests B and/or C have not been generated yet -- run the generation step first.\n")
} else if (isTRUE(B$identical_avg) && isTRUE(B$identical_med)) {
  cat("  GENERATOR FAULT: reducing one draw set at 0.99 and 0.999 gave identical\n")
  cat("  output. The quant argument is not reaching the reduction. Every NN table\n")
  cat("  produced by this path is at an unknown level.\n")
} else if (isTRUE(C$identical_avg)) {
  cat("  ANOMALY: same-draws reductions differ but two independent runs agree.\n")
  cat("  Investigate before drawing any conclusion.\n")
} else {
  cat("  GENERATOR SOUND: same-draws reductions differ (max|d| = ")
  cat(sprintf("%.6g", B$max_abs_diff_avg))
  cat("), and two independent\n  runs differ (max|d| = ")
  cat(sprintf("%.6g", C$max_abs_diff_avg))
  cat(").\n  The shipped d=21 pair being byte-identical therefore cannot have come from\n")
  cat("  this code path -- the duplication is a FILE-LEVEL accident (one run saved\n")
  cat("  or copied under two names), not a systematic generator defect.\n")
  cat("  Consequence: tables written by single-level runs are sound; only the\n")
  cat("  affected filenames misstate their level.\n")
}
cat(sprintf("\nNote: regenerated at n=1000, while the shipped d=21 tables are n=5000.\n"))
cat("These are not drop-in replacements; the identity question does not depend on n.\n")
cat(sprintf("\nwrote %s\ndone\n", OUT))
