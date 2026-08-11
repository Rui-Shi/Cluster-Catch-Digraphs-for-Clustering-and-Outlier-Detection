#!/usr/bin/env Rscript
# 71_regen_nnd_alpha001.R -- generate GENUINE alpha = 0.1% NND quantile tables
# at the five dimensions where SUN-MCCD's published real-data results were
# produced from a table that actually holds 1%.
#
# WHY THIS RUN EXISTS
# -------------------
# WP2_RESULTS section 9 recovered the alpha level behind every shipped NN
# quantile table. At d = 12, 16, 18, 19 and 21 the file named "999" (0.1%) is
# byte-identical to its "99" (1%) sibling -- one Monte Carlo run saved under two
# names -- and the recovered content is 1%. The paper states SUN-MCCD uses
# alpha = 0.1% for d >= 10, so five real-data rows (vowels, PenDigits,
# lymphography, hepatitis, waveform) plus UN-MCCD on waveform were computed at a
# level the manuscript does not claim.
#
# The chosen fix keeps the manuscript's stated schedule and repairs the tables,
# rather than restating the schedule to match the tables. That leaves the
# simulation study, every figure, and Table tab:alpha_real untouched; only the
# affected real-data rows are recomputed afterwards.
#
# WHY n IS THE DATA SET'S OWN n, NOT 5000
# ---------------------------------------
# The shipped tables all carry 5000 entries, but nothing in the generator makes
# column k depend on the total n. Both engines build a FRESH k-point sample for
# every subsample size k (01i_nn_multiquant_table.R:112-119 for orig_list,
# :132-138 for fast_stream), so column k is a function of k and d alone and the
# total n only decides how many columns get produced. A table generated to n is
# therefore identical in distribution to an n=5000 table truncated at n.
#
# That matters because the cost is O(sum_k k^2 d) ~ n^3 d / 3. Generating all
# five at n=5000 is roughly 209 h on 22 cores; generating each only as far as
# its data set actually needs is about 28 h. The lookup can never run off the
# end, because a table covering n entries covers every possible cluster size in
# an n-point data set.
#
#   d=12  vowels        n=1452
#   d=16  PenDigits     n=3200
#   d=18  lymphography  n=148
#   d=19  hepatitis     n=74
#   d=21  waveform      n=3443
#
# niter = 10000 matches the original generators (R/NN-test_quantile/10d_99%.R,
# 20d_999%.R). This is not a detail to economise on: at niter = 250 -- the depth
# the earlier d=21 attempt reached -- quantile(x, 0.001) IS the sample minimum,
# whose expected CDF position is 1/251 ~ 0.004, so it would be a 0.4% table
# wearing a 0.1% name. Reproducing the very defect being repaired.
#
# DURABILITY. CLAUDE.md records that this volume is an external USB-NVMe
# enclosure that has dropped off the bus under sustained write load, and
# requires long runs to persist per cell rather than buffer to the end. Each
# dimension is split into chunks sized to land every ~20 minutes; every chunk's
# raw draws hit disk before the next starts, and completed chunks are skipped on
# re-invocation. Draws are i.i.d. across iterations, so chunking and pooling is
# statistically identical to one long run. The reduction to quantiles happens
# ONCE on the pooled draws -- pooling quantiles would not be equivalent, since
# quantile(rbind(A,B)) is not the mean of the two quantiles.
#
# NOTHING IS INSTALLED BY THIS SCRIPT. Output goes to a staging directory. The
# tables only reach R/NN-test_quantile after 72_verify_regen999.R confirms each
# one reads back at its nominal level, using the same estimator that recovered
# the mislabelling in the first place.
#
# MEASURED COST (22 cores, this machine, 2026-08-11). An earlier estimate here
# was calibrated on a comment in 58_gen_d21_chunked.R quoting "roughly sixteen
# hours" for d=21/n=5000/niter=2000 -- that figure was an a-priori guess written
# before the run, and the run only ever completed 250 iterations, so it was
# never a measurement. Timing six iterations per dimension gives the real
# serial cost of one iteration:
#
#   d=12  n=1452    33 s/iter  ->  10000 iters / 22 cores =   4.2 h
#   d=16  n=3200   396 s/iter  ->                             50   h
#   d=21  n=3443  ~615 s/iter  ->                             78   h
#   d=18  n=148, d=19 n=74                                  < 5 min
#
# There is no shortcut on n: UN_CCD.R:241 slices simul$average[1:n] with
# n = nrow(dx), the FULL data matrix handed to the CCD, so the table needs one
# entry per data point. A shorter table yields NA past its end, not a clamp.
#
# Usage:
#   Rscript revision_experiments/71_regen_nnd_alpha001.R [cores] [niter] [dims]
#   defaults: 22 cores, 10000 iterations, dims = 19,18,12 (the cheap block)
#   e.g.      Rscript ... 22 10000 16      # the 50 h dimension
#             Rscript ... 22 10000 21      # the 78 h dimension
#
# Resume: safe to re-invoke at any time; finished chunks are skipped.

suppressMessages(library(here))
source(here::here("revision_experiments", "01i_nn_multiquant_table.R"))

args   <- commandArgs(trailingOnly = TRUE)
CORES  <- if (length(args) >= 1) as.integer(args[[1]]) else 22L
NITER  <- if (length(args) >= 2) as.integer(args[[2]]) else 10000L
DIMS   <- if (length(args) >= 3) as.integer(strsplit(args[[3]], ",")[[1]]) else c(19L, 18L, 12L)

ENGINE  <- "fast_stream"
QUANTS  <- c(0.95, 0.99, 0.999)   # same draws reduced at every level we may cite
SEEDBAS <- 20260811L
TARGET_CHUNK_SEC <- 1200          # aim for a chunk landing every ~20 min

OUTROOT <- here::here("R/NN-test_quantile_regen999")
PROG    <- file.path(OUTROOT, "PROGRESS.tsv")
dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)

# d, n, and the data set that fixes n. Ordered cheapest first so the early
# dimensions are on disk within minutes and a misconfiguration surfaces then
# rather than sixteen hours in.
JOBS <- data.frame(
  d       = c(19L,          18L,             12L,      16L,         21L),
  n       = c(74L,          148L,            1452L,    3200L,       3443L),
  dataset = c("hepatitis",  "lymphography",  "vowels", "PenDigits", "waveform"),
  stringsAsFactors = FALSE
)
JOBS <- JOBS[JOBS$d %in% DIMS, , drop = FALSE]
stopifnot(nrow(JOBS) > 0)

# ---------------------------------------------------------------------------
# progress log -- append-only, one line per event, machine-readable.
# This file is the monitor's entire interface to the run.
# ---------------------------------------------------------------------------
say <- function(fmt, ...) {
  line <- sprintf("%s\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(fmt, ...))
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = PROG, append = TRUE)
  flush.console()
}

#' Core-seconds for ONE iteration, from six-iteration timings taken on this
#' machine on 2026-08-11 (22 cores, so six iterations run concurrently and the
#' wall time of that batch is one iteration's serial cost). Measured values are
#' used directly; anything else falls back to the n^3 * d model anchored on the
#' d=16 measurement, which reproduces the d=12 point to within 20%.
MEASURED <- c("12" = 33, "16" = 396, "18" = 0.7, "19" = 0.15)
est_core_sec_per_iter <- function(d, n) {
  k <- as.character(d)
  if (k %in% names(MEASURED)) return(unname(MEASURED[k]))
  396 * (n / 3200)^3 * (d / 16)
}

say("RUNSTART\tcores=%d\tniter=%d\tengine=%s\tdims=%s", CORES, NITER, ENGINE,
    paste(JOBS$d, collapse = ","))

total_est <- sum(vapply(seq_len(nrow(JOBS)), function(i)
  est_core_sec_per_iter(JOBS$d[i], JOBS$n[i]) * NITER / CORES, numeric(1)))
say("ESTIMATE\ttotal_wall_hours=%.1f", total_est / 3600)

for (j in seq_len(nrow(JOBS))) {
  D <- JOBS$d[j]; N <- JOBS$n[j]; DS <- JOBS$dataset[j]

  outdir <- file.path(OUTROOT, sprintf("%dd", D))
  chkdir <- file.path(outdir, "chunks")
  dir.create(chkdir, recursive = TRUE, showWarnings = FALSE)

  final <- file.path(outdir, sprintf("NN-test-simul_%dd_999%%.RData", D))
  if (file.exists(final)) { say("DIMSKIP\td=%d\talready complete", D); next }

  est_wall <- est_core_sec_per_iter(D, N) * NITER / CORES
  nchunk   <- max(1L, min(40L, as.integer(round(est_wall / TARGET_CHUNK_SEC))))
  perchk   <- as.integer(ceiling(NITER / nchunk))
  nchunk   <- as.integer(ceiling(NITER / perchk))   # re-derive after rounding

  say("DIMSTART\td=%d\tn=%d\tdataset=%s\tchunks=%dx%d\test_wall_h=%.2f",
      D, N, DS, nchunk, perchk, est_wall / 3600)

  chunk_file <- function(i) file.path(chkdir, sprintf("draws_chunk%03d.RData", i))

  for (i in seq_len(nchunk)) {
    if (file.exists(chunk_file(i))) { say("CHUNKSKIP\td=%d\tchunk=%d/%d", D, i, nchunk); next }
    # last chunk absorbs the remainder so the pooled total is exactly NITER
    this_iter <- if (i == nchunk) NITER - perchk * (nchunk - 1L) else perchk
    if (this_iter <= 0L) next

    t0 <- Sys.time()
    draws <- nn_mc_draws(n = N, d = D, niter = this_iter, cores = CORES,
                         seed = SEEDBAS + 1000L * D + i, engine = ENGINE)
    NN.dist.ave.mat <- draws$ave
    NN.dist.med.mat <- draws$med
    nn_meta <- draws$meta
    save(NN.dist.ave.mat, NN.dist.med.mat, nn_meta, file = chunk_file(i))
    rm(draws, NN.dist.ave.mat, NN.dist.med.mat); gc(verbose = FALSE)

    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    say("CHUNK\td=%d\tchunk=%d/%d\titer=%d\tsecs=%.1f\teta_h=%.2f",
        D, i, nchunk, this_iter, el, el * (nchunk - i) / 3600)
  }

  # ---- pool, reduce once, write -------------------------------------------
  say("POOL\td=%d\tloading %d chunks", D, nchunk)
  ave <- matrix(0, nrow = NITER, ncol = N)
  med <- matrix(0, nrow = NITER, ncol = N)
  row <- 0L
  for (i in seq_len(nchunk)) {
    e <- new.env(); load(chunk_file(i), envir = e)
    A <- get("NN.dist.ave.mat", envir = e); M <- get("NN.dist.med.mat", envir = e)
    idx <- (row + 1L):(row + nrow(A))
    ave[idx, ] <- A; med[idx, ] <- M
    row <- row + nrow(A)
    rm(e, A, M); gc(verbose = FALSE)
  }
  if (row != NITER) { say("ERROR\td=%d\tpooled %d rows, expected %d", D, row, NITER); quit(status = 1) }

  pooled <- list(ave = ave, med = med,
                 meta = list(n = N, d = D, niter = NITER, cores = CORES,
                             seed = SEEDBAS + 1000L * D, engine = ENGINE,
                             chunks = nchunk, dataset = DS,
                             generated = format(Sys.time()),
                             R_version = R.version.string,
                             provenance = "71_regen_nnd_alpha001.R"))
  w <- write_nn_tables(pooled, QUANTS, outdir, overwrite = FALSE)
  rm(ave, med, pooled); gc(verbose = FALSE)

  z <- w$zero_frac
  say("DIMDONE\td=%d\tn=%d\tfiles=%d\tzero_frac_999=%.4f",
      D, N, length(w$files), z$frac_zero_average_excl_pos1[z$token == "999"])
}

say("ALLDONE\tstaging=%s", OUTROOT)
say("NEXT\tRscript revision_experiments/72_verify_regen999.R  (verify before installing)")
