#!/usr/bin/env Rscript
# 58_gen_d21_chunked.R -- durable, resumable generation of the d=21 n=5000 NN
# quantile pair.
#
# WHY THIS EXISTS. 01i_nn_multiquant_table.R runs the whole Monte Carlo in one
# parLapply and calls write_nn_tables() only after every iteration finishes. At
# d=21, n=5000, niter=2000 that is roughly sixteen hours during which nothing at
# all is on disk. Two problems:
#
#   1. Durability. The volume this repo sits on is an external USB-NVMe
#      enclosure that dropped off the bus mid-write earlier tonight, logging
#      paging I/O errors and a Delayed Write Failed against $Mft. A sixteen-hour
#      all-or-nothing run on that hardware is a bad bet, and it contradicts the
#      rule recorded in CLAUDE.md that long runs must persist incrementally.
#   2. Observability. With no intermediate output, a wedged job and a working
#      job look identical for sixteen hours.
#
# Chunking fixes both. The draws are i.i.d. across iterations, so splitting
# niter into chunks with distinct seeds and pooling them afterwards is
# statistically identical to one long run -- this is a durability change, not a
# methodological one. The reduction to quantiles happens ONCE, on the pooled
# draw matrices, which is what keeps the two alpha levels a genuine same-draws
# pair (the whole point of the test in 56_verify_d21_duplication.R).
#
# Pooling quantiles would NOT be equivalent: quantile(rbind(A,B)) is not the
# mean of quantile(A) and quantile(B). Hence chunks store raw draws, never
# reduced tables.
#
# 01i is sourced, not modified. Its .nn_mq_is_main() guard keys on the --file=
# basename, so sourcing does not trigger its CLI.
#
# Usage:
#   Rscript revision_experiments/58_gen_d21_chunked.R [n_chunks] [iter_per_chunk] [cores]
#   defaults: 16 chunks x 125 iterations = 2000, 14 cores
#
# Resume: completed chunk files are skipped. Safe to re-invoke after any failure.

suppressMessages(library(here))
source(here::here("revision_experiments", "01i_nn_multiquant_table.R"))

args    <- commandArgs(trailingOnly = TRUE)
NCHUNK  <- if (length(args) >= 1) as.integer(args[[1]]) else 16L
PERCHK  <- if (length(args) >= 2) as.integer(args[[2]]) else 125L
CORES   <- if (length(args) >= 3) as.integer(args[[3]]) else 14L

D       <- 21L
N       <- 5000L
ENGINE  <- "fast_stream"
QUANTS  <- c(0.99, 0.999)
SEEDBAS <- 20260810L

OUTDIR  <- here::here("R/NN-test_quantile_d21_regen/samedraws")
CHKDIR  <- file.path(OUTDIR, "chunks")
dir.create(CHKDIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== 58_gen_d21_chunked ===\nd=%d n=%d chunks=%d x %d = %d iterations  cores=%d engine=%s\noutdir=%s\n",
            D, N, NCHUNK, PERCHK, NCHUNK * PERCHK, CORES, ENGINE, OUTDIR))
flush.console()

chunk_file <- function(i) file.path(CHKDIR, sprintf("draws_chunk%02d.RData", i))

# ---- generate any missing chunks ------------------------------------------
for (i in seq_len(NCHUNK)) {
  f <- chunk_file(i)
  if (file.exists(f)) {
    cat(sprintf("CHUNK_SKIP %d/%d (already on disk)\n", i, NCHUNK)); flush.console()
    next
  }
  t0 <- Sys.time()
  dr <- nn_mc_draws(n = N, d = D, niter = PERCHK, cores = CORES,
                    seed = SEEDBAS + i, engine = ENGINE)
  ave <- dr$ave; med <- dr$med; meta <- dr$meta
  # write to a temporary name first, then rename: a torn write leaves a .part
  # file rather than a corrupt chunk that the resume logic would trust.
  tmp <- paste0(f, ".part")
  save(ave, med, meta, file = tmp)
  file.rename(tmp, f)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("CHUNK_DONE %d/%d iters=%d sec=%.1f  (%.1f s/iter)\n",
              i, NCHUNK, PERCHK, el, el / PERCHK)); flush.console()
  rm(dr, ave, med); gc(verbose = FALSE)
}

# ---- pool and reduce -------------------------------------------------------
cat("POOLING chunks\n"); flush.console()
ave_l <- vector("list", NCHUNK); med_l <- vector("list", NCHUNK)
for (i in seq_len(NCHUNK)) {
  e <- new.env(); load(chunk_file(i), envir = e)
  ave_l[[i]] <- get("ave", envir = e); med_l[[i]] <- get("med", envir = e)
}
ave_all <- do.call(rbind, ave_l); med_all <- do.call(rbind, med_l)
rm(ave_l, med_l); gc(verbose = FALSE)

stopifnot("pooled row count must equal total iterations" =
            nrow(ave_all) == NCHUNK * PERCHK,
          "pooled column count must equal n" = ncol(ave_all) == N,
          "ave and med must have the same shape" =
            identical(dim(ave_all), dim(med_all)))
cat(sprintf("pooled draws: %d x %d\n", nrow(ave_all), ncol(ave_all)))

draws <- list(ave = ave_all, med = med_all,
              meta = list(n = N, d = D, niter = NCHUNK * PERCHK,
                          engine = ENGINE, seed_base = SEEDBAS,
                          chunks = NCHUNK, iter_per_chunk = PERCHK,
                          pooled = TRUE))

res <- write_nn_tables(draws, QUANTS, OUTDIR, overwrite = FALSE)
cat("WROTE:\n"); for (f in res$files) cat(sprintf("  %s\n", f))

# ---- the question this whole run exists to answer --------------------------
tabs <- res$tables
nm   <- names(tabs)
if (length(nm) >= 2) {
  a <- tabs[[nm[1]]]; b <- tabs[[nm[2]]]
  same <- identical(a$average, b$average) && identical(a$median, b$median)
  cat(sprintf("\nSAME_DRAWS_TWO_LEVELS identical=%s  max|d avg|=%.6g\n",
              same, max(abs(a$average - b$average))))
  cat(sprintf("VERDICT: %s\n", if (same)
    "GENERATOR CONFLATES LEVELS -- quant is not reaching the reduction" else
    "GENERATOR SOUND -- the shipped duplicate pair cannot have come from this path"))
}
cat("ALL_CHUNKS_COMPLETE\n")
