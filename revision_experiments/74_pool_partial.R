#!/usr/bin/env Rscript
# 74_pool_partial.R -- build the quantile tables from however many chunks of a
# 71 run have landed so far, instead of waiting for all of them.
#
# WHY THIS TURNS A DECISION INTO A NON-DECISION. 71 writes raw draws per chunk
# and reduces only at the end, so its niter is fixed when the run starts.
# Choosing niter up front means choosing between a table that is ready in a day
# and one that is ready in six -- with a 2026-08-24 deadline, that is a real
# bet. But the chunks are i.i.d. blocks of draws: whatever has landed at any
# moment IS a valid Monte Carlo run of that size. Pooling on demand means the
# run can simply be stopped when the table is good enough, and the choice made
# with the evidence in hand rather than in advance.
#
# 73 measured what the choice costs at d=12: two INDEPENDENT 2500-iteration
# runs give 0.1% tables differing by a median 0.11%, and two 5000-iteration
# runs by 0.08%, against the 2% tolerance at which 72's T1 accepts a table --
# and against the 0.04-0.18% by which the freshly generated 1% table differs
# from the shipped one. The tail of this statistic is well enough behaved that
# a quarter of the iterations costs almost nothing.
#
# The pooled niter is recorded in the table's own metadata and in the filename
# of the draws file, so a table can never be mistaken for a deeper run than it
# was. That is the same discipline the whole work package exists to restore.
#
# Usage:
#   Rscript revision_experiments/74_pool_partial.R <d> [max_chunks]
#   e.g.  Rscript revision_experiments/74_pool_partial.R 21
#         Rscript revision_experiments/74_pool_partial.R 21 10   # cap at 10

suppressMessages(library(here))
source(here::here("revision_experiments", "01i_nn_multiquant_table.R"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript revision_experiments/74_pool_partial.R <d> [max_chunks]")
D    <- as.integer(args[[1]])
MAXC <- if (length(args) >= 2) as.integer(args[[2]]) else NA_integer_

QUANTS <- c(0.95, 0.99, 0.999)
OUTDIR <- here::here("R/NN-test_quantile_regen999", sprintf("%dd", D))
CHKDIR <- file.path(OUTDIR, "chunks")

fs <- sort(list.files(CHKDIR, pattern = "^draws_chunk[0-9]+\\.RData$", full.names = TRUE))
if (!length(fs)) stop("no chunks on disk for d=", D, " in ", CHKDIR)
if (!is.na(MAXC)) fs <- head(fs, MAXC)
cat(sprintf("d=%d: pooling %d chunk file(s)\n", D, length(fs)))

# A chunk still being written would deserialise short or fail; load defensively
# and drop anything unreadable rather than silently pooling a truncated block.
blocks <- list()
for (p in fs) {
  ok <- tryCatch({
    e <- new.env(); load(p, envir = e)
    list(ave = get("NN.dist.ave.mat", envir = e),
         med = get("NN.dist.med.mat", envir = e))
  }, error = function(err) { cat(sprintf("  SKIP %s (%s)\n", basename(p), conditionMessage(err))); NULL })
  if (!is.null(ok)) blocks[[length(blocks) + 1]] <- ok
}
if (!length(blocks)) stop("no readable chunks")

ave <- do.call(rbind, lapply(blocks, `[[`, "ave"))
med <- do.call(rbind, lapply(blocks, `[[`, "med"))
rm(blocks); gc(verbose = FALSE)
NITER <- nrow(ave); N <- ncol(ave)
cat(sprintf("  pooled %d iterations x %d subsample sizes\n", NITER, N))

pooled <- list(ave = ave, med = med,
               meta = list(n = N, d = D, niter = NITER, engine = "fast_stream",
                           chunks_pooled = length(fs), partial = TRUE,
                           generated = format(Sys.time()),
                           R_version = R.version.string,
                           provenance = "74_pool_partial.R over 71's chunks"))

# write_nn_tables refuses to clobber, so a re-pool at a deeper niter goes to a
# fresh directory keyed by the iteration count -- previous tables stay put and
# stay comparable.
dest <- file.path(OUTDIR, sprintf("pooled_n%d", NITER))
w <- write_nn_tables(pooled, QUANTS, dest, overwrite = FALSE)
cat(sprintf("  wrote %d files to %s\n", length(w$files), dest))
print(w$zero_frac[, c("token", "alpha", "frac_zero_average_excl_pos1", "min_nonzero_average")],
      row.names = FALSE, digits = 4)
cat("\nnext: point 72_verify_regen999.R at this directory before installing\n")
