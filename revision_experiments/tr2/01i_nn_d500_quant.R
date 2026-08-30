#!/usr/bin/env Rscript
# revision_experiments/01h_nn_d500_quant.R -- NN(UN-CCD) quantile table for
# the new WP4 grid-2 cell d=500 (n=500). See revision_experiments/
# 04_wp4_runtime.R and its FLAGGED/SKIPPED_NO_TABLE docs: no
# R/NN-test_quantile/NN-test-simul_500d_999%.RData table exists on disk
# (confirmed by a throwaway inventory read on 2026-07-19 -- every other
# table the redesigned grids need already exists: NN d=5/10 at n<=2000,
# n<=5000-row tables; NN/RK d=50/100 at n<=500, 1000-row tables; RK d=5/10
# with >=5000 rows. Only NN d=500 is missing.)
#
# WHY A FULL BAND, NOT KNOTS: 01g_nn_sizes_quant.R (which extended the
# Musk/Speech NN tables past their existing n=1000/5000 entries) could get
# away with 6 interpolation knots because it was SPLICING onto an already-
# dense existing table and the envelope is documented near-flat over that
# high-n band. Here there is no existing d=500 table to splice onto and no
# evidence the envelope is flat over the LOW end of the band (small k is
# exactly where NN-distance quantiles change fastest, k^(-1/d) is far from
# its asymptote). n=500 is also small enough that computing every integer
# size is tractable (see runtime estimate below), so this script generates
# the FULL band 1..500, no knots, no interpolation.
#
# MACHINERY: reuses 01g_nn_sizes_quant.R's own approach verbatim --
# rep_one_size() draws ONE realization of k points via the ORIGINAL
# rpoisball.unit (R/ccds/NN_Dist_Est.R, NOT the 01e_nn_fast.R override; the
# task brief for this script specifically calls for reusing rpoisball.unit),
# and NNDestP.sizes.lower.quant()'s parallel-over-reps-per-size pattern with
# clusterSetRNGStream(cl, seed) for reproducibility is copied unmodified
# (only the `sizes` argument differs: 2:500 here instead of a knot vector).
# k=1 is NOT simulated (a lone point has no defined NN distance) -- entry 1
# of both average/median is hardcoded to 0, exactly mirroring 01e_nn_fast.R's
# `c(0, NN.dist.temp[1, ])` convention and NNDestP.simpois.lower.quant's own
# (R/ccds/NN_Dist_Est.R) treatment of the same edge case.
#
# HONEST RUNTIME ESTIMATE (measured 2026-07-19, single-threaded benchmark,
# NOT this script -- see task report): rep_one_size() reuses the ORIGINAL
# rpoisball.unit, whose MASS::mvrnorm call pays a fresh eigen(500x500)
# decomposition on EVERY call (~0.05-0.1 s, roughly constant in k) on top of
# the O(k^2*d) dist() cost. One full serial sweep over k=2..500 (one
# replicate at every size) measured at 88.4 s. At niter=10000 and cores=20:
#   serial total    ~ 245 h  (88.4 s * 10000 / 3600)
#   parallel wall   ~ 12.3 h (88.4 s * 10000 / 20 / 3600)
# at cores=22: ~11.2 h wall. This is markedly slower than a rpoisball.unit.
# fast()-based rewrite would be (that sampler has no per-call eigen cost;
# an earlier throwaway benchmark of the fast sampler projected ~4 h wall at
# cores=20 for the same band) -- but the task brief for this file explicitly
# specifies reusing rpoisball.unit and 01g's pattern, so the slower, more
# conservative path is what is implemented here. If ~12 h is judged too
# long when this is actually queued, swapping rep_one_size()'s call to
# rpoisball.unit.fast() (01e_nn_fast.R) would cut it roughly 2.5-3x with no
# distributional change (rnorm-equivalence proof in 01e's header) --
# flagged here as an option, not applied.
#
# THIS SCRIPT IS WRITTEN BUT NOT RUN (per task instructions). It is invoked
# by the orchestrator's queue script (queue_wp4_rerun.ps1) as the first
# production step, before any 04_wp4_runtime.R grid cells run.
#
# MODES (PowerShell; Rscript under Bash segfaults):
#   Rscript "revision_experiments/01h_nn_d500_quant.R" smoke
#   Rscript "revision_experiments/01h_nn_d500_quant.R" gen [niter] [cores]
#
# gen writes R/NN-test_quantile/NN-test-simul_500d_999%.RData (object
# `simul` = list(average, median), length 500 each -- the exact shape/name
# harness.R's get_simul("NN", 500) and nn_quant_for_d(500) == "999" expect).

suppressPackageStartupMessages({ library(parallel); library(here) })

source(here::here("R/ccds/NN_Dist_Est.R"))   # canonical (ORIGINAL, non-fast) rpoisball.unit

D          <- 500L
NMAX       <- 500L
QUANT      <- 0.999
QUANT_STR  <- "999"
TABLE_DIR  <- here::here("R/NN-test_quantile")
OUT_FILE   <- file.path(TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", D, QUANT_STR))
SEED       <- 20260719L

# One Monte Carlo rep for one size: simulate k points via the ORIGINAL
# rpoisball.unit, return c(mean, median) of the per-point NN distances.
# Identical statistic to NNDest(P).simpois and to 01g_nn_sizes_quant.R's
# rep_one_size() (copied verbatim; d is fixed at 500 here so the signature
# drops the d argument that 01g carries generically).
rep_one_size <- function(k, d = D) {
  x <- rpoisball.unit(k, d)
  dd <- as.matrix(dist(x))
  diag(dd) <- Inf
  nn <- apply(dd, 1, min)
  c(ave = mean(nn), med = median(nn))
}

# Full-band null quantiles at EVERY size 2..nmax (no knots, no
# interpolation -- see header). Parallel over reps within each size;
# reproducible via clusterSetRNGStream(seed). Structurally identical to
# 01g_nn_sizes_quant.R's NNDestP.sizes.lower.quant(), just renamed and with
# `sizes` always the full 2:nmax band.
NNDestP.fullband.lower.quant <- function(nmax, d = D, quant = QUANT, niter = 10000,
                                         cores = 20, seed = SEED) {
  sizes <- 2:nmax
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, c("rep_one_size", "rpoisball.unit"), envir = environment())
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS))))

  out <- list(average = numeric(length(sizes)), median = numeric(length(sizes)))
  t_start <- Sys.time()
  for (idx in seq_along(sizes)) {
    k <- sizes[idx]
    m <- parSapply(cl, seq_len(niter), function(i, k, d) rep_one_size(k, d), k = k, d = d)
    out$average[idx] <- unname(quantile(m["ave", ], 1 - quant))
    out$median[idx]  <- unname(quantile(m["med", ], 1 - quant))
    if (idx %% 25 == 0 || idx == length(sizes)) {
      cat(sprintf("  size %5d/%d done (%.1f min elapsed)\n", k, nmax,
                  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
    }
  }
  out
}

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[[1]] else "smoke"

if (MODE == "smoke") {
  # Cheap sanity check ONLY (d=5, two tiny sizes, niter=40, 4 cores) --
  # proves the machinery runs and the k=1 convention is sound. Does NOT
  # touch D=500 or write any table file.
  cat("=== 01h smoke: d=5, sizes {50,100}, niter=40, cores=4 ===\n")
  r <- NNDestP.fullband.lower.quant(nmax = 100, d = 5, niter = 40, cores = 4, seed = 1L)
  # nmax=100 here means sizes 2:100; just check shape/finiteness, not d=500 specifics.
  stopifnot(length(r$average) == 99, length(r$median) == 99,
            all(is.finite(r$average)), all(is.finite(r$median)),
            all(r$average > 0), all(r$median > 0))
  cat("SMOKE PASS\n")

} else if (MODE == "gen") {
  niter <- if (length(args) >= 2) as.integer(args[[2]]) else 10000L
  cores <- if (length(args) >= 3) as.integer(args[[3]]) else 20L
  cat(sprintf("=== 01h gen: d=%d nmax=%d niter=%d cores=%d quant=%.3f seed=%d ===\nStart: %s\n",
              D, NMAX, niter, cores, QUANT, SEED, format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  band <- NNDestP.fullband.lower.quant(nmax = NMAX, d = D, quant = QUANT, niter = niter,
                                       cores = cores, seed = SEED)

  # Prepend the k=1 convention (no defined NN distance for a single point;
  # matches 01e_nn_fast.R's `c(0, NN.dist.temp[1, ])` and is never actually
  # indexed by nnccd.radi's ascend/descend loops, which start at
  # low.num - 1 = 2 at the earliest -- see harness.R/04_wp4_runtime.R).
  simul <- list(
    average = c(0, band$average),
    median  = c(0, band$median)
  )

  # ---- SAVE FIRST -----------------------------------------------------------
  # 2026-07-22: the first production run (niter=5000, 22 cores, 339 min) was
  # LOST because the monotonicity stopifnot below fired before save() and R
  # halted. Never discard hours of Monte Carlo on a verification failure:
  # persist the artifact, then verify and report loudly. A suspect table on
  # disk can be inspected or deleted; a halted run cannot be recovered.
  save(simul, file = OUT_FILE)
  cat(sprintf("Saved: %s\n", OUT_FILE))

  # ---- Verification block ---------------------------------------------------
  # Structural checks stay hard (these indicate a real bug, and the file is
  # already on disk so nothing is lost by stopping).
  stopifnot(
    "length must be exactly NMAX"        = length(simul$average) == NMAX,
    "length must be exactly NMAX (med)"  = length(simul$median) == NMAX,
    "no NAs in average"                  = !anyNA(simul$average),
    "no NAs in median"                   = !anyNA(simul$median),
    "entry 1 (k=1) must be 0 by convention" = simul$average[1] == 0 && simul$median[1] == 0,
    "entries 2.. must be positive"       = all(simul$average[-1] > 0) && all(simul$median[-1] > 0)
  )

  # Trend check, REGIME-AWARE. Mean NN distance scales as k^(-1/d), so the
  # expected decline across the whole band is 1 - (NMAX/2)^(-1/d): ~29% at
  # d=10 but only ~1.1% at d=500. Where the true trend is that flat, MC noise
  # dominates and per-step monotonicity is a coin flip -- the original
  # ">70% of steps non-increasing" test is then guaranteed to fail on a
  # perfectly good table (the 2026-07-22 run measured 0.540). So require
  # per-step monotonicity only when the trend is big enough to see, and
  # otherwise check the SMOOTHED endpoints plus overall flatness.
  expected_decline   <- 1 - (NMAX / 2)^(-1 / D)
  frac_nonincreasing <- mean(diff(simul$average[-1]) <= 0)
  cat(sprintf("[verify] expected band decline under k^(-1/d): %.2f%% | frac of steps non-increasing: %.3f\n",
              100 * expected_decline, frac_nonincreasing))
  if (expected_decline > 0.05) {
    stopifnot("average envelope should be broadly decreasing in k (allow MC noise)" =
                frac_nonincreasing > 0.7)
  } else {
    w    <- max(10L, NMAX %/% 20L)
    head_m <- mean(head(simul$average[-1], w)); tail_m <- mean(tail(simul$average[-1], w))
    cat(sprintf("[verify] near-flat regime: mean(first %d)=%.6f mean(last %d)=%.6f ratio=%.4f\n",
                w, head_m, w, tail_m, tail_m / head_m))
    stopifnot("smoothed envelope must not INCREASE materially" = tail_m <= head_m * 1.02)
    stopifnot("envelope must be near-flat at high d (|ratio-1| < 10%)" = abs(tail_m / head_m - 1) < 0.10)
  }
  cat(sprintf("[verify] head(average, 5) = %s\n", paste(round(head(simul$average, 5), 5), collapse = ", ")))
  cat(sprintf("[verify] tail(average, 5) = %s\n", paste(round(tail(simul$average, 5), 5), collapse = ", ")))
  cat("[verify] PASS\nDONE gen (elapsed shown per-size above)\n")

} else {
  stop("Unknown mode. Use smoke | gen [niter] [cores]", call. = FALSE)
}
