# 01g_nn_sizes_quant.R -- size-targeted NN null-quantile generation
#
# WHY (FINDINGS "Design upgrade 2026-07-17"): table entry k is the null
# (1-quant) quantile of the mean/median NN distance of a k-point Poisson
# process in the unit d-ball -- dataset-independent. Existing production
# tables carry entries 1..1000; full-data Musk (n=3062, d=166) and Speech
# (n=3686, d=400) need the band 1001..n. Generating every size costs
# Sum k^2*d ~ n^3*d/3 per iteration; generating a handful of KNOT sizes and
# interpolating costs ~k^2*d per size -- ~1000x less. The envelope is nearly
# flat over the band (mean NN dist ~ k^(-1/d): 0.7% total variation over
# 1000->3062 at d=166; 0.3% over 1000->3686 at d=400), so linear
# interpolation between knots is essentially exact. The interpolated band is
# offset-anchored at the existing table's entry 1000 so there is no seam at
# the splice boundary; the raw knot-vs-table discrepancy at k=1000 is
# reported as a Monte Carlo cross-validation of the whole procedure.
#
# The sampler is the ORIGINAL rpoisball.unit sourced from R/ccds/
# NN_Dist_Est.R (no reimplementation -- zero distributional-drift risk).
#
# MODES (PowerShell; Rscript under Bash segfaults):
#   Rscript "revision_experiments/tr2/01g_nn_sizes_quant.R" smoke
#   Rscript "revision_experiments/tr2/01g_nn_sizes_quant.R" probe <musk|speech> [cores]
#   Rscript "revision_experiments/tr2/01g_nn_sizes_quant.R" gen   <musk|speech> <niter> [cores]
#
# gen writes:
#   results/nn_sizes_knots_<dataset>.RData        (raw knot quantiles + meta)
#   R/NN-test_quantile/NN-test-simul_<d>d_999%_n<nmax>_spliced.RData
#     (object `simul` = list(average, median), length nmax -- same shape the
#      construction code consumes; entries 1..1000 are byte-identical to the
#      production table, entries 1001..nmax are anchored interpolation)

suppressPackageStartupMessages({ library(parallel); library(here) })

source(here::here("R/ccds/NN_Dist_Est.R"))   # canonical rpoisball.unit (+MASS)

SPECS <- list(
  musk   = list(d = 166, nmax = 3062, knots = c(1000, 1400, 1900, 2400, 2800, 3062)),
  speech = list(d = 400, nmax = 3686, knots = c(1000, 1500, 2100, 2700, 3300, 3686))
)
QUANT      <- 0.999
TABLE_DIR  <- here::here("R/NN-test_quantile")
RESULTS    <- here::here("revision_experiments/results/tr2")

# One Monte Carlo rep for one size: simulate k points, return c(mean, median)
# of the per-point NN distances. Identical statistic to NNDest(P).simpois.
rep_one_size <- function(k, d) {
  x <- rpoisball.unit(k, d)
  dd <- as.matrix(dist(x))
  diag(dd) <- Inf
  nn <- apply(dd, 1, min)
  c(ave = mean(nn), med = median(nn))
}

# Null quantiles at the requested sizes only. Parallel over reps; reproducible
# via clusterSetRNGStream(seed).
NNDestP.sizes.lower.quant <- function(sizes, d, quant = QUANT, niter = 1000,
                                      cores = 20, seed = 20260717L) {
  cl <- makeCluster(cores)
  on.exit(stopCluster(cl), add = TRUE)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, c("rep_one_size", "rpoisball.unit"), envir = environment())
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS))))

  out <- list(average = numeric(0), median = numeric(0))
  for (k in sizes) {
    t0 <- Sys.time()
    m <- parSapply(cl, seq_len(niter), function(i, k, d) rep_one_size(k, d), k = k, d = d)
    q_ave <- unname(quantile(m["ave", ], 1 - quant))
    q_med <- unname(quantile(m["med", ], 1 - quant))
    out$average[as.character(k)] <- q_ave
    out$median[as.character(k)]  <- q_med
    cat(sprintf("  size %5d: q_ave=%.6f q_med=%.6f  (%.1f s, niter=%d)\n",
                k, q_ave, q_med, as.numeric(difftime(Sys.time(), t0, units = "secs")), niter))
  }
  out
}

# Anchored linear interpolation of one envelope component over 1001..nmax,
# spliced onto the production entries 1..1000.
splice_component <- function(existing_1k, knots, knot_vals, nmax) {
  offset <- existing_1k[1000] - knot_vals[knots == 1000]
  band_x <- 1001:nmax
  band   <- approx(knots, knot_vals + offset, xout = band_x, rule = 2)$y
  c(existing_1k[1:1000], band)
}

load_production <- function(d) {
  f <- file.path(TABLE_DIR, sprintf("NN-test-simul_%dd_999%%.RData", d))
  stopifnot(file.exists(f))
  e <- new.env(); load(f, envir = e)
  e$simul
}

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[[1]] else "smoke"

if (MODE == "smoke") {
  cat("=== 01g smoke: d=5, sizes {50,100}, niter=40, cores=4 ===\n")
  r <- NNDestP.sizes.lower.quant(c(50, 100), d = 5, niter = 40, cores = 4, seed = 1L)
  stopifnot(length(r$average) == 2, length(r$median) == 2,
            all(is.finite(r$average)), all(is.finite(r$median)),
            all(r$average > 0), all(r$median > 0),
            r$average["100"] < r$average["50"])   # NN dist shrinks with k
  cat("SMOKE PASS\n")

} else if (MODE == "probe") {
  ds <- SPECS[[args[[2]]]]
  cores <- if (length(args) >= 3) as.integer(args[[3]]) else 20L
  kmax <- max(ds$knots); np <- 2L * cores    # two scheduling waves
  cat(sprintf("=== 01g probe: %s (d=%d) size=%d niter=%d cores=%d ===\n",
              args[[2]], ds$d, kmax, np, cores))
  t0 <- Sys.time()
  invisible(NNDestP.sizes.lower.quant(kmax, d = ds$d, niter = np, cores = cores))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  per_rep_wall <- wall / np                       # amortized wall per rep at this core count
  cat(sprintf("probe wall %.1f s -> per-rep (amortized) %.2f s\n", wall, per_rep_wall))
  for (nit in c(1000, 2000, 10000)) {
    est_knots <- per_rep_wall * nit * sum((ds$knots / kmax)^2)  # k^2 cost scaling across knots
    cat(sprintf("  projected gen at niter=%5d: ~%.1f min for all %d knots\n",
                nit, est_knots / 60, length(ds$knots)))
  }

} else if (MODE == "gen") {
  ds <- SPECS[[args[[2]]]]
  niter <- as.integer(args[[3]])
  cores <- if (length(args) >= 4) as.integer(args[[4]]) else 20L
  cat(sprintf("=== 01g gen: %s d=%d nmax=%d knots={%s} niter=%d cores=%d quant=%.3f ===\nStart: %s\n",
              args[[2]], ds$d, ds$nmax, paste(ds$knots, collapse = ","), niter, cores, QUANT,
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  knots_res <- NNDestP.sizes.lower.quant(ds$knots, d = ds$d, niter = niter, cores = cores)

  prod <- load_production(ds$d)
  rel <- function(a, b) abs(a - b) / abs(b)
  v_ave <- rel(knots_res$average["1000"], prod$average[1000])
  v_med <- rel(knots_res$median["1000"],  prod$median[1000])
  cat(sprintf("cross-check at k=1000 (fresh knot vs production table): ave rel.diff %.3f%%, med rel.diff %.3f%%\n",
              100 * v_ave, 100 * v_med))

  simul <- list(
    average = splice_component(prod$average, ds$knots, knots_res$average, ds$nmax),
    median  = splice_component(prod$median,  ds$knots, knots_res$median,  ds$nmax)
  )
  stopifnot(length(simul$average) == ds$nmax, !anyNA(simul$average), !anyNA(simul$median))

  kf <- file.path(RESULTS, sprintf("nn_sizes_knots_%s.RData", args[[2]]))
  meta <- list(dataset = args[[2]], d = ds$d, nmax = ds$nmax, knots = ds$knots,
               niter = niter, quant = QUANT, seed = 20260717L,
               crosscheck_rel_ave = v_ave, crosscheck_rel_med = v_med,
               timestamp = Sys.time())
  save(knots_res, meta, file = kf)
  sf <- file.path(TABLE_DIR, sprintf("NN-test-simul_%dd_999%%_n%d_spliced.RData", ds$d, ds$nmax))
  save(simul, file = sf)
  cat(sprintf("Saved knots: %s\nSaved spliced table: %s\nDONE gen %s\n", kf, sf, args[[2]]))

} else {
  stop("Unknown mode. Use smoke | probe <musk|speech> [cores] | gen <musk|speech> <niter> [cores]")
}
