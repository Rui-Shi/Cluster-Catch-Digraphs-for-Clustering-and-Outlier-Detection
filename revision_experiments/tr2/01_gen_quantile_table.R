#!/usr/bin/env Rscript
# 01_gen_quantile_table.R
# Task T3 Phase A: quantile-table generator for the PR-D-26-05767 revision
# experiments pipeline.
#
# Generates a per-dimension Monte-Carlo quantile lookup table for CCD
# construction, reusing the ORIGINAL generator functions
# (KestP.simpois.edge.quantile in R/ccds/Kest.R, and
#  NNDestP.simpois.lower.quant in R/ccds/NN_Dist_Est.R) unmodified.
# Exception 1: for RK at d >= 342 the original edge-correction weight overflows
# (gamma() -> Inf/NaN); a log-space-stable override defined in THIS file is
# selected automatically there (see "Numerically stable RK override" below;
# validated by 01c_validate_rk_stable.R). Original files are never touched.
# Exception 2: for NN (every d), a performance override defined in
# revision_experiments/tr2/01e_nn_fast.R is used in place of the original
# NNDestP.simpois.lower.quant -- distributionally identical output (proven,
# not approximated; see that file's header), same statistical definition and
# object shape, but avoids an O(d^3) eigendecomposition per draw and cuts
# peak per-worker RAM by ~2-3 orders of magnitude by streaming the
# per-iteration workload instead of materializing it as a list first.
# Validated by 01f_validate_nn_fast.R (statistical equivalence + timing +
# RAM). Needed to make d=400 at niter=10000/cores=20 tractable at all
# (original path RAM-caps at ~16 cores there); originals are never touched.
#
# --------------------------------------------------------------------------
# Reverse-engineered production parameters (see revision_experiments notes /
# final report for full derivation). Summary:
#
# Both KestP.simpois.edge.quantile() and NNDestP.simpois.lower.quant() are
# shape-invariant in d: the SAVED object's component shapes depend only on
# the per-iteration sample size (m for RK, n for NN) and, for RK, the radius
# count rn (r <- seq(1/rn, 1, 1/rn)), NEVER on d itself. d only changes the
# cost per iteration (dist()/mvrnorm() cost scales with d; the RK edge-
# correction integrate() calls may also get costlier/less stable at high d).
#
# Inspecting the recovered production .RData files directly
# (R/RK-test_quantile/*.RData, R/NN-test_quantile/*.RData) shows TWO tiers,
# not one fixed parameter set:
#
#   Tier "low/mid d" (observed at RK d=10,20,35, q=0.999; NN d=10,20, q=0.999):
#     RK: Kest.m is 2000 x 50000  -> niter=2000, m*rn=50000, rn=10 -> m=5000
#     NN: average/median length 5000                              -> n=5000
#     (This tier's shapes do NOT match the currently committed low-d driver
#      scripts, e.g. R/RK-test_quantile/10d_999%.R has n=1000, niter=10000 --
#      those committed drivers are stale/draft versions; the shapes on disk
#      are ground truth per the task brief.)
#
#   Tier "high d" (observed at RK d=50,100, q=0.999; NN d=50,100, q=0.999):
#     RK: Kest.m is 10000 x 10000 -> niter=10000, m*rn=10000, rn=10 -> m=1000
#     NN: average/median length 1000                               -> n=1000
#     This tier EXACTLY matches the one still-committed high-d driver,
#     R/NN-test_quantile/50100d_999%.R (d=50 then d=100, iteN=10000, n=1000,
#     quant=0.999), giving strong, direct (not just shape-inferred)
#     confirmation for n=1000/niter=10000 at d>=50. No committed RK driver
#     exists for d=50/100, so m=1000/niter=10000 there is shape-inferred, but
#     it aligns with NN's confirmed value for the same d's and quantile,
#     which is the best evidence available.
#
# The revision's new target dimensions (166, 274, 400, 500, 1000, 1555) are
# all well above the d=35/d=50 tier boundary, and we have no evidence of a
# third, even-more-reduced tier beyond d=100. We therefore adopt the "high d"
# tier verbatim as the production parameter set for all new tables:
#   RK: m = 1000, rn = 10, niter = 10000, quan = 0.999
#   NN: n = 1000,          niter = 10000, quant = 0.999, shape = "sphere"
# rn is fixed at 10 in every recovered file (r always seq(0.1, 1, by=0.1)).
# --------------------------------------------------------------------------

PROD_M   <- 1000L   # RK per-iteration sample size (== NN's n)
PROD_RN  <- 10L      # RK radius-bin count; r <- seq(1/rn, 1, 1/rn)
PROD_N   <- 1000L    # NN per-iteration sample size

suppressPackageStartupMessages({
  library(here)
  library(parallel)
  library(doParallel)
  library(MASS)
  library(igraph)
  library(cluster)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript \"revision_experiments/tr2/01_gen_quantile_table.R\" ",
    "<variant RK|NN> <d> <quant e.g. 0.999> <niter> <cores> [outdir]",
    call. = FALSE
  )
}

variant <- toupper(args[[1]])
d       <- as.integer(args[[2]])
quant   <- as.numeric(args[[3]])
niter   <- as.integer(args[[4]])
cores   <- as.integer(args[[5]])
outdir  <- if (length(args) >= 6) args[[6]] else NA_character_

if (!variant %in% c("RK", "NN")) {
  stop("variant must be RK or NN, got: ", variant, call. = FALSE)
}
stopifnot(
  "d must be a positive integer"      = is.finite(d) && d > 0,
  "quant must be in (0,1)"            = is.finite(quant) && quant > 0 && quant < 1,
  "niter must be a positive integer"  = is.finite(niter) && niter > 0,
  "cores must be a positive integer"  = is.finite(cores) && cores > 0
)

# Quantile filename token: strip the leading "0." from the decimal (matches
# the convention already used on disk, e.g. 0.999 -> "999", 0.99 -> "99").
quant_str <- sub("^0\\.", "", format(quant, trim = TRUE))

cat(sprintf(
  "=== 01_gen_quantile_table.R ===\nvariant=%s d=%d quant=%s (token '%s') niter=%d cores=%d outdir=%s\n",
  variant, d, quant, quant_str, niter, cores,
  if (is.na(outdir)) "<production default>" else outdir
))

repo_root <- here::here()

# ---- Resolve output path -------------------------------------------------
if (is.na(outdir)) {
  out_subdir <- if (variant == "RK") "R/RK-test_quantile" else "R/NN-test_quantile"
  out_dir_path <- file.path(repo_root, out_subdir)
} else {
  out_dir_path <- outdir
}
if (!dir.exists(out_dir_path)) {
  dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)
}
out_file <- file.path(
  out_dir_path,
  sprintf("%s-test-simul_%dd_%s%%.RData", variant, d, quant_str)
)
cat("Output file:", out_file, "\n")

# ---- Source the ORIGINAL generator functions (unmodified) ----------------
if (variant == "RK") {
  source(here::here("R/ccds/Kest.R"))
} else {
  source(here::here("R/ccds/NN_Dist_Est.R"))
  # Optimized NN override (rnorm-equivalence + streaming) -- see file header
  # and "Exception 2" note above. Selected for every NN run, all d.
  source(here::here("revision_experiments/tr2/01e_nn_fast.R"))
}

# ---- Numerically stable RK override for very high d -----------------------
# The ORIGINAL RK edge-correction weight in KestP.simpois.edge.quantile is
#     cons  <- (sqrt(pi)*gamma((d+1)/2)) / (2*gamma(d/2+1))
#     ftemp <- cons * (1 / integrate(sin(t)^d, 0, acos(dist/2))$value)
# gamma() overflows to Inf at arguments > ~171.6, so:
#     d <= 341 : cons finite (correct)
#     d  = 342 : gamma(d/2+1) = Inf, cons silently 0  -> all-zero weights
#     d >= 343 : cons = Inf/Inf = NaN                 -> all-NaN Kest.m and
#                quantile() aborts ("missing values and NaN's not allowed")
# Additionally sin(t)^d underflows to exactly 0 for far pairs at very high d,
# giving 1/0 = Inf weights and 0*Inf = NaN under the original masking.
#
# The override below computes the SAME quantity in log space:
#     log(cons) = log(sqrt(pi)) + lgamma((d+1)/2) - log(2) - lgamma(d/2+1)
#     log(I)    = d*log(sin(b)) + log( integrate(exp(d*(log sin u - log sin b)), 0, b) )
#                 with b = acos(dist/2)   [exact factorization, no approximation]
# and masks BEFORE multiplying, so Inf weights on masked (far) pairs never
# meet a 0. Validated against the original at d where both work by
# revision_experiments/tr2/01c_validate_rk_stable.R. Selected automatically for
# d >= RK_STABLE_MIN_D; the untouched original code path is used below that,
# so regeneration at d = 166, 274 keeps the original behavior bit-for-bit.
RK_STABLE_MIN_D <- 342L

KestP.simpois.edge.quantile.stable <- function(m, d, rn, quan, niter) {
  r <- seq(1 / rn, 1, 1 / rn)
  SimuOnce <- function() {
    rpoisball.unit <- function(n, d) {
      r1 <- runif(n, 0, 1)^(1 / d)
      norm.data <- matrix(mvrnorm(n, rep(0, d), diag(d)), ncol = d, byrow = T)
      data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
      data1 <- apply(data1, 1, function(x) x * r1)
      return(data1)
    }
    temp <- rpoisball.unit(m, d)
    temp.dist <- as.matrix(dist(temp))

    # log-space edge-correction weights (see block comment above)
    log_cons <- log(sqrt(pi)) + lgamma((d + 1) / 2) - log(2) - lgamma(d / 2 + 1)
    log_ftemp_one <- function(t) {
      b <- acos(min(max(t / 2, -1), 1))
      if (b <= 0) return(Inf)          # zero-length integral -> Inf weight (always masked)
      lsb <- log(sin(b))
      I_scaled <- integrate(function(u) exp(d * (log(sin(u)) - lsb)), 0, b)$value
      logI <- d * lsb + log(I_scaled)
      return(log_cons - logI)
    }
    ftemp <- sapply(temp.dist, log_ftemp_one, simplify = TRUE)
    ftemp <- matrix(exp(ftemp), nrow = nrow(temp.dist), byrow = FALSE)

    diag(temp.dist) <- Inf
    result <- sapply(r, function(x) {
      Mtemp <- (temp.dist < x)
      Mtemp[lower.tri(Mtemp)] <- 0
      W <- ftemp
      W[Mtemp == 0] <- 0               # mask FIRST: Inf weights on far pairs never meet a 0
      sumM <- cumsum(2 * colSums(W))
      return(sumM / ((1:m) * (1:m)))
    }, simplify = TRUE)
    return(as.vector(result))
  }
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  Kest.m <- foreach(1:niter, .packages = c("MASS", "cluster", "igraph")) %dopar% SimuOnce()
  stopCluster(cl)

  Kest.m <- do.call(rbind, Kest.m)

  Kest.quan <- list()
  for (cur_quan in quan) {
    temp <- apply(Kest.m, 2, quantile, probs = as.numeric(cur_quan))
    Kest.quan[[as.character(cur_quan)]] <- matrix(temp, nrow = m)
  }

  return(list(Kest.m = Kest.m, quan = Kest.quan, r = r))
}

# ---- Honor the `cores` CLI argument -------------------------------------
# KestP.simpois.edge.quantile() and NNDestP.simpois.lower.quant() both call
# a hard-coded `cores <- detectCores()` internally (no cores parameter
# exists in the original signatures) and then `makeCluster(cores)`. Since
# these functions are sourced at top level, their closure environment is
# .GlobalEnv; R resolves the free symbol `detectCores` by looking in the
# function's enclosing environment (.GlobalEnv) BEFORE walking the search
# path where the real parallel::detectCores lives. Shadowing it here lets us
# honor the requested core count without touching R/ccds/Kest.R or
# R/ccds/NN_Dist_Est.R.
requested_cores <- cores
detectCores <- function(...) requested_cores  # nolint: intentional shadow

# ---- Run ------------------------------------------------------------------
t_start <- Sys.time()
cat("Start:", format(t_start), "\n")

if (variant == "RK") {
  if (d >= RK_STABLE_MIN_D) {
    cat(sprintf(
      "NOTE: d=%d >= %d -- original RK weights overflow (gamma() -> Inf/NaN).\n      Using log-space stable override KestP.simpois.edge.quantile.stable\n      (validated vs original at lower d by 01c_validate_rk_stable.R).\n",
      d, RK_STABLE_MIN_D
    ))
    cat(sprintf(
      "Calling KestP.simpois.edge.quantile.stable(m=%d, d=%d, rn=%d, quan=%s, niter=%d)\n",
      PROD_M, d, PROD_RN, quant, niter
    ))
    simul <- KestP.simpois.edge.quantile.stable(PROD_M, d, PROD_RN, quant, niter)
  } else {
    cat(sprintf(
      "Calling KestP.simpois.edge.quantile(m=%d, d=%d, rn=%d, quan=%s, niter=%d)\n",
      PROD_M, d, PROD_RN, quant, niter
    ))
    simul <- KestP.simpois.edge.quantile(PROD_M, d, PROD_RN, quant, niter)
  }
} else {
  cat(sprintf(
    "Calling NNDestP.simpois.lower.quant.fast(n=%d, d=%d, quant=%s, niter=%d, shape='sphere') [optimized override, see 01e_nn_fast.R]\n",
    PROD_N, d, quant, niter
  ))
  simul <- NNDestP.simpois.lower.quant.fast(PROD_N, d, quant = quant, niter = niter, shape = "sphere")
}

t_end <- Sys.time()
elapsed_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
cat("End:", format(t_end), "\n")
cat(sprintf("Elapsed seconds: %.3f\n", elapsed_sec))

# ---- Report shape actually produced (sanity echo) -------------------------
if (variant == "RK") {
  quant_key <- as.character(quant)
  cat(sprintf(
    "Result shapes: Kest.m %s | quan[['%s']] %s | r length %d\n",
    paste(dim(simul$Kest.m), collapse = "x"),
    quant_key,
    paste(dim(simul$quan[[quant_key]]), collapse = "x"),
    length(simul$r)
  ))
} else {
  cat(sprintf(
    "Result shapes: average length %d | median length %d\n",
    length(simul$average), length(simul$median)
  ))
}

save(simul, file = out_file)
cat("Saved:", out_file, "\n")
cat(sprintf("TOTAL_ELAPSED_SECONDS=%.3f\n", elapsed_sec))
cat("=== DONE ===\n")
