# 01e_nn_fast.R  (T3 Phase B-prep)
# Optimized override of NNDestP.simpois.lower.quant() (R/ccds/NN_Dist_Est.R,
# untouched) for production-scale NN quantile-table generation at high d
# (166, 274, 400). Sourced by 01_gen_quantile_table.R for every NN run.
#
# Two independent optimizations, both validated by 01f_validate_nn_fast.R:
#
# --------------------------------------------------------------------------
# (1) rnorm-equivalence: rpoisball.unit.fast() replaces
#         norm.data <- matrix(mvrnorm(n, rep(0,d), diag(d)), ncol=d, byrow=T)
#     with
#         norm.data <- matrix(rnorm(n * d), nrow = n, ncol = d)
#
#     MASS::mvrnorm(n, mu, Sigma) internally draws X <- matrix(rnorm(n*p), n, p)
#     and rotates it by the eigendecomposition of Sigma:
#         eS <- eigen(Sigma, symmetric = TRUE)
#         X <- X %*% t(eS$vectors %*% diag(sqrt(pmax(eS$values, 0)), p))  [+ mu]
#     For Sigma = diag(d): every eigenvalue is 1, so the scaling is a no-op,
#     and the rotation matrix eS$vectors is SOME fixed orthogonal matrix Q
#     (LAPACK's arbitrary-but-deterministic choice for a degenerate spectrum
#     -- empirically NOT the identity; verified in 01f). Rotating n iid
#     N(0,1) vectors by any fixed orthogonal Q preserves N(0, I_d):
#     if Z ~ N(0,I) then QZ ~ N(0, Q Q^T) = N(0, I). So the swap is
#     DISTRIBUTIONALLY IDENTICAL to the original call, not an approximation
#     -- it just skips the O(d^3) eigen(diag(d)) LAPACK call that mvrnorm()
#     pays on every single invocation (once per x = 1..n inside SimuOnce,
#     i.e. n eigendecompositions of an already-diagonal d x d matrix per
#     Monte-Carlo iteration).
#
#     NOTE: the RNG call sequence differs from the original (rnorm() is not
#     rotated, mvrnorm()'s draws are) so outputs are NOT bitwise reproducible
#     from the original under the same seed -- validated statistically
#     instead in 01f_validate_nn_fast.R.
#
# (2) Streaming SimuOnce: the original SimuOnce() builds
#         data.simu.list <- lapply(1:n, rpoisball.unit, d = d)
#     BEFORE any distance computation, holding n matrices of sizes
#     1*d, 2*d, ..., n*d simultaneously resident: sum_{x=1}^n x*d elements =
#     d*n*(n+1)/2. At n=1000, d=400 that's 400*500500*8 bytes ~= 1.6 GB
#     PER WORKER just for this list (see revision_experiments/results/probes
#     /probe_report.csv: NN d=400 was capped at 16 of 20 cores by RAM).
#     data.simu.list[[x]] is used exactly once, immediately, inside the
#     x-th term of the sapply(2:n, ...) that follows (list element x=1 is
#     never used for a distance computation at all -- only the "NN distance
#     of a single point is 0" placeholder is used, reproduced below via the
#     same hand-added leading `c(0, ...)` as the original). Generating
#     data.temp on demand inside the sapply loop -- instead of pre-building
#     the whole list -- keeps at most ONE matrix resident at a time
#     (<= n*d elements, e.g. 1000*400*8 bytes ~= 3.2 MB), cutting peak
#     per-worker RAM by roughly 2-3 orders of magnitude. Output is
#     unaffected: same draws (in expectation, same distribution), same
#     summary statistics, same object shape.
#
# Output contract preserved exactly: NNDestP.simpois.lower.quant.fast()
# returns list(average = <length n>, median = <length n>), identical to
# the original NNDestP.simpois.lower.quant().

# Top-level copy: used by serial validation code (01f) and any interactive
# caller. NOT visible to %dopar% workers -- foreach serializes SimuOnce's
# closure but not the master's .GlobalEnv, so NNDestP.simpois.lower.quant.fast
# below defines its own copy INSIDE SimuOnce, exactly as the original
# NNDestP.simpois.lower.quant does with rpoisball.unit. Keep the two
# definitions textually identical.
rpoisball.unit.fast <- function(n, d) {
  r1 <- runif(n, 0, 1)^(1 / d)
  norm.data <- matrix(rnorm(n * d), nrow = n, ncol = d)
  data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
  data1 <- apply(data1, 1, function(x) x * r1)
  return(data1)
}

NNDestP.simpois.lower.quant.fast <- function(n, d, quant = 0.99, niter = 100, shape = "sphere") {
  NN.dist.ave.mat <- NULL
  NN.dist.med.mat <- NULL

  if (shape == "sphere") {
    SimuOnce <- function() {
      # Defined INSIDE SimuOnce (mirroring the original) so %dopar% workers
      # -- which receive SimuOnce's closure but not the master's .GlobalEnv --
      # can resolve it. Keep textually identical to the top-level copy above.
      rpoisball.unit.fast <- function(n, d) {
        r1 <- runif(n, 0, 1)^(1 / d)
        norm.data <- matrix(rnorm(n * d), nrow = n, ncol = d)
        data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
        data1 <- apply(data1, 1, function(x) x * r1)
        return(data1)
      }

      # Streaming: data.temp is generated fresh for each x inside the
      # sapply, never pre-materialized as a full list (see (2) above).
      NN.dist.temp <- sapply(2:n, function(x) {
        data.temp <- rpoisball.unit.fast(x, d)
        data.dist <- as.matrix(dist(data.temp))
        diag(data.dist) <- Inf
        NN.dist.ttemp <- apply(data.dist, 1, min) # the Nearest Neighbor distance for each point
        NN.dist.ttemp.ave <- mean(NN.dist.ttemp)
        NN.dist.ttemp.med <- median(NN.dist.ttemp)
        return(c(NN.dist.ttemp.ave, NN.dist.ttemp.med))
      })
      NN.dist.temp.ave <- c(0, NN.dist.temp[1, ])
      NN.dist.temp.med <- c(0, NN.dist.temp[2, ])
      return(c(NN.dist.temp.ave = list(NN.dist.temp.ave), NN.dist.temp.med = list(NN.dist.temp.med)))
    }

    cores <- detectCores()
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    NNest.m <- foreach(1:niter, .packages = c("MASS", "cluster", "igraph")) %dopar% SimuOnce()
    stopCluster(cl)

    for (i in 1:niter) {
      NN.dist.ave.mat <- rbind(NN.dist.ave.mat, NNest.m[[i]]$NN.dist.temp.ave)
      NN.dist.med.mat <- rbind(NN.dist.med.mat, NNest.m[[i]]$NN.dist.temp.med)
    }
  }

  quant.ave.lower <- sapply(1:n, function(x) {
    quant <- quantile(NN.dist.ave.mat[, x], 1 - quant)
  })
  names(quant.ave.lower) <- NULL

  quant.med.lower <- sapply(1:n, function(x) {
    quant <- quantile(NN.dist.med.mat[, x], 1 - quant)
  })
  names(quant.med.lower) <- NULL

  return(list(average = quant.ave.lower, median = quant.med.lower))
}
