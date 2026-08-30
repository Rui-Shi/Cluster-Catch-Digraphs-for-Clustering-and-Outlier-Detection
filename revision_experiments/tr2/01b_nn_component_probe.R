#!/usr/bin/env Rscript
# 01b_nn_component_probe.R  (T3 Phase A support)
# Times the per-iteration cost COMPONENTS of NNDestP.simpois.lower.quant()
# (R/ccds/NN_Dist_Est.R) at a given dimension d, without running a full
# Monte-Carlo iteration. Used to size the NN probes and to extrapolate
# production cost where a single full iteration would exceed the probe
# budget (relevant at d=1555).
#
# One SimuOnce iteration in NNDestP.simpois.lower.quant does, for n=1000:
#   for x in 1..n:  rpoisball.unit(x, d)   [mvrnorm -> eigen(diag(d)) O(d^3)
#                                            per CALL, + x*d normal draws]
#   for x in 2..n:  dist(x points, d dims) [~x^2*d/2 flops] + row mins
# Total dist work = sum_{x=2}^{n} x^2*d/2 ~ (n^3/6)*d, i.e. ~334x the cost
# of the single largest (x=1000) dist call.
#
# Usage: Rscript "revision_experiments/01b_nn_component_probe.R" <d> [n=1000] [reps=3]

suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: 01b_nn_component_probe.R <d> [n] [reps]")
d    <- as.integer(args[[1]])
n    <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
reps <- if (length(args) >= 3) as.integer(args[[3]]) else 3L

cat(sprintf("=== NN component probe: d=%d n=%d reps=%d ===\n", d, n, reps))

tmed <- function(expr_fun, reps) {
  ts <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- Sys.time(); expr_fun(); ts[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  median(ts)
}

# (1) mvrnorm fixed overhead per call (eigen of diag(d)): draw only 2 points
t_eigen <- tmed(function() invisible(mvrnorm(2, rep(0, d), diag(d))), reps)
cat(sprintf("mvrnorm(2, d=%d) per-call [~eigen overhead]: %.4f s\n", d, t_eigen))

# (2) full rpoisball.unit at x = n (as defined inside NNDestP's SimuOnce)
rpoisball.unit <- function(n, d) {
  r1 <- runif(n, 0, 1)^(1 / d)
  norm.data <- matrix(mvrnorm(n, rep(0, d), diag(d)), ncol = d, byrow = TRUE)
  data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
  data1 <- apply(data1, 1, function(x) x * r1)
  return(data1)
}
t_gen_n <- tmed(function() invisible(rpoisball.unit(n, d)), reps)
cat(sprintf("rpoisball.unit(%d, d=%d): %.4f s\n", n, d, t_gen_n))

# (3) dist + row-min step at x = n (the largest single dist call)
dat <- rpoisball.unit(n, d)
t_dist_n <- tmed(function() {
  dm <- as.matrix(dist(dat)); diag(dm) <- Inf; invisible(apply(dm, 1, min))
}, reps)
cat(sprintf("dist+rowmin at x=%d, d=%d: %.4f s\n", n, d, t_dist_n))

# (4) model one full SimuOnce iteration:
#     eigen overhead: n calls; generation ~ sum_x x = n(n+1)/2 draws vs n at
#     x=n -> gen scales ~ (n+1)/2 relative to t_gen_n's non-eigen part;
#     dist work: sum_{x=2}^n x^2 relative to n^2.
gen_noneigen <- max(t_gen_n - t_eigen, 0)
sum_x  <- n * (n + 1) / 2
sum_x2 <- n * (n + 1) * (2 * n + 1) / 6
iter_est <- n * t_eigen + gen_noneigen * (sum_x / n) + t_dist_n * (sum_x2 / n^2)
cat(sprintf("Modeled ONE SimuOnce iteration at d=%d: %.1f s  (eigen part %.1f s, gen part %.1f s, dist part %.1f s)\n",
            d, iter_est, n * t_eigen, gen_noneigen * (sum_x / n), t_dist_n * (sum_x2 / n^2)))
cat(sprintf("MODELED_ITER_SECONDS=%.3f\n", iter_est))
