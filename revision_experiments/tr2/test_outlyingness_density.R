# test_outlyingness_density.R -- regression test for the Vic_Den high-d overflow
#
# Bug (found 2026-07-16, WP5 Arrhythmia d=274): Vic_Den computed
#   (size / R^d)^(1/d)
# which evaluates R^d FIRST. For d=274 that overflows double precision for
# any radius R > exp(709/274) ~= 13.3 (Inf -> Den = 0) and underflows for
# R < exp(-708/274) ~= 0.076 (0 -> Den = Inf). Standardized distances at
# d=274 concentrate near sqrt(2d) ~= 23, so covering radii routinely exceed
# 13.3. Downstream: OOS = 0/0 = NaN (crashes the tie-break loop with
# "missing value where TRUE/FALSE needed"), IOS = 1/0 = Inf (fails the
# score-sanity assert; 207/452 non-finite on Arrhythmia).
#
# Fix: algebraically identical reordering, size^(1/d) / R -- no intermediate
# R^d. This test asserts (a) finiteness at d=274 with realistic radii and
# (b) numerical agreement with the old formula at low d where the old
# formula is stable (guards the published low-d results).
#
# Run: Rscript "revision_experiments/tr2/test_outlyingness_density.R"

source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))

fail <- 0
check <- function(ok, label) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail <<- fail + 1
}

# --- (a) high-d finiteness: n=40, d=274, distances ~ sqrt(2d) ~ 23 --------
set.seed(42)
n <- 40; d <- 274
X <- scale(matrix(rnorm(n * d), n, d))
dd <- as.matrix(dist(X))
cat(sprintf("median pairwise distance at d=%d: %.2f (overflow threshold %.2f)\n",
            d, median(dd[upper.tri(dd)]), exp(709 / d)))
# radii chosen like covering radii: each point's 5th-NN distance (all > 13.3 here)
R <- apply(dd, 1, function(r) sort(r)[6])
stopifnot(all(R > exp(709 / d)))   # test precondition: radii in the overflow regime

den <- Vic_Den(X, R, d)
check(all(is.finite(den)) && all(den > 0), "Vic_Den finite and positive at d=274, R > 13.3")

oos <- OOS(X, R, d)
check(all(is.finite(oos)), "OOS finite at d=274 (no 0/0 = NaN)")

ios <- IOS(X, R, d)
check(all(is.finite(ios)), "IOS finite at d=274 (no 1/0 = Inf)")

# --- (b) low-d agreement with the original formula (stable regime) --------
set.seed(7)
n2 <- 60; d2 <- 5
X2 <- scale(matrix(rnorm(n2 * d2), n2, d2))
dd2 <- as.matrix(dist(X2))
R2 <- apply(dd2, 1, function(r) sort(r)[4])
den_new <- Vic_Den(X2, R2, d2)
den_old <- sapply(1:n2, function(x) {           # original line-12 formula
  size <- length(which(dd2[x, ] <= R2[x]))
  (size / R2[x]^d2)^(1 / d2)
})
rel <- max(abs(den_new - den_old) / den_old)
check(rel < 1e-12, sprintf("low-d (d=5) agreement with original formula (max rel diff %.2e)", rel))

cat(if (fail == 0) "\nALL PASS\n" else sprintf("\n%d FAILURE(S)\n", fail))
if (fail > 0) quit(status = 1)
