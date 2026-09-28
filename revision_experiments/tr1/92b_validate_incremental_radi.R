#!/usr/bin/env Rscript
# 92_validate_incremental_radi.R -- the incremental nearest-neighbour radius
# search in R/ccds/UN_CCD.R (2026-09-23) must return EXACTLY the radii of the
# earlier per-step recomputation (R/ccds/UN_CCD_radi_recompute_reference.R).
# Compares identical() radii, both search directions, on synthetic two-cluster
# data and on real data sets, with the production quantile tables, and reports
# the speed-up. Appends one row per case to results/tr1/wp6_incr/92_validate.csv.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))
source(here::here("R", "ccds", "UN_CCD_radi_recompute_reference.R"))

OUTDIR <- here::here("revision_experiments/results/tr1/wp6_incr")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, "92_validate.csv")

run_case <- function(tag, X, d, dir, scores = FALSE, low.num = 3) {
  n <- nrow(X)
  tab <- get_simul("NN", d, quant = nn_quant_label_paper_SUN(d), n = n)
  t0 <- proc.time()[["elapsed"]]
  a <- nnccd.radi.recompute(X, method = dir, low.num = low.num, simul = tab$simul, scores = scores)$R
  t_old <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  b <- nnccd.radi(X, method = dir, low.num = low.num, simul = tab$simul, scores = scores)$R
  t_new <- proc.time()[["elapsed"]] - t0
  row <- data.frame(case = tag, n = n, d = d, direction = dir, scores = scores,
                    identical = identical(a, b), n_diff = sum(a != b),
                    zero_radii = sum(b == 0), t_old = t_old, t_new = t_new,
                    stringsAsFactors = FALSE)
  write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
              col.names = !file.exists(OUT))
  cat(sprintf("%-22s n=%4d d=%3d %-7s scores=%-5s identical=%s  old=%.1fs new=%.1fs\n",
              tag, n, d, dir, scores, identical(a, b), t_old, t_new))
  flush.console()
  identical(a, b)
}

gen2 <- function(seed, n, d, gauss = FALSE) {
  set.seed(seed)
  n1 <- floor(n * 0.475); n2 <- floor(n * 0.475); n0 <- n - n1 - n2
  g <- function(m, c) if (gauss) matrix(rnorm(m * d, c, 0.3), m) else matrix(runif(m * d, c - 1, c + 1), m)
  rbind(g(n1, 0), g(n2, 4), matrix(runif(n0 * d, -3, 7), n0))
}

ok <- TRUE
for (d in c(2, 3, 5, 10, 20)) for (gauss in c(FALSE, TRUE)) for (dir in c("ascend", "descend")) {
  X <- gen2(100 + d, 200, d, gauss)
  ok <- run_case(sprintf("%s_d%d", if (gauss) "gauss" else "unif", d), X, d, dir) && ok
}
X <- gen2(7, 200, 5); ok <- run_case("unif_d5_scores", X, 5, "ascend", scores = TRUE) && ok
ok <- run_case("unif_d5_scores", X, 5, "descend", scores = TRUE) && ok
ok <- run_case("unif_d5_lownum2", X, 5, "ascend", low.num = 4) && ok
X <- rbind(X, X[1:5, ]); ok <- run_case("unif_d5_duplicates", X, 5, "ascend") && ok
ok <- run_case("unif_d5_duplicates", X, 5, "descend") && ok
for (ds in c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli", "stamps", "WDBC")) {
  dat <- load_real_dataset(ds)
  for (dir in c("ascend", "descend")) ok <- run_case(ds, dat$X, dat$d, dir) && ok
}
X <- gen2(11, 1000, 10); ok <- run_case("unif_d10_n1000", X, 10, "ascend") && ok
cat(if (ok) "ALL_IDENTICAL\n" else "MISMATCH_FOUND\n")
