#!/usr/bin/env Rscript
# 96_wp6_direction_timing.R -- ascending vs descending radius search of UN-MCCD
# and SUN-MCCD, timed on the same data sets: the WP6 grid point n = 500, d = 10
# (cell 3 of tr1/91_wp6_runtime.R: same generator and seeds), 10 replicates,
# alternating order, single thread, idle machine, incremental search.
# Why: the runtime table uses the ascending search (the setting of every real
# data set); the d >= 10 simulation scripts used the descending one, and the
# manuscript states how much slower that is.
# Output (appended per measurement): results/tr1/wp6_incr/96_direction_raw.csv

suppressPackageStartupMessages(library(here))
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))

local({                                   # get_simul memoizer, as in 91/94
  .orig <- get_simul
  .cache <- new.env(parent = emptyenv())
  get_simul <<- function(variant = c("RK", "NN"), d, quant = NULL, n = NULL) {
    variant <- match.arg(variant)
    key <- paste(variant, d, quant, sep = "|")
    if (!exists(key, envir = .cache, inherits = FALSE)) {
      rm(list = ls(.cache, all.names = TRUE), envir = .cache)
      assign(key, .orig(variant, d, quant, n = NULL), envir = .cache)
    }
    tab <- get(key, envir = .cache, inherits = FALSE)
    check_simul_extent(tab$simul, variant, d, n, tab$file)
    tab
  }
})

gen_uniform_wp6 <- function(seed, n, d, cont = 0.05) {       # verbatim from 91
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3
  mu1 <- rep(3, d); mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu  <- (mu1 + mu2) / 2
  n1 <- round(n * (1 - cont) * 0.5)
  n0 <- round(n * cont)
  n2 <- n - n1 - n0
  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) +
             matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) +
             matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i <- 0L; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis && sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier <- rbind(outlier, temp); i <- i + 1L
    }
  }
  if (!is.null(outlier)) rownames(outlier) <- NULL
  X <- rbind(data1, data2, outlier)
  colnames(X) <- paste0("V", seq_len(d))
  list(X = X, Y = c(rep(1L, n1 + n2), rep(0L, n0)))
}
seed_for <- function(cell_index, rep) 60000L + cell_index * 100L + rep
N <- 500; D <- 10; CELL <- 3L; REPS <- 10L

OUT <- here::here("revision_experiments/results/tr1/wp6_incr/96_direction_raw.csv")
for (method in c("UN-MCCD", "SUN-MCCD")) {
  q <- if (method == "UN-MCCD") nn_quant_label_paper_UN(D) else nn_quant_label_paper_SUN(D)
  invisible(get_simul("NN", D, quant = q, n = N))
  for (r in seq_len(REPS)) {
    dat <- gen_uniform_wp6(seed_for(CELL, r), N, D)
    for (dir in if (r %% 2 == 1) c("ascend", "descend") else c("descend", "ascend")) {
      gc_before <- gc(reset = TRUE)
      t0 <- proc.time()[["elapsed"]]
      res <- if (method == "UN-MCCD") unmccd_method(X = dat$X, d = D, method = dir)
             else sunmccd_method(X = dat$X, d = D, method = dir, min.cls = 0.05)
      el <- proc.time()[["elapsed"]] - t0
      gc_after <- gc(reset = FALSE)
      row <- data.frame(method = method, n = N, d = D, rep = r, direction = dir,
                        time_s = round(el, 4),
                        mem_delta_mb = round(sum(gc_after[, 6]) - sum(gc_before[, 2]), 3))
      write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
                  col.names = !file.exists(OUT))
      cat(sprintf("%-8s rep %2d %-7s %.2fs\n", method, r, dir, el)); flush.console()
    }
  }
}
cat("DIRECTION_DONE\n")
