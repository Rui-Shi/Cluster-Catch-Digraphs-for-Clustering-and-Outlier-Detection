#!/usr/bin/env Rscript
# 94_wp6_ab_old_vs_new.R -- paired timing of the earlier (per-step recompute)
# and the incremental nearest-neighbour radius search inside UN-MCCD and
# SUN-MCCD, on the WP6 runtime grid (same generator, cells and seeds as
# tr1/91_wp6_runtime.R), ascending search, single thread, idle machine.
#
# Why: the re-measurement in results/tr1/wp6_incr (new code) came out 6-27%
# slower than results/tr1/wp6 (old code, measured on another day). A paired
# design removes day-to-day machine drift: every replicate draws ONE data set
# and times BOTH implementations on it, in an order that alternates between
# replicates. It also times the radius search alone, to show its share.
#
# Output (appended per measurement): results/tr1/wp6_incr/94_ab_raw.csv
# Usage: Rscript 94_wp6_ab_old_vs_new.R [--reps=10]

suppressPackageStartupMessages(library(here))
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))
source(here::here("R/ccds/UN_CCD_radi_recompute_reference.R"))
nnccd.radi.new <- nnccd.radi

args <- commandArgs(trailingOnly = TRUE)
REPS <- { h <- grep("^--reps=", args, value = TRUE); if (length(h)) as.integer(sub("^--reps=", "", h)) else 10L }

# get_simul memoizer, identical to 91_wp6_runtime.R (keeps table loads out of the timed span)
local({
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

# generator, cells and seeds copied verbatim from 91_wp6_runtime.R
gen_uniform_wp6 <- function(seed, n, d, cont = 0.05) {
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
  Y <- c(rep(1L, n1 + n2), rep(0L, n0))
  list(X = X, Y = Y)
}
CELLS <- data.frame(cell_index = 1:8,
                    n = c(100, 250, 500, 1000, 2000, 500, 500, 500),
                    d = c(10, 10, 10, 10, 10, 5, 50, 100))
seed_for <- function(cell_index, rep) 60000L + cell_index * 100L + rep

OUT <- here::here("revision_experiments/results/tr1/wp6_incr/94_ab_raw.csv")
done <- if (file.exists(OUT)) read.csv(OUT, stringsAsFactors = FALSE) else NULL
append_row <- function(row) {
  write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
              col.names = !file.exists(OUT), qmethod = "double")
}
is_done <- function(ci, method, impl, part, r) {
  !is.null(done) && any(done$cell_index == ci & done$method == method &
                        done$impl == impl & done$part == part & done$rep == r)
}

call_det <- function(method, X, d) {
  if (method == "UN-MCCD") unmccd_method(X = X, d = d, method = "ascend")
  else sunmccd_method(X = X, d = d, method = "ascend", min.cls = 0.05)
}
measure <- function(fn) {
  gc_before <- gc(reset = TRUE)
  t0 <- proc.time()[["elapsed"]]
  v <- fn()
  el <- proc.time()[["elapsed"]] - t0
  gc_after <- gc(reset = FALSE)
  list(v = v, time = el, mem_delta = sum(gc_after[, 6]) - sum(gc_before[, 2]))
}

t_start <- Sys.time()
for (k in seq_len(nrow(CELLS))) {
  cell <- CELLS[k, ]
  for (method in c("UN-MCCD", "SUN-MCCD")) {
    q <- if (method == "UN-MCCD") nn_quant_label_paper_UN(cell$d) else nn_quant_label_paper_SUN(cell$d)
    invisible(get_simul("NN", cell$d, quant = q, n = cell$n))          # warm the cache
    for (r in seq_len(REPS)) {
      dat <- gen_uniform_wp6(seed_for(cell$cell_index, r), cell$n, cell$d)
      order_impl <- if (r %% 2 == 1) c("old", "new") else c("new", "old")
      res <- list()
      for (impl in order_impl) {
        assign("nnccd.radi", if (impl == "old") nnccd.radi.recompute else nnccd.radi.new, envir = .GlobalEnv)
        tab <- get_simul("NN", cell$d, quant = q, n = cell$n)
        for (part in c("detector", "radius_search")) {
          if (is_done(cell$cell_index, method, impl, part, r)) next
          m <- if (part == "detector") measure(function() call_det(method, dat$X, cell$d))
               else measure(function() nnccd.radi(dat$X, method = "ascend", low.num = 3, simul = tab$simul))
          if (part == "detector") res[[impl]] <- m$v$score
          append_row(data.frame(cell_index = cell$cell_index, n = cell$n, d = cell$d, method = method,
                                rep = r, impl = impl, order = paste(order_impl, collapse = ">"),
                                part = part, time_s = round(m$time, 4),
                                mem_delta_mb = round(m$mem_delta, 3), stringsAsFactors = FALSE))
        }
      }
      same <- if (length(res) == 2) identical(res$old, res$new) else NA
      cat(sprintf("n=%4d d=%3d %-8s rep %2d  scores identical=%s\n", cell$n, cell$d, method, r, same))
      flush.console()
    }
  }
}
assign("nnccd.radi", nnccd.radi.new, envir = .GlobalEnv)
cat(sprintf("AB_DONE in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
