#!/usr/bin/env Rscript
# 93_validate_incremental_extended.R -- extended check that the incremental
# nearest-neighbour radius search (R/ccds/UN_CCD.R, 2026-09-23) changes no
# result. Extends 92b_validate_incremental_radi.R in two ways:
#
#  Part R (radii): nnccd.radi() new vs the verbatim earlier version
#    (R/ccds/UN_CCD_radi_recompute_reference.R) over 5 seeds x 7 dimensions
#    (2..100) x 6 generators (uniform, Gaussian, rounded values with many ties,
#    one cluster, five clusters, half the points duplicated) x both directions.
#    If the earlier version errors on a case, the new one must raise the same error.
#
#  Part E (end to end): the full UN-MCCD and SUN-MCCD outputs (score, cluster,
#    radii, connectivity) with the earlier radius search swapped in vs the new
#    one, on all 16 real data sets and the four above d = 21, ascending search
#    (the production setting); descending search as well on the real data sets
#    with n <= 350 and on simulation-style settings at d = 10, 20, 50.
#
# Appends one row per case to results/tr1/wp6_incr/93_validate_extended.csv.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))
source(here::here("R", "ccds", "UN_CCD_radi_recompute_reference.R"))
nnccd.radi.new <- nnccd.radi

OUTDIR <- here::here("revision_experiments/results/tr1/wp6_incr")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
OUT <- file.path(OUTDIR, "93_validate_extended.csv")
args <- commandArgs(trailingOnly = TRUE)
PARTS <- if (length(args)) strsplit(args[1], ",")[[1]] else c("R", "E")

log_row <- function(part, case, n, d, dir, method, same, note, t_old, t_new) {
  row <- data.frame(part = part, case = case, n = n, d = d, direction = dir, method = method,
                    identical = same, note = note, t_old = round(t_old, 2), t_new = round(t_new, 2),
                    stringsAsFactors = FALSE)
  write.table(row, OUT, sep = ",", row.names = FALSE, append = file.exists(OUT),
              col.names = !file.exists(OUT), qmethod = "double")
  cat(sprintf("[%s] %-28s n=%4d d=%3d %-7s %-8s identical=%s %s old=%.1fs new=%.1fs\n",
              part, case, n, d, dir, method, same, note, t_old, t_new))
  flush.console()
  isTRUE(same)
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "err"))
  list(v = v, t = proc.time()[["elapsed"]] - t0)
}

same_result <- function(a, b) {
  if (inherits(a, "err") || inherits(b, "err")) {
    return(list(same = inherits(a, "err") && inherits(b, "err") && identical(a$msg, b$msg),
                note = if (inherits(a, "err")) paste0("both error: ", substr(a$msg, 1, 40)) else "one side errored"))
  }
  list(same = identical(a, b), note = "")
}

ok <- TRUE

# ---------------------------------------------------------------- Part R
gen <- function(kind, seed, n, d) {
  set.seed(seed)
  n0 <- max(1, round(0.05 * n))
  box <- function(m, c) matrix(runif(m * d, c - 1, c + 1), m)
  out <- matrix(runif(n0 * d, -3, 7), n0)
  X <- switch(kind,
    uniform  = rbind(box((n - n0) %/% 2, 0), box(n - n0 - (n - n0) %/% 2, 4), out),
    gaussian = rbind(matrix(rnorm((n - n0) %/% 2 * d, 0, 0.4), ncol = d),
                     matrix(rnorm((n - n0 - (n - n0) %/% 2) * d, 4, 0.4), ncol = d), out),
    ties     = round(rbind(box((n - n0) %/% 2, 0), box(n - n0 - (n - n0) %/% 2, 4), out), 1),
    one      = rbind(box(n - n0, 0), out),
    five     = rbind(do.call(rbind, lapply(0:4, function(k) box((n - n0) %/% 5, 3 * k))),
                     matrix(runif((n - n0 - 5 * ((n - n0) %/% 5) + n0) * d, -2, 14), ncol = d)),
    dups     = { Z <- rbind(box((n - n0) %/% 4, 0), box((n - n0) %/% 4, 4), out)
                 rbind(Z, Z[seq_len(n - nrow(Z)), , drop = FALSE]) })
  X[seq_len(n), , drop = FALSE]
}

if ("R" %in% PARTS) {
  for (d in c(2, 3, 5, 10, 20, 50, 100)) {
    tab <- get_simul("NN", d, quant = nn_quant_label_paper_SUN(d), n = 150)
    for (kind in c("uniform", "gaussian", "ties", "one", "five", "dups")) for (seed in 1:5) {
      X <- gen(kind, 1000 * d + seed, 150, d)
      for (dir in c("ascend", "descend")) {
        a <- timed(nnccd.radi.recompute(X, method = dir, low.num = 3, simul = tab$simul)$R)
        b <- timed(nnccd.radi.new(X, method = dir, low.num = 3, simul = tab$simul)$R)
        s <- same_result(a$v, b$v)
        ok <- log_row("R", sprintf("%s_s%d", kind, seed), 150, d, dir, "radii", s$same, s$note, a$t, b$t) && ok
      }
    }
  }
}

# ---------------------------------------------------------------- Part E
set_impl <- function(f) assign("nnccd.radi", f, envir = .GlobalEnv)

# Sanity: the detectors must actually reach the swapped-in function.
set_impl(function(...) stop("OVERRIDE_HIT"))
hit <- tryCatch({ unmccd_method(X = gen("uniform", 1, 60, 3), d = 3); FALSE },
                error = function(e) grepl("OVERRIDE_HIT", conditionMessage(e)))
set_impl(nnccd.radi.new)
stopifnot("detectors do not call the global nnccd.radi" = hit)

run_det <- function(method, X, d, dir) {
  f <- if (method == "UN-MCCD") function() unmccd_method(X = X, d = d, method = dir)
       else function() sunmccd_method(X = X, d = d, method = dir, min.cls = 0.05)
  r <- f()
  r[c("score", "cluster", "radii", "connectivity", "unassigned_rows", "singleton_lost_rows", "quant_label")]
}

e2e <- function(case, X, d, dir, Y = NULL) {
  for (method in c("UN-MCCD", "SUN-MCCD")) {
    set_impl(nnccd.radi.recompute); a <- timed(run_det(method, X, d, dir))
    set_impl(nnccd.radi.new);       b <- timed(run_det(method, X, d, dir))
    s <- same_result(a$v, b$v)
    note <- s$note
    if (s$same && !inherits(b$v, "err") && !is.null(Y)) {
      m <- evaluate(Y, b$v$score, REAL_DATA_THRESHOLDS[[method]])
      note <- sprintf("BA=%.3f F2=%.3f", m[["BA"]], m[["F2"]])
    }
    ok <<- log_row("E", case, nrow(X), d, dir, method, s$same, note, a$t, b$t) && ok
  }
}

load_wp5 <- function(name) {
  df <- read.csv(here::here("revision_experiments/results/tr1/wp5/data", paste0(name, ".csv")))
  d <- ncol(df) - 1
  list(X = as.matrix(df[, seq_len(d)]), Y = df$label, d = d)
}

if ("E" %in% PARTS) {
  REAL <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli", "stamps",
            "WDBC", "pima", "shuffle", "vowels", "PenDigits", "waveform", "thyroid",
            "pageblocks", "wilt")
  for (ds in REAL) {
    dat <- load_real_dataset(ds)
    e2e(ds, dat$X, dat$d, "ascend", dat$Y)
    if (nrow(dat$X) <= 350) e2e(ds, dat$X, dat$d, "descend", dat$Y)
  }
  for (ds in c("letter", "mnist", "musk", "arrhythmia")) {
    dat <- load_wp5(ds)
    e2e(ds, dat$X, dat$d, "ascend", dat$Y)
  }
  for (d in c(10, 20, 50)) for (kind in c("uniform", "gaussian")) {
    X <- gen(kind, 77 + d, 200, d)
    e2e(sprintf("sim_%s", kind), X, d, "descend")
  }
}

set_impl(nnccd.radi.new)
cat(if (ok) "ALL_IDENTICAL\n" else "MISMATCH_FOUND\n")
