#!/usr/bin/env Rscript
# revision_experiments/49_wp2a_rk_saturation.R
#
# WP2a (Neurocomputing revision): does the RK spatial-randomness test's
# quantile envelope saturate with dimension, and does that saturation explain
# a loss of alpha-sensitivity in the U-MCCD / SU-MCCD detectors?  This is the
# RK counterpart to a saturation claim being re-derived on the NN side by a
# separate agent (see FINDINGS.md / 42_wp2a_alpha_sweep.R for the NN-side
# duplicate-table problem that motivated this cross-check).
#
# Part 1 (table analysis, zero detector runs): for every RK dimension on
#   disk, derive the quantile vector at alpha in {0.10, 0.05, 0.01, 0.001}
#   (i.e. probability q in {0.90, 0.95, 0.99, 0.999}) straight from the raw
#   Monte-Carlo draws (Kest.m) already stored in each RK-test-simul_*.RData
#   file, and report the zero-quantile fraction, near-zero fraction, and the
#   median/IQR of the non-zero entries; also report how much the vector moves
#   between alpha=0.01 and alpha=0.001, normalised by the median non-zero
#   value.
#
# Part 2 (behavioural): run U-MCCD and SU-MCCD at n=200 over d in
#   {2,3,5,10,20,35,50,100} and alpha in {0.10,0.05,0.01,0.001}, on both the
#   uniform-cluster and Gaussian-cluster generators from
#   revision_experiments/09_wp3_synthetic.R (gen_uniform / gen_gaussian),
#   generalised from their original fixed d=10,n=500 to arbitrary (d, n) --
#   every other constant (cont=0.05, cls_dis=3, otl_dis=2, r_min=0.7,
#   r_max=1.3, noise_level=0.01 for the Gaussian setting) is kept byte-for-
#   byte identical to the original script. min.cls is passed to SU-MCCD as a
#   PROPORTION of n (0.0625), never a count -- see CLAUDE.md / harness.R WP0
#   notes on that unit trap.
#
# RK-alpha derivation machinery (Part 2's get_simul() shadow) follows the
# same design as 42_wp2a_alpha_sweep.R: get_simul() is overridden AFTER
# harness.R / wp0_mccd_methods.R are sourced (both source() into globalenv,
# so an earlier override would be clobbered); the shadow loads ONE base RK
# file per d (the paper-default token from rk_quant_label_paper()), derives
# the requested probability's $quan entry directly from $Kest.m, and -- before
# trusting the derivation -- reproduces the file's OWN stored $quan entry
# from $Kest.m and aborts if that reproduction is not bit-identical. This
# script only READS harness.R / wp0_mccd_methods.R / the .RData tables; it
# never edits them, per the hard read-only rule for this task.
#
# Usage:
#   Rscript 49_wp2a_rk_saturation.R --part1                  # table analysis
#   Rscript 49_wp2a_rk_saturation.R --part2 --control         # low-d control (d=2,3)
#   Rscript 49_wp2a_rk_saturation.R --part2 [d1,d2,...]        # behavioural grid, chunked by d
#   Rscript 49_wp2a_rk_saturation.R --summary                  # Part 2 aggregation + join with Part 1
#
# Checkpointed on (setting, d, alpha_token, rep): safe to re-invoke after a
# PowerShell 600s timeout.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- only now is it safe to define and install the get_simul shadow --------

RK_DIR <- RK_QUANT_TABLE_DIR
RESULTS_DIR <- here::here("revision_experiments/results/tr1")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# SAT_SUFFIX lets a second, parallelised run write to a physically different
# file so it never collides with an earlier single-threaded run still
# appending to the canonical name (Windows file writes from two independent
# processes are not safely interleaved). Set WP2A_SAT_SUFFIX=_v2 etc.
# Merge/rename to the canonical name once the earlier writer has exited.
SAT_SUFFIX   <- Sys.getenv("WP2A_SAT_SUFFIX", unset = "")
TABLE_CSV    <- file.path(RESULTS_DIR, "wp2a_rk_table_degeneracy.csv")
SAT_CSV      <- file.path(RESULTS_DIR, sprintf("wp2a_rk_saturation%s.csv", SAT_SUFFIX))
SAT_SUM_CSV  <- file.path(RESULTS_DIR, sprintf("wp2a_rk_saturation_summary%s.csv", SAT_SUFFIX))
FINDINGS_MD  <- file.path(RESULTS_DIR, "wp2a_rk_saturation_findings.md")

ALPHAS <- c(0.10, 0.05, 0.01, 0.001)          # alpha levels requested
PROBS  <- 1 - ALPHAS                           # {0.90, 0.95, 0.99, 0.999}
names(PROBS) <- sprintf("a%s", ALPHAS)
# Exact alpha value per name key, used everywhere instead of recomputing
# 1 - (1 - alpha) -- that round trip is not guaranteed bit-identical to the
# literal alpha value in IEEE double arithmetic, which silently breaks any
# `== 0.01`-style comparison downstream (confirmed empirically: a direct test
# run produced NA move_norm_* columns from exactly this bug before this fix).
ALPHA_OF_NM <- setNames(ALPHAS, names(PROBS))

# The full set of RK dimensions on disk (2..35, 50, 100), read directly from
# filenames rather than hardcoded, so a missing/extra file is visible.
disk_rk_dims <- function() {
  f <- list.files(RK_DIR, pattern = "^RK-test-simul_[0-9]+d_[0-9]+%\\.RData$")
  as.integer(sub("^RK-test-simul_([0-9]+)d_.*$", "\\1", f))
}
ALL_RK_D <- sort(unique(disk_rk_dims()))

rk_files_for_d <- function(d) {
  list.files(RK_DIR, pattern = sprintf("^RK-test-simul_%dd_[0-9]+%%\\.RData$", d))
}

# ---------------------------------------------------------------------------
# Part 1: pure table analysis. One base file per d (paper-default token from
# wp0_mccd_methods.R's rk_quant_label_paper(); falls back to whatever token
# is on disk if the default isn't present -- none of the 36 dims hit that
# fallback, verified below).
# ---------------------------------------------------------------------------
rk_derive_quan_vec <- function(Kest.m, q, chunk = 5000L) {
  # Same formula as Kest.R:406-411 / 42_wp2a_alpha_sweep.R's rk_derive_quan(),
  # but returns the FLAT length-ncol(Kest.m) vector (no reshape to matrix) --
  # that flat vector is exactly what Part 1 needs to characterise.
  nc <- ncol(Kest.m)
  out <- numeric(nc)
  i <- 1L
  while (i <= nc) {
    j <- min(i + chunk - 1L, nc)
    out[i:j] <- apply(Kest.m[, i:j, drop = FALSE], 2, quantile, probs = q)
    i <- j + 1L
  }
  out
}

part1_one_dim <- function(d) {
  files <- rk_files_for_d(d)
  base_label <- rk_quant_label_paper(d)
  base_file <- sprintf("RK-test-simul_%dd_%s%%.RData", d, base_label)
  if (!(base_file %in% files)) {
    base_file <- files[1]
    base_label <- sub(sprintf("^RK-test-simul_%dd_([0-9]+)%%\\.RData$", d), "\\1", base_file)
    cat(sprintf("  [d=%d] paper-default token not on disk; falling back to %s\n", d, base_file))
  }
  path <- file.path(RK_DIR, base_file)
  t0 <- Sys.time()
  e <- new.env(); load(path, envir = e)
  s <- get("simul", envir = e)
  Kest.m <- s[["Kest.m"]]
  stopifnot(!is.null(Kest.m))
  nc <- ncol(Kest.m)

  # Guard: reproduce the file's OWN stored quan entry from Kest.m before
  # trusting any derived probability at this d.
  own_q <- as.numeric(paste0("0.", base_label))
  own_key <- as.character(own_q)
  stored <- s[["quan"]][[own_key]]
  stopifnot(!is.null(stored))
  m_rows <- if (is.null(dim(stored))) NA_integer_ else nrow(stored)
  stored_flat <- as.numeric(stored)  # matrix or vector, flattened column-major
  redone_flat <- rk_derive_quan_vec(Kest.m, own_q)
  # stored may be reshaped to a matrix(nrow=m) from the same flat vector;
  # matrix() with byrow=FALSE preserves column-major order, so flattening
  # both back to vectors must match exactly regardless of shape.
  guard_ok <- isTRUE(all.equal(redone_flat, stored_flat, tolerance = 0))
  if (!guard_ok) {
    stop(sprintf("part1_one_dim(): DERIVATION GUARD FAILED at d=%d (%s) -- recomputed quan[['%s']] does not reproduce the stored vector",
                 d, base_file, own_key))
  }

  rows <- list()
  prev_vec <- NULL
  for (nm in names(PROBS)) {
    q <- PROBS[[nm]]
    v <- if (isTRUE(all.equal(q, own_q, tolerance = 0))) stored_flat else rk_derive_quan_vec(Kest.m, q)
    n_zero <- sum(v == 0)
    n_tiny <- sum(v < 1e-12)
    nz <- v[v > 1e-12]
    med_nz <- if (length(nz)) median(nz) else NA_real_
    q1 <- if (length(nz)) quantile(nz, 0.25, names = FALSE) else NA_real_
    q3 <- if (length(nz)) quantile(nz, 0.75, names = FALSE) else NA_real_
    rows[[nm]] <- data.frame(
      d = d, alpha = ALPHA_OF_NM[[nm]], prob = q, base_file = base_file, base_token = base_label,
      n_entries = nc, frac_zero = n_zero / nc, frac_below_1e12 = n_tiny / nc,
      median_nonzero = med_nz, iqr_nonzero = q3 - q1, q1_nonzero = q1, q3_nonzero = q3,
      stringsAsFactors = FALSE)
    if (nm == "a0.01") a01 <- v
    if (nm == "a0.001") a0001 <- v
  }
  out <- do.call(rbind, rows)

  # Movement between alpha=0.01 and alpha=0.001, normalised by the median
  # non-zero entry of the alpha=0.01 vector (the tighter of the two commonly-
  # used levels; using the looser one as the normaliser is an equally
  # defensible choice and is reported too for robustness).
  mad_diff <- mean(abs(a01 - a0001))
  med01 <- out$median_nonzero[out$alpha == 0.01]
  med0001 <- out$median_nonzero[out$alpha == 0.001]
  out$move_mad_a01_a0001 <- mad_diff
  out$move_norm_by_med01 <- if (length(med01) && !is.na(med01) && med01 > 0) mad_diff / med01 else NA_real_
  out$move_norm_by_med0001 <- if (length(med0001) && !is.na(med0001) && med0001 > 0) mad_diff / med0001 else NA_real_
  out$guard_ok <- guard_ok
  out$m_rows_stored <- m_rows
  out$wall_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rm(Kest.m, s, e, a01, a0001); gc(FALSE)
  out
}

run_part1 <- function(dims = ALL_RK_D) {
  cat(sprintf("Part 1: %d RK dimensions on disk: %s\n", length(dims), paste(dims, collapse = ",")))
  done <- if (file.exists(TABLE_CSV)) unique(read.csv(TABLE_CSV)$d) else integer(0)
  pend <- setdiff(dims, done)
  cat(sprintf("Part 1: %d already done, %d pending: %s\n", length(done), length(pend), paste(pend, collapse = ",")))
  for (d in pend) {
    cat(sprintf("\n[Part1] d=%d ...\n", d))
    res <- part1_one_dim(d)
    write.table(res, TABLE_CSV, sep = ",", row.names = FALSE,
                col.names = !file.exists(TABLE_CSV), append = file.exists(TABLE_CSV))
    cat(sprintf("  d=%-3d frac_zero(a=0.10..0.001)=%s  wall=%.1fs\n",
                d, paste(sprintf("%.3f", res$frac_zero), collapse = "/"), res$wall_sec[1]))
    flush.console()
  }
  cat("\nPart 1 done.\n")
}

# ---------------------------------------------------------------------------
# Part 2: behavioural grid. get_simul() shadow -- derives the RK quan entry
# for an arbitrary probability from one base file's Kest.m, guarded exactly
# as in Part 1 / 42_wp2a_alpha_sweep.R.
# ---------------------------------------------------------------------------
WP2A <- new.env(parent = emptyenv())
WP2A$cache_d <- NA_integer_
WP2A$cache_label <- NA_character_
WP2A$cache_simul <- NULL
WP2A$cache_file <- NA_character_

orig_get_simul <- get_simul

rk_base_table <- function(d) {
  base_label <- rk_quant_label_paper(d)
  files <- rk_files_for_d(d)
  base_file <- sprintf("RK-test-simul_%dd_%s%%.RData", d, base_label)
  if (!(base_file %in% files)) { base_file <- files[1]; base_label <- sub(sprintf("^RK-test-simul_%dd_([0-9]+)%%\\.RData$", d), "\\1", base_file) }
  if (!identical(WP2A$cache_d, d) || !identical(WP2A$cache_label, base_label)) {
    path <- file.path(RK_DIR, base_file)
    t0 <- Sys.time()
    e <- new.env(); load(path, envir = e)
    s <- get("simul", envir = e)
    if (is.null(s$Kest.m)) stop("rk_base_table(): no $Kest.m in ", path)
    own_q <- as.numeric(paste0("0.", base_label)); own_key <- as.character(own_q)
    stored <- s$quan[[own_key]]
    if (is.null(stored)) stop("rk_base_table(): stored $quan has no key ", own_key)
    m <- nrow(stored)
    redone <- rk_derive_quan_vec(s$Kest.m, own_q)
    if (!isTRUE(all.equal(redone, as.numeric(stored), tolerance = 0))) {
      stop(sprintf("rk_base_table(): DERIVATION GUARD FAILED for %s", path))
    }
    WP2A$cache_simul <- s; WP2A$cache_d <- d; WP2A$cache_label <- base_label
    WP2A$cache_file <- path; WP2A$cache_m <- m
    cat(sprintf("    [rk base] loaded %s in %.1fs (guard OK)\n", basename(path),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  WP2A$cache_simul
}

get_simul_shadow <- function(variant = c("RK", "NN"), d, quant = NULL) {
  variant <- match.arg(variant)
  if (variant == "RK" && !is.null(quant) && is.numeric(quant)) {
    # quant is a bare probability (e.g. 0.95) rather than a file-label token
    # -> derive it from the cached base file for this d.
    s <- rk_base_table(d)
    key <- as.character(quant)
    if (is.null(s$quan[[key]])) {
      s$quan[[key]] <- matrix(rk_derive_quan_vec(s$Kest.m, quant), nrow = WP2A$cache_m)
      WP2A$cache_simul <- s
    }
    return(list(simul = s, quant = quant, quant_label = sub("^0\\.", "", format(quant, trim = TRUE, scientific = FALSE)),
                file = WP2A$cache_file))
  }
  orig_get_simul(variant, d, quant)
}
assign("get_simul", get_simul_shadow, envir = globalenv())

# umccd_method()/sumccd_method() call get_simul("RK", d, quant = q_label)
# where q_label defaults to a STRING token. We instead call the underlying
# detector directly with a bare probability via the shadow, by passing a
# NUMERIC quant through a thin re-implementation of the two wrappers (the
# originals coerce quant to a file-label string internally via
# rk_quant_label_paper()/q_label logic when NULL, but pass a non-NULL
# argument straight to get_simul() unchanged -- confirmed by reading
# umccd_method/sumccd_method above: `tab <- get_simul("RK", d, quant = q_label)`
# with q_label = quant when quant is supplied. A numeric quant is therefore
# accepted unchanged by the wrapper and only interpreted by our shadow.)
umccd_alpha <- function(X, d, prob) {
  tr <- umccd_method(X = X, d = d, quant = prob)
  tr
}
sumccd_alpha <- function(X, d, prob, min.cls) {
  tr <- sumccd_method(X = X, d = d, quant = prob, min.cls = min.cls)
  tr
}

# ---------------------------------------------------------------------------
# Synthetic generators -- generalised from 09_wp3_synthetic.R's gen_gaussian /
# gen_uniform (which hardcode d=10, n=500). All other constants kept
# identical to that script; only d and n become parameters. Verified
# rpoisball.unit()/mvrnorm() are already in globalenv via harness.R's source
# chain (RKCCD_OOS_IOS.R -> ... -> Kest.R defines rpoisball.unit; MASS is
# attached by harness.R for mvrnorm).
# ---------------------------------------------------------------------------
gen_gaussian_dn <- function(seed, d, n) {
  cont <- 0.05
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3

  mu1 <- rep(3, d)
  mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu <- rbind(mu1, mu2); mu <- apply(mu, 2, mean)

  n1 <- round(n * (1 - cont) * 0.5)
  n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)

  noise_level <- 0.01
  sigma <- 1 / sqrt(qchisq(1 - noise_level, d))

  set.seed(seed)
  data1 <- mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 <- mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    r1 <- sqrt(sum((temp - mu1)^2)); r2 <- sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) { outlier <- rbind(outlier, temp); i <- i + 1 }
  }
  rownames(outlier) <- NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

gen_uniform_dn <- function(seed, d, n) {
  cont <- 0.05
  cls_dis <- 3; otl_dis <- 2; r_min <- 0.7; r_max <- 1.3

  mu1 <- rep(3, d)
  mu2 <- c(3 + cls_dis, rep(3, d - 1))
  mu <- rbind(mu1, mu2); mu <- apply(mu, 2, mean)

  n1 <- round(n * (1 - cont) * 0.5)
  n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)

  set.seed(seed)
  data1 <- rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 <- rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i <- 0; outlier <- NULL
  while (i < n0) {
    temp <- rpoisball.unit(1, d) * 5 + mu
    r1 <- sqrt(sum((temp - mu1)^2)); r2 <- sqrt(sum((temp - mu2)^2))
    if (r1 > otl_dis & r2 > otl_dis) { outlier <- rbind(outlier, temp); i <- i + 1 }
  }
  rownames(outlier) <- NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

GENERATORS <- list(gaussian = gen_gaussian_dn, uniform = gen_uniform_dn)

S_MIN <- 0.0625
N_SYN <- 200L

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
parse_idx <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else as.integer(strsplit(s, ";")[[1]])

has_result2 <- function(csv_path, keys) {
  if (!file.exists(csv_path)) return(FALSE)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  for (k in names(keys)) if (!(k %in% names(df))) return(FALSE)
  match_all <- rep(TRUE, nrow(df))
  for (k in names(keys)) match_all <- match_all & (as.character(df[[k]]) == as.character(keys[[k]]))
  any(match_all)
}

# Pure computation, no I/O: returns a data.frame of up to 8 rows (2 methods x
# 4 alpha levels) for one (setting, d, rep). Split out from the writer so a
# parallel worker can return rows to the master instead of each worker
# appending to SAT_CSV directly -- concurrent small appends from independent
# OS processes are not guaranteed atomic/non-interleaved on Windows, so only
# the master process ever writes.
run_part2_one_rows <- function(setting, d, rep_id, base_seed, skip_existing = TRUE) {
  gen <- GENERATORS[[setting]]
  seed <- base_seed + rep_id
  dat <- gen(seed, d, N_SYN)
  X <- dat$X; n <- dat$n; n0 <- dat$n0
  Y <- c(rep(1, n - n0), rep(0, n0))

  out_rows <- list()
  for (nm in names(PROBS)) {
    prob <- unname(PROBS[[nm]])
    alpha <- unname(ALPHA_OF_NM[[nm]])
    keys <- c(setting = setting, d = as.character(d), rep = as.character(rep_id), alpha = as.character(alpha))
    if (skip_existing && isTRUE(has_result2(SAT_CSV, keys))) next

    t0 <- Sys.time()
    out_u <- tryCatch({
      r <- umccd_alpha(X, d, prob)
      list(score = r$score, status = "ok", note = NA_character_)
    }, error = function(e) list(score = rep(NA_real_, n), status = "error", note = conditionMessage(e)))

    out_s <- tryCatch({
      r <- sumccd_alpha(X, d, prob, min.cls = S_MIN)
      list(score = r$score, status = "ok", note = NA_character_)
    }, error = function(e) list(score = rep(NA_real_, n), status = "error", note = conditionMessage(e)))
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    for (meth_nm in c("U-MCCD", "SU-MCCD")) {
      out <- if (meth_nm == "U-MCCD") out_u else out_s
      if (out$status == "ok") {
        m <- evaluate(Y, out$score, 0.5)
        flagged <- which(out$score > 0.5)
      } else {
        m <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
        flagged <- integer(0)
      }
      out_rows[[length(out_rows) + 1]] <- data.frame(
        setting = setting, method = meth_nm, d = d, n = n, n0 = n0,
        rep = rep_id, seed = seed, alpha = alpha, prob = prob,
        TPR = unname(m[["TPR"]]), TNR = unname(m[["TNR"]]), BA = unname(m[["BA"]]), F2 = unname(m[["F2"]]),
        n_flagged = length(flagged), flagged_idx = collapse_idx(flagged),
        status = out$status, note = out$note, wall_sec = wall,
        timestamp = format(Sys.time()), stringsAsFactors = FALSE)
    }
  }
  if (!length(out_rows)) return(NULL)
  do.call(rbind, out_rows)
}

write_rows <- function(rows) {
  if (is.null(rows) || !nrow(rows)) return(invisible(NULL))
  write.table(rows, SAT_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(SAT_CSV), append = file.exists(SAT_CSV))
}

# Backward-compatible single-process entry point (writes directly).
run_part2_one <- function(setting, d, rep_id, base_seed) {
  write_rows(run_part2_one_rows(setting, d, rep_id, base_seed))
}

run_part2 <- function(dims, n_reps = 30L, settings = c("gaussian", "uniform")) {
  base_seeds <- c(gaussian = 123L, uniform = 123L)
  for (d in dims) {
    cat(sprintf("\n[Part2] d=%d ...\n", d))
    for (setting in settings) {
      for (r in seq_len(n_reps)) {
        t0 <- Sys.time()
        run_part2_one(setting, d, r, base_seeds[[setting]])
        cat(sprintf("  [%s] d=%-3d rep=%-3d done in %.1fs\n", setting, d, r,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        flush.console()
      }
    }
    WP2A$cache_d <- NA_integer_; WP2A$cache_simul <- NULL
    invisible(gc(FALSE))
  }
}

# ---------------------------------------------------------------------------
# Parallel driver. Single-call detector cost is 6-40s per alpha level
# (measured directly: umccd_alpha at d=2,n=200 took 6.3/5.7/14.8/37.9s at
# alpha=0.10/0.05/0.01/0.001; sumccd_method adds a comparable amount on top).
# A sequential run of the requested 20-50 reps x 8 dims x 2 settings would
# take hours; this driver farms REPLICATES (not alpha levels -- all 4 alpha
# levels for one rep share one cached base RK table, so splitting alpha
# across workers would multiply the expensive file-load+derive step instead
# of amortising it) out to a bounded-size PSOCK cluster, one dimension at a
# time so the per-dimension RK cache reset in run_part2() is preserved.
# Workers each independently source harness.R/wp0_mccd_methods.R and install
# their own get_simul() shadow -- the exact pattern 09_wp3_synthetic.R uses
# for its own cluster workers.
run_part2_parallel <- function(dims, n_reps = 20L, settings = c("gaussian", "uniform"),
                                cores = 5L) {
  cores <- max(1L, min(cores, n_reps))
  repo_root <- here::here()
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, "repo_root", envir = environment())
  invisible(parallel::clusterEvalQ(cl, {
    setwd(repo_root)
    suppressMessages(library(here))
    # Re-sourcing this same file on each worker: commandArgs(trailingOnly=TRUE)
    # is empty in a parallel worker (its args are MASTER=/PORT=/... from
    # parallel's own PSOCK bootstrap, not --part1/--part2), so the dispatch
    # block at the bottom of this file falls through to the harmless "Usage:"
    # print and every function/variable/shadow-install above it still runs.
    suppressMessages(source(here::here("revision_experiments", "49_wp2a_rk_saturation.R")))
    TRUE
  }))
  base_seeds <- c(gaussian = 123L, uniform = 123L)
  for (d in dims) {
    cat(sprintf("\n[Part2-parallel] d=%d (cores=%d) ...\n", d, cores))
    for (setting in settings) {
      t0 <- Sys.time()
      # Only ship pending reps -- has_result2() against the MASTER's view of
      # SAT_CSV, taken once per (d,setting) before dispatch (a worker crash
      # mid-dimension just means its reps get skipped as "pending" next call).
      done_reps <- integer(0)
      if (file.exists(SAT_CSV)) {
        existing <- tryCatch(read.csv(SAT_CSV, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(existing) && nrow(existing)) {
          es <- existing[existing$setting == setting & existing$d == d, ]
          for (r in unique(es$rep)) {
            if (length(unique(es$alpha[es$rep == r])) >= length(PROBS) &&
                length(unique(es$method[es$rep == r])) >= 2) done_reps <- c(done_reps, r)
          }
        }
      }
      pend <- setdiff(seq_len(n_reps), done_reps)
      if (!length(pend)) { cat(sprintf("  [%s] d=%d: all %d reps already recorded\n", setting, d, n_reps)); next }
      cat(sprintf("  [%s] d=%d: dispatching %d/%d pending reps to %d workers\n", setting, d, length(pend), n_reps, cores))
      # Workers COMPUTE and return rows; only the master WRITES (see
      # run_part2_one_rows()'s docstring for why -- concurrent appends from
      # independent OS processes are not safely interleaved on Windows).
      res_list <- parallel::parLapply(cl, pend, function(r, setting, d, seed0) {
        run_part2_one_rows(setting, d, r, seed0, skip_existing = FALSE)
      }, setting = setting, d = d, seed0 = base_seeds[[setting]])
      for (rows in res_list) write_rows(rows)
      cat(sprintf("  [%s] d=%-3d %d reps done in %.1fs\n", setting, d, length(pend),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      flush.console()
    }
    invisible(parallel::clusterEvalQ(cl, { WP2A$cache_d <- NA_integer_; WP2A$cache_simul <- NULL; gc(FALSE); TRUE }))
  }
}

run_control <- function(n_reps = 10L, parallel_cores = 5L) {
  cat("==== CONTROL: is the alpha effect clearly non-zero at d=2,3? ====\n")
  if (parallel_cores > 1) run_part2_parallel(c(2, 3), n_reps = n_reps, cores = parallel_cores) else run_part2(c(2, 3), n_reps = n_reps)
  df <- read.csv(SAT_CSV, stringsAsFactors = FALSE)
  df <- df[df$d %in% c(2, 3) & df$status == "ok", ]
  for (meth in unique(df$method)) for (d in c(2, 3)) for (setting in unique(df$setting)) {
    sub <- df[df$method == meth & df$d == d & df$setting == setting, ]
    if (!nrow(sub)) next
    agg <- aggregate(BA ~ rep, sub, function(x) max(x) - min(x))
    cat(sprintf("  %-8s %-9s d=%d: mean BA-range over alpha = %.4f (n_reps=%d), max=%.4f\n",
                setting, meth, d, mean(agg$BA), nrow(agg), max(agg$BA)))
  }
}

# ---------------------------------------------------------------------------
# Summary / join
# ---------------------------------------------------------------------------
build_summary <- function() {
  tbl <- read.csv(TABLE_CSV, stringsAsFactors = FALSE)
  sat <- read.csv(SAT_CSV, stringsAsFactors = FALSE,
                  colClasses = c(flagged_idx = "character"))
  sat$flagged_idx[is.na(sat$flagged_idx)] <- ""
  ok <- sat[sat$status == "ok", ]

  rows <- list()
  for (setting in unique(ok$setting)) for (meth in unique(ok$method)) for (d in unique(ok$d[ok$setting == setting & ok$method == meth])) {
    sub <- ok[ok$setting == setting & ok$method == meth & ok$d == d, ]
    reps <- unique(sub$rep)
    ba_ranges <- numeric(0); min_jacs <- numeric(0); all_ident <- logical(0)
    for (r in reps) {
      s2 <- sub[sub$rep == r, ]
      s2 <- s2[order(-s2$alpha), ]  # 0.10 -> 0.001
      if (nrow(s2) < 2) next
      ba_ranges <- c(ba_ranges, max(s2$BA) - min(s2$BA))
      sets <- lapply(s2$flagged_idx, parse_idx)
      nlev <- length(sets); ident <- TRUE; mj <- 1
      for (i in 1:(nlev - 1)) for (j in (i + 1):nlev) {
        a <- sets[[i]]; b <- sets[[j]]
        if (!identical(a, b)) ident <- FALSE
        u <- length(union(a, b)); mj <- min(mj, if (u == 0) 1 else length(intersect(a, b)) / u)
      }
      min_jacs <- c(min_jacs, mj); all_ident <- c(all_ident, ident)
    }
    if (!length(ba_ranges)) next
    rows[[length(rows) + 1]] <- data.frame(
      setting = setting, method = meth, d = d, n_reps = length(ba_ranges),
      BA_range_mean = mean(ba_ranges), BA_range_sd = sd(ba_ranges),
      min_jaccard_mean = mean(min_jacs), min_jaccard_min = min(min_jacs),
      frac_identical_all_alpha = mean(all_ident),
      stringsAsFactors = FALSE)
  }
  behav <- do.call(rbind, rows)
  behav <- behav[order(behav$setting, behav$method, behav$d), ]
  write.csv(behav, SAT_SUM_CSV, row.names = FALSE)
  cat(sprintf("Wrote %s (%d rows)\n", SAT_SUM_CSV, nrow(behav)))

  # Join with Part 1's frac_zero at alpha=0.001 (the finest level) per d.
  t1 <- tbl[abs(tbl$alpha - 0.001) < 1e-9, c("d", "frac_zero", "move_norm_by_med01")]
  merged <- merge(behav, t1, by = "d", all.x = TRUE)
  merged <- merged[order(merged$setting, merged$method, merged$d), ]
  cat("\n==== JOIN: behavioural alpha-effect vs. table-level zero fraction ====\n")
  print(merged[, c("setting", "method", "d", "BA_range_mean", "min_jaccard_mean",
                    "frac_identical_all_alpha", "frac_zero", "move_norm_by_med01")], row.names = FALSE)
  invisible(merged)
}

# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (any(args == "--part1")) {
  run_part1()
} else if (any(args == "--control")) {
  nrep <- if (length(args) >= 2) as.integer(args[2]) else 10L
  run_control(nrep)
} else if (any(args == "--part2")) {
  rest <- args[args != "--part2"]
  dims <- if (length(rest) >= 1 && nzchar(rest[1])) as.integer(strsplit(rest[1], ",")[[1]]) else c(2, 3, 5, 10, 20, 35, 50, 100)
  nrep <- if (length(rest) >= 2) as.integer(rest[2]) else 30L
  cores <- if (length(rest) >= 3) as.integer(rest[3]) else 5L
  if (cores > 1) run_part2_parallel(dims, n_reps = nrep, cores = cores) else run_part2(dims, n_reps = nrep)
} else if (any(args == "--summary")) {
  build_summary()
} else {
  cat("Usage: --part1 | --control [n_reps] | --part2 [d1,d2,...] [n_reps] [cores] | --summary\n")
}
