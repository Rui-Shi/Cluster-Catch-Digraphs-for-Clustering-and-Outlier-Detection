#!/usr/bin/env Rscript
# revision_experiments/tr1/86_wp8_outlier_types.R
#
# WP8 experiment 2 -- local / bridged / collective outliers (R3.7). See
# WP8_PROTOCOL.md for the pre-declared design; this file implements it.
#
# R3.7 argues the mutual-catch-graph rule (outliers lack connectivity to a
# cluster core) is tuned for GLOBAL outliers and may be restrictive for:
#   - local outliers embedded within a cluster's own extent,
#   - anomalies bridged to a cluster core through a chain of connecting
#     points,
#   - internally well-connected collective groups of anomalies near a
#     cluster boundary.
# Three generators, one per case, each replacing only the outlier placement
# in the same two-cluster base (n1=95, n2=94 at cont=0.05, same geometry
# constants as 85_wp8_boundary_fp.R / 55_wp2c_simulation_arm.R). n0 = 10
# throughout, so evaluate()/count_scores2 is safe (n0 > 0), unlike experiment
# 4 (88), which has no outliers at all.
#
# Settings: type in {local_interior, local_shell, bridge, collective} x d in
# {3, 10}, n = 200, n0 = 10. 8 settings, 100 reps, 9 methods (4 MCCD + 5
# baselines) per rep.
#
# Single row per cell -- has_result() on the main CSV gates skip/restart
# directly (no done-marker needed, unlike 85).
#
# 2026-09-05 (Opus verification WP8_VERIFICATION.md, FAIL -> fixed): the
# original `local` and `bridge` generators did not produce the named
# phenomenon (measured local outlier/regular NN ratio 0.61 at d=10; measured
# 5.2-6.2 of 10 `bridge` points landed INSIDE a realised cluster ball). See
# the dated notes appended to WP8_PROTOCOL.md for the replacement
# construction and the measured acceptance-check values.
#
# 2026-09-05, round 2 (Opus re-verification WP8_REVERIFICATION.md): `bridge`
# still collapsed to a compact, `collective`-like blob under the round-1 fix
# (clearance 0.45 with 3-D N(0,0.15) jitter eroded ~10% of the realised gap);
# rebuilt with clearance 0.05 and jitter projected onto the orthogonal
# complement of the cluster axis, so the axial span is analytic (jitter can
# only increase, never decrease, distance from either centre) and no redraw
# is needed. `local` (round-1 shell construction) is kept as `local_shell`;
# the ORIGINAL pre-rewrite interior construction (git 5ed690c) is restored
# as `local_interior` so both placements are reported side by side -- the
# interior arm's near-zero TPR at d=10 (an interior point with doubled NN
# distance is not realisable at this n; see WP8_PROTOCOL.md) is itself part
# of the R3.7 answer, not a defect to fix away. See the round-2 dated note in
# WP8_PROTOCOL.md.
#
# Usage:
#   Rscript 86_wp8_outlier_types.R --smoke
#     One rep of EVERY (type, d) setting -- all 8 canonical settings, all 9
#     default methods -- written to
#     results/tr1/wp8/smoke/86_wp8_outlier_types.csv.
#   Rscript 86_wp8_outlier_types.R --summarize
#     Means and SEs of TPR/TNR/BA/F2 per (type, d, method); writes
#     results/tr1/wp8/86_outlier_types_summary.csv.
#   Rscript 86_wp8_outlier_types.R [settings] [reps] [budget] [methods]
#     settings: comma list of setting_id (from the fixed canonical grid
#               type x d in {3,10}, e.g. "local_interior_d3,bridge_d10"), or
#               "ALL" (default). The canonical grid is fixed regardless of
#               this selection, so the seed for a given setting never
#               depends on which subset of settings a chunk happens to
#               cover.
#     reps:     target replicate count (default 100)
#     budget:   stop starting new reps after this many seconds (default 480)
#     methods:  comma list, default the 9-method WP8 default (see below)

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))

ALL_THRESHOLDS <- REAL_DATA_THRESHOLDS
DEFAULT_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
                      "LOF", "DBSCAN", "MST", "ODIN", "iForest")
MCCD_METHODS <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD")
S_MIN <- 0.05

BASE_SEED <- 8601L
CONT <- 0.05
N_NOMINAL <- 200L
CLS_DIS <- 3; R_MIN <- 0.7; R_MAX <- 1.3

RESDIR    <- here::here("revision_experiments/results/tr1/wp8")
MAIN_CSV  <- file.path(RESDIR, "86_outlier_types.csv")
CLUST_CSV <- file.path(RESDIR, "86_outlier_types_clusters.csv")
SUMM_CSV  <- file.path(RESDIR, "86_outlier_types_summary.csv")
SMOKE_DIR <- file.path(RESDIR, "smoke")
SMOKE_MAIN  <- file.path(SMOKE_DIR, "86_wp8_outlier_types.csv")
SMOKE_CLUST <- file.path(SMOKE_DIR, "86_wp8_outlier_types_clusters.csv")

# ---------------------------------------------------------------------------
# Base two-cluster regular points, shared by all three generators.
# `attr(pts, "scale")` records the realised jitter multiplier
# (runif(1, R_MIN, R_MAX)) so callers can compute geometry (cluster surface,
# standoff) relative to what was ACTUALLY drawn, not the nominal range --
# needed by `bridge` (clearance from each cluster's realised surface) and
# `collective` (standoff past cluster 1's realised surface).
# `core_frac` (default 1) restricts the regular draw to a ball of radius
# `core_frac * scale` while `scale` itself still records the FULL jitter --
# used only by `local`, whose regular points must occupy a dense CORE
# strictly inside the nominal support so its outliers (drawn in the outer
# shell) are genuinely isolated from them.
# ---------------------------------------------------------------------------
draw_cluster <- function(n_k, d, mu_k, core_frac = 1) {
  scale <- runif(1, R_MIN, R_MAX)
  pts <- rpoisball.unit(n_k, d) * (core_frac * scale) +
    matrix(rep(mu_k, n_k), ncol = d, byrow = TRUE)
  attr(pts, "scale") <- scale
  pts
}

# ---------------------------------------------------------------------------
# Three outlier-type generators (WP8_PROTOCOL.md experiment 2, as amended by
# the 2026-09-05 dated notes below the original declarations)
# ---------------------------------------------------------------------------

# local_shell -- an outlier that is deep inside a cluster's spatial extent
# but locally isolated. The host cluster's own regular points are confined
# to a dense CORE of radius `core_frac * scale` = 0.4*scale; each local
# outlier sits at `runif(1, 0.8, 1.0) * scale` (shell 0.8-1.0*scale) along a
# uniformly random direction from its host centre -- inside the nominal
# support (<= scale) but strictly outside the core, so no rejection sampling
# is needed to keep the outlier isolated. (Comment corrected 2026-09-05,
# round 2, WP8_REVERIFICATION.md: this had drifted to describe the FIRST
# attempt's core 0.5*scale / shell 0.75-1.0*scale, which measured ratio 1.89
# at d=10 -- below the required 2.0. The code below, and the values in this
# comment, are the SECOND, adopted attempt: core 0.4*scale / shell
# 0.8-1.0*scale, measured ratio 6.032 at d=3 and 2.540 at d=10 -- see the
# round-1 acceptance-check note in WP8_PROTOCOL.md.) Per
# WP8_REVERIFICATION.md's two-placement decision, this shell construction is
# reported alongside the restored original interior construction
# (`local_interior`, below) rather than replacing it.
gen_local_shell <- function(seed, n, d, cont) {
  mu1 <- rep(3, d); mu2 <- c(3 + CLS_DIS, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- draw_cluster(n1, d, mu1, core_frac = 0.4)
  data2 <- draw_cluster(n2, d, mu2, core_frac = 0.4)
  scale1 <- attr(data1, "scale"); scale2 <- attr(data2, "scale")
  n_k1 <- ceiling(n0 / 2)
  outlier <- do.call(rbind, lapply(seq_len(n0), function(i) {
    host <- if (i <= n_k1) 1L else 2L
    mu_k <- if (host == 1L) mu1 else mu2
    scale_k <- if (host == 1L) scale1 else scale2
    dirv <- rnorm(d); dirv <- dirv / sqrt(sum(dirv^2))
    r <- runif(1, 0.8, 1.0) * scale_k
    mu_k + dirv * r
  }))
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)), n_inside = -1L)
}

# local_interior -- the ORIGINAL pre-rewrite construction (recovered
# verbatim from git 5ed690c, 2026-09-05 round 2, per WP8_REVERIFICATION.md's
# two-placement decision): a fixed 0.4*scale interior offset from the host
# centre, with the host cluster's own regular points kept OUT of a 0.3-radius
# ball around each outlier via rejection sampling (up to 200 attempts, kept
# anyway past that and logged via `n_forced`). This is the construction
# WP8_VERIFICATION.md originally FAILED (measured outlier-NN/regular-NN
# ratio 1.28 at d=3, 0.61 at d=10) -- it is restored here not because the
# ratio problem was wrong, but because the failure mode IS the R3.7 answer:
# at d=10, n=95, an interior point with doubled NN distance is not
# realisable at all (median within-cluster NN 0.699 vs realised radius
# 1.009), so an "embedded, isolated" local outlier at that dimension is a
# geometric impossibility, not a construction bug, and every method's
# TPR = 0 there is expected to be reported as such, not engineered away.
draw_cluster_avoid <- function(n_k, d, mu_k, avoid = NULL) {
  scale <- runif(1, R_MIN, R_MAX)
  pts <- matrix(NA_real_, nrow = n_k, ncol = d)
  n_forced <- 0L
  for (i in seq_len(n_k)) {
    ok <- FALSE
    for (attempt in 1:200) {
      cand <- as.numeric(rpoisball.unit(1, d)) * scale + mu_k
      clash <- FALSE
      if (!is.null(avoid)) {
        for (a in avoid) if (sqrt(sum((cand - a$center)^2)) < a$radius) { clash <- TRUE; break }
      }
      if (!clash) { ok <- TRUE; break }
    }
    if (!ok) n_forced <- n_forced + 1L
    pts[i, ] <- cand
  }
  attr(pts, "n_forced") <- n_forced
  attr(pts, "scale") <- scale
  pts
}

gen_local_interior <- function(seed, n, d, cont) {
  mu1 <- rep(3, d); mu2 <- c(3 + CLS_DIS, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  n_k1 <- ceiling(n0 / 2)
  gaps <- vector("list", n0)
  for (i in seq_len(n0)) {
    host <- if (i <= n_k1) 1L else 2L
    mu_k <- if (host == 1L) mu1 else mu2
    g <- mu_k + as.numeric(rpoisball.unit(1, d)) * 0.4
    gaps[[i]] <- list(center = g, radius = 0.3, host = host)
  }
  avoid1 <- Filter(function(a) a$host == 1L, gaps)
  avoid2 <- Filter(function(a) a$host == 2L, gaps)
  data1 <- draw_cluster_avoid(n1, d, mu1, avoid1)
  data2 <- draw_cluster_avoid(n2, d, mu2, avoid2)
  outlier <- do.call(rbind, lapply(gaps, function(a) a$center))
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)), n_inside = -1L)
}

# bridge -- a chain of n0 points spanning ONLY the gap between the two
# clusters' REALISED surfaces (not the nominal centre-to-centre segment).
# ROUND 1 FIX (2026-09-05): with s1, s2 the realised jitter scales,
# t_lo = (s1+0.25)/CLS_DIS, t_hi = 1-(s2+0.25)/CLS_DIS place the chain's
# endpoints just past each cluster's own realised radius; t_i interpolates
# n0 points evenly inside (t_lo, t_hi). The earlier construction spanned the
# full centre-to-centre segment regardless of the realised cluster sizes, so
# 5-6 of 10 points per rep landed inside a cluster's own realised ball.
# ROUND 2 FIX (2026-09-05, WP8_REVERIFICATION.md): round 1's 0.45 clearance
# with isotropic 3-D N(0,0.15) jitter still ate ~10% of the realised gap
# (t_lo > t_hi in 37-41% of reps), collapsing the chain to a compact,
# `collective`-like blob with TPR = 0 everywhere. Fixed by (a) shrinking
# clearance back to 0.05 -- it only needs to clear the realised surface, not
# absorb any jitter -- and (b) restricting jitter to the PERPENDICULAR
# subspace: draw a raw d-dim N(0,1) vector, subtract its projection onto the
# unit cluster-axis vector, scale the remainder to sd 0.03 per coordinate.
# Because the axis component of the jitter is now exactly zero, a
# perpendicular displacement strictly INCREASES Euclidean distance from both
# mu1 and mu2 relative to the (already-clear) axial placement (Pythagoras:
# distance^2 = axial_distance^2 + perp_distance^2) -- n_inside = 0 is
# therefore analytic, not merely likely, and the bounded-redraw loop from
# round 1 is removed entirely (nothing left for it to catch). See the
# round-2 dated note in WP8_PROTOCOL.md for the 100-rep measured values.
gen_bridge <- function(seed, n, d, cont) {
  mu1 <- rep(3, d); mu2 <- c(3 + CLS_DIS, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- draw_cluster(n1, d, mu1)
  data2 <- draw_cluster(n2, d, mu2)
  s1 <- attr(data1, "scale"); s2 <- attr(data2, "scale")
  clearance <- 0.05
  axis <- (mu2 - mu1) / CLS_DIS  # unit vector along the cluster axis
  t_lo <- (s1 + clearance) / CLS_DIS
  t_hi <- 1 - (s2 + clearance) / CLS_DIS
  outlier <- do.call(rbind, lapply(seq_len(n0), function(i) {
    t_i <- t_lo + (i - 0.5) / n0 * (t_hi - t_lo)
    base <- mu1 + t_i * (mu2 - mu1)
    g <- rnorm(d)
    g_perp <- g - sum(g * axis) * axis  # strip the axial component
    base + 0.03 * g_perp
  }))
  d1 <- sqrt(rowSums(sweep(outlier, 2, mu1)^2))
  d2 <- sqrt(rowSums(sweep(outlier, 2, mu2)^2))
  n_inside <- sum(d1 < s1 | d2 < s2)
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)), n_inside = n_inside)
}

# collective -- a tight minority group just past cluster 1's REALISED
# surface. FIX (2026-09-05): standoff is `s1 + 0.5` with s1 the realised
# jitter scale (was the fixed nominal r_max = 1.3), so the group's distance
# from the boundary no longer depends on which cluster's jitter happened to
# be drawn; correct the protocol's "just past the cluster's typical edge" to
# this measured stand-off (0.5 past the REALISED radius, not the nominal
# max).
gen_collective <- function(seed, n, d, cont) {
  mu1 <- rep(3, d); mu2 <- c(3 + CLS_DIS, rep(3, d - 1))
  n1 <- round(n * (1 - cont) * 0.5); n2 <- round(n * (1 - cont) * 0.5) - 1
  n0 <- round(n * cont)
  set.seed(seed)
  data1 <- draw_cluster(n1, d, mu1)
  data2 <- draw_cluster(n2, d, mu2)
  s1 <- attr(data1, "scale")
  u <- rnorm(d); u <- u / sqrt(sum(u^2))
  center <- mu1 + u * (s1 + 0.5)
  outlier <- matrix(rep(center, n0), nrow = n0, byrow = TRUE) +
    matrix(rnorm(n0 * d, 0, 0.15), nrow = n0)
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0, n1 = n1, n2 = n2,
       true_cluster = c(rep(1L, n1), rep(2L, n2)), n_inside = -1L)
}

GENERATORS <- list(local_interior = gen_local_interior, local_shell = gen_local_shell,
                    bridge = gen_bridge, collective = gen_collective)

build_settings <- function(dims) {
  s <- do.call(rbind, lapply(names(GENERATORS), function(g)
    do.call(rbind, lapply(dims, function(d)
      data.frame(type = g, d = d, stringsAsFactors = FALSE)))))
  s$setting_id <- sprintf("%s_d%d", s$type, s$d)
  s
}

# Fixed canonical grid (WP8_REVERIFICATION.md, round 2: type x d in {3,10},
# type in {local_interior, local_shell, bridge, collective}). CANON never
# changes with how a run is chunked, which is what fixes the seed bug: s_i
# is always this table's row index, never the row index of a re-built,
# possibly-DIMS-restricted table.
FULL_DIMS <- c(3L, 10L)
CANON <- build_settings(FULL_DIMS)

resolve_settings <- function(sel) {
  if (is.null(sel) || (length(sel) == 1 && toupper(sel) == "ALL")) return(CANON)
  out <- CANON[match(sel, CANON$setting_id), ]
  if (anyNA(out$setting_id)) stop(sprintf("unknown setting_id(s): %s",
                                          paste(sel[is.na(match(sel, CANON$setting_id))], collapse = ", ")))
  out
}

# ---------------------------------------------------------------------------
# One (setting, d, rep, method) cell
# ---------------------------------------------------------------------------
run_cell <- function(setting_id, type, d, rep_id, seed, method, main_csv, clust_csv) {
  keys <- c(setting_id = setting_id, d = d, rep = rep_id, method = method)
  if (isTRUE(has_result(main_csv, keys))) {
    cat(sprintf("[skip] %s d=%d rep=%d method=%s\n", setting_id, d, rep_id, method))
    return(invisible("skip"))
  }

  dat <- GENERATORS[[type]](seed, N_NOMINAL, d, CONT)
  X <- dat$X; n <- dat$n; n0 <- dat$n0; n_reg <- dat$n1 + dat$n2
  Y <- c(rep(1L, n_reg), rep(0L, n0))

  extra <- if (method %in% c("SU-MCCD", "SUN-MCCD")) list(min.cls = S_MIN) else list()
  out <- tryCatch({
    res <- do.call(METHOD_REGISTRY[[method]], c(list(X = X, d = d, Y = Y), extra))
    stopifnot(length(res$score) == n, !anyNA(res$score))
    m <- evaluate(Y, res$score, ALL_THRESHOLDS[[method]])
    # note = "-" (never "" or NA) on success. has_result() treats ANY NA
    # in a non-key column as an incomplete/partial row (its own
    # truncation guard) -- and an empty string doesn't dodge this: a CSV
    # column that is "" on EVERY row round-trips through read.csv's
    # automatic type.convert() as logical NA (an all-empty character
    # column has no type to infer, so it collapses to `logical NA`),
    # which is exactly the NA has_result is watching for. A non-empty
    # placeholder avoids the type.convert collapse. The same reasoning is
    # why `n_inside` uses the sentinel -1L (not NA_integer_) for the two
    # types (local, collective) where it is not applicable -- see the NA
    # convention note in WP8_PROTOCOL.md.
    list(m = m, res = res, status = "ok", note = "-")
  }, error = function(e) list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
                              res = NULL, status = "error",
                              note = substr(conditionMessage(e), 1, 200)))

  row <- list(setting_id = setting_id, type = type, d = d, rep = rep_id, seed = seed,
              method = method,
              TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
              BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
              n_flagged = if (is.null(out$res)) NA_integer_ else sum(out$res$score >= ALL_THRESHOLDS[[method]]),
              n_inside = as.integer(dat$n_inside),
              status = out$status, note = out$note)
  append_result(main_csv, row)

  if (identical(out$status, "ok") && method %in% MCCD_METHODS && !is.null(out$res$cluster)) {
    # Block write: all n_reg cluster rows for this cell in ONE append_result
    # call (append_result accepts a list of equal-length vectors and
    # constructs a multi-row data.frame from it) -- 6.4 ms per 189-row block
    # vs 773 ms row by row (WP8_VERIFICATION.md).
    idx <- seq_len(n_reg)
    dc <- out$res$cluster[idx]
    append_result(clust_csv, list(
      setting_id = rep(setting_id, n_reg), d = rep(d, n_reg), rep = rep(rep_id, n_reg),
      seed = rep(seed, n_reg), method = rep(method, n_reg), row_index = idx,
      true_cluster = dat$true_cluster[idx],
      detected_cluster = ifelse(is.na(dc), NA_integer_, dc)))
  }

  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE %s d=%d rep=%d method=%s TPR=%.3f TNR=%.3f\n",
                setting_id, d, rep_id, method, out$m[["TPR"]], out$m[["TNR"]]))
  } else {
    cat(sprintf("CELL_FAIL %s d=%d rep=%d method=%s note=%s\n",
                setting_id, d, rep_id, method, out$note))
  }
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# --summarize: means and SEs over reps, status=="ok" only, de-duplicated on
# the key (setting_id, d, method) keeping the LAST occurrence -- a retried
# cell leaves an error row and a later ok row under the same key.
# ---------------------------------------------------------------------------
do_summarize <- function() {
  if (!file.exists(MAIN_CSV)) { cat("no results yet\n"); return(invisible(NULL)) }
  df_all <- read.csv(MAIN_CSV, stringsAsFactors = FALSE)
  key <- paste(df_all$setting_id, df_all$d, df_all$rep, df_all$method, sep = "\r")
  df_all <- df_all[!duplicated(key, fromLast = TRUE), ]
  df <- df_all[df_all$status == "ok", ]

  rows <- list()
  for (sid in unique(df_all$setting_id)) for (meth in unique(df$method[df$setting_id == sid])) {
    sub <- df[df$setting_id == sid & df$method == meth, ]
    if (!nrow(sub)) next
    n_err <- sum(df_all$setting_id == sid & df_all$method == meth & df_all$status != "ok")
    rows[[length(rows) + 1L]] <- data.frame(
      setting_id = sid, type = sub$type[1], d = sub$d[1], method = meth, n_reps = nrow(sub),
      TPR_mean = mean(sub$TPR), TPR_se = sd(sub$TPR) / sqrt(nrow(sub)),
      TNR_mean = mean(sub$TNR), TNR_se = sd(sub$TNR) / sqrt(nrow(sub)),
      BA_mean  = mean(sub$BA),  BA_se  = sd(sub$BA)  / sqrt(nrow(sub)),
      F2_mean  = mean(sub$F2),  F2_se  = sd(sub$F2)  / sqrt(nrow(sub)),
      n_flagged_mean = mean(sub$n_flagged),
      n_inside_mean = mean(sub$n_inside[sub$n_inside >= 0]),
      n_error = n_err, stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  S <- S[order(S$type, S$d, S$method), ]
  write.csv(S, SUMM_CSV, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMM_CSV, nrow(S)))
  print(S, row.names = FALSE)
  invisible(S)
}

# ---------------------------------------------------------------------------
# --run: full (subset) grid, checkpointed by budget
# ---------------------------------------------------------------------------
do_run <- function(sel, n_reps, budget, methods) {
  st <- resolve_settings(sel)
  cat(sprintf("86_wp8_outlier_types: %d settings x %d reps x %d methods (budget %.0fs)\n",
              nrow(st), n_reps, length(methods), budget))
  t_start <- Sys.time()
  for (row in seq_len(nrow(st))) {
    s <- st[row, ]
    s_i <- match(s$setting_id, CANON$setting_id)
    for (rep_id in seq_len(n_reps)) {
      if (as.numeric(difftime(Sys.time(), t_start, units = "secs")) > budget) {
        cat(sprintf("[budget %gs exhausted -- rerun to continue] next up: %s rep %d\n",
                    budget, s$setting_id, rep_id))
        return(invisible("budget_stop"))
      }
      seed <- BASE_SEED + 100000L * s_i + rep_id
      for (m in methods) run_cell(s$setting_id, s$type, s$d, rep_id, seed, m, MAIN_CSV, CLUST_CSV)
    }
  }
  cat("86_wp8_outlier_types: run complete.\n")
}

# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_SMOKE <- "--smoke" %in% args
MODE_SUMM  <- "--summarize" %in% args
args <- args[!args %in% c("--smoke", "--summarize")]

if (!isTRUE(getOption("wp8.no_main", FALSE))) {
  if (MODE_SMOKE) {
    dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("==== 86_wp8_outlier_types --smoke: all %d settings (both d, all 4 types), 1 rep, %d methods ====\n",
                nrow(CANON), length(DEFAULT_METHODS)))
    for (row in seq_len(nrow(CANON))) {
      s <- CANON[row, ]
      seed <- BASE_SEED + 100000L * row + 9001L
      for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$type, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
    }
    cat("---- re-invoking the same cells to demonstrate resume/skip ----\n")
    for (row in seq_len(nrow(CANON))) {
      s <- CANON[row, ]
      seed <- BASE_SEED + 100000L * row + 9001L
      for (m in DEFAULT_METHODS) run_cell(s$setting_id, s$type, s$d, 9001L, seed, m, SMOKE_MAIN, SMOKE_CLUST)
    }
    cat("ALL_CELLS_COMPLETE\n")
    quit(save = "no", status = 0)
  } else if (MODE_SUMM) {
    do_summarize()
    quit(save = "no", status = 0)
  } else {
    sel     <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else NULL
    n_reps  <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 100L
    budget  <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 480
    methods <- if (length(args) >= 4 && nzchar(args[4])) strsplit(args[4], ",")[[1]] else DEFAULT_METHODS
    do_run(sel, n_reps, budget, methods)
  }
}
