#!/usr/bin/env Rscript
# revision_experiments/tr3/33_type1.R
#
# Type I error of the UN-CCD radius search (KBS review: Reviewer 2 says the
# MC-SRT, applied at many radii with no multiplicity correction, inflates false
# positives; Reviewer 3.1 asks for finite-sample Type I behaviour). The paper
# says the size of the test at one radius lies between alpha and 2 * alpha.
#
# The test, as nnccd.radi.ascend.one() / nnccd.radi.descend.one() apply it
# (R/ccds/UN_CCD.R): for a ball with centre i and m members (the points of the
# ball other than the centre; r = distance from the centre to the farthest
# member), take the nearest-neighbour distances of the members among the
# members, and reject complete spatial randomness (CSR) if
#     mean(nnd) / r < q_mean[m]   OR   median(nnd) / r < q_median[m],
# q_*[m] being entry m of the NN table (simul$average, simul$median). The
# descending search accepts a ball when both statistics are above their entries,
# which is the same rule negated.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/33_type1.R <part>
#     test     size of the test at ONE radius (one row per d and m)
#     ccd      the whole procedure on n uniform points, no clusters present
#              (one row per replicate; one done row per cell d, n)
#     summary  type1_test_summary.csv and type1_ccd_summary.csv from the two above
#     check    (extra, not part of the experiment) replays the radius search of
#              the method with this script's own statistic and table index and
#              compares the radii; also compares t1_one() with ccd_cluster()
#   TR3_CORES overrides the number of workers (default detectCores()).
#   Reruns skip finished cells (done files). Move a result file aside to recompute.
#
# SETTINGS (environment variables; the defaults are the full experiment)
#   TR3_OUT_DIR        output folder (default TR3_RES of tr3_lib.R)
#   TR3_TYPE1_B        part test: independent CSR samples per (d, m), default 10000
#                      (at most 20000)
#   TR3_TYPE1_D        part test: dimensions, default 2,3,5,10,20
#   TR3_TYPE1_M        part test: numbers of members, default 5,10,20,50,100,200,500
#   TR3_TYPE1_CCD_D    part ccd: dimensions, default 2,3,5,10,20
#   TR3_TYPE1_CCD_N    part ccd: sample sizes, default 50,100,200,500
#   TR3_TYPE1_REPS     part ccd: replicates per cell, default 500
#
# TABLES AND LEVELS
#   For each d the table and the search direction are read from the UN-CCDs rows
#   of settings_simulation.csv (one table per d across all rows; the script
#   stops otherwise): d = 2 NN-test-simul_2d_85%, 3 _3d_90%, 5 _5d_95%,
#   10 _10d_99%, 20 _20d_999% (read on 2026-10-25). alpha = 1 - level.
#
# PART test
#   B independent samples under CSR "exactly as the table generator does"
#   (NNDestP.simpois.lower.quant / .ccd_qt_nn_once): m points iid uniform in the
#   unit ball, drawn by a copy of rpoisball.unit(). The statistics are the
#   clustering's with r = 1: the generator's ball has radius 1, so the
#   clustering's division by r is a division by 1. Entry m of the table is the
#   lower (1 - level) quantile of the generator's statistic for m points, which
#   is the clustering's entry for a ball with m members (env.*[j - 1] with
#   j - 1 members in the ascending search, env.*.rev[j + 1] = entry n - j with
#   n - j members in the descending one). Part check verifies this by replay.
#   SUPPLEMENTARY null (columns ball_*): in the clustering the farthest member
#   lies ON the boundary of the ball (r is defined by it). For a Poisson process,
#   given r, the other m - 1 members are iid uniform in the ball of radius r and
#   the farthest is uniform on its boundary. The columns ball_* repeat the test
#   with m - 1 points iid uniform in the unit ball plus one point uniform on the
#   unit sphere. The generator does not condition on this, so the two nulls can
#   differ; the columns are there to show by how much. rej_* is the size of the test
#   as tabulated; ball_* is the null the clustering applies when the ball lies inside
#   the window (verifier, 2026-10-25). Both are reported.
#   Random numbers: RNGkind("L'Ecuyer-CMRG"), set.seed(TEST_SEED) on the main
#   process, then one stream per chunk of 100 draws via parallel::nextRNGStream,
#   indexed by the position of (d, m, chunk) in the FULL default grids; each task
#   assigns its stream to .Random.seed. The results therefore depend on neither
#   the number of workers nor the scheduling, nor on the subset of the grid or
#   on B (B = 200 is the first two chunks of B = 10000). This replaces
#   clusterSetRNGStream, which ties streams to workers and, with dynamic
#   scheduling, to the order in which tasks finish.
#   TEST_SEED = 3310001. The table generator uses 20260927 + d (quantile_table.R)
#   and the shipped tables came from NNDestP.simpois.lower.quant, which sets no
#   seed; the simulation drivers use 2024.
#
# PART ccd
#   UN-CCDs with the published settings of the UNIFORM rows of
#   settings_simulation.csv (table and direction for d and n; low.num 3, niter
#   1000, greedy2: all inside UNCCD_clustering()) on n points iid uniform in
#   [0, 1]^d, REPS replicates per cell. Data: set.seed(seed) on the main
#   process, then REPS matrices matrix(runif(n * d), n, d); seed =
#   33000000 + 1000 * d + n (distinct for every cell; the first r replicates of a
#   cell do not depend on REPS). The clustering is deterministic.
#   The method is called as ccd_cluster() calls it, UNCCD_clustering(datax,
#   simul, method = direction), but directly, to read Int.D from the result.
#   Per replicate: k_hat = si.ind (what the simulations average), si, n_dom2 =
#   length(Int.D) (the second dominating set; k_hat can only be 1 or at most
#   n_dom2), n_dom1 = length(D) (D is dominate.mat.greedy2()'s ordering of all n
#   points, so n_dom1 = n always; kept as a check), seconds (elapsed time
#   of the call), data_fp (fingerprint of the data, tr3_fp()).
#   k_hat = 1 is what nnccd.silhouette() returns when the second dominating set
#   has fewer than 2 members or no k >= 2 has a positive silhouette (with more
#   than one point in the last cluster), so it is reachable.
#
# OUTPUT (OUT_DIR)
#   type1_test.csv, type1_test_done.csv
#   type1_ccd_replicates.csv, type1_ccd_done.csv
#   type1_test_summary.csv, type1_ccd_summary.csv
# Numbers are written with fmt17(); a cell is appended as soon as it finishes.

args <- commandArgs(trailingOnly = TRUE)
PARTS <- c("test", "ccd", "summary", "check")
if (length(args) < 1 || !(args[1] %in% PARTS)) {
  stop("usage: Rscript revision_experiments/tr3/33_type1.R <test|ccd|summary|check>")
}
PART <- args[1]
source(here::here("revision_experiments/tr3/tr3_lib.R"))

OUT_DIR <- Sys.getenv("TR3_OUT_DIR", TR3_RES)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- function(name) file.path(OUT_DIR, name)
TEST_CSV <- out_path("type1_test.csv");           TEST_DONE <- out_path("type1_test_done.csv")
REP_CSV  <- out_path("type1_ccd_replicates.csv"); CCD_DONE  <- out_path("type1_ccd_done.csv")
TEST_SUM <- out_path("type1_test_summary.csv");   CCD_SUM   <- out_path("type1_ccd_summary.csv")

env_int <- function(name, default) {
  v <- Sys.getenv(name, "")
  as.integer(if (nzchar(v)) strsplit(v, ",")[[1]] else default)
}
FULL_D <- c(2L, 3L, 5L, 10L, 20L)
FULL_M <- c(5L, 10L, 20L, 50L, 100L, 200L, 500L)
B        <- env_int("TR3_TYPE1_B", 10000L)
T_D      <- env_int("TR3_TYPE1_D", FULL_D)
T_M      <- env_int("TR3_TYPE1_M", FULL_M)
CCD_D    <- env_int("TR3_TYPE1_CCD_D", FULL_D)
CCD_N    <- env_int("TR3_TYPE1_CCD_N", c(50L, 100L, 200L, 500L))
REPS     <- env_int("TR3_TYPE1_REPS", 500L)
TEST_SEED <- 3310001L
CHUNK    <- 100L        # draws per random-number stream
MAX_B    <- 20000L

# ---------------------------------------------------------------------------
# Tables and levels from settings_simulation.csv
# ---------------------------------------------------------------------------

#' Table file, search direction, level and alpha of UN-CCDs at dimension d,
#' over all its rows of settings_simulation.csv (both distributions, all n).
t1_table_for <- function(d) {
  s <- tr3_settings("simulation")
  u <- s[s$method == "UN-CCDs" & s$d == as.character(d), ]
  if (!nrow(u)) stop("t1_table_for(): no UN-CCDs rows for d = ", d)
  if (length(unique(u$table)) != 1 || length(unique(u$direction)) != 1) {
    stop("t1_table_for(): d = ", d, " has several tables or directions: ",
         paste(unique(paste(u$table, u$direction)), collapse = "; "))
  }
  lab <- sub("^NN-test-simul_[0-9]+d_([0-9]+)%\\.RData$", "\\1", u$table[1])
  if (identical(lab, u$table[1])) stop("t1_table_for(): cannot read the level of ", u$table[1])
  level <- as.numeric(paste0("0.", lab))
  list(file = u$table[1], direction = u$direction[1], level = level, alpha = round(1 - level, 10))
}

# ---------------------------------------------------------------------------
# The draw and the statistics (global functions: foreach exports them)
# ---------------------------------------------------------------------------

# rpoisball.unit() of R/ccds/NN_Dist_Est.R, as .ccd_qt_nn_once() copies it:
# n points iid uniform in the unit ball of R^d, returned as an n x d matrix.
t1_rpoisball <- function(n, d) {
  r1 <- runif(n, 0, 1)^(1/d)
  norm.data <- matrix(MASS::mvrnorm(n, rep(0, d), diag(d)), ncol = d, byrow = T)
  data1 <- apply(norm.data, 1, function(x) x/sqrt(sum(x^2)))
  data1 <- apply(data1, 1, function(x) x*r1)
  return(data1)
}

# one point uniform on the unit sphere, as a 1 x d matrix
t1_sphere <- function(d) {
  z <- matrix(MASS::mvrnorm(1, rep(0, d), diag(d)), ncol = d, byrow = T)
  z / sqrt(sum(z^2))
}

# (mean, median) of the nearest-neighbour distances among the rows of a
# distance matrix, divided by r (the clustering's statistic; NNDest.dist.f())
t1_stat_D <- function(D, r = 1) {
  diag(D) <- Inf
  nn <- D[cbind(seq_len(nrow(D)), max.col(-D, ties.method = "first"))]
  c(mean(nn) / r, median(nn) / r)
}
t1_stat <- function(X, r = 1) t1_stat_D(as.matrix(dist(X)), r)

# one chunk of nb draws for one (d, m): counts of rejections under the
# generator's null (4 numbers) and the supplementary null (4 numbers)
t1_chunk <- function(stream, nb, d, m, q_mean, q_med) {
  assign(".Random.seed", stream, envir = globalenv())
  cnt <- integer(8)
  for (b in seq_len(nb)) {
    st <- t1_stat(t1_rpoisball(m, d))
    a <- st[1] < q_mean; z <- st[2] < q_med
    cnt[1:4] <- cnt[1:4] + c(a, z, a & z, a | z)
    st <- t1_stat(rbind(t1_rpoisball(m - 1, d), t1_sphere(d)))
    a <- st[1] < q_mean; z <- st[2] < q_med
    cnt[5:8] <- cnt[5:8] + c(a, z, a & z, a | z)
  }
  cnt
}

# one replicate of part ccd: scalars only
t1_one <- function(x, par) {
  t0 <- proc.time()[["elapsed"]]
  res <- UNCCD_clustering(datax = x, simul = par$simul, method = par$direction)
  secs <- proc.time()[["elapsed"]] - t0
  if (length(res$si) != 1 || length(res$si.ind) != 1) stop("t1_one(): si or si.ind is not a scalar")
  c(k_hat = res$si.ind, si = res$si, n_dom2 = length(res$Int.D), n_dom1 = length(res$D),
    seconds = secs)
}

binom_se <- function(p, n) sqrt(p * (1 - p) / n)

# ---------------------------------------------------------------------------
# Part test
# ---------------------------------------------------------------------------
if (PART == "test") {
  if (B > MAX_B || B < 1) stop("33_type1.R: TR3_TYPE1_B must be between 1 and ", MAX_B)
  if (!all(T_D %in% FULL_D) || !all(T_M %in% FULL_M)) {
    stop("33_type1.R: d and m must come from the default grids (the random-number streams are indexed by them)")
  }
  source(here::here("R/ccds/quantile_table.R"))   # load_quantile_table() for tr3_load_table(); functions only
  tr3_library(c("parallel", "doParallel", "MASS"))
  nch_max <- MAX_B %/% CHUNK
  nch <- ceiling(B / CHUNK)
  sizes <- c(rep(CHUNK, nch - 1), B - CHUNK * (nch - 1))

  # one stream per (d, m, chunk) of the full grids, in a fixed order
  RNGkind("L'Ecuyer-CMRG")
  set.seed(TEST_SEED)
  total <- length(FULL_D) * length(FULL_M) * nch_max
  streams <- vector("list", total)
  sd0 <- .Random.seed
  for (i in seq_len(total)) { streams[[i]] <- sd0; sd0 <- parallel::nextRNGStream(sd0) }

  cores <- tr3_cores()
  cat(sprintf("33_type1.R test: B=%d d=%s m=%s cores=%d seed=%d out=%s\n", B,
              paste(T_D, collapse = ","), paste(T_M, collapse = ","), cores, TEST_SEED, OUT_DIR))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  t_start <- Sys.time()

  for (d in T_D) {
    tb <- t1_table_for(d)
    simul <- tr3_load_table(tb$file)
    for (m in T_M) {
      key <- list(d = d, m = m, B = B)
      if (has_result(TEST_DONE, key)) {
        cat(sprintf("  d=%d m=%d: done, skipped\n", d, m))
        next
      }
      if (length(simul$average) < m || length(simul$median) < m) stop("table ", tb$file, " is shorter than m = ", m)
      q_mean <- simul$average[m]
      q_med <- simul$median[m]
      base <- ((match(d, FULL_D) - 1) * length(FULL_M) + (match(m, FULL_M) - 1)) * nch_max
      streams_cell <- streams[base + seq_len(nch)]
      t0 <- Sys.time()
      r <- foreach(stream = streams_cell, nb = sizes, .packages = "MASS") %dopar%
        t1_chunk(stream, nb, d, m, q_mean, q_med)
      cnt <- colSums(do.call(rbind, r))
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      p <- cnt / B
      rows <- data.frame(d = d, m = m, B = B, table = tb$file, direction = tb$direction,
                         level = fmt17(tb$level), alpha = fmt17(tb$alpha), seed = TEST_SEED,
                         q_mean = fmt17(q_mean), q_median = fmt17(q_med),
                         n_mean = cnt[1], n_median = cnt[2], n_both = cnt[3], n_rule = cnt[4],
                         rej_mean = fmt17(p[1]), se_mean = fmt17(binom_se(p[1], B)),
                         rej_median = fmt17(p[2]), se_median = fmt17(binom_se(p[2], B)),
                         rej_rule = fmt17(p[4]), se_rule = fmt17(binom_se(p[4], B)),
                         ball_n_mean = cnt[5], ball_n_median = cnt[6], ball_n_both = cnt[7],
                         ball_n_rule = cnt[8],
                         ball_rej_mean = fmt17(p[5]), ball_se_mean = fmt17(binom_se(p[5], B)),
                         ball_rej_median = fmt17(p[6]), ball_se_median = fmt17(binom_se(p[6], B)),
                         ball_rej_rule = fmt17(p[8]), ball_se_rule = fmt17(binom_se(p[8], B)),
                         stringsAsFactors = FALSE)
      tr3_append(TEST_CSV, rows)
      tr3_append(TEST_DONE, data.frame(key, seconds = fmt17(secs), cores = cores,
                                       finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                       stringsAsFactors = FALSE))
      cat(sprintf("  d=%2d m=%3d alpha=%.3f  mean %.4f  median %.4f  rule %.4f (se %.4f)  | ball: rule %.4f  [%.1f s]\n",
                  d, m, tb$alpha, p[1], p[2], p[4], binom_se(p[4], B), p[8], secs))
    }
  }
  print(Sys.time() - t_start)
}

# ---------------------------------------------------------------------------
# Part ccd
# ---------------------------------------------------------------------------
if (PART == "ccd") {
  tr3_source_method("UN-CCDs", "simulation")
  cores <- tr3_cores()
  cat(sprintf("33_type1.R ccd: d=%s n=%s reps=%d cores=%d out=%s\n", paste(CCD_D, collapse = ","),
              paste(CCD_N, collapse = ","), REPS, cores, OUT_DIR))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  t_start <- Sys.time()

  # foreach is called at top level, as in 10_sim.R (it exports only what the
  # calling environment holds, following global functions into their bodies)
  for (d in CCD_D) {
    for (n in CCD_N) {
      key <- list(d = d, n = n, reps = REPS)
      if (has_result(CCD_DONE, key)) {
        cat(sprintf("  d=%d n=%d: done, skipped\n", d, n))
        next
      }
      s <- tr3_sim_setting("UN-CCDs", "uniform", d, n)
      p <- list(simul = tr3_load_table(s$table), direction = s$direction)
      seed <- 33000000L + 1000L * d + n
      set.seed(seed)
      data_list <- lapply(seq_len(REPS), function(i) matrix(runif(n * d), nrow = n, ncol = d))
      t0 <- Sys.time()
      r <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% t1_one(x, p)
      r <- do.call(rbind, r)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rows <- data.frame(key, rep = seq_len(nrow(r)), seed = seed, table = s$table,
                         direction = s$direction,
                         k_hat = fmt17(r[, "k_hat"]), si = fmt17(r[, "si"]),
                         n_dom2 = fmt17(r[, "n_dom2"]), n_dom1 = fmt17(r[, "n_dom1"]),
                         seconds = fmt17(r[, "seconds"]),
                         data_fp = fmt17(vapply(data_list, tr3_fp, numeric(1))),
                         stringsAsFactors = FALSE)
      tr3_append(REP_CSV, rows)
      tr3_append(CCD_DONE, data.frame(key, n_rows = nrow(rows), seed = seed, cores = cores,
                                      wall_secs = fmt17(secs),
                                      finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      stringsAsFactors = FALSE))
      cat(sprintf("  d=%2d n=%3d (%s, %s): P(k_hat = 1) = %.3f, mean k_hat = %.3f, mean n_dom2 = %.2f, mean %.2f s/rep, wall %.1f s\n",
                  d, n, s$table, s$direction, mean(r[, "k_hat"] == 1), mean(r[, "k_hat"]),
                  mean(r[, "n_dom2"]), mean(r[, "seconds"]), secs))
    }
  }
  print(Sys.time() - t_start)
}

# ---------------------------------------------------------------------------
# Part summary
# ---------------------------------------------------------------------------
if (PART == "summary") {
  num <- function(df, cols) { for (k in cols) df[[k]] <- as.numeric(df[[k]]); df }

  tt <- tr3_read_done(TEST_CSV, TEST_DONE, keys = c("d", "m", "B"))
  if (is.null(tt)) stop("33_type1.R summary: no finished cells in ", TEST_CSV)
  tt <- num(tt, setdiff(names(tt), c("table", "direction")))
  tt <- tt[order(tt$d, tt$m), ]
  z <- function(p, a, B) (p - a) / sqrt(a * (1 - a) / B)
  ts <- data.frame(d = tt$d, m = tt$m, B = tt$B, alpha = tt$alpha,
                   rej_mean = tt$rej_mean, se_mean = tt$se_mean,
                   rej_median = tt$rej_median, se_median = tt$se_median,
                   rej_rule = tt$rej_rule, se_rule = tt$se_rule,
                   rule_over_alpha = tt$rej_rule / tt$alpha,
                   z_mean_vs_alpha = z(tt$rej_mean, tt$alpha, tt$B),
                   z_median_vs_alpha = z(tt$rej_median, tt$alpha, tt$B),
                   ball_rej_mean = tt$ball_rej_mean, ball_rej_median = tt$ball_rej_median,
                   ball_rej_rule = tt$ball_rej_rule, ball_se_rule = tt$ball_se_rule,
                   ball_rule_over_alpha = tt$ball_rej_rule / tt$alpha)
  out <- ts
  for (k in names(out)) if (!(k %in% c("d", "m", "B"))) out[[k]] <- fmt17(out[[k]])
  write.csv(out, TEST_SUM, row.names = FALSE, quote = FALSE)
  cat("wrote", TEST_SUM, "(", nrow(ts), "rows )\n")
  print(data.frame(d = ts$d, m = ts$m, alpha = ts$alpha, mean = round(ts$rej_mean, 4),
                   median = round(ts$rej_median, 4), rule = round(ts$rej_rule, 4),
                   se_rule = round(ts$se_rule, 4), ratio = round(ts$rule_over_alpha, 2),
                   ball_rule = round(ts$ball_rej_rule, 4)), row.names = FALSE)

  cc <- tr3_read_done(REP_CSV, CCD_DONE, keys = c("d", "n", "reps"), row_key = "rep")
  if (is.null(cc)) stop("33_type1.R summary: no finished cells in ", REP_CSV)
  cc <- num(cc, c("d", "n", "reps", "k_hat", "si", "n_dom2", "n_dom1", "seconds"))
  cells <- unique(cc[, c("d", "n", "reps")])
  cells <- cells[order(cells$d, cells$n), ]
  cs <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    g <- cc[cc$d == cells$d[i] & cc$n == cells$n[i] & cc$reps == cells$reps[i], ]
    p1 <- mean(g$k_hat == 1)
    data.frame(d = cells$d[i], n = cells$n[i], reps = cells$reps[i], n_reps = nrow(g),
               prop_k1 = p1, se_prop_k1 = binom_se(p1, nrow(g)), mean_k_hat = mean(g$k_hat),
               mean_n_dom2 = mean(g$n_dom2), mean_n_dom1 = mean(g$n_dom1),
               mean_si = mean(g$si), mean_seconds = mean(g$seconds))
  }))
  out <- cs
  for (k in setdiff(names(out), c("d", "n", "reps", "n_reps"))) out[[k]] <- fmt17(out[[k]])
  write.csv(out, CCD_SUM, row.names = FALSE, quote = FALSE)
  cat("wrote", CCD_SUM, "(", nrow(cs), "rows )\n")
  print(data.frame(d = cs$d, n = cs$n, reps = cs$n_reps, prop_k1 = round(cs$prop_k1, 3),
                   se = round(cs$se_prop_k1, 3), mean_k = round(cs$mean_k_hat, 3),
                   n_dom2 = round(cs$mean_n_dom2, 2), secs = round(cs$mean_seconds, 3)),
        row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Part check (extra): is this script's statistic and table index the clustering's?
# ---------------------------------------------------------------------------
# Replays the radius search of nnccd.radi() for every centre of small uniform
# data sets with t1_stat_D() and the table entry indexed by the NUMBER OF
# MEMBERS (length of the member set), and compares the radii with the method's.
# Ascending: first ball with a rejection; the radius is that of the ball before.
# Descending: first ball (from the largest) with no rejection; its radius.
if (PART == "check") {
  tr3_source_method("UN-CCDs", "simulation")
  replay_radius <- function(ddx, i, simul, direction, low.num = 3) {
    n <- nrow(ddx)
    rejects <- function(mem, r) {
      st <- t1_stat_D(ddx[mem, mem, drop = FALSE], r)
      st[1] < simul$average[length(mem)] | st[2] < simul$median[length(mem)]
    }
    if (direction == "ascend") {
      o <- order(ddx[i, ])
      for (j in low.num:n) {
        if (rejects(o[2:j], ddx[i, o[j]])) return(if (j == low.num) 0 else ddx[i, o[j - 1]])
      }
      return(0)
    }
    o <- order(ddx[i, ], decreasing = TRUE)
    for (j in 1:(n - low.num)) {
      if (!rejects(o[j:(n - 1)], ddx[i, o[j]])) return(ddx[i, o[j]])
    }
    0
  }
  set.seed(33000001)
  for (d in c(2L, 10L)) {
    tb <- t1_table_for(d)
    simul <- tr3_load_table(tb$file)
    n <- 40L
    ok <- 0L; tot <- 0L; mx <- 0
    for (rep in 1:3) {
      X <- matrix(runif(n * d), n, d)
      ddx <- as.matrix(dist(X))
      R_m <- nnccd.radi(dx = X, quantile = "lower", method = tb$direction, low.num = 3,
                        quant = 0.9, simul = simul, niter = 1000)$R
      R_me <- vapply(seq_len(n), function(i) replay_radius(ddx, i, simul, tb$direction), numeric(1))
      ok <- ok + sum(abs(R_m - R_me) < 1e-12); tot <- tot + n; mx <- max(mx, abs(R_m - R_me))
    }
    cat(sprintf("check d=%d (%s, %s): radii equal for %d of %d centres, max |diff| = %g, nonzero radii in the last data set: %d\n",
                d, tb$file, tb$direction, ok, tot, mx, sum(R_m > 0)))
  }
  s <- tr3_sim_setting("UN-CCDs", "uniform", 2, 50)
  p <- list(simul = tr3_load_table(s$table), direction = s$direction)
  same <- vapply(1:3, function(i) {
    X <- matrix(runif(100), 50, 2)
    a <- ccd_cluster("UN-CCDs", X, p); b <- t1_one(X, p)
    a$si.ind == b[["k_hat"]] && a$si == b[["si"]]
  }, logical(1))
  cat(sprintf("check t1_one() vs ccd_cluster(): k_hat and si equal in %d of 3 data sets\n", sum(same)))
}
