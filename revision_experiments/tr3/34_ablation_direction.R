#!/usr/bin/env Rscript
# revision_experiments/tr3/34_ablation_direction.R
#
# Ablation of the radius-search direction of UN-CCDs (reviewer comment 3.3).
# The published method searches each covering-ball radius ascending for d < 10
# and descending for d >= 10 (paper Section 3.2, Algorithm 4). Here the d >= 10
# cells are also run ascending, on the paper's own data sets, and compared with
# the descending (published) replicates already in sim_replicates.csv.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/34_ablation_direction.R run [reps] [cells]
#   Rscript revision_experiments/tr3/34_ablation_direction.R summary
#     [reps]   replicates per cell to run (default 500 = all). A smaller number is a
#              smoke test: it is accepted only when TR3_OUT_DIR points elsewhere than
#              results/tr3, because a short cell would otherwise be marked done.
#              The data are generated in full (500 per cluster number, in the
#              scripts' order) and the first [reps] of each list are used.
#     [cells]  comma-separated subset of <dist>:<d>:<clusters>, e.g.
#              uniform:10:2,gaussian:20:5 (default: all 12 cells)
#   TR3_CORES overrides the number of workers (default detectCores()). The results
#   do not depend on it: data are generated on the main process and the clustering
#   is deterministic.
#   TR3_OUT_DIR overrides the folder of every output (default results/tr3). Inputs
#   (sim_replicates.csv, sim_done.csv, the quantile tables) are always read from the
#   usual places.
#
# Cells (part run): UN-CCDs, d in {10, 20}, n = 200, dist in {uniform, gaussian},
# 2, 3 and 5 clusters, 500 replicates (12 cells). Data: tr3_sim_setting("UN-CCDs",
# dist, d, 200) and tr3_sim_data(): the published scripts' seed (2024), random-number
# order and quantile table (NN-test-simul_10d_99%, NN-test-simul_20d_999%).
# The only change against 10_sim.R is the search direction, forced to "ascend" in
# the parameters handed to tr3_one(); the same table, the same parallel setup.
#
# Before a cell is run (and skipped when the cell is done)
#   1. Fingerprint check. tr3_fp() of every generated data set must equal the
#      data_fp of the stored `published` row of the same replicate (compared as
#      fmt17 text, which is how both were written). Otherwise stop.
#   2. Reproduction check. The DESCENDING arm is rerun on the first 3 replicates
#      and its ARI, Sil and k_hat must equal the stored `published` rows, compared as
#      fmt17 text (17 significant digits, i.e. bit-identical doubles). Otherwise
#      stop. The largest absolute difference is printed and recorded.
#
# The descending arm is not recomputed: its replicates are the `published` rows of
# results/tr3/sim_replicates.csv (cells listed in sim_done.csv; they were rerun on
# 2026-10-23 with the corrected descending index, see README.md).
#
# Part summary: per cell and arm the mean ARI, mean Sil, % of replicates with
# k_hat = true number of clusters, mean seconds per replicate; a paired two-sided
# Wilcoxon signed-rank test on the per-replicate ARI differences (ascend minus
# descend), wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = TRUE), with
# the number of nonzero differences, and Holm-adjusted p-values over the 12 cells
# (p.adjust(..., n = 12), so a partial run is adjusted as if all 12 were present).
# If all differences of a cell are zero the test is undefined and its p-value is NA.
# Seconds: the two arms were timed in different runs (the descending rows on
# 2026-10-23, the ascending rows now), so the comparison of seconds is indicative.
#
# OUTPUT (OUT_DIR, default revision_experiments/results/tr3/)
#   ablation_direction_replicates.csv  one row per ASCEND replicate, the columns of
#                          sim_replicates.csv: variant ("ascend"), method, dist, d, n,
#                          n_clusters, iteN, rep, ari, sil, k_hat, seconds, param ("NA"),
#                          data_fp
#   ablation_direction_done.csv        one row per finished cell; reruns skip these cells
#   ablation_direction_repro.csv       the reproduction check: per cell and replicate
#                          (1 to 3) the stored and the rerun descending values
#   ablation_direction_summary.csv     part summary (one row per cell, both arms)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !(args[1] %in% c("run", "summary"))) {
  stop("usage: Rscript revision_experiments/tr3/34_ablation_direction.R run [reps] [cells] | summary")
}
source(here::here("revision_experiments/tr3/tr3_lib.R"))

PART    <- args[1]
OUT_DIR <- Sys.getenv("TR3_OUT_DIR", TR3_RES)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- normalizePath(OUT_DIR, winslash = "/", mustWork = TRUE)

REP_CSV     <- file.path(OUT_DIR, "ablation_direction_replicates.csv")
DONE_CSV    <- file.path(OUT_DIR, "ablation_direction_done.csv")
REPRO_CSV   <- file.path(OUT_DIR, "ablation_direction_repro.csv")
SUMMARY_CSV <- file.path(OUT_DIR, "ablation_direction_summary.csv")

METHOD <- "UN-CCDs"
NN     <- "200"
KEYS   <- c("variant", "method", "dist", "d", "n", "n_clusters")

# the 12 cells, grouped by (dist, d): one data generation and one table per group
GROUPS <- expand.grid(d = c("10", "20"), dist = c("uniform", "gaussian"), stringsAsFactors = FALSE)
GROUPS <- GROUPS[, c("dist", "d")]
CLUSTERS <- c(2, 3, 5)
ALL_CELLS <- do.call(rbind, lapply(seq_len(nrow(GROUPS)), function(i)
  data.frame(dist = GROUPS$dist[i], d = GROUPS$d[i], n_clusters = as.character(CLUSTERS),
             stringsAsFactors = FALSE)))
stopifnot(nrow(ALL_CELLS) == 12)

# the stored descending rows (variant published) of the 12 cells
published_rows <- function() {
  pub <- tr3_read_done(tr3_path("sim_replicates.csv"), tr3_path("sim_done.csv"), KEYS, "rep")
  if (is.null(pub)) stop("34_ablation_direction.R: no results/tr3/sim_replicates.csv or sim_done.csv")
  pub <- pub[pub$variant == "published" & pub$method == METHOD & pub$n == NN &
               pub$d %in% c("10", "20"), , drop = FALSE]
  stopifnot(nrow(pub) > 0, all(pub$iteN == "500"))
  pub$rep <- as.integer(pub$rep)
  pub[order(pub$dist, as.integer(pub$d), as.integer(pub$n_clusters), pub$rep), , drop = FALSE]
}

# ===========================================================================
if (PART == "run") {
  REPS_ALL <- 500L
  REPS <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else REPS_ALL
  stopifnot(!is.na(REPS), REPS >= 3, REPS <= REPS_ALL)
  if (REPS < REPS_ALL && identical(OUT_DIR, normalizePath(TR3_RES, winslash = "/", mustWork = FALSE))) {
    stop("34_ablation_direction.R: [reps] < 500 writes short cells; set TR3_OUT_DIR to a scratch folder")
  }
  want <- ALL_CELLS
  if (length(args) >= 3 && nzchar(args[3])) {
    sel <- strsplit(args[3], ",")[[1]]
    ok <- sel %in% paste(want$dist, want$d, want$n_clusters, sep = ":")
    if (!all(ok)) stop("34_ablation_direction.R: unknown cell(s): ", paste(sel[!ok], collapse = ", "))
    want <- want[paste(want$dist, want$d, want$n_clusters, sep = ":") %in% sel, , drop = FALSE]
  }

  pub <- published_rows()
  tr3_source_method(METHOD, "simulation")

  cores <- tr3_cores()
  cat(sprintf("34_ablation_direction.R run: %s n=%s reps=%d cells=%d cores=%d out=%s\n",
              METHOD, NN, REPS, nrow(want), cores, OUT_DIR))
  t_start <- Sys.time()

  cl <- makeCluster(cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  cell_key <- function(DIST, D, N) list(variant = "ascend", method = METHOD, dist = DIST, d = D, n = NN,
                                        n_clusters = as.character(N))

  # Every foreach() below is called at top level, as in 10_sim.R: foreach exports to
  # the workers only what the calling environment holds, following global functions
  # into their bodies (tr3_one -> ccd_cluster -> the method).
  for (g in seq_len(nrow(GROUPS))) {
    DIST <- GROUPS$dist[g]; D <- GROUPS$d[g]
    Ns <- want$n_clusters[want$dist == DIST & want$d == D]
    Ns <- Ns[!vapply(Ns, function(N) has_result(DONE_CSV, cell_key(DIST, D, N)), logical(1))]
    if (!length(Ns)) {
      cat(sprintf("%s d=%s: done or not selected, skipped\n", DIST, D))
      next
    }
    s <- tr3_sim_setting(METHOD, DIST, D, NN)
    stopifnot(s$direction == "descend", s$iteN == as.character(REPS_ALL), s$n == NN)
    par_desc <- list(simul = tr3_load_table(s$table), direction = "descend")
    par_asc  <- list(simul = par_desc$simul, direction = "ascend")   # the only difference
    data <- tr3_sim_data(s)       # always all three cluster numbers, in the scripts' order
    cat(sprintf("%s d=%s: table=%s, published direction=%s, running direction=%s\n",
                DIST, D, s$table, s$direction, par_asc$direction))

    for (N in Ns) {
      key <- cell_key(DIST, D, N)
      g_dat <- data[[N]]
      data_list <- g_dat$data[seq_len(REPS)]
      label_T <- g_dat$label_T
      st <- pub[pub$dist == DIST & pub$d == D & pub$n_clusters == N, , drop = FALSE]
      if (!identical(st$rep[st$rep <= REPS], seq_len(REPS))) {
        stop(sprintf("34_ablation_direction.R: stored published rows of %s d=%s %s clusters do not cover replicates 1..%d",
                     DIST, D, N, REPS))
      }
      st <- st[seq_len(REPS), , drop = FALSE]
      t0 <- Sys.time()

      # 1. fingerprints of the generated data against the stored published rows
      fp <- fmt17(vapply(data_list, tr3_fp, numeric(1)))
      bad <- which(fp != st$data_fp)
      if (length(bad)) {
        stop(sprintf("34_ablation_direction.R: data_fp differs from the stored published rows in %s d=%s %s clusters at replicate(s) %s",
                     DIST, D, N, paste(head(bad, 10), collapse = ", ")))
      }
      cat(sprintf("  %s clusters: fingerprints of %d data sets equal the stored rows\n", N, REPS))

      # 2. reproduction of the descending arm on the first 3 replicates
      NR <- 3L
      p <- par_desc
      rr <- foreach(x = data_list[seq_len(NR)], .packages = c("MASS", "cluster")) %dopar% tr3_one(x, label_T, METHOD, p)
      rr <- do.call(rbind, rr)
      chk <- data.frame(variant = "published", method = METHOD, dist = DIST, d = D, n = NN, n_clusters = N,
                        rep = seq_len(NR),
                        ari_stored = st$ari[seq_len(NR)], ari_rerun = fmt17(rr[, "ari"]),
                        sil_stored = st$sil[seq_len(NR)], sil_rerun = fmt17(rr[, "sil"]),
                        k_hat_stored = st$k_hat[seq_len(NR)], k_hat_rerun = fmt17(rr[, "k_hat"]),
                        max_abs_diff = fmt17(pmax(abs(rr[, "ari"] - as.numeric(st$ari[seq_len(NR)])),
                                                  abs(rr[, "sil"] - as.numeric(st$sil[seq_len(NR)])),
                                                  abs(rr[, "k_hat"] - as.numeric(st$k_hat[seq_len(NR)])))),
                        seconds_stored = st$seconds[seq_len(NR)], seconds_rerun = fmt17(rr[, "seconds"]),
                        stringsAsFactors = FALSE)
      chk$equal <- chk$ari_stored == chk$ari_rerun & chk$sil_stored == chk$sil_rerun &
        chk$k_hat_stored == chk$k_hat_rerun
      tr3_append(REPRO_CSV, chk)
      if (!all(chk$equal)) {
        stop(sprintf("34_ablation_direction.R: the descending rerun differs from the stored rows in %s d=%s %s clusters (replicate(s) %s; largest abs difference %s); see %s",
                     DIST, D, N, paste(chk$rep[!chk$equal], collapse = ", "),
                     format(max(as.numeric(chk$max_abs_diff)), digits = 3), REPRO_CSV))
      }
      cat(sprintf("  %s clusters: descending rerun on replicates 1-%d equals the stored rows (fmt17 text; largest abs difference %s)\n",
                  N, NR, format(max(as.numeric(chk$max_abs_diff)), digits = 3)))

      # the ascending arm
      p <- par_asc
      stopifnot(identical(p$direction, "ascend"))
      r <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% tr3_one(x, label_T, METHOD, p)
      r <- do.call(rbind, r)

      rows <- data.frame(key, iteN = s$iteN, rep = seq_len(nrow(r)),
                         ari = fmt17(r[, "ari"]), sil = fmt17(r[, "sil"]),
                         k_hat = fmt17(r[, "k_hat"]), seconds = fmt17(r[, "seconds"]),
                         param = "NA", data_fp = fp, stringsAsFactors = FALSE)
      tr3_append(REP_CSV, rows)
      tr3_append(DONE_CSV, data.frame(key, iteN = s$iteN, n_rows = nrow(rows), cores = cores,
                                      wall_secs = fmt17(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                                      finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      stringsAsFactors = FALSE))
      cat(sprintf("  %s clusters, ascend, %s d=%s: ", N, DIST, D))
      tr3_print_summary(METHOD, r[, "ari"], r[, "sil"], r[, "k_hat"], as.numeric(N))
    }
  }
  print(Sys.time() - t_start)
}

# ===========================================================================
if (PART == "summary") {
  asc <- tr3_read_done(REP_CSV, DONE_CSV, KEYS, "rep")
  if (is.null(asc)) stop("34_ablation_direction.R summary: no finished ascend cells in ", OUT_DIR)
  asc <- asc[asc$variant == "ascend", , drop = FALSE]
  pub <- published_rows()
  num <- function(x) as.numeric(x)

  out <- list()
  for (i in seq_len(nrow(ALL_CELLS))) {
    C <- ALL_CELLS[i, ]
    a <- asc[asc$dist == C$dist & asc$d == C$d & asc$n_clusters == C$n_clusters, , drop = FALSE]
    if (!nrow(a)) next
    a$rep <- as.integer(a$rep)
    a <- a[order(a$rep), , drop = FALSE]
    p <- pub[pub$dist == C$dist & pub$d == C$d & pub$n_clusters == C$n_clusters, , drop = FALSE]
    p <- p[match(a$rep, p$rep), , drop = FALSE]
    if (anyNA(p$rep)) stop("34_ablation_direction.R summary: a stored published replicate is missing in ", C$dist, " d=", C$d, " ", C$n_clusters)
    if (!identical(a$data_fp, p$data_fp)) stop("34_ablation_direction.R summary: data_fp differs between the arms in ", C$dist, " d=", C$d, " ", C$n_clusters)
    N <- as.numeric(C$n_clusters)
    diff <- num(a$ari) - num(p$ari)                        # ascend minus descend
    nz <- sum(diff != 0)
    pw <- if (nz > 0) {
      stats::wilcox.test(num(a$ari), num(p$ari), paired = TRUE, exact = FALSE, correct = TRUE)$p.value
    } else NA_real_
    out[[length(out) + 1]] <- data.frame(
      dist = C$dist, d = C$d, n = NN, n_clusters = C$n_clusters, reps = nrow(a),
      descend_mean_ari = fmt17(mean(num(p$ari))), ascend_mean_ari = fmt17(mean(num(a$ari))),
      descend_mean_sil = fmt17(mean(num(p$sil))), ascend_mean_sil = fmt17(mean(num(a$sil))),
      descend_pct_correct_k = fmt17(mean(num(p$k_hat) == N) * 100),
      ascend_pct_correct_k = fmt17(mean(num(a$k_hat) == N) * 100),
      descend_mean_seconds = fmt17(mean(num(p$seconds))), ascend_mean_seconds = fmt17(mean(num(a$seconds))),
      mean_ari_diff = fmt17(mean(diff)), n_nonzero_diff = nz,
      p_wilcoxon = fmt17(pw), stringsAsFactors = FALSE)
  }
  if (!length(out)) stop("34_ablation_direction.R summary: nothing to summarise")
  S <- do.call(rbind, out)
  S$p_holm <- fmt17(stats::p.adjust(suppressWarnings(num(S$p_wilcoxon)), method = "holm", n = 12))
  utils::write.csv(S, SUMMARY_CSV, row.names = FALSE, quote = TRUE)

  show <- S
  for (k in setdiff(names(S), c("dist", "d", "n", "n_clusters", "reps", "n_nonzero_diff"))) {
    v <- suppressWarnings(num(S[[k]]))
    show[[k]] <- ifelse(is.na(v), "NA", sprintf(if (k %in% c("p_wilcoxon", "p_holm")) "%.4g" else "%.3f", v))
  }
  print(show, row.names = FALSE)
  cat(sprintf("wrote %s (%d of 12 cells)\n", SUMMARY_CSV, nrow(S)))
}
