#!/usr/bin/env Rscript
# revision_experiments/54_wp2c_smin_rule.R
#
# WP2(c): a LABEL-FREE rule for S_min, and what it costs.
#
# Reviewer points R3.2 / R5.2: S_min is set from the true contamination rate,
# i.e. from the labels the method is supposed to be discovering. That charge
# is valid for the SIMULATION study (main text L605, "with S_min set to half
# the contamination level") and only there -- the real-data S_min is never
# specified in the manuscript at all. This script handles the REAL-DATA side
# (Steps 1-2 of the WP2c brief: plateau structure + candidate rules);
# 55_wp2c_simulation_arm.R handles the simulation side, which is where the
# oracle rule actually lives.
#
# --------------------------------------------------------------------------
# THE ONE STRUCTURAL FACT THIS SCRIPT IS BUILT ON
# --------------------------------------------------------------------------
# min.cls (== S_min) enters the detectors through EXACTLY ONE expression:
#
#   R/ccds/RK_CCD_New.R:113   lenD = length(which(graph$catch >= round(min.cls*n)))   [SU-MCCD]
#   R/ccds/UN_CCD.R:102       lenD = length(which(graph$catch >= round(min.cls*n)))   [SUN-MCCD]
#
# (SU-MCCDs.R:11 -> RKCCD_correct_quant -> rccd.silhouette;
#  SUN-MCCD.R:14 -> nnccd.silhouette. Verified by grep over the whole repo:
#  no other min.cls consumer is on either call path.)
#
# Therefore  k = round(S_min * n)  is a SUFFICIENT STATISTIC for S_min: two
# S_min values that round to the same k give bit-identical output, by
# construction, with no experiment needed. Everything below exploits this:
#   * plateaus are reported in both S_min and k units;
#   * a query S_min (e.g. the contamination rate) that shares its k with an
#     already-measured grid point needs NO new run -- it is the same run;
#   * only genuinely new k values are measured (--measure mode).
# R's round() is half-to-even; we call the same round(), so the k we compute
# is the k the detector computes.
#
# --------------------------------------------------------------------------
# INPUTS (all READ-ONLY -- this script never writes to any of them)
# --------------------------------------------------------------------------
#   results/tr1/wp2a_smin_sweep16.csv   227 ok cells, 12 datasets x 2 methods
#                                       x 10-point grid (PenDigits partial).
#                                       Has flagged_idx -> plateau-capable.
#   results/tr1/wp2a_smin_bigfour.csv   large-set cells, WRITTEN BY A JOB THAT
#                                       MAY STILL BE RUNNING. Read as partial;
#                                       never written.
#   results/tr1/regen_smin_grid.csv     5 datasets, 8-point grid, metrics only
#                                       (no flagged_idx) -> used as a metric
#                                       cross-check, not for plateaus.
#   results/tr1/dataset_inventory.csv   n, d, n0 per dataset.
#
# OUTPUTS (results/tr1/)
#   wp2c_extra_cells.csv        new (dataset, method, S_min) cells measured
#                               here, same schema as wp2a_smin_sweep16.csv
#   wp2c_adaptive_real.csv      adaptive-rule fixed-point traces on real data
#   wp2c_plateau_table.csv      one row per (dataset, method, plateau)
#   wp2c_rule_candidates.csv    one row per (dataset, method, rule)
#
# --------------------------------------------------------------------------
# USAGE (each invocation must finish inside the 600 s shell timeout; every
# result row is appended to disk the moment it is computed)
# --------------------------------------------------------------------------
#   Rscript 54_wp2c_smin_rule.R analyze
#       Plateaus + query resolution from whatever is on disk. Prints the list
#       of (dataset, method, S_min) cells still needed, and writes
#       wp2c_plateau_table.csv.
#
#   Rscript 54_wp2c_smin_rule.R measure [datasets] [budget_sec]
#       Measures the cells "analyze" says are missing (query points whose k
#       is not covered by any measured grid point). Appends per cell.
#
#   Rscript 54_wp2c_smin_rule.R adaptive [datasets] [budget_sec]
#       Runs candidate rule 3 (fraction of the median CCD cluster size) as a
#       fixed-point iteration, for each rho in RHO_GRID. Appends per cell.
#
#   Rscript 54_wp2c_smin_rule.R candidates
#       Ranks the candidate rules; writes wp2c_rule_candidates.csv.

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments", "harness.R")))
suppressMessages(source(here::here("revision_experiments", "wp0_mccd_methods.R")))

REPO   <- here::here()
RESDIR <- file.path(REPO, "revision_experiments/results/tr1")

IN_SWEEP16 <- file.path(RESDIR, "wp2a_smin_sweep16.csv")
IN_BIGFOUR <- file.path(RESDIR, "wp2a_smin_bigfour.csv")
IN_REGEN   <- file.path(RESDIR, "regen_smin_grid.csv")
IN_INVENT  <- file.path(RESDIR, "dataset_inventory.csv")

OUT_EXTRA    <- file.path(RESDIR, "wp2c_extra_cells.csv")
OUT_ADAPT    <- file.path(RESDIR, "wp2c_adaptive_real.csv")
OUT_PLATEAU  <- file.path(RESDIR, "wp2c_plateau_table.csv")
OUT_CANDID   <- file.path(RESDIR, "wp2c_rule_candidates.csv")

METHODS <- c("SU-MCCD", "SUN-MCCD")

# Candidate rule 2 searches this grid of fixed constants for the one that
# maximises plateau membership. Chosen BEFORE any result was looked at:
# the 10-point sweep grid plus the two published/near-published values.
FIXED_GRID <- c(0.005, 0.01, 0.02, 0.03, 0.05, 0.0625, 0.10, 0.15, 0.20, 0.30)

# Candidate rule 3: S_min = rho * (median cluster size) / n.
# rho values fixed a priori: 1/4, 1/2, 1 ("a candidate cluster must hold at
# least a quarter / half / all of what a typical already-found cluster holds").
RHO_GRID <- c(0.25, 0.5, 1.0)

MAX_FIXPOINT_ITER <- 8L
CELL_TIMEOUT_SEC  <- 240

# PenDigits costs ~200 s/cell (measured, wp2a_smin_sweep16.csv elapsed_sec);
# excluded from all NEW measurement here. Its existing 7 cells are still read.
NO_NEW_MEASUREMENT <- c("PenDigits", "pageblocks", "thyroid", "waveform", "wilt")

ROW_COLS <- c("dataset", "method", "s_min", "n", "d",
              "TPR", "TNR", "BA", "F2",
              "n_flagged", "n_clusters", "median_cluster_size", "cluster_sizes",
              "flagged_idx", "elapsed_sec", "status", "note", "timestamp")

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
collapse_num <- function(v) if (length(v) == 0) "" else paste(v, collapse = ";")

read_csv_safe <- function(p) {
  if (!file.exists(p)) return(NULL)
  tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
}

inventory <- function() {
  inv <- read_csv_safe(IN_INVENT)
  stopifnot(!is.null(inv))
  inv
}

# ---------------------------------------------------------------------------
# Load every measured (dataset, method, S_min) cell that carries a flagged set
# ---------------------------------------------------------------------------
load_grid <- function() {
  parts <- list()
  for (nm in c("sweep16", "bigfour", "extra")) {
    p <- switch(nm, sweep16 = IN_SWEEP16, bigfour = IN_BIGFOUR, extra = OUT_EXTRA)
    df <- read_csv_safe(p)
    if (is.null(df) || nrow(df) == 0) next
    df <- df[df$status == "ok", , drop = FALSE]
    if (nrow(df) == 0) next
    df$source <- nm
    parts[[nm]] <- df[, c(ROW_COLS, "source")]
  }
  if (length(parts) == 0) stop("no measured cells found")
  g <- do.call(rbind, parts)
  g$k <- round(g$s_min * g$n)
  # If the same k was measured more than once (different S_min, same k, or the
  # same cell in two files), keep the first and assert the outputs agree --
  # that assertion is itself the empirical check on the sufficient-statistic
  # claim above.
  key <- paste(g$dataset, g$method, g$k, sep = "|")
  dup <- key[duplicated(key)]
  for (kk in unique(dup)) {
    sub <- g[key == kk, ]
    if (length(unique(sub$flagged_idx)) != 1L) {
      stop(sprintf("SUFFICIENT-STATISTIC VIOLATION: %s has %d distinct flagged sets at the same k",
                   kk, length(unique(sub$flagged_idx))))
    }
  }
  g <- g[!duplicated(key), , drop = FALSE]
  g[order(g$dataset, g$method, g$s_min), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Plateaus: maximal contiguous runs of measured grid points (ordered by S_min)
# with a bit-identical flagged set.
# ---------------------------------------------------------------------------
plateaus_for <- function(sub) {
  sub <- sub[order(sub$s_min), , drop = FALSE]
  fi  <- sub$flagged_idx
  run <- cumsum(c(1L, as.integer(fi[-1] != fi[-length(fi)])))
  out <- lapply(unique(run), function(r) {
    ix <- which(run == r)
    data.frame(
      dataset = sub$dataset[1], method = sub$method[1],
      plateau_id = r,
      s_min_lo = min(sub$s_min[ix]), s_min_hi = max(sub$s_min[ix]),
      k_lo = min(sub$k[ix]), k_hi = max(sub$k[ix]),
      n_grid_points = length(ix),
      s_min_points = paste(sub$s_min[ix], collapse = ";"),
      n = sub$n[1], d = sub$d[1],
      TPR = sub$TPR[ix[1]], TNR = sub$TNR[ix[1]],
      BA = sub$BA[ix[1]], F2 = sub$F2[ix[1]],
      n_flagged = sub$n_flagged[ix[1]],
      n_clusters = sub$n_clusters[ix[1]],
      median_cluster_size = sub$median_cluster_size[ix[1]],
      stringsAsFactors = FALSE)
  })
  p <- do.call(rbind, out)
  p$width_s_min <- p$s_min_hi - p$s_min_lo
  p$is_widest_by_points <- p$n_grid_points == max(p$n_grid_points)
  p
}

# Which plateau does a measured k belong to?
plateau_of_k <- function(pl, sub, k) {
  hit <- sub$s_min[abs(sub$k - k) < 0.5]
  if (length(hit) == 0) return(NA_integer_)
  s <- hit[1]
  pid <- pl$plateau_id[sapply(seq_len(nrow(pl)), function(i)
    s %in% as.numeric(strsplit(pl$s_min_points[i], ";")[[1]]))]
  if (length(pid) == 0) NA_integer_ else pid[1]
}

QUERY_NAMES <- c("fixed_0.0625", "contamination", "half_contamination")

query_smins <- function(n, n0) {
  c(fixed_0.0625 = 0.0625, contamination = n0 / n, half_contamination = 0.5 * n0 / n)
}

# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------
do_analyze <- function() {
  g   <- load_grid()
  inv <- inventory()

  pl_all <- list(); q_all <- list(); missing <- list()

  for (ds in sort(unique(g$dataset))) {
    ivr <- inv[inv$dataset == ds, ]
    if (nrow(ivr) == 0) { cat(sprintf("[warn] %s not in inventory; skipped\n", ds)); next }
    n <- ivr$n[1]; n0 <- ivr$n0[1]
    qs <- query_smins(n, n0)
    qk <- round(qs * n)

    for (m in METHODS) {
      sub <- g[g$dataset == ds & g$method == m, , drop = FALSE]
      if (nrow(sub) < 2) next
      pl <- plateaus_for(sub)
      pl_all[[paste(ds, m)]] <- pl

      pids <- sapply(qk, function(k) plateau_of_k(pl, sub, k))
      for (qi in seq_along(qs)) {
        if (is.na(pids[qi])) {
          # k not covered by any measured grid point -> needs a run.
          # Use S_min = k/n exactly, which round()s back to the same k.
          missing[[length(missing) + 1L]] <- data.frame(
            dataset = ds, method = m, query = QUERY_NAMES[qi],
            s_min_needed = qk[qi] / n, k = qk[qi], n = n,
            stringsAsFactors = FALSE)
        }
      }
      q_all[[paste(ds, m)]] <- data.frame(
        dataset = ds, method = m, n = n, n0 = n0, contam = n0 / n,
        k_fixed = qk[1], k_contam = qk[2], k_half = qk[3],
        pid_fixed = pids[1], pid_contam = pids[2], pid_half = pids[3],
        all_three_same_plateau =
          !anyNA(pids) && length(unique(pids)) == 1L,
        fixed_eq_half = !anyNA(pids[c(1, 3)]) && pids[1] == pids[3],
        fixed_eq_contam = !anyNA(pids[c(1, 2)]) && pids[1] == pids[2],
        stringsAsFactors = FALSE)
    }
  }

  PL <- do.call(rbind, pl_all); rownames(PL) <- NULL
  QT <- do.call(rbind, q_all);  rownames(QT) <- NULL

  # merge the query resolution onto the plateau table as a companion block
  PL$flag_fixed  <- FALSE; PL$flag_contam <- FALSE; PL$flag_half <- FALSE
  for (i in seq_len(nrow(QT))) {
    sel <- PL$dataset == QT$dataset[i] & PL$method == QT$method[i]
    if (!is.na(QT$pid_fixed[i]))  PL$flag_fixed[sel  & PL$plateau_id == QT$pid_fixed[i]]  <- TRUE
    if (!is.na(QT$pid_contam[i])) PL$flag_contam[sel & PL$plateau_id == QT$pid_contam[i]] <- TRUE
    if (!is.na(QT$pid_half[i]))   PL$flag_half[sel   & PL$plateau_id == QT$pid_half[i]]   <- TRUE
  }
  write.csv(PL, OUT_PLATEAU, row.names = FALSE)
  cat(sprintf("wrote %s (%d plateau rows over %d dataset x method pairs)\n",
              OUT_PLATEAU, nrow(PL), nrow(QT)))

  cat("\n=== per (dataset, method): do 0.0625, contamination and half-contamination share a plateau? ===\n")
  print(QT[, c("dataset", "method", "contam", "k_fixed", "k_contam", "k_half",
               "pid_fixed", "pid_contam", "pid_half", "all_three_same_plateau")],
        row.names = FALSE)
  cat(sprintf("\nall three in the SAME plateau: %d of %d (dataset, method) pairs\n",
              sum(QT$all_three_same_plateau), nrow(QT)))
  bad <- QT[!QT$all_three_same_plateau, ]
  if (nrow(bad)) {
    cat("pairs where they do NOT coincide:\n")
    print(bad[, c("dataset", "method", "pid_fixed", "pid_contam", "pid_half")], row.names = FALSE)
  }

  # dataset-level: all three in the same plateau for BOTH methods
  byds <- tapply(QT$all_three_same_plateau, QT$dataset, all)
  cat(sprintf("\ndatasets where it holds for BOTH methods: %d of %d -- %s\n",
              sum(byds), length(byds), paste(names(byds)[byds], collapse = ", ")))
  cat(sprintf("datasets where it fails for at least one method: %s\n",
              paste(names(byds)[!byds], collapse = ", ")))

  cat("\n=== cluster-count degeneracy at each query point ===\n")
  deg <- do.call(rbind, lapply(seq_len(nrow(QT)), function(i) {
    sel <- PL$dataset == QT$dataset[i] & PL$method == QT$method[i]
    gv <- function(fl) { r <- PL[sel & PL[[fl]], ]; if (nrow(r)) r$n_clusters[1] else NA_integer_ }
    data.frame(dataset = QT$dataset[i], method = QT$method[i],
               ncls_fixed = gv("flag_fixed"), ncls_contam = gv("flag_contam"),
               ncls_half = gv("flag_half"), stringsAsFactors = FALSE)
  }))
  print(deg, row.names = FALSE)
  cat(sprintf("\nsingle-macro-cluster (n_clusters == 1) at S_min = 0.0625: %d of %d pairs\n",
              sum(deg$ncls_fixed == 1, na.rm = TRUE), sum(!is.na(deg$ncls_fixed))))
  cat(sprintf("single-macro-cluster at S_min = half-contamination:        %d of %d pairs\n",
              sum(deg$ncls_half == 1, na.rm = TRUE), sum(!is.na(deg$ncls_half))))

  if (length(missing)) {
    MS <- do.call(rbind, missing)
    MS <- MS[!(MS$dataset %in% NO_NEW_MEASUREMENT), , drop = FALSE]
    MS <- MS[!duplicated(paste(MS$dataset, MS$method, MS$k)), , drop = FALSE]
    cat(sprintf("\n=== %d measurable cell(s) missing (query k not on the measured grid) ===\n", nrow(MS)))
    if (nrow(MS)) print(MS, row.names = FALSE)
  } else {
    cat("\nno missing cells: every query point's k is already on the measured grid.\n")
  }
  invisible(list(PL = PL, QT = QT, grid = g))
}

# ---------------------------------------------------------------------------
# measure: run the cells analyze reports as missing
# ---------------------------------------------------------------------------
run_cell <- function(dataset, method, s_min, out_csv, note = NA_character_) {
  ex <- read_csv_safe(out_csv)
  if (!is.null(ex) && any(ex$dataset == dataset & ex$method == method &
                          abs(ex$s_min - s_min) < 1e-12)) {
    cat(sprintf("[skip] %s x %s s_min=%.6f already recorded\n", dataset, method, s_min))
    return(invisible("skip"))
  }
  dat <- load_real_dataset(dataset)
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT_SEC, transient = TRUE)
    res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = s_min)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    stopifnot(length(res$score) == n, !anyNA(res$score))
    m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[[method]])
    cl <- res$cluster; cl <- cl[!is.na(cl)]
    tab <- table(cl); cs <- as.integer(tab)
    list(m = m, score = res$score, n_clusters = length(tab),
         mcs = if (length(cs)) median(cs) else NA_real_, cs = cs,
         status = "ok", note = note)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")), score = NULL,
         n_clusters = NA_integer_, mcs = NA_real_, cs = integer(0),
         status = "error", note = substr(conditionMessage(e), 1, 300))
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fl <- if (is.null(out$score)) integer(0) else which(out$score == 1)
  row <- as.data.frame(list(
    dataset = dataset, method = method, s_min = s_min, n = n, d = d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
    n_flagged = if (is.null(out$score)) NA_integer_ else length(fl),
    n_clusters = out$n_clusters, median_cluster_size = out$mcs,
    cluster_sizes = collapse_num(out$cs), flagged_idx = collapse_idx(fl),
    elapsed_sec = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stringsAsFactors = FALSE)[ROW_COLS]
  append_result(out_csv, row)   # <- on disk immediately
  cat(sprintf("  %-12s x %-9s s_min=%-9.6f k=%-4d BA=%s F2=%s ncls=%-3s %s wall=%.1fs\n",
              dataset, method, s_min, round(s_min * n),
              if (is.na(out$m[["BA"]])) "NA" else sprintf("%.3f", out$m[["BA"]]),
              if (is.na(out$m[["F2"]])) "NA" else sprintf("%.3f", out$m[["F2"]]),
              ifelse(is.na(out$n_clusters), "NA", out$n_clusters), out$status, wall))
  flush.console()
  invisible(out$status)
}

do_measure <- function(ds_filter, budget) {
  a <- do_analyze()
  g <- a$grid; inv <- inventory()
  todo <- list()
  for (ds in sort(unique(g$dataset))) {
    if (ds %in% NO_NEW_MEASUREMENT) next
    if (!is.null(ds_filter) && !(ds %in% ds_filter)) next
    ivr <- inv[inv$dataset == ds, ]; if (!nrow(ivr)) next
    n <- ivr$n[1]; n0 <- ivr$n0[1]
    qk <- round(query_smins(n, n0) * n)
    for (m in METHODS) {
      sub <- g[g$dataset == ds & g$method == m, , drop = FALSE]
      if (!nrow(sub)) next
      for (k in unique(qk)) {
        if (!any(abs(sub$k - k) < 0.5)) {
          todo[[length(todo) + 1L]] <- list(ds = ds, m = m, s = k / n)
        }
      }
    }
  }
  cat(sprintf("\n54_wp2c_smin_rule.R measure: %d cell(s) to run\n", length(todo)))
  t0 <- Sys.time()
  for (job in todo) {
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > budget) {
      cat("[budget exhausted -- rerun to continue]\n"); break
    }
    run_cell(job$ds, job$m, job$s, OUT_EXTRA, note = "wp2c query point: k not on the wp2a grid")
  }
}

# ---------------------------------------------------------------------------
# adaptive: candidate rule 3, S_min = rho * median(cluster size) / n
# ---------------------------------------------------------------------------
#
# The rule is circular (the partition depends on S_min). Resolution, fixed in
# advance and reported as-is:
#
#   iterate 0:  S_min^(0) = 0   (the DETECTORS' OWN DEFAULT, min.cls = 0 --
#               no cluster-size floor at all, so lenD counts every dominating
#               ball. This is the neutral seed: it uses no label, no tuning
#               constant, and it is the value the code ships with.)
#   iterate t:  run the detector at S_min^(t); take m_t = median size of the
#               macro-clusters it returns; set S_min^(t+1) = rho * m_t / n.
#   stop when   k^(t+1) == k^(t)  (FIXED POINT), or a cycle repeats, or
#               MAX_FIXPOINT_ITER iterations elapse.
#
# Because k = round(S_min*n) is a sufficient statistic, convergence is tested
# on k, not on the real number S_min, and a cycle in k is detectable exactly.
# Whether it converges is REPORTED, not assumed.
adaptive_trace <- function(dataset, method, rho, budget_left) {
  dat <- load_real_dataset(dataset)
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  s <- 0; seen_k <- c(); converged <- FALSE; cyc <- FALSE
  for (it in 0:MAX_FIXPOINT_ITER) {
    k <- round(s * n)
    t0 <- Sys.time()
    out <- tryCatch({
      setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT_SEC, transient = TRUE)
      res <- METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, min.cls = s)
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      cl <- res$cluster; cl <- cl[!is.na(cl)]
      tab <- table(cl); cs <- as.integer(tab)
      list(m = evaluate(Y, res$score, REAL_DATA_THRESHOLDS[[method]]),
           score = res$score, ncls = length(tab),
           mcs = if (length(cs)) median(cs) else NA_real_, cs = cs, status = "ok")
    }, error = function(e) {
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")), score = NULL,
           ncls = NA_integer_, mcs = NA_real_, cs = integer(0),
           status = paste0("error: ", substr(conditionMessage(e), 1, 200)))
    })
    wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fl <- if (is.null(out$score)) integer(0) else which(out$score == 1)
    s_next <- if (is.na(out$mcs)) NA_real_ else rho * out$mcs / n
    k_next <- if (is.na(s_next)) NA_integer_ else round(s_next * n)
    if (!is.na(k_next) && k_next == k) converged <- TRUE
    if (!is.na(k_next) && k_next %in% seen_k && !converged) cyc <- TRUE
    seen_k <- c(seen_k, k)

    append_result(OUT_ADAPT, as.data.frame(list(
      dataset = dataset, method = method, rho = rho, iter = it,
      s_min = s, k = k, n = n, d = d,
      TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
      BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
      n_flagged = length(fl), n_clusters = out$ncls,
      median_cluster_size = out$mcs, cluster_sizes = collapse_num(out$cs),
      flagged_idx = collapse_idx(fl),
      s_min_next = s_next, k_next = k_next,
      converged = converged, cycled = cyc,
      elapsed_sec = wall, status = out$status,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), stringsAsFactors = FALSE))
    cat(sprintf("  %-12s %-9s rho=%.2f it=%d s=%.6f k=%-4d ncls=%-3s med=%-6s BA=%s -> k_next=%s%s\n",
                dataset, method, rho, it, s, k,
                ifelse(is.na(out$ncls), "NA", out$ncls),
                ifelse(is.na(out$mcs), "NA", out$mcs),
                if (is.na(out$m[["BA"]])) "NA" else sprintf("%.3f", out$m[["BA"]]),
                ifelse(is.na(k_next), "NA", k_next),
                if (converged) "  [FIXED POINT]" else if (cyc) "  [CYCLE]" else ""))
    flush.console()
    if (converged || cyc || is.na(k_next) || out$status != "ok") break
    s <- s_next
  }
  invisible(NULL)
}

adaptive_done <- function(dataset, method, rho) {
  ex <- read_csv_safe(OUT_ADAPT)
  if (is.null(ex)) return(FALSE)
  sub <- ex[ex$dataset == dataset & ex$method == method & abs(ex$rho - rho) < 1e-12, ]
  if (!nrow(sub)) return(FALSE)
  any(sub$converged | sub$cycled | sub$iter == MAX_FIXPOINT_ITER | sub$status != "ok")
}

do_adaptive <- function(ds_filter, budget) {
  inv <- inventory()
  g   <- load_grid()
  dss <- if (is.null(ds_filter)) setdiff(sort(unique(g$dataset)), NO_NEW_MEASUREMENT) else ds_filter
  t0 <- Sys.time()
  for (ds in dss) for (m in METHODS) for (rho in RHO_GRID) {
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > budget) {
      cat("[budget exhausted -- rerun to continue]\n"); return(invisible(NULL))
    }
    if (adaptive_done(ds, m, rho)) {
      cat(sprintf("[skip] %s x %s rho=%.2f already traced\n", ds, m, rho)); next
    }
    adaptive_trace(ds, m, rho, budget)
  }
}

# ---------------------------------------------------------------------------
# candidates: rank the rules
# ---------------------------------------------------------------------------
#
# HIT CRITERION. The brief asks how many datasets each rule "keeps inside the
# measured plateau". Plateau membership is bookkeeping over a particular grid;
# the thing it is a proxy for is simply "does this rule produce the SAME
# FLAGGED SET as the reference?". We therefore score on flagged-set identity,
# which is exact, grid-independent, and coincides with plateau membership
# whenever both S_min values sit on the reference grid. plateau_id is reported
# alongside so the two can be cross-checked.
#
# Two references, neither silently preferred:
#   hit_vs_half   -- same flagged set as the oracle half-contamination rule
#                    (the rule the manuscript's simulation study uses)
#   hit_vs_fixed  -- same flagged set as the published fixed 0.0625 (the value
#                    the real-data runs already used)
#
# CANDIDATE RULE 3 IS REPORTED IN TWO VARIANTS, because the circularity
# resolution matters enormously:
#   adaptive_one_step   S_min = rho * median(cluster size at min.cls = 0) / n,
#                       evaluated ONCE and stopped. Genuinely adaptive.
#   adaptive_fixed_pt   the same map iterated to a fixed point in
#                       k = round(S_min*n).
# Both come from the same trace: iter 0 is the neutral seed, iter 1 is the
# one-step rule, the last iterate is the fixed point.
do_candidates <- function() {
  g <- load_grid(); inv <- inventory()
  ad <- read_csv_safe(OUT_ADAPT)

  # Pool of measured cells for metric/flagged-set lookup: the reference grid
  # PLUS every adaptive-trace iterate (also real runs, at k values the grid
  # does not cover). Plateau ids still come from the reference grid only.
  pool <- g[, c("dataset", "method", "k", "n", "TPR", "TNR", "BA", "F2",
                "n_flagged", "n_clusters", "flagged_idx")]
  if (!is.null(ad) && nrow(ad)) {
    a2 <- ad[ad$status == "ok", ]
    if (nrow(a2)) {
      a2$flagged_idx[is.na(a2$flagged_idx)] <- ""
      pool <- rbind(pool, a2[, c("dataset", "method", "k", "n", "TPR", "TNR",
                                 "BA", "F2", "n_flagged", "n_clusters", "flagged_idx")])
    }
  }
  pool <- pool[!duplicated(paste(pool$dataset, pool$method, pool$k)), ]

  rows <- list()
  for (ds in sort(unique(g$dataset))) {
    ivr <- inv[inv$dataset == ds, ]; if (!nrow(ivr)) next
    n <- ivr$n[1]; n0 <- ivr$n0[1]
    for (m in METHODS) {
      sub <- g[g$dataset == ds & g$method == m, , drop = FALSE]
      if (nrow(sub) < 2) next
      pl  <- plateaus_for(sub)
      pol <- pool[pool$dataset == ds & pool$method == m, , drop = FALSE]
      pid <- function(k) plateau_of_k(pl, sub, k)
      row_at <- function(k) {
        r <- pol[abs(pol$k - k) < 0.5, , drop = FALSE]
        if (!nrow(r)) NULL else r[1, ]
      }
      k_half  <- round(0.5 * n0 / n * n)
      k_fixed <- round(0.0625 * n)
      r_half  <- row_at(k_half); r_fixed <- row_at(k_fixed)
      p_half  <- pid(k_half);    p_fixed <- pid(k_fixed)

      add <- function(rule, param, s_val, extra = "") {
        if (is.na(s_val)) return(invisible(NULL))
        k <- round(s_val * n); r <- row_at(k)
        rows[[length(rows) + 1L]] <<- data.frame(
          dataset = ds, method = m, rule = rule, rule_param = param,
          label_free = !startsWith(rule, "oracle"),
          s_min = s_val, k = k, n = n, n0 = n0,
          measured = !is.null(r),
          plateau_id = pid(k),
          plateau_id_oracle_half = p_half, plateau_id_fixed0625 = p_fixed,
          hit_vs_half  = !is.null(r) && !is.null(r_half)  && r$flagged_idx == r_half$flagged_idx,
          hit_vs_fixed = !is.null(r) && !is.null(r_fixed) && r$flagged_idx == r_fixed$flagged_idx,
          BA  = if (is.null(r)) NA_real_ else r$BA,
          F2  = if (is.null(r)) NA_real_ else r$F2,
          TPR = if (is.null(r)) NA_real_ else r$TPR,
          TNR = if (is.null(r)) NA_real_ else r$TNR,
          n_clusters = if (is.null(r)) NA_integer_ else r$n_clusters,
          n_flagged  = if (is.null(r)) NA_integer_ else r$n_flagged,
          BA_minus_oracle = if (is.null(r) || is.null(r_half)) NA_real_ else r$BA - r_half$BA,
          F2_minus_oracle = if (is.null(r) || is.null(r_half)) NA_real_ else r$F2 - r_half$F2,
          note = extra, stringsAsFactors = FALSE)
      }

      add("oracle_half_contamination", "0.5*n0/n", 0.5 * n0 / n)
      add("oracle_full_contamination", "n0/n",     n0 / n)
      for (fc in FIXED_GRID) add("fixed_constant", sprintf("%.4f", fc), fc)

      if (!is.null(ad)) {
        for (rho in RHO_GRID) {
          tr <- ad[ad$dataset == ds & ad$method == m & abs(ad$rho - rho) < 1e-12, ]
          if (!nrow(tr)) next
          tr <- tr[order(tr$iter), ]
          s_one  <- tr$s_min_next[tr$iter == 0]           # one-step: rho * m0 / n
          last   <- tr[nrow(tr), ]
          s_fixp <- if (!is.na(last$s_min_next)) last$s_min_next else last$s_min
          add("adaptive_one_step", sprintf("rho=%.2f", rho),
              if (length(s_one)) s_one[1] else NA_real_,
              sprintf("seed=min.cls_0 median_cluster_size=%s", tr$median_cluster_size[tr$iter == 0][1]))
          add("adaptive_fixed_pt", sprintf("rho=%.2f", rho), s_fixp,
              sprintf("iters=%d converged=%s cycled=%s k_path=%s",
                      nrow(tr), last$converged, last$cycled,
                      paste(tr$k, collapse = ">")))
        }
      }
    }
  }
  CD <- do.call(rbind, rows); rownames(CD) <- NULL
  write.csv(CD, OUT_CANDID, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", OUT_CANDID, nrow(CD)))

  lf  <- CD[CD$label_free, ]
  key <- paste(lf$rule, lf$rule_param)
  rank <- data.frame(
    rule = tapply(lf$rule, key, function(x) x[1]),
    param = tapply(lf$rule_param, key, function(x) x[1]),
    n_pairs = tapply(lf$hit_vs_half, key, length),
    n_measured = tapply(lf$measured, key, sum),
    hits_vs_oracle_half = tapply(lf$hit_vs_half, key, sum),
    hits_vs_fixed0625 = tapply(lf$hit_vs_fixed, key, sum),
    mean_BA = round(tapply(lf$BA, key, function(x) mean(x, na.rm = TRUE)), 4),
    mean_F2 = round(tapply(lf$F2, key, function(x) mean(x, na.rm = TRUE)), 4),
    mean_dBA_vs_oracle = round(tapply(lf$BA_minus_oracle, key, function(x) mean(x, na.rm = TRUE)), 4),
    mean_dF2_vs_oracle = round(tapply(lf$F2_minus_oracle, key, function(x) mean(x, na.rm = TRUE)), 4),
    frac_single_cluster = round(tapply(lf$n_clusters, key,
                                       function(x) mean(x == 1, na.rm = TRUE)), 3),
    stringsAsFactors = FALSE)
  rank <- rank[order(-rank$hits_vs_oracle_half, -rank$mean_BA), ]
  cat("\n=== candidate ranking (label-free rules); hit = flagged set identical to the reference ===\n")
  print(rank, row.names = FALSE)

  # Restrict to the pairs where the oracle itself is measured, so the
  # comparison is apples-to-apples.
  or <- CD[CD$rule == "oracle_half_contamination", ]
  cat(sprintf("\noracle (half-contamination) reference over %d (dataset, method) pairs: mean BA = %.4f, mean F2 = %.4f; single-cluster in %d of %d\n",
              nrow(or), mean(or$BA, na.rm = TRUE), mean(or$F2, na.rm = TRUE),
              sum(or$n_clusters == 1, na.rm = TRUE), sum(!is.na(or$n_clusters))))

  cat("\n=== adaptive fixed points: where does the iteration land? ===\n")
  fp <- CD[CD$rule == "adaptive_fixed_pt", c("dataset", "method", "rule_param",
                                             "s_min", "k", "n_clusters", "BA", "note")]
  print(fp, row.names = FALSE)
  invisible(CD)
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "analyze"
dsf  <- if (length(args) >= 2 && nzchar(args[2]) && toupper(args[2]) != "ALL")
          strsplit(args[2], ",")[[1]] else NULL
bud  <- if (length(args) >= 3) as.numeric(args[3]) else 480

cat(sprintf("54_wp2c_smin_rule.R: mode=%s datasets=%s budget=%.0fs\n",
            mode, if (is.null(dsf)) "ALL" else paste(dsf, collapse = ","), bud))

switch(mode,
  analyze    = invisible(do_analyze()),
  measure    = do_measure(dsf, bud),
  adaptive   = do_adaptive(dsf, bud),
  candidates = invisible(do_candidates()),
  stop(sprintf("unknown mode '%s'", mode))
)
cat("54_wp2c_smin_rule.R: done.\n")
