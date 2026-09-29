#!/usr/bin/env Rscript
# revision_experiments/tr1/90_wp7_clustering_quality.R
#
# WP7 -- clustering quality (R3.6). See WP7_PROTOCOL.md for the pre-declared
# design. SCORING ONLY: this script runs no detector. It reads the
# cluster-assignment files WP8's four scripts (85-88) already write, scores
# the MCCD methods' partitions against ground truth (ARI, NMI, AMI, k-hat,
# unassigned fraction), and separately regenerates each cell's data from its
# recorded seed to run DBSCAN and HDBSCAN as a density-clustering comparison
# on the SAME data.
#
# Usage:
#   Rscript 90_wp7_clustering_quality.R --smoke
#     Processes results/tr1/wp8/smoke/ (85 and 88 have no smoke file as of
#     2026-09-05 -- they were launched straight to production -- so this
#     covers 86 and 87). Writes to results/tr1/wp7/smoke/.
#   Rscript 90_wp7_clustering_quality.R [path] [budget] [--scripts=85,86,87,88]
#     path:   directory holding the four *_clusters.csv files
#             (default results/tr1/wp8/)
#     budget: stop starting new cells after this many seconds (default 1800)
#   Rscript 90_wp7_clustering_quality.R --summarize [--smoke]
#     Reads wp7_metrics_long.csv and writes wp7_summary_metrics.csv and
#     wp7_khat_table.csv (in the smoke output dir if --smoke is also given).
#
# Never edits shared/harness.R, wp0_mccd_methods.R, methods/, R/, or any of
# 85-88 (those are sourced read-only, see load_wp8_env() below).
#
# DO NOT START A PRODUCTION RUN BEFORE ALL FOUR WP8 GRIDS (85-88) HAVE
# FINISHED. An "incomplete"/partial cell written mid-flight (WP8 still
# appending rows to a cluster file WP7 is reading) is never revisited by
# WP7's resumable has_result() gating -- it is scored once, from whatever
# rows existed at read time, and stays that way (WP7_VERIFICATION.md
# change 5). See WP7_PROTOCOL.md's dated 2026-09-06 note under section 1.

suppressMessages(library(here))
suppressPackageStartupMessages(library(mclust))
source(here::here("revision_experiments", "shared", "harness.R"))

# ---------------------------------------------------------------------------
# 1. Metric implementations (WP7_PROTOCOL.md section 3).
#    ARI: mclust::adjustedRandIndex (installed; aricode is not, so NMI/AMI
#    are implemented here from the contingency table and validated against
#    sklearn -- see the protocol and the task report for the validation
#    table).
# ---------------------------------------------------------------------------
wp7_entropy <- function(counts) {
  N <- sum(counts); p <- counts[counts > 0] / N
  -sum(p * log(p))
}

wp7_mutual_information <- function(cm) {
  N <- sum(cm); a <- rowSums(cm); b <- colSums(cm)
  mi <- 0
  for (i in seq_len(nrow(cm))) for (j in seq_len(ncol(cm))) {
    nij <- cm[i, j]
    if (nij > 0) mi <- mi + (nij / N) * log((N * nij) / (a[i] * b[j]))
  }
  mi
}

# Expected mutual information under the hypergeometric null (Vinh, Epps &
# Bailey, 2010), needed for AMI. O(R*C*max_range) -- trivial at N=200.
wp7_expected_mi <- function(a, b, N) {
  la <- lfactorial(a); lNa <- lfactorial(N - a)
  lb <- lfactorial(b); lNb <- lfactorial(N - b)
  lN <- lfactorial(N)
  emi <- 0
  for (i in seq_along(a)) for (j in seq_along(b)) {
    ai <- a[i]; bj <- b[j]
    lo <- max(1, ai + bj - N); hi <- min(ai, bj)
    if (lo > hi) next
    for (nij in lo:hi) {
      term_log <- la[i] + lb[j] + lNa[i] + lNb[j] - lN -
        lfactorial(nij) - lfactorial(ai - nij) - lfactorial(bj - nij) -
        lfactorial(N - ai - bj + nij)
      term <- exp(term_log)
      mi_term <- (nij / N) * log((N * nij) / (ai * bj))
      emi <- emi + term * mi_term
    }
  }
  emi
}

#' ARI/NMI/AMI for one pair of label vectors (no NA allowed -- resolve
#' unassigned points before calling this, see wp7_score_pair()).
#' Returns c(ari=, nmi=, ami=).
wp7_clustering_metrics <- function(true_lab, pred_lab) {
  ut <- sort(unique(true_lab)); up <- sort(unique(pred_lab))
  cm <- table(factor(true_lab, levels = ut), factor(pred_lab, levels = up))
  # BUG FOUND WHILE WRITING THIS (see WP7_PROTOCOL.md section 3): table()'s
  # dimnames propagate through `cm[i,j]` into named scalars, which corrupts
  # c(nmi=nmi, ami=ami) into names "nmi.<level>" -- r["nmi"] then silently
  # returns NA. Strip names immediately, before any indexing.
  cm <- unname(as.matrix(cm))
  a <- unname(rowSums(cm)); b <- unname(colSums(cm)); N <- sum(cm)
  Hu <- wp7_entropy(a); Hv <- wp7_entropy(b)
  mi <- wp7_mutual_information(cm)

  if (Hu == 0 && Hv == 0) {
    nmi <- 1
  } else if (Hu == 0 || Hv == 0) {
    nmi <- 0
  } else {
    nmi <- 2 * mi / (Hu + Hv)
  }

  if (Hu == 0 && Hv == 0) {
    ami <- 1
  } else {
    emi <- wp7_expected_mi(a, b, N)
    denom <- (Hu + Hv) / 2 - emi
    ami <- if (abs(denom) < 1e-12) 0 else (mi - emi) / denom
  }

  ari <- tryCatch(mclust::adjustedRandIndex(true_lab, pred_lab),
                  error = function(e) NA_real_)
  if (is.na(ari) || is.nan(ari)) ari <- 0

  c(ari = unname(ari), nmi = unname(nmi), ami = unname(ami))
}

# ---------------------------------------------------------------------------
# 2. Per-cell scoring: true_cluster/detected_cluster (both length n_reg,
#    detected_cluster may contain NA = unassigned) -> both unassigned
#    treatments (WP7_PROTOCOL.md section 4) + k-hat (section 5).
# ---------------------------------------------------------------------------
wp7_score_pair <- function(row_index, true_cluster, detected_cluster) {
  n_reg <- length(true_cluster)
  is_na <- is.na(detected_cluster)
  n_unassigned <- sum(is_na)
  unassigned_frac <- n_unassigned / n_reg
  true_k <- length(unique(true_cluster))
  k_hat <- length(unique(detected_cluster[!is_na]))

  score_one <- function(t, p) {
    if (length(t) < 2 || length(unique(t)) < 1 || length(unique(p)) < 1) {
      return(list(ari = 0, nmi = 0, ami = 0, status = "degenerate",
                  note = sprintf("n_scored=%d: too few points to score", length(t))))
    }
    m <- tryCatch(wp7_clustering_metrics(t, p),
                  error = function(e) c(ari = NA_real_, nmi = NA_real_, ami = NA_real_))
    if (anyNA(m)) {
      list(ari = 0, nmi = 0, ami = 0, status = "degenerate",
           note = "clustering_metrics() returned NA; sentinel 0 written")
    } else {
      list(ari = unname(m["ari"]), nmi = unname(m["nmi"]), ami = unname(m["ami"]),
           status = "ok", note = "-")
    }
  }

  # singleton: unassigned points get unique labels (never agree with anything)
  pred_singleton <- ifelse(is_na, paste0("UNASSN_", row_index), as.character(detected_cluster))
  m_singleton <- score_one(as.character(true_cluster), pred_singleton)

  # excluded: unassigned points dropped entirely
  keep <- !is_na
  m_excluded <- if (sum(keep) < 2) {
    list(ari = 0, nmi = 0, ami = 0, status = "degenerate",
         note = sprintf("only %d of %d points assigned; excluded treatment has <2 to score",
                         sum(keep), n_reg))
  } else {
    score_one(as.character(true_cluster[keep]), as.character(detected_cluster[keep]))
  }

  list(n_reg = n_reg, n_unassigned = n_unassigned, unassigned_frac = unassigned_frac,
       true_k = true_k, k_hat = k_hat, singleton = m_singleton, excluded = m_excluded)
}

# ---------------------------------------------------------------------------
# 3. WP8 script access -- read-only sys.source() of 85/86/87/88 for their
#    generator functions (WP7_PROTOCOL.md section 6). Each is sourced into
#    its OWN environment: all four define GENERATORS/CANON/BASE_SEED/etc.
#    under the same names, so sharing one environment would clobber them.
#    options(wp8.no_main=TRUE) prevents any of the four from launching a
#    live run or writing to results/tr1/wp8/.
# ---------------------------------------------------------------------------
WP8_SCRIPT_FILE <- c(
  "85" = "85_wp8_boundary_fp.R",
  "86" = "86_wp8_outlier_types.R",
  "87" = "87_wp8_small_cluster.R",
  "88" = "88_wp8_csr_violation.R"
)
.wp8_env_cache <- new.env()
load_wp8_env <- function(script_id) {
  key <- script_id
  if (!is.null(.wp8_env_cache[[key]])) return(.wp8_env_cache[[key]])
  e <- new.env(parent = globalenv())
  op <- options(wp8.no_main = TRUE)
  on.exit(options(op))
  sys.source(here::here("revision_experiments", "tr1", WP8_SCRIPT_FILE[[script_id]]), envir = e)
  .wp8_env_cache[[key]] <- e
  e
}

#' Regenerate one cell's data exactly as the WP8 driver would have, using
#' the recorded seed (never recomputed from CANON's row order). Returns
#' list(X =, true_cluster =, n0 =, n_reg =).
#'
#' WP7_VERIFICATION.md change 1 (blocking): for 85/86, X is the FULL
#' regenerated matrix (n_reg regular rows FOLLOWED BY n0 outlier rows, per
#' gen_*_boundary()/gen_bridge() etc.'s own `rbind(data1, data2, outlier)`
#' convention) -- exactly what the WP8 drivers themselves pass to
#' METHOD_REGISTRY[[method]] (85_wp8_boundary_fp.R / 86_wp8_outlier_types.R,
#' `X <- dat$X`, no truncation). true_cluster stays regular-points-only
#' (length n_reg) since that is all the cluster file ever records; the
#' caller must run DBSCAN/HDBSCAN on the FULL X and subset the returned
#' labels to seq_len(n_reg) before scoring, mirroring the drivers'
#' `res$cluster[seq_len(n_reg)]`. 87/88 already generate n_reg-only X (no
#' contamination in either generator), so this distinction is a no-op there.
regenerate_cell <- function(script_id, env, generator_or_m, d, seed) {
  if (script_id == "85" || script_id == "86") {
    dat <- env$GENERATORS[[generator_or_m]](seed, env$N_NOMINAL, d, env$CONT)
    n_reg <- dat$n1 + dat$n2
    list(X = dat$X,
         true_cluster = dat$true_cluster[seq_len(n_reg)],
         n0 = dat$n0, n_reg = n_reg)
  } else if (script_id == "87") {
    m <- as.integer(generator_or_m)
    dat <- env$gen_small_cluster(seed, env$N_TOTAL, d, m)
    list(X = dat$X, true_cluster = dat$true_cluster, n0 = 0L, n_reg = dat$n)
  } else if (script_id == "88") {
    X <- env$GENERATORS[[generator_or_m]](seed, env$N_NOMINAL, d)
    list(X = X, true_cluster = rep(1L, nrow(X)), n0 = 0L, n_reg = nrow(X))
  } else stop("unknown script_id: ", script_id)
}

# ---------------------------------------------------------------------------
# 4. DBSCAN (native labels) and HDBSCAN (python round trip).
#    WP7_PROTOCOL.md section 6 for the cont= rule per script and the
#    python-not-R rationale for HDBSCAN.
# ---------------------------------------------------------------------------
run_dbscan_native <- function(X, cont) {
  labels <- DBSCAN(X, k = 4, quant = cont)  # from harness.R's sourced DBSCAN.R; 0 = noise
  ifelse(labels == 0, NA_integer_, as.integer(labels))
}

WP7_TMP_DIR <- here::here("revision_experiments/results/tr1/wp7/tmp")
PYTHON_EXE  <- here::here("revision_experiments/.venv/python.exe")
HDBSCAN_PY  <- here::here("revision_experiments/tr1/90b_wp7_hdbscan.py")

run_hdbscan_python <- function(X, tag) {
  dir.create(WP7_TMP_DIR, recursive = TRUE, showWarnings = FALSE)
  in_csv  <- file.path(WP7_TMP_DIR, sprintf("in_%s.csv", tag))
  out_csv <- file.path(WP7_TMP_DIR, sprintf("out_%s.csv", tag))
  on.exit({ suppressWarnings(file.remove(in_csv)); suppressWarnings(file.remove(out_csv)) }, add = TRUE)
  write.csv(as.data.frame(X), in_csv, row.names = FALSE)
  status <- tryCatch(
    system2(PYTHON_EXE, args = shQuote(c(HDBSCAN_PY, in_csv, out_csv)),
            stdout = TRUE, stderr = TRUE),
    error = function(e) { attr(e, "wp7_error") <- TRUE; e })
  ok <- file.exists(out_csv)
  if (!ok) {
    msg <- if (is.list(status) || is.character(status)) paste(status, collapse = "; ") else "no output file"
    return(list(labels = NULL, status = "error", note = substr(paste("hdbscan python call failed:", msg), 1, 200)))
  }
  lab <- tryCatch(read.csv(out_csv, stringsAsFactors = FALSE)$label,
                  error = function(e) NULL)
  if (is.null(lab) || length(lab) != nrow(X)) {
    return(list(labels = NULL, status = "error",
                note = sprintf("hdbscan output has %s rows, expected %d", length(lab), nrow(X))))
  }
  list(labels = ifelse(lab == -1L, NA_integer_, as.integer(lab)), status = "ok", note = "-")
}

# ---------------------------------------------------------------------------
# 5. Cluster-file reading, de-dup, completeness (WP7_PROTOCOL.md section 2)
# ---------------------------------------------------------------------------
CLUSTER_FILE <- list(
  "85" = list(prod = "85_boundary_fp_clusters.csv",  smoke = "85_wp8_boundary_fp_clusters.csv",
              setting_col = "setting_id", expected_n = 189L, true_k = 2L),
  "86" = list(prod = "86_outlier_types_clusters.csv", smoke = "86_wp8_outlier_types_clusters.csv",
              setting_col = "setting_id", expected_n = 189L, true_k = 2L),
  "87" = list(prod = "87_small_cluster_clusters.csv", smoke = "87_wp8_small_cluster_clusters.csv",
              setting_col = "m",         expected_n = 200L, true_k = 3L),
  "88" = list(prod = "88_csr_violation_clusters.csv", smoke = "88_wp8_csr_violation_clusters.csv",
              setting_col = "setting_id", expected_n = 200L, true_k = 1L)
)

#' setting_key: the value used everywhere downstream to identify a
#' (generator/type/m) x d combination. For 85/86/88 this is setting_id
#' itself; for 87 it is m (kept as the bare integer, joined with d later).
wp7_setting_key <- function(script_id, df) {
  if (script_id == "87") as.character(df$m) else as.character(df$setting_id)
}

#' generator/type key to feed regenerate_cell()'s GENERATORS[[...]] lookup:
#' for 85/86/88, setting_id is "<generator>_d<d>" -- strip the suffix. 87
#' has no such column; m is passed through unchanged.
wp7_generator_key <- function(script_id, setting_id, d) {
  if (script_id == "87") return(NA_character_)
  sub(sprintf("_d%d$", d), "", setting_id)
}

read_cluster_file <- function(script_id, dir, smoke) {
  meta <- CLUSTER_FILE[[script_id]]
  fname <- if (smoke) meta$smoke else meta$prod
  path <- file.path(dir, fname)
  if (!file.exists(path)) {
    message(sprintf("[%s] no cluster file at %s -- skipping", script_id, path))
    return(NULL)
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  setting_col <- meta$setting_col
  key <- paste(df[[setting_col]], df$d, df$rep, df$method, df$row_index, sep = "\r")
  n_before <- nrow(df)
  df <- df[!duplicated(key, fromLast = TRUE), ]
  n_dup <- n_before - nrow(df)
  if (n_dup > 0) message(sprintf("[%s] dropped %d duplicate row(s) on (setting,d,rep,method,row_index)", script_id, n_dup))

  df$setting_key <- wp7_setting_key(script_id, df)
  df$cell_key <- paste(df$setting_key, df$d, df$rep, df$method, sep = "\r")

  # Completeness assertion per cell (section 2): exactly expected_n distinct
  # row_index values. A failing cell is flagged, not dropped.
  agg <- stats::aggregate(row_index ~ cell_key, data = df, FUN = function(x) length(unique(x)))
  names(agg) <- c("cell_key", "n_rows")
  bad <- agg$cell_key[agg$n_rows != meta$expected_n]
  if (length(bad) > 0) {
    message(sprintf("[%s] %d cell(s) with row count != %d (incomplete cluster block, e.g. crash between metrics row and cluster block)",
                    script_id, length(bad), meta$expected_n))
  }
  list(df = df, meta = meta, incomplete_keys = bad)
}

# ---------------------------------------------------------------------------
# 6. Output writer -- one 2-row block per (script,setting_key,d,rep,method),
#    gated by has_result() on the "excluded" row as the completeness marker
#    (WP7_PROTOCOL.md section 7).
# ---------------------------------------------------------------------------
WP7_COLS <- c("script", "setting_key", "d", "rep", "seed", "method", "unassigned_treatment",
              "n_reg", "n_unassigned", "unassigned_frac", "true_k", "k_hat",
              "ari", "nmi", "ami", "status", "note")

wp7_done <- function(out_csv, script, setting_key, d, rep_id, method) {
  isTRUE(has_result(out_csv, c(script = script, setting_key = setting_key, d = d,
                               rep = rep_id, method = method, unassigned_treatment = "excluded")))
}

wp7_write_cell <- function(out_csv, script, setting_key, d, rep_id, seed, method, scored) {
  mk_row <- function(treat, m) {
    list(script = script, setting_key = setting_key, d = d, rep = rep_id, seed = seed,
         method = method, unassigned_treatment = treat,
         n_reg = scored$n_reg, n_unassigned = scored$n_unassigned,
         unassigned_frac = scored$unassigned_frac, true_k = scored$true_k, k_hat = scored$k_hat,
         ari = m$ari, nmi = m$nmi, ami = m$ami, status = m$status, note = m$note)
  }
  r1 <- mk_row("singleton", scored$singleton)
  r2 <- mk_row("excluded", scored$excluded)
  block <- Map(function(...) c(...), r1, r2)
  append_result(out_csv, block)
}

# ---------------------------------------------------------------------------
# 7. Main run: MCCD-derived rows (from the cluster file directly) + the
#    DBSCAN/HDBSCAN regeneration pass (one per distinct (setting,d,rep)
#    cell, shared across whichever MCCD methods are present).
# ---------------------------------------------------------------------------
do_run <- function(wp8_dir, out_csv, budget, scripts_sel, smoke) {
  t_start <- Sys.time()
  time_left <- function() budget - as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  for (script_id in scripts_sel) {
    if (time_left() <= 0) { cat(sprintf("[budget exhausted before %s]\n", script_id)); break }
    loaded <- read_cluster_file(script_id, wp8_dir, smoke)
    if (is.null(loaded)) next
    df <- loaded$df; meta <- loaded$meta
    methods_seen <- sort(unique(df$method))
    cat(sprintf("[%s] %d rows, %d cells, methods: %s\n", script_id,
                nrow(df), length(unique(df$cell_key)), paste(methods_seen, collapse = ", ")))

    env <- load_wp8_env(script_id)

    # --- 7a. MCCD-derived rows, straight from the cluster file ---
    for (cell_key in unique(df$cell_key)) {
      if (time_left() <= 0) { cat(sprintf("[budget exhausted mid-%s, MCCD pass]\n", script_id)); break }
      sub <- df[df$cell_key == cell_key, ]
      setting_key <- sub$setting_key[1]; d <- sub$d[1]; rep_id <- sub$rep[1]
      method <- sub$method[1]; seed <- sub$seed[1]

      if (wp7_done(out_csv, script_id, setting_key, d, rep_id, method)) next

      if (cell_key %in% loaded$incomplete_keys) {
        scored <- list(n_reg = nrow(sub), n_unassigned = sum(is.na(sub$detected_cluster)),
                        unassigned_frac = NA_real_, true_k = meta$true_k,
                        k_hat = length(unique(sub$detected_cluster[!is.na(sub$detected_cluster)])),
                        singleton = list(ari = 0, nmi = 0, ami = 0, status = "incomplete",
                                          note = sprintf("row count %d != expected %d", nrow(sub), meta$expected_n)),
                        excluded = list(ari = 0, nmi = 0, ami = 0, status = "incomplete",
                                         note = sprintf("row count %d != expected %d", nrow(sub), meta$expected_n)))
        scored$unassigned_frac <- scored$n_unassigned / scored$n_reg
        wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, method, scored)
        next
      }

      ord <- order(sub$row_index)
      scored <- wp7_score_pair(sub$row_index[ord], sub$true_cluster[ord], sub$detected_cluster[ord])
      wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, method, scored)
    }

    # --- 7b. DBSCAN/HDBSCAN, one regeneration per distinct (setting,d,rep) cell ---
    cell_tab <- unique(df[, c("setting_key", "d", "rep", "seed")])
    for (ci in seq_len(nrow(cell_tab))) {
      if (time_left() <= 0) { cat(sprintf("[budget exhausted mid-%s, DBSCAN/HDBSCAN pass]\n", script_id)); break }
      setting_key <- cell_tab$setting_key[ci]; d <- cell_tab$d[ci]
      rep_id <- cell_tab$rep[ci]; seed <- cell_tab$seed[ci]

      need_dbscan  <- !wp7_done(out_csv, script_id, setting_key, d, rep_id, "DBSCAN")
      need_hdbscan <- !wp7_done(out_csv, script_id, setting_key, d, rep_id, "HDBSCAN")
      if (!need_dbscan && !need_hdbscan) next

      # WP7_VERIFICATION.md change 2 (blocking): a done cell can have a
      # short MCCD block (missing cluster rows after a crash between the
      # metrics row and the cluster block, WP8_REVERIFICATION.md) that is
      # flagged in loaded$incomplete_keys (section 2). Such a method's rows
      # are excluded from BOTH the reference choice and the cross-method
      # true_cluster check below -- an incomplete block's row_index set is
      # a strict subset, and comparing it directly against a complete
      # block's true_cluster vector fails identical() on length alone, which
      # previously wrote DBSCAN/HDBSCAN as status="error" for the WHOLE cell
      # even though the cell's other (complete) methods verify cleanly.
      ref_rows_all <- df[df$setting_key == setting_key & df$d == d & df$rep == rep_id, ]
      ref_rows <- ref_rows_all[!(ref_rows_all$cell_key %in% loaded$incomplete_keys), ]
      excluded_incomplete_methods <- setdiff(unique(ref_rows_all$method), unique(ref_rows$method))
      if (nrow(ref_rows) == 0) {
        # Every method present for this (setting,d,rep) is incomplete -- no
        # clean reference exists. Fall back to the unfiltered rows (best
        # available information) rather than crash on ref_methods[1] being
        # NA; this is expected to be unreached in practice (WP7_VERIFICATION
        # .md records exactly one short block across all of WP8's output,
        # and other methods for that same cell are complete).
        ref_rows <- ref_rows_all
      }
      ref_methods <- unique(ref_rows$method)
      ref_one <- ref_rows[ref_rows$method == ref_methods[1], ]
      ref_one <- ref_one[order(ref_one$row_index), ]

      gen_key <- wp7_generator_key(script_id, setting_key, d)
      regen <- tryCatch(regenerate_cell(script_id, env, if (script_id == "87") setting_key else gen_key, d, seed),
                         error = function(e) NULL)

      verify_ok <- TRUE; verify_notes <- character(0)
      if (is.null(regen)) {
        verify_ok <- FALSE; verify_notes <- c(verify_notes, "regeneration failed (see console)")
      } else {
        if (length(regen$true_cluster) != nrow(ref_one) ||
            !identical(as.integer(regen$true_cluster), as.integer(ref_one$true_cluster))) {
          verify_ok <- FALSE
          verify_notes <- c(verify_notes, sprintf("regenerated true_cluster/row-count mismatch: n=%d vs %d",
                                                   length(regen$true_cluster), nrow(ref_one)))
        }
        if (length(ref_methods) > 1) {
          for (mm in ref_methods[-1]) {
            other <- ref_rows[ref_rows$method == mm, ]; other <- other[order(other$row_index), ]
            if (!identical(as.integer(other$true_cluster), as.integer(ref_one$true_cluster))) {
              verify_ok <- FALSE
              verify_notes <- c(verify_notes, sprintf("cross-method true_cluster mismatch (%s vs %s)", mm, ref_methods[1]))
            }
          }
        }
      }
      if (length(excluded_incomplete_methods)) {
        verify_notes <- c(verify_notes, sprintf("excluded incomplete method(s) from reference/cross-check: %s",
                                                 paste(excluded_incomplete_methods, collapse = ",")))
      }
      verify_note <- if (length(verify_notes)) paste(verify_notes, collapse = "; ") else "-"

      if (!verify_ok) {
        bad <- list(n_reg = nrow(ref_one), n_unassigned = 0L, unassigned_frac = 0,
                    true_k = meta$true_k, k_hat = 0L,
                    singleton = list(ari = 0, nmi = 0, ami = 0, status = "error", note = substr(verify_note, 1, 200)),
                    excluded = list(ari = 0, nmi = 0, ami = 0, status = "error", note = substr(verify_note, 1, 200)))
        if (need_dbscan)  wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, "DBSCAN", bad)
        if (need_hdbscan) wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, "HDBSCAN", bad)
        next
      }

      # WP7_VERIFICATION.md change 1 (blocking): X is now the FULL
      # regenerated matrix for 85/86 (n_reg regular rows + n0 outlier rows,
      # regular rows first -- see regenerate_cell()'s comment above). DBSCAN
      # and HDBSCAN are run on the FULL matrix, exactly as the WP8 drivers
      # ran the MCCD methods on it, and the returned labels are subset to
      # seq_len(n_reg) before scoring -- mirroring the drivers' own
      # `res$cluster[seq_len(n_reg)]`. For 87/88, X already has exactly
      # n_reg rows (no contamination in either generator), so the subset is
      # a no-op there.
      X <- regen$X; true_cluster <- regen$true_cluster; row_index <- seq_len(regen$n_reg)
      cont <- if (script_id %in% c("85", "86")) regen$n0 / (regen$n0 + regen$n_reg) else 0

      if (need_dbscan) {
        dlab_full <- tryCatch(run_dbscan_native(X, cont), error = function(e) NULL)
        scored <- if (is.null(dlab_full) || length(dlab_full) != nrow(X)) {
          m <- list(ari = 0, nmi = 0, ami = 0, status = "error", note = "DBSCAN() call failed")
          list(n_reg = length(true_cluster), n_unassigned = 0L, unassigned_frac = 0,
               true_k = meta$true_k, k_hat = 0L, singleton = m, excluded = m)
        } else {
          dlab <- dlab_full[seq_len(regen$n_reg)]
          wp7_score_pair(row_index, true_cluster, dlab)
        }
        wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, "DBSCAN", scored)
      }

      if (need_hdbscan) {
        tag <- paste(script_id, setting_key, d, rep_id, sep = "_")
        hres <- run_hdbscan_python(X, tag)
        scored <- if (identical(hres$status, "error") || length(hres$labels) != nrow(X)) {
          m <- list(ari = 0, nmi = 0, ami = 0, status = "error",
                     note = if (identical(hres$status, "error")) hres$note else
                       sprintf("hdbscan label count %d != nrow(X) %d", length(hres$labels), nrow(X)))
          list(n_reg = length(true_cluster), n_unassigned = 0L, unassigned_frac = 0,
               true_k = meta$true_k, k_hat = 0L, singleton = m, excluded = m)
        } else {
          hlab <- hres$labels[seq_len(regen$n_reg)]
          wp7_score_pair(row_index, true_cluster, hlab)
        }
        wp7_write_cell(out_csv, script_id, setting_key, d, rep_id, seed, "HDBSCAN", scored)
      }
    }
  }
  cat("90_wp7_clustering_quality: run complete (or budget-stopped; rerun to continue).\n")
}

# ---------------------------------------------------------------------------
# 8. --summarize
# ---------------------------------------------------------------------------
do_summarize <- function(out_csv, summ_metrics_csv, khat_csv) {
  if (!file.exists(out_csv)) { cat("no WP7 results yet\n"); return(invisible(NULL)) }
  df <- utils::read.csv(out_csv, stringsAsFactors = FALSE)
  key <- paste(df$script, df$setting_key, df$d, df$rep, df$method, df$unassigned_treatment, sep = "\r")
  df <- df[!duplicated(key, fromLast = TRUE), ]
  df_ok <- df[df$status == "ok", ]
  # WP7_VERIFICATION.md change 5 (non-blocking): unassigned_frac_mean is
  # averaged over ok + degenerate rows, not ok alone -- unassigned_frac is a
  # real, meaningful count for a degenerate cell (e.g. <2 points assigned)
  # even though that cell's ari/nmi/ami are sentinel 0s under status =
  # "degenerate". "incomplete"/"error" rows are still excluded: their
  # unassigned_frac is either not comparable (incomplete: computed from a
  # partial row set) or a bare 0 sentinel (error), not a real observation.
  df_frac <- df[df$status %in% c("ok", "degenerate"), ]

  se <- function(x) if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_

  rows <- list()
  for (grp in split(df_ok, list(df_ok$script, df_ok$setting_key, df_ok$d, df_ok$method, df_ok$unassigned_treatment), drop = TRUE)) {
    frac_grp <- df_frac[df_frac$script == grp$script[1] & df_frac$setting_key == grp$setting_key[1] &
                         df_frac$d == grp$d[1] & df_frac$method == grp$method[1] &
                         df_frac$unassigned_treatment == grp$unassigned_treatment[1], ]
    rows[[length(rows) + 1L]] <- data.frame(
      script = grp$script[1], setting_key = grp$setting_key[1], d = grp$d[1],
      method = grp$method[1], unassigned_treatment = grp$unassigned_treatment[1],
      n_reps = nrow(grp),
      ari_mean = mean(grp$ari), ari_se = se(grp$ari),
      nmi_mean = mean(grp$nmi), nmi_se = se(grp$nmi),
      ami_mean = mean(grp$ami), ami_se = se(grp$ami),
      unassigned_frac_mean = mean(frac_grp$unassigned_frac, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, rows); rownames(S) <- NULL
  if (!is.null(S)) S <- S[order(S$script, S$setting_key, S$d, S$method, S$unassigned_treatment), ]
  write.csv(S, summ_metrics_csv, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", summ_metrics_csv, if (is.null(S)) 0 else nrow(S)))

  # k-hat count table: independent of unassigned_treatment (dedupe first).
  # WP7_VERIFICATION.md change 3 (blocking): filter on status == "ok" --
  # otherwise the sentinel k_hat = 0 written for "error" rows (regeneration/
  # verification failure) and the partial k_hat from "incomplete" rows (a
  # done cell missing cluster rows) enter the k-hat distribution as if they
  # were real observations.
  khat_df <- df[df$unassigned_treatment == "excluded" & df$status == "ok", ]
  krows <- list()
  for (grp in split(khat_df, list(khat_df$script, khat_df$setting_key, khat_df$d, khat_df$method), drop = TRUE)) {
    tab <- table(grp$k_hat)
    for (kv in names(tab)) {
      krows[[length(krows) + 1L]] <- data.frame(
        script = grp$script[1], setting_key = grp$setting_key[1], d = grp$d[1],
        method = grp$method[1], true_k = grp$true_k[1],
        k_hat = as.integer(kv), n_reps_at_this_khat = as.integer(tab[[kv]]),
        stringsAsFactors = FALSE)
    }
  }
  K <- do.call(rbind, krows); rownames(K) <- NULL
  if (!is.null(K)) K <- K[order(K$script, K$setting_key, K$d, K$method, K$k_hat), ]
  write.csv(K, khat_csv, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", khat_csv, if (is.null(K)) 0 else nrow(K)))
  invisible(list(summary = S, khat = K))
}

# ---------------------------------------------------------------------------
# 9. CLI
# ---------------------------------------------------------------------------
if (!isTRUE(getOption("wp7.no_main", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  MODE_SMOKE <- "--smoke" %in% args
  MODE_SUMM  <- "--summarize" %in% args
  scripts_arg <- args[grepl("^--scripts=", args)]
  scripts_sel <- if (length(scripts_arg)) strsplit(sub("^--scripts=", "", scripts_arg[1]), ",")[[1]] else names(CLUSTER_FILE)
  args <- args[!grepl("^--", args)]

  WP8_DIR       <- here::here("revision_experiments/results/tr1/wp8")
  WP8_SMOKE_DIR <- file.path(WP8_DIR, "smoke")
  OUT_DIR       <- here::here("revision_experiments/results/tr1/wp7")
  OUT_SMOKE_DIR <- file.path(OUT_DIR, "smoke")

  in_dir  <- if (MODE_SMOKE) WP8_SMOKE_DIR else (if (length(args) >= 1 && nzchar(args[1])) args[1] else WP8_DIR)
  out_dir <- if (MODE_SMOKE) OUT_SMOKE_DIR else OUT_DIR
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv         <- file.path(out_dir, "wp7_metrics_long.csv")
  summ_metrics_csv <- file.path(out_dir, "wp7_summary_metrics.csv")
  khat_csv        <- file.path(out_dir, "wp7_khat_table.csv")

  if (MODE_SUMM) {
    do_summarize(out_csv, summ_metrics_csv, khat_csv)
  } else {
    budget <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(args[2]) else 1800
    cat(sprintf("90_wp7_clustering_quality: input=%s output=%s scripts=%s budget=%gs\n",
                in_dir, out_dir, paste(scripts_sel, collapse = ","), budget))
    do_run(in_dir, out_csv, budget, scripts_sel, MODE_SMOKE)
  }
}
