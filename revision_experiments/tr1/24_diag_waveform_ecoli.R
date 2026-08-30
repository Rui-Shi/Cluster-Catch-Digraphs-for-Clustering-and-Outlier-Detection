#!/usr/bin/env Rscript
# revision_experiments/tr1/24_diag_waveform_ecoli.R
#
# Agent D diagnostic script for the NEUCOM-D-26-15191 real-data regeneration
# (see REGEN_SPEC.md). DIAGNOSIS ONLY -- this script never modifies harness.R,
# wp0_mccd_methods.R, 13_wp0_gate.R, published_realdata_truth.csv, or anything
# under methods/, R/, data/, simulations/. It writes only to its own CSV shard,
# results/diag_waveform_ecoli.csv.
#
# Investigates why waveform and ecoli fail to reproduce Section 6 published
# TPR/TNR/BA/F2 for all four proposed detectors (U-MCCD, SU-MCCD, UN-MCCD,
# SUN-MCCD):
#   waveform (n=3443, d=21, n0=100): all four under-flag TPR badly.
#   ecoli    (n=336,  d=7,  n0=8):   all four are TPR +1/8 = +0.125 vs published.
#
# Usage (run ONE task per invocation -- some waveform tasks run 5-10 minutes):
#   Rscript 24_diag_waveform_ecoli.R meta
#   Rscript 24_diag_waveform_ecoli.R waveform_umccd_default
#   Rscript 24_diag_waveform_ecoli.R waveform_unmccd_default
#   Rscript 24_diag_waveform_ecoli.R waveform_unmccd_q99
#   Rscript 24_diag_waveform_ecoli.R waveform_sunmccd_q99
#   Rscript 24_diag_waveform_ecoli.R ecoli_full

suppressMessages(library(here))
source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))

OUT_CSV   <- here::here("revision_experiments/results/tr1/diag_waveform_ecoli.csv")
TRUTH_CSV <- here::here("revision_experiments/tr1/published_realdata_truth.csv")
truth <- read.csv(TRUTH_CSV, stringsAsFactors = FALSE)

get_truth <- function(dataset, method) {
  v <- setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2"))
  if (is.na(method)) return(v)
  sub <- truth[truth$dataset == dataset & truth$method == method, ]
  if (nrow(sub) == 0) return(v)
  for (i in seq_len(nrow(sub))) if (sub$metric[i] %in% names(v)) v[sub$metric[i]] <- sub$value[i]
  v
}

with_timeout <- function(fn, timeout_sec) {
  setTimeLimit(cpu = Inf, elapsed = timeout_sec, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  fn()
}

ROW_COLS <- c("dataset", "method", "hypothesis", "variant", "quant_label", "min_cls",
              "n", "d", "n0", "TPR", "TNR", "BA", "F2",
              "published_TPR", "published_TNR", "published_BA", "published_F2",
              "max_abs_diff", "match_3dp", "note", "t_total", "status", "timestamp")

mk_row <- function(dataset, method = NA_character_, hypothesis, variant, quant_label = NA_character_,
                    min_cls = NA_real_, n = NA_integer_, d = NA_integer_, n0 = NA_integer_,
                    TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
                    note = NA_character_, t_total = NA_real_, status = "ok") {
  pub <- get_truth(dataset, method)
  vals <- c(TPR = TPR, TNR = TNR, BA = BA, F2 = F2)
  if (all(!is.na(vals)) && all(!is.na(pub))) {
    max_diff <- max(abs(vals - pub))
    match3 <- all(round(vals, 3) == round(pub, 3))
  } else {
    max_diff <- NA_real_; match3 <- NA
  }
  data.frame(
    dataset = dataset, method = method, hypothesis = hypothesis, variant = variant,
    quant_label = quant_label, min_cls = min_cls, n = n, d = d, n0 = n0,
    TPR = TPR, TNR = TNR, BA = BA, F2 = F2,
    published_TPR = unname(pub["TPR"]), published_TNR = unname(pub["TNR"]),
    published_BA = unname(pub["BA"]), published_F2 = unname(pub["F2"]),
    max_abs_diff = max_diff, match_3dp = match3, note = note, t_total = t_total, status = status,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE
  )[ROW_COLS]
}

write_row <- function(row) append_result(OUT_CSV, row)

args <- commandArgs(trailingOnly = TRUE)
TASK <- if (length(args) >= 1) args[1] else stop("Usage: Rscript 24_diag_waveform_ecoli.R <task>")
cat(sprintf("=== 24_diag_waveform_ecoli.R: task = %s ===\n", TASK))

# ---------------------------------------------------------------------------
if (TASK == "meta") {

  for (ds in c("waveform", "ecoli")) {
    dat <- load_real_dataset(ds)
    n0 <- sum(dat$Y == 0)
    cat(sprintf("[meta] %-9s n=%d d=%d n0=%d\n", ds, dat$n, dat$d, n0))
  }

  # ---- H3 (waveform): quantile-table existence / metadata, and probe whether
  # the RK/NN tables needed for a stricter alpha at d=21 exist at all ----
  cat("\n--- H3: quantile table inspection (waveform, d=21) ---\n")
  inspect_table <- function(variant, d, q) {
    fname <- sprintf("%s-test-simul_%dd_%s%%.RData",
                      ifelse(variant == "RK", "RK", "NN"), d, q)
    dir <- ifelse(variant == "RK", "R/RK-test_quantile", "R/NN-test_quantile")
    path <- here::here(dir, fname)
    exists_f <- file.exists(path)
    cat(sprintf("[H3] %s : exists=%s\n", path, exists_f))
    if (!exists_f) return(invisible(NULL))
    e <- new.env(); load(path, envir = e)
    cat(sprintf("     objects: %s\n", paste(ls(e), collapse = ", ")))
    simul <- get("simul", envir = e)
    cat(sprintf("     names(simul) = %s\n", paste(names(simul), collapse = ", ")))
    if (!is.null(simul$quan)) {
      for (qn in names(simul$quan)) {
        qt <- simul$quan[[qn]]
        cat(sprintf("     quan[['%s']]: dim = %s\n", qn, paste(dim(qt), collapse = " x ")))
      }
    }
    if (!is.null(simul$average)) cat(sprintf("     length(average)=%d length(median)=%d\n", length(simul$average), length(simul$median)))
    if (!is.null(simul$r)) cat(sprintf("     length(r)=%d\n", length(simul$r)))
    invisible(simul)
  }
  inspect_table("RK", 21, "999")
  inspect_table("NN", 21, "999")
  inspect_table("NN", 21, "99")

  cat("\n--- H1 pre-check: do stricter-alpha tables even exist at d=21? ---\n")
  probe_quant <- function(variant, d, q) {
    res <- tryCatch({ get_simul(variant, d, quant = q); "EXISTS" },
                     error = function(e) conditionMessage(e))
    cat(sprintf("[probe] get_simul(%s, d=%d, quant=%s) -> %s\n", variant, d, q, res))
    res
  }
  rk99  <- probe_quant("RK", 21, "99")
  rk95  <- probe_quant("RK", 21, "95")
  nn99  <- probe_quant("NN", 21, "99")   # expected to EXIST
  nn95  <- probe_quant("NN", 21, "95")   # expected to be missing
  write_row(mk_row("waveform", "U-MCCD/SU-MCCD", "H1-alpha-sweep", "rk-table-probe-99", quant_label = "99",
                    note = paste0("get_simul(RK,21,99) -> ", rk99), status = ifelse(rk99 == "EXISTS", "ok", "missing_table")))
  write_row(mk_row("waveform", "U-MCCD/SU-MCCD", "H1-alpha-sweep", "rk-table-probe-95", quant_label = "95",
                    note = paste0("get_simul(RK,21,95) -> ", rk95), status = ifelse(rk95 == "EXISTS", "ok", "missing_table")))
  write_row(mk_row("waveform", "UN-MCCD/SUN-MCCD", "H1-alpha-sweep", "nn-table-probe-99", quant_label = "99",
                    note = paste0("get_simul(NN,21,99) -> ", nn99), status = ifelse(nn99 == "EXISTS", "ok", "missing_table")))
  write_row(mk_row("waveform", "UN-MCCD/SUN-MCCD", "H1-alpha-sweep", "nn-table-probe-95", quant_label = "95",
                    note = paste0("get_simul(NN,21,95) -> ", nn95), status = ifelse(nn95 == "EXISTS", "ok", "missing_table")))

  # ---- H4 (waveform): preprocessing / MADN scaling ----
  cat("\n--- H4: waveform preprocessing (MADN scale_R) ---\n")
  dat_wf <- load_real_dataset("waveform")
  Xs <- as.matrix(dat_wf$X)  # already scale_R()'d inside RealData_Collection.R
  col_sd  <- apply(Xs, 2, sd)
  col_min <- apply(Xs, 2, min)
  col_max <- apply(Xs, 2, max)
  cat("[H4] per-column sd (post MADN-scale):\n");  print(round(col_sd, 4))
  cat("[H4] per-column min (post scale):\n");      print(round(col_min, 4))
  cat("[H4] per-column max (post scale):\n");      print(round(col_max, 4))
  cat(sprintf("[H4] any NA/NaN/Inf in scaled X: NA=%s NaN=%s Inf=%s\n",
              anyNA(Xs), any(is.nan(Xs)), any(is.infinite(Xs))))
  cat(sprintf("[H4] sd ratio max/min across the 21 columns: %.4f\n", max(col_sd) / min(col_sd)))

  raw <- foreign::read.arff(here::here("data/outlier_detection/Waveform_withoutdupl_v01.arff"))
  raw$id <- NULL
  raw_X <- raw[, setdiff(names(raw), "outlier"), drop = FALSE]
  raw_mad <- sapply(raw_X, mad)
  cat("[H4] raw (pre-scale) per-column MAD:\n"); print(round(raw_mad, 6))
  cat(sprintf("[H4] any raw column MAD == 0: %s ; min raw MAD = %.6f\n", any(raw_mad == 0), min(raw_mad)))
  note_h4 <- sprintf("sd_ratio=%.4f; min_raw_mad=%.6f; any_raw_mad_zero=%s; any_scaled_inf_or_nan=%s",
                      max(col_sd) / min(col_sd), min(raw_mad), any(raw_mad == 0),
                      any(is.nan(Xs)) || any(is.infinite(Xs)))
  write_row(mk_row("waveform", NA_character_, "H4-preprocessing", "madn-scale-check",
                    n = dat_wf$n, d = dat_wf$d, n0 = sum(dat_wf$Y == 0), note = note_h4, status = "ok"))

  # ---- H3 (ecoli): label definition / source classes ----
  cat("\n--- H3: ecoli label definition ---\n")
  dat_ec <- load_real_dataset("ecoli")
  n0_ec <- sum(dat_ec$Y == 0)
  cat(sprintf("[H3-ecoli] n0 = %d (loader)\n", n0_ec))
  ecoli_csv <- read.csv(here::here("data/outlier_detection/ecoli.csv"), stringsAsFactors = FALSE)
  cat("[H3-ecoli] ecoli.csv 'class' column distribution:\n"); print(table(ecoli_csv$class))
  mat_path <- here::here("data/outlier_detection/ecoli.mat")
  note_mat <- "ecoli.mat not present"
  if (file.exists(mat_path)) {
    m <- R.matlab::readMat(mat_path)
    cat(sprintf("[H3-ecoli] ecoli.mat top-level names: %s\n", paste(names(m), collapse = ", ")))
    for (nm in names(m)) {
      obj <- m[[nm]]
      cat(sprintf("   %s: class=%s dim=%s\n", nm, class(obj)[1], paste(dim(obj), collapse = "x")))
      if (is.numeric(obj) && length(obj) < 400) {
        cat(sprintf("     unique values: %s\n", paste(sort(unique(as.vector(obj))), collapse = ", ")))
      }
    }
    note_mat <- sprintf("ecoli.mat names: %s", paste(names(m), collapse = ","))
  } else {
    cat("[H3-ecoli] ecoli.mat not found -- cannot cross-check against a richer multi-class label.\n")
  }
  note_h3 <- sprintf("n0=%d; csv class table: %s; %s", n0_ec,
                      paste(capture.output(print(table(ecoli_csv$class))), collapse = " / "), note_mat)
  write_row(mk_row("ecoli", NA_character_, "H3-label-definition", "class-source-check",
                    n = dat_ec$n, d = dat_ec$d, n0 = n0_ec, note = note_h3, status = "ok"))

  cat("\n=== meta task done ===\n")

# ---------------------------------------------------------------------------
} else if (TASK == "waveform_umccd_default") {

  dat <- load_real_dataset("waveform")
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  n0 <- sum(Y == 0)
  cat(sprintf("[waveform] n=%d d=%d n0=%d\n", n, d, n0))

  cell <- tryCatch({
    with_timeout(function() {
      t0 <- Sys.time()
      res <- umccd_method(X = X, d = d, Y = Y)  # default quant -> rk_quant_label_paper(21) = "999"
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[["U-MCCD"]])
      radii <- res$radii[[1]]
      cat(sprintf("[U-MCCD default] quant_label=%s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f wall=%.1fs\n",
                  res$quant_label, m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
      cat(sprintf("[H2-RK radii] n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f unassigned_rows=%d\n",
                  length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows))
      note <- sprintf("RK radii: n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f; unassigned_rows=%d",
                       length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows)
      write_row(mk_row("waveform", "U-MCCD", "H2-RK-radii", "default-quant999", quant_label = res$quant_label,
                        n = n, d = d, n0 = n0, TPR = unname(m["TPR"]), TNR = unname(m["TNR"]),
                        BA = unname(m["BA"]), F2 = unname(m["F2"]), note = note, t_total = wall, status = "ok"))
    }, 560)
  }, error = function(e) {
    cat(sprintf("[U-MCCD default] ERROR/TIMEOUT: %s\n", conditionMessage(e)))
    write_row(mk_row("waveform", "U-MCCD", "H2-RK-radii", "default-quant999",
                      n = n, d = d, n0 = n0, note = conditionMessage(e), status = "error_or_timeout"))
  })

# ---------------------------------------------------------------------------
} else if (TASK == "waveform_unmccd_default") {

  dat <- load_real_dataset("waveform")
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  n0 <- sum(Y == 0)
  cat(sprintf("[waveform] n=%d d=%d n0=%d\n", n, d, n0))

  cell <- tryCatch({
    with_timeout(function() {
      t0 <- Sys.time()
      res <- unmccd_method(X = X, d = d, Y = Y)  # default -> nn_quant_label_paper_UN(21) = "999"
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[["UN-MCCD"]])
      radii <- res$radii[[1]]
      cat(sprintf("[UN-MCCD default] quant_label=%s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f wall=%.1fs\n",
                  res$quant_label, m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
      cat(sprintf("[H2-NND radii] n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f unassigned_rows=%d\n",
                  length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows))
      note <- sprintf("NND radii: n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f; unassigned_rows=%d",
                       length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows)
      write_row(mk_row("waveform", "UN-MCCD", "H2-NND-radii", "default-quant999", quant_label = res$quant_label,
                        n = n, d = d, n0 = n0, TPR = unname(m["TPR"]), TNR = unname(m["TNR"]),
                        BA = unname(m["BA"]), F2 = unname(m["F2"]), note = note, t_total = wall, status = "ok"))
    }, 560)
  }, error = function(e) {
    cat(sprintf("[UN-MCCD default] ERROR/TIMEOUT: %s\n", conditionMessage(e)))
    write_row(mk_row("waveform", "UN-MCCD", "H2-NND-radii", "default-quant999",
                      n = n, d = d, n0 = n0, note = conditionMessage(e), status = "error_or_timeout"))
  })

# ---------------------------------------------------------------------------
} else if (TASK == "waveform_unmccd_q99") {

  dat <- load_real_dataset("waveform")
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  n0 <- sum(Y == 0)
  cat(sprintf("[waveform] n=%d d=%d n0=%d\n", n, d, n0))

  cell <- tryCatch({
    with_timeout(function() {
      t0 <- Sys.time()
      res <- unmccd_method(X = X, d = d, Y = Y, quant = "99")  # forced stricter alpha
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[["UN-MCCD"]])
      radii <- res$radii[[1]]
      cat(sprintf("[UN-MCCD quant=99] TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f wall=%.1fs\n",
                  m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
      cat(sprintf("[H1/H2 radii @99] n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f unassigned_rows=%d\n",
                  length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows))
      note <- sprintf("forced quant=99; radii zero_frac=%.4f min=%.6f median=%.6f max=%.6f; unassigned_rows=%d",
                       mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows)
      write_row(mk_row("waveform", "UN-MCCD", "H1-alpha-sweep", "forced-quant99", quant_label = "99",
                        n = n, d = d, n0 = n0, TPR = unname(m["TPR"]), TNR = unname(m["TNR"]),
                        BA = unname(m["BA"]), F2 = unname(m["F2"]), note = note, t_total = wall, status = "ok"))
    }, 560)
  }, error = function(e) {
    cat(sprintf("[UN-MCCD quant=99] ERROR/TIMEOUT: %s\n", conditionMessage(e)))
    write_row(mk_row("waveform", "UN-MCCD", "H1-alpha-sweep", "forced-quant99", quant_label = "99",
                      n = n, d = d, n0 = n0, note = conditionMessage(e), status = "error_or_timeout"))
  })

# ---------------------------------------------------------------------------
} else if (TASK == "waveform_sunmccd_q99") {

  dat <- load_real_dataset("waveform")
  X <- dat$X; Y <- dat$Y; d <- dat$d; n <- dat$n
  n0 <- sum(Y == 0)
  min_cls_val <- 0.5 * (n0 / n)  # half_contam reading, same convention as wp0_gate_v3.csv
  cat(sprintf("[waveform] n=%d d=%d n0=%d min.cls(half_contam)=%.6f\n", n, d, n0, min_cls_val))

  cell <- tryCatch({
    with_timeout(function() {
      t0 <- Sys.time()
      res <- sunmccd_method(X = X, d = d, Y = Y, quant = "99", min.cls = min_cls_val)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[["SUN-MCCD"]])
      radii <- res$radii[[1]]
      cat(sprintf("[SUN-MCCD quant=99] TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f wall=%.1fs\n",
                  m["TPR"], m["TNR"], m["BA"], m["F2"], wall))
      cat(sprintf("[H1/H2 radii @99] n=%d zero_frac=%.4f min=%.6f median=%.6f max=%.6f unassigned_rows=%d\n",
                  length(radii), mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows))
      note <- sprintf("forced quant=99, min.cls=half_contam=%.6f; radii zero_frac=%.4f min=%.6f median=%.6f max=%.6f; unassigned_rows=%d",
                       min_cls_val, mean(radii == 0), min(radii), median(radii), max(radii), res$unassigned_rows)
      write_row(mk_row("waveform", "SUN-MCCD", "H1-alpha-sweep", "forced-quant99", quant_label = "99",
                        min_cls = min_cls_val, n = n, d = d, n0 = n0,
                        TPR = unname(m["TPR"]), TNR = unname(m["TNR"]), BA = unname(m["BA"]), F2 = unname(m["F2"]),
                        note = note, t_total = wall, status = "ok"))
    }, 560)
  }, error = function(e) {
    cat(sprintf("[SUN-MCCD quant=99] ERROR/TIMEOUT: %s\n", conditionMessage(e)))
    write_row(mk_row("waveform", "SUN-MCCD", "H1-alpha-sweep", "forced-quant99", quant_label = "99",
                      min_cls = min_cls_val, n = n, d = d, n0 = n0,
                      note = conditionMessage(e), status = "error_or_timeout"))
  })

# ---------------------------------------------------------------------------
} else if (TASK == "ecoli_full") {

  dat <- load_real_dataset("ecoli")
  X0 <- as.matrix(dat$X); Y <- dat$Y; d <- dat$d; n <- dat$n
  n0 <- sum(Y == 0)
  outlier_idx <- which(Y == 0)
  cat(sprintf("[ecoli] n=%d d=%d n0=%d outlier_idx=%s\n", n, d, n0, paste(outlier_idx, collapse = ",")))

  madn_scale   <- function(x) { M <- median(x); s <- mad(x); (x - M) / s }
  zscore_scale <- function(x) as.numeric(scale(x))
  variants <- list(
    "as-supplied" = X0,
    "madn"        = apply(X0, 2, madn_scale),
    "zscore"      = apply(X0, 2, zscore_scale)
  )

  min_cls_val <- 0.5 * (n0 / n)  # half_contam, matches wp0_gate_v3.csv convention

  methods_list <- list(
    "U-MCCD"   = function(X) umccd_method(X = X, d = d, Y = Y),
    "SU-MCCD"  = function(X) sumccd_method(X = X, d = d, Y = Y, min.cls = min_cls_val),
    "UN-MCCD"  = function(X) unmccd_method(X = X, d = d, Y = Y),
    "SUN-MCCD" = function(X) sunmccd_method(X = X, d = d, Y = Y, min.cls = min_cls_val)
  )

  flagged_by_method_asis <- list()

  for (vname in names(variants)) {
    Xv <- variants[[vname]]
    for (mname in names(methods_list)) {
      t0 <- Sys.time()
      res <- tryCatch(methods_list[[mname]](Xv), error = function(e) e)
      wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (inherits(res, "error")) {
        cat(sprintf("[ecoli %s / %-9s] ERROR: %s\n", vname, mname, conditionMessage(res)))
        write_row(mk_row("ecoli", mname, "H2-preprocessing", vname, n = n, d = d, n0 = n0,
                          note = paste("error:", conditionMessage(res)), t_total = wall, status = "error"))
        next
      }
      thr <- REAL_DATA_THRESHOLDS[[mname]]
      m <- evaluate(Y, res$score, thr)
      flagged_outliers <- outlier_idx[res$score[outlier_idx] == 1]
      missed_outliers  <- outlier_idx[res$score[outlier_idx] == 0]
      cat(sprintf("[ecoli %-11s / %-9s] TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f | flagged=%s | missed=%s | wall=%.2fs\n",
                  vname, mname, m["TPR"], m["TNR"], m["BA"], m["F2"],
                  paste(flagged_outliers, collapse = ","), paste(missed_outliers, collapse = ","), wall))
      if (vname == "as-supplied") {
        flagged_by_method_asis[[mname]] <- list(flagged = flagged_outliers, missed = missed_outliers)
      }
      note <- sprintf("flagged_outlier_idx=[%s]; missed_outlier_idx=[%s]; unassigned_rows=%d",
                       paste(flagged_outliers, collapse = ","), paste(missed_outliers, collapse = ","), res$unassigned_rows)
      write_row(mk_row("ecoli", mname, "H1-H2-ecoli", vname, n = n, d = d, n0 = n0,
                        min_cls = if (mname %in% c("SU-MCCD", "SUN-MCCD")) min_cls_val else NA_real_,
                        TPR = unname(m["TPR"]), TNR = unname(m["TNR"]), BA = unname(m["BA"]), F2 = unname(m["F2"]),
                        note = note, t_total = wall, status = "ok"))
    }
  }

  cat("\n[H1 summary] missed true-outlier index by method (as-supplied preprocessing):\n")
  for (mname in names(flagged_by_method_asis)) {
    cat(sprintf("  %-9s missed=%s\n", mname, paste(flagged_by_method_asis[[mname]]$missed, collapse = ",")))
  }
  all_missed_sorted <- lapply(flagged_by_method_asis, function(x) sort(x$missed))
  same_across <- length(unique(all_missed_sorted)) == 1
  cat(sprintf("[H1] identical missed-outlier index across all 4 methods: %s\n", same_across))
  if (same_across && length(all_missed_sorted[[1]]) == 1) {
    idx <- all_missed_sorted[[1]]
    cat(sprintf("[H1] the single common missed outlier is row index %d; feature values:\n", idx))
    print(X0[idx, , drop = FALSE])
  }
  write_row(mk_row("ecoli", NA_character_, "H1-summary", "missed-index-consistency",
                    n = n, d = d, n0 = n0,
                    note = sprintf("missed_by_method=%s; identical_across_methods=%s",
                                    paste(sprintf("%s:[%s]", names(all_missed_sorted),
                                                  sapply(all_missed_sorted, paste, collapse = ",")), collapse = "; "),
                                    same_across),
                    status = "ok"))

} else {
  stop(sprintf("Unknown task: %s", TASK))
}

cat("\n=== 24_diag_waveform_ecoli.R: task complete ===\n")
