#!/usr/bin/env Rscript
# 42_wp2a_alpha_sweep.R -- WP2(a): sensitivity of the four MCCD detectors to
# alpha, the significance level of the spatial-randomness test, on the 16 real
# data sets.
#
#   Rscript 42_wp2a_alpha_sweep.R --inventory          # Steps 1 & 2
#   Rscript 42_wp2a_alpha_sweep.R [datasets] [methods] # Step 3 (grid)
#   Rscript 42_wp2a_alpha_sweep.R --summary            # Step 4
#   Rscript 42_wp2a_alpha_sweep.R --cost [d] [niter]   # Step 5 (NN gen cost)
#
# WHY.  Reviewers R1.4, R3.1, R3.2, R5.2 and AE.2 all object that alpha in the
# submitted manuscript varies both with d AND with method, with no stated
# justification. This script measures what alpha actually buys.
#
# THE KNOB.  alpha enters only through which Monte-Carlo quantile table the
# CCD construction consults. wp0_mccd_methods.R's four wrappers already expose
# it as `quant`, a FILENAME TOKEN string ("90", "95", "99", "999"), which they
# hand to harness.R's get_simul(). U-MCCD / SU-MCCD read the RK tables,
# UN-MCCD / SUN-MCCD read the NN tables. No override is needed to sweep the
# tokens that exist on disk.
#
# THE ASYMMETRY (established empirically by --inventory; see the findings md).
#   RK: the saved object is list(Kest.m =, quan =, r =) -- Kest.m is the FULL
#       niter x (m*rn) matrix of raw Monte-Carlo draws (Kest.R:412). The
#       reduced envelope in $quan is just apply(Kest.m, 2, quantile, probs=q)
#       reshaped to nrow = m (Kest.R:406-411). So every alpha is recoverable
#       from a file already on disk at zero simulation cost. That is why the
#       RK files are 280-620 MB and the NN files are 70 KB.
#   NN: NNDestP.simpois.lower.quant returns list(average =, median =) only
#       (R/ccds/NN_Dist_Est.R:111) -- the per-iteration draws are discarded
#       inside the generator. A new alpha therefore requires a new Monte-Carlo
#       run. The NN half of this sweep is limited to the tokens on disk.
#
# HOW THE RK DERIVATION IS WIRED.  get_simul() is shadowed by a global binding
# installed AFTER every source() (harness.R and wp0_mccd_methods.R both define
# into globalenv, so an override installed earlier would be silently clobbered
# -- the same trap documented in 40_wp2a_tolerance_sweep.R). When
# WP2A$rk_derive is TRUE the shadow loads ONE base file for that d, recomputes
# the $quan entry for the requested probability straight from $Kest.m, and
# hands the detector a table whose $quan key matches the numeric quantile it
# will index with. Deriving every alpha at a given d from the SAME base file
# is deliberate: it holds the Monte-Carlo replicate fixed so the only thing
# moving across the sweep is alpha. Where a second RK file exists on disk for
# the same d, that file is ALSO run at its own token (alpha_source = "disk"),
# which gives a free independent-replicate check on the derivation.
#
# DERIVATION GUARD.  Before any derived cell runs, rk_derive_quan() is asked to
# reproduce the base file's OWN stored $quan entry from its $Kest.m. If the
# recomputation is not identical to the stored matrix, the script aborts. A
# silent mismatch would make every "derived" row a different estimator rather
# than a different alpha.
#
# Checkpointed on (dataset, method, alpha_token, alpha_source): a restart
# skips completed cells, so a PowerShell 600 s timeout costs at most one cell.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- only now is it safe to define and install the get_simul shadow ---------

args <- commandArgs(trailingOnly = TRUE)
MODE_INVENTORY <- any(args == "--inventory")
MODE_SUMMARY   <- any(args == "--summary")
MODE_COST      <- any(args == "--cost")
MODE_NNCMP     <- any(args == "--nncmp")
args <- args[!args %in% c("--inventory", "--summary", "--cost", "--nncmp")]

S_MIN           <- 0.0625   # proportion of n, not a count (manuscript value)
MIN_CLS_METHODS <- c("SU-MCCD", "SUN-MCCD")
CELL_TIMEOUT    <- 900      # seconds; a cell over this is recorded as timed out

OUT_CSV      <- here::here("revision_experiments/results/tr1/wp2a_alpha_sweep.csv")
SUMMARY_CSV  <- here::here("revision_experiments/results/tr1/wp2a_alpha_summary.csv")
FINDINGS_MD  <- here::here("revision_experiments/results/tr1/wp2a_alpha_findings.md")

# The 12 distinct dimensions carried by the 16 real data sets.
REAL_D <- c(5, 6, 7, 8, 9, 10, 12, 16, 18, 19, 21, 30)

# Dataset order: the eight small sets first (they are the affordable grid),
# then the eight large ones in increasing n.
SMALL_SETS <- c("hepatitis", "lymphography", "glass", "WBC",
                "vertebral", "ecoli", "stamps", "WDBC")
LARGE_SETS <- c("pima", "shuffle", "vowels", "PenDigits",
                "waveform", "thyroid", "pageblocks", "wilt")
ALL_SETS   <- c(SMALL_SETS, LARGE_SETS)

RK_METHODS <- c("U-MCCD", "SU-MCCD")
NN_METHODS <- c("UN-MCCD", "SUN-MCCD")

# alpha levels to DERIVE from the RK raw draws (probabilities, not alphas).
RK_DERIVED_Q <- c(0.90, 0.95, 0.99, 0.999)

token_of_q <- function(q) sub("^0\\.", "", format(q, trim = TRUE, scientific = FALSE))
q_of_token <- function(tok) as.numeric(paste0("0.", tok))
alpha_of_q <- function(q) 1 - q

# ---------------------------------------------------------------------------
# Step 1: what is on disk
# ---------------------------------------------------------------------------
disk_tokens <- function(variant, d) {
  dir <- if (variant == "RK") RK_QUANT_TABLE_DIR else NN_QUANT_TABLE_DIR
  pat <- sprintf("^%s-test-simul_%dd_([0-9]+)%%\\.RData$", variant, d)
  f <- list.files(dir, pattern = pat)
  toks <- sub(pat, "\\1", f)
  toks[order(q_of_token(toks))]
}

inventory_table <- function() {
  rows <- lapply(REAL_D, function(d) data.frame(
    d = d,
    RK_tokens = paste(disk_tokens("RK", d), collapse = " "),
    RK_n = length(disk_tokens("RK", d)),
    NN_tokens = paste(disk_tokens("NN", d), collapse = " "),
    NN_n = length(disk_tokens("NN", d)),
    stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Step 2: what is actually inside a stored `simul`
# ---------------------------------------------------------------------------
inspect_one <- function(path) {
  e <- new.env(); load(path, envir = e)
  cat(sprintf("\n----- %s  (%.1f MB on disk) -----\n", basename(path),
              file.size(path) / 1024^2))
  cat("  objects in file :", paste(ls(e), collapse = ", "), "\n")
  if (!exists("simul", envir = e)) { cat("  NO OBJECT NAMED simul\n"); return(invisible(NULL)) }
  s <- get("simul", envir = e)
  cat("  class(simul)    :", paste(class(s), collapse = "/"), "\n")
  cat("  names(simul)    :", paste(names(s), collapse = ", "), "\n")
  for (nm in names(s)) {
    x <- s[[nm]]
    cat(sprintf("    $%-8s %-9s dim=%-14s len=%-8s\n", nm, class(x)[1],
                if (is.null(dim(x))) "-" else paste(dim(x), collapse = "x"),
                length(x)))
    if (is.list(x)) cat("               keys:", paste(names(x), collapse = ", "), "\n")
  }
  invisible(s)
}

# ---------------------------------------------------------------------------
# NN table diff: do the two NN tables at a given d actually differ numerically?
#
# The grid finds the NN detectors returning identical labels at both alpha
# levels for every data set. Two very different explanations: (i) the envelopes
# differ but the decision does not -- genuine saturation of the test; (ii) the
# two files hold (nearly) the same numbers, in which case "alpha" was never
# varied at all. This separates them.
# ---------------------------------------------------------------------------
nn_table_compare <- function() {
  cat("\n==== NN quantile vectors: do the two tables at each d differ? ====\n")
  for (dd in REAL_D) {
    toks <- disk_tokens("NN", dd)
    if (length(toks) < 2) { cat(sprintf("  d=%-3d only one table (%s)\n", dd, toks)); next }
    g <- lapply(toks, function(tk) {
      e <- new.env()
      load(file.path(NN_QUANT_TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", dd, tk)), envir = e)
      get("simul", envir = e)
    })
    a1 <- g[[1]]$average; a2 <- g[[2]]$average
    m1 <- g[[1]]$median;  m2 <- g[[2]]$median
    cat(sprintf("  d=%-3d %s vs %s :: identical(avg)=%-5s max|dAvg|=%.4g max rel=%.4g ; identical(med)=%-5s max|dMed|=%.4g ; mean(a1)=%.4g mean(a2)=%.4g\n",
                dd, toks[1], toks[2], identical(a1, a2), max(abs(a1 - a2)),
                max(abs(a1 - a2) / pmax(abs(a1), 1e-30)), identical(m1, m2),
                max(abs(m1 - m2)), mean(a1), mean(a2)))
  }
}

# ---------------------------------------------------------------------------
# RK derivation: rebuild the $quan entry for an arbitrary probability from the
# stored raw draws, exactly as Kest.R:406-411 does it.
#   temp <- apply(Kest.m, 2, quantile, probs = q);  matrix(temp, nrow = m)
# Chunked over columns so the transient copy apply() makes stays small (the
# d<10 tables are 2000 x 50000 doubles, ~800 MB resident).
# ---------------------------------------------------------------------------
rk_derive_quan <- function(Kest.m, m, q, chunk = 5000L) {
  nc <- ncol(Kest.m)
  out <- numeric(nc)
  i <- 1L
  while (i <= nc) {
    j <- min(i + chunk - 1L, nc)
    out[i:j] <- apply(Kest.m[, i:j, drop = FALSE], 2, quantile, probs = q)
    i <- j + 1L
  }
  matrix(out, nrow = m)
}

# Control environment for the shadow (an environment, not bare globals, so the
# detectors' recursive source() calls cannot stomp it).
WP2A <- new.env(parent = emptyenv())
WP2A$rk_derive     <- FALSE
WP2A$rk_base_label <- NULL   # token of the base file to derive FROM
WP2A$rk_target_q   <- NULL   # probability to derive
WP2A$rk_cache_d    <- NA_integer_
WP2A$rk_cache      <- NULL   # the one loaded base simul, with derived keys added
WP2A$rk_cache_file <- NA_character_
WP2A$guard_done    <- new.env(parent = emptyenv())

orig_get_simul <- get_simul   # capture BEFORE shadowing

rk_base_table <- function(d, base_label) {
  if (!identical(WP2A$rk_cache_d, d) || !identical(WP2A$rk_cache_label, base_label)) {
    t0 <- Sys.time()
    base <- orig_get_simul("RK", d, quant = base_label)
    WP2A$rk_cache       <- base$simul
    WP2A$rk_cache_d     <- d
    WP2A$rk_cache_label <- base_label
    WP2A$rk_cache_file  <- base$file
    cat(sprintf("    [rk base] loaded %s in %.1fs\n", basename(base$file),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    # ---- DERIVATION GUARD: reproduce the file's own stored envelope --------
    s <- WP2A$rk_cache
    if (is.null(s$Kest.m)) stop("rk_base_table(): no $Kest.m in ", base$file,
                                " -- RK alpha cannot be derived from this file")
    own_q   <- base$quant
    own_key <- as.character(own_q)
    stored  <- s$quan[[own_key]]
    if (is.null(stored)) stop("rk_base_table(): stored $quan has no key ", own_key)
    m <- nrow(stored)
    redone <- rk_derive_quan(s$Kest.m, m, own_q)
    if (!identical(dim(redone), dim(stored)) || !isTRUE(all.equal(redone, stored, tolerance = 0))) {
      stop(sprintf("rk_base_table(): DERIVATION GUARD FAILED for %s -- recomputed $quan[['%s']] does not reproduce the stored matrix",
                   base$file, own_key))
    }
    cat(sprintf("    [rk base] derivation guard OK (reproduced stored $quan[['%s']], %dx%d)\n",
                own_key, nrow(stored), ncol(stored)))
    WP2A$rk_cache_m <- m
  }
  WP2A$rk_cache
}

get_simul_shadow <- function(variant = c("RK", "NN"), d, quant = NULL) {
  variant <- match.arg(variant)
  if (variant == "RK" && isTRUE(WP2A$rk_derive)) {
    s   <- rk_base_table(d, WP2A$rk_base_label)
    q   <- WP2A$rk_target_q
    key <- as.character(q)
    if (is.null(s$quan[[key]])) {
      t0 <- Sys.time()
      s$quan[[key]] <- rk_derive_quan(s$Kest.m, WP2A$rk_cache_m, q)
      WP2A$rk_cache <- s
      cat(sprintf("    [rk derive] q=%s from %s in %.1fs\n", key,
                  basename(WP2A$rk_cache_file),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    return(list(simul = s, quant = q, quant_label = token_of_q(q),
                file = WP2A$rk_cache_file))
  }
  orig_get_simul(variant, d, quant)
}
assign("get_simul", get_simul_shadow, envir = globalenv())

# ---------------------------------------------------------------------------
# The per-dataset cell plan
# ---------------------------------------------------------------------------
# For an RK method at dimension d:
#   base file = the paper's default token (rk_quant_label_paper: d<10 -> "99",
#   d>=10 -> "999"); all of RK_DERIVED_Q are derived from it (source
#   "derived"), and every OTHER token on disk at that d is additionally run
#   from its own file (source "disk") as an independent-replicate check.
# For an NN method: every token on disk at that d, source "disk". Nothing can
#   be derived.
cell_plan <- function(d, method) {
  if (method %in% RK_METHODS) {
    base <- rk_quant_label_paper(d)
    if (!base %in% disk_tokens("RK", d)) base <- disk_tokens("RK", d)[1]
    derived <- data.frame(alpha_token = token_of_q(RK_DERIVED_Q),
                          alpha_q = RK_DERIVED_Q, alpha_source = "derived",
                          base_token = base, stringsAsFactors = FALSE)
    others <- setdiff(disk_tokens("RK", d), base)
    disk <- if (length(others)) data.frame(
      alpha_token = others, alpha_q = q_of_token(others),
      alpha_source = "disk", base_token = others, stringsAsFactors = FALSE)
      else NULL
    rbind(derived, disk)
  } else {
    toks <- disk_tokens("NN", d)
    data.frame(alpha_token = toks, alpha_q = q_of_token(toks),
               alpha_source = "disk", base_token = toks, stringsAsFactors = FALSE)
  }
}

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
parse_idx    <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else as.integer(strsplit(s, ";")[[1]])

# ---------------------------------------------------------------------------
# Step 3: the grid
# ---------------------------------------------------------------------------
run_grid <- function(datasets, methods) {
  cat(sprintf("42_wp2a: datasets = %s\n", paste(datasets, collapse = ", ")))
  cat(sprintf("42_wp2a: methods  = %s\n", paste(methods, collapse = ", ")))
  cat(sprintf("42_wp2a: S_min    = %g   out = %s\n", S_MIN, OUT_CSV))

  for (ds in datasets) {
    dat <- load_real_dataset(ds)
    d <- dat$d
    cat(sprintf("\n== %s  n=%d d=%d n_out=%d ==\n", ds, dat$n, d, sum(dat$Y == 0)))

    # RK methods first while the (large) base table for this d is hot.
    for (meth in intersect(c(RK_METHODS, NN_METHODS), methods)) {
      plan <- cell_plan(d, meth)
      for (k in seq_len(nrow(plan))) {
        tok <- plan$alpha_token[k]; qq <- plan$alpha_q[k]
        src <- plan$alpha_source[k]; base <- plan$base_token[k]
        keys <- c(dataset = ds, method = meth, alpha_token = tok, alpha_source = src)
        if (isTRUE(has_result(OUT_CSV, keys))) {
          cat(sprintf("  [skip] %s x %s x a=%s (%s)\n", ds, meth, tok, src)); next
        }

        if (meth %in% RK_METHODS && src == "derived") {
          WP2A$rk_derive <- TRUE; WP2A$rk_base_label <- base; WP2A$rk_target_q <- qq
        } else {
          WP2A$rk_derive <- FALSE
        }

        t0 <- Sys.time()
        out <- tryCatch({
          setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
          res <- if (meth %in% MIN_CLS_METHODS) {
            METHOD_REGISTRY[[meth]](X = dat$X, d = d, Y = dat$Y, quant = tok, min.cls = S_MIN)
          } else {
            METHOD_REGISTRY[[meth]](X = dat$X, d = d, Y = dat$Y, quant = tok)
          }
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[meth]]),
               res = res, status = "ok", note = NA_character_)
        }, error = function(e) {
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          msg <- conditionMessage(e)
          list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
               res = NULL,
               status = if (grepl("elapsed time limit|reached elapsed", msg)) "timeout" else "error",
               note = msg)
        })
        wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        flagged <- if (is.null(out$res)) integer(0) else which(out$res$score > 0.5)
        ncls <- if (is.null(out$res)) NA_integer_ else
          length(unique(out$res$cluster[!is.na(out$res$cluster)]))

        append_result(OUT_CSV, list(
          dataset = ds, method = meth,
          alpha_token = tok, alpha_q = qq, alpha = alpha_of_q(qq),
          alpha_source = src, base_token = base,
          variant = if (meth %in% RK_METHODS) "RK" else "NN",
          n = dat$n, d = d, n_outliers = sum(dat$Y == 0),
          TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
          BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
          n_flagged  = if (is.null(out$res)) NA_integer_ else length(flagged),
          n_clusters = ncls,
          flagged_idx = collapse_idx(flagged),
          elapsed_sec = wall,
          s_min = if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_,
          status = out$status, note = out$note, timestamp = format(Sys.time())
        ))

        cat(sprintf("  %-12s %-9s a=%-5s(%-7s) BA=%s F2=%s nflag=%-5s ncls=%-3s %-7s %.1fs\n",
                    ds, meth, tok, src,
                    if (is.na(out$m[["BA"]])) "  NA " else sprintf("%.3f", out$m[["BA"]]),
                    if (is.na(out$m[["F2"]])) "  NA " else sprintf("%.3f", out$m[["F2"]]),
                    if (is.null(out$res)) "NA" else length(flagged),
                    if (is.na(ncls)) "NA" else ncls, out$status, wall))
        flush.console()
      }
    }
    WP2A$rk_derive <- FALSE
    WP2A$rk_cache <- NULL; WP2A$rk_cache_d <- NA_integer_
    WP2A$rk_cache_label <- NULL
    invisible(gc(FALSE))
  }
  cat("\n42_wp2a: grid done\n")
}

# Record cells we deliberately did not run, with the cost that justified it.
mark_not_run <- function(datasets, methods, reason) {
  for (ds in datasets) {
    d <- switch(ds, hepatitis = 19, lymphography = 18, glass = 9, WBC = 9,
                vertebral = 6, ecoli = 7, stamps = 9, WDBC = 30, pima = 8,
                shuffle = 9, vowels = 12, PenDigits = 16, waveform = 21,
                thyroid = 6, pageblocks = 10, wilt = 5)
    nn <- switch(ds, hepatitis = 74, lymphography = 148, glass = 213, WBC = 223,
                 vertebral = 240, ecoli = 336, stamps = 340, WDBC = 367,
                 pima = 555, shuffle = 1013, vowels = 1452, PenDigits = 3200,
                 waveform = 3443, thyroid = 3656, pageblocks = 4795, wilt = 4819)
    for (meth in methods) {
      plan <- cell_plan(d, meth)
      for (k in seq_len(nrow(plan))) {
        keys <- c(dataset = ds, method = meth, alpha_token = plan$alpha_token[k],
                  alpha_source = plan$alpha_source[k])
        if (isTRUE(has_result(OUT_CSV, keys))) next
        append_result(OUT_CSV, list(
          dataset = ds, method = meth, alpha_token = plan$alpha_token[k],
          alpha_q = plan$alpha_q[k], alpha = alpha_of_q(plan$alpha_q[k]),
          alpha_source = plan$alpha_source[k], base_token = plan$base_token[k],
          variant = if (meth %in% RK_METHODS) "RK" else "NN",
          n = nn, d = d, n_outliers = NA_integer_,
          TPR = NA_real_, TNR = NA_real_, BA = NA_real_, F2 = NA_real_,
          n_flagged = NA_integer_, n_clusters = NA_integer_, flagged_idx = "",
          elapsed_sec = NA_real_,
          s_min = if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_,
          status = "not_run", note = reason, timestamp = format(Sys.time())))
      }
    }
  }
  cat("42_wp2a: not_run rows written\n")
}

# ---------------------------------------------------------------------------
# Step 4: analysis
# ---------------------------------------------------------------------------
build_summary <- function() {
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE,
                 colClasses = c(alpha_token = "character", base_token = "character",
                                flagged_idx = "character"))
  df$flagged_idx[is.na(df$flagged_idx)] <- ""
  ok <- df[df$status == "ok", ]

  rows <- list()
  for (ds in unique(ok$dataset)) for (meth in unique(ok$method[ok$dataset == ds])) {
    for (srcgrp in unique(ok$alpha_source[ok$dataset == ds & ok$method == meth])) {
      sub <- ok[ok$dataset == ds & ok$method == meth & ok$alpha_source == srcgrp, ]
      sub <- sub[order(sub$alpha_q), ]
      if (nrow(sub) < 1) next
      sets <- lapply(sub$flagged_idx, parse_idx)
      # pairwise across the sorted alpha ladder
      n_levels <- nrow(sub)
      all_identical <- TRUE; min_jac <- 1
      if (n_levels > 1) for (i in 1:(n_levels - 1)) for (j in (i + 1):n_levels) {
        a <- sets[[i]]; b <- sets[[j]]
        if (!identical(a, b)) all_identical <- FALSE
        u <- length(union(a, b))
        min_jac <- min(min_jac, if (u == 0) 1 else length(intersect(a, b)) / u)
      }
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, method = meth, variant = sub$variant[1],
        alpha_source = srcgrp, n = sub$n[1], d = sub$d[1],
        n_levels = n_levels,
        alpha_tokens = paste(sub$alpha_token, collapse = ";"),
        labels_identical_across_alpha = all_identical,
        min_pairwise_jaccard = min_jac,
        n_flagged_min = min(sub$n_flagged), n_flagged_max = max(sub$n_flagged),
        BA_min = min(sub$BA), BA_max = max(sub$BA), BA_range = max(sub$BA) - min(sub$BA),
        F2_min = min(sub$F2), F2_max = max(sub$F2), F2_range = max(sub$F2) - min(sub$F2),
        TPR_range = max(sub$TPR) - min(sub$TPR), TNR_range = max(sub$TNR) - min(sub$TNR),
        BA_by_alpha = paste(sprintf("%s=%.4f", sub$alpha_token, sub$BA), collapse = ";"),
        nflag_by_alpha = paste(sprintf("%s=%d", sub$alpha_token, sub$n_flagged), collapse = ";"),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$d, out$dataset, out$method, out$alpha_source), ]
  dir.create(dirname(SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("42_wp2a: wrote %s (%d rows)\n", SUMMARY_CSV, nrow(out)))

  cat("\n==== CELL COUNTS ====\n")
  print(table(df$status))
  cat(sprintf("  ok cells: %d over %d datasets\n", nrow(ok), length(unique(ok$dataset))))

  cat("\n==== DOES THE LABEL SET MOVE WITH alpha? (>=2 levels only) ====\n")
  m <- out[out$n_levels > 1, ]
  cat(sprintf("  (dataset,method,source) groups with >=2 alpha levels: %d\n", nrow(m)))
  cat(sprintf("  groups where the flagged set is IDENTICAL at every alpha: %d\n",
              sum(m$labels_identical_across_alpha)))
  cat(sprintf("  groups where it CHANGES: %d\n", sum(!m$labels_identical_across_alpha)))
  mv <- m[!m$labels_identical_across_alpha, ]
  if (nrow(mv)) {
    mv <- mv[order(-mv$BA_range), ]
    for (i in seq_len(nrow(mv))) cat(sprintf(
      "    %-12s %-9s d=%-3d %-8s jac_min=%.3f  BA %.3f..%.3f (range %.3f)  nflag %d..%d  [%s]\n",
      mv$dataset[i], mv$method[i], mv$d[i], mv$alpha_source[i], mv$min_pairwise_jaccard[i],
      mv$BA_min[i], mv$BA_max[i], mv$BA_range[i], mv$n_flagged_min[i], mv$n_flagged_max[i],
      mv$BA_by_alpha[i]))
  }

  cat("\n==== LARGEST MOVES ====\n")
  if (nrow(m)) {
    k <- which.max(m$BA_range)
    cat(sprintf("  largest BA range: %.4f  (%s / %s / %s, d=%d, %s)\n",
                m$BA_range[k], m$dataset[k], m$method[k], m$alpha_source[k], m$d[k],
                m$BA_by_alpha[k]))
    k2 <- which.max(m$F2_range)
    cat(sprintf("  largest F2 range: %.4f  (%s / %s / %s, d=%d)\n",
                m$F2_range[k2], m$dataset[k2], m$method[k2], m$alpha_source[k2], m$d[k2]))
  }

  cat("\n==== SATURATION: does the flagged set move across the alpha ladder? ====\n")
  cat("  tables_distinct = FALSE means the two NN files compared are byte-identical\n")
  cat("  duplicates, so identical output is guaranteed by the input and is NOT\n")
  cat("  evidence that alpha stopped binding. Only d=5 and d=10 hold genuinely\n")
  cat("  different NN tables (see --nncmp).\n")
  sat <- m[, c("dataset", "method", "variant", "d", "alpha_source",
               "labels_identical_across_alpha", "min_pairwise_jaccard", "n_levels")]
  sat <- sat[order(sat$variant, sat$d), ]
  sat$tables_distinct <- mapply(function(v, dd) {
    if (v == "RK") return(TRUE)
    toks <- disk_tokens("NN", dd)
    if (length(toks) < 2) return(NA)
    g <- lapply(toks, function(tk) {
      e <- new.env()
      load(file.path(NN_QUANT_TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", dd, tk)), envir = e)
      get("simul", envir = e)$average
    })
    !identical(g[[1]], g[[2]])
  }, sat$variant, sat$d)

  for (i in seq_len(nrow(sat))) cat(sprintf(
    "  %-3s d=%-3d %-12s %-9s %-8s levels=%d tables_distinct=%-5s  %s\n",
    sat$variant[i], sat$d[i], sat$dataset[i], sat$method[i], sat$alpha_source[i],
    sat$n_levels[i], sat$tables_distinct[i],
    if (sat$labels_identical_across_alpha[i]) "labels IDENTICAL" else
      sprintf("labels move (jac_min=%.3f)", sat$min_pairwise_jaccard[i])))

  for (v in c("NN", "RK")) {
    sv <- sat[sat$variant == v & !is.na(sat$tables_distinct) & sat$tables_distinct, ]
    if (!nrow(sv)) { cat(sprintf("\n  %s: no group has genuinely distinct tables at two alpha levels\n", v)); next }
    movers <- sv$d[!sv$labels_identical_across_alpha]
    satd   <- sv$d[sv$labels_identical_across_alpha]
    cat(sprintf("\n  %s (INFORMATIVE groups only -- tables actually differ): alpha binds at d in {%s}; labels identical at d in {%s}\n",
                v, paste(sort(unique(movers)), collapse = ","),
                paste(sort(unique(satd)), collapse = ",")))
    if (length(movers) && length(satd) && min(satd) > max(movers))
      cat(sprintf("  %s: saturation sets in at d = %d (no informative group at d >= %d moves)\n",
                  v, min(satd), min(satd)))
  }
  nun <- sum(!is.na(sat$tables_distinct) & !sat$tables_distinct)
  cat(sprintf("\n  UNINFORMATIVE groups (duplicate NN tables, identical by construction): %d of %d\n",
              nun, nrow(sat)))

  cat("\n==== DERIVED vs DISK at the same alpha (RK MC-replicate check) ====\n")
  rk <- ok[ok$variant == "RK", ]
  for (ds in unique(rk$dataset)) for (meth in unique(rk$method[rk$dataset == ds])) {
    s <- rk[rk$dataset == ds & rk$method == meth, ]
    for (tok in unique(s$alpha_token)) {
      t2 <- s[s$alpha_token == tok, ]
      if (length(unique(t2$alpha_source)) < 2) next
      a <- parse_idx(t2$flagged_idx[t2$alpha_source == "derived"][1])
      b <- parse_idx(t2$flagged_idx[t2$alpha_source == "disk"][1])
      u <- length(union(a, b))
      cat(sprintf("  %-12s %-9s a=%-4s derived vs disk: %s (jac=%.3f, nflag %d vs %d)\n",
                  ds, meth, tok, if (identical(a, b)) "IDENTICAL" else "DIFFERS",
                  if (u == 0) 1 else length(intersect(a, b)) / u, length(a), length(b)))
    }
  }

  bad <- df[!df$status %in% c("ok"), ]
  cat(sprintf("\n==== NON-OK CELLS: %d ====\n", nrow(bad)))
  if (nrow(bad)) print(table(bad$status, bad$dataset))
  invisible(out)
}

# ---------------------------------------------------------------------------
# Step 5: NN generation cost probe (measure only; never launch production)
# ---------------------------------------------------------------------------
run_cost <- function(d = 5L, niter = 200L, cores = 4L) {
  outdir <- here::here("revision_experiments/results/tr1/cost_probe")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  script <- here::here("revision_experiments/01_gen_quantile_table.R")
  cat(sprintf("42_wp2a --cost: NN d=%d niter=%d cores=%d -> %s\n", d, niter, cores, outdir))
  t0 <- Sys.time()
  st <- system2(file.path(R.home("bin"), "Rscript"),
                args = c(shQuote(script), "NN", d, "0.95", niter, cores, shQuote(outdir)),
                stdout = TRUE, stderr = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(paste(st, collapse = "\n"), "\n")
  cat(sprintf("\n42_wp2a --cost: wall = %.1f s for niter=%d at n=%d(default PROD_N) cores=%d\n",
              el, niter, 1000L, cores))
  cat(sprintf("  per-iteration = %.4f s (contended: other agents are running)\n", el / niter))
  cat(sprintf("  projected niter=10000 (production) = %.1f s = %.2f h\n",
              el / niter * 10000, el / niter * 10000 / 3600))
  invisible(el)
}

# ---------------------------------------------------------------------------
if (MODE_INVENTORY) {
  cat("==== STEP 1: quantile tokens on disk, per variant x real-data dimension ====\n\n")
  inv <- inventory_table()
  print(inv, row.names = FALSE)
  cat("\n==== STEP 2: what is inside the stored `simul` object ====\n")
  inspect_one(file.path(RK_QUANT_TABLE_DIR, "RK-test-simul_9d_99%.RData"))
  inspect_one(file.path(NN_QUANT_TABLE_DIR, "NN-test-simul_9d_95%.RData"))
  inspect_one(file.path(NN_QUANT_TABLE_DIR, "NN-test-simul_9d_99%.RData"))
  nn_table_compare()
} else if (MODE_NNCMP) {
  nn_table_compare()
} else if (MODE_SUMMARY) {
  build_summary()
} else if (MODE_COST) {
  dd <- if (length(args) >= 1) as.integer(args[1]) else 5L
  ni <- if (length(args) >= 2) as.integer(args[2]) else 200L
  run_cost(dd, ni, 4L)
} else {
  DATASETS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else SMALL_SETS
  METHODS  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else
                c(RK_METHODS, NN_METHODS)
  if (length(args) >= 3 && args[3] == "not_run") {
    mark_not_run(DATASETS, METHODS, if (length(args) >= 4) args[4] else "not run: cost")
  } else {
    run_grid(DATASETS, METHODS)
  }
}
