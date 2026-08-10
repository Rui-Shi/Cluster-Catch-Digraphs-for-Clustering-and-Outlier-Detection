#!/usr/bin/env Rscript
# 50_wp2a_nn_alpha_bigsets.R -- WP2(a): NND-arm alpha sensitivity on the two
# real data sets whose NN quantile tables are genuinely distinct at two
# tokens, including SUN-MCCD, the paper's headline method.
#
# WHY THESE TWO DATASETS ONLY. Every other real data set (see
# 42_wp2a_alpha_sweep.R's --nncmp step) sits at a dimension where the NN
# quantile table generator produced two byte-identical files under different
# filenames -- comparing "99%" vs "99.9%" there compares a file with itself,
# so any "alpha does not matter" finding at those d would be an artifact of
# duplicate inputs, not a property of the test. Only two real-data dimensions
# carry genuinely different NN tables:
#   wilt        n=4819 d=5   tokens "90" / "95"
#   pageblocks  n=4795 d=10  tokens "99" / "999"
# Grid: 2 datasets x {UN-MCCD, SUN-MCCD} x 2 tokens each = 8 cells.
#
# THE KNOB. wp0_mccd_methods.R's unmccd_method/sunmccd_method both take
# `quant` as an override: when non-NULL it is treated as the exact NN
# filename-token string ("90"/"95"/"99"/"999") and handed straight to
# get_simul("NN", d, quant = <token>), which is what selects which table
# unmccd/sunmccd's own internals see (UNMCCD_outlier/SUNMCCD_outlier never
# receive `quant` as a bare number at all -- they take `simul` directly, see
# wp0_mccd_methods.R L316-348). So sweeping alpha here needs nothing beyond
# `METHOD_REGISTRY[["UN-MCCD"]](X=, d=, Y=, quant=tok)` -- no RK-style
# derivation shadow is needed (unlike 42_wp2a_alpha_sweep.R's RK arm), because
# both tokens at both dimensions are already on disk.
#
# get_simul SHADOW (installed AFTER every source(), same discipline as
# 42/43/45 -- sourcing wp0_mccd_methods.R defines closures whose free
# reference to `get_simul` is resolved by CURRENT globalenv binding at call
# time, so an override installed before sourcing would be clobbered, and one
# installed after is honoured by every already-defined wrapper). Here the
# shadow does no derivation; it is a thin provenance tap that records which
# file path and file size get_simul actually returned for the call the
# detector just made, so every output row states exactly which .RData file
# produced it (columns table_file, table_size_bytes).
#
# DISTINCT-TABLE GUARD (run before ANY cell, in both --guard mode and the
# real grid). This is exactly the failure this sweep exists to rule out: if
# either dataset's two on-disk tokens turned out to be duplicates after all,
# every "alpha doesn't matter" cell here would silently repeat the same
# mistake documented in 42's comments. Both files are loaded directly (not
# through get_simul, so this check is independent of the shadow above) and
# compared with identical() on $average and $median; the script stop()s
# loudly if any pair matches.
#
# Usage:
#   Rscript 50_wp2a_nn_alpha_bigsets.R --guard
#     Runs only the distinct-table guard (no cells), prints the comparison,
#     exits. Cheap (loads 4 small .RData files, <1s).
#
#   Rscript 50_wp2a_nn_alpha_bigsets.R --smoke
#     Runs the exact same per-cell code path (run_cell()) against a small
#     SYNTHETIC matrix (n=150, d=5) instead of loading wilt/pageblocks, using
#     the real NN d=5 tokens ("90","95") and UN-MCCD, writing to a separate
#     results/tr1/wp2a_nn_alpha_bigsets_smoke.csv so the production file is
#     never touched. Proves the full plumbing (get_simul shadow, method call,
#     evaluate(), CSV row, resume) cheaply. Safe to run twice to see resume
#     skip the second time.
#
#   Rscript 50_wp2a_nn_alpha_bigsets.R [datasets] [methods] [alpha_tokens] [min_cls]
#     The real grid. All four positional args are optional:
#       datasets     comma list, default "wilt,pageblocks"
#       methods      comma list, default "UN-MCCD,SUN-MCCD"
#       alpha_tokens comma list, or "ALL" (default) to auto-discover every
#                    token on disk for that dataset's dimension (the
#                    intended production use -- "ALL" resolves to exactly
#                    {"90","95"} for wilt and {"99","999"} for pageblocks;
#                    an explicit list is an advanced override applied
#                    UNIFORMLY to every dataset in the run, for debugging
#                    only -- a token missing on disk for some dataset's d
#                    is a hard error via get_simul(), not silently skipped)
#       min_cls      numeric, proportion of n, default 0.0625 (manuscript
#                    value); passed to SUN-MCCD only -- UN-MCCD's wrapper
#                    takes no min.cls argument at all, so it is never passed
#                    to it, per the task brief.
#
# Checkpointed on (dataset, method, alpha_token): a restart skips completed
# cells via harness.R's has_result(). Per-cell timeout: 2400s, recorded as
# status="timeout" with the measured elapsed time -- never silently dropped.
# CSV is appended and flushed (write.csv/write.table, not a held-open
# connection) after every single cell, so a killed process loses at most the
# cell in flight.

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# --- get_simul provenance shadow: installed AFTER every source() -----------
orig_get_simul <- get_simul   # capture BEFORE shadowing

PROV <- new.env(parent = emptyenv())
PROV$last_file  <- NA_character_
PROV$last_size  <- NA_real_
PROV$last_quant <- NA_character_

get_simul_shadow <- function(variant = c("RK", "NN"), d, quant = NULL) {
  res <- orig_get_simul(variant, d, quant)
  PROV$last_file  <- res$file
  PROV$last_size  <- suppressWarnings(file.size(res$file))
  PROV$last_quant <- res$quant_label
  res
}
assign("get_simul", get_simul_shadow, envir = globalenv())

# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_GUARD <- any(args == "--guard")
MODE_SMOKE <- any(args == "--smoke")
args <- args[!args %in% c("--guard", "--smoke")]

CELL_TIMEOUT <- 2400   # seconds; per task spec
NN_METHODS   <- c("UN-MCCD", "SUN-MCCD")
MIN_CLS_METHODS <- c("SUN-MCCD")   # UN-MCCD's wrapper takes no min.cls arg

OUT_CSV   <- here::here("revision_experiments/results/tr1/wp2a_nn_alpha_bigsets.csv")
SMOKE_CSV <- here::here("revision_experiments/results/tr1/wp2a_nn_alpha_bigsets_smoke.csv")

# The two real-data dimensions this sweep exists for, and what is on disk.
BIGSET_D <- list(wilt = 5L, pageblocks = 10L)

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
token_of_q   <- function(q) sub("^0\\.", "", format(q, trim = TRUE, scientific = FALSE))
q_of_token   <- function(tok) as.numeric(paste0("0.", tok))
alpha_of_q   <- function(q) 1 - q

disk_tokens <- function(d) {
  pat <- sprintf("^NN-test-simul_%dd_([0-9]+)%%\\.RData$", d)
  f <- list.files(NN_QUANT_TABLE_DIR, pattern = pat)
  toks <- sub(pat, "\\1", f)
  toks[order(q_of_token(toks))]
}

# ---------------------------------------------------------------------------
# Distinct-table guard: load every on-disk token's raw table directly
# (bypassing get_simul/the shadow entirely) and confirm no two are identical.
# ---------------------------------------------------------------------------
run_guard <- function(datasets = names(BIGSET_D)) {
  cat("==== 50_wp2a: distinct-table guard ====\n")
  ok <- TRUE
  for (ds in datasets) {
    d <- BIGSET_D[[ds]]
    toks <- disk_tokens(d)
    cat(sprintf("  %-12s d=%-3d tokens on disk: %s\n", ds, d, paste(toks, collapse = ", ")))
    if (length(toks) < 2) {
      cat(sprintf("    *** ABORT: fewer than 2 tokens on disk for d=%d -- this dataset cannot support an alpha sweep\n", d))
      ok <- FALSE
      next
    }
    tabs <- lapply(toks, function(tk) {
      e <- new.env()
      load(file.path(NN_QUANT_TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tk)), envir = e)
      get("simul", envir = e)
    })
    for (i in 1:(length(toks) - 1)) for (j in (i + 1):length(toks)) {
      same_avg <- identical(tabs[[i]]$average, tabs[[j]]$average)
      same_med <- identical(tabs[[i]]$median,  tabs[[j]]$median)
      max_d_avg <- max(abs(tabs[[i]]$average - tabs[[j]]$average))
      cat(sprintf("    %s (%s) vs %s (%s): identical(average)=%-5s identical(median)=%-5s max|d avg|=%.4g\n",
                  toks[i], format(file.size(file.path(NN_QUANT_TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, toks[i]))), big.mark = ","),
                  toks[j], format(file.size(file.path(NN_QUANT_TABLE_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, toks[j]))), big.mark = ","),
                  same_avg, same_med, max_d_avg))
      if (same_avg && same_med) {
        cat(sprintf("    *** GUARD FAILED: %s and %s are byte-identical at d=%d -- aborting, this sweep would be comparing a file with itself\n",
                    toks[i], toks[j], d))
        ok <- FALSE
      }
    }
  }
  if (!ok) stop("50_wp2a: distinct-table guard FAILED -- see log above. Refusing to run the grid.")
  cat("==== 50_wp2a: distinct-table guard PASSED (all on-disk token pairs are genuinely different) ====\n\n")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Core cell: generic over (X, Y, d) so the exact same code path serves both
# the real grid (real datasets) and --smoke (synthetic matrix).
# ---------------------------------------------------------------------------
run_cell <- function(dataset_label, X, Y, d, method, alpha_token, min_cls, out_csv) {
  keys <- c(dataset = dataset_label, method = method, alpha_token = alpha_token)
  if (isTRUE(has_result(out_csv, keys))) {
    cat(sprintf("[skip, already recorded] %s x %s x alpha_token=%s\n", dataset_label, method, alpha_token))
    return(invisible("skip"))
  }

  PROV$last_file <- NA_character_; PROV$last_size <- NA_real_; PROV$last_quant <- NA_character_
  n <- nrow(X)
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
    res <- if (method %in% MIN_CLS_METHODS) {
      METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, quant = alpha_token, min.cls = min_cls)
    } else {
      METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, quant = alpha_token)
    }
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    if (length(res$score) != n || anyNA(res$score)) {
      stop(sprintf("score sanity check failed: length=%d (n=%d), anyNA=%s", length(res$score), n, anyNA(res$score)))
    }
    m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[[method]])
    list(m = m, res = res, status = "ok", note = NA_character_)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")), res = NULL,
         status = if (grepl("elapsed time limit|reached elapsed", msg)) "timeout" else "error",
         note = msg)
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  flagged <- if (is.null(out$res)) integer(0) else which(out$res$score > 0.5)
  ncls <- if (is.null(out$res)) NA_integer_ else
    length(unique(out$res$cluster[!is.na(out$res$cluster)]))

  row <- list(
    dataset = dataset_label, method = method,
    alpha_token = alpha_token, alpha = alpha_of_q(q_of_token(alpha_token)),
    n = n, d = d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
    n_flagged  = if (is.null(out$res)) NA_integer_ else length(flagged),
    n_clusters = ncls,
    flagged_idx = collapse_idx(flagged),
    table_file = PROV$last_file, table_size_bytes = PROV$last_size,
    elapsed_sec = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time())
  )
  append_result(out_csv, row)

  note_1line <- if (is.na(out$note)) "" else substr(gsub("[\r\n]+", " ", out$note), 1, 150)
  fields <- sprintf("dataset=%s method=%s param=alpha_token:%s status=%s sec=%.1f",
                     dataset_label, method, alpha_token, out$status, wall)
  if (identical(out$status, "ok")) {
    cat(sprintf("CELL_DONE %s BA=%.3f F2=%.3f nflag=%s ncls=%s table=%s\n",
                fields,
                if (is.na(out$m[["BA"]])) NA_real_ else out$m[["BA"]],
                if (is.na(out$m[["F2"]])) NA_real_ else out$m[["F2"]],
                if (is.null(out$res)) "NA" else length(flagged),
                ifelse(is.na(ncls), "NA", ncls),
                basename(ifelse(is.na(PROV$last_file), "NA", PROV$last_file))))
  } else {
    cat(sprintf("CELL_FAIL %s note=%s\n", fields, note_1line))
  }
  flush.console()
  invisible(out$status)
}

# ---------------------------------------------------------------------------
# --guard mode
# ---------------------------------------------------------------------------
if (MODE_GUARD) {
  run_guard()
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# --smoke mode: exact same run_cell() path, synthetic n=150 d=5 matrix, the
# real on-disk d=5 tokens ("90","95"), method UN-MCCD, separate output CSV.
# ---------------------------------------------------------------------------
if (MODE_SMOKE) {
  run_guard(datasets = "wilt")   # d=5 guard, cheap, also exercised here
  set.seed(42)
  n_s <- 150L; d_s <- 5L
  X_s <- matrix(rnorm(n_s * d_s), nrow = n_s, ncol = d_s)
  # push the last 10 rows away from the bulk so there is a genuine minority
  # component for the mutual-catch-graph translation to find.
  X_s[(n_s - 9):n_s, ] <- X_s[(n_s - 9):n_s, ] + 6
  Y_s <- c(rep(1L, n_s - 10L), rep(0L, 10L))   # 1 = regular, 0 = outlier

  cat("==== 50_wp2a --smoke: synthetic n=150, d=5, method=UN-MCCD, tokens 90 & 95 ====\n")
  for (tok in disk_tokens(d_s)) {
    run_cell("SMOKE_SYNTH_D5", X_s, Y_s, d_s, "UN-MCCD", tok, min_cls = NA_real_, out_csv = SMOKE_CSV)
  }
  cat("\n==== 50_wp2a --smoke: re-invoking the SAME cells to demonstrate resume/skip ====\n")
  for (tok in disk_tokens(d_s)) {
    run_cell("SMOKE_SYNTH_D5", X_s, Y_s, d_s, "UN-MCCD", tok, min_cls = NA_real_, out_csv = SMOKE_CSV)
  }
  cat("ALL_CELLS_COMPLETE\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Real grid
# ---------------------------------------------------------------------------
resolve_list <- function(s, default) if (is.null(s) || !nzchar(s)) default else strsplit(s, ",")[[1]]

DATASETS     <- resolve_list(if (length(args) >= 1) args[1] else NULL, names(BIGSET_D))
METHODS      <- resolve_list(if (length(args) >= 2) args[2] else NULL, NN_METHODS)
ALPHA_ARG    <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "ALL"
MIN_CLS      <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 0.0625

stopifnot(all(DATASETS %in% names(BIGSET_D)))
stopifnot(all(METHODS %in% NN_METHODS))

cat(sprintf("50_wp2a_nn_alpha_bigsets: datasets = %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("50_wp2a_nn_alpha_bigsets: methods  = %s\n", paste(METHODS, collapse = ", ")))
cat(sprintf("50_wp2a_nn_alpha_bigsets: alpha_tokens arg = %s ; min_cls = %g ; out = %s\n",
            ALPHA_ARG, MIN_CLS, OUT_CSV))

run_guard(DATASETS)

for (ds in DATASETS) {
  d <- BIGSET_D[[ds]]
  toks <- if (identical(toupper(ALPHA_ARG), "ALL")) disk_tokens(d) else strsplit(ALPHA_ARG, ",")[[1]]
  cat(sprintf("\n== %s  d=%d  alpha tokens = %s ==\n", ds, d, paste(toks, collapse = ", ")))
  dat <- load_real_dataset(ds)
  stopifnot(dat$d == d)   # loader's own dimension must agree with BIGSET_D

  for (meth in METHODS) {
    for (tok in toks) {
      run_cell(ds, dat$X, dat$Y, d, meth, tok, min_cls = MIN_CLS, out_csv = OUT_CSV)
    }
  }
}

cat("ALL_CELLS_COMPLETE\n")
