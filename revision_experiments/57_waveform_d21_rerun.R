#!/usr/bin/env Rscript
# 57_waveform_d21_rerun.R -- waveform (n=3443, d=21, 100 outliers) rerun on
# UN-MCCD and SUN-MCCD against a FRESH, genuinely-distinct d=21 NN quantile
# pair at 99% / 99.9%, plus a shipped-table baseline for comparison against
# what the manuscript currently reports.
#
# WHY THIS EXISTS. The manuscript claimed the NND-based spatial-randomness
# test "ceases to discriminate" at d=21, evidenced by waveform returning
# identical output at the 99% and 99.9% NN quantile tables. Both shipped
# files (R/NN-test_quantile/NN-test-simul_21d_99%.RData and
# ..._999%.RData) are byte-identical -- 21 such duplicate pairs exist on
# disk (d = 6-9, 11-19, 21-28; see 56_verify_d21_duplication.R) -- so that
# observation compared a file with itself. The claim has been withdrawn.
# A genuine pair is being generated at d=21, n=5000, niter=2000 and will
# land in R/NN-test_quantile_d21_regen/samedraws/ as
# NN-test-simul_21d_99%.RData and NN-test-simul_21d_999%.RData. That
# generation takes ~16h and is NOT done as this file is written -- this
# script only prepares the rerun; see --dryrun below for validation that
# does not require the fresh files to exist yet.
#
# GRID (6 cells, dataset fixed to waveform):
#   {UN-MCCD, SUN-MCCD} x {fresh:99, fresh:999, shipped:999}
# The shipped baseline uses only token "999" (not also "99") because the
# shipped 99%/999% pair is asserted identical by guard 3 below -- running
# both would be redundant, and "999" is what the manuscript's published
# waveform rows for UN-MCCD/SUN-MCCD actually used (wp0_gate_v2.csv,
# regen_final_waveform.csv: quant_label="999" for both methods at d=21,
# via nn_quant_label_paper_UN(21)/nn_quant_label_paper_SUN(21), both of
# which resolve to "999" for d>=20 / d>=10 respectively).
# SUN-MCCD gets min.cls = 0.0625 (the manuscript's S_min); UN-MCCD's
# wrapper (unmccd_method(), wp0_mccd_methods.R:316) takes no min.cls
# argument at all, so none is passed to it, per task instruction.
#
# TABLE SELECTION / get_simul SHADOW. harness.R's get_simul(variant, d,
# quant) (~line 283) resolves NN tables from the single fixed directory
# NN_QUANT_TABLE_DIR = R/NN-test_quantile. The fresh pair lives in a
# DIFFERENT directory (R/NN-test_quantile_d21_regen/samedraws), so it
# cannot be selected by token alone -- get_simul("NN", 21, quant="999")
# would always resolve to the shipped file no matter what. Following the
# pattern in 50_wp2a_nn_alpha_bigsets.R (source order: harness.R then
# wp0_mccd_methods.R, THEN install the shadow -- sourcing wp0_mccd_methods.R
# defines unmccd_method/sunmccd_method closures whose free reference to
# `get_simul` is resolved by the CURRENT .GlobalEnv binding at CALL time,
# so a shadow installed before sourcing would be clobbered by nothing [it
# isn't reassigned by sourcing], but installed after is the documented-safe
# order and is what 50/42/43/45 all use), this script installs
# get_simul_shadow() after both sources.
#
# The shadow is table-source-aware, not merely d-aware: a package-private
# flag TABLE_SOURCE_STATE$active in {"fresh","shipped",NA} is set by the
# driver immediately before each detector call and cleared immediately
# after (on.exit). When active and variant=="NN", the shadow builds the
# path directly from FRESH_NN_DIR or NN_QUANT_TABLE_DIR (whichever the
# active state names) and the caller's own `quant` token -- it does NOT
# fall through to nn_quant_for_d()/nn_quant_label_paper_*() bucketing,
# because both fresh and shipped requests here always pass an explicit
# alpha_token. When TABLE_SOURCE_STATE$active is NA (not set), the shadow
# defers unconditionally to the captured original get_simul(), so nothing
# outside this script's own cells is affected. Every cell records the
# absolute path, file size, and table vector length actually used (columns
# table_file, table_size_bytes, table_len) via a provenance tap (PROV),
# same discipline as 50_wp2a_nn_alpha_bigsets.R's PROV/get_simul_shadow.
#
# THREE GUARDS, run before any cell (see run_all_guards()):
#   Guard 1 (distinctness): fresh 99% vs fresh 999% must NOT be identical()
#     on $average or $median. Skipped gracefully (not aborted) if the fresh
#     files do not exist yet; aborts loudly if they exist AND are identical
#     (that would mean the new pair reproduces the original defect).
#   Guard 2 (coverage): each fresh table's $average/$median must have
#     length >= 3443 (waveform's n), or nnccd.radi's pmin(1:n, length(...))
#     clamp (harness.R:154-158) would silently reuse the table's last entry
#     for the remaining points with no warning. Reports actual lengths;
#     skipped gracefully if fresh files do not exist yet; aborts if a fresh
#     file exists and is short.
#   Guard 3 (shipped baseline): asserts (not merely assumes) that the
#     shipped d=21 99%/999% pair IS identical on both $average and $median,
#     documenting the comparison's baseline rather than taking it on faith.
#     Aborts if this ever turns out false (it always runs; nothing to wait
#     for).
#
# MECHANICS: command-line driven (methods, table_sources as args), per-cell
# CSV append via harness.R's append_result(), resume via has_result() keyed
# on (dataset, method, table_source, alpha_token), 3600s per-cell timeout
# via setTimeLimit() recorded as status="timeout", CELL_DONE/CELL_FAIL/
# ALL_CELLS_COMPLETE stdout lines -- all copied from
# 50_wp2a_nn_alpha_bigsets.R's run_cell()/run_guard() conventions.
#
# MODES:
#   --guard    Run only the three guards, print, exit. No cells, no methods
#              sourced beyond harness.R/wp0_mccd_methods.R. Cheap.
#   --dryrun   Run the three guards, then for each of the 6 cells resolve
#              and print exactly which table file would be used, its size,
#              and its vector length -- WITHOUT calling any detector.
#              Handles missing fresh files gracefully (prints "NOT YET
#              AVAILABLE", does not error).
#   --proof    Demonstrate the get_simul shadow's table-source routing on a
#              small case that does not depend on the (not-yet-existent)
#              fresh d=21 files: see run_proof() below for the exact route
#              taken (d=10, tokens "99"/"999", both genuinely distinct on
#              disk -- confirmed by file size/mtime, not assumed).
#   (default)  Real grid. Positional args:
#                Rscript 57_waveform_d21_rerun.R [methods] [table_sources]
#              methods:       comma list, default "UN-MCCD,SUN-MCCD"
#              table_sources: comma list, default "fresh,shipped"
#              (alpha_token is NOT a free argument here -- fresh always
#              runs both "99" and "999"; shipped always runs "999" only;
#              this fixes the grid at exactly the 6 cells specified.)

suppressMessages(library(here))
source(here::here("revision_experiments", "harness.R"))
source(here::here("revision_experiments", "wp0_mccd_methods.R"))

# ---------------------------------------------------------------------------
# get_simul shadow: table-source-aware, installed AFTER every source().
# ---------------------------------------------------------------------------
orig_get_simul <- get_simul   # capture BEFORE shadowing

FRESH_NN_DIR <- here::here("R/NN-test_quantile_d21_regen/samedraws")
FRESH_D      <- 21L

TABLE_SOURCE_STATE <- new.env(parent = emptyenv())
TABLE_SOURCE_STATE$active <- NA_character_   # "fresh" | "shipped" | NA

PROV <- new.env(parent = emptyenv())
PROV$last_file  <- NA_character_
PROV$last_size  <- NA_real_
PROV$last_quant <- NA_character_
PROV$last_avg_len <- NA_integer_
PROV$last_med_len <- NA_integer_

#' Load a raw NN .RData file directly (bypassing get_simul), for guards and
#' the --dryrun path resolver. Returns NULL (not an error) if missing.
load_nn_table_or_null <- function(dir, d, tok) {
  path <- file.path(dir, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))
  if (!file.exists(path)) return(NULL)
  e <- new.env()
  load(path, envir = e)
  if (!exists("simul", envir = e)) stop(sprintf("load_nn_table_or_null(): %s has no object named 'simul'", path))
  list(simul = get("simul", envir = e), file = path, size = file.size(path))
}

#' get_simul_shadow(): when TABLE_SOURCE_STATE$active is "fresh" or
#' "shipped" AND variant=="NN", resolve the path directly from FRESH_NN_DIR
#' or NN_QUANT_TABLE_DIR using the caller's explicit `quant` token (no
#' bucketing). Otherwise defer unconditionally to the captured original.
get_simul_shadow <- function(variant = c("RK", "NN"), d, quant = NULL) {
  active <- TABLE_SOURCE_STATE$active
  if (identical(variant, "NN") && !is.na(active)) {
    if (is.null(quant)) stop("get_simul_shadow(): TABLE_SOURCE_STATE active but quant is NULL -- every cell in this script must pass an explicit alpha_token")
    dir <- if (identical(active, "fresh")) FRESH_NN_DIR else NN_QUANT_TABLE_DIR
    fname <- sprintf("NN-test-simul_%dd_%s%%.RData", d, quant)
    path <- file.path(dir, fname)
    if (!file.exists(path)) {
      stop(sprintf("get_simul_shadow(): %s table missing for d=%d, quant=%s.\nExpected file: %s",
                    active, d, quant, path))
    }
    e <- new.env()
    load(path, envir = e)
    if (!exists("simul", envir = e)) stop(sprintf("get_simul_shadow(): file %s does not contain an object named 'simul'", path))
    res <- list(simul = get("simul", envir = e), quant = as.numeric(paste0("0.", quant)),
                quant_label = quant, file = path)
  } else {
    res <- orig_get_simul(variant, d, quant)
  }
  PROV$last_file    <- res$file
  PROV$last_size    <- suppressWarnings(file.size(res$file))
  PROV$last_quant   <- res$quant_label
  PROV$last_avg_len <- length(res$simul$average)
  PROV$last_med_len <- length(res$simul$median)
  res
}
assign("get_simul", get_simul_shadow, envir = globalenv())

# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------

#' Guard 3: shipped d=21 pair IS expected identical. Always runs (files
#' always exist). Aborts if the assumption is violated.
guard3_shipped_identical <- function(abort_on_fail = TRUE) {
  cat("---- Guard 3: shipped d=21 99%/999% pair should be identical (documented baseline) ----\n")
  a <- load_nn_table_or_null(NN_QUANT_TABLE_DIR, FRESH_D, "99")
  b <- load_nn_table_or_null(NN_QUANT_TABLE_DIR, FRESH_D, "999")
  if (is.null(a) || is.null(b)) {
    msg <- sprintf("Guard 3 FAILED: shipped d=%d table(s) missing (99: %s, 999: %s) -- cannot even establish the documented baseline.",
                    FRESH_D, !is.null(a), !is.null(b))
    cat("  ", msg, "\n")
    if (abort_on_fail) stop(msg)
    return(list(ok = FALSE, identical_avg = NA, identical_med = NA))
  }
  ia <- identical(a$simul$average, b$simul$average)
  im <- identical(a$simul$median,  b$simul$median)
  cat(sprintf("  shipped 99%%  file=%s size=%d bytes len(avg)=%d len(med)=%d\n",
              a$file, a$size, length(a$simul$average), length(a$simul$median)))
  cat(sprintf("  shipped 999%% file=%s size=%d bytes len(avg)=%d len(med)=%d\n",
              b$file, b$size, length(b$simul$average), length(b$simul$median)))
  cat(sprintf("  identical(average)=%s  identical(median)=%s\n", ia, im))
  if (!(ia && im)) {
    msg <- "Guard 3 FAILED: shipped d=21 99% vs 999% are NOT identical -- the documented baseline assumption (byte-identical shipped pair) no longer holds; the fresh-vs-shipped comparison this script exists to run would be comparing against a moving target. Investigate before proceeding."
    cat("  ***", msg, "\n")
    if (abort_on_fail) stop(msg)
    return(list(ok = FALSE, identical_avg = ia, identical_med = im))
  }
  cat("  Guard 3 PASSED: shipped d=21 pair confirmed byte-identical, as expected.\n\n")
  list(ok = TRUE, identical_avg = ia, identical_med = im)
}

#' Guard 1: fresh 99% vs fresh 999% must NOT be identical. Skips gracefully
#' (returns ok=NA) if either fresh file does not exist yet.
guard1_fresh_distinct <- function(abort_on_fail = TRUE) {
  cat("---- Guard 1: fresh d=21 99%/999% pair must be genuinely distinct ----\n")
  a <- load_nn_table_or_null(FRESH_NN_DIR, FRESH_D, "99")
  b <- load_nn_table_or_null(FRESH_NN_DIR, FRESH_D, "999")
  if (is.null(a) || is.null(b)) {
    cat(sprintf("  fresh 99%%: %s   fresh 999%%: %s\n",
                if (is.null(a)) "NOT YET AVAILABLE" else "present", if (is.null(b)) "NOT YET AVAILABLE" else "present"))
    cat("  Guard 1 SKIPPED (not a failure): fresh table generation has not landed yet.\n")
    cat(sprintf("  Expected location once done: %s\n\n", FRESH_NN_DIR))
    return(list(ok = NA, identical_avg = NA, identical_med = NA))
  }
  ia <- identical(a$simul$average, b$simul$average)
  im <- identical(a$simul$median,  b$simul$median)
  mx <- tryCatch(max(abs(a$simul$average - b$simul$average)), error = function(e) NA_real_)
  cat(sprintf("  fresh 99%%  file=%s size=%d bytes\n", a$file, a$size))
  cat(sprintf("  fresh 999%% file=%s size=%d bytes\n", b$file, b$size))
  cat(sprintf("  identical(average)=%s  identical(median)=%s  max|d avg|=%s\n", ia, im, format(mx)))
  if (ia && im) {
    msg <- "Guard 1 FAILED: fresh d=21 99% vs 999% are byte-identical -- the new pair reproduces the original defect this script exists to fix. Refusing to run the grid."
    cat("  ***", msg, "\n")
    if (abort_on_fail) stop(msg)
    return(list(ok = FALSE, identical_avg = ia, identical_med = im))
  }
  cat("  Guard 1 PASSED: fresh pair is genuinely distinct.\n\n")
  list(ok = TRUE, identical_avg = ia, identical_med = im)
}

#' Guard 2: each fresh table's average/median length must be >= 3443
#' (waveform's n). Skips gracefully if fresh files are absent.
guard2_fresh_length <- function(abort_on_fail = TRUE, min_len = 3443L) {
  cat(sprintf("---- Guard 2: fresh table length >= %d (waveform n) ----\n", min_len))
  toks <- c("99", "999")
  ok_all <- NA
  results <- list()
  any_present <- FALSE
  for (tok in toks) {
    t <- load_nn_table_or_null(FRESH_NN_DIR, FRESH_D, tok)
    if (is.null(t)) {
      cat(sprintf("  fresh %-4s%%: NOT YET AVAILABLE\n", tok))
      results[[tok]] <- list(avg_len = NA_integer_, med_len = NA_integer_, ok = NA)
      next
    }
    any_present <- TRUE
    al <- length(t$simul$average); ml <- length(t$simul$median)
    ok <- (al >= min_len) && (ml >= min_len)
    cat(sprintf("  fresh %-4s%%: len(average)=%d  len(median)=%d  %s\n",
                tok, al, ml, if (ok) "OK" else "*** TOO SHORT ***"))
    results[[tok]] <- list(avg_len = al, med_len = ml, ok = ok)
    if (isFALSE(ok_all) || is.na(ok_all)) ok_all <- ok else ok_all <- ok_all && ok
  }
  if (!any_present) {
    cat("  Guard 2 SKIPPED (not a failure): fresh tables not yet available.\n\n")
    return(list(ok = NA, results = results))
  }
  if (!isTRUE(ok_all)) {
    msg <- sprintf("Guard 2 FAILED: at least one fresh table has average/median length < %d -- nnccd.radi's pmin(1:n, length(...)) clamp (harness.R:154-158) would silently reuse the table's last entry for the rest of waveform's points. Refusing to run the grid.", min_len)
    cat("  ***", msg, "\n")
    if (abort_on_fail) stop(msg)
    return(list(ok = FALSE, results = results))
  }
  cat("  Guard 2 PASSED (for all fresh tables present).\n\n")
  list(ok = TRUE, results = results)
}

run_all_guards <- function(abort_on_fail = TRUE) {
  cat("==== 57_waveform_d21_rerun: guards ====\n\n")
  g3 <- guard3_shipped_identical(abort_on_fail)
  g1 <- guard1_fresh_distinct(abort_on_fail)
  g2 <- guard2_fresh_length(abort_on_fail)
  cat("==== guards complete ====\n\n")
  list(guard1 = g1, guard2 = g2, guard3 = g3)
}

# ---------------------------------------------------------------------------
# The 6-cell design table (dataset fixed to "waveform").
# ---------------------------------------------------------------------------
ALL_CELLS <- list(
  list(method = "UN-MCCD",  table_source = "fresh",   alpha_token = "99"),
  list(method = "UN-MCCD",  table_source = "fresh",   alpha_token = "999"),
  list(method = "SUN-MCCD", table_source = "fresh",   alpha_token = "99"),
  list(method = "SUN-MCCD", table_source = "fresh",   alpha_token = "999"),
  list(method = "UN-MCCD",  table_source = "shipped", alpha_token = "999"),
  list(method = "SUN-MCCD", table_source = "shipped", alpha_token = "999")
)

DATASET_LABEL <- "waveform"
CELL_TIMEOUT  <- 3600   # seconds, per task spec
MIN_CLS_SUN   <- 0.0625 # manuscript S_min
OUT_CSV       <- here::here("revision_experiments/results/tr1/wp2a_waveform_d21_rerun.csv")

collapse_idx <- function(v) if (length(v) == 0) "" else paste(sort(as.integer(v)), collapse = ";")
alpha_of_tok <- function(tok) 1 - as.numeric(paste0("0.", tok))

# ---------------------------------------------------------------------------
# resolve_cell_table(): resolve (without calling a detector) exactly which
# file a given cell would use, via the SAME shadow the real grid uses.
# Returns a list with ok=TRUE/FALSE(missing)/error, plus file/size/len.
# ---------------------------------------------------------------------------
resolve_cell_table <- function(table_source, alpha_token) {
  TABLE_SOURCE_STATE$active <- table_source
  on.exit(TABLE_SOURCE_STATE$active <- NA_character_, add = TRUE)
  out <- tryCatch({
    res <- get_simul("NN", FRESH_D, quant = alpha_token)
    list(ok = TRUE, file = res$file, size = suppressWarnings(file.size(res$file)),
         avg_len = length(res$simul$average), med_len = length(res$simul$median))
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  out
}

# ---------------------------------------------------------------------------
# run_cell(): the real per-cell driver, run_cell()'s conventions copied from
# 50_wp2a_nn_alpha_bigsets.R (setTimeLimit-based 3600s timeout, PROV tap,
# CELL_DONE/CELL_FAIL stdout lines, has_result()/append_result() resume).
# ---------------------------------------------------------------------------
run_cell <- function(method, table_source, alpha_token, X, Y, d, out_csv = OUT_CSV) {
  keys <- c(dataset = DATASET_LABEL, method = method, table_source = table_source, alpha_token = alpha_token)
  if (isTRUE(has_result(out_csv, keys))) {
    cat(sprintf("[skip, already recorded] %s x %s x table_source=%s x alpha_token=%s\n",
                DATASET_LABEL, method, table_source, alpha_token))
    return(invisible("skip"))
  }

  PROV$last_file <- NA_character_; PROV$last_size <- NA_real_; PROV$last_quant <- NA_character_
  PROV$last_avg_len <- NA_integer_; PROV$last_med_len <- NA_integer_

  TABLE_SOURCE_STATE$active <- table_source
  on.exit(TABLE_SOURCE_STATE$active <- NA_character_, add = TRUE)

  n <- nrow(X)
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
    res <- if (identical(method, "SUN-MCCD")) {
      METHOD_REGISTRY[[method]](X = X, d = d, Y = Y, quant = alpha_token, min.cls = MIN_CLS_SUN)
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

  flagged <- if (is.null(out$res)) integer(0) else which(out$res$score > REAL_DATA_THRESHOLDS[[method]])
  ncls <- if (is.null(out$res)) NA_integer_ else length(unique(out$res$cluster[!is.na(out$res$cluster)]))

  row <- list(
    dataset = DATASET_LABEL, method = method, table_source = table_source,
    alpha_token = alpha_token, alpha = alpha_of_tok(alpha_token),
    n = n, d = d,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA  = unname(out$m[["BA"]]),  F2  = unname(out$m[["F2"]]),
    n_flagged  = if (is.null(out$res)) NA_integer_ else length(flagged),
    n_clusters = ncls,
    flagged_idx = collapse_idx(flagged),
    table_file = PROV$last_file, table_size_bytes = PROV$last_size,
    table_len  = PROV$last_avg_len,
    elapsed_sec = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time())
  )
  append_result(out_csv, row)

  note_1line <- if (is.na(out$note)) "" else substr(gsub("[\r\n]+", " ", out$note), 1, 150)
  fields <- sprintf("dataset=%s method=%s table_source=%s alpha_token=%s status=%s sec=%.1f",
                     DATASET_LABEL, method, table_source, alpha_token, out$status, wall)
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
# --proof: demonstrate the get_simul shadow's routing without needing the
# not-yet-existent fresh d=21 files.
#
# ROUTE TAKEN: the fresh directory for d=21 does not exist yet, so per task
# instruction we fall back to "two different shipped tables at some other
# dimension where a genuine distinct pair exists". d=10 and d=20 both carry
# genuinely distinct 99%/999% pairs (confirmed by file size AND mtime,
# not assumed): 10d_99% is 70717 bytes / 2024-04-30, 10d_999% is 70940
# bytes / 2024-05-04; 20d_99% is 70151 bytes / 2024-05-09, 20d_999% is
# 70384 bytes / 2024-05-15. This script uses d=10.
#
# What is demonstrated:
#   (a) TABLE_SOURCE_STATE$active="shipped" with quant="99" returns exactly
#       the object in NN_QUANT_TABLE_DIR/NN-test-simul_10d_99%.RData.
#   (b) TABLE_SOURCE_STATE$active="fresh" with quant="999", but with
#       FRESH_NN_DIR temporarily REPOINTED at NN_QUANT_TABLE_DIR (since no
#       separate fresh directory exists for d=10 -- there is nothing else
#       it could physically point at), returns exactly the object in
#       NN_QUANT_TABLE_DIR/NN-test-simul_10d_999%.RData -- a DIFFERENT file
#       from (a), proving the shadow does not silently collapse "fresh" and
#       "shipped" requests onto the same table when their tokens differ.
#   (c) Directory-isolation check: FRESH_NN_DIR is pointed at a path that
#       does not exist, and a "fresh" request is shown to error out loudly
#       rather than silently falling back to the shipped directory -- this
#       is the actual failure mode the fresh-vs-shipped design must avoid,
#       and it is checked directly on the real FRESH_NN_DIR/get_simul_shadow
#       codepath (not a re-implementation).
run_proof <- function() {
  cat("==== 57_waveform_d21_rerun: --proof (get_simul shadow routing) ====\n")
  cat("Fresh d=21 directory not required to exist for this proof; using the\n")
  cat("task's documented fallback route: d=10, tokens 99%/999%, both present\n")
  cat("on disk and genuinely distinct (different size AND mtime).\n\n")

  d10_99  <- load_nn_table_or_null(NN_QUANT_TABLE_DIR, 10L, "99")
  d10_999 <- load_nn_table_or_null(NN_QUANT_TABLE_DIR, 10L, "999")
  if (is.null(d10_99) || is.null(d10_999)) stop("--proof: d=10 shipped tables unexpectedly missing")
  cat(sprintf("Reference: 10d_99%%  file=%s size=%d\n", d10_99$file, d10_99$size))
  cat(sprintf("Reference: 10d_999%% file=%s size=%d\n", d10_999$file, d10_999$size))
  cat(sprintf("Reference distinctness: identical(average)=%s identical(median)=%s\n\n",
              identical(d10_99$simul$average, d10_999$simul$average),
              identical(d10_99$simul$median,  d10_999$simul$median)))

  # (a) shipped route
  TABLE_SOURCE_STATE$active <- "shipped"
  res_shipped <- get_simul("NN", 10L, quant = "99")
  TABLE_SOURCE_STATE$active <- NA_character_
  ok_a <- identical(res_shipped$file, d10_99$file) && identical(res_shipped$simul$average, d10_99$simul$average)
  cat(sprintf("(a) shipped/quant=99  -> file=%s  matches reference 99%% object: %s\n", res_shipped$file, ok_a))

  # (b) fresh route, FRESH_NN_DIR temporarily repointed at the shipped dir
  #     (no separate physical fresh dir exists for d=10; see comment above)
  old_fresh_dir <- FRESH_NN_DIR
  FRESH_NN_DIR <<- NN_QUANT_TABLE_DIR
  TABLE_SOURCE_STATE$active <- "fresh"
  res_fresh <- get_simul("NN", 10L, quant = "999")
  TABLE_SOURCE_STATE$active <- NA_character_
  FRESH_NN_DIR <<- old_fresh_dir
  ok_b1 <- identical(res_fresh$file, d10_999$file) && identical(res_fresh$simul$average, d10_999$simul$average)
  ok_b2 <- !identical(res_fresh$simul$average, res_shipped$simul$average)
  cat(sprintf("(b) fresh/quant=999   -> file=%s  matches reference 999%% object: %s\n", res_fresh$file, ok_b1))
  cat(sprintf("    fresh result differs from (a)'s shipped/99%% result: %s (must be TRUE)\n\n", ok_b2))

  # (c) directory isolation: FRESH_NN_DIR points nowhere; "fresh" must error,
  #     not silently fall back to NN_QUANT_TABLE_DIR.
  FRESH_NN_DIR <<- here::here("revision_experiments/results/tr1/__does_not_exist__")
  TABLE_SOURCE_STATE$active <- "fresh"
  c_result <- tryCatch({
    get_simul("NN", 10L, quant = "99")
    "NO ERROR RAISED -- FAIL"
  }, error = function(e) conditionMessage(e))
  TABLE_SOURCE_STATE$active <- NA_character_
  FRESH_NN_DIR <<- old_fresh_dir
  ok_c <- grepl("missing for d=10", c_result)
  cat(sprintf("(c) fresh request against a nonexistent fresh dir raised an error: %s\n", ok_c))
  cat(sprintf("    error message: %s\n\n", c_result))

  all_ok <- ok_a && ok_b1 && ok_b2 && ok_c
  cat(sprintf("==== --proof verdict: %s ====\n", if (all_ok) "PASSED (shadow routes fresh vs shipped correctly and does not silently fall back)" else "*** FAILED -- see above ***"))
  if (!all_ok) stop("--proof: shadow routing proof failed")
  invisible(all_ok)
}

# ---------------------------------------------------------------------------
# --dryrun: guards + per-cell path resolution, no detector calls.
# ---------------------------------------------------------------------------
run_dryrun <- function() {
  run_all_guards(abort_on_fail = TRUE)
  cat("==== 57_waveform_d21_rerun: --dryrun (resolving all 6 cells, no detector calls) ====\n\n")
  for (cell in ALL_CELLS) {
    r <- resolve_cell_table(cell$table_source, cell$alpha_token)
    if (isTRUE(r$ok)) {
      cat(sprintf("dataset=%s method=%-9s table_source=%-8s alpha_token=%-4s -> file=%s size=%d bytes len(avg)=%d len(med)=%d\n",
                  DATASET_LABEL, cell$method, cell$table_source, cell$alpha_token, r$file, r$size, r$avg_len, r$med_len))
    } else {
      cat(sprintf("dataset=%s method=%-9s table_source=%-8s alpha_token=%-4s -> NOT YET AVAILABLE (%s)\n",
                  DATASET_LABEL, cell$method, cell$table_source, cell$alpha_token, r$error))
    }
  }
  cat("\n==== --dryrun complete: no detector was called, no CSV row was written ====\n")
}

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
MODE_GUARD  <- any(args == "--guard")
MODE_DRYRUN <- any(args == "--dryrun")
MODE_PROOF  <- any(args == "--proof")
args <- args[!args %in% c("--guard", "--dryrun", "--proof")]

if (MODE_GUARD) {
  run_all_guards(abort_on_fail = TRUE)
  quit(save = "no", status = 0)
}

if (MODE_PROOF) {
  run_proof()
  quit(save = "no", status = 0)
}

if (MODE_DRYRUN) {
  run_dryrun()
  quit(save = "no", status = 0)
}

# --- Real grid ---------------------------------------------------------
resolve_list <- function(s, default) if (is.null(s) || !nzchar(s)) default else strsplit(s, ",")[[1]]
REQ_METHODS  <- resolve_list(if (length(args) >= 1) args[1] else NULL, c("UN-MCCD", "SUN-MCCD"))
REQ_SOURCES  <- resolve_list(if (length(args) >= 2) args[2] else NULL, c("fresh", "shipped"))
stopifnot(all(REQ_METHODS %in% c("UN-MCCD", "SUN-MCCD")))
stopifnot(all(REQ_SOURCES %in% c("fresh", "shipped")))

run_all_guards(abort_on_fail = TRUE)

cells_to_run <- Filter(function(c) c$method %in% REQ_METHODS && c$table_source %in% REQ_SOURCES, ALL_CELLS)
cat(sprintf("57_waveform_d21_rerun: %d of %d cells selected (methods=%s, table_sources=%s)\n",
            length(cells_to_run), length(ALL_CELLS),
            paste(REQ_METHODS, collapse = ","), paste(REQ_SOURCES, collapse = ",")))

dat <- load_real_dataset(DATASET_LABEL)
stopifnot(dat$d == FRESH_D)
cat(sprintf("waveform loaded: n=%d d=%d n0(outliers)=%d\n\n", dat$n, dat$d, sum(dat$Y == 0)))

for (cell in cells_to_run) {
  run_cell(cell$method, cell$table_source, cell$alpha_token, dat$X, dat$Y, dat$d, out_csv = OUT_CSV)
}

cat("ALL_CELLS_COMPLETE\n")
