#!/usr/bin/env Rscript
# 99_table_seed_stability.R -- how much do the NND-based detectors (UN-MCCD,
# SUN-MCCD) change when the Monte Carlo quantile table is rebuilt under an
# independent seed?  Protocol: TABLE_SEED_STABILITY_PROTOCOL.md (same folder).
#
#   Rscript revision_experiments/tr1/99_table_seed_stability.R --anchor
#   Rscript ... --generate <d> <seeds,csv> [--dry-run] [--cores 22] [--niter 10000]
#   Rscript ... --run <datasets,csv> [--methods UN-MCCD,SUN-MCCD]
#   Rscript ... --summary
#
# Run from the repository root (here::here() finds it from the working directory).
#
# Modes
#   --anchor    run the production-table cells (3 data sets x 2 methods) and compare
#               TPR/TNR/BA/F2 with the published values to three decimals; exits
#               with status 2 on any mismatch.
#   --generate  build NND quantile tables (99% and 99.9% from the same 10,000 draws)
#               under R/NN-test_quantile_seedstab/<d>d/seed<k>/ ; chunked, resumable.
#   --run       run every available table for the named data sets.
#   --summary   comparisons, per-(data set, method) summary and the md5 manifest.
#
# Nothing here writes to R/NN-test_quantile/ or R/NN-test_quantile_regen999/.

suppressMessages(library(here))

args <- commandArgs(trailingOnly = TRUE)
has  <- function(flag) any(args == flag)
opt  <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}

MODE <- if (has("--anchor")) "anchor" else if (has("--generate")) "generate" else
        if (has("--run")) "run" else if (has("--summary")) "summary" else
        stop("give one of --anchor, --generate, --run, --summary")

REPO     <- here::here()
OUT_DIR  <- here::here("revision_experiments/results/tr1/table_seed")
CELLS    <- file.path(OUT_DIR, "table_seed_cells.csv")
COMPAR   <- file.path(OUT_DIR, "table_seed_comparisons.csv")
SUMMARY  <- file.path(OUT_DIR, "table_seed_summary.csv")
MANIFEST <- file.path(OUT_DIR, "tables.csv")

PROD_DIR  <- here::here("R/NN-test_quantile")
REGEN_DIR <- here::here("R/NN-test_quantile_regen999")
SEED_ROOT <- here::here("R/NN-test_quantile_seedstab")

S_MIN        <- 0.05          # label-free constant used for every real data set
CELL_TIMEOUT <- 7200          # seconds; generous, a timed-out cell would be silently lost
MIN_CLS_METHODS <- "SUN-MCCD"
METHODS_ALL  <- c("UN-MCCD", "SUN-MCCD")
DATASETS_ALL <- c("hepatitis", "lymphography", "vowels")
N_OF <- c("12" = 1452L, "18" = 148L, "19" = 74L)   # n of the data set that fixes the table extent

# Published values (supplement real-data tables), TPR TNR BA F2.
PUBLISHED <- rbind(
  data.frame(dataset = "hepatitis",    method = "UN-MCCD",  TPR = 0.714, TNR = 0.657, BA = 0.686, F2 = 0.446),
  data.frame(dataset = "hepatitis",    method = "SUN-MCCD", TPR = 0.857, TNR = 0.552, BA = 0.705, F2 = 0.469),
  data.frame(dataset = "lymphography", method = "UN-MCCD",  TPR = 1.000, TNR = 0.859, BA = 0.930, F2 = 0.600),
  data.frame(dataset = "lymphography", method = "SUN-MCCD", TPR = 1.000, TNR = 0.880, BA = 0.940, F2 = 0.638),
  data.frame(dataset = "vowels",       method = "UN-MCCD",  TPR = 1.000, TNR = 0.541, BA = 0.771, F2 = 0.263),
  data.frame(dataset = "vowels",       method = "SUN-MCCD", TPR = 0.978, TNR = 0.338, BA = 0.658, F2 = 0.194))

# ---------------------------------------------------------------------------
# Generic helpers (no harness needed)
# ---------------------------------------------------------------------------

rel_path <- function(p) {
  root <- paste0(normalizePath(REPO, winslash = "/"), "/")
  q <- normalizePath(p, winslash = "/", mustWork = FALSE)
  if (startsWith(q, root)) substring(q, nchar(root) + 1L) else q
}
fmt_num <- function(v) if (length(v) == 0) "" else paste(sprintf("%.15g", v), collapse = ";")
fmt_int <- function(v) if (length(v) == 0) "" else paste(as.integer(v), collapse = ";")
parse_num <- function(s) if (is.na(s) || !nzchar(s)) numeric(0) else suppressWarnings(as.numeric(strsplit(s, ";")[[1]]))
parse_int <- function(s) if (is.na(s) || !nzchar(s)) integer(0) else suppressWarnings(as.integer(strsplit(s, ";")[[1]]))

adj_rand <- function(a, b) {
  keep <- !is.na(a) & !is.na(b); a <- a[keep]; b <- b[keep]
  n <- length(a); if (n < 2) return(NA_real_)
  if (length(unique(a)) == 1 && length(unique(b)) == 1) return(1)
  tab <- table(a, b); cc <- function(x) sum(x * (x - 1) / 2)
  sij <- cc(as.vector(tab)); sa <- cc(rowSums(tab)); sb <- cc(colSums(tab))
  tot <- n * (n - 1) / 2; expected <- sa * sb / tot; maxi <- (sa + sb) / 2
  if (isTRUE(all.equal(maxi, expected))) return(1)
  (sij - expected) / (maxi - expected)
}
jaccard <- function(a, b) {
  u <- length(union(a, b)); if (u == 0) return(1)
  length(intersect(a, b)) / u
}

table_file <- function(dir, d, label) file.path(dir, sprintf("NN-test-simul_%dd_%s%%.RData", d, label))
md5_of <- function(p) unname(tools::md5sum(p))

# complete seed-stability table sets (COMPLETE marker written after both tables)
seed_dirs <- function(d) {
  base <- file.path(SEED_ROOT, sprintf("%dd", d))
  if (!dir.exists(base)) return(character(0))
  ds <- list.dirs(base, recursive = FALSE, full.names = TRUE)
  ds <- ds[grepl("seed[0-9]+$", ds)]
  ds <- ds[file.exists(file.path(ds, "COMPLETE"))]
  ds[order(as.integer(sub("^.*seed", "", ds)))]
}

# ---------------------------------------------------------------------------
# --generate
# ---------------------------------------------------------------------------

chunk_seed <- function(k, d, i) as.integer(20261010 + 100000 * k + 1000 * d + i)

#' Seeds already used elsewhere in the code base (for the collision check).
seeds_in_use <- function() {
  s71 <- as.vector(outer(c(12, 16, 18, 19, 21), 1:60, function(d, i) 20260811 + 1000 * d + i))
  s71b <- 20260811 + 1000 * c(12, 16, 18, 19, 21)                  # the "seed" stored in 71's table meta
  s64 <- as.vector(outer(1:100, 1:600, function(r, d) 20260810 + 1000 * r + d))
  sqt <- 20260927 + 1:600                                           # R/ccds/quantile_table.R default
  s01 <- 20260809 + 1:600                                           # 01i_nn_multiquant_table.R default
  list("71 chunks" = s71, "71 table meta" = s71b, "64" = s64, "quantile_table.R" = sqt, "01i default" = s01)
}

# est. core-seconds per iteration at 10,000 draws on 22 cores (user-supplied
# timings for these tables: 12 s at d=19, 45 s at d=18, 7.7 h at d=12)
EST_CORE_SEC <- c("19" = 12 * 22 / 10000, "18" = 45 * 22 / 10000, "12" = 7.7 * 3600 * 22 / 10000)
TARGET_CHUNK_SEC <- 1200

chunk_plan <- function(d, cores, niter) {
  est_wall <- unname(EST_CORE_SEC[as.character(d)]) * niter / cores
  nchunk <- max(1L, min(40L, as.integer(round(est_wall / TARGET_CHUNK_SEC))))
  perchk <- as.integer(ceiling(niter / nchunk))
  nchunk <- as.integer(ceiling(niter / perchk))
  list(est_wall = est_wall, nchunk = nchunk, perchk = perchk,
       iters = c(rep(perchk, nchunk - 1L), niter - perchk * (nchunk - 1L)))
}

read_plan <- function(f) {
  l <- readLines(f, warn = FALSE); kv <- strsplit(l, "=", fixed = TRUE)
  setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
}

do_generate <- function() {
  gi <- which(args == "--generate")
  d  <- as.integer(args[gi + 1])
  ks <- as.integer(strsplit(args[gi + 2], ",")[[1]])
  stopifnot(length(d) == 1, !is.na(d), !anyNA(ks), as.character(d) %in% names(N_OF))
  dry   <- has("--dry-run")
  cores <- as.integer(opt("--cores", 22L))
  niter <- as.integer(opt("--niter", 10000L))
  n     <- N_OF[[as.character(d)]]

  # ---- seed collision check ------------------------------------------------
  plan0 <- chunk_plan(d, cores, niter)
  new_seeds <- unlist(lapply(ks, function(k) chunk_seed(k, d, seq_len(plan0$nchunk))))
  stopifnot(!anyDuplicated(new_seeds))
  cat("seed collision check (new seeds vs seeds in use elsewhere):\n")
  for (nm in names(seeds_in_use())) {
    ov <- intersect(new_seeds, seeds_in_use()[[nm]])
    cat(sprintf("  %-18s %7d seeds in use, overlap = %d\n", nm, length(seeds_in_use()[[nm]]), length(ov)))
    if (length(ov) > 0) stop("seed collision with ", nm)
  }
  cat(sprintf("  new seeds: %d, range %d..%d, all distinct\n", length(new_seeds), min(new_seeds), max(new_seeds)))

  for (k in ks) {
    outdir <- file.path(SEED_ROOT, sprintf("%dd", d), sprintf("seed%d", k))
    chkdir <- file.path(outdir, "chunks")
    seeds_k <- chunk_seed(k, d, seq_len(plan0$nchunk))
    cat(sprintf("\n=== seed index %d, d=%d, n=%d ===\n  outdir   : %s\n  chunks   : %d x %d iterations (last %d), est. wall %.2f h at %d cores\n  chunk seeds: %s\n",
                k, d, n, rel_path(outdir), plan0$nchunk, plan0$perchk, tail(plan0$iters, 1),
                plan0$est_wall / 3600, cores, paste(seeds_k, collapse = ",")))
    cat(sprintf("  tables   : %s ; %s\n", rel_path(table_file(outdir, d, "99")), rel_path(table_file(outdir, d, "999"))))
    if (file.exists(file.path(outdir, "COMPLETE"))) { cat("  already COMPLETE, skipped\n"); next }
    if (dry) { cat("  [dry run: nothing created, nothing computed]\n"); next }

    suppressMessages(source(here::here("revision_experiments/tr1/01i_nn_multiquant_table.R")))
    dir.create(chkdir, recursive = TRUE, showWarnings = FALSE)
    planf <- file.path(outdir, "PLAN.txt")
    plan_lines <- c(sprintf("n=%d", n), sprintf("d=%d", d), sprintf("niter=%d", niter), sprintf("cores=%d", cores),
                    sprintf("seed_index=%d", k), sprintf("nchunk=%d", plan0$nchunk), sprintf("perchk=%d", plan0$perchk),
                    sprintf("engine=fast_stream"), sprintf("chunk_seeds=%s", paste(seeds_k, collapse = ",")))
    if (file.exists(planf)) {
      old <- read_plan(planf)
      for (key in c("n", "d", "niter", "cores", "nchunk", "perchk", "chunk_seeds")) {
        new <- sub(paste0("^", key, "="), "", plan_lines[startsWith(plan_lines, paste0(key, "="))])
        if (!identical(unname(old[key]), new)) stop(sprintf("PLAN.txt mismatch on %s (on disk %s, now %s): refusing to resume with a different plan", key, old[key], new))
      }
    } else writeLines(plan_lines, planf)

    log_f <- file.path(outdir, "PROGRESS.tsv")
    say <- function(fmt, ...) {
      line <- sprintf("%s\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), sprintf(fmt, ...))
      cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = log_f, append = TRUE); flush.console()
    }
    chunk_file <- function(i) file.path(chkdir, sprintf("draws_chunk%03d.RData", i))
    chunk_ok <- function(i) {
      f <- chunk_file(i); if (!file.exists(f)) return(FALSE)
      e <- new.env(); ok <- tryCatch({ load(f, envir = e)
        is.matrix(e$NN.dist.ave.mat) && nrow(e$NN.dist.ave.mat) == plan0$iters[i] && ncol(e$NN.dist.ave.mat) == n &&
        identical(dim(e$NN.dist.ave.mat), dim(e$NN.dist.med.mat)) && identical(as.integer(e$nn_meta$seed), seeds_k[i]) },
        error = function(e) FALSE)
      isTRUE(ok)
    }
    say("RUNSTART\td=%d\tk=%d\tcores=%d\tniter=%d", d, k, cores, niter)
    for (i in seq_len(plan0$nchunk)) {
      if (chunk_ok(i)) { say("CHUNKSKIP\tchunk=%d/%d", i, plan0$nchunk); next }
      if (file.exists(chunk_file(i))) say("CHUNKBAD\tchunk=%d/%d unreadable or wrong shape, recomputing", i, plan0$nchunk)
      t0 <- Sys.time()
      dr <- nn_mc_draws(n = n, d = d, niter = plan0$iters[i], cores = cores, seed = seeds_k[i], engine = "fast_stream")
      NN.dist.ave.mat <- dr$ave; NN.dist.med.mat <- dr$med; nn_meta <- dr$meta
      tmp <- paste0(chunk_file(i), ".tmp")
      save(NN.dist.ave.mat, NN.dist.med.mat, nn_meta, file = tmp)
      if (!file.rename(tmp, chunk_file(i))) stop("could not move chunk into place: ", chunk_file(i))
      rm(dr, NN.dist.ave.mat, NN.dist.med.mat); gc(verbose = FALSE)
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      say("CHUNK\tchunk=%d/%d\titer=%d\tsecs=%.1f\teta_h=%.2f", i, plan0$nchunk, plan0$iters[i], el, el * (plan0$nchunk - i) / 3600)
    }
    say("POOL\tloading %d chunks", plan0$nchunk)
    ave <- matrix(0, nrow = niter, ncol = n); med <- matrix(0, nrow = niter, ncol = n); row <- 0L
    for (i in seq_len(plan0$nchunk)) {
      e <- new.env(); load(chunk_file(i), envir = e)
      idx <- (row + 1L):(row + nrow(e$NN.dist.ave.mat))
      ave[idx, ] <- e$NN.dist.ave.mat; med[idx, ] <- e$NN.dist.med.mat; row <- row + nrow(e$NN.dist.ave.mat)
      rm(e); gc(verbose = FALSE)
    }
    stopifnot(row == niter)
    pooled <- list(ave = ave, med = med,
                   meta = list(n = n, d = d, niter = niter, cores = cores, seed = seeds_k[1],
                               seed_index = k, chunk_seeds = seeds_k, engine = "fast_stream",
                               chunks = plan0$nchunk, perchk = plan0$perchk,
                               generated = format(Sys.time()), R_version = R.version.string,
                               provenance = "99_table_seed_stability.R"))
    w <- write_nn_tables(pooled, c(0.99, 0.999), outdir, overwrite = TRUE)   # outdir has no COMPLETE: safe to redo
    m99 <- md5_of(table_file(outdir, d, "99")); m999 <- md5_of(table_file(outdir, d, "999"))
    writeLines(c(sprintf("finished=%s", format(Sys.time())), sprintf("md5_99=%s", m99), sprintf("md5_999=%s", m999)),
               file.path(outdir, "COMPLETE"))
    z <- w$zero_frac
    say("DONE\td=%d\tk=%d\tmd5_99=%s\tmd5_999=%s\tzero_frac_999=%.4f", d, k, m99, m999,
        z$frac_zero_average_excl_pos1[z$token == "999"])
    rm(ave, med, pooled, w); gc(verbose = FALSE)
  }
  cat("\n99_table_seed_stability: generate finished\n")
}

# ---------------------------------------------------------------------------
# --anchor / --run
# ---------------------------------------------------------------------------

setup_detectors <- function() {
  # generation OFF, before any source()
  options(ccd.quantile.generate = FALSE)
  Sys.setenv(CCD_QUANTILE_GENERATE = "FALSE")
  suppressMessages({
    source(here::here("revision_experiments", "shared", "harness.R"))
    source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))
  })
  stopifnot("generation must be switched off" = !.ccd_qt_generate_allowed())

  TAP <<- new.env(parent = emptyenv())
  TAP$radii <- NULL; TAP$radii_n <- 0L
  TAP$intended <- NA_character_; TAP$files <- character(0); TAP$sums <- NULL

  # taps installed only now, after every source() has run
  local({
    orig_nnrad <- get("nnccd.radi", envir = globalenv())
    assign("nnccd.radi", function(...) {
      v <- orig_nnrad(...)
      TAP$radii_n <- TAP$radii_n + 1L
      if (is.null(TAP$radii)) TAP$radii <- v$R
      v
    }, envir = globalenv())

    orig_get_simul <- get("get_simul", envir = globalenv())
    viol <- function(msg) stop(structure(class = c("tableTapViolation", "error", "condition"),
                                         list(message = msg, call = NULL)))
    assign("get_simul", function(...) {
      res <- orig_get_simul(...)
      f <- normalizePath(res$file, winslash = "/", mustWork = TRUE)
      TAP$files <- c(TAP$files, f)
      if (length(TAP$files) != 1L) viol(sprintf("get_simul called %d times in one cell", length(TAP$files)))
      if (!identical(f, TAP$intended)) viol(sprintf("table mismatch: resolved %s, intended %s", f, TAP$intended))
      if (!is.null(attr(res$simul, "ccd_table")) || !is.null(attr(res$simul, "generated")))
        viol("loaded table carries a ccd_table/generated attribute (generated on the spot?)")
      TAP$sums <- c(avg = sum(res$simul$average), med = sum(res$simul$median),
                    len = min(length(res$simul$average), length(res$simul$median)))
      res
    }, envir = globalenv())
  })
}

tap_reset <- function(intended) {
  TAP$radii <- NULL; TAP$radii_n <- 0L
  TAP$intended <- normalizePath(intended, winslash = "/", mustWork = TRUE)
  TAP$files <- character(0); TAP$sums <- NULL
}

run_cell <- function(meth, X, Y, d, intended_file, intended_dir) {
  assign("NN_QUANT_TABLE_DIR", intended_dir, envir = globalenv())   # get_simul() reads this global
  tap_reset(intended_file)
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
    res <- if (meth %in% MIN_CLS_METHODS) {
      METHOD_REGISTRY[[meth]](X = X, d = d, Y = Y, min.cls = S_MIN)
    } else METHOD_REGISTRY[[meth]](X = X, d = d, Y = Y)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    list(res = res, status = "ok", note = NA_character_)
  }, tableTapViolation = function(e) {          # never turned into an error row
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    stop(e)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    list(res = NULL, status = if (grepl("elapsed time limit", msg)) "timeout" else "error", note = msg)
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(out$res)) return(list(status = out$status, note = out$note, elapsed = wall, score = NULL,
                                    metrics = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
                                    radii = numeric(0), cluster = integer(0), unassigned = NA_integer_))
  if (length(TAP$files) != 1L) stop("table tap did not fire exactly once")
  m <- evaluate(Y, out$res$score, REAL_DATA_THRESHOLDS[[meth]])
  list(status = "ok", note = NA_character_, elapsed = wall, score = out$res$score, metrics = m,
       radii = if (is.null(TAP$radii)) numeric(0) else TAP$radii,
       cluster = out$res$cluster, unassigned = out$res$unassigned_rows)
}

tables_for <- function(d, meth) {
  lab <- if (meth == "UN-MCCD") nn_quant_label_paper_UN(d) else nn_quant_label_paper_SUN(d)
  rows <- list(data.frame(table_id = "prod", role = "production", seed_index = NA_integer_, dir = PROD_DIR, label = lab,
                          stringsAsFactors = FALSE))
  if (meth == "UN-MCCD" && lab == "99")
    rows[[length(rows) + 1]] <- data.frame(table_id = "regen", role = "regen999_fresh_1pct", seed_index = NA_integer_,
                                           dir = file.path(REGEN_DIR, sprintf("%dd", d)), label = lab, stringsAsFactors = FALSE)
  for (sd in seed_dirs(d))
    rows[[length(rows) + 1]] <- data.frame(table_id = basename(sd), role = "new_seed", seed_index = as.integer(sub("^seed", "", basename(sd))),
                                           dir = sd, label = lab, stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

done_keys <- function() {
  if (!file.exists(CELLS)) return(character(0))
  df <- tryCatch(read.csv(CELLS, stringsAsFactors = FALSE, colClasses = "character"), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(character(0))
  df <- df[df$status == "ok", , drop = FALSE]
  paste(df$dataset, df$method, df$table_id, sep = "|")
}

run_grid <- function(datasets, methods, only_prod = FALSE) {
  setup_detectors()
  for (ds in datasets) {
    dat <- load_real_dataset(ds)
    n_exp <- N_OF[[as.character(dat$d)]]
    stopifnot("data set n does not match the n the tables are built for" = dat$n == n_exp)
    for (meth in methods) {
      tabs <- tables_for(dat$d, meth)
      if (only_prod) tabs <- tabs[tabs$table_id == "prod", , drop = FALSE]
      prod_md5 <- md5_of(table_file(PROD_DIR, dat$d, tabs$label[1]))
      for (j in seq_len(nrow(tabs))) {
        tb <- tabs[j, ]
        key <- paste(ds, meth, tb$table_id, sep = "|")
        f <- table_file(tb$dir, dat$d, tb$label)
        if (key %in% done_keys()) { cat(sprintf("[skip done] %s\n", key)); next }
        if (!file.exists(f)) { cat(sprintf("[skip missing table] %s : %s\n", key, rel_path(f))); next }
        md5 <- md5_of(f)
        if (tb$table_id != "prod" && identical(md5, prod_md5))
          stop(sprintf("%s: alternative table is byte-identical to the production table (%s); refusing to record it", key, md5))
        r <- run_cell(meth, dat$X, dat$Y, dat$d, f, tb$dir)
        fl <- if (is.null(r$score)) integer(0) else which(r$score > 0.5)
        mt <- r$metrics
        row <- list(
          dataset = ds, method = meth, table_id = tb$table_id, table_role = tb$role,
          table_seed_index = tb$seed_index, table_label = tb$label,
          table_path = rel_path(f), table_md5 = md5,
          tap_path = if (length(TAP$files)) rel_path(TAP$files[1]) else NA_character_,
          tap_calls = length(TAP$files),
          table_len = if (is.null(TAP$sums)) NA_real_ else unname(TAP$sums["len"]),
          table_sum_average = if (is.null(TAP$sums)) NA_real_ else unname(TAP$sums["avg"]),
          table_sum_median  = if (is.null(TAP$sums)) NA_real_ else unname(TAP$sums["med"]),
          n = dat$n, d = dat$d, n_outliers = sum(dat$Y == 0),
          TPR = unname(mt[["TPR"]]), TNR = unname(mt[["TNR"]]), BA = unname(mt[["BA"]]), F2 = unname(mt[["F2"]]),
          n_flagged = if (is.null(r$score)) NA_integer_ else length(fl),
          n_clusters = if (length(r$cluster)) length(unique(r$cluster[!is.na(r$cluster)])) else NA_integer_,
          unassigned_rows = r$unassigned, radii_len = length(r$radii), radii_calls = TAP$radii_n,
          flagged_idx = fmt_int(fl), cluster_vec = fmt_int(r$cluster), radii_vec = fmt_num(r$radii),
          s_min = if (meth %in% MIN_CLS_METHODS) S_MIN else NA_real_,
          elapsed_sec = r$elapsed, status = r$status, note = r$note, timestamp = format(Sys.time()))
        append_result(CELLS, row)
        cat(sprintf("  %-12s %-9s %-6s TPR=%.3f TNR=%.3f BA=%.3f F2=%.3f nflag=%s ncl=%s %s %.1fs md5=%s\n",
                    ds, meth, tb$table_id, mt[["TPR"]], mt[["TNR"]], mt[["BA"]], mt[["F2"]],
                    if (is.null(r$score)) "NA" else length(fl), row$n_clusters, r$status, r$elapsed, substr(md5, 1, 8)))
        flush.console()
      }
    }
  }
}

do_anchor <- function() {
  cat("--anchor: production-table cells vs published values (3 decimals)\n")
  run_grid(DATASETS_ALL, METHODS_ALL, only_prod = TRUE)
  df <- read.csv(CELLS, stringsAsFactors = FALSE)
  df <- df[df$table_id == "prod" & df$status == "ok", ]
  allok <- TRUE
  for (i in seq_len(nrow(PUBLISHED))) {
    p <- PUBLISHED[i, ]; r <- df[df$dataset == p$dataset & df$method == p$method, ]
    if (nrow(r) != 1) { cat(sprintf("  MISSING %s %s\n", p$dataset, p$method)); allok <- FALSE; next }
    got <- round(c(r$TPR, r$TNR, r$BA, r$F2), 3); pub <- c(p$TPR, p$TNR, p$BA, p$F2)
    ok <- all(abs(got - pub) < 1e-9); allok <- allok && ok
    cat(sprintf("  %-12s %-9s got %s | published %s | %s\n", p$dataset, p$method,
                paste(sprintf("%.3f", got), collapse = " "), paste(sprintf("%.3f", pub), collapse = " "), if (ok) "MATCH" else "MISMATCH"))
  }
  if (!allok) { cat("ANCHOR FAILED\n"); quit(status = 2) }
  cat("ANCHOR OK: all six production cells match\n")
}

do_run <- function() {
  ds <- strsplit(args[which(args == "--run") + 1], ",")[[1]]
  me <- strsplit(opt("--methods", paste(METHODS_ALL, collapse = ",")), ",")[[1]]
  stopifnot(all(ds %in% DATASETS_ALL), all(me %in% METHODS_ALL))
  run_grid(ds, me)
  cat("99_table_seed_stability: run finished\n")
}

# ---------------------------------------------------------------------------
# --summary
# ---------------------------------------------------------------------------

chunk_meta <- function(chdir) {
  fs <- sort(list.files(chdir, pattern = "^draws_chunk[0-9]+\\.RData$", full.names = TRUE))
  if (!length(fs)) return(NULL)
  ms <- lapply(fs, function(f) { e <- new.env(); load(f, envir = e); e$nn_meta })
  list(seeds = vapply(ms, function(m) as.integer(m$seed), 1L), cores = unique(vapply(ms, function(m) as.integer(m$cores), 1L)),
       engine = unique(vapply(ms, function(m) as.character(m$engine), "")), niter = sum(vapply(ms, function(m) as.integer(m$niter), 1L)),
       nchunk = length(fs))
}

build_manifest <- function() {
  cells <- if (file.exists(CELLS)) read.csv(CELLS, stringsAsFactors = FALSE) else NULL
  rows <- list()
  add <- function(role, table_id, d, label, f, chdir = NULL, note = "") {
    if (!file.exists(f)) return()
    e <- new.env(); load(f, envir = e); s <- e$simul
    cm <- if (!is.null(chdir) && dir.exists(chdir)) chunk_meta(chdir) else NULL
    md5 <- md5_of(f)
    rows[[length(rows) + 1]] <<- data.frame(
      role = role, table_id = table_id, d = d, label = label, path = rel_path(f), md5 = md5,
      n_entries = length(s$average), niter = if (!is.null(cm)) cm$niter else NA, cores = if (!is.null(cm)) paste(cm$cores, collapse = ";") else NA,
      engine = if (!is.null(cm)) paste(cm$engine, collapse = ";") else NA, n_chunks = if (!is.null(cm)) cm$nchunk else NA,
      chunk_seeds = if (!is.null(cm)) paste(cm$seeds, collapse = ";") else "legacy/unknown",
      used_in_cells = if (is.null(cells)) 0L else sum(cells$table_md5 == md5 & cells$status == "ok"),
      note = note, stringsAsFactors = FALSE)
  }
  for (d in c(12, 18, 19)) {
    for (lab in c("99", "999")) {
      fp <- table_file(PROD_DIR, d, lab); fr <- table_file(file.path(REGEN_DIR, sprintf("%dd", d)), d, lab)
      same <- file.exists(fp) && file.exists(fr) && identical(md5_of(fp), md5_of(fr))
      add("production", "prod", d, lab, fp, if (same) file.path(REGEN_DIR, sprintf("%dd", d), "chunks") else NULL,
          if (same) "byte-identical to the regen999 table (installed copy of 71_regen_nnd_alpha001.R output)" else
                    "legacy table (generator R/NN-test_quantile/*.R has no set.seed); seed/niter/engine not recoverable from the file")
      if (!same) add("regen999_fresh", "regen", d, lab, fr, file.path(REGEN_DIR, sprintf("%dd", d), "chunks"), "71_regen_nnd_alpha001.R")
    }
    for (sd in seed_dirs(d)) for (lab in c("99", "999"))
      add("new_seed", basename(sd), d, lab, table_file(sd, d, lab), file.path(sd, "chunks"), "99_table_seed_stability.R")
  }
  m <- do.call(rbind, rows)
  write.csv(m, MANIFEST, row.names = FALSE)
  cat(sprintf("wrote %s (%d tables)\n", rel_path(MANIFEST), nrow(m)))
  m
}

do_summary <- function() {
  stopifnot(file.exists(CELLS))
  chr <- c(flagged_idx = "character", cluster_vec = "character", radii_vec = "character")
  df <- read.csv(CELLS, stringsAsFactors = FALSE, colClasses = chr)
  for (k in names(chr)) df[[k]][is.na(df[[k]])] <- ""
  df <- df[df$status == "ok", ]
  comp <- list(); summ <- list()
  for (ds in unique(df$dataset)) for (meth in unique(df$method[df$dataset == ds])) {
    sub <- df[df$dataset == ds & df$method == meth, ]
    pr <- sub[sub$table_id == "prod", ]
    if (nrow(pr) != 1) { cat(sprintf("no unique production row for %s %s, skipped\n", ds, meth)); next }
    n <- pr$n; fp <- parse_int(pr$flagged_idx); cp <- parse_int(pr$cluster_vec); rp <- parse_num(pr$radii_vec)
    alt <- sub[sub$table_id != "prod", ]
    cl <- list()
    for (j in seq_len(nrow(alt))) {
      a <- alt[j, ]; fa <- parse_int(a$flagged_idx); ca <- parse_int(a$cluster_vec); ra <- parse_num(a$radii_vec)
      pos <- length(ra) == length(rp) && length(rp) > 0
      rel <- if (pos) { k <- rp > 0; abs(ra[k] - rp[k]) / rp[k] } else numeric(0)
      cl[[j]] <- data.frame(
        dataset = ds, method = meth, table_id = a$table_id, table_role = a$table_role, n = n, d = a$d,
        md5_differs_from_prod = !identical(a$table_md5, pr$table_md5),
        share_labels_changed = (length(setdiff(fa, fp)) + length(setdiff(fp, fa))) / n,
        n_labels_changed = length(setdiff(fa, fp)) + length(setdiff(fp, fa)),
        flagged_prod = length(fp), flagged_alt = length(fa), jaccard = jaccard(fa, fp),
        ari = adj_rand(ca, cp), clusters_prod = pr$n_clusters, clusters_alt = a$n_clusters,
        TPR_prod = pr$TPR, TNR_prod = pr$TNR, BA_prod = pr$BA, F2_prod = pr$F2,
        TPR_alt = a$TPR, TNR_alt = a$TNR, BA_alt = a$BA, F2_alt = a$F2,
        dBA = a$BA - pr$BA, dF2 = a$F2 - pr$F2,
        radius_median_rel_change = if (length(rel)) median(rel) else NA_real_,
        radius_max_rel_change = if (length(rel)) max(rel) else NA_real_,
        radius_n_points_compared = length(rel),
        radius_zero_flip = if (pos) sum((ra == 0) != (rp == 0)) else NA_integer_,
        stringsAsFactors = FALSE)
    }
    cdf <- if (length(cl)) do.call(rbind, cl) else NULL
    if (!is.null(cdf)) comp[[length(comp) + 1]] <- cdf
    allr <- sub
    summ[[length(summ) + 1]] <- data.frame(
      dataset = ds, method = meth, n = n, d = pr$d,
      n_tables_total = nrow(allr), n_tables_compared = nrow(alt),
      tables = paste(allr$table_id, collapse = ";"),
      all_md5_differ_from_prod = if (is.null(cdf)) NA else all(cdf$md5_differs_from_prod),
      F2_min = min(allr$F2), F2_max = max(allr$F2), BA_min = min(allr$BA), BA_max = max(allr$BA),
      F2_prod = pr$F2, BA_prod = pr$BA,
      jaccard_mean = if (is.null(cdf)) NA_real_ else mean(cdf$jaccard), jaccard_min = if (is.null(cdf)) NA_real_ else min(cdf$jaccard),
      share_labels_changed_max = if (is.null(cdf)) NA_real_ else max(cdf$share_labels_changed),
      n_labels_changed_max = if (is.null(cdf)) NA_integer_ else max(cdf$n_labels_changed),
      n_tables_labels_identical = if (is.null(cdf)) NA_integer_ else sum(cdf$n_labels_changed == 0),
      ari_min = if (is.null(cdf)) NA_real_ else min(cdf$ari, na.rm = TRUE),
      clusters_prod = pr$n_clusters, clusters_min = min(allr$n_clusters), clusters_max = max(allr$n_clusters),
      radius_median_rel_change_max_over_tables = if (is.null(cdf)) NA_real_ else max(cdf$radius_median_rel_change, na.rm = TRUE),
      radius_max_rel_change_max_over_tables = if (is.null(cdf)) NA_real_ else max(cdf$radius_max_rel_change, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  comp <- if (length(comp)) do.call(rbind, comp) else data.frame()
  summ <- do.call(rbind, summ)
  write.csv(comp, COMPAR, row.names = FALSE); write.csv(summ, SUMMARY, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows), %s (%d rows)\n", rel_path(COMPAR), nrow(comp), rel_path(SUMMARY), nrow(summ)))
  build_manifest()
  options(width = 220)
  print(summ[, c("dataset", "method", "n_tables_total", "all_md5_differ_from_prod", "F2_min", "F2_max", "BA_min", "BA_max",
                 "jaccard_mean", "jaccard_min", "share_labels_changed_max", "ari_min", "clusters_prod", "clusters_min", "clusters_max",
                 "radius_median_rel_change_max_over_tables", "radius_max_rel_change_max_over_tables")], row.names = FALSE, digits = 4)
}

switch(MODE, anchor = do_anchor(), generate = do_generate(), run = do_run(), summary = do_summary())
