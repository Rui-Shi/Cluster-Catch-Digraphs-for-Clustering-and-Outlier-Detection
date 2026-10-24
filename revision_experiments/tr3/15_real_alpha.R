#!/usr/bin/env Rscript
# revision_experiments/tr3/15_real_alpha.R
#
# Step 2 of ALPHA_PROTOCOL.md: every (config, level, data set) cell of the
# significance-level experiment for UN-CCDs on the nine real data sets.
#
# USAGE (from the clone root, in PowerShell, CCD_QUANTILE_GENERATE=FALSE)
#   Rscript revision_experiments/tr3/15_real_alpha.R
#
# Configs   A  table d = d of the loaded data; shipped table if the file exists,
#              else the generated one (14_alpha_tables.R); all 9 data sets
#           P  the script's d (user_knowledge 6, wholesale 8); same source rule
#           G  generated table at every level; all 9 data sets
#           S  as A, with ascend for seeds and breast_cancer
# Levels    20, 15, 10, 5, 1 (percent) = tables 80, 85, 90, 95, 99
# Direction the script's, from settings_realdata.csv (S: ascend)
# Metrics   as 11_real.R: ARI = adjustedRandIndex(label, Y), Sil =
#           mean(silhouette(label, dist(X))[, 3]) (NA for one cluster),
#           k_hat = length(unique(label))
# Tables    generated tables are loaded and their `ccd_table` and `generated`
#           attributes are removed in memory, so nnccd.radi() uses the table as
#           it is and never grows it. Shipped tables go through tr3_load_table().
# Workers   one PSOCK cluster of 20, parLapplyLB, D31 tasks first. Each task is
#           wrapped in tryCatch; its row (with the error text, if any) goes to
#           <RES>/real_alpha_parts/<config>_<label>_<dataset>.csv, and a task
#           whose part file exists is skipped. At the end the parts are merged
#           into <RES>/real_alpha.csv.
# A missing or unusable table, or a data set of the wrong size, is not a cell
# result: the script stops before any task runs.
#
# Test switches (the full run uses none of them; defaults in brackets)
#   TR3_ALPHA_RES       base folder for outputs [revision_experiments/results/tr3]
#   TR3_ALPHA_NITER     the niter the generated tables must record [10000]
#   TR3_ALPHA_DATASETS  comma-separated subset of data sets [all nine]
#   TR3_ALPHA_LEVELS    comma-separated table labels, e.g. 90,95 [80,85,90,95,99]
#   TR3_ALPHA_CONFIGS   comma-separated subset of A,P,G,S [A,P,G,S]
#   TR3_ALPHA_WORKERS   number of workers [20]

source(here::here("revision_experiments/tr3/tr3_lib.R"))

RES      <- Sys.getenv("TR3_ALPHA_RES", TR3_RES)
TABDIR   <- file.path(RES, "alpha_tables")
PARTDIR  <- file.path(RES, "real_alpha_parts")
OUT_CSV  <- file.path(RES, "real_alpha.csv")
LOG      <- file.path(RES, "_logs", "alpha_runs.log")
NITER    <- as.integer(Sys.getenv("TR3_ALPHA_NITER", "10000"))
WORKERS  <- as.integer(Sys.getenv("TR3_ALPHA_WORKERS", "20"))
ROOT     <- normalizePath(here::here(), winslash = "/")

dir.create(PARTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG, append = TRUE)
}

# ---- the protocol's design ---------------------------------------------------
LABELS <- c("80", "85", "90", "95", "99")
ALPHA  <- c("80" = "20", "85" = "15", "90" = "10", "95" = "5", "99" = "1")
# loaded size, true k and table d of the script (protocol table), in protocol order
DS <- data.frame(
  dataset  = c("iris", "seeds", "user_knowledge", "wholesale", "breast_cancer",
               "aggregation", "asymmetric", "R15", "D31"),
  n        = c(147, 210, 258, 440, 449, 788, 1000, 600, 3100),
  d        = c(4, 7, 5, 6, 9, 2, 2, 2, 2),
  script_d = c(4, 7, 6, 8, 9, 2, 2, 2, 2),
  stringsAsFactors = FALSE)
EXTENT <- c("2" = 1000, "4" = 147, "5" = 258, "6" = 440, "7" = 210, "8" = 440, "9" = 449)

sel <- function(env, all) {
  v <- Sys.getenv(env, "")
  if (nzchar(v)) strsplit(v, ",")[[1]] else all
}
DATASETS <- sel("TR3_ALPHA_DATASETS", DS$dataset)
LEVELSEL <- sel("TR3_ALPHA_LEVELS", LABELS)
CONFIGS  <- sel("TR3_ALPHA_CONFIGS", c("A", "P", "G", "S"))
stopifnot(all(DATASETS %in% DS$dataset), all(LEVELSEL %in% LABELS), all(CONFIGS %in% c("A", "P", "G", "S")))

ship_name <- function(d, label) paste0("NN-test-simul_", d, "d_", label, "%.RData")
gen_name  <- function(d, label) paste0("NN-test-simul_", d, "d_", label, "%_gen_n", EXTENT[[as.character(d)]], ".RData")
ship_path <- function(d, label) here::here("R/NN-test_quantile", ship_name(d, label))

# ---- the task list -----------------------------------------------------------
settings <- lapply(setNames(DS$dataset, DS$dataset), function(ds) tr3_real_setting("UN-CCDs", ds))
tasks <- list()
for (cfg in CONFIGS) {
  for (ds in DATASETS) {
    if (cfg == "P" && !(ds %in% c("user_knowledge", "wholesale"))) next
    if (cfg == "S" && !(ds %in% c("seeds", "breast_cancer"))) next
    i <- match(ds, DS$dataset)
    d_tab <- if (cfg == "P") DS$script_d[i] else DS$d[i]
    direction <- if (cfg == "S") "ascend" else settings[[ds]]$direction
    for (lab in LEVELSEL) {
      shipped <- cfg != "G" && file.exists(ship_path(d_tab, lab))
      tasks[[length(tasks) + 1]] <- list(
        config = cfg, level = ALPHA[[lab]], label = lab, dataset = ds,
        d_table = d_tab, direction = direction,
        source = if (shipped) "shipped" else "generated",
        file = if (shipped) ship_name(d_tab, lab) else gen_name(d_tab, lab),
        path = if (shipped) ship_path(d_tab, lab) else file.path(TABDIR, gen_name(d_tab, lab)),
        n_expected = DS$n[i], d_expected = DS$d[i],
        part = file.path(PARTDIR, sprintf("%s_%s_%s.csv", cfg, lab, ds)))
    }
  }
}
logmsg("15_real_alpha.R start: %d tasks (configs %s; data sets %s; levels %s), %d workers, RES=%s",
       length(tasks), paste(CONFIGS, collapse = ","), paste(DATASETS, collapse = ","),
       paste(LEVELSEL, collapse = ","), WORKERS, RES)

# ---- pre-flight: infrastructure problems stop the run, they are not cell results
real_master <- tr3_load_real(DATASETS)
for (ds in DATASETS) {
  i <- match(ds, DS$dataset)
  D <- real_master[[ds]]
  if (nrow(D) != DS$n[i] || ncol(D) - 1 != DS$d[i]) {
    stop(sprintf("pre-flight: %s loads as %d x %d, protocol says %d x %d", ds, nrow(D), ncol(D) - 1,
                 DS$n[i], DS$d[i]))
  }
}
rm(real_master)
checked <- character(0)
for (tk in tasks) {
  if (tk$path %in% checked) next
  checked <- c(checked, tk$path)
  if (!file.exists(tk$path)) stop("pre-flight: missing table ", tk$path)
  if (tk$source == "generated") {
    e <- new.env(); load(tk$path, envir = e)
    g <- attr(e$simul, "generated")
    ext <- EXTENT[[as.character(tk$d_table)]]
    ok <- !is.null(g) && identical(as.integer(g$niter), NITER) && identical(as.integer(g$cores), 22L) &&
      identical(as.integer(g$seed), as.integer(20260927 + tk$d_table)) &&
      identical(as.integer(g$d), as.integer(tk$d_table)) && identical(as.integer(g$extent), as.integer(ext)) &&
      isTRUE(all.equal(g$level, as.numeric(tk$label) / 100)) &&
      length(e$simul$average) == ext && length(e$simul$median) == ext &&
      !anyNA(e$simul$average) && !anyNA(e$simul$median)
    if (!ok) stop("pre-flight: generated table does not match the protocol (niter, cores, seed, d, level, extent): ", tk$path)
  }
}
logmsg("pre-flight passed: %d table files, data set sizes as in the protocol", length(checked))

# ---- one task, on a worker ---------------------------------------------------
run_task <- function(tk) {
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch({
    if (tk$source == "shipped") {
      simul <- tr3_load_table(tk$file)
    } else {
      e <- new.env(); load(tk$path, envir = e)
      simul <- e$simul
      attr(simul, "ccd_table") <- NULL
      attr(simul, "generated") <- NULL
    }
    D <- real[[tk$dataset]]
    n_col <- ncol(D)
    X <- D[, -n_col]
    Y <- D[, n_col]
    res <- ccd_cluster("UN-CCDs", X, list(simul = simul, direction = tk$direction))
    label <- res$label
    n_clusters <- length(unique(label))
    ari <- adjustedRandIndex(label, Y)
    sil <- if (n_clusters > 1) mean(silhouette(label, dist(X))[, 3]) else NA_real_
    list(n = nrow(X), d = ncol(X), ari = ari, sil = sil, k_hat = n_clusters, error = "")
  }, error = function(e) {
    list(n = NA, d = NA, ari = NA_real_, sil = NA_real_, k_hat = NA,
         error = gsub("[\r\n]+", " ", conditionMessage(e)))
  })
  secs <- proc.time()[["elapsed"]] - t0
  row <- data.frame(config = tk$config, level = tk$level, dataset = tk$dataset,
                    table = tk$file, source = tk$source,
                    n = if (is.na(out$n)) "NA" else as.character(out$n),
                    d = if (is.na(out$d)) "NA" else as.character(out$d),
                    ari = fmt17(out$ari), sil = fmt17(out$sil),
                    k_hat = if (is.na(out$k_hat)) "NA" else as.character(out$k_hat),
                    seconds = fmt17(secs), error = out$error, stringsAsFactors = FALSE)
  tmp <- paste0(tk$part, ".tmp")
  utils::write.csv(row, tmp, row.names = FALSE)
  if (!file.rename(tmp, tk$part)) stop("could not rename ", tmp)
  cat(sprintf("[%s] %s %s%% %s: %.1f s%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tk$config,
              tk$level, tk$dataset, secs, if (nzchar(out$error)) paste0(" ERROR: ", out$error) else ""),
      file = tk$logfile, append = TRUE)
  invisible(NULL)
}

todo <- Filter(function(tk) !file.exists(tk$part), tasks)
logmsg("%d tasks to run, %d part files already exist", length(todo), length(tasks) - length(todo))
if (length(todo)) {
  ord <- order(vapply(todo, function(tk) tk$n_expected, numeric(1)), decreasing = TRUE)   # D31 first
  todo <- lapply(todo[ord], function(tk) { tk$logfile <- LOG; tk })
  cl <- parallel::makeCluster(min(WORKERS, length(todo)))
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterCall(cl, function(root, ds) {
    setwd(root)
    source(file.path(root, "revision_experiments/tr3/tr3_lib.R"))
    tr3_source_method("UN-CCDs", "realdata")
    assign("real", tr3_load_real(ds), envir = globalenv())
    invisible(NULL)
  }, ROOT, DATASETS)
  logmsg("workers ready")
  parallel::clusterExport(cl, "run_task")
  t0 <- Sys.time()
  invisible(parallel::parLapplyLB(cl, todo, run_task, chunk.size = 1))
  logmsg("tasks done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  parallel::stopCluster(cl)
}

# ---- merge -------------------------------------------------------------------
parts <- list.files(PARTDIR, pattern = "\\.csv$", full.names = TRUE)
rows <- lapply(parts, function(p) utils::read.csv(p, stringsAsFactors = FALSE, colClasses = "character"))
merged <- do.call(rbind, rows)
merged <- merged[order(match(merged$config, c("A", "P", "G", "S")), -as.numeric(merged$level),
                 match(merged$dataset, DS$dataset)), ]
merged$error[is.na(merged$error)] <- ""
tmp <- paste0(OUT_CSV, ".tmp")
utils::write.csv(merged, tmp, row.names = FALSE)
file.rename(tmp, OUT_CSV)
nerr <- sum(nzchar(merged$error))
logmsg("merged %d part files into %s; %d cells with an error", length(parts), OUT_CSV, nerr)
if (nerr) for (i in which(nzchar(merged$error))) logmsg("  failed: %s %s%% %s: %s", merged$config[i], merged$level[i], merged$dataset[i], merged$error[i])
logmsg("15_real_alpha.R done")
