#!/usr/bin/env Rscript
# revision_experiments/tr3/14_alpha_tables.R
#
# Step 1 of ALPHA_PROTOCOL.md: the seven NN-test quantile tables that the
# significance-level experiment needs, one Monte Carlo run per (d, extent),
# reduced at levels 0.80, 0.85, 0.90, 0.95 and 0.99.
#
# USAGE (from the clone root, in PowerShell, CCD_QUANTILE_GENERATE=FALSE)
#   Rscript revision_experiments/tr3/14_alpha_tables.R
#
# Draws:     .ccd_qt_draws(.ccd_qt_nn_once, list(n = extent, d = d),
#            niter = 10000, cores = 22, seed = 20260927 + d)   (22 workers fixed)
# Reduction: as the NN branch of ccd_generate_quantile_table(): the lower
#            quantile at 1 - level of the mean and of the median NND, per number
#            of points.
# Output:    <RES>/alpha_tables/NN-test-simul_<d>d_<label>%_gen_n<extent>.RData
#            holding `simul` (average, median) with a `generated` attribute.
#            <RES>/alpha_table_check.csv   generated against shipped tables
#            <RES>/_logs/alpha_tables.log
# A table file that already exists is skipped.
#
# Test switches (the full run uses none of them; defaults in brackets)
#   TR3_ALPHA_RES      base folder for outputs [revision_experiments/results/tr3]
#   TR3_ALPHA_NITER    Monte Carlo iterations [10000]; must be >= 22
#   TR3_ALPHA_ONLY_D   comma-separated d to generate [all seven]

source(here::here("R/ccds/quantile_table.R"))   # .ccd_qt_draws, .ccd_qt_nn_once only
suppressPackageStartupMessages(library(MASS))

RES    <- Sys.getenv("TR3_ALPHA_RES", here::here("revision_experiments/results/tr3"))
TABDIR <- file.path(RES, "alpha_tables")
LOG    <- file.path(RES, "_logs", "alpha_tables.log")
CHECK  <- file.path(RES, "alpha_table_check.csv")
NITER  <- as.integer(Sys.getenv("TR3_ALPHA_NITER", "10000"))
CORES  <- 22L
LEVELS <- c("80" = 0.80, "85" = 0.85, "90" = 0.90, "95" = 0.95, "99" = 0.99)
# (d, extent), small tables first and d = 2 last
TABLES <- data.frame(d = c(4, 7, 5, 6, 8, 9, 2), extent = c(147, 210, 258, 440, 440, 449, 1000))
only_d <- Sys.getenv("TR3_ALPHA_ONLY_D", "")
if (nzchar(only_d)) TABLES <- TABLES[TABLES$d %in% as.integer(strsplit(only_d, ",")[[1]]), ]
stopifnot(NITER >= CORES)

dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(LOG), recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = LOG, append = TRUE)
}

gen_file <- function(d, label, extent) {
  file.path(TABDIR, paste0("NN-test-simul_", d, "d_", label, "%_gen_n", extent, ".RData"))
}
ship_file <- function(d, label) {
  here::here(paste0("R/NN-test_quantile/NN-test-simul_", d, "d_", label, "%.RData"))
}

# generated against shipped: largest relative difference over entries 2..min(extent, shipped length)
write_check <- function() {
  rows <- list()
  for (i in seq_len(nrow(TABLES))) {
    d <- TABLES$d[i]; ext <- TABLES$extent[i]
    for (lab in names(LEVELS)) {
      gf <- gen_file(d, lab, ext); sf <- ship_file(d, lab)
      if (!file.exists(gf) || !file.exists(sf)) next
      e1 <- new.env(); load(gf, envir = e1)
      e2 <- new.env(); load(sf, envir = e2)
      m <- min(ext, length(e2$simul$average), length(e2$simul$median))
      idx <- 2:m
      rel <- function(g, s) max(abs(g[idx] - s[idx]) / abs(s[idx]))
      rows[[length(rows) + 1]] <- data.frame(
        d = d, level = paste0(lab, "%"), extent = ext, shipped_length = length(e2$simul$average),
        entries_compared = length(idx),
        max_rel_diff_average = sprintf("%.17g", rel(e1$simul$average, e2$simul$average)),
        max_rel_diff_median  = sprintf("%.17g", rel(e1$simul$median, e2$simul$median)),
        stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) {
    out <- do.call(rbind, rows)
    tmp <- paste0(CHECK, ".tmp")
    utils::write.csv(out, tmp, row.names = FALSE)
    file.rename(tmp, CHECK)
    logmsg("table check written: %d rows", nrow(out))
  }
}

logmsg("14_alpha_tables.R start: niter=%d cores=%d RES=%s", NITER, CORES, RES)
t_all <- Sys.time()
for (i in seq_len(nrow(TABLES))) {
  d <- TABLES$d[i]; ext <- TABLES$extent[i]
  files <- vapply(names(LEVELS), function(lab) gen_file(d, lab, ext), character(1))
  if (all(file.exists(files))) {
    logmsg("d=%d extent=%d: all five files exist, skipped", d, ext)
    next
  }
  seed <- 20260927 + d
  logmsg("d=%d extent=%d: draws start (niter=%d, cores=%d, seed=%d)", d, ext, NITER, CORES, seed)
  t0 <- Sys.time()
  out <- .ccd_qt_draws(.ccd_qt_nn_once, list(n = as.integer(ext), d = as.integer(d)),
                       NITER, CORES, seed)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ave <- do.call(rbind, lapply(out, `[[`, "ave"))
  med <- do.call(rbind, lapply(out, `[[`, "med"))
  rm(out)
  logmsg("d=%d extent=%d: draws done in %.0f s (matrix %d x %d)", d, ext, secs, nrow(ave), ncol(ave))
  for (lab in names(LEVELS)) {
    f <- gen_file(d, lab, ext)
    if (file.exists(f)) { logmsg("  %s%%: exists, skipped", lab); next }
    quant <- LEVELS[[lab]]
    # reduction as in ccd_generate_quantile_table() (lower tail, 1 - quant)
    q.ave <- sapply(1:ext, function(x) quantile(ave[, x], 1 - quant))
    names(q.ave) <- NULL
    q.med <- sapply(1:ext, function(x) quantile(med[, x], 1 - quant))
    names(q.med) <- NULL
    simul <- list(average = q.ave, median = q.med)
    attr(simul, "generated") <- list(d = as.integer(d), level = quant, extent = as.integer(ext),
                                     niter = NITER, cores = CORES, seed = as.integer(seed),
                                     seconds = secs, date = format(Sys.time()),
                                     R = R.version.string)
    tmp <- paste0(f, ".tmp")
    save(simul, file = tmp)
    if (!file.rename(tmp, f)) stop("could not rename ", tmp)
    logmsg("  %s%%: saved %s", lab, basename(f))
  }
  rm(ave, med); gc(verbose = FALSE)
  write_check()
}
write_check()
logmsg("14_alpha_tables.R done in %.1f min", as.numeric(difftime(Sys.time(), t_all, units = "mins")))
