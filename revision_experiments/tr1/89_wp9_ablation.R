#!/usr/bin/env Rscript
# revision_experiments/89_wp9_ablation.R
#
# WP9 (R3.10a): SUN-MCCD component ablation driver. See
# revision_experiments/tr1/WP9_PROTOCOL.md for the full design and
# revision_experiments/tr1/wp9_sun_variants.R for the implementation
# (nnccd.radi.ablate, the environment-surgery substitution, the 5-variant
# grid, and the synthetic generators).
#
# Usage:
#   Rscript 89_wp9_ablation.R --smoke
#       1 setting (uniform_d3), all 5 variants, 1 rep, plus a bit-identity
#       assertion that the `stock` variant reproduces stock SUN-MCCD
#       (wp0_mccd_methods.R::sunmccd_method) exactly. Output:
#       results/tr1/wp9/smoke/wp9_ablation_smoke.csv
#
#   Rscript 89_wp9_ablation.R --run [reps] [budget_sec]
#       Full grid: 4 settings x 5 variants x `reps` (default 100) replicates.
#       Per-cell append (append_result/has_result from harness.R), so a
#       restart after a G: drive drop skips completed (setting, variant, rep)
#       cells. `budget_sec` (default 480) stops early so one invocation fits
#       the 600s tool timeout; re-invoke (same command) to continue. Output:
#       results/tr1/wp9/wp9_ablation.csv
#
#   Rscript 89_wp9_ablation.R --summarize
#       Builds results/tr1/wp9/wp9_ablation_summary.csv (mean/sd over reps
#       per setting x variant) from wp9_ablation.csv.
#
# NOTE: this driver was written and smoke-tested by the WP9 author; the FULL
# (--run) grid is NOT executed in this session -- only --smoke and parse()
# checks, per the task's "you do NOT run the real experiment" instruction.
# ---------------------------------------------------------------------------

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))
source(here::here("revision_experiments", "tr1", "wp0_mccd_methods.R"))
source(here::here("revision_experiments", "tr1", "wp9_sun_variants.R"))

RESDIR      <- here::here("revision_experiments/results/tr1/wp9")
SMOKE_DIR   <- file.path(RESDIR, "smoke")
OUT_CSV     <- file.path(RESDIR, "wp9_ablation.csv")
SMOKE_CSV   <- file.path(SMOKE_DIR, "wp9_ablation_smoke.csv")
SUMMARY_CSV <- file.path(RESDIR, "wp9_ablation_summary.csv")

MIN_CLS  <- 0.05     # S_min, CLAUDE.md's current constant (WP2c)
LOW_NUM  <- 3         # SUN-MCCD's own default (SUN-MCCD.R)
CELL_TIMEOUT <- 300   # per-cell hard cap (s)

# ---------------------------------------------------------------------------
# One (setting, variant, rep) cell
# ---------------------------------------------------------------------------
run_one_cell <- function(setting_row, variant_id, rep_id) {
  seed <- WP9_BASE_SEED + rep_id
  s <- setting_row
  dat <- if (s$generator == "uniform") gen_uniform(seed, s$n_nominal, s$d, s$contam)
         else                          gen_gaussian(seed, s$n_nominal, s$d, s$contam)
  X <- dat$X; n <- dat$n; n0 <- dat$n0
  Y <- c(rep(1, n - n0), rep(0, n0))   # outliers are the trailing n0 rows

  v <- WP9_VARIANTS[[variant_id]]
  t0 <- Sys.time()
  out <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = CELL_TIMEOUT, transient = TRUE)
    res <- sunmccd_ablate_method(X = X, d = s$d, Y = Y, method = v$method,
                                  min.cls = MIN_CLS, low.num = LOW_NUM,
                                  stat = v$stat, remove_centre = v$remove_centre)
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    stopifnot(length(res$score) == n, !anyNA(res$score))
    m <- evaluate(Y, res$score, REAL_DATA_THRESHOLDS[["SUN-MCCD"]])
    list(m = m, score = res$score, cluster = res$cluster, radii = res$radii,
         status = "ok", note = NA_character_)
  }, error = function(e) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    msg <- conditionMessage(e)
    list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
         score = NULL, cluster = NULL, radii = NULL,
         status = if (grepl("elapsed time limit|reached elapsed time limit", msg)) "timeout" else "error",
         note = msg)
  })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  n_clusters <- if (is.null(out$cluster)) NA_integer_ else length(unique(stats::na.omit(out$cluster)))
  radii_vec  <- if (is.null(out$radii)) NULL else out$radii[[1]]
  mean_radius <- if (is.null(radii_vec)) NA_real_ else mean(radii_vec)

  list(
    setting_id = s$setting_id, generator = s$generator, d = s$d, n = n, n0 = n0,
    rep = rep_id, seed = seed, variant_id = variant_id,
    stat = v$stat, remove_centre = v$remove_centre, method = v$method,
    TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
    BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
    n_flagged = if (is.null(out$score)) NA_integer_ else sum(out$score == 1),
    n_clusters = n_clusters, mean_radius = mean_radius,
    elapsed_sec = wall, status = out$status, note = out$note,
    timestamp = format(Sys.time()),
    score = out$score, radii_vec = radii_vec   # not written to CSV; used by --smoke's bit-identity check
  )
}

row_for_csv <- function(cell) {
  cell[!(names(cell) %in% c("score", "radii_vec"))]
}

# ---------------------------------------------------------------------------
# --smoke: 1 setting, all 5 variants, 1 rep, + bit-identity assertion
# ---------------------------------------------------------------------------
do_smoke <- function() {
  dir.create(SMOKE_DIR, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(SMOKE_CSV)) file.remove(SMOKE_CSV)  # smoke is always a clean rerun

  setting_row <- WP9_SETTINGS[WP9_SETTINGS$setting_id == "uniform_d3", ][1, ]
  cat(sprintf("wp9 smoke: setting=%s, variants=%s, rep=1\n",
              setting_row$setting_id, paste(names(WP9_VARIANTS), collapse = ", ")))

  cells <- list()
  for (vid in names(WP9_VARIANTS)) {
    cat(sprintf("  running variant %-20s ... ", vid)); flush.console()
    cell <- run_one_cell(setting_row, vid, rep_id = 1L)
    cells[[vid]] <- cell
    append_result(SMOKE_CSV, row_for_csv(cell))
    cat(sprintf("status=%s BA=%s F2=%s ncls=%s meanR=%s (%.2fs)\n",
                cell$status,
                ifelse(is.na(cell$BA), "NA", sprintf("%.3f", cell$BA)),
                ifelse(is.na(cell$F2), "NA", sprintf("%.3f", cell$F2)),
                cell$n_clusters,
                ifelse(is.na(cell$mean_radius), "NA", sprintf("%.3f", cell$mean_radius)),
                cell$elapsed_sec))
  }

  # ---- bit-identity assertion: stock variant vs stock SUN-MCCD ----
  cat("\nwp9 smoke: bit-identity check (stock ablate variant vs stock sunmccd_method) ...\n")
  seed <- WP9_BASE_SEED + 1L
  dat <- gen_uniform(seed, setting_row$n_nominal, setting_row$d, setting_row$contam)
  X <- dat$X; n <- dat$n; n0 <- dat$n0
  Y <- c(rep(1, n - n0), rep(0, n0))

  stock_ablate <- sunmccd_ablate_method(X = X, d = setting_row$d, Y = Y, method = "ascend",
                                         min.cls = MIN_CLS, low.num = LOW_NUM,
                                         stat = "both", remove_centre = TRUE)
  stock_plain  <- sunmccd_method(X = X, d = setting_row$d, Y = Y, method = "ascend",
                                  min.cls = MIN_CLS, low.num = LOW_NUM)

  score_ok <- identical(stock_ablate$score, stock_plain$score)
  radii_ok <- identical(stock_ablate$radii, stock_plain$radii)
  cluster_ok <- identical(stock_ablate$cluster, stock_plain$cluster)

  cat(sprintf("  score identical()   : %s\n", score_ok))
  cat(sprintf("  radii identical()   : %s\n", radii_ok))
  cat(sprintf("  cluster identical() : %s\n", cluster_ok))

  stopifnot(
    "BIT-IDENTITY FAILED: stock ablate variant's score differs from stock SUN-MCCD -- substitution mechanism is broken" = score_ok,
    "BIT-IDENTITY FAILED: stock ablate variant's radii differ from stock SUN-MCCD" = radii_ok,
    "BIT-IDENTITY FAILED: stock ablate variant's cluster assignment differs from stock SUN-MCCD" = cluster_ok
  )
  cat("wp9 smoke: BIT-IDENTITY OK -- stock variant reproduces stock SUN-MCCD exactly.\n")

  # cross-check: the "stock" row already written to SMOKE_CSV (rep=1, same
  # seed/setting) must ALSO match this independent stock_plain call.
  stopifnot(
    "wp9 smoke: the stock row's score (as run inside the variant loop) does not match the independent stock_plain re-check" =
      identical(cells[["stock"]]$score, stock_plain$score)
  )
  cat("wp9 smoke: cross-check OK -- variant-loop 'stock' row matches the independent re-check too.\n")

  invisible(cells)
}

# ---------------------------------------------------------------------------
# --run: full grid, checkpointed
# ---------------------------------------------------------------------------
do_run <- function(n_reps = 100L, budget = 480) {
  t_start <- Sys.time()
  for (i in seq_len(nrow(WP9_SETTINGS))) {
    s <- WP9_SETTINGS[i, ]
    for (vid in names(WP9_VARIANTS)) {
      for (rr in seq_len(n_reps)) {
        keys <- c(setting_id = s$setting_id, variant_id = vid, rep = as.character(rr))
        if (isTRUE(has_result(OUT_CSV, keys))) next
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        if (elapsed > budget) {
          cat(sprintf("[budget %gs reached after %.1fs -- rerun to continue] next up: %s x %s x rep %d\n",
                      budget, elapsed, s$setting_id, vid, rr))
          return(invisible("budget_stop"))
        }
        cell <- run_one_cell(s, vid, rr)
        append_result(OUT_CSV, row_for_csv(cell))
        cat(sprintf("  %-12s x %-20s x rep %3d  status=%-7s BA=%s F2=%s (%.2fs)\n",
                    s$setting_id, vid, rr, cell$status,
                    ifelse(is.na(cell$BA), "NA", sprintf("%.3f", cell$BA)),
                    ifelse(is.na(cell$F2), "NA", sprintf("%.3f", cell$F2)),
                    cell$elapsed_sec))
        flush.console()
      }
    }
  }
  cat("89_wp9_ablation: full grid complete.\n")
}

# ---------------------------------------------------------------------------
# --summarize
# ---------------------------------------------------------------------------
do_summarize <- function() {
  if (!file.exists(OUT_CSV)) { cat("no results yet\n"); return(invisible(NULL)) }
  df <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
  df <- df[df$status == "ok", ]
  rows <- list()
  for (sid in unique(df$setting_id)) for (vid in unique(df$variant_id[df$setting_id == sid])) {
    sub <- df[df$setting_id == sid & df$variant_id == vid, ]
    rows[[length(rows) + 1]] <- data.frame(
      setting_id = sid, variant_id = vid,
      stat = sub$stat[1], remove_centre = sub$remove_centre[1], method = sub$method[1],
      n_reps = nrow(sub),
      TPR_mean = mean(sub$TPR), TPR_sd = sd(sub$TPR),
      TNR_mean = mean(sub$TNR), TNR_sd = sd(sub$TNR),
      BA_mean  = mean(sub$BA),  BA_sd  = sd(sub$BA),
      F2_mean  = mean(sub$F2),  F2_sd  = sd(sub$F2),
      n_clusters_mean = mean(sub$n_clusters),
      mean_radius_mean = mean(sub$mean_radius),
      stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, rows)
  write.csv(S, SUMMARY_CSV, row.names = FALSE)
  cat(sprintf("wrote %s (%d rows)\n", SUMMARY_CSV, nrow(S)))
  invisible(S)
}

# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if ("--smoke" %in% args) {
  do_smoke()
} else if ("--run" %in% args) {
  rest <- args[args != "--run"]
  n_reps <- if (length(rest) >= 1 && nzchar(rest[1])) as.integer(rest[1]) else 100L
  budget <- if (length(rest) >= 2 && nzchar(rest[2])) as.numeric(rest[2]) else 480
  do_run(n_reps, budget)
} else if ("--summarize" %in% args) {
  do_summarize()
} else {
  cat("Usage: Rscript 89_wp9_ablation.R --smoke | --run [reps] [budget_sec] | --summarize\n")
}
