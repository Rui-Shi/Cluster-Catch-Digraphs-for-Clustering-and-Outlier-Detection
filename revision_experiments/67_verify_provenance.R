#!/usr/bin/env Rscript
# 67_verify_provenance.R -- independent re-derivation of the alpha-provenance
# verdicts produced by 64/65.
#
# WHAT 65 CONCLUDED, AND WHY IT NEEDS AN INDEPENDENT CHECK
# --------------------------------------------------------
# 65 reports that the duplicated NN quantile tables hold the LESS EXTREME of
# their two candidate levels: 5% at d = 6,7,8,9 and 1% at d = 18,19. Those
# verdicts decide whether Table tab:alpha_real in the manuscript is true, so
# they are re-derived here from the raw Monte-Carlo draws with an estimator
# written independently of 65's code. Nothing is taken from 65 except the
# draws themselves.
#
# THE ESTIMATOR. A shipped table entry S_k is claimed to be the alpha-level
# lower-tail quantile of the NND statistic at subsample size k. Pool the fresh
# draws at that k across all replicate runs and read off the empirical CDF
# position of S_k:
#
#     level_k = mean(pooled_draws_k <= S_k)
#
# If the claim is right this concentrates on alpha. Taking the median over the
# several hundred available k averages out both the shipped table's own Monte
# Carlo error and the pooled draws'. This is a DIRECT estimate of the level --
# it does not ask "which of two candidates is nearer", so unlike a two-sample
# comparison it can also return a value matching NEITHER candidate, which is
# the outcome that would signal a deeper problem.
#
# WHY THE CONTROLS CARRY THE ARGUMENT. The estimator is only as good as the
# assumption that the fresh generator samples the same distribution as
# whatever produced the shipped tables. That is not assumed here, it is
# tested: d = 5, 10 and 20 have genuinely distinct pairs whose levels ARE
# known from their filenames, and the estimator must return those levels. Six
# such files are available spanning alpha = 10%, 5%, 1%, 0.1%. If the controls
# do not come back on their nominal values the unknown verdicts are void.
#
# COLUMN ALIGNMENT IS ALSO TESTED, NOT ASSUMED. The draw matrix has one column
# per subsample size and the shipped vector has 5000 entries; an off-by-one
# between them would silently bias every estimate. Both candidate alignments
# are evaluated on the controls and the one that reproduces the nominal levels
# is reported, together with what the other one would have given.
#
# Read-only apart from its own CSV.

suppressMessages(library(here))
PDIR  <- here::here("R/NN-test_quantile_provenance")
SDIR  <- here::here("R/NN-test_quantile")
OUT   <- here::here("revision_experiments/results/tr1/wp2a_provenance_verify.csv")

reps <- sort(list.files(PDIR, pattern = "^rep[0-9]+$"))
cat(sprintf("replicate directories found: %d (%s)\n", length(reps),
            paste(range(reps), collapse = " .. ")))

load1 <- function(p, obj) { e <- new.env(); load(p, envir = e); get(obj, envir = e) }

#' Pool the per-iteration draw matrices for one (d, n) across every replicate.
pool_draws <- function(d) {
  fs <- unlist(lapply(reps, function(r)
    list.files(file.path(PDIR, r), pattern = sprintf("^NN-draws_%dd_n[0-9]+\\.RData$", d),
               full.names = TRUE)))
  if (!length(fs)) return(NULL)
  ave <- do.call(rbind, lapply(fs, load1, obj = "NN.dist.ave.mat"))
  med <- do.call(rbind, lapply(fs, load1, obj = "NN.dist.med.mat"))
  list(ave = ave, med = med, n_files = length(fs))
}

#' Empirical CDF position of each shipped entry within the pooled draws.
#' shift = 0 pairs shipped[k] with column k; shift = 1 pairs it with column k+1.
level_curve <- function(shipped, draws, shift) {
  K <- min(length(shipped), ncol(draws) - shift)
  vapply(seq_len(K), function(k)
    mean(draws[, k + shift] <= shipped[k]), numeric(1))
}

# ---- targets ---------------------------------------------------------------
# truth = NA for the duplicated dimensions, the known level for the controls.
targets <- rbind(
  data.frame(d =  5, tok = "90",  truth = 0.10),
  data.frame(d =  5, tok = "95",  truth = 0.05),
  data.frame(d = 10, tok = "99",  truth = 0.01),
  data.frame(d = 10, tok = "999", truth = 0.001),
  data.frame(d = 20, tok = "99",  truth = 0.01),
  data.frame(d = 20, tok = "999", truth = 0.001),
  data.frame(d =  6, tok = "95",  truth = NA),
  data.frame(d =  7, tok = "95",  truth = NA),
  data.frame(d =  8, tok = "95",  truth = NA),
  data.frame(d =  9, tok = "95",  truth = NA),
  data.frame(d = 18, tok = "99",  truth = NA),
  data.frame(d = 19, tok = "99",  truth = NA)
)

rows <- list()
for (i in seq_len(nrow(targets))) {
  d <- targets$d[i]; tok <- targets$tok[i]; truth <- targets$truth[i]
  f <- file.path(SDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))
  if (!file.exists(f)) { cat(sprintf("  d=%d tok=%s: shipped file absent, skipped\n", d, tok)); next }
  s <- load1(f, "simul")
  P <- pool_draws(d)
  if (is.null(P)) { cat(sprintf("  d=%d: no pooled draws, skipped\n", d)); next }
  for (stat in c("average", "median")) {
    D <- if (stat == "average") P$ave else P$med
    for (shift in c(0L, 1L)) {
      lv <- level_curve(s[[stat]], D, shift)
      rows[[length(rows) + 1]] <- data.frame(
        d = d, shipped_token = tok, statistic = stat, shift = shift,
        pooled_draws = nrow(D), n_k = length(lv),
        level_median = median(lv), level_mean = mean(lv),
        level_q25 = unname(quantile(lv, 0.25)), level_q75 = unname(quantile(lv, 0.75)),
        n_zero = sum(lv == 0), truth = truth, stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, rows)

# ---- which alignment reproduces the controls? ------------------------------
cat("\n=== alignment test on the known-level files ===\n")
ctl <- res[!is.na(res$truth), ]
for (sh in c(0L, 1L)) {
  cc <- ctl[ctl$shift == sh, ]
  err <- abs(log10(pmax(cc$level_median, 1e-6)) - log10(cc$truth))
  cat(sprintf("  shift=%d : median |log10 error| over %d control cells = %.4f  (max %.4f)\n",
              sh, nrow(cc), median(err), max(err)))
}
BEST <- if (median(abs(log10(pmax(ctl$level_median[ctl$shift == 0], 1e-6)) -
                       log10(ctl$truth[ctl$shift == 0]))) <=
            median(abs(log10(pmax(ctl$level_median[ctl$shift == 1], 1e-6)) -
                       log10(ctl$truth[ctl$shift == 1])))) 0L else 1L
cat(sprintf("  -> adopting shift = %d\n", BEST))

b <- res[res$shift == BEST, ]
cat("\n=== CONTROLS (level must land on `truth`) ===\n")
cb <- b[!is.na(b$truth), ]
cb$ratio <- cb$level_median / cb$truth
print(cb[, c("d","shipped_token","statistic","truth","level_median","ratio","n_k","pooled_draws")],
      row.names = FALSE, digits = 4)
ok <- all(cb$ratio > 0.5 & cb$ratio < 2.5)
cat(sprintf("  all %d control cells within a factor 2.5 of nominal: %s\n", nrow(cb), ok))

cat("\n=== UNKNOWNS (duplicated pairs) ===\n")
ub <- b[is.na(b$truth), ]
cand <- function(d) if (d <= 9) c(0.05, 0.01) else c(0.01, 0.001)
ub$nearest <- vapply(seq_len(nrow(ub)), function(i) {
  cs <- cand(ub$d[i]); cs[which.min(abs(log10(ub$level_median[i]) - log10(cs)))] }, numeric(1))
ub$ratio_to_nearest <- ub$level_median / ub$nearest
ub$ratio_to_other <- vapply(seq_len(nrow(ub)), function(i) {
  cs <- cand(ub$d[i]); ub$level_median[i] / cs[cs != ub$nearest[i]] }, numeric(1))
print(ub[, c("d","shipped_token","statistic","level_median","level_q25","level_q75",
             "nearest","ratio_to_nearest","ratio_to_other")], row.names = FALSE, digits = 4)

cat("\n=== VERDICT ===\n")
for (d in sort(unique(ub$d))) {
  sub <- ub[ub$d == d, ]
  agree <- length(unique(sub$nearest)) == 1
  cat(sprintf("  d=%2d  level = %.4f / %.4f (ave/med)  -> alpha = %s%s\n",
              d, sub$level_median[sub$statistic == "average"],
              sub$level_median[sub$statistic == "median"],
              if (agree) sprintf("%g", sub$nearest[1]) else "AMBIGUOUS",
              if (agree) "" else "  *** statistics disagree ***"))
}

if (!ok) cat("\n*** CONTROLS FAILED -- the unknown verdicts above are VOID ***\n")
write.csv(res, OUT, row.names = FALSE)
cat(sprintf("\nwrote %s\ndone\n", OUT))
