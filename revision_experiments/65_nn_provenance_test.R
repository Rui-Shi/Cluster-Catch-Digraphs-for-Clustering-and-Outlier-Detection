#!/usr/bin/env Rscript
# 65_nn_provenance_test.R -- which alpha actually generated each shipped NN
# quantile table? Discriminator + positive control + behavioural check.
#
# WHY THIS FILE EXISTS
# --------------------
# 21 pairs of shipped NN quantile tables are identical in content (one Monte
# Carlo run saved under two level names; 62_audit_nn_duplicates.R). For nine
# of the paper's real data sets in the small-n tier the alpha behind the
# published number is therefore not readable off the filename, and the paper's
# Table tab:alpha_real asserts one anyway. 64_nn_provenance_gen.R generated
# fresh replicate clouds at both candidate levels for each affected dimension;
# this file tests the shipped content against them.
#
# THE DESIGN DECISION THIS FILE EMBODIES
# --------------------------------------
# 46_verify_wp2a_alpha.R compared a duplicated dimension's table against a
# NEIGHBOURING dimension's known pair and failed, because the per-unit-
# dimension effect (0.01734) is 9.4x the alpha effect (0.00185): the
# comparison measured dimension, not alpha. Everything here compares WITHIN a
# dimension, against a replicate cloud that supplies the missing noise scale.
#
# THREE TESTS, IN DECREASING ORDER OF TRUST
# -----------------------------------------
# TEST 1  Pointwise sign test against the replicate clouds (the commissioned
#   primary test). A lower-tail quantile is monotone in alpha, so at every
#   subsample size k the level-A cloud sits above the level-B cloud. Ask
#   whether the shipped S_k is nearer the A-cloud median or the B-cloud
#   median and count over k. Columns of the generator's draw matrix come from
#   independent samples (R/ccds/NN_Dist_Est.R:24 draws a fresh sample for
#   every subsample size), so the k are independent trials and the count is a
#   genuine binomial. k where the two cloud medians coincide carry no
#   information and are excluded from the denominator rather than counted as
#   agreement. "Inside the A range" / "inside the B range" are reported
#   separately so "consistent with both" and "consistent with neither" stay
#   visible instead of being forced into the binary.
#
# TEST 1L  Empirical-level estimator (ADDED, not commissioned -- see the
#   caveat it exists to answer). Test 1 compares a shipped quantile estimate
#   against fresh quantile estimates computed at niter = 250. A lower-tail
#   quantile estimate is biased UPWARD at small niter -- quantile(x, 0.001)
#   from 250 draws is essentially the sample minimum, whose expected rank is
#   1/251 = 0.004, not 0.001 -- and the bias is level-dependent, so at the
#   99/999 dimensions the fresh 999-cloud does not sit where the true 0.001
#   quantile sits. The shipped tables' niter is not recorded (the generation
#   scripts left in R/NN-test_quantile use n=1000, but the shipped tables are
#   length 5000, so those scripts are not what produced them), so the bias
#   cannot be matched away. Test 1L sidesteps it: pool the raw draws across
#   all replicates (8 x 250 = 2000 per (d,k)) and report the empirical level
#   of the shipped value, p_k = Fhat_k(S_k). E[p_k] is the alpha that
#   generated S_k, whatever niter either side used. This is bias-free where
#   Test 1 is not, so where the two disagree, Test 1L is the one to believe --
#   and the positive control below says whether they disagree at all.
#
# TEST 1B  POSITIVE CONTROL, mandatory. Three dimensions carry genuinely
#   distinct shipped pairs, so the true label is known from the filename:
#   d=10 and d=20 (99%/999%) control the extreme regime that d=18 and d=19
#   live in, and d=5 (90%/95%) controls the less-extreme regime that d=6..9
#   live in -- without it the 95-vs-99 verdicts would rest on a discriminator
#   validated only at 99-vs-999. Tests 1 and 1L are run on BOTH shipped files
#   at all three (12 known-answer cases, spanning alpha = 10%, 5%, 1%, 0.1%)
#   and scored. If the discriminator does not recover the right label there,
#   it has no power and no verdict is reported at the unknown dimensions. The
#   reverse sanity check -- that the fresh A and B tables at a given d are
#   themselves distinguishable by the same statistic -- is the n_informative /
#   n_clean_sep columns: if those collapse, there is nothing to detect at that
#   d and every data set there is NO POWER.
#
# TEST 2  Behavioural (secondary, and strictly noisier than Test 1). Run
#   UN-MCCD and SUN-MCCD on each data set with the shipped table and with each
#   fresh level-A and level-B replicate, and ask whether the published metric
#   sits in the A-distribution, the B-distribution, both, or neither. It is
#   noisier because the published numbers were produced BY the shipped table,
#   so it can only ask which fresh level reproduces the shipped table's
#   BEHAVIOUR -- the question Test 1 asks directly of the table itself. Where
#   the A and B runs give identical metrics (the detector is insensitive to
#   the level on that data set -- e.g. SUN-MCCD collapsing to one cluster)
#   the cell is reported NO POWER, which is a result, not a gap.
#
# TABLE SOURCING. Follows the shadow pattern documented in
# 57_waveform_d21_rerun.R: source harness.R, then wp0_mccd_methods.R, THEN
# install a get_simul() shadow that, when TABLE_STATE$path is set and
# variant == "NN", loads exactly that file, and otherwise defers
# unconditionally to the captured original. methods/, R/ccds/, harness.R and
# wp0_mccd_methods.R are not modified. Metrics always go through the
# harness's evaluate() (label polarity: Y == 1 regular, Y == 0 outlier).
# SUN-MCCD uses min.cls = 0.05, the current S_min operating point, expressed
# as a PROPORTION of n; UN-MCCD's wrapper takes no min.cls argument.
#
# OUTPUTS (revision_experiments/results/tr1/)
#   wp2a_provenance_signtest.csv     Test 1 + Test 1L, per (d, shipped file, statistic)
#   wp2a_provenance_control.csv      the same rows restricted to d=10/20, scored
#   wp2a_provenance_behavioural.csv  Test 2, per (dataset, method, table)
#   wp2a_provenance_findings.md      written by hand from these
#
# CLI
#   Rscript revision_experiments/65_nn_provenance_test.R --test1
#   Rscript revision_experiments/65_nn_provenance_test.R --test2 [dataset,csv]
#   (--test2 appends per cell and skips cells already recorded, so it can be
#    chunked by data set and resumed after an interruption.)

suppressPackageStartupMessages(library(here))

PROV_ROOT <- here::here("R/NN-test_quantile_provenance")
SHIP_DIR  <- here::here("R/NN-test_quantile")
RESDIR    <- here::here("revision_experiments/results/tr1")
REPS      <- 1:8

# d -> generation n, candidate tokens (A = less extreme = larger quantile
# value; B = more extreme), role, and the data sets served.
DIMSPEC <- list(
  list(d =  6, n = 555,  tokA = "95", tokB = "99",  role = "unknown", datasets = c("vertebral")),
  list(d =  7, n = 555,  tokA = "95", tokB = "99",  role = "unknown", datasets = c("ecoli")),
  list(d =  8, n = 555,  tokA = "95", tokB = "99",  role = "unknown", datasets = c("pima")),
  list(d =  9, n = 1013, tokA = "95", tokB = "99",  role = "unknown", datasets = c("glass", "WBC", "stamps", "shuffle")),
  list(d = 18, n = 555,  tokA = "99", tokB = "999", role = "unknown", datasets = c("lymphography")),
  list(d = 19, n = 555,  tokA = "99", tokB = "999", role = "unknown", datasets = c("hepatitis")),
  # POSITIVE CONTROLS: shipped pairs that are genuinely distinct, so the true
  # label is known from the filename. d=10 and d=20 control the 99/999 regime
  # (d=18, d=19). d=5 controls the LESS EXTREME regime (d=6..9): its shipped
  # 90%/95% pair is distinct, and although the pair is 90/95 rather than
  # 95/99, it is the only known-label case available below d=10 and it tests
  # the same thing -- whether the estimator reads the right alpha when both
  # candidates are in the easy, well-estimated part of the tail.
  list(d =  5, n = 555,  tokA = "90", tokB = "95",  role = "control", datasets = character(0)),
  list(d = 10, n = 555,  tokA = "99", tokB = "999", role = "control", datasets = character(0)),
  list(d = 20, n = 555,  tokA = "99", tokB = "999", role = "control", datasets = character(0))
)
dimspec_of <- function(d) Filter(function(s) s$d == d, DIMSPEC)[[1]]
alpha_of_tok <- function(tok) 1 - as.numeric(paste0("0.", tok))

load_simul <- function(path) {
  if (!file.exists(path)) return(NULL)
  e <- new.env(); load(path, envir = e); get("simul", envir = e)
}
fresh_tab  <- function(r, d, tok) load_simul(file.path(PROV_ROOT, sprintf("rep%02d", r),
                                                       sprintf("NN-test-simul_%dd_%s%%.RData", d, tok)))
ship_tab   <- function(d, tok)    load_simul(file.path(SHIP_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok)))
fresh_draws <- function(r, d, n) {
  p <- file.path(PROV_ROOT, sprintf("rep%02d", r), sprintf("NN-draws_%dd_n%d.RData", d, n))
  if (!file.exists(p)) return(NULL)
  e <- new.env(); load(p, envir = e)
  list(ave = get("NN.dist.ave.mat", envir = e), med = get("NN.dist.med.mat", envir = e))
}

# ===========================================================================
# TEST 1 + TEST 1L
# ===========================================================================

#' @param S numeric shipped vector (length >= K)
#' @param A,B  matrices reps x K of the fresh level-A / level-B tables
#' @param pool matrix (reps*niter) x K of the raw draws
discriminate <- function(S, A, B, pool, alphaA, alphaB) {
  K <- ncol(A)
  idx <- 2:K                       # position 1 is the generator's structural 0
  Sk    <- S[idx]
  medA  <- apply(A[, idx, drop = FALSE], 2, median)
  medB  <- apply(B[, idx, drop = FALSE], 2, median)
  minA  <- apply(A[, idx, drop = FALSE], 2, min); maxA <- apply(A[, idx, drop = FALSE], 2, max)
  minB  <- apply(B[, idx, drop = FALSE], 2, min); maxB <- apply(B[, idx, drop = FALSE], 2, max)

  informative <- medA != medB
  clean_sep   <- minA > maxB                     # clouds do not overlap at all
  nearer_A    <- abs(Sk - medA) < abs(Sk - medB)
  tie         <- abs(Sk - medA) == abs(Sk - medB)

  nI  <- sum(informative & !tie)
  nA  <- sum(informative & !tie & nearer_A)
  nC  <- sum(clean_sep)
  nCI <- sum(clean_sep & !tie)
  nCA <- sum(clean_sep & !tie & nearer_A)

  bt  <- if (nI > 0) stats::binom.test(max(nA, nI - nA), nI, 0.5)$p.value else NA_real_
  btC <- if (nCI > 0) stats::binom.test(max(nCA, nCI - nCA), nCI, 0.5)$p.value else NA_real_

  inA <- Sk >= minA & Sk <= maxA
  inB <- Sk >= minB & Sk <= maxB

  # ---- Test 1L: empirical level of the shipped value in the pooled draws ----
  M   <- nrow(pool)
  phat <- vapply(seq_along(idx), function(j) mean(pool[, idx[j]] < Sk[j]), 0.0)
  # p = 0 means S is below all M pooled draws; the informative floor is 1/(2M).
  phat_f <- pmax(phat, 1 / (2 * M))
  vote_A <- abs(log(phat_f) - log(alphaA)) < abs(log(phat_f) - log(alphaB))
  nLA <- sum(vote_A); nL <- length(vote_A)
  btL <- if (nL > 0) stats::binom.test(max(nLA, nL - nLA), nL, 0.5)$p.value else NA_real_

  list(
    n_k = length(idx),
    n_informative = nI, n_favour_A = nA, frac_A = if (nI > 0) nA / nI else NA_real_, p_sign = bt,
    n_clean_sep = nC, n_clean_informative = nCI, frac_A_clean = if (nCI > 0) nCA / nCI else NA_real_,
    p_sign_clean = btC,
    n_inside_A_only = sum(inA & !inB), n_inside_B_only = sum(inB & !inA),
    n_inside_both = sum(inA & inB), n_inside_neither = sum(!inA & !inB),
    mean_sep = mean(medA - medB), median_sep = median(medA - medB),
    mean_cloud_sd_A = mean(apply(A[, idx, drop = FALSE], 2, sd)),
    mean_cloud_sd_B = mean(apply(B[, idx, drop = FALSE], 2, sd)),
    level_mean = mean(phat), level_median = median(phat),
    level_geomean = exp(mean(log(phat_f))), n_level_zero = sum(phat == 0),
    n_level_vote_A = nLA, frac_level_A = nLA / nL, p_level = btL,
    pool_M = M
  )
}

verdict_of <- function(frac_A, p, tokA, tokB, p_cut = 1e-3) {
  if (is.na(frac_A) || is.na(p)) return("no power")
  if (p >= p_cut) return("inconclusive")
  if (frac_A > 0.5) tokA else tokB
}

run_test1 <- function() {
  out <- list()
  for (spec in DIMSPEC) {
    d <- spec$d; n <- spec$n; tokA <- spec$tokA; tokB <- spec$tokB
    aA <- alpha_of_tok(tokA); aB <- alpha_of_tok(tokB)

    A_ave <- do.call(rbind, lapply(REPS, function(r) fresh_tab(r, d, tokA)$average))
    B_ave <- do.call(rbind, lapply(REPS, function(r) fresh_tab(r, d, tokB)$average))
    A_med <- do.call(rbind, lapply(REPS, function(r) fresh_tab(r, d, tokA)$median))
    B_med <- do.call(rbind, lapply(REPS, function(r) fresh_tab(r, d, tokB)$median))
    if (is.null(A_ave) || nrow(A_ave) != length(REPS)) stop("missing fresh tables at d=", d)

    dl <- lapply(REPS, function(r) fresh_draws(r, d, n))
    pool_ave <- do.call(rbind, lapply(dl, `[[`, "ave"))
    pool_med <- do.call(rbind, lapply(dl, `[[`, "med"))

    ship_toks <- if (identical(spec$role, "control")) {
      # A control is only a control if its shipped pair really is distinct.
      # Assert it rather than trusting 62's global count.
      sa <- ship_tab(d, tokA); sb <- ship_tab(d, tokB)
      if (identical(sa$average, sb$average) && identical(sa$median, sb$median))
        stop(sprintf("d=%d is marked 'control' but its shipped pair IS identical -- it has no known label", d))
      c(tokA, tokB)
    } else {
      # unknown dimensions: the two shipped files are identical in content, so
      # one evaluation covers both. Verified, not assumed, below.
      sa <- ship_tab(d, tokA); sb <- ship_tab(d, tokB)
      if (!identical(sa$average, sb$average) || !identical(sa$median, sb$median))
        stop(sprintf("d=%d: shipped pair is NOT identical -- DIMSPEC role is wrong", d))
      tokA
    }

    for (stok in ship_toks) {
      S <- ship_tab(d, stok)
      for (stat in c("average", "median")) {
        r <- discriminate(S[[stat]],
                          if (stat == "average") A_ave else A_med,
                          if (stat == "average") B_ave else B_med,
                          if (stat == "average") pool_ave else pool_med,
                          aA, aB)
        out[[length(out) + 1]] <- data.frame(
          d = d, role = spec$role, n_gen = n, niter_per_rep = nrow(pool_ave) / length(REPS),
          n_reps = length(REPS), statistic = stat,
          shipped_file = sprintf("NN-test-simul_%dd_%s%%.RData", d, stok),
          shipped_token = stok,
          shipped_pair_identical = !identical(spec$role, "control"),
          tokA = tokA, tokB = tokB, alphaA = aA, alphaB = aB,
          as.data.frame(r),
          verdict_sign  = verdict_of(r$frac_A, r$p_sign, tokA, tokB),
          verdict_level = verdict_of(r$frac_level_A, r$p_level, tokA, tokB),
          verdict_level_direct = if (abs(log(max(r$level_mean, 1e-6)) - log(aA)) <
                                     abs(log(max(r$level_mean, 1e-6)) - log(aB))) tokA else tokB,
          truth = if (identical(spec$role, "control")) stok else NA_character_,
          datasets = paste(spec$datasets, collapse = ";"),
          stringsAsFactors = FALSE)
      }
    }
    cat(sprintf("test1 d=%d done\n", d)); flush.console()
  }
  res <- do.call(rbind, out)
  dir.create(RESDIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(res, file.path(RESDIR, "wp2a_provenance_signtest.csv"), row.names = FALSE)

  ctl <- res[res$role == "control", ]
  ctl$sign_correct  <- ctl$verdict_sign  == ctl$truth
  ctl$level_correct <- ctl$verdict_level == ctl$truth
  ctl$level_direct_correct <- ctl$verdict_level_direct == ctl$truth
  write.csv(ctl, file.path(RESDIR, "wp2a_provenance_control.csv"), row.names = FALSE)

  cat("\n==== POSITIVE CONTROL (d=5, 10, 20: true label known from filename) ====\n")
  print(ctl[, c("d", "statistic", "shipped_token", "truth", "n_informative", "frac_A", "p_sign",
                "verdict_sign", "sign_correct", "level_mean", "verdict_level_direct",
                "level_direct_correct")], row.names = FALSE, digits = 4)
  cat("\n==== UNKNOWN DIMENSIONS ====\n")
  unk <- res[res$role == "unknown", ]
  print(unk[, c("d", "statistic", "tokA", "tokB", "n_k", "n_informative", "n_clean_sep",
                "frac_A", "p_sign", "verdict_sign", "level_mean", "verdict_level_direct",
                "n_inside_A_only", "n_inside_B_only", "n_inside_both", "n_inside_neither")],
        row.names = FALSE, digits = 4)
  invisible(res)
}

# ===========================================================================
# TEST 2 -- behavioural
# ===========================================================================
BEHAV_CSV <- file.path(RESDIR, "wp2a_provenance_behavioural.csv")
MIN_CLS_SUN <- 0.05     # current S_min operating point, a PROPORTION of n
DATASET_D <- c(hepatitis = 19, lymphography = 18, glass = 9, WBC = 9, vertebral = 6,
               ecoli = 7, stamps = 9, pima = 8, shuffle = 9)

run_test2 <- function(which_datasets) {
  source(here::here("revision_experiments", "harness.R"))
  source(here::here("revision_experiments", "wp0_mccd_methods.R"))

  orig_get_simul <- get_simul                       # capture BEFORE shadowing
  TABLE_STATE <- new.env(parent = emptyenv()); TABLE_STATE$path <- NA_character_
  PROV <- new.env(parent = emptyenv()); PROV$file <- NA_character_; PROV$len <- NA_integer_

  get_simul_shadow <- function(variant = c("RK", "NN"), d, quant = NULL) {
    p <- TABLE_STATE$path
    if (identical(variant, "NN") && !is.na(p)) {
      if (!file.exists(p)) stop("get_simul_shadow(): table missing: ", p)
      e <- new.env(); load(p, envir = e)
      if (!exists("simul", envir = e)) stop("get_simul_shadow(): no 'simul' in ", p)
      res <- list(simul = get("simul", envir = e),
                  quant = as.numeric(paste0("0.", sub(".*_([0-9]+)%\\.RData$", "\\1", basename(p)))),
                  quant_label = sub(".*_([0-9]+)%\\.RData$", "\\1", basename(p)), file = p)
    } else {
      res <- orig_get_simul(variant, d, quant)
    }
    PROV$file <- res$file; PROV$len <- length(res$simul$average)
    res
  }
  assign("get_simul", get_simul_shadow, envir = globalenv())

  pub <- read.csv(file.path(RESDIR, "final_comparison.csv"), stringsAsFactors = FALSE)

  for (ds in which_datasets) {
    d <- DATASET_D[[ds]]; spec <- dimspec_of(d)
    dat <- load_real_dataset(ds)
    stopifnot(dat$d == d)
    cat(sprintf("\n--- %s (n=%d d=%d n_out=%d) ---\n", ds, dat$n, dat$d, sum(dat$Y == 0)))

    tables <- list(list(id = "shipped", token = spec$tokA, rep = NA_integer_,
                        path = file.path(SHIP_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, spec$tokA))))
    for (tok in c(spec$tokA, spec$tokB)) for (r in REPS)
      tables[[length(tables) + 1]] <- list(id = sprintf("fresh_%s_rep%02d", tok, r), token = tok, rep = r,
                                           path = file.path(PROV_ROOT, sprintf("rep%02d", r),
                                                            sprintf("NN-test-simul_%dd_%s%%.RData", d, tok)))

    for (method in c("UN-MCCD", "SUN-MCCD")) {
      pr <- pub[pub$dataset == ds & pub$method == method, ]
      for (tb in tables) {
        keys <- c(dataset = ds, method = method, table_id = tb$id)
        if (isTRUE(has_result(BEHAV_CSV, keys))) next
        TABLE_STATE$path <- tb$path
        t0 <- Sys.time()
        out <- tryCatch({
          res <- if (identical(method, "SUN-MCCD"))
            METHOD_REGISTRY[[method]](X = dat$X, d = d, Y = dat$Y, quant = tb$token, min.cls = MIN_CLS_SUN)
          else
            METHOD_REGISTRY[[method]](X = dat$X, d = d, Y = dat$Y, quant = tb$token)
          list(m = evaluate(dat$Y, res$score, REAL_DATA_THRESHOLDS[[method]]), res = res,
               status = "ok", note = NA_character_)
        }, error = function(e) list(m = setNames(rep(NA_real_, 4), c("TPR", "TNR", "BA", "F2")),
                                    res = NULL, status = "error", note = conditionMessage(e)))
        TABLE_STATE$path <- NA_character_
        ncls <- if (is.null(out$res)) NA_integer_ else length(unique(out$res$cluster[!is.na(out$res$cluster)]))
        nflag <- if (is.null(out$res)) NA_integer_ else sum(out$res$score > REAL_DATA_THRESHOLDS[[method]])
        append_result(BEHAV_CSV, list(
          dataset = ds, method = method, table_id = tb$id, table_source = if (tb$id == "shipped") "shipped" else "fresh",
          token = tb$token, alpha = alpha_of_tok(tb$token), replicate = tb$rep, n = dat$n, d = d,
          TPR = unname(out$m[["TPR"]]), TNR = unname(out$m[["TNR"]]),
          BA = unname(out$m[["BA"]]), F2 = unname(out$m[["F2"]]),
          published_TPR = if (nrow(pr)) pr$TPR[1] else NA_real_,
          published_TNR = if (nrow(pr)) pr$TNR[1] else NA_real_,
          published_BA  = if (nrow(pr)) pr$BA[1]  else NA_real_,
          published_F2  = if (nrow(pr)) pr$F2[1]  else NA_real_,
          n_flagged = nflag, n_clusters = ncls,
          table_file = PROV$file, table_len = PROV$len,
          elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
          status = out$status, note = out$note, timestamp = format(Sys.time())))
        cat(sprintf("  %-9s %-18s BA=%.4f F2=%.4f ncls=%s %s\n", method, tb$id,
                    out$m[["BA"]], out$m[["F2"]], ncls, out$status)); flush.console()
      }
    }
  }
  cat("\nTEST2_COMPLETE\n")
}

# ---------------------------------------------------------------------------
# Test 2 reduction. For each (dataset, method) ask, per metric, whether the
# PUBLISHED value falls inside the level-A replicate range, the level-B range,
# both, or neither -- and, first, whether the two ranges are the same set, in
# which case the cell has NO POWER and no verdict may be read off it. The
# published numbers are rounded to 3 dp, so ranges are widened by 5e-4 before
# the containment test; anything tighter would reject on rounding alone.
# ---------------------------------------------------------------------------
summarize_test2 <- function() {
  df <- read.csv(BEHAV_CSV, stringsAsFactors = FALSE)
  df <- df[df$status == "ok", ]
  TOL <- 5e-4
  rows <- list()
  for (ds in unique(df$dataset)) for (mth in unique(df$method[df$dataset == ds])) {
    sub <- df[df$dataset == ds & df$method == mth, ]
    spec <- dimspec_of(sub$d[1])
    A <- sub[sub$table_source == "fresh" & sub$token == spec$tokA, ]
    B <- sub[sub$table_source == "fresh" & sub$token == spec$tokB, ]
    SH <- sub[sub$table_source == "shipped", ]
    for (met in c("BA", "F2", "TPR", "TNR")) {
      a <- A[[met]]; b <- B[[met]]; pubv <- A[[paste0("published_", met)]][1]
      coincide <- identical(sort(round(a, 10)), sort(round(b, 10)))
      inA <- pubv >= min(a) - TOL && pubv <= max(a) + TOL
      inB <- pubv >= min(b) - TOL && pubv <= max(b) + TOL
      rows[[length(rows) + 1]] <- data.frame(
        dataset = ds, method = mth, d = sub$d[1], metric = met,
        tokA = spec$tokA, tokB = spec$tokB, published = pubv,
        shipped = if (nrow(SH)) SH[[met]][1] else NA_real_,
        shipped_matches_published = if (nrow(SH)) abs(SH[[met]][1] - pubv) <= TOL else NA,
        A_min = min(a), A_max = max(a), A_median = median(a), A_n_distinct = length(unique(round(a, 10))),
        B_min = min(b), B_max = max(b), B_median = median(b), B_n_distinct = length(unique(round(b, 10))),
        distributions_coincide = coincide,
        published_in_A = inA, published_in_B = inB,
        verdict = if (coincide) "NO POWER (A and B give the same values)"
                  else if (inA && inB) "consistent with both"
                  else if (inA) spec$tokA else if (inB) spec$tokB else "consistent with neither",
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, rows)
  write.csv(res, file.path(RESDIR, "wp2a_provenance_behavioural_summary.csv"), row.names = FALSE)
  print(res[res$metric %in% c("BA", "F2"),
            c("dataset", "method", "d", "metric", "published", "shipped", "shipped_matches_published",
              "A_min", "A_max", "B_min", "B_max", "distributions_coincide", "verdict")],
        row.names = FALSE, digits = 4)
  invisible(res)
}

# ===========================================================================
args <- commandArgs(trailingOnly = TRUE)
if (any(args == "--summarize2")) {
  summarize_test2()
} else if (any(args == "--test1")) {
  run_test1()
} else if (any(args == "--test2")) {
  rest <- args[args != "--test2"]
  dss <- if (length(rest)) strsplit(rest[1], ",")[[1]] else names(DATASET_D)
  run_test2(dss)
} else {
  run_test1()
  run_test2(names(DATASET_D))
}
