#!/usr/bin/env Rscript
# 66_alpha_claim_audit.R -- mechanical cross-check of every alpha-level claim
# the manuscript makes (main text L907-910 + Table tab:alpha_real L1098-1111;
# SupplementaryMaterial.tex L718-740) against what the code actually loads and
# whether the loaded file's content is distinguishable from its sibling at the
# same dimension.
#
# Context: 62_audit_nn_duplicates.R established that 21 (d, token-pair)
# combinations in R/NN-test_quantile/ are byte-identical files saved under two
# alpha labels -- one Monte-Carlo run, two names. That fact alone does not
# tell us which manuscript sentences it invalidates: it depends on (a) which
# token each detector actually resolves to at each dimension used by the real
# data sets and the simulation grid, and (b) whether the manuscript claims a
# level for a file that turns out to be a duplicate, or claims two DIFFERENT
# levels for UN-MCCD/SUN-MCCD that turn out to load the SAME file.
#
# This script answers that mechanically instead of by inspection:
#   1. Recomputes, for every d used by the 16 real data sets and the
#      simulation grid {2,3,5,10,20,50,100}, the exact filename token each of
#      the four detectors resolves to. The resolvers are transcribed verbatim
#      from the live code (cited by file:line below) rather than re-derived,
#      so a change to the real resolver is a change this script's assertions
#      would need to be re-synced against -- the citations make that visible.
#   2. For every (variant, d) with 2+ candidate token files on disk, loads
#      them and compares content directly (NN: $average/$median; RK: $quan,
#      the sub-list actually indexed by the detector -- see
#      wp0_mccd_methods.R:162-166). File size is used as a fast pre-filter
#      (unequal size already proves inequality without loading), and full
#      content loading is the deciding check whenever sizes match.
#   3. Encodes the manuscript's claimed levels (Table tab:alpha_real's three
#      row buckets, and the two SM piecewise schedules) as literal lookup
#      tables transcribed from the .tex sources, and compares them against
#      what step 1 resolved.
#   4. Classifies every claim as VERIFIED / UNVERIFIABLE / FALSE using the
#      rule stated in the task brief: VERIFIED = distinct file, code loads it;
#      UNVERIFIABLE = file exists but is a content-duplicate of another
#      labeled file at the same d, so its true generating alpha is unknown;
#      FALSE = the code demonstrably loads something other than the claim
#      (label mismatch, missing file, or a claimed UN/SUN distinction that
#      collapses to one file).
#
# Read-only: only reads R/RK-test_quantile/, R/NN-test_quantile/, and
# results/tr1/dataset_inventory.csv. Writes exactly one file:
#   results/tr1/wp2a_alpha_claim_audit.csv

suppressMessages(library(here))

RK_DIR <- here::here("R/RK-test_quantile")
NN_DIR <- here::here("R/NN-test_quantile")
INV_CSV <- here::here("revision_experiments/results/tr1/dataset_inventory.csv")
OUT_CSV <- here::here("revision_experiments/results/tr1/wp2a_alpha_claim_audit.csv")

# ---------------------------------------------------------------------------
# 1. Resolver functions, transcribed verbatim from the live code (not
#    sourced, to avoid pulling in the full detector dependency chain for what
#    is meant to be a purely documentary audit -- sourcing harness.R crashed
#    the R session on this machine with a segfault, confirmed separately).
# ---------------------------------------------------------------------------

# revision_experiments/wp0_mccd_methods.R:177
#   rk_quant_label_paper <- function(d) if (d < 10) "99" else "999"
# (main text line ~880: alpha=1% for d<10, alpha=0.1% for d>=10.)
rk_quant_label_paper <- function(d) if (d < 10) "99" else "999"

# revision_experiments/harness.R:265-271 (nn_quant_for_d), the shared step
# function both NN wrappers below build on.
nn_quant_for_d <- function(d) {
  if (d <= 2) "85"
  else if (d <= 4) "90"
  else if (d <= 9) "95"
  else if (d <= 19) "99"
  else "999"
}

# revision_experiments/wp0_mccd_methods.R:185
#   nn_quant_label_paper_UN <- function(d) nn_quant_for_d(d)
nn_quant_label_paper_UN <- function(d) nn_quant_for_d(d)

# revision_experiments/wp0_mccd_methods.R:190
#   nn_quant_label_paper_SUN <- function(d) if (d < 10) nn_quant_for_d(d) else "999"
nn_quant_label_paper_SUN <- function(d) if (d < 10) nn_quant_for_d(d) else "999"

tok_to_pct <- function(tok) {
  # Filename tokens are the QUANTILE level with the decimal point removed:
  # "99" = the 99% quantile, "999" = the 99.9% quantile, "95" = 95%, "85" =
  # 85%. The manuscript reports the complementary SIGNIFICANCE LEVEL alpha =
  # 100% - quantile (e.g. the 99% quantile table is used at alpha = 1%, the
  # 99.9% quantile table at alpha = 0.1%) -- this is what Table tab:alpha_real
  # and the SM piecewise schedules both state. A 2-digit token is read as
  # XX.0%; a 3-digit token as XX.X% (one implied decimal place).
  quant_pct <- if (nchar(tok) <= 2) as.numeric(tok) else as.numeric(tok) / 10
  alpha_pct <- 100 - quant_pct
  paste0(format(alpha_pct, trim = TRUE), "%")
}

rk_path <- function(d, tok) file.path(RK_DIR, sprintf("RK-test-simul_%dd_%s%%.RData", d, tok))
nn_path <- function(d, tok) file.path(NN_DIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))

# ---------------------------------------------------------------------------
# 2. Manuscript claims, transcribed literally from the .tex sources.
# ---------------------------------------------------------------------------

# Table tab:alpha_real, CCD_OutlierDetection_Neurocomputing.tex L1098-1111.
# Row bucketing exactly as printed (dataset lists, d range, RK / UN / SUN %).
alpha_real_rows <- list(
  list(datasets = c("wilt","vertebral","thyroid","ecoli","pima","glass","WBC","Shuttle","stamps"),
       d_lo = 5, d_hi = 9, rk = "1%", un = "5%", sun = "5%"),
  list(datasets = c("pageblocks","vowels","PenDigits","lymphography","hepatitis"),
       d_lo = 10, d_hi = 19, rk = "0.1%", un = "1%", sun = "0.1%"),
  list(datasets = c("waveform","WDBC"),
       d_lo = 20, d_hi = Inf, rk = "0.1%", un = "0.1%", sun = "0.1%")
)

# SupplementaryMaterial.tex L718-740, the two piecewise schedules for the
# simulation grid d in {2,3,5,10,20,50,100}, plus the RK rule at L719-720.
sim_rk_claim <- function(d) if (d < 10) "1%" else "0.1%"
sim_alpha_UN_claim <- c(`2` = "15%", `3` = "10%", `5` = "5%", `10` = "1%",
                         `20` = "0.1%", `50` = "0.1%", `100` = "0.1%")
sim_alpha_SUN_claim <- c(`2` = "15%", `3` = "10%", `5` = "5%", `10` = "0.1%",
                          `20` = "0.1%", `50` = "0.1%", `100` = "0.1%")

# ---------------------------------------------------------------------------
# 3. Content comparators.
# ---------------------------------------------------------------------------

load_simul <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

# Returns "identical" / "distinct" / "single" (no sibling at this d) /
# "missing" (the token file itself is absent).
# `all_tokens_at_d`: character vector of every token that has a file at this
# d (for the given variant), used to find the sibling(s) to compare against.
dup_check <- function(variant, d, tok, all_tokens_at_d) {
  siblings <- setdiff(all_tokens_at_d, tok)
  target_path <- if (variant == "NN") nn_path(d, tok) else rk_path(d, tok)
  if (!file.exists(target_path)) return(list(status = "missing", partner = NA_character_))
  if (length(siblings) == 0) return(list(status = "single", partner = NA_character_))

  target_size <- file.size(target_path)
  for (stok in siblings) {
    spath <- if (variant == "NN") nn_path(d, stok) else rk_path(d, stok)
    if (!file.exists(spath)) next
    if (file.size(spath) != target_size) next  # different size => provably distinct, skip load
    # sizes match: only a full content load settles it
    a <- load_simul(target_path)
    b <- load_simul(spath)
    same <- if (variant == "NN") {
      identical(a$average, b$average) && identical(a$median, b$median)
    } else {
      identical(a$quan, b$quan)
    }
    if (same) return(list(status = "identical", partner = stok))
  }
  list(status = "distinct", partner = NA_character_)
}

# ---------------------------------------------------------------------------
# 4. Part A -- the 16 real data sets.
# ---------------------------------------------------------------------------

inv <- read.csv(INV_CSV, stringsAsFactors = FALSE)
stopifnot(nrow(inv) == 16)
inv$display_name <- ifelse(inv$dataset == "shuffle", "Shuttle", inv$dataset)

# every d actually present among the RK / NN directories, restricted to the
# dims we need (real-data dims + simulation dims), so the "all tokens at this
# d" lookups below are exhaustive and not hand-curated.
sim_dims <- c(2, 3, 5, 10, 20, 50, 100)
real_dims <- unique(inv$d)
all_dims_needed <- sort(unique(c(real_dims, sim_dims)))

list_tokens_at_d <- function(variant, d) {
  dir <- if (variant == "NN") NN_DIR else RK_DIR
  pat <- sprintf("^%s-test-simul_%dd_([0-9]+)%%%%\\.RData$",
                 if (variant == "NN") "NN" else "RK", d)
  fs <- list.files(dir, full.names = FALSE)
  fs <- fs[grepl(sprintf("^%s-test-simul_%dd_[0-9]+%%\\.RData$",
                          if (variant == "NN") "NN" else "RK", d), fs)]
  sub(sprintf("^%s-test-simul_%dd_([0-9]+)%%\\.RData$",
              if (variant == "NN") "NN" else "RK", d), "\\1", fs)
}

# lookup for the manuscript's claimed row, given d
lookup_row <- function(d) {
  for (r in alpha_real_rows) if (d >= r$d_lo && d <= r$d_hi) return(r)
  NULL
}

real_rows <- list()
for (i in seq_len(nrow(inv))) {
  ds <- inv$display_name[i]; d <- inv$d[i]; n <- inv$n[i]
  row <- lookup_row(d)
  claimed_rk <- if (!is.null(row) && ds %in% row$datasets) row$rk else NA
  claimed_un <- if (!is.null(row) && ds %in% row$datasets) row$un else NA
  claimed_sun <- if (!is.null(row) && ds %in% row$datasets) row$sun else NA
  row_ok <- !is.null(row) && ds %in% row$datasets  # dataset actually falls in the row it's printed under

  # -- RK --
  rk_tok <- rk_quant_label_paper(d)
  rk_tokens_here <- list_tokens_at_d("RK", d)
  rk_dup <- dup_check("RK", d, rk_tok, rk_tokens_here)
  rk_computed_pct <- tok_to_pct(rk_tok)

  # -- NN, UN-MCCD --
  un_tok <- nn_quant_label_paper_UN(d)
  nn_tokens_here <- list_tokens_at_d("NN", d)
  un_dup <- dup_check("NN", d, un_tok, nn_tokens_here)
  un_computed_pct <- tok_to_pct(un_tok)

  # -- NN, SUN-MCCD --
  sun_tok <- nn_quant_label_paper_SUN(d)
  sun_dup <- dup_check("NN", d, sun_tok, nn_tokens_here)
  sun_computed_pct <- tok_to_pct(sun_tok)

  # UN vs SUN: same file trivially, or different files -- and if different
  # files, are they nonetheless content-identical? (case (b) in the brief)
  if (un_tok == sun_tok) {
    un_sun_relationship <- "same token (paper does not claim a distinction)"
  } else {
    up <- nn_path(d, un_tok); sp <- nn_path(d, sun_tok)
    if (file.exists(up) && file.exists(sp) && file.size(up) == file.size(sp)) {
      a <- load_simul(up); b <- load_simul(sp)
      same <- identical(a$average, b$average) && identical(a$median, b$median)
      un_sun_relationship <- if (same)
        "DIFFERENT tokens claimed, IDENTICAL content (case b: claimed distinction did not exist)"
      else "different tokens, content genuinely differs (distinction real)"
    } else if (file.exists(up) && file.exists(sp)) {
      un_sun_relationship <- "different tokens, different file sizes (distinction real)"
    } else {
      un_sun_relationship <- "cannot compare: a token file is missing"
    }
  }

  un_sun_case_b <- grepl("^DIFFERENT tokens claimed, IDENTICAL", un_sun_relationship)

  verdict_of <- function(computed_pct, claimed_pct, dupstat) {
    if (is.na(claimed_pct)) return("NO CLAIM (dataset not placed in any table row)")
    if (dupstat$status == "missing") return("FALSE: claimed level has no file on disk")
    if (computed_pct != claimed_pct) return("FALSE: code resolves a different token than claimed")
    if (dupstat$status == "identical") return(sprintf("UNVERIFIABLE: file is content-duplicate of token %s%%", dupstat$partner))
    if (dupstat$status == "single") return("VERIFIED (single file at this d, no sibling to cross-check against)")
    "VERIFIED (distinct file, content differs from every sibling token at this d)"
  }

  real_rows[[length(real_rows) + 1]] <- data.frame(
    section = "real_data", dataset = ds, n = n, d = d,
    row_claimed_d_range = if (!is.null(row)) (if (is.infinite(row$d_hi)) sprintf(">=%s", row$d_lo) else sprintf("%s-%s", row$d_lo, row$d_hi)) else NA,
    dataset_in_correct_row = row_ok,
    rk_token = rk_tok, rk_claimed_pct = claimed_rk, rk_computed_pct = rk_computed_pct,
    rk_file = basename(rk_path(d, rk_tok)), rk_dup_status = rk_dup$status, rk_dup_partner = rk_dup$partner,
    verdict_rk = verdict_of(rk_computed_pct, claimed_rk, rk_dup),
    un_token = un_tok, un_claimed_pct = claimed_un, un_computed_pct = un_computed_pct,
    un_file = basename(nn_path(d, un_tok)), un_dup_status = un_dup$status, un_dup_partner = un_dup$partner,
    verdict_un = verdict_of(un_computed_pct, claimed_un, un_dup),
    sun_token = sun_tok, sun_claimed_pct = claimed_sun, sun_computed_pct = sun_computed_pct,
    sun_file = basename(nn_path(d, sun_tok)), sun_dup_status = sun_dup$status, sun_dup_partner = sun_dup$partner,
    verdict_sun = verdict_of(sun_computed_pct, claimed_sun, sun_dup),
    un_sun_relationship = un_sun_relationship,
    un_sun_distinction_is_false = un_sun_case_b,
    stringsAsFactors = FALSE
  )
}
real_df <- do.call(rbind, real_rows)

# ---------------------------------------------------------------------------
# 5. Part B -- the simulation grid, d in {2,3,5,10,20,50,100}.
# ---------------------------------------------------------------------------

sim_rows <- list()
for (d in sim_dims) {
  rk_tok <- rk_quant_label_paper(d)
  rk_tokens_here <- list_tokens_at_d("RK", d)
  rk_dup <- dup_check("RK", d, rk_tok, rk_tokens_here)
  rk_claim <- sim_rk_claim(d)
  rk_computed <- tok_to_pct(rk_tok)
  rk_verdict <- if (rk_dup$status == "missing") "FALSE: claimed level has no file on disk" else
    if (rk_computed != rk_claim) "FALSE: code resolves a different token than claimed" else
    if (rk_dup$status == "identical") sprintf("UNVERIFIABLE: duplicate of %s%%", rk_dup$partner) else
    if (rk_dup$status == "single") "VERIFIED (single file, no sibling)" else "VERIFIED (distinct)"

  un_tok <- nn_quant_label_paper_UN(d)
  nn_tokens_here <- list_tokens_at_d("NN", d)
  un_dup <- dup_check("NN", d, un_tok, nn_tokens_here)
  un_claim <- unname(sim_alpha_UN_claim[as.character(d)])
  un_computed <- tok_to_pct(un_tok)
  un_verdict <- if (un_dup$status == "missing") "FALSE: claimed level has no file on disk" else
    if (un_computed != un_claim) "FALSE: code resolves a different token than claimed" else
    if (un_dup$status == "identical") sprintf("UNVERIFIABLE: duplicate of %s%%", un_dup$partner) else
    if (un_dup$status == "single") "VERIFIED (single file, no sibling)" else "VERIFIED (distinct)"

  sun_tok <- nn_quant_label_paper_SUN(d)
  sun_dup <- dup_check("NN", d, sun_tok, nn_tokens_here)
  sun_claim <- unname(sim_alpha_SUN_claim[as.character(d)])
  sun_computed <- tok_to_pct(sun_tok)
  sun_verdict <- if (sun_dup$status == "missing") "FALSE: claimed level has no file on disk" else
    if (sun_computed != sun_claim) "FALSE: code resolves a different token than claimed" else
    if (sun_dup$status == "identical") sprintf("UNVERIFIABLE: duplicate of %s%%", sun_dup$partner) else
    if (sun_dup$status == "single") "VERIFIED (single file, no sibling)" else "VERIFIED (distinct)"

  if (un_tok == sun_tok) {
    un_sun_rel <- "same token (schedules coincide at this d)"
  } else {
    up <- nn_path(d, un_tok); sp <- nn_path(d, sun_tok)
    if (file.exists(up) && file.exists(sp) && file.size(up) == file.size(sp)) {
      a <- load_simul(up); b <- load_simul(sp)
      same <- identical(a$average, b$average) && identical(a$median, b$median)
      un_sun_rel <- if (same) "DIFFERENT tokens claimed, IDENTICAL content (case b)" else "content genuinely differs"
    } else if (file.exists(up) && file.exists(sp)) {
      un_sun_rel <- "different file sizes (distinction real)"
    } else un_sun_rel <- "cannot compare: a token file is missing"
  }

  sim_rows[[length(sim_rows) + 1]] <- data.frame(
    section = "simulation", dataset = NA, n = NA, d = d,
    row_claimed_d_range = NA, dataset_in_correct_row = NA,
    rk_token = rk_tok, rk_claimed_pct = rk_claim, rk_computed_pct = rk_computed,
    rk_file = basename(rk_path(d, rk_tok)), rk_dup_status = rk_dup$status, rk_dup_partner = rk_dup$partner,
    verdict_rk = rk_verdict,
    un_token = un_tok, un_claimed_pct = un_claim, un_computed_pct = un_computed,
    un_file = basename(nn_path(d, un_tok)), un_dup_status = un_dup$status, un_dup_partner = un_dup$partner,
    verdict_un = un_verdict,
    sun_token = sun_tok, sun_claimed_pct = sun_claim, sun_computed_pct = sun_computed,
    sun_file = basename(nn_path(d, sun_tok)), sun_dup_status = sun_dup$status, sun_dup_partner = sun_dup$partner,
    verdict_sun = sun_verdict,
    un_sun_relationship = un_sun_rel,
    un_sun_distinction_is_false = grepl("^DIFFERENT tokens claimed, IDENTICAL", un_sun_rel),
    stringsAsFactors = FALSE
  )
}
sim_df <- do.call(rbind, sim_rows)

# ---------------------------------------------------------------------------
# 6. Write + console summary.
# ---------------------------------------------------------------------------

out <- rbind(real_df, sim_df)
write.csv(out, OUT_CSV, row.names = FALSE)
cat(sprintf("Written: %s (%d rows)\n", OUT_CSV, nrow(out)))

cat("\n=== Real data set verdicts ===\n")
print(real_df[, c("dataset","d","dataset_in_correct_row",
                   "verdict_rk","verdict_un","verdict_sun","un_sun_relationship")],
      row.names = FALSE)

cat("\n=== Simulation grid verdicts ===\n")
print(sim_df[, c("d","verdict_rk","verdict_un","verdict_sun","un_sun_relationship")],
      row.names = FALSE)

cat("\n=== Tallies (real data, 16 x 3 = 48 claims) ===\n")
vs <- c(real_df$verdict_rk, real_df$verdict_un, real_df$verdict_sun)
bucket <- ifelse(grepl("^FALSE", vs), "FALSE",
           ifelse(grepl("^UNVERIFIABLE", vs), "UNVERIFIABLE",
           ifelse(grepl("^VERIFIED", vs), "VERIFIED", "NO CLAIM")))
print(table(bucket))

cat("\n=== Tallies (simulation grid, 7 x 3 = 21 claims) ===\n")
vs2 <- c(sim_df$verdict_rk, sim_df$verdict_un, sim_df$verdict_sun)
bucket2 <- ifelse(grepl("^FALSE", vs2), "FALSE",
            ifelse(grepl("^UNVERIFIABLE", vs2), "UNVERIFIABLE",
            ifelse(grepl("^VERIFIED", vs2), "VERIFIED", "NO CLAIM")))
print(table(bucket2))

cat("\ndone\n")
