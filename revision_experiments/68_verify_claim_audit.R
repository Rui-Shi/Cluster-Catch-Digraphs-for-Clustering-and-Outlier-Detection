#!/usr/bin/env Rscript
# 68_verify_claim_audit.R -- independent check of 66's documentary audit.
#
# 66 concluded that the damage is confined to the real-data NN column: the RK
# column of tab:alpha_real and the whole simulation schedule were reported
# sound. Those are exculpatory findings, and exculpatory findings deserve the
# same scrutiny as incriminating ones -- if wrong, they leave a false claim
# standing in the manuscript. Two specific things are re-derived here:
#
#   1. RK DISTINCTNESS at the dimensions the sixteen real data sets actually
#      use. 66 asserts every RK pair is distinct. Checked by CONTENT (loading
#      and comparing the objects), not by file size, since two Monte Carlo
#      runs can easily produce files of equal size.
#   2. THE SIMULATION GRID d in {2,3,5,10,20,50,100}: for each method, does a
#      file exist at the claimed token, and is it distinct from its siblings
#      at that d?
#
# 66 also reported that sourcing harness.R segfaulted its R session, so it
# transcribed the four resolver functions rather than calling them. This
# script calls the LIVE functions, so any divergence between the transcription
# and the real code shows up as a disagreement with 66's CSV -- which is
# checked directly at the end.
#
# Read-only apart from its own CSV.

suppressMessages(library(here))
source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

RKDIR <- here::here("R/RK-test_quantile")
NNDIR <- here::here("R/NN-test_quantile")
load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

#' Content comparison that works for both table families: RK tables carry the
#' raw draw matrix Kest.m as well as the reduced quantile, NN tables only the
#' reduction. Compare every numeric field the two share.
same_content <- function(a, b) {
  ka <- intersect(names(a), names(b))
  if (!length(ka)) return(NA)
  all(vapply(ka, function(k) identical(a[[k]], b[[k]]), logical(1)))
}
sibs <- function(dir, prefix, d) {
  fs <- list.files(dir, pattern = sprintf("^%s_%dd_[0-9]+%%\\.RData$", prefix, d),
                   full.names = TRUE)
  fs[order(as.numeric(paste0("0.", sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs)))))]
}

inv <- read.csv(here::here("revision_experiments/results/tr1/dataset_inventory.csv"),
                stringsAsFactors = FALSE)

cat("=== 1. RK tables at the dimensions the real data sets use ===\n")
rk_rows <- list()
for (d in sort(unique(inv$d))) {
  fs <- sibs(RKDIR, "RK-test-simul", d)
  tok_claimed <- rk_quant_for_d(d)
  f_used <- file.path(RKDIR, sprintf("RK-test-simul_%dd_%s%%.RData", d, tok_claimed))
  dup <- NA
  if (length(fs) >= 2) {
    a <- load1(fs[1]); b <- load1(fs[2]); dup <- same_content(a, b); rm(a, b); gc(verbose = FALSE)
  }
  rk_rows[[length(rk_rows) + 1]] <- data.frame(
    d = d, n_files = length(fs),
    tokens = paste(sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs)), collapse = "/"),
    token_used = tok_claimed, file_exists = file.exists(f_used),
    pair_identical = dup,
    datasets = paste(inv$dataset[inv$d == d], collapse = ";"), stringsAsFactors = FALSE)
  cat(sprintf("  d=%3d  files=%d [%s]  used=%s  exists=%s  pair_identical=%s\n",
              d, length(fs), rk_rows[[length(rk_rows)]]$tokens, tok_claimed,
              file.exists(f_used), dup))
}
rk <- do.call(rbind, rk_rows)
cat(sprintf("  -> RK pairs identical anywhere: %d ; RK files missing: %d\n",
            sum(rk$pair_identical %in% TRUE), sum(!rk$file_exists)))

cat("\n=== 2. simulation grid d in {2,3,5,10,20,50,100} ===\n")
sim_rows <- list()
for (d in c(2, 3, 5, 10, 20, 50, 100)) {
  fs <- sibs(NNDIR, "NN-test-simul", d)
  toks <- sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs))
  un <- nn_quant_label_paper_UN(d); sn <- nn_quant_label_paper_SUN(d)
  dup <- NA
  if (length(fs) >= 2) {
    # is the file UN uses identical to ANY sibling at this d?
    fu <- file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, un))
    if (file.exists(fu)) {
      a <- load1(fu)
      dup <- any(vapply(setdiff(fs, fu), function(g) isTRUE(same_content(a, load1(g))), logical(1)))
    }
  }
  ok_un <- file.exists(file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, un)))
  ok_sn <- file.exists(file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, sn)))
  sim_rows[[length(sim_rows) + 1]] <- data.frame(
    d = d, tokens = paste(toks, collapse = "/"), UN_token = un, SUN_token = sn,
    UN_file_exists = ok_un, SUN_file_exists = ok_sn,
    UN_file_duplicated = dup, stringsAsFactors = FALSE)
  cat(sprintf("  d=%3d  on disk [%s]  UN->%s (%s)  SUN->%s (%s)  UN dup=%s\n",
              d, paste(toks, collapse = "/"), un, ok_un, sn, ok_sn, dup))
}
sim <- do.call(rbind, sim_rows)

cat("\n=== 3. live resolvers vs 66's transcription ===\n")
f66 <- here::here("revision_experiments/results/tr1/wp2a_alpha_claim_audit.csv")
if (file.exists(f66)) {
  a66 <- read.csv(f66, stringsAsFactors = FALSE)
  cat(sprintf("  66's CSV: %d rows, columns: %s\n", nrow(a66),
              paste(names(a66), collapse = ", ")))
  # find whichever columns hold the resolved tokens and compare to live calls
  cand <- grep("token|quant|alpha", names(a66), ignore.case = TRUE, value = TRUE)
  cat(sprintf("  candidate token columns: %s\n", paste(cand, collapse = ", ")))
  if ("dataset" %in% names(a66)) {
    m <- merge(a66, inv[, c("dataset", "d")], by = "dataset", suffixes = c("", "_inv"))
    dcol <- if ("d" %in% names(m)) "d" else "d_inv"
    m$live_UN  <- vapply(m[[dcol]], nn_quant_label_paper_UN, character(1))
    m$live_SUN <- vapply(m[[dcol]], nn_quant_label_paper_SUN, character(1))
    for (cc in cand) {
      if (all(as.character(m[[cc]]) == m$live_UN))  cat(sprintf("  column %-24s == live UN resolver\n", cc))
      if (all(as.character(m[[cc]]) == m$live_SUN)) cat(sprintf("  column %-24s == live SUN resolver\n", cc))
    }
  }
} else cat("  66's CSV not found -- cannot cross-check\n")

# also confirm the live resolvers agree with what 63 recorded (which used them)
t63 <- read.csv(here::here("revision_experiments/results/tr1/wp2a_provenance_targets.csv"),
                stringsAsFactors = FALSE)
t63$live_UN  <- vapply(t63$d, nn_quant_label_paper_UN, character(1))
t63$live_SUN <- vapply(t63$d, nn_quant_label_paper_SUN, character(1))
cat(sprintf("  63's recorded tokens still reproduce from live code: UN %s, SUN %s\n",
            all(as.character(t63$UN_token) == t63$live_UN),
            all(as.character(t63$SUN_token) == t63$live_SUN)))

write.csv(rk,  here::here("revision_experiments/results/tr1/wp2a_verify_rk_distinct.csv"), row.names = FALSE)
write.csv(sim, here::here("revision_experiments/results/tr1/wp2a_verify_sim_grid.csv"),   row.names = FALSE)
cat("\nwrote wp2a_verify_rk_distinct.csv and wp2a_verify_sim_grid.csv\ndone\n")
