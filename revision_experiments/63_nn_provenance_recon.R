#!/usr/bin/env Rscript
# 63_nn_provenance_recon.R -- reconnaissance for the alpha-provenance study.
#
# Thirteen of the sixteen real data sets sit at a dimension where the two
# shipped NN quantile files are identical in content, so the alpha level the
# published result was computed at is not readable off the filename. Before
# any Monte Carlo is commissioned we need three facts per affected dimension:
#
#   1. WHICH TOKEN each method actually loads (nn_quant_label_paper_UN/SUN),
#      and therefore which claimed alpha is at risk for which data set.
#   2. THE SAMPLE SIZE the shipped table was generated at -- length(average)
#      is (n - 1), since the generator tabulates subsample sizes 2..n. This
#      bounds the n a fresh replicate has to reach to be comparable, and the
#      NN generator costs O(n^3 d) per iteration, so it is the whole budget.
#   3. WHETHER THE TWO SHIPPED FILES DIFFER at that d -- re-derived here per
#      dimension rather than taken from 62's global count, so the per-data-set
#      table stands on its own.
#
# Read-only. Emits results/tr1/wp2a_provenance_targets.csv.

suppressMessages(library(here))
source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

inv <- read.csv(here::here("revision_experiments/results/tr1/dataset_inventory.csv"),
                stringsAsFactors = FALSE)
NNDIR <- here::here("R/NN-test_quantile")

tabs_at <- function(d) {
  fs <- list.files(NNDIR, pattern = sprintf("^NN-test-simul_%dd_[0-9]+%%\\.RData$", d),
                   full.names = TRUE)
  toks <- sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs))
  o <- order(as.numeric(paste0("0.", toks))); fs[o]
}
load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

rows <- list()
for (i in seq_len(nrow(inv))) {
  d <- inv$d[i]; n <- inv$n[i]
  fs <- tabs_at(d)
  toks <- sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs))
  un  <- nn_quant_label_paper_UN(d)
  sun <- nn_quant_label_paper_SUN(d)
  dup <- NA; len_un <- NA_integer_; len_sun <- NA_integer_
  if (length(fs) >= 2) {
    a <- load1(fs[1]); b <- load1(fs[2])
    dup <- identical(a$average, b$average) && identical(a$median, b$median)
  }
  fu <- file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, un))
  fsn <- file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, sun))
  if (file.exists(fu))  len_un  <- length(load1(fu)$average)
  if (file.exists(fsn)) len_sun <- length(load1(fsn)$average)
  rows[[i]] <- data.frame(
    dataset = inv$dataset[i], n = n, d = d,
    tokens_on_disk = paste(toks, collapse = "/"),
    UN_token = un, SUN_token = sun,
    same_table_for_UN_and_SUN = (un == sun) || isTRUE(dup),
    pair_identical = dup,
    shipped_n_UN  = if (is.na(len_un))  NA_integer_ else len_un + 1L,
    shipped_n_SUN = if (is.na(len_sun)) NA_integer_ else len_sun + 1L,
    covers_dataset_n = if (is.na(len_un)) NA else (len_un + 1L) >= n,
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
res <- res[order(res$n), ]
print(res, row.names = FALSE)

cat("\n=== affected data sets (pair identical at their d) ===\n")
aff <- res[which(res$pair_identical), ]
cat(sprintf("  %d of %d data sets\n", nrow(aff), nrow(res)))
cat(sprintf("  dimensions involved: %s\n", paste(sort(unique(aff$d)), collapse = ",")))

cat("\n=== data sets where the paper claims DIFFERENT alpha for UN vs SUN\n")
cat("    but both methods load the SAME file ===\n")
nom <- res[which(res$pair_identical & res$UN_token != res$SUN_token), ]
if (nrow(nom)) print(nom[, c("dataset","n","d","UN_token","SUN_token")], row.names = FALSE) else
  cat("  none\n")

# generation budget: one table per dimension at the largest n needed there
cat("\n=== per-dimension generation budget (one table serves every data set at that d) ===\n")
# 145 s/iteration measured at n=5000, d=21 (58_gen_d21_chunked.R); cost is O(n^3 d)
cost <- function(n, d) 145 * (n / 5000)^3 * (d / 21)
for (d in sort(unique(aff$d))) {
  nn <- max(aff$n[aff$d == d]); ds <- aff$dataset[aff$d == d]
  cat(sprintf("  d=%2d  n_max=%4d  (%s)  ~%.3f s/iter  ->  250 iters = %s\n",
              d, nn, paste(ds, collapse = ", "), cost(nn, d),
              format(round(cost(nn, d) * 250 / 60, 1), nsmall = 1)))
}

write.csv(res, here::here("revision_experiments/results/tr1/wp2a_provenance_targets.csv"),
          row.names = FALSE)
cat("\nwrote results/tr1/wp2a_provenance_targets.csv\ndone\n")
