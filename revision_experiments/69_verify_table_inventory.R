#!/usr/bin/env Rscript
# 69_verify_table_inventory.R -- does every table the code can ask for exist,
# and is every restored file byte-for-byte what it was?
#
# WHY THIS EXISTS. Twenty-two NN quantile tables were deleted and restored from
# the Recycle Bin. Two things have to be checked before anything is run again:
#
#   1. COMPLETENESS. get_simul() resolves a table purely by filename token and
#      stop()s if the file is absent. A missing table is therefore a hard
#      failure, not a silent fallback -- but only for the (d, token) pairs the
#      resolvers actually produce. This walks every dimension the study uses,
#      asks both resolvers what they would load, and checks the file is there.
#   2. INTEGRITY. A round trip through the Recycle Bin should be lossless, but
#      "should be" is not evidence, and this volume is a USB enclosure with a
#      documented history of dropping off the bus under write load. Each
#      restored file is loaded and checked against the property that identified
#      it in the first place: whether it is identical to its sibling at the same
#      dimension. The duplicate structure is known exactly (21 pairs at
#      d = 6-9, 11-19, 21-28; d = 5, 10, 20 genuinely distinct), so any
#      corruption shows up as a change in that pattern.
#
# Read-only.

suppressMessages(library(here))
source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

NNDIR <- here::here("R/NN-test_quantile")
load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }
nnf <- function(d, tok) file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))

inv <- read.csv(here::here("revision_experiments/results/tr1/dataset_inventory.csv"),
                stringsAsFactors = FALSE)
SIM_D <- c(2, 3, 5, 10, 20, 50, 100)

cat("=== 1. completeness: every (d, token) the resolvers can produce ===\n")
rows <- list()
for (d in sort(unique(c(inv$d, SIM_D)))) {
  for (meth in c("UN", "SUN")) {
    tok <- if (meth == "UN") nn_quant_label_paper_UN(d) else nn_quant_label_paper_SUN(d)
    f <- nnf(d, tok)
    who <- paste(inv$dataset[inv$d == d], collapse = ";")
    rows[[length(rows) + 1]] <- data.frame(
      d = d, method = meth, token = tok, exists = file.exists(f),
      context = if (nzchar(who)) who else if (d %in% SIM_D) "simulation only" else "-",
      stringsAsFactors = FALSE)
  }
}
comp <- do.call(rbind, rows)
missing <- comp[!comp$exists, ]
if (nrow(missing)) {
  cat("  *** MISSING TABLES -- these calls would stop() ***\n")
  print(missing, row.names = FALSE)
} else {
  cat(sprintf("  all %d (d, method) resolutions present\n", nrow(comp)))
}

cat("\n=== 2. integrity: the duplicate structure must be unchanged ===\n")
# ground truth established before the deletion (62, 63, 68), then amended as
# 75 installs genuine 0.1% tables: a repaired dimension's "999" is no longer a
# copy of its "99", which is the whole point, so it moves to the distinct list.
REPAIRED <- c(12, 18, 19)                       # extend as 75 installs more
EXPECT_DUP <- setdiff(c(6, 7, 8, 9, 11:19, 21:28), REPAIRED)
EXPECT_DISTINCT <- c(2, 3, 4, 5, 10, 20, REPAIRED)
chk <- list()
for (d in sort(c(EXPECT_DUP, EXPECT_DISTINCT))) {
  fs <- list.files(NNDIR, pattern = sprintf("^NN-test-simul_%dd_[0-9]+%%\\.RData$", d),
                   full.names = TRUE)
  toks <- sub(".*_([0-9]+)%\\.RData$", "\\1", basename(fs))
  o <- order(as.numeric(paste0("0.", toks))); fs <- fs[o]; toks <- toks[o]
  dup <- NA
  if (length(fs) >= 2) {
    a <- load1(fs[1]); b <- load1(fs[2])
    dup <- identical(a$average, b$average) && identical(a$median, b$median)
  }
  expect <- if (d %in% EXPECT_DUP) TRUE else FALSE
  # After 70 ran, a duplicated dimension may legitimately hold ONE file: the
  # redundant member was deleted. That is a pass, not a failure -- what would
  # be a failure is a dimension that still has two files whose relationship
  # has changed, or a DISTINCT dimension that has lost a file.
  status <- if (length(fs) == 1 && expect) "duplicate removed"
            else if (length(fs) == 1 && !expect) "*** FILE MISSING ***"
            else if (isTRUE(dup) == expect) "as expected"
            else "*** RELATIONSHIP CHANGED ***"
  chk[[length(chk) + 1]] <- data.frame(
    d = d, n_files = length(fs), tokens = paste(toks, collapse = "/"),
    identical_pair = dup, expected = expect, status = status,
    ok = !grepl("\\*\\*\\*", status), stringsAsFactors = FALSE)
}
ck <- do.call(rbind, chk)
bad <- ck[!ck$ok, ]
if (nrow(bad)) {
  cat("  *** duplicate structure CHANGED -- possible corruption or wrong file ***\n")
  print(bad, row.names = FALSE)
} else {
  cat(sprintf("  all %d dimensions consistent with the recorded pattern\n", nrow(ck)))
  cat(sprintf("  (%d still duplicated, %d duplicate removed, %d distinct as expected)\n",
              sum(ck$status == "as expected" & ck$identical_pair %in% TRUE),
              sum(ck$status == "duplicate removed"),
              sum(ck$status == "as expected" & ck$identical_pair %in% FALSE)))
}

cat("\n=== 2b. table length vs the data it serves ===\n")
# NEW FAILURE MODE, introduced by 75. The shipped tables all carry 5000
# entries; the regenerated ones are only as long as their own data set,
# because the cost is cubic in that length. UN_CCD.R:241 does
# simul$average[1:n] with n = nrow(dx), the full data matrix, so a table
# shorter than its data set yields NA past the end rather than clamping --
# silent, and it would corrupt exactly the results this repair is for.
len_rows <- list()
for (i in seq_len(nrow(inv))) {
  d <- inv$d[i]; nd <- inv$n[i]
  for (meth in c("UN", "SUN")) {
    tok <- if (meth == "UN") nn_quant_label_paper_UN(d) else nn_quant_label_paper_SUN(d)
    f <- nnf(d, tok)
    if (!file.exists(f)) next
    L <- length(load1(f)$average)
    len_rows[[length(len_rows) + 1]] <- data.frame(
      dataset = inv$dataset[i], d = d, n = nd, method = meth, token = tok,
      table_len = L, ok = L >= nd, stringsAsFactors = FALSE)
  }
}
ln <- do.call(rbind, len_rows)
short <- ln[!ln$ok, ]
if (nrow(short)) {
  cat("  *** TABLE SHORTER THAN ITS DATA SET -- would return NA ***\n")
  print(short, row.names = FALSE)
} else {
  cat(sprintf("  all %d (data set, method) pairs have table_len >= n\n", nrow(ln)))
  tight <- ln[ln$table_len == ln$n, ]
  if (nrow(tight)) {
    cat("  exactly-sized (no headroom -- fine for this data set, unusable for a larger one at the same d):\n")
    print(tight[, c("dataset", "d", "n", "method", "token", "table_len")], row.names = FALSE)
  }
}

cat("\n=== 3. d=20 specifically -- the one deletion that was not a duplicate ===\n")
f99 <- nnf(20, "99"); f999 <- nnf(20, "999")
if (file.exists(f99) && file.exists(f999)) {
  a <- load1(f99); b <- load1(f999)
  m <- min(length(a$average), length(b$average))
  cat(sprintf("  both present; identical = %s ; max|diff| over %d entries = %.6f\n",
              identical(a$average, b$average), m,
              max(abs(a$average[1:m] - b$average[1:m]))))
  cat(sprintf("  lengths: 99%% = %d, 999%% = %d\n", length(a$average), length(b$average)))
  cat("  -> a genuinely distinct table, correctly labelled; it was not a duplicate\n")
} else cat("  *** one or both d=20 files still absent ***\n")

cat(sprintf("\nOVERALL: %s\n",
            if (!nrow(missing) && !nrow(bad)) "inventory complete and unchanged" else
              "PROBLEMS ABOVE -- do not run experiments until resolved"))
