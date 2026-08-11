#!/usr/bin/env Rscript
# 62_audit_nn_duplicates.R -- exhaustive audit of the NN quantile table
# duplication, per dimension.
#
# The claim "21 pairs are byte-identical" came from one MD5 sweep, and it now
# underpins a disclosure in the manuscript about which alpha levels can be
# stated. It deserves a stricter derivation than file hashing:
#
#   - Hashing proves FILE identity. What actually matters is CONTENT identity:
#     two files could differ in a serialization detail (timestamp, compression)
#     while holding the same numbers, or -- far less likely -- collide. Compare
#     the loaded $average and $median vectors directly.
#   - A per-dimension view catches what a global "list the duplicate groups"
#     sweep hides: dimensions carrying only one file (nothing to compare),
#     dimensions with three or more files, and dimensions whose files are
#     genuinely distinct.
#   - Report every dimension, not only the duplicated ones, so the denominator
#     is visible and the count can be checked rather than trusted.
#
# Read-only.

suppressMessages(library(here))
DIR <- here::here("R/NN-test_quantile")

fs <- list.files(DIR, pattern = "^NN-test-simul_.*\\.RData$", full.names = TRUE)
parse1 <- function(p) {
  b <- basename(p)
  d <- suppressWarnings(as.integer(sub("^NN-test-simul_(\\d+)d_.*$", "\\1", b)))
  tk <- sub("^NN-test-simul_\\d+d_([0-9]+)%.*$", "\\1", b)
  spliced <- grepl("spliced", b)
  data.frame(path = p, base = b, d = d, tok = tk, spliced = spliced,
             size = file.size(p), stringsAsFactors = FALSE)
}
meta <- do.call(rbind, lapply(fs, parse1))
meta <- meta[!is.na(meta$d) & grepl("^[0-9]+$", meta$tok), ]
meta <- meta[order(meta$d, as.numeric(paste0("0.", meta$tok))), ]

cat(sprintf("NN table files parsed: %d across %d dimensions\n",
            nrow(meta), length(unique(meta$d))))
cat(sprintf("(%d spliced variants included and flagged separately)\n\n",
            sum(meta$spliced)))

load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

dup_pairs <- 0L; distinct_pairs <- 0L; single <- 0L
rows <- list()
for (dd in sort(unique(meta$d))) {
  sub <- meta[meta$d == dd, ]
  if (nrow(sub) == 1) {
    single <- single + 1L
    rows[[length(rows) + 1]] <- data.frame(
      d = dd, n_files = 1L, tokens = sub$tok[1], verdict = "single file",
      max_abs_diff = NA_real_, stringsAsFactors = FALSE)
    next
  }
  tabs <- lapply(sub$path, load1)
  for (i in 1:(nrow(sub) - 1)) for (j in (i + 1):nrow(sub)) {
    a <- tabs[[i]]; b <- tabs[[j]]
    same <- identical(a$average, b$average) && identical(a$median, b$median)
    m <- min(length(a$average), length(b$average))
    md <- max(abs(a$average[1:m] - b$average[1:m]), na.rm = TRUE)
    if (same) dup_pairs <- dup_pairs + 1L else distinct_pairs <- distinct_pairs + 1L
    rows[[length(rows) + 1]] <- data.frame(
      d = dd, n_files = nrow(sub),
      tokens = sprintf("%s vs %s", sub$tok[i], sub$tok[j]),
      verdict = if (same) "IDENTICAL CONTENT" else "distinct",
      max_abs_diff = md, stringsAsFactors = FALSE)
  }
  rm(tabs); gc(verbose = FALSE)
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE, digits = 5)

cat(sprintf("\n=== TOTALS ===\n"))
cat(sprintf("  dimensions with a single file (nothing to compare): %d\n", single))
cat(sprintf("  pairs compared                                    : %d\n",
            dup_pairs + distinct_pairs))
cat(sprintf("  pairs IDENTICAL in content                        : %d\n", dup_pairs))
cat(sprintf("  pairs distinct                                    : %d\n", distinct_pairs))
d_dup <- sort(unique(res$d[res$verdict == "IDENTICAL CONTENT"]))
cat(sprintf("  dimensions with at least one identical pair       : %d  (%s)\n",
            length(d_dup), paste(d_dup, collapse = ",")))
d_dis <- sort(unique(res$d[res$verdict == "distinct"]))
cat(sprintf("  dimensions with a genuinely distinct pair         : %d  (%s)\n",
            length(d_dis), paste(d_dis, collapse = ",")))

# does file-level hashing agree with content comparison?
cat("\n=== hash vs content cross-check ===\n")
h <- tools::md5sum(meta$path)
hd <- 0L
for (dd in sort(unique(meta$d))) {
  sub <- meta[meta$d == dd, ]
  if (nrow(sub) < 2) next
  for (i in 1:(nrow(sub) - 1)) for (j in (i + 1):nrow(sub))
    if (h[[sub$path[i]]] == h[[sub$path[j]]]) hd <- hd + 1L
}
cat(sprintf("  pairs with identical MD5      : %d\n", hd))
cat(sprintf("  pairs with identical content  : %d\n", dup_pairs))
cat(sprintf("  -> %s\n", if (hd == dup_pairs)
  "AGREE (byte identity and content identity coincide)" else
  "DISAGREE -- some pair matches on one measure but not the other"))
cat("\ndone\n")
