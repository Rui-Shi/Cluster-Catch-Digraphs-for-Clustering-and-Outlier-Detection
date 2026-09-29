#!/usr/bin/env Rscript
# 47_alpha_operating_range.R -- how big is the alpha effect INSIDE the range the
# paper actually operates in?
#
# 42's sweep spans alpha in {0.10, 0.05, 0.01, 0.001}. The manuscript never uses
# 0.10 or 0.05 for the RK methods: its schedule is 1% for d < 10 and 0.1% for
# d >= 10. Reporting a sensitivity range that includes settings nobody would
# choose overstates the problem, exactly as reporting only the flat pairs would
# understate it. So split the two:
#
#   FULL      all swept levels                   -- the honest worst case
#   OPERATING alpha in {0.01, 0.001} only        -- what a user of this paper faces
#
# Read-only. Writes one summary CSV.

suppressMessages(library(here))
SWEEP <- here::here("revision_experiments/results/tr1/wp2a_alpha_sweep.csv")
OUT   <- here::here("revision_experiments/results/tr1/wp2a_alpha_operating_range.csv")

df <- read.csv(SWEEP, stringsAsFactors = FALSE, colClasses = c(flagged_idx = "character"))
df$flagged_idx[is.na(df$flagged_idx)] <- ""
ok <- df[df$status == "ok", ]

# duplicate-table dimensions: NN groups there compare a file with itself
dupd <- c(6,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28)
ok$informative <- (ok$variant == "RK") | !(ok$d %in% dupd)
ok <- ok[ok$informative, ]

pidx <- function(s) if (!nzchar(s)) integer(0) else sort(as.integer(strsplit(s, ";")[[1]]))
jac  <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

cat(sprintf("alpha levels present in the informative rows: %s\n",
            paste(sort(unique(ok$alpha)), collapse = ", ")))

summarise <- function(sub) {
  if (nrow(sub) < 2) return(NULL)
  sets <- lapply(sub$flagged_idx, pidx)
  jmin <- 1
  for (a in seq_along(sets)) for (b in seq_along(sets)) if (a < b) jmin <- min(jmin, jac(sets[[a]], sets[[b]]))
  data.frame(n_levels = nrow(sub),
             identical = all(vapply(sets[-1], function(s) identical(s, sets[[1]]), logical(1))),
             BA_range = max(sub$BA) - min(sub$BA),
             F2_range = max(sub$F2) - min(sub$F2),
             jaccard_min = jmin, stringsAsFactors = FALSE)
}

grp <- unique(ok[, c("dataset", "method", "d")])
rows <- list()
for (i in seq_len(nrow(grp))) {
  g   <- grp[i, ]
  sub <- ok[ok$dataset == g$dataset & ok$method == g$method, ]
  full <- summarise(sub)
  oper <- summarise(sub[sub$alpha %in% c(0.01, 0.001), ])
  if (is.null(full)) next
  rows[[length(rows) + 1]] <- data.frame(
    dataset = g$dataset, method = g$method, d = g$d,
    full_n = full$n_levels, full_BA_range = full$BA_range,
    full_F2_range = full$F2_range, full_jaccard_min = full$jaccard_min,
    oper_n  = if (is.null(oper)) NA_integer_ else oper$n_levels,
    oper_identical = if (is.null(oper)) NA else oper$identical,
    oper_BA_range = if (is.null(oper)) NA_real_ else oper$BA_range,
    oper_F2_range = if (is.null(oper)) NA_real_ else oper$F2_range,
    oper_jaccard_min = if (is.null(oper)) NA_real_ else oper$jaccard_min,
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
res <- res[order(-res$oper_BA_range), ]
write.csv(res, OUT, row.names = FALSE)

cat("\n=== FULL swept range (alpha 0.10 / 0.05 / 0.01 / 0.001) ===\n")
cat(sprintf("  groups=%d  median dBA=%.4f  mean dBA=%.4f  max dBA=%.4f  median jaccard_min=%.4f\n",
            nrow(res), median(res$full_BA_range), mean(res$full_BA_range),
            max(res$full_BA_range), median(res$full_jaccard_min)))

o <- res[!is.na(res$oper_BA_range), ]
cat("\n=== OPERATING range only (alpha 0.01 vs 0.001 -- the paper's schedule) ===\n")
cat(sprintf("  groups=%d  labels identical: %d/%d\n", nrow(o), sum(o$oper_identical), nrow(o)))
cat(sprintf("  median dBA=%.4f  mean dBA=%.4f  max dBA=%.4f\n",
            median(o$oper_BA_range), mean(o$oper_BA_range), max(o$oper_BA_range)))
cat(sprintf("  median dF2=%.4f  median jaccard_min=%.4f  min jaccard_min=%.4f\n",
            median(o$oper_F2_range), median(o$oper_jaccard_min), min(o$oper_jaccard_min)))
k <- which.max(o$oper_BA_range)
cat(sprintf("  largest operating dBA: %.4f at %s / %s (d=%d)\n",
            o$oper_BA_range[k], o$dataset[k], o$method[k], o$d[k]))
cat(sprintf("\n  shrinkage: median dBA falls %.4f -> %.4f (%.0f%% of the full-range effect)\n",
            median(res$full_BA_range), median(o$oper_BA_range),
            100 * median(o$oper_BA_range) / median(res$full_BA_range)))

cat("\n  per group:\n")
print(o[, c("dataset","method","d","oper_identical","oper_BA_range","oper_F2_range",
            "oper_jaccard_min","full_BA_range")], row.names = FALSE, digits = 4)
cat(sprintf("\nwrote %s\n", OUT))
