#!/usr/bin/env Rscript
# 46_verify_wp2a_alpha.R -- independent verification of 42_wp2a_alpha_sweep.R.
#
# 42 returned three claims that would change the manuscript, so each is
# re-derived here from the primary artefacts rather than from 42's summary CSV.
#
# (A) NN DUPLICATE TABLES. 42 reports that the NN quantile files at many d are
#     byte-identical pairs, one MC run saved under two names. File hashing has
#     already confirmed 21 such pairs. What hashing CANNOT say is *which* alpha
#     the shared content represents -- and that matters, because the published
#     real-data runs consumed one of these files believing it to be a specific
#     level. d = 5, 10 and 20 carry genuinely distinct pairs, so they calibrate
#     the answer: the lower-tail quantile must be monotone in alpha (the 0.1%
#     quantile lies below the 1% quantile pointwise). Compare a duplicated d
#     against its distinct neighbours and report which level its content tracks.
#
# (B) RK RAW DRAWS. 42 claims RK's `simul` retains the full per-iteration draw
#     matrix Kest.m, so every alpha is recoverable from disk at zero cost, and
#     that $quan is reproducible as apply(Kest.m, 2, quantile, probs = q).
#     Re-derive that equality here on a file 42 did not special-case.
#
# (C) EFFECT SIZE. Recompute the per-group BA/F2 ranges and Jaccard floors from
#     the raw per-cell CSV, ignoring wp2a_alpha_summary.csv entirely. 42 reports
#     a largest dBA of 0.482 and "no plateau"; if its summary builder is wrong,
#     (C) and the summary disagree. Also splits informative (RK, and the NN d
#     with distinct files) from uninformative (duplicate-table NN) groups --
#     a group whose two inputs are the same bytes cannot evidence anything, and
#     must not be counted as "stable".
#
# Read-only with respect to every result file and every table. Writes nothing.

suppressMessages(library(here))

NNDIR <- here::here("R/NN-test_quantile")
RKDIR <- here::here("R/RK-test_quantile")
SWEEP <- here::here("revision_experiments/results/tr1/wp2a_alpha_sweep.csv")
SUMM  <- here::here("revision_experiments/results/tr1/wp2a_alpha_summary.csv")

load1 <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("simul", envir = e)) stop("no object named 'simul' in ", basename(path))
  get("simul", envir = e)
}
nnf <- function(d, tok) file.path(NNDIR, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))

# ---- (A) which alpha do the duplicated NN files actually hold? -------------
cat("=== (A) NN duplicate tables: identifying the surviving alpha ===\n")

# Calibrate on the dimensions that carry genuinely distinct files.
cat("\n-- reference dimensions with DISTINCT 99 / 999 files --\n")
for (d in c(10, 20)) {
  f99 <- nnf(d, "99"); f999 <- nnf(d, "999")
  if (!file.exists(f99) || !file.exists(f999)) { cat(sprintf("  d=%d: missing\n", d)); next }
  a <- load1(f99); b <- load1(f999)
  idx <- c(50, 100, 200, 500, 1000)
  idx <- idx[idx <= length(a$average)]
  cat(sprintf("  d=%2d  identical=%s\n", d, identical(a, b)))
  cat(sprintf("        99%%  average[%s] = %s\n", paste(idx, collapse = ","),
              paste(sprintf("%.5f", a$average[idx]), collapse = " ")))
  cat(sprintf("        999%% average[%s] = %s\n", paste(idx, collapse = ","),
              paste(sprintf("%.5f", b$average[idx]), collapse = " ")))
  cat(sprintf("        999 < 99 pointwise (expected for a lower-tail quantile): %s\n",
              all(b$average <= a$average, na.rm = TRUE)))
}

# Now the duplicated dimensions, compared against the distinct d=20 pair.
cat("\n-- duplicated dimensions, content compared against the d=20 reference --\n")
ref99  <- load1(nnf(20, "99"))
ref999 <- load1(nnf(20, "999"))
for (d in c(19, 21, 22)) {
  f99 <- nnf(d, "99"); f999 <- nnf(d, "999")
  if (!file.exists(f99) || !file.exists(f999)) { cat(sprintf("  d=%d: missing a file\n", d)); next }
  a <- load1(f99); b <- load1(f999)
  same <- identical(a, b)
  idx <- c(50, 100, 200, 500, 1000); idx <- idx[idx <= length(a$average)]
  # distance of this file's content to each d=20 reference level
  m  <- min(length(a$average), length(ref99$average))
  d99  <- mean(abs(a$average[1:m] - ref99$average[1:m]),  na.rm = TRUE)
  d999 <- mean(abs(a$average[1:m] - ref999$average[1:m]), na.rm = TRUE)
  cat(sprintf("  d=%2d  files identical=%-5s  mean|content - d20@99%%|=%.6f   mean|content - d20@999%%|=%.6f  -> nearer %s\n",
              d, same, d99, d999, if (d99 < d999) "99%" else "999%"))
  cat(sprintf("        content average[%s] = %s\n", paste(idx, collapse = ","),
              paste(sprintf("%.5f", a$average[idx]), collapse = " ")))
}

# Is that "nearer" verdict actually decisive? Only if the gap BETWEEN the two
# alpha levels is large relative to the gap BETWEEN adjacent dimensions. If the
# dimension effect dominates, the discriminator is measuring d, not alpha, and
# the question of which level survived cannot be settled this way.
m <- min(length(ref99$average), length(ref999$average))
gap_level <- mean(abs(ref99$average[1:m] - ref999$average[1:m]), na.rm = TRUE)
a19 <- load1(nnf(19, "99")); a21 <- load1(nnf(21, "99"))
m2  <- min(length(a19$average), length(a21$average))
gap_dim <- mean(abs(a19$average[1:m2] - a21$average[1:m2]), na.rm = TRUE) / 2
cat(sprintf("\n  DISCRIMINATOR POWER:\n"))
cat(sprintf("    mean |99%% - 999%%| at fixed d=20 (the alpha effect)      = %.6f\n", gap_level))
cat(sprintf("    mean per-unit-d change around d=20 (the dimension effect) = %.6f\n", gap_dim))
cat(sprintf("    ratio dimension/alpha = %.2f  -> verdict is %s\n", gap_dim / gap_level,
            if (gap_dim > gap_level) "NOT DECISIVE (dimension effect swamps alpha effect)" else "usable"))

# ---- (B) RK stores raw draws? ----------------------------------------------
cat("\n=== (B) RK simul structure and quantile reproducibility ===\n")
rk_probe <- c("RK-test-simul_9d_99%.RData", "RK-test-simul_20d_999%.RData")
for (fn in rk_probe) {
  p <- file.path(RKDIR, fn)
  if (!file.exists(p)) { cat(sprintf("  %s: missing\n", fn)); next }
  s <- load1(p)
  # NB: use [[ ]] throughout. R's $ does PARTIAL matching on lists, so s$q
  # silently resolves to s$quan and dumps a 50,000-element vector.
  km <- s[["Kest.m"]]; qn <- s[["quan"]]
  cat(sprintf("\n  %s\n    names: %s\n", fn, paste(names(s), collapse = ", ")))
  if (!is.null(km)) {
    cat(sprintf("    dim(Kest.m) = %s   class=%s   storage=%s\n",
                paste(dim(km), collapse = " x "), class(km)[1], typeof(km)))
  }
  if (!is.null(qn)) {
    qv <- if (is.list(qn)) qn[[1]] else qn
    cat(sprintf("    quan: is.list=%s  outer length=%d  inner length=%d\n",
                is.list(qn), length(qn), length(qv)))
  }
  if (!is.null(km) && !is.null(qn)) {
    # infer the level from the filename token, then reproduce $quan from the draws
    tok <- sub("^.*_([0-9]+)%\\.RData$", "\\1", fn)
    q   <- as.numeric(paste0("0.", tok))
    qv  <- if (is.list(qn)) qn[[1]] else qn
    rec <- apply(km, 2, stats::quantile, probs = q, names = FALSE, na.rm = TRUE)
    tgt <- as.numeric(qv)
    if (length(rec) != length(tgt)) {
      cat(sprintf("    length mismatch: recomputed %d vs stored %d -- cannot compare\n",
                  length(rec), length(tgt)))
    } else {
      md <- max(abs(rec - tgt), na.rm = TRUE)
      cat(sprintf("    recomputed quantile(Kest.m, %.3f) vs stored quan: max|diff|=%.3e  equal=%s\n",
                  q, md, md < 1e-9))
    }
    cat(sprintf("    -> RAW PER-ITERATION DRAWS RETAINED: %s (%d rows)\n",
                !is.null(dim(km)) && nrow(km) > 1, if (is.null(dim(km))) 0L else nrow(km)))
  }
}

# ---- (C) recompute effect sizes from the raw per-cell CSV ------------------
cat("\n=== (C) effect size recomputed from wp2a_alpha_sweep.csv ===\n")
df <- read.csv(SWEEP, stringsAsFactors = FALSE, colClasses = c(flagged_idx = "character"))
df$flagged_idx[is.na(df$flagged_idx)] <- ""
ok <- df[df$status == "ok", ]
cat(sprintf("  rows=%d  ok=%d  timeout=%d  not_run=%d\n", nrow(df), nrow(ok),
            sum(df$status == "timeout"), sum(df$status == "not_run")))

pidx <- function(s) if (!nzchar(s)) integer(0) else sort(as.integer(strsplit(s, ";")[[1]]))
jac  <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

# a group is informative only if its alpha levels came from genuinely different inputs
dupd <- c(6,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28)
grp <- unique(ok[, c("dataset", "method", "variant", "d")])
res <- list()
for (i in seq_len(nrow(grp))) {
  g <- grp[i, ]
  sub <- ok[ok$dataset == g$dataset & ok$method == g$method, ]
  if (nrow(sub) < 2) next
  sets <- lapply(sub$flagged_idx, pidx)
  allsame <- all(vapply(sets[-1], function(s) identical(s, sets[[1]]), logical(1)))
  jmin <- 1
  for (a in seq_along(sets)) for (b in seq_along(sets)) if (a < b) jmin <- min(jmin, jac(sets[[a]], sets[[b]]))
  informative <- (g$variant == "RK") || !(g$d %in% dupd)
  res[[length(res) + 1]] <- data.frame(
    dataset = g$dataset, method = g$method, variant = g$variant, d = g$d,
    n_alpha = nrow(sub), informative = informative, labels_identical = allsame,
    BA_range = max(sub$BA) - min(sub$BA), F2_range = max(sub$F2) - min(sub$F2),
    jaccard_min = jmin, stringsAsFactors = FALSE)
}
res <- do.call(rbind, res)

inf <- res[res$informative, ]; uninf <- res[!res$informative, ]
cat(sprintf("\n  groups total=%d   informative=%d   uninformative(duplicate-table NN)=%d\n",
            nrow(res), nrow(inf), nrow(uninf)))
cat(sprintf("  labels identical across alpha:  informative %d/%d   uninformative %d/%d\n",
            sum(inf$labels_identical), nrow(inf), sum(uninf$labels_identical), nrow(uninf)))
cat(sprintf("  -> every uninformative group identical? %s  (expected TRUE: same bytes in, same labels out)\n",
            all(uninf$labels_identical)))
cat(sprintf("\n  informative groups: median dBA=%.4f  mean dBA=%.4f  max dBA=%.4f\n",
            median(inf$BA_range), mean(inf$BA_range), max(inf$BA_range)))
cat(sprintf("                      median dF2=%.4f  median jaccard_min=%.4f  min jaccard_min=%.4f\n",
            median(inf$F2_range), median(inf$jaccard_min), min(inf$jaccard_min)))
k <- which.max(inf$BA_range)
cat(sprintf("  largest dBA: %.4f at %s / %s (d=%d, %s)\n", inf$BA_range[k],
            inf$dataset[k], inf$method[k], inf$d[k], inf$variant[k]))
j <- which.min(inf$jaccard_min)
cat(sprintf("  lowest jaccard: %.4f at %s / %s (d=%d)\n", inf$jaccard_min[j],
            inf$dataset[j], inf$method[j], inf$d[j]))

cat("\n  informative groups, sorted by BA_range:\n")
o <- inf[order(-inf$BA_range), c("dataset","method","variant","d","n_alpha","BA_range","F2_range","jaccard_min")]
print(o, row.names = FALSE, digits = 4)

# cross-check against 42's own summary
if (file.exists(SUMM)) {
  sm <- read.csv(SUMM, stringsAsFactors = FALSE)
  cat(sprintf("\n  42's summary rows=%d; columns: %s\n", nrow(sm), paste(names(sm), collapse = ", ")))
  if ("BA_range" %in% names(sm)) {
    cat(sprintf("  42 max BA_range=%.4f vs recomputed %.4f -> %s\n",
                max(sm$BA_range, na.rm = TRUE), max(res$BA_range),
                if (abs(max(sm$BA_range, na.rm = TRUE) - max(res$BA_range)) < 1e-6) "AGREES" else "DISAGREES"))
  }
}

cat("\ndone\n")
