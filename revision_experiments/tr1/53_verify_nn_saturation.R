#!/usr/bin/env Rscript
# 53_verify_nn_saturation.R -- independent verification of 48_wp2a_nn_saturation.R.
#
# 48 reports that the manuscript's high-dimensional saturation claim is REFUTED
# and, further, that the alpha effect GROWS with dimension. That would remove a
# claim from the paper and replace it with its opposite, so it is re-derived
# here from the raw per-cell CSV and from the generated tables themselves,
# never from 48's summary.
#
# (A) TABLE DISTINCTNESS. The original claim died because two quantile files
#     were byte-identical. If 48's 28 new tables had any duplicate pair, the
#     replacement would inherit the same defect. 48 reports 28 distinct
#     SHA-256 hashes; hashing is necessary but not sufficient -- two tables
#     could differ in a timestamp attribute yet hold identical numbers. So
#     compare the CONTENT vectors pairwise within each d.
#
# (B) DESIGN COMPLETENESS. A "square" grid is what makes the d-comparison
#     valid. Verify every (setting, d, method, alpha) has the same replicate
#     count, that all statuses are ok, and that replicate seeds line up across
#     alpha within a d -- if alpha levels were run on DIFFERENT data, the
#     within-replicate range is meaningless.
#
# (C) EFFECT SIZES recomputed from the 4480 raw cells, ignoring the summary.
#
# (D) NOISE FLOOR. 48's own caveat is that BA degrades at high d, so a growing
#     BA range might reflect a fragile detector rather than a sharper test.
#     The design is paired (same replicate, same data, only alpha changes) and
#     the detectors are deterministic, so an alpha difference is genuinely
#     caused by alpha. But magnitude still needs context: compare the
#     within-replicate across-alpha range against the across-replicate spread
#     at FIXED alpha. If the alpha effect is small next to ordinary
#     replicate-to-replicate variation, "alpha matters more at high d" is
#     technically true and practically minor. Report both.
#
# Read-only. Writes nothing.

suppressMessages(library(here))

RAW  <- here::here("revision_experiments/results/tr1/wp2a_nn_saturation.csv")
SUMM <- here::here("revision_experiments/results/tr1/wp2a_nn_saturation_summary.csv")
TDIR <- here::here("R/NN-test_quantile_n200")

# ---- (A) are the 28 generated tables genuinely distinct in CONTENT? --------
cat("=== (A) distinctness of the generated n=200 tables ===\n")
if (!dir.exists(TDIR)) {
  cat(sprintf("  directory %s missing -- cannot check\n", TDIR))
} else {
  fs <- list.files(TDIR, pattern = "\\.RData$", full.names = TRUE)
  cat(sprintf("  %d .RData files found\n", length(fs)))
  meta <- do.call(rbind, lapply(fs, function(p) {
    b <- basename(p)
    dd <- suppressWarnings(as.integer(sub("^.*_(\\d+)d_.*$", "\\1", b)))
    tk <- sub("^.*_(\\d+)%\\.RData$", "\\1", b)
    data.frame(path = p, base = b, d = dd, tok = tk, stringsAsFactors = FALSE)
  }))
  meta <- meta[!is.na(meta$d) & grepl("^[0-9]+$", meta$tok), ]
  cat(sprintf("  parsed %d as (d, token) tables across d = %s\n", nrow(meta),
              paste(sort(unique(meta$d)), collapse = ",")))
  bad <- 0L; cmp <- 0L
  for (dd in sort(unique(meta$d))) {
    sub <- meta[meta$d == dd, ]
    if (nrow(sub) < 2) next
    tabs <- lapply(sub$path, function(p) { e <- new.env(); load(p, envir = e)
                                           get("simul", envir = e) })
    for (i in 1:(nrow(sub) - 1)) for (j in (i + 1):nrow(sub)) {
      cmp <- cmp + 1L
      sa <- identical(tabs[[i]]$average, tabs[[j]]$average)
      sm <- identical(tabs[[i]]$median,  tabs[[j]]$median)
      if (sa && sm) {
        bad <- bad + 1L
        cat(sprintf("  *** IDENTICAL CONTENT: d=%d %s vs %s\n", dd, sub$tok[i], sub$tok[j]))
      }
    }
    rm(tabs); gc(verbose = FALSE)
  }
  cat(sprintf("  %d within-d pairs compared, %d identical -> %s\n", cmp, bad,
              if (bad == 0) "ALL DISTINCT" else "DUPLICATES PRESENT -- RESULT INVALID"))
}

# ---- (B) design completeness ----------------------------------------------
cat("\n=== (B) design completeness ===\n")
df <- read.csv(RAW, stringsAsFactors = FALSE, colClasses = c(flagged_idx = "character"))
df$flagged_idx[is.na(df$flagged_idx)] <- ""
cat(sprintf("  rows=%d  statuses: %s\n", nrow(df),
            paste(sprintf("%s=%d", names(table(df$status)), table(df$status)), collapse = " ")))

acol <- if ("alpha" %in% names(df)) "alpha" else "alpha_token"
scol <- if ("setting" %in% names(df)) "setting" else "gen"
cat(sprintf("  grouping columns: setting='%s' alpha='%s'\n", scol, acol))
cat(sprintf("  d levels: %s\n", paste(sort(unique(df$d)), collapse = ",")))
cat(sprintf("  alpha levels: %s\n", paste(sort(unique(df[[acol]])), collapse = ",")))
cat(sprintf("  methods: %s\n", paste(sort(unique(df$method)), collapse = ",")))

tb <- table(df[[scol]], df$d, df$method, df[[acol]])
cat(sprintf("  cells per (setting,d,method,alpha): min=%d max=%d -> %s\n",
            min(tb), max(tb), if (min(tb) == max(tb)) "SQUARE" else "RAGGED"))

# do the alpha levels share replicate seeds within a (setting,d,method)?
if ("rep" %in% names(df)) {
  key <- unique(df[, c(scol, "d", "method")])
  mism <- 0L
  for (i in seq_len(nrow(key))) {
    s <- df[df[[scol]] == key[i, scol] & df$d == key[i, "d"] & df$method == key[i, "method"], ]
    reps <- split(sort(s$rep), s[[acol]])
    if (length(unique(lapply(reps, identity))) > 1) mism <- mism + 1L
  }
  cat(sprintf("  groups whose alpha levels use different replicate sets: %d -> %s\n",
              mism, if (mism == 0) "PAIRED CORRECTLY" else "NOT PAIRED -- range is meaningless"))
}

# ---- (C) effect sizes recomputed ------------------------------------------
cat("\n=== (C) alpha effect recomputed from raw cells ===\n")
pidx <- function(s) if (!nzchar(s)) integer(0) else sort(as.integer(strsplit(s, ";")[[1]]))
jacm <- function(sets) { m <- 1
  for (a in seq_along(sets)) for (b in seq_along(sets)) if (a < b) {
    u <- length(union(sets[[a]], sets[[b]]))
    m <- min(m, if (u == 0) 1 else length(intersect(sets[[a]], sets[[b]])) / u) }
  m }

key <- unique(df[, c(scol, "d", "method")])
out <- list()
for (i in seq_len(nrow(key))) {
  s <- df[df[[scol]] == key[i, scol] & df$d == key[i, "d"] & df$method == key[i, "method"], ]
  reps <- unique(s$rep)
  rng <- numeric(0); jm <- numeric(0); ident <- logical(0)
  orng <- numeric(0); oident <- logical(0)
  for (rp in reps) {
    r <- s[s$rep == rp, ]
    if (nrow(r) < 2) next
    rng   <- c(rng, max(r$BA) - min(r$BA))
    sets  <- lapply(r$flagged_idx, pidx)
    jm    <- c(jm, jacm(sets))
    ident <- c(ident, all(vapply(sets[-1], function(x) identical(x, sets[[1]]), logical(1))))
    o <- r[r[[acol]] %in% c(0.01, 0.001), ]
    if (nrow(o) == 2) {
      orng   <- c(orng, abs(diff(o$BA)))
      os     <- lapply(o$flagged_idx, pidx)
      oident <- c(oident, identical(os[[1]], os[[2]]))
    }
  }
  # (D) noise floor: spread of BA across replicates at a FIXED alpha
  fixed <- s[s[[acol]] == min(s[[acol]]), ]
  out[[length(out) + 1]] <- data.frame(
    setting = key[i, scol], d = key[i, "d"], method = key[i, "method"],
    n_reps = length(rng),
    mean_BA_range = mean(rng), mean_min_jaccard = mean(jm),
    frac_identical = mean(ident),
    oper_mean_BA_range = if (length(orng)) mean(orng) else NA_real_,
    oper_frac_identical = if (length(oident)) mean(oident) else NA_real_,
    BA_sd_across_reps_fixed_alpha = sd(fixed$BA),
    mean_BA = mean(s$BA),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
res$alpha_effect_vs_noise <- res$mean_BA_range / res$BA_sd_across_reps_fixed_alpha

cat("\n  recomputed (compare mean_BA_range against 48's summary):\n")
print(res[order(res$setting, res$method, res$d),
          c("setting","d","method","n_reps","mean_BA_range","mean_min_jaccard",
            "frac_identical","oper_mean_BA_range","oper_frac_identical")],
      row.names = FALSE, digits = 4)

# saturation criterion, re-applied
res$indist <- res$mean_BA_range < 0.01 & res$mean_min_jaccard > 0.95
cat(sprintf("\n  cells meeting the saturation criterion (range<0.01 AND jaccard>0.95): %d / %d\n",
            sum(res$indist), nrow(res)))
cat(sprintf("  -> %s\n", if (sum(res$indist) == 0)
            "NO SATURATION AT ANY DIMENSION -- claim refuted" else "some cells saturate"))

# trend direction
lowd  <- res[res$d <= 5,  ]
highd <- res[res$d >= 20, ]
cat(sprintf("\n  mean_BA_range: d<=5 mean=%.4f   d>=20 mean=%.4f   ratio=%.2f\n",
            mean(lowd$mean_BA_range), mean(highd$mean_BA_range),
            mean(highd$mean_BA_range) / mean(lowd$mean_BA_range)))
cat(sprintf("  mean_min_jaccard: d<=5 mean=%.4f  d>=20 mean=%.4f\n",
            mean(lowd$mean_min_jaccard), mean(highd$mean_min_jaccard)))
cat(sprintf("  -> alpha effect %s with dimension\n",
            if (mean(highd$mean_BA_range) > mean(lowd$mean_BA_range)) "GROWS" else "SHRINKS"))

cat("\n=== (D) is the effect large relative to replicate noise? ===\n")
print(res[order(res$d), c("setting","d","method","mean_BA","mean_BA_range",
                          "BA_sd_across_reps_fixed_alpha","alpha_effect_vs_noise")],
      row.names = FALSE, digits = 4)
cat(sprintf("\n  median alpha_effect_vs_noise: d<=5 %.2f   d>=20 %.2f\n",
            median(lowd$mean_BA_range / lowd$BA_sd_across_reps_fixed_alpha),
            median(highd$mean_BA_range / highd$BA_sd_across_reps_fixed_alpha)))

# cross-check against 48's summary
if (file.exists(SUMM)) {
  sm <- read.csv(SUMM, stringsAsFactors = FALSE)
  if ("mean_BA_range" %in% names(sm)) {
    cat(sprintf("\n  48 summary rows=%d; max mean_BA_range 48=%.4f vs recomputed=%.4f -> %s\n",
                nrow(sm), max(sm$mean_BA_range), max(res$mean_BA_range),
                if (abs(max(sm$mean_BA_range) - max(res$mean_BA_range)) < 1e-6) "AGREES" else "DISAGREES"))
  }
}
cat("\ndone\n")
