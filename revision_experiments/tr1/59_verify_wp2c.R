#!/usr/bin/env Rscript
# 59_verify_wp2c.R -- independent verification of WP2(c).
#
# WP2(c) recommends replacing the oracle S_min rule with a fixed constant, and
# reports that one published number moves. That is the only sanctioned
# methodological change in this revision cycle, so the numbers behind it are
# re-derived here from the raw per-cell simulation file, not from
# wp2c_simulation_summary.csv.
#
# Claims under test:
#   (A) Pairing is real. The oracle-vs-label-free delta is only meaningful if
#       both rules ran on the SAME replicate (same setting, method, rep/seed).
#       If the rules were run on different draws the "paired SD" is wrong and
#       every conclusion softens. Check the join before trusting any delta.
#   (B) Per-rule mean paired dBA against the oracle, its SD, and the count of
#       cells where the rule is worse.
#   (C) frac_ncls_1 by rule -- the claim that the single-cluster collapse is
#       caused by the S_min VALUE rather than by the method.
#   (D) The split by method at 0.0625: SUN-MCCD reported 0/884, SU-MCCD 26.9%.
#       This is what decides whether WP7 survives for the headline method.
#   (E) The recommended rule's worst cell.
#
# Read-only. Writes nothing.

suppressMessages(library(here))
ARM <- here::here("revision_experiments/results/tr1/wp2c_simulation_arm.csv")
SUM <- here::here("revision_experiments/results/tr1/wp2c_simulation_summary.csv")

d <- read.csv(ARM, stringsAsFactors = FALSE)
d <- d[d$status == "ok", ]
cat(sprintf("rows ok = %d | settings = %d | rules = %d | methods = %s\n",
            nrow(d), length(unique(d$setting_id)), length(unique(d$rule)),
            paste(sort(unique(d$method)), collapse = ",")))

ORACLE <- "oracle_half_contamination"
stopifnot(ORACLE %in% d$rule)

# ---- (A) is the design genuinely paired? ----------------------------------
cat("\n=== (A) pairing check ===\n")
key <- c("setting_id", "method", "rep")
orc <- d[d$rule == ORACLE, ]
cat(sprintf("  oracle rows: %d ; unique (setting,method,rep): %d -> %s\n",
            nrow(orc), nrow(unique(orc[key])),
            if (nrow(orc) == nrow(unique(orc[key]))) "no duplicate keys" else "DUPLICATE KEYS"))
# do the other rules use the same replicate seeds?
mism <- 0L
for (rl in setdiff(unique(d$rule), ORACLE)) {
  s <- d[d$rule == rl, ]
  m <- merge(s[c(key, "seed")], orc[c(key, "seed")], by = key, suffixes = c("_r", "_o"))
  if (nrow(m) == 0) { cat(sprintf("  %-28s NO OVERLAP with oracle\n", rl)); mism <- mism + 1L; next }
  bad <- sum(m$seed_r != m$seed_o)
  cat(sprintf("  %-28s matched %4d/%4d cells, seed mismatches %d\n",
              rl, nrow(m), nrow(s), bad))
  if (bad > 0) mism <- mism + 1L
}
cat(sprintf("  -> %s\n", if (mism == 0) "PAIRED ON IDENTICAL SEEDS" else "PAIRING SUSPECT"))

# ---- (B) paired deltas ------------------------------------------------------
cat("\n=== (B) paired delta vs oracle, recomputed ===\n")
res <- list()
for (rl in setdiff(unique(d$rule), ORACLE)) {
  s <- d[d$rule == rl, ]
  m <- merge(s[c(key, "BA", "F2")], orc[c(key, "BA", "F2")], by = key,
             suffixes = c("_r", "_o"))
  m$dBA <- m$BA_r - m$BA_o; m$dF2 <- m$F2_r - m$F2_o
  # per (setting, method) cell means, then summarise across the 18 cells
  cellmean <- aggregate(cbind(dBA, dF2) ~ setting_id + method, data = m, FUN = mean)
  res[[length(res) + 1]] <- data.frame(
    rule = rl, n_paired = nrow(m), n_cells = nrow(cellmean),
    mean_dBA = mean(m$dBA), sd_dBA = sd(m$dBA),
    mean_dF2 = mean(m$dF2),
    cells_worse = sum(cellmean$dBA < -1e-9),
    worst_cell_dBA = min(cellmean$dBA),
    worst_cell = cellmean$setting_id[which.min(cellmean$dBA)],
    worst_method = cellmean$method[which.min(cellmean$dBA)],
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, res)
res <- res[order(-res$mean_dBA), ]
print(res, row.names = FALSE, digits = 4)

# ---- (C) collapse fraction by rule -----------------------------------------
cat("\n=== (C) fraction of runs returning a single cluster, by rule ===\n")
cl <- aggregate(cbind(frac_ncls_1 = n_clusters == 1) ~ rule, data = d, FUN = mean)
cl$n <- aggregate(n_clusters ~ rule, data = d, FUN = length)$n_clusters
cl$k_hat_correct <- aggregate(cbind(ok = n_clusters == 2) ~ rule, data = d, FUN = mean)$ok
print(cl[order(cl$frac_ncls_1), ], row.names = FALSE, digits = 4)

# ---- (D) the split that decides WP7 ----------------------------------------
cat("\n=== (D) collapse at fixed_0.0625, by method ===\n")
f <- d[d$rule == "fixed_0.0625", ]
if (nrow(f)) {
  by <- aggregate(cbind(frac1 = n_clusters == 1) ~ method, data = f, FUN = mean)
  by$n <- aggregate(n_clusters ~ method, data = f, FUN = length)$n_clusters
  print(by, row.names = FALSE, digits = 4)
  cat("\n  SU-MCCD by setting:\n")
  su <- f[f$method == "SU-MCCD", ]
  bs <- aggregate(cbind(frac1 = n_clusters == 1) ~ setting_id, data = su, FUN = mean)
  bs$n <- aggregate(n_clusters ~ setting_id, data = su, FUN = length)$n_clusters
  bs$k <- aggregate(k ~ setting_id, data = su, FUN = function(x) x[1])$k
  print(bs[order(-bs$frac1), ], row.names = FALSE, digits = 4)
}

# ---- (E) the recommended rule's worst cell ---------------------------------
cat("\n=== (E) recommended rule fixed_0.0500, worst cell detail ===\n")
r5 <- res[res$rule == "fixed_0.0500", ]
if (nrow(r5)) {
  cat(sprintf("  mean dBA=%.4f  cells worse=%d/%d  worst=%.4f at %s / %s\n",
              r5$mean_dBA, r5$cells_worse, r5$n_cells, r5$worst_cell_dBA,
              r5$worst_cell, r5$worst_method))
  s <- d[d$rule == "fixed_0.0500", ]
  m <- merge(s[c(key, "BA")], orc[c(key, "BA")], by = key, suffixes = c("_r", "_o"))
  w <- m[m$setting_id == r5$worst_cell & m$method == r5$worst_method, ]
  cat(sprintf("  worst cell: n_paired=%d  mean dBA=%.4f  sd=%.4f  se=%.4f\n",
              nrow(w), mean(w$BA_r - w$BA_o), sd(w$BA_r - w$BA_o),
              sd(w$BA_r - w$BA_o) / sqrt(nrow(w))))
  sw <- d[d$rule == "fixed_0.0500" & d$setting_id == r5$worst_cell &
            d$method == r5$worst_method, ]
  cat(sprintf("  and its single-cluster fraction: %.3f\n", mean(sw$n_clusters == 1)))
}

# cross-check against the agent's own summary
if (file.exists(SUM)) {
  sm <- read.csv(SUM, stringsAsFactors = FALSE)
  cat(sprintf("\n  summary file rows=%d, columns: %s\n", nrow(sm),
              paste(head(names(sm), 12), collapse = ", ")))
  if ("dBA_mean" %in% names(sm) && "rule" %in% names(sm)) {
    a <- aggregate(dBA_mean ~ rule, data = sm, FUN = mean)
    cmp <- merge(a, res[c("rule", "mean_dBA")], by = "rule")
    cmp$agrees <- abs(cmp$dBA_mean - cmp$mean_dBA) < 5e-3
    print(cmp, row.names = FALSE, digits = 4)
  }
}
cat("\ndone\n")
