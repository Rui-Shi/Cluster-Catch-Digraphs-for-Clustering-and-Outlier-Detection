#!/usr/bin/env Rscript
# 95_wp6_incr_summary.R -- summaries of the paired UN-/SUN-MCCD timings.
#
# The UN-MCCD and SUN-MCCD rows of the runtime tables (main text: the runtime
# table; supplement: the n- and d-sweep tables and the log-log slope table)
# report the incremental nearest-neighbour radius search, measured in the
# paired run of 94_wp6_ab_old_vs_new.R (arm impl == "new", part == "detector").
# The U-MCCD, SU-MCCD and baseline rows come from 91_wp6_runtime.R
# --summarize (results/tr1/wp6/). This script derives the UN/SUN rows with the
# same definitions 91 uses: median and CV (sd/mean) of time_s and of
# mem_delta_mb per cell, and lm(log(y) ~ log(n)) over every replicate of the
# d = 10 sweep for the slopes. It also summarizes the ascending-vs-descending
# timing of 96_wp6_direction_timing.R.
#
# Reads:  results/tr1/wp6_incr/94_ab_raw.csv, 96_direction_raw.csv,
#         and (table-load column only) wp6_incr/91_wp6_runtime_{n,d}.csv
# Writes: results/tr1/wp6_incr/95_ab_runtime_n.csv   (d = 10 sweep)
#         results/tr1/wp6_incr/95_ab_runtime_d.csv   (n = 500 sweep)
#         results/tr1/wp6_incr/95_ab_slope.csv
#         results/tr1/wp6_incr/95_direction_summary.csv
# Usage:  Rscript revision_experiments/tr1/95_wp6_incr_summary.R

suppressPackageStartupMessages(library(here))
DIR <- here::here("revision_experiments/results/tr1/wp6_incr")

ab <- read.csv(file.path(DIR, "94_ab_raw.csv"), stringsAsFactors = FALSE)
ab <- ab[ab$impl == "new" & ab$part == "detector", ]
stopifnot(nrow(ab) == 160L)

cv <- function(x) if (length(x) > 1) stats::sd(x) / mean(x) else NA_real_

summ <- function(sub) {
  keys <- unique(sub[, c("n", "d", "method")])
  keys <- keys[order(keys$n, keys$d, keys$method), ]
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    g <- sub[sub$n == k$n & sub$d == k$d & sub$method == k$method, ]
    data.frame(n = k$n, d = k$d, method = k$method, n_reps = nrow(g),
               median_time_s = stats::median(g$time_s), mean_time_s = mean(g$time_s),
               cv_time = cv(g$time_s),
               median_mem_delta_mb = stats::median(g$mem_delta_mb),
               cv_mem_delta = cv(g$mem_delta_mb))
  })
  do.call(rbind, rows)
}

# table-load medians: not timed in the paired run; taken from the 91 re-timing
# of the same code (wp6_incr/91_*), where they are recorded per cell
add_table_load <- function(agg, f) {
  p <- file.path(DIR, f)
  if (!file.exists(p)) { agg$median_table_load_s <- NA_real_; return(agg) }
  t <- read.csv(p, stringsAsFactors = FALSE)[, c("n", "d", "method", "median_table_load_s")]
  merge(agg, t, by = c("n", "d", "method"), all.x = TRUE, sort = FALSE)
}

agg_n <- add_table_load(summ(ab[ab$d == 10, ]), "91_wp6_runtime_n.csv")
agg_d <- add_table_load(summ(ab[ab$n == 500, ]), "91_wp6_runtime_d.csv")
agg_n <- agg_n[order(agg_n$n, agg_n$method), ]
agg_d <- agg_d[order(agg_d$d, agg_d$method), ]
write.csv(agg_n, file.path(DIR, "95_ab_runtime_n.csv"), row.names = FALSE)
write.csv(agg_d, file.path(DIR, "95_ab_runtime_d.csv"), row.names = FALSE)

fit <- function(m, metric, sub, y) {
  f <- stats::lm(log(y) ~ log(n), data = sub)
  cf <- summary(f)$coefficients
  data.frame(method = m, metric = metric, n_points = nrow(sub),
             slope = cf["log(n)", "Estimate"], se = cf["log(n)", "Std. Error"],
             r_squared = summary(f)$r.squared,
             stated_bound = if (metric == "time") "n^3" else "n^2")
}
slope <- do.call(rbind, lapply(c("UN-MCCD", "SUN-MCCD"), function(m) {
  s <- ab[ab$d == 10 & ab$method == m, ]
  sm <- s[is.finite(s$mem_delta_mb) & s$mem_delta_mb > 0, ]
  rbind(fit(m, "time", s, s$time_s), fit(m, "mem", sm, sm$mem_delta_mb))
}))
write.csv(slope, file.path(DIR, "95_ab_slope.csv"), row.names = FALSE)

dr <- read.csv(file.path(DIR, "96_direction_raw.csv"), stringsAsFactors = FALSE)
keys <- unique(dr[, c("method", "n", "d", "direction")])
dsum <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  g <- dr[dr$method == k$method & dr$n == k$n & dr$d == k$d & dr$direction == k$direction, ]
  data.frame(k, n_reps = nrow(g), median_time_s = stats::median(g$time_s),
             cv_time = cv(g$time_s), median_mem_delta_mb = stats::median(g$mem_delta_mb))
}))
write.csv(dsum, file.path(DIR, "95_direction_summary.csv"), row.names = FALSE)

fmt <- function(df) print(format(df, digits = 4), row.names = FALSE)
cat("== d = 10 sweep\n"); fmt(agg_n)
cat("== n = 500 sweep\n"); fmt(agg_d)
cat("== slopes\n"); fmt(slope)
cat("== search direction\n"); fmt(dsum)
