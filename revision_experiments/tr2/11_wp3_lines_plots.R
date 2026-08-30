#!/usr/bin/env Rscript
# revision_experiments/tr2/11_wp3_lines_plots.R
#
# WP3 threshold-sensitivity plots (Task T7 Phase A).
# Reads whichever of the two sensitivity CSVs exist and renders 300-dpi PNG
# line plots into revision_experiments/results/tr2/figures/:
#   - wp3_sensitivity_synthetic.csv -> wp3_synthetic_<setting>.png
#       x = cutoff multiplier, y = BA / F2 (mean over reps), one line per
#       OS method, dashed vertical line at multiplier = 1 (the calibrated
#       cutoff). If the CSV holds only probe-scale reps (< 100), the
#       subtitle flags the figure as a placeholder pending the production
#       sweep.
#   - wp3_sensitivity_real.csv -> wp3_real_<dataset>.png
#       x = absolute cutoff, dashed vertical line at the calibrated value 2.
# Rerun after the synthetic production run and/or Musk scores land; it
# regenerates every figure it has data for and skips the rest with a log
# line.
#
# Visual conventions: fixed method -> color assignment (colorblind-checked
# palette, Okabe-Ito derived; validated for CVD separation and surface
# contrast), plus linetype and point shape as secondary encoding so method
# identity is never carried by color alone (also survives grayscale print).
#
# Run:  Rscript "revision_experiments/tr2/11_wp3_lines_plots.R"

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

RESULTS_DIR <- here::here("revision_experiments/results/tr2")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
SYNTH_CSV   <- file.path(RESULTS_DIR, "wp3_sensitivity_synthetic.csv")
REAL_CSV    <- file.path(RESULTS_DIR, "wp3_sensitivity_real.csv")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
METHOD_COLORS <- c("RKCCD-OOS" = "#0072B2", "RKCCD-IOS" = "#D55E00",
                   "UNCCD-OOS" = "#009E73", "UNCCD-IOS" = "#B0567F")
METHOD_LINETYPES <- c("RKCCD-OOS" = "solid", "RKCCD-IOS" = "longdash",
                      "UNCCD-OOS" = "dotdash", "UNCCD-IOS" = "twodash")
METHOD_SHAPES <- c("RKCCD-OOS" = 16, "RKCCD-IOS" = 17,
                   "UNCCD-OOS" = 15, "UNCCD-IOS" = 18)

theme_wp3 <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey30", size = 9),
      axis.title = element_text(size = 10)
    )
}

scales_wp3 <- function() {
  list(
    scale_color_manual(values = METHOD_COLORS, breaks = OS_METHODS),
    scale_linetype_manual(values = METHOD_LINETYPES, breaks = OS_METHODS),
    scale_shape_manual(values = METHOD_SHAPES, breaks = OS_METHODS)
  )
}

save_png <- function(p, filename, width = 7.5, height = 4.6) {
  path <- file.path(FIG_DIR, filename)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  cat(sprintf("  wrote %s\n", path))
}

# Reshape one wide (BA_*, F2_*) frame into long metric rows for faceting.
to_long <- function(df, value_cols, sd_cols = NULL) {
  out <- do.call(rbind, lapply(seq_along(value_cols), function(i) {
    metric <- names(value_cols)[i]
    d <- df
    d$metric <- metric
    d$value <- df[[value_cols[i]]]
    d$value_sd <- if (!is.null(sd_cols)) df[[sd_cols[i]]] else NA_real_
    d
  }))
  out$metric <- factor(out$metric, levels = names(value_cols))
  out
}

# ---------------------------------------------------------------------------
# Synthetic figures
# ---------------------------------------------------------------------------
if (file.exists(SYNTH_CSV)) {
  syn <- read.csv(SYNTH_CSV, stringsAsFactors = FALSE)
  syn$method <- factor(syn$method, levels = OS_METHODS)
  max_reps <- max(syn$n_reps)
  probe_note <- if (max_reps < 100) {
    sprintf("PROBE DATA (n_reps = %d) -- placeholder pending the production sweep", max_reps)
  } else {
    sprintf("Monte Carlo means over %d replicates", max_reps)
  }

  for (s in unique(syn$setting)) {
    d <- syn[syn$setting == s, ]
    dl <- to_long(d, c(BA = "BA_mean", F2 = "F2_mean"), c("BA_sd", "F2_sd"))
    p <- ggplot(dl, aes(x = multiplier, y = value, color = method,
                        linetype = method, shape = method)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.4) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.9) +
      facet_wrap(~metric, nrow = 1) +
      scales_wp3() +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = sprintf("Cutoff sensitivity, synthetic %s setting (d = 10)", s),
        subtitle = paste0(probe_note,
                          "\ndashed line = calibrated cutoff (multiplier 1); calibrated values differ per method"),
        x = "Cutoff multiplier (x calibrated value)", y = NULL
      ) +
      theme_wp3()
    save_png(p, sprintf("wp3_synthetic_%s.png", s))
  }
} else {
  cat(sprintf("[skip] %s not found -- no synthetic figures.\n", basename(SYNTH_CSV)))
}

# ---------------------------------------------------------------------------
# Real-data figures
# ---------------------------------------------------------------------------
if (file.exists(REAL_CSV)) {
  rea <- read.csv(REAL_CSV, stringsAsFactors = FALSE)
  rea$method <- factor(rea$method, levels = OS_METHODS)

  for (ds in unique(rea$dataset)) {
    d <- rea[rea$dataset == ds, ]
    calib <- d$calibrated[1]
    dl <- to_long(d, c(BA = "BA", F2 = "F2"))
    p <- ggplot(dl, aes(x = cutoff, y = value, color = method,
                        linetype = method, shape = method)) +
      geom_vline(xintercept = calib, linetype = "dashed", color = "grey45", linewidth = 0.4) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.9) +
      facet_wrap(~metric, nrow = 1) +
      scales_wp3() +
      coord_cartesian(ylim = c(0, 1)) +
      labs(
        title = sprintf("Cutoff sensitivity, %s (n = %d, d = %d)", ds, d$n[1], d$d[1]),
        subtitle = sprintf("Dashed line = calibrated real-data cutoff (%.0f)", calib),
        x = "Score cutoff", y = NULL
      ) +
      theme_wp3()
    save_png(p, sprintf("wp3_real_%s.png", ds))
  }
} else {
  cat(sprintf("[skip] %s not found -- no real-data figures.\n", basename(REAL_CSV)))
}

cat("11_wp3_lines_plots.R done.\n")
