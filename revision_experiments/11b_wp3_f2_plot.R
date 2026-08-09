#!/usr/bin/env Rscript
# revision_experiments/11b_wp3_f2_plot.R
#
# Manuscript figure for Appendix C (WP3 cutoff sensitivity).
#
# 11_plots.R renders six figures: BA and F2 side by side, one file per
# synthetic setting and one per real dataset. The appendix now shows only the
# F2 curves for the three synthetic settings, in a single three-panel figure,
# so this script writes that one file. 11_plots.R is left alone; its six
# figures remain the full record behind the appendix.
#
# Colors, linetypes and shapes are copied from 11_plots.R so the appendix
# stays visually consistent with the rest of the WP3 material.
#
# Run:  Rscript "revision_experiments/11b_wp3_f2_plot.R"

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

RESULTS_DIR <- here::here("revision_experiments/results/tr2")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
SYNTH_CSV   <- file.path(RESULTS_DIR, "wp3_sensitivity_synthetic.csv")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
METHOD_COLORS <- c("RKCCD-OOS" = "#0072B2", "RKCCD-IOS" = "#D55E00",
                   "UNCCD-OOS" = "#009E73", "UNCCD-IOS" = "#B0567F")
METHOD_LINETYPES <- c("RKCCD-OOS" = "solid", "RKCCD-IOS" = "longdash",
                      "UNCCD-OOS" = "dotdash", "UNCCD-IOS" = "twodash")
METHOD_SHAPES <- c("RKCCD-OOS" = 16, "RKCCD-IOS" = 17,
                   "UNCCD-OOS" = 15, "UNCCD-IOS" = 18)

SETTING_LABELS <- c(uniform  = "Uniform clusters",
                    gaussian = "Gaussian clusters",
                    matern   = "Random cluster process")

stopifnot(file.exists(SYNTH_CSV))
syn <- read.csv(SYNTH_CSV, stringsAsFactors = FALSE)
stopifnot(all(names(SETTING_LABELS) %in% syn$setting))

syn$method  <- factor(syn$method, levels = OS_METHODS)
syn$setting <- factor(SETTING_LABELS[syn$setting], levels = SETTING_LABELS)

p <- ggplot(syn, aes(x = multiplier, y = F2_mean, color = method,
                     linetype = method, shape = method)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_wrap(~setting, nrow = 1) +
  scale_color_manual(values = METHOD_COLORS, breaks = OS_METHODS) +
  scale_linetype_manual(values = METHOD_LINETYPES, breaks = OS_METHODS) +
  scale_shape_manual(values = METHOD_SHAPES, breaks = OS_METHODS) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Cutoff multiplier (x calibrated value)", y = expression(F[2])) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 10)
  )

out <- file.path(FIG_DIR, "wp3_synthetic_F2.png")
ggsave(out, p, width = 9.5, height = 3.5, dpi = 300, bg = "white")
cat(sprintf("wrote %s (reps = %d per grid point)\n", out, max(syn$n_reps)))
