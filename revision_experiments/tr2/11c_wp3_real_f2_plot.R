#!/usr/bin/env Rscript
# revision_experiments/tr2/11c_wp3_real_f2_plot.R
#
# Companion to 11b_wp3_synthetic_f2_plot.R: the real-data half of the Appendix C
# cutoff-sensitivity figure, as a single three-panel F2 figure placed directly
# below the synthetic one in the manuscript.
#
# Panels: WBC, Thyroid, vowels. vowels replaced Arrhythmia on 2026-08-16 --
# Arrhythmia (d=274) can only be scored by the two UN-CCD methods, so its panel
# carried 2 of 4 curves. See 10_wp3_real.R for why vowels and not WDBC.
#
# Unlike the synthetic figure, the horizontal axis is the cutoff itself, not a
# multiplier: all four scores share the single adopted real-data cutoff of 2
# (manuscript Section 4), so there is nothing to normalize away.
#
# Colors, linetypes and shapes are copied from 11_wp3_lines_plots.R / 11b so all of the
# WP3 material reads the same.
#
# Run:  Rscript "revision_experiments/tr2/11c_wp3_real_f2_plot.R"

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

RESULTS_DIR <- here::here("revision_experiments/results/tr2")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
REAL_CSV    <- file.path(RESULTS_DIR, "wp3_sensitivity_real.csv")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OS_METHODS <- c("RKCCD-OOS", "RKCCD-IOS", "UNCCD-OOS", "UNCCD-IOS")
METHOD_COLORS <- c("RKCCD-OOS" = "#0072B2", "RKCCD-IOS" = "#D55E00",
                   "UNCCD-OOS" = "#009E73", "UNCCD-IOS" = "#B0567F")
METHOD_LINETYPES <- c("RKCCD-OOS" = "solid", "RKCCD-IOS" = "longdash",
                      "UNCCD-OOS" = "dotdash", "UNCCD-IOS" = "twodash")
METHOD_SHAPES <- c("RKCCD-OOS" = 16, "RKCCD-IOS" = 17,
                   "UNCCD-OOS" = 15, "UNCCD-IOS" = 18)

PANELS <- c("WBC", "Thyroid", "vowels")
CALIBRATED <- 2

stopifnot(file.exists(REAL_CSV))
rea <- read.csv(REAL_CSV, stringsAsFactors = FALSE)
missing <- setdiff(PANELS, unique(rea$dataset))
if (length(missing) > 0) {
  stop(sprintf("missing from %s: %s -- rerun 10_wp3_real.R",
               basename(REAL_CSV), paste(missing, collapse = ", ")))
}

rea <- rea[rea$dataset %in% PANELS, ]

# Every panel must carry all four scores; the whole point of dropping
# Arrhythmia was to stop showing a panel with two.
for (ds in PANELS) {
  got <- sort(unique(rea$method[rea$dataset == ds]))
  if (!setequal(got, OS_METHODS)) {
    stop(sprintf("%s carries %d of 4 scores (%s)", ds, length(got),
                 paste(got, collapse = ", ")))
  }
}

# Facet label carries n and d, so the figure states its own settings.
lab <- vapply(PANELS, function(ds) {
  r <- rea[rea$dataset == ds, ][1, ]
  sprintf("%s (n = %d, d = %d)", ds, r$n, r$d)
}, character(1))

rea$method  <- factor(rea$method, levels = OS_METHODS)
rea$dataset <- factor(lab[rea$dataset], levels = unname(lab))

p <- ggplot(rea, aes(x = cutoff, y = F2, color = method,
                     linetype = method, shape = method)) +
  geom_vline(xintercept = CALIBRATED, linetype = "dashed", color = "grey45",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_manual(values = METHOD_COLORS, breaks = OS_METHODS) +
  scale_linetype_manual(values = METHOD_LINETYPES, breaks = OS_METHODS) +
  scale_shape_manual(values = METHOD_SHAPES, breaks = OS_METHODS) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Cutoff", y = expression(F[2])) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 10)
  )

out <- file.path(FIG_DIR, "wp3_real_F2.png")
ggsave(out, p, width = 9.5, height = 3.5, dpi = 300, bg = "white")
cat(sprintf("wrote %s (%d cutoffs per curve, panels: %s)\n", out,
            length(unique(rea$cutoff)), paste(PANELS, collapse = ", ")))
