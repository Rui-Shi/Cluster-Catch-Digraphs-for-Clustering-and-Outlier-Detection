# 25_smin_match.R -- which S_min value reproduces the PUBLISHED real-data rows?
#
# The manuscript specifies S_min only for the simulations (main text L605,
# "half the contamination level", scoped to the simulation settings). It is
# never specified for the real data. Agent A's sweep (21_regen_smin_grid.R)
# gives us the metric as a step function of S_min on five datasets, so we can
# invert the question: for each (dataset, shape-adaptive method), which grid
# values of S_min land on the published row?
#
# If a single constant reproduces several published cells while the
# contamination-derived reading does not, that constant is what the real-data
# runs actually used, and it is what the revision must document.
#
# Read-only. Writes nothing.

suppressMessages(library(here))

GRID  <- read.csv(here::here("revision_experiments/results/tr1/regen_smin_grid.csv"),
                  stringsAsFactors = FALSE)
TRUTH <- read.csv(here::here("revision_experiments/published_realdata_truth.csv"),
                  stringsAsFactors = FALSE)

published <- function(ds, meth) {
  s <- TRUTH[tolower(TRUTH$dataset) == tolower(ds) & TRUTH$method == meth, ]
  if (nrow(s) == 0) return(NULL)
  setNames(as.numeric(s$value[match(c("TPR", "TNR", "BA", "F2"), s$metric)]),
           c("TPR", "TNR", "BA", "F2"))
}

# tolerance: 3 d.p. exact, plus a 1-unit-in-the-last-place band, because
# several published cells are one true negative away from ours and we want to
# see near-misses rather than silently call them mismatches.
cat(sprintf("\n%-11s %-9s %-9s %7s %7s %7s %7s   %s\n",
            "dataset", "method", "s_min", "TPR", "TNR", "BA", "F2", "vs published"))
cat(strrep("-", 88), "\n")

hits <- list()
for (ds in sort(unique(GRID$dataset))) {
  for (meth in sort(unique(GRID$method))) {
    sub <- GRID[GRID$dataset == ds & GRID$method == meth, ]
    if (nrow(sub) == 0) next
    pub <- published(ds, meth)
    if (is.null(pub)) next
    for (i in order(sub$variant)) {
      r <- sub[i, ]
      got <- c(TPR = r$TPR, TNR = r$TNR, BA = r$BA, F2 = r$F2)
      dif <- max(abs(round(got, 3) - pub))
      tag <- if (dif < 5e-4) "EXACT" else if (dif <= 1.5e-3) "within 0.0015" else ""
      if (nzchar(tag)) {
        cat(sprintf("%-11s %-9s %-9s %7.3f %7.3f %7.3f %7.3f   %s (max diff %.4f)\n",
                    ds, meth, r$variant, got[1], got[2], got[3], got[4], tag, dif))
        hits[[length(hits) + 1]] <- data.frame(dataset = ds, method = meth,
                                               variant = r$variant, diff = dif,
                                               stringsAsFactors = FALSE)
      }
    }
  }
}

if (length(hits) == 0) {
  cat("  (no grid point reproduces any published row)\n")
} else {
  h <- do.call(rbind, hits)
  cat("\n--- how many published cells each S_min value reproduces ---\n")
  print(sort(table(h$variant), decreasing = TRUE))
}

# And the converse: which (dataset, method) pairs are reproduced by NO grid point?
cat("\n--- published cells reproduced by NO S_min in the grid ---\n")
for (ds in sort(unique(GRID$dataset))) {
  for (meth in sort(unique(GRID$method))) {
    sub <- GRID[GRID$dataset == ds & GRID$method == meth, ]
    if (nrow(sub) == 0) next
    pub <- published(ds, meth)
    if (is.null(pub)) next
    best <- min(apply(sub[, c("TPR", "TNR", "BA", "F2")], 1,
                      function(g) max(abs(round(as.numeric(g), 3) - pub))))
    if (best > 1.5e-3) {
      cat(sprintf("  %-11s %-9s   best any grid point can do: max diff %.3f\n",
                  ds, meth, best))
    }
  }
}
cat("\ndone\n")
