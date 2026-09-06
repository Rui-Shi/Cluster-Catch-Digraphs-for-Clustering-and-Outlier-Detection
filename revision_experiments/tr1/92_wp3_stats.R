#!/usr/bin/env Rscript
# 92_wp3_stats.R -- WP3: paired real-data significance tests, simulation
# variability (to the extent per-replicate data exist), and the parameter
# provenance table, following tr1/WP3_PROTOCOL.md.
#
# This script computes NO new metrics. Part (a) reuses the exact per-cell
# real-data table tr1/82_wp4_metrics.R already assembled (fx = the corrected
# nine-method table; combined("T1") = fx plus the eight WP4 competitors at
# the contamination-0.1 threshold) -- the same numbers already printed in
# tab:Real_Data_Aggregate and tab:Real_Data_Result_Summary. Part (b) reads
# whatever per-replicate simulation output exists (see WP3_PROTOCOL.md (b)
# for the search that found none for the manuscript's own simulation grid).
# Part (c) writes out the parameter-provenance table drafted in the protocol.
#
# 82_wp4_metrics.R is sourced, not edited or re-derived, per instruction.

suppressMessages(library(here))
options(width = 200)

TR1 <- here::here("revision_experiments/tr1")
WP3 <- here::here("revision_experiments/results/tr1/wp3")
dir.create(WP3, recursive = TRUE, showWarnings = FALSE)

cat("=== WP3: sourcing 82_wp4_metrics.R for the shared per-cell assembly ===\n")
E82 <- new.env()
invisible(capture.output(source(file.path(TR1, "82_wp4_metrics.R"), local = E82)))
cat("  82_wp4_metrics.R sourced OK (its own internal PRINTED gate already passed,\n")
cat("  or source() would have stopped with an error).\n")

# ---------------------------------------------------------------------------
# GATE (re-asserted here, not just trusted from 82's internal stop())
#
# NOTE: 82_wp4_metrics.R's own internal PRINTED constant and its a_pub =
# agg(fc) check what was printed in the manuscript BEFORE the alpha-repair
# values (CLAUDE.md, "tab:alpha_real needs correcting") were folded into the
# main text. The manuscript has since been edited (\revised{} spans visible
# today in tab:Real_Data_Aggregate for SUN-MCCD and UN-MCCD's mean/median F2
# and BA), so 82's own PRINTED is now STALE relative to the actual current
# tab:Real_Data_Aggregate. WP3's gate must check against the manuscript as it
# stands today, not against 82's authoring-time snapshot, so the 9 values
# below were copied by hand from CCD_OutlierDetection_Neurocomputing.tex
# (both blocks of tab:Real_Data_Aggregate) on 2026-09-05, and the comparison
# is against E82$fx (alpha-corrected), which is what those printed values
# actually equal -- NOT against E82$fc (pre-alpha-repair), which is what
# 82's own internal gate uses. This discrepancy is reported to the verifier.
# ---------------------------------------------------------------------------
cat("\n=== GATE (re-asserted in 92): nine-method aggregate vs the CURRENT tab:Real_Data_Aggregate ===\n")
PRINTED_92 <- data.frame(
  method    = c("SUN-MCCD","SU-MCCD","U-MCCD","iForest","LOF","MST","UN-MCCD","ODIN","DBSCAN"),
  mean_F2   = c(0.317, 0.302, 0.291, 0.260, 0.254, 0.237, 0.229, 0.229, 0.234),
  median_F2 = c(0.329, 0.301, 0.280, 0.115, 0.195, 0.237, 0.216, 0.224, 0.133),
  mean_BA   = c(0.718, 0.696, 0.709, 0.632, 0.654, 0.606, 0.659, 0.619, 0.599),
  median_BA = c(0.713, 0.692, 0.708, 0.526, 0.577, 0.599, 0.641, 0.585, 0.546),
  stringsAsFactors = FALSE)
a_pub_92 <- E82$agg(E82$fx)
mm92 <- merge(PRINTED_92, a_pub_92, by = "method", suffixes = c("_printed", "_computed"))
stopifnot(nrow(mm92) == 9)
bad92 <- 0L
for (col in c("mean_F2", "median_F2", "mean_BA", "median_BA")) {
  dd <- abs(mm92[[paste0(col, "_printed")]] - mm92[[paste0(col, "_computed")]])
  for (i in which(dd > 0.0005)) {
    bad92 <- bad92 + 1L
    cat(sprintf("  MISMATCH %-9s %-10s printed %.3f, computed %.3f\n",
                mm92$method[i], col, mm92[[paste0(col, "_printed")]][i],
                mm92[[paste0(col, "_computed")]][i]))
  }
}
if (bad92 == 0L) {
  cat("  all 36 printed values (as they stand in the manuscript TODAY) reproduced exactly -- gate PASSED\n")
} else {
  stop(sprintf("%d printed values do not reproduce -- WP3 stats are not trustworthy, stopping", bad92))
}
cat("\n  NOTE: 82_wp4_metrics.R's OWN internal PRINTED gate (agg(fc) vs its hardcoded\n")
cat("  PRINTED constant: SUN-MCCD mean_F2=0.321, mean_BA=0.722, etc.) is STALE -- it\n")
cat("  reproduces an earlier, pre-alpha-repair snapshot of the manuscript, not the\n")
cat("  \\revised{} values now in tab:Real_Data_Aggregate. This does not affect WP3's\n")
cat("  own results (which gate against the current text, above), but is flagged here\n")
cat("  because 82 is read-only for WP3 and was not touched to fix it.\n")

# ---------------------------------------------------------------------------
# (a) paired real-data tests
# ---------------------------------------------------------------------------
ORDER <- E82$ORDER
PROPOSED <- E82$PROPOSED
c1 <- E82$combined("T1")   # 17 methods x 16 data sets, F2/BA already at published rounding
stopifnot(nrow(c1) == 17 * 16)
METHODS <- sort(unique(c1$method))
stopifnot(length(METHODS) == 17, "SUN-MCCD" %in% METHODS)

wide <- function(metric) {
  m <- matrix(NA_real_, nrow = length(ORDER), ncol = length(METHODS),
              dimnames = list(ORDER, METHODS))
  for (i in seq_len(nrow(c1))) m[c1$dataset[i], c1$method[i]] <- c1[[metric]][i]
  stopifnot(!anyNA(m))
  m
}
F2W <- wide("F2"); BAW <- wide("BA")

one_test <- function(primary, opponent, metric, W) {
  x <- W[, primary]; y <- W[, opponent]
  d <- x - y
  n_ahead  <- sum(d > 0)
  n_behind <- sum(d < 0)
  n_tied   <- sum(d == 0)
  wt <- tryCatch(
    suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = NULL)),
    error = function(e) NULL
  )
  if (is.null(wt)) {
    # every difference is zero -- wilcox.test errors ("cannot compute exact
    # p-value with ties" / "not enough finite observations"); report as a
    # degenerate all-tied comparison rather than crash.
    return(data.frame(primary = primary, opponent = opponent, metric = metric,
                       n = length(x), n_ahead = n_ahead, n_behind = n_behind,
                       n_tied = n_tied, V = NA_real_, p_raw = NA_real_,
                       exact_used = NA, note = "all 16 differences are zero; no test possible",
                       stringsAsFactors = FALSE))
  }
  exact_used <- grepl("exact", wt$method, ignore.case = TRUE)
  data.frame(primary = primary, opponent = opponent, metric = metric,
             n = length(x), n_ahead = n_ahead, n_behind = n_behind,
             n_tied = n_tied, V = unname(wt$statistic), p_raw = wt$p.value,
             exact_used = exact_used, note = "", stringsAsFactors = FALSE)
}

run_block <- function(primary, metric, W) {
  opponents <- setdiff(METHODS, primary)
  rows <- do.call(rbind, lapply(opponents, one_test, primary = primary, metric = metric, W = W))
  rows$p_holm <- p.adjust(rows$p_raw, method = "holm")
  rows[order(rows$p_raw), ]
}

PRIMARY_MAIN <- "SUN-MCCD"
PRIMARY_SECONDARY <- c("U-MCCD", "SU-MCCD", "UN-MCCD")

cat("\n=== (a) paired Wilcoxon: SUN-MCCD vs the other 16 methods, F2 ===\n")
sun_f2 <- run_block(PRIMARY_MAIN, "F2", F2W)
print(sun_f2[, c("primary","opponent","n","n_ahead","n_behind","n_tied","V","p_raw","p_holm","exact_used")], row.names = FALSE)

cat("\n=== (a) paired Wilcoxon: SUN-MCCD vs the other 16 methods, BA ===\n")
sun_ba <- run_block(PRIMARY_MAIN, "BA", BAW)
print(sun_ba[, c("primary","opponent","n","n_ahead","n_behind","n_tied","V","p_raw","p_holm","exact_used")], row.names = FALSE)

cat("\n=== (a) secondary: U-MCCD / SU-MCCD / UN-MCCD vs the other 16 methods ===\n")
secondary <- do.call(rbind, lapply(PRIMARY_SECONDARY, function(pm) {
  rbind(run_block(pm, "F2", F2W), run_block(pm, "BA", BAW))
}))
for (pm in PRIMARY_SECONDARY) for (met in c("F2","BA")) {
  cat(sprintf("\n--- %s vs 16, %s ---\n", pm, met))
  g <- secondary[secondary$primary == pm & secondary$metric == met, ]
  print(g[, c("opponent","n","n_ahead","n_behind","n_tied","V","p_raw","p_holm","exact_used")], row.names = FALSE)
}

wilcoxon_real <- rbind(sun_f2, sun_ba, secondary)
wilcoxon_real <- wilcoxon_real[, c("primary","metric","opponent","n","n_ahead","n_behind",
                                    "n_tied","V","p_raw","p_holm","exact_used","note")]
write.csv(wilcoxon_real, file.path(WP3, "wilcoxon_real.csv"), row.names = FALSE)

sig_holm <- function(rows) sum(rows$p_holm < 0.05, na.rm = TRUE)
cat(sprintf("\nSUN-MCCD: Holm-significant (p_holm<0.05) on %d/16 F2 comparisons, %d/16 BA comparisons\n",
            sig_holm(sun_f2), sig_holm(sun_ba)))

# ---------------------------------------------------------------------------
# (b) simulation variability
# ---------------------------------------------------------------------------
cat("\n=== (b) simulation-variability: searching for per-replicate output ===\n")
SIM_ROOT <- here::here("simulations/outlier_detection")
all_files <- list.files(SIM_ROOT, recursive = TRUE, full.names = FALSE)
ext <- function(pat) sum(grepl(pat, all_files, ignore.case = TRUE))
n_R    <- ext("\\.R$")
n_out  <- ext("\\.out$")
n_csv  <- ext("\\.csv$")
n_rds  <- ext("\\.rds$")
n_rdat <- ext("\\.RData$")
cat(sprintf("  %s : %d files total\n", SIM_ROOT, length(all_files)))
cat(sprintf("  .R scripts: %d | .out SLURM logs: %d | .csv: %d | .rds: %d | .RData: %d\n",
            n_R, n_out, n_csv, n_rds, n_rdat))

wp2c_path <- here::here("revision_experiments/results/tr1/wp2c_simulation_arm.csv")
wp2c_exists <- file.exists(wp2c_path)
inv_lines <- c(
  "# WP3 simulation-variability inventory",
  "",
  "Generated by `tr1/92_wp3_stats.R`, following the search declared in",
  "`tr1/WP3_PROTOCOL.md` (b).",
  "",
  sprintf("Searched `%s` recursively: **%d files**.", SIM_ROOT, length(all_files)),
  sprintf("- `.R` driver scripts: %d", n_R),
  sprintf("- `.out` SLURM stdout logs: %d", n_out),
  sprintf("- `.csv` result files: %d", n_csv),
  sprintf("- `.rds` result files: %d", n_rds),
  sprintf("- `.RData` result files: %d", n_rdat),
  "",
  "**Conclusion.** Every driver ends with `save.image()` to an `.RData` file",
  "in the same directory as the script -- that is where the raw",
  "per-replicate TPR/TNR/BA/F2 vectors behind the manuscript's Section 5",
  "simulation tables and figures (uniform clusters, Gaussian clusters, the",
  "six robustness checks, and the Neyman-Scott experiments) would live.",
  "`*.RData` is listed in this repo's `.gitignore`; the file count above",
  "confirms none survive in this checkout. The SLURM `.out` logs print only",
  "a single aggregate line per job (e.g. \"The mean success rate is ...\"),",
  "not a per-replicate vector, so no standard deviation is recoverable from",
  "them either.",
  "",
  "**No standard deviation or Monte Carlo standard error can be computed for",
  "any of the manuscript's own reported simulation means** (main synthetic",
  "grid, six robustness checks, or Neyman-Scott settings) from data present",
  "in this repository.",
  ""
)
if (wp2c_exists) {
  wp2c <- read.csv(wp2c_path, stringsAsFactors = FALSE)
  wp2c <- wp2c[wp2c$status == "ok", ]
  # NOTE: s_min is EXCLUDED from the grouping key. For the fixed-value rules
  # (fixed_0.03/0.05/0.0625, oracle_half_contamination, min_cls_zero) s_min is
  # constant within a (setting_id, method, rule) cell anyway. For
  # "adaptive_one_step_rho0.5" it is a per-replicate COMPUTED quantity that
  # varies slightly rep to rep (e.g. gaussian_d10_n200_c02/SU-MCCD: 95 of 100
  # reps land on s_min=0.25, the rest scattered near it) -- grouping on its
  # exact value would shatter a single 100-replicate design cell into dozens
  # of 1-2-row pseudo-cells. s_min's mean/range is kept as an informational
  # column instead.
  cell_cols <- c("setting_id","generator","d","contam","n_nominal","method","rule","label_free")
  agg_wp2c <- do.call(rbind, lapply(split(wp2c, interaction(wp2c[cell_cols], drop = TRUE)), function(g) {
    if (nrow(g) < 2) return(NULL)
    data.frame(g[1, cell_cols], n_rep = nrow(g),
               mean_s_min = mean(g$s_min), min_s_min = min(g$s_min), max_s_min = max(g$s_min),
               mean_F2 = mean(g$F2), sd_F2 = sd(g$F2), mcse95_F2 = 1.96 * sd(g$F2) / sqrt(nrow(g)),
               mean_TPR = mean(g$TPR), sd_TPR = sd(g$TPR),
               mean_TNR = mean(g$TNR), sd_TNR = sd(g$TNR),
               mean_BA = mean(g$BA), sd_BA = sd(g$BA), mcse95_BA = 1.96 * sd(g$BA) / sqrt(nrow(g)),
               matches_manuscript_cell = FALSE,
               note = "WP2(c) S_min-sensitivity cell; design (d, n, contam, method, rule) does not match any manuscript-reported simulation cell -- see WP3_PROTOCOL.md (b)",
               stringsAsFactors = FALSE)
  }))
  agg_wp2c <- agg_wp2c[order(agg_wp2c$generator, agg_wp2c$d, agg_wp2c$n_nominal, agg_wp2c$method, agg_wp2c$rule), ]
  write.csv(agg_wp2c, file.path(WP3, "sim_variability.csv"), row.names = FALSE)
  cat(sprintf("\n  wp2c_simulation_arm.csv found: %d ok rows -> %d per-cell SD rows written to sim_variability.csv\n",
              nrow(wp2c), nrow(agg_wp2c)))
  cat("  (labeled matches_manuscript_cell = FALSE throughout -- see WP3_PROTOCOL.md (b))\n")
  rng_f2 <- range(agg_wp2c$sd_F2); rng_ba <- range(agg_wp2c$sd_BA)
  inv_lines <- c(inv_lines,
    "**What does exist.** `55_wp2c_simulation_arm.R` (a WP2(c) script) wrote",
    sprintf("genuine per-replicate rows to `%s`", basename(wp2c_path)),
    sprintf("(%d rows with status \"ok\", collapsing to %d per-cell groups).", nrow(wp2c), nrow(agg_wp2c)),
    "It answers a narrower question -- S_min rule sensitivity, at d=10 only,",
    "n in {200,500}, SU-MCCD and SUN-MCCD only, 100 replicates per cell -- and",
    "its cells do not correspond to any manuscript-reported figure or",
    "in-text value (a direct check, gaussian d=10 contam=5% n=200 SUN-MCCD,",
    "gives mean F2 = 0.787 against the manuscript's F2 = 0.830 for the same",
    "(d, generator) combination, confirming the mismatch).",
    sprintf("`92_wp3_stats.R` still computes and writes its per-cell SD/MC-SE to"),
    "`sim_variability.csv`, every row flagged `matches_manuscript_cell = FALSE`,",
    sprintf("as the only real evidence in this repo of Monte Carlo noise magnitude"),
    sprintf("for these methods (sd(F2) ranges %.4f-%.4f, sd(BA) ranges %.4f-%.4f",
            rng_f2[1], rng_f2[2], rng_ba[1], rng_ba[2]),
    "across its cells, at 100 replicates -- an indicative, not a substitute,",
    "figure for the manuscript's own 1000-replicate cells).",
    "")
} else {
  cat("\n  wp2c_simulation_arm.csv not found -- writing inventory only, no sim_variability.csv\n")
  inv_lines <- c(inv_lines, sprintf("`%s` was not found either.", wp2c_path), "")
}
writeLines(inv_lines, file.path(WP3, "sim_variability_INVENTORY.md"))

# ---------------------------------------------------------------------------
# (c) parameter provenance table
# ---------------------------------------------------------------------------
cat("\n=== (c) parameter provenance table ===\n")
prov <- read.csv(textConnection('method,parameter,value,how_set,source
LOF,k,"{11,...,30}; largest LOF value retained per point",fixed constant (literature-recommended envelope),"main text sec:rand-clust-proc, reused unchanged for real data"
LOF,outlier threshold,1.5,fixed constant (matches Breunig et al. reported optimum),main text sec:rand-clust-proc
DBSCAN,MinPts,4,fixed constant,"main text sec:rand-clust-proc and sec:Real-Data-Examples"
DBSCAN,Eps,4-distance heuristic evaluated at each data set own true contamination rate,USES TRUE CONTAMINATION,"main text sec:Real-Data-Examples"
MST,inconsistent-edge threshold,"1.2 (real data); 1.7 to 1.1 decreasing with d (Neyman-Scott simulation)",fixed constant (held across all 16 real data sets),"main text sec:Real-Data-Examples and sec:rand-clust-proc"
MST,minimum cluster size to avoid outlier flag,2% of n,fixed constant,main text sec:Real-Data-Examples
ODIN,k,round(sqrt(n)),fixed rule (function of n; label-free),main text sec:rand-clust-proc
ODIN,in-degree threshold T,round(n^(1/3)),fixed rule (function of n; label-free),main text sec:rand-clust-proc
iForest,number of trees,1000,fixed constant (library-recommended),main text sec:rand-clust-proc
iForest,sub-sample size,256,fixed constant (library-recommended),main text sec:rand-clust-proc
iForest,outlier-score threshold,0.55,fixed constant,main text sec:rand-clust-proc
"U-MCCD; SU-MCCD (RK-based)",MC-SRT significance level alpha,"step function of d: 1% (d<10), 0.1% (d>=10)",fixed rule (function of d; pre-declared; label-free),"main text sec:Uniform_clusters; tab:alpha_real for realized values"
"U-MCCD; SU-MCCD",KS-CCD density-search tolerance tau,0.01,fixed constant (sensitivity checked over 1e-3 to 1; WP2a),Supplementary Material algorithm boxes
"U-MCCD; SU-MCCD",Monte Carlo replications behind RK quantile tables,5000 (99% tables) / 10000 (99.9% tables),"fixed constant, computed once and stored",R/RK-test_quantile/*.R
SU-MCCD,S_min (minimum cluster fraction),0.05,"fixed constant, label-free (same value for every real data set and the simulations)",main text sec:Real-Data-Examples
"UN-MCCD; SUN-MCCD (NND-based)",MC-SRT significance level alpha,"step function of d, own schedule per method",fixed rule (function of d; pre-declared; label-free),"main text sec:Uniform_clusters; tab:alpha_real for realized values"
"UN-MCCD; SUN-MCCD",Monte Carlo replications behind NN quantile tables,5000 / 10000 by level,fixed constant,R/NN-test_quantile/*.R
SUN-MCCD,S_min,0.05,"fixed constant, label-free",main text sec:Real-Data-Examples
ECOD,none,"- (parameter-free by construction)",library default,tr1/WP4_PROTOCOL.md section 2
COPOD,none,"- (parameter-free by construction)",library default,tr1/WP4_PROTOCOL.md section 2
DIF,architecture/training defaults,"pyod 3.6.1 get_params() defaults (hidden_neurons=[500,100] etc.)",library default,tr1/WP4_PROTOCOL.md section 5
DIF,random seed,"{1,...,5}; mean over seeds reported",fixed constant (5 seeds; pre-declared),tr1/WP4_PROTOCOL.md section 2
LUNAR,architecture/training defaults,"pyod 3.6.1 get_params() defaults (n_neighbours=5 etc.)",library default,tr1/WP4_PROTOCOL.md section 5
LUNAR,random seed,"{1,...,5}; mean over seeds reported",fixed constant (5 seeds; pre-declared),tr1/WP4_PROTOCOL.md section 2
GLOSH (HDBSCAN),min_cluster_size,5,library default,tr1/WP4_PROTOCOL.md sections 2 and 5
GLOSH (HDBSCAN),min_samples,None (falls back to min_cluster_size),library default,tr1/WP4_PROTOCOL.md section 5
OPTICS,min_samples,5,library default,tr1/WP4_PROTOCOL.md section 5
OPTICS,xi,0.05,library default,tr1/WP4_PROTOCOL.md section 5
MutualKNN,k,"best of {5,10,15,20,30} by F2 under true labels, per data set",ORACLE-K (uses true labels),"tr1/WP4_PROTOCOL.md section 2, oracle concession on k"
SNN,k,"best of {5,10,15,20,30} by F2 under true labels, per data set",ORACLE-K (uses true labels),tr1/WP4_PROTOCOL.md section 2
"all 8 WP4 competitors (T1 thresholding)",contamination fraction to threshold the score,"0.1, identical for every data set",library default; label-free,"tr1/WP4_PROTOCOL.md section 3; main text sec:Real-Data-Examples"
'), stringsAsFactors = FALSE)
write.csv(prov, file.path(WP3, "parameter_provenance.csv"), row.names = FALSE)
cat(sprintf("  wrote %d rows to parameter_provenance.csv\n", nrow(prov)))

cat("\n=== done ===\n")
cat(sprintf("wrote %s\n", file.path(WP3, "wilcoxon_real.csv")))
cat(sprintf("wrote %s\n", file.path(WP3, "sim_variability_INVENTORY.md")))
if (wp2c_exists) cat(sprintf("wrote %s\n", file.path(WP3, "sim_variability.csv")))
cat(sprintf("wrote %s\n", file.path(WP3, "parameter_provenance.csv")))
