#!/usr/bin/env Rscript
# revision_experiments/tr3/35_paired_tests.R
#
# Paired significance tests between the CCD methods (Reviewer 1: the
# conclusions rest on observed means without formal statistical tests).
# Reads only; no clustering is run.
#
# Per-replicate results of all three CCD methods exist at n = 50: 30 cells per
# method (uniform / gaussian x d in 2, 3, 5, 10, 20 x 2, 3, 5 clusters), 500
# replicates each, on the same data sets across methods (variant "published").
# For each cell, UN-CCDs are paired with RK-CCDs and with KS-CCDs by replicate
# number, and the per-replicate ARI differences (UN minus other) go to a
# two-sided Wilcoxon signed-rank test,
#   wilcox.test(x, y, paired = TRUE, exact = FALSE, correct = TRUE).
# Zero differences are dropped by the test (n_nonzero counts the rest); a cell
# whose differences are all zero gets p = 1 and V = NA. The 60 p-values are
# adjusted together with Holm's method. Verdict at 0.05 on the adjusted p:
# "UN higher" / "UN lower" (sign of the mean difference) or "no difference".
# Cells where the signs of the mean and the median difference differ are listed
# on the console (a median of 0 against a nonzero mean counts as a difference).
#
# Second check: the per-cell mean ARIs of the CSV, at 3 decimals, against the
# n = 50 ARI values of the six simulation tables of Supplementary_Material.tex
# (tab:{uni,Gau}_{UN,RK,KS}-CCDs_Apdx; columns ARI/Sil/SR per d = 2, 3, 5, 10, 20;
# the first three rows of each table are n = 50 with 2, 3, 5 clusters).
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/35_paired_tests.R
#   TR3_OUT_DIR   output folder (default: results/tr3, i.e. TR3_RES)
#   SUPP_TEX      path of Supplementary_Material.tex
#                 (default G:/Submissions/TR3/TR3_10182026/Supplementary_Material.tex)
#
# INPUTS (revision_experiments/results/tr3/, read only)
#   sim_replicates.csv, sim_done.csv   variant "published", n = 50, the three
#                                      CCD methods; cells kept as listed in sim_done.csv
# OUTPUTS (TR3_OUT_DIR)
#   paired_tests.csv        one row per cell and comparison (60): dist, d, n,
#                           n_clusters, comparison, mean_UN, mean_other, mean_diff,
#                           median_diff, n_nonzero, V, p, p_holm, verdict
#   paired_tests_check.csv  one row per cell and method (90): cell, method,
#                           csv_mean_3dp, tex_value, match

source(here::here("revision_experiments/tr3/tr3_lib.R"))

OUT_DIR <- sub("[/\\\\]+$", "", Sys.getenv("TR3_OUT_DIR", TR3_RES))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SUPP_TEX <- Sys.getenv("SUPP_TEX", "G:/Submissions/TR3/TR3_10182026/Supplementary_Material.tex")

VARIANT <- "published"
N_FIX   <- "50"
IT_N    <- 500L
DISTS   <- c("uniform", "gaussian")
DS      <- c("2", "3", "5", "10", "20")
KS_     <- c("2", "3", "5")
ALPHA   <- 0.05

# ---------------------------------------------------------------------------
# 1. Rows of the 90 cells
# ---------------------------------------------------------------------------
keys <- c("variant", "method", "dist", "d", "n", "n_clusters")
rep_df <- tr3_read_done(tr3_path("sim_replicates.csv"), tr3_path("sim_done.csv"), keys, "rep")
if (is.null(rep_df)) stop("35_paired_tests.R: sim_replicates.csv or sim_done.csv is missing or empty")
rep_df <- rep_df[rep_df$variant == VARIANT & rep_df$n == N_FIX & rep_df$method %in% TR3_METHODS, , drop = FALSE]
rep_df$rep <- as.integer(rep_df$rep)
rep_df$ari <- as.numeric(rep_df$ari)
if (anyNA(rep_df$ari)) stop("35_paired_tests.R: NA in ari")
cat(sprintf("35_paired_tests.R: %d rows (variant %s, n = %s, %s); out = %s\n",
            nrow(rep_df), VARIANT, N_FIX, paste(TR3_METHODS, collapse = ", "), OUT_DIR))

get_cell <- function(method, dist, d, k) {
  r <- rep_df[rep_df$method == method & rep_df$dist == dist & rep_df$d == d & rep_df$n_clusters == k, , drop = FALSE]
  cellname <- sprintf("%s %s d=%s n=%s k=%s", method, dist, d, N_FIX, k)
  if (nrow(r) != IT_N) stop(sprintf("35_paired_tests.R: %s has %d replicates, not %d", cellname, nrow(r), IT_N))
  if (!identical(sort(r$rep), seq_len(IT_N))) stop("35_paired_tests.R: replicate numbers are not 1..", IT_N, " in ", cellname)
  r[order(r$rep), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 2. Paired Wilcoxon tests
# ---------------------------------------------------------------------------
res <- list()
means_csv <- list()   # for the check against the supplement
n_warn <- 0L
for (dist in DISTS) for (d in DS) for (k in KS_) {
  un <- get_cell("UN-CCDs", dist, d, k)
  oth <- list("RK-CCDs" = get_cell("RK-CCDs", dist, d, k), "KS-CCDs" = get_cell("KS-CCDs", dist, d, k))
  for (m in names(oth)) {
    if (!identical(un$data_fp, oth[[m]]$data_fp)) {
      bad <- which(un$data_fp != oth[[m]]$data_fp)
      stop(sprintf("35_paired_tests.R: data_fp differs between UN-CCDs and %s in %s d=%s k=%s (%d replicates, first rep %d)",
                   m, dist, d, k, length(bad), un$rep[bad[1]]))
    }
  }
  means_csv[[paste(dist, d, k, "UN-CCDs")]] <- mean(un$ari)
  for (m in names(oth)) {
    o <- oth[[m]]
    means_csv[[paste(dist, d, k, m)]] <- mean(o$ari)
    diff <- un$ari - o$ari
    nz <- sum(diff != 0)
    if (nz == 0) {
      V <- NA_real_; p <- 1
      cat(sprintf("  all differences are zero: %s d=%s k=%s UN vs %s (p = 1, V = NA)\n", dist, d, k, m))
    } else {
      w <- withCallingHandlers(
        wilcox.test(un$ari, o$ari, paired = TRUE, exact = FALSE, correct = TRUE),
        warning = function(w) { n_warn <<- n_warn + 1L; message("  warning: ", conditionMessage(w),
                                                                sprintf(" (%s d=%s k=%s UN vs %s)", dist, d, k, m)); invokeRestart("muffleWarning") })
      V <- unname(w$statistic); p <- w$p.value
    }
    res[[length(res) + 1L]] <- data.frame(
      dist = dist, d = d, n = N_FIX, n_clusters = k,
      comparison = paste0("UN vs ", sub("-CCDs", "", m)),
      mean_UN = mean(un$ari), mean_other = mean(o$ari), mean_diff = mean(diff),
      median_diff = median(diff), n_nonzero = nz, V = V, p = p,
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, res)
stopifnot(nrow(tab) == 60)
tab$p_holm <- p.adjust(tab$p, method = "holm")
tab$verdict <- ifelse(tab$p_holm < ALPHA & tab$mean_diff > 0, "UN higher",
               ifelse(tab$p_holm < ALPHA & tab$mean_diff < 0, "UN lower", "no difference"))
# p_holm < alpha with a mean difference of exactly 0 cannot be labelled by sign
if (any(tab$p_holm < ALPHA & tab$mean_diff == 0)) stop("35_paired_tests.R: significant cell with mean_diff = 0")

out <- data.frame(dist = tab$dist, d = tab$d, n = tab$n, n_clusters = tab$n_clusters,
                  comparison = tab$comparison,
                  mean_UN = fmt17(tab$mean_UN), mean_other = fmt17(tab$mean_other),
                  mean_diff = fmt17(tab$mean_diff), median_diff = fmt17(tab$median_diff),
                  n_nonzero = tab$n_nonzero, V = fmt17(tab$V), p = fmt17(tab$p),
                  p_holm = fmt17(tab$p_holm), verdict = tab$verdict, stringsAsFactors = FALSE)
utils::write.csv(out, file.path(OUT_DIR, "paired_tests.csv"), row.names = FALSE, quote = TRUE)
cat(sprintf("wrote %s (%d rows); %d wilcox.test warning(s)\n",
            file.path(OUT_DIR, "paired_tests.csv"), nrow(out), n_warn))

tab$d <- factor(tab$d, levels = DS)
cat("\nVerdicts per comparison (Holm over 60 tests, alpha = 0.05)\n")
print(table(comparison = tab$comparison, verdict = tab$verdict))
cat("\nVerdicts per comparison and d\n")
print(table(comparison = tab$comparison, d = tab$d, verdict = tab$verdict))
cat("\nVerdicts per comparison and distribution\n")
print(table(comparison = tab$comparison, dist = tab$dist, verdict = tab$verdict))
cat("\nVerdicts per comparison and number of clusters\n")
print(table(comparison = tab$comparison, k = tab$n_clusters, verdict = tab$verdict))

dis <- tab[sign(tab$mean_diff) != sign(tab$median_diff), , drop = FALSE]
opp <- sign(tab$mean_diff) * sign(tab$median_diff) < 0
cat(sprintf("\nCells where the sign of the mean difference differs from the sign of the median difference: %d of %d\n",
            nrow(dis), nrow(tab)))
cat(sprintf("  of these, median exactly 0 with a nonzero mean: %d; opposite nonzero signs: %d\n",
            sum(dis$median_diff == 0), sum(opp)))
cat(sprintf("  significant cells (Holm) among the %d: %d\n", nrow(dis), sum(dis$verdict != "no difference")))
if (nrow(dis)) {
  print(data.frame(dist = dis$dist, d = dis$d, k = dis$n_clusters, comparison = dis$comparison,
                   mean_diff = signif(dis$mean_diff, 4), median_diff = signif(dis$median_diff, 4),
                   n_nonzero = dis$n_nonzero, p_holm = signif(dis$p_holm, 3), verdict = dis$verdict),
        row.names = FALSE)
}
cat(sprintf("\nUnadjusted p < 0.05: %d of 60; Holm-adjusted p < 0.05: %d of 60\n",
            sum(tab$p < ALPHA), sum(tab$p_holm < ALPHA)))
cat("Verdict 'no difference' with raw p < 0.05 (lost to Holm):\n")
print(data.frame(tab[tab$p < ALPHA & tab$p_holm >= ALPHA, c("dist", "d", "n_clusters", "comparison", "mean_diff", "p", "p_holm")],
                 row.names = NULL))

# ---------------------------------------------------------------------------
# 3. Means of the CSV against the n = 50 ARI values of the supplement
# ---------------------------------------------------------------------------
strip_red <- function(x) {
  repeat {
    y <- gsub("\\\\textcolor\\{red\\}\\{([^{}]*)\\}", "\\1", x)
    if (identical(y, x)) break
    x <- y
  }
  x
}
tex <- readLines(SUPP_TEX, warn = FALSE)
TABLES <- expand.grid(dist = c("uniform", "gaussian"), method = TR3_METHODS, stringsAsFactors = FALSE)
chk <- list()
for (i in seq_len(nrow(TABLES))) {
  dist <- TABLES$dist[i]; m <- TABLES$method[i]
  lab <- sprintf("\\label{tab:%s_%s_Apdx}", if (dist == "uniform") "uni" else "Gau", m)
  li <- which(grepl(lab, tex, fixed = TRUE))
  if (length(li) != 1) stop("35_paired_tests.R: label ", lab, " found ", length(li), " times in ", SUPP_TEX)
  # the table's own n=50 row is the closest one above its caption line
  r50 <- tail(grep("n=50\\$", tex[seq_len(li)]), 1)
  # header: the method name must be in this table's header (checked in the 8 lines above the n=50 row)
  hdr <- strip_red(paste(tex[(r50 - 8):(r50 - 1)], collapse = " "))
  if (!grepl(m, hdr, fixed = TRUE)) stop("35_paired_tests.R: header of ", lab, " does not name ", m)
  rows <- tex[r50 + 0:2]
  for (j in 1:3) {
    f <- strsplit(strip_red(sub("\\\\\\\\.*$", "", rows[j])), "&", fixed = TRUE)[[1]]
    f <- trimws(f)
    if (length(f) != 17) stop(sprintf("35_paired_tests.R: %s row %d has %d fields, not 17", lab, j, length(f)))
    if (!grepl(sprintf("%s clusters", KS_[j]), f[2], fixed = TRUE)) stop("35_paired_tests.R: unexpected row label '", f[2], "' in ", lab)
    for (jd in seq_along(DS)) {
      tv <- f[3 + 3 * (jd - 1)]               # ARI, Sil, SR per d
      mv <- sprintf("%.3f", means_csv[[paste(dist, DS[jd], KS_[j], m)]])
      chk[[length(chk) + 1L]] <- data.frame(
        cell = sprintf("%s d=%s n=%s k=%s", dist, DS[jd], N_FIX, KS_[j]), method = m,
        csv_mean_3dp = mv, tex_value = tv, match = identical(mv, tv), stringsAsFactors = FALSE)
    }
  }
}
chk <- do.call(rbind, chk)
stopifnot(nrow(chk) == 90)
utils::write.csv(chk, file.path(OUT_DIR, "paired_tests_check.csv"), row.names = FALSE, quote = TRUE)
cat(sprintf("\nwrote %s (%d rows)\n", file.path(OUT_DIR, "paired_tests_check.csv"), nrow(chk)))
cat(sprintf("Means of the CSV against %s: %d of %d match at 3 decimals\n", basename(SUPP_TEX), sum(chk$match), nrow(chk)))
if (any(!chk$match)) {
  cat("Mismatches:\n")
  mm <- chk[!chk$match, ]
  mm$mean_17 <- vapply(seq_len(nrow(mm)), function(i) {
    p <- strsplit(sub(" n=50", "", mm$cell[i]), " ")[[1]]
    sprintf("%.6f", means_csv[[paste(p[1], sub("d=", "", p[2]), sub("k=", "", p[3]), mm$method[i])]])
  }, character(1))
  print(mm, row.names = FALSE)
}
