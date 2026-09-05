#!/usr/bin/env Rscript
# 18_wp2_assumption_checks.R -- the two WP2 checklist items that cannot be
# settled by reading the manuscript (Ceyhan, 2026-08-22).
#
# CHECK A -- empty outbound neighbourhood.
#   The manuscript states: "If N_O(x_i) = empty, ... by convention we set
#   OOS(x_i) := infinity so that such points receive the maximum outlyingness
#   rank." The code in Outlyingness_Score.R does NOT implement that convention:
#       out_nei.index = which(ddatax[x,] <= R[x])   # diag set to Inf
#       score = mean(Den[out_nei.index]) / Den[x]
#   with an empty index this is mean(numeric(0)) = NaN, not +Inf. The question
#   is therefore whether the empty case can arise at all. Both constructions
#   have a scores-mode fallback that should forbid it:
#     UN-CCD  nnccd.radi : if (scores && R_i == 0) R_i <- sort(ddx[i,])[2]
#     RK-CCD  see RK_CCD_New.R
#   This check counts, per dataset and per construction, how many points have a
#   covering ball containing no other point.
#
# CHECK B -- similarity equivariance of the radius selection.
#   Proposition 2.1 assumes the CCD construction is equivariant under
#   T(x) = cQx + b: radii scale by c and the neighbour relations and cluster
#   assignments are preserved. That is an assumption about the RADIUS SEARCH,
#   not about the score algebra, so it has to be measured. We apply a random
#   rotation, a translation and a scaling, rebuild, and compare radii (should be
#   exactly c times the originals, up to floating point), the adjacency matrix
#   M = [d(i,j) <= r_i] (should be identical), and the resulting scores.
#
# Runs on small datasets only -- these are structural properties, not
# performance numbers, so glass/vertebral/stamps suffice.
#
# CLI (PowerShell; Rscript segfaults under Bash):
#   Rscript "revision_experiments/tr2/18_wp2_assumption_checks.R" [dataset ...]

suppressMessages(library(here))
suppressMessages(source(here::here("revision_experiments", "shared", "harness.R")))

OUT_CSV <- here::here("revision_experiments/results/tr2/wp2_assumption_checks.csv")
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
SETS <- if (length(args) > 0) args else c("glass", "vertebral", "hepatitis", "stamps")

# Deterministic rotation + translation + scaling, no Math.random equivalent
# needed: a fixed seed keeps the check reproducible.
set.seed(20260822)

#' Radii and adjacency for one construction, in the point order the scorers use.
build <- function(kind, X, d) {
  if (kind == "RK") {
    tab <- get_simul("RK", d)
    g <- RKCCD_correct_quant(X, r.seq = 10, dom.method = "greedy2",
                             quan = tab$quant, simul = tab$simul,
                             niter = 1000, scores = TRUE, min.cls = 0)
  } else {
    tab <- get_simul("NN", d)
    g <- nnccd_clustering_quantile(X, low.num = 3, quantile = "lower",
                                   method = "ascend", dom.method = "greedy2",
                                   simul = tab$simul, niter = 1000, scores = TRUE)
  }
  R <- g$R[order(g$D)]
  D <- as.matrix(dist(X))
  list(R = R, M = (D <= R), label = g$label)
}

rows <- list()

for (ds in SETS) {
  dat <- load_real_dataset(ds)
  X <- as.matrix(dat$X); d <- dat$d; n <- dat$n
  cat(sprintf("\n===== %s: n=%d, d=%d =====\n", ds, n, d))

  for (kind in c("RK", "NN")) {
    b <- tryCatch(build(kind, X, d),
                  error = function(e) { cat(sprintf("  %s build ERROR: %s\n", kind,
                                                    conditionMessage(e))); NULL })
    if (is.null(b)) next

    # ---- CHECK A -------------------------------------------------------
    # Outbound count EXCLUDING self, exactly as OOS() does (diag -> Inf).
    Mo <- b$M; diag(Mo) <- FALSE
    out_deg <- rowSums(Mo)
    n_empty <- sum(out_deg == 0)
    n_zero_r <- sum(b$R <= 0)

    # Does the shipped OOS actually produce non-finite scores here?
    sc <- OOS(X, b$R, d)
    n_nan <- sum(is.nan(sc)); n_inf <- sum(is.infinite(sc))

    cat(sprintf("  [A] %s: min outbound degree %d, points with EMPTY outbound set %d, r<=0 %d | OOS NaN %d, Inf %d\n",
                kind, min(out_deg), n_empty, n_zero_r, n_nan, n_inf))

    # ---- CHECK B -------------------------------------------------------
    Q <- qr.Q(qr(matrix(rnorm(d * d), d, d)))   # random orthogonal
    cc <- 3.7                                    # scaling
    bb <- rnorm(d) * 5                           # translation
    Xt <- X %*% t(Q) * cc + matrix(bb, n, d, byrow = TRUE)

    bt <- tryCatch(build(kind, Xt, d),
                   error = function(e) { cat(sprintf("  %s transformed build ERROR: %s\n",
                                                     kind, conditionMessage(e))); NULL })
    if (is.null(bt)) next

    rel_r <- max(abs(bt$R - cc * b$R) / pmax(cc * b$R, .Machine$double.eps))
    adj_same <- identical(unname(b$M), unname(bt$M))
    n_adj_diff <- sum(b$M != bt$M)
    lab_same <- if (is.null(b$label) || is.null(bt$label)) NA
                else length(unique(paste(b$label, bt$label))) == length(unique(b$label))

    sc_t <- OOS(Xt, bt$R, d)
    finite_both <- is.finite(sc) & is.finite(sc_t)
    oos_rel <- if (any(finite_both))
      max(abs(sc_t[finite_both] - sc[finite_both]) /
            pmax(abs(sc[finite_both]), .Machine$double.eps)) else NA_real_

    cat(sprintf("  [B] %s: max rel radius error vs c*r = %.3e | adjacency identical: %s (%d cells differ) | clusters preserved: %s | max rel OOS change %.3e\n",
                kind, rel_r, adj_same, n_adj_diff, lab_same, oos_rel))

    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, n = n, d = d, construction = kind,
      min_outbound_degree = min(out_deg), n_empty_outbound = n_empty,
      n_radius_le_zero = n_zero_r, n_oos_nan = n_nan, n_oos_inf = n_inf,
      max_rel_radius_error = rel_r, adjacency_identical = adj_same,
      n_adjacency_cells_differ = n_adj_diff, clusters_preserved = lab_same,
      max_rel_oos_change = oos_rel, stringsAsFactors = FALSE)
  }
}

out <- do.call(rbind, rows)
write.csv(out, OUT_CSV, row.names = FALSE)

cat("\n==== verdicts ====\n")
cat(sprintf("  A. empty outbound set observed in %d of %d (dataset, construction) cells; total points affected %d\n",
            sum(out$n_empty_outbound > 0), nrow(out), sum(out$n_empty_outbound)))
cat(sprintf("     OOS returned NaN for %d points and +Inf for %d points overall\n",
            sum(out$n_oos_nan), sum(out$n_oos_inf)))
cat(sprintf("  B. radii scale by c to within %.2e; adjacency identical in %d of %d cells\n",
            max(out$max_rel_radius_error), sum(out$adjacency_identical), nrow(out)))
cat(sprintf("\nWrote %s\n18_wp2_assumption_checks.R done.\n", OUT_CSV))
