#!/usr/bin/env Rscript
# revision_experiments/tr1/wp0_mccd_methods.R
#
# WP0: wires the paper's four MCCD detectors (U-MCCD, SU-MCCD, UN-MCCD,
# SUN-MCCD) into the revision_experiments harness's METHOD_REGISTRY /
# REAL_DATA_THRESHOLDS, by APPENDING to those lists. harness.R itself is
# never modified (see CLAUDE.md hard rule); this file only reads it.
#
# PARAMETER FIX (this revision): the wrapper defaults below used to fall
# back silently to get_simul()'s rk_quant_for_d()/nn_quant_for_d() buckets
# and to min.cls=0 whenever no override was supplied. Those buckets were
# reverse-engineered from a DIFFERENT paper's RKCCD-OOS/UNCCD-OOS simulation
# drivers (see harness.R's comment above get_simul()) and are wrong for the
# RK-based methods here (and for SUN-MCCD's d>=10 case); min.cls=0 disables
# SU-MCCD/SUN-MCCD's shape-adaptive step entirely. See the "Paper-faithful
# parameter resolver" section below (added by this fix) for the real rule,
# sourced from the manuscript itself (main text L880, L605; Supplementary
# Material L655-677), and exposed as separately-callable, separately-named
# functions so a reviewer can audit each one against the cited line.
#
# Assumes revision_experiments/shared/harness.R has already been source()'d (for
# get_simul(), evaluate(), METHOD_REGISTRY, REAL_DATA_THRESHOLDS, and the
# `here` package already being attached). Source order matters:
#
#   source(here::here("revision_experiments/shared/harness.R"))
#   source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))
#
# ---------------------------------------------------------------------------
# Detector contract (verified by recon + a throwaway probe run against
# hepatitis/glass/ecoli before this file was written; see WP0 report):
#
#   RUMCCD_outlier (methods/outlier_detection/RU-MCCDs.R)   -> U-MCCD
#     (datax, simul=NULL, quant=0.99, cls=NULL, ind=NULL, niter=1000)
#   SUMCCD_outlier (methods/outlier_detection/SU-MCCDs.R)   -> SU-MCCD
#     (datax, simul=NULL, min.cls=0, low.num=2, quant=0.99)
#   UNMCCD_outlier (methods/outlier_detection/UN-MCCD.R)    -> UN-MCCD
#     (datax, simul, method="ascend", niter=1000)
#   SUNMCCD_outlier (methods/outlier_detection/SUN-MCCD.R)  -> SUN-MCCD
#     (datax, simul=NULL, min.cls=0, method="ascend", low.num=3)
#
# All four return c(clusters=list(...), label=list(...), radii=list(...))
# -- clusters[[i]] is a row-subset of datax for macro-cluster i, label[[i]]
# is the per-row mutual-catch-graph connected-component id WITHIN cluster i
# (same order as clusters[[i]]'s rows), and radii[[1]] is whatever the
# detector's internal RK/NN radius object was at return time (see the "radii"
# section of the WP0 report for which are real and which are NULL).
#
# NEW SURFACE AREA: all four source R/ccds/mKNN_CCD_functions.R (for
# connected.ksccd.m()), which is NOT sourced anywhere in harness.R today.
# It defines exactly two functions (connected.ksccd.m, connected.ksccd.m.old)
# that do not collide with anything harness.R already loads. Loading it here
# is therefore additive only.
# ---------------------------------------------------------------------------

if (!exists("RUMCCD_outlier"))  source(here::here("methods/outlier_detection/RU-MCCDs.R"))
if (!exists("SUMCCD_outlier"))  source(here::here("methods/outlier_detection/SU-MCCDs.R"))
if (!exists("UNMCCD_outlier"))  source(here::here("methods/outlier_detection/UN-MCCD.R"))
if (!exists("SUNMCCD_outlier")) source(here::here("methods/outlier_detection/SUN-MCCD.R"))

# ---------------------------------------------------------------------------
# Translation: detector's (clusters, label, radii) -> harness score interface
# ---------------------------------------------------------------------------
#
# Row identity is tracked by ROWNAME, not by raw column-1 value (the legacy
# R/general_functions/count.R postprocessor matches on raw values, which the
# task brief flags as unsafe under ties on real data -- confirmed unsafe in
# principle, not reused here).
#
# Outlier rule (matches count()'s semantics, reimplemented safely): within
# each macro-cluster's mutual-catch-graph component labeling, the LARGEST
# component (by row count; ties broken by the lowest component id, i.e.
# which.max()'s default tie rule, same as count()) is the "regular" majority;
# every row in a smaller (minority) component is scored 1 (outlier). A
# macro-cluster that is a single fully-connected component (table(label) has
# one level) contributes zero outliers, by construction.
#
# Two failure modes are tracked and MUST be reported, not silently absorbed:
#
#   1. unassigned_rows: rows never claimed by any clusters[[i]]. Defaulted to
#      score 0 (regular) per the task spec. Should be 0 in practice -- if
#      nonzero, that's a real finding about the detector's row coverage.
#
#   2. singleton_lost_rows: matrix row-subsetting with a *single* TRUE index
#      drops to a bare vector in R, and the returned vector's `names()` are
#      the DATASET'S COLUMN NAMES, not the row's rowname -- row identity is
#      unrecoverable from the object the detector actually returns. This is
#      a real, verified R semantics gap (confirmed with a direct probe: see
#      report), not a hypothetical. Since clusters[[i]] is produced by the
#      detector's OWN internal subsetting (which we are not allowed to
#      touch -- append-only rule), a singleton macro-cluster is structurally
#      unrecoverable here. Empirically, no singleton macro-cluster occurred
#      for hepatitis/glass/ecoli across any of the 4 methods (verified before
#      writing this file), so singleton_lost_rows was 0 in the smoke test,
#      but the code path is real and is exercised defensively (not dead
#      code) should a future dataset produce one. Rows lost this way are
#      folded into the same "unassigned" bucket (scored 0, counted, and
#      reported) since we cannot do anything more precise with them.
mccd_translate <- function(detector_result, n) {
  clusters <- detector_result$clusters
  labels   <- detector_result$label
  radii    <- detector_result$radii

  score        <- rep(0, n)
  cluster_id   <- rep(NA_integer_, n)
  connectivity <- rep(NA_integer_, n)
  claimed      <- rep(FALSE, n)
  singleton_lost_rows <- 0L

  for (i in seq_along(clusters)) {
    cl <- clusters[[i]]
    lb <- labels[[i]]

    if (is.null(dim(cl))) {
      # Row-subsetting inside the detector dropped this (size-1) cluster to
      # a bare vector; its rowname is gone. Cannot map back to an index.
      singleton_lost_rows <- singleton_lost_rows + length(lb)
      next
    }

    idx <- suppressWarnings(as.integer(rownames(cl)))
    if (length(idx) != length(lb) || any(is.na(idx))) {
      stop(sprintf("mccd_translate(): cluster %d rowname/label mismatch (length %d vs %d, any NA idx: %s) -- rowname propagation assumption violated",
                   i, length(idx), length(lb), any(is.na(idx))))
    }

    tab <- table(lb)
    majority <- as.integer(names(tab)[which.max(tab)])
    score[idx]        <- ifelse(lb != majority, 1, 0)
    cluster_id[idx]   <- i
    connectivity[idx] <- lb
    claimed[idx]      <- TRUE
  }

  unassigned_rows <- sum(!claimed)  # includes singleton_lost_rows, by construction

  list(score = score, cluster = cluster_id, connectivity = connectivity,
       radii = radii, unassigned_rows = unassigned_rows,
       singleton_lost_rows = singleton_lost_rows)
}

# ---------------------------------------------------------------------------
# Paper-faithful parameter resolver (WP0 fix)
# ---------------------------------------------------------------------------
#
# Do NOT reuse harness.R's rk_quant_for_d()/nn_quant_for_d() buckets for the
# RK-based methods (U-MCCD, SU-MCCD) or for SUN-MCCD's d>=10 case -- those
# buckets were written for a different paper's methods (see harness.R's own
# comment immediately above get_simul()). This paper's alpha schedule is
# stated explicitly in:
#   - Main text, line 880 (RK-based methods, U-MCCD & SU-MCCD)
#   - Main text, line 605 (S_min "set to half the contamination level")
#   - Supplementary Material, lines 655-677 (RK schedule restated; the full
#     NND alpha(d) tables for UN-MCCD and SUN-MCCD, including the one point
#     where they diverge)
#
# get_simul()'s own `quant` argument is a FILE-LABEL STRING ("99", "999",
# "95", ...) matching the RData filename, not a bare probability. The
# resolvers below return that label, so callers do
#   get_simul(variant, d, quant = <label>)
# and are guaranteed the loaded $simul table and the numeric quantile handed
# to the RK-CCD internals come from the SAME file. This is not a style
# preference: RUMCCD_outlier/SUMCCD_outlier pass `quant` straight through to
# ccd.Kest.edge.quantile(), which indexes `simul$quan[[as.character(quant)]]`
# -- and each RK-test-simul_*d_*%.RData file's $quan sub-list has EXACTLY
# ONE key, matching only its own filename (verified directly: the 9d_99%
# file's $quan has key "0.99" only; the 9d_999% file has "0.999" only).
# Loading one file's table while telling the RK internals "this is quantile
# X" for a different X returns NULL from that index, `any(NULL > y)` is
# FALSE, and every point silently gets radius 0 -- a real failure mode this
# resolver design rules out by construction, not a hypothetical. The NND
# detectors (UNMCCD_outlier/SUNMCCD_outlier) have no analogous risk: they
# take no separate `quant` scalar at all, so the quantile is controlled
# entirely by which file's $average/$median vectors get loaded.

#' RK-based methods (U-MCCD, SU-MCCD): main text line 880 / SM line 656.
#' alpha = 1% for d<10, alpha = 0.1% for d>=10  =>  quant label "99"/"999".
rk_quant_label_paper <- function(d) if (d < 10) "99" else "999"

#' NND-based UN-MCCD: SM lines 657-667 (alpha = 15%,10%,5%,1%,0.1% at
#' d=2,3,5,10,{20,50,100}). This coincides with harness.R's nn_quant_for_d()
#' step rule -- confirmed correct for this paper, unlike the RK bucket -- so
#' it is kept as an explicit, separately-named pass-through (not
#' reimplemented independently) so the two definitions cannot silently
#' drift apart.
nn_quant_label_paper_UN <- function(d) nn_quant_for_d(d)

#' NND-based SUN-MCCD: SM lines 668-677. Identical to UN-MCCD for d<10;
#' forces alpha=0.1% ("999") already at d=10, one step earlier than
#' UN-MCCD (which is still at alpha=1%, "99", at d=10..19).
nn_quant_label_paper_SUN <- function(d) if (d < 10) nn_quant_for_d(d) else "999"

#' S_min / min.cls ("set to half the contamination level", main text line
#' 605) is genuinely ambiguous about which half-of-what to round and which
#' n0/contamination source to trust (see task brief re: Table 5's published
#' outlier counts being wrong for glass and ecoli). This function returns
#' ALL defensible readings, deduplicated by resulting integer value -- it
#' deliberately does not pick a winner; the caller (13_wp0_gate.R) runs and
#' reports every distinct value.
#'
#' @param n0 loader's actual outlier count for the dataset (sum(Y==0)).
#' @param n loader's actual sample size.
#' @param contamination_pct_table5 the manuscript's published Table 5
#'   contamination fraction (e.g. 0.095 for hepatitis) -- may disagree with
#'   n0/n for datasets where Table 5's own outlier count is wrong.
#' @return named list: integer min.cls value (as character key) -> character
#'   vector of reading labels ("a","b","c","d") that produced it.
mccd_min_cls_readings <- function(n0, n, contamination_pct_table5) {
  readings <- c(
    a = floor(n0 / 2),                                  # floor(n0/2)
    b = ceiling(n0 / 2),                                 # ceiling(n0/2)
    c = round(0.5 * (n0 / n) * n),                        # 0.5 * contamination(loader) * n, as stated
    d = floor(0.5 * contamination_pct_table5 * n)         # 0.5 * contamination(Table 5) * n
  )
  split(names(readings), unname(readings))
}

# ---------------------------------------------------------------------------
# UNITS FIX (WP0 gate re-run): mccd_min_cls_readings() above returns RAW
# COUNTS (e.g. floor(n0/2) = 15 for vertebral). Those counts were being
# handed straight to the `min.cls` argument of SUMCCD_outlier/SUNMCCD_outlier.
# Confirmed directly against source: R/ccds/RK_CCD_New.R line 113
# (rccd.silhouette, used by SU-MCCD via RKCCD_correct_quant) and
# R/ccds/UN_CCD.R line 102 (nnccd.silhouette, used by SUN-MCCD) both compute
#   lenD = length(which(graph$catch >= round(min.cls * n)))
# i.e. `min.cls` IS A PROPORTION OF n (the line-110/-101 comment above each
# says so explicitly: "the minimum percentage accepted as a cluster"), not a
# count. Handing it the raw count 15 against n=240 makes the threshold
# round(15*240)=3600 -- unreachable by any cluster's catch count, so lenD
# collapses to 0 and the shape-adaptive step is disabled. This is a real,
# verified unit mismatch, not a hypothesis.
#
# This function is intentionally NOT a replacement for
# mccd_min_cls_readings() above (left untouched, as the record of what the
# pre-fix wiring computed) -- it is the corrected, unit-consistent resolver,
# separately named so the two are never confused. `n` here must be the
# WHOLE-DATASET n (verified: SUMCCD_outlier/SUNMCCD_outlier build ddatax /
# call nnccd.silhouette on the full `datax` argument, not a per-cluster
# subset, so round(min.cls*n) is always relative to the whole dataset).
#
# Two readings only, per the task brief (main text L605, S_min "set to half
# the contamination level"), both expressed as proportions of n so that the
# detector's own round(min.cls*n) recovers the intended COUNT:
#   half_contam = 0.5 * (n0/n)   -- the literal reading of L605
#   full_contam = n0/n           -- comparison reading (full, not half,
#                                    contamination); reported alongside, not
#                                    chosen as a winner.
#
# @param n0 loader's actual outlier count (sum(Y==0)).
# @param n  loader's actual sample size.
# @return named list: half_contam=<proportion>, full_contam=<proportion>.
mccd_min_cls_proportion_readings <- function(n0, n) {
  list(half_contam = 0.5 * (n0 / n), full_contam = n0 / n)
}

# ---------------------------------------------------------------------------
# Method wrappers: function(X, d, Y = NULL, ...) -> list(score=, t_construct=,
# t_total=, ...). Timing note: all four detectors are single monolithic calls
# (digraph construction + mutual-catch-graph component search happen inside
# one function body with no exposed sub-timing hook); splitting construction
# from total would require a second, duplicate, redundant call into the
# detector's internals (as e.g. rkccd_oos_method() does for RKCCD-OOS by
# calling RKCCD_correct_quant() once for t_construct and RKCCD_OOS() again
# for t_total). We deliberately do NOT do that here, to avoid doubling
# reported runtime for a distinction that would be synthetic anyway (there's
# no natural "scoring-only" phase for these clustering-based detectors, only
# clustering itself). Per the task's explicit fallback: t_construct and
# t_total are the SAME single measured wall-clock time for the one call.
# ---------------------------------------------------------------------------

# NOTE on `quant` arguments below: unlike the pre-fix code, `quant` here is
# the FILE-LABEL STRING get_simul() expects ("99"/"999"/"95"/...), NOT a
# bare numeric probability -- this is what makes it safe to override without
# desynchronizing the loaded table from the value handed to the detector
# (see the resolver section above). NULL (the default) resolves via the
# paper-faithful rule for that method; this is the actual bug fix. `min.cls`
# remains an explicit, un-resolved argument (default 0 = shape-adaptive step
# disabled, a legitimate reading in its own right, not a stand-in for "the
# paper's value") -- callers must pass a value from mccd_min_cls_readings()
# explicitly, by design (see that function's docstring: the ambiguity is not
# ours to resolve silently).

umccd_method <- function(X, d, Y = NULL, quant = NULL, ...) {
  X <- as.matrix(X)
  rownames(X) <- as.character(seq_len(nrow(X)))
  quant_source <- if (is.null(quant)) "paper default (main text L880): d<10 -> 99%, d>=10 -> 999%" else "override"
  q_label <- if (is.null(quant)) rk_quant_label_paper(d) else quant
  tab <- get_simul("RK", d, quant = q_label)
  t0 <- Sys.time()
  res <- RUMCCD_outlier(datax = X, simul = tab$simul, quant = tab$quant)
  t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- mccd_translate(res, nrow(X))
  list(score = tr$score, t_construct = t, t_total = t,
       cluster = tr$cluster, connectivity = tr$connectivity, radii = tr$radii,
       unassigned_rows = tr$unassigned_rows, singleton_lost_rows = tr$singleton_lost_rows,
       quant_used = tab$quant, quant_label = tab$quant_label, quant_source = quant_source,
       min_cls_used = NA_integer_)
}

sumccd_method <- function(X, d, Y = NULL, quant = NULL, min.cls = 0, low.num = 2, ...) {
  X <- as.matrix(X)
  rownames(X) <- as.character(seq_len(nrow(X)))
  quant_source <- if (is.null(quant)) "paper default (main text L880): d<10 -> 99%, d>=10 -> 999%" else "override"
  q_label <- if (is.null(quant)) rk_quant_label_paper(d) else quant
  tab <- get_simul("RK", d, quant = q_label)
  t0 <- Sys.time()
  res <- SUMCCD_outlier(datax = X, simul = tab$simul, min.cls = min.cls, low.num = low.num, quant = tab$quant)
  t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- mccd_translate(res, nrow(X))
  list(score = tr$score, t_construct = t, t_total = t,
       cluster = tr$cluster, connectivity = tr$connectivity, radii = tr$radii,
       unassigned_rows = tr$unassigned_rows, singleton_lost_rows = tr$singleton_lost_rows,
       quant_used = tab$quant, quant_label = tab$quant_label, quant_source = quant_source,
       min_cls_used = min.cls)
}

unmccd_method <- function(X, d, Y = NULL, method = "ascend", quant = NULL, ...) {
  X <- as.matrix(X)
  rownames(X) <- as.character(seq_len(nrow(X)))
  quant_source <- if (is.null(quant)) "paper default (SM L657-667; == harness nn_quant_for_d(), unchanged)" else "override"
  q_label <- if (is.null(quant)) nn_quant_label_paper_UN(d) else quant
  tab <- get_simul("NN", d, quant = q_label)
  t0 <- Sys.time()
  res <- UNMCCD_outlier(datax = X, simul = tab$simul, method = method)
  t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- mccd_translate(res, nrow(X))
  list(score = tr$score, t_construct = t, t_total = t,
       cluster = tr$cluster, connectivity = tr$connectivity, radii = tr$radii,
       unassigned_rows = tr$unassigned_rows, singleton_lost_rows = tr$singleton_lost_rows,
       quant_used = tab$quant, quant_label = tab$quant_label, quant_source = quant_source,
       min_cls_used = NA_integer_)
}

sunmccd_method <- function(X, d, Y = NULL, method = "ascend", min.cls = 0, low.num = 3, quant = NULL, ...) {
  X <- as.matrix(X)
  rownames(X) <- as.character(seq_len(nrow(X)))
  quant_source <- if (is.null(quant)) "paper default (SM L668-677): forces 999% at d>=10, unlike UN-MCCD" else "override"
  q_label <- if (is.null(quant)) nn_quant_label_paper_SUN(d) else quant
  tab <- get_simul("NN", d, quant = q_label)
  t0 <- Sys.time()
  res <- SUNMCCD_outlier(datax = X, simul = tab$simul, min.cls = min.cls, method = method, low.num = low.num)
  t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- mccd_translate(res, nrow(X))
  list(score = tr$score, t_construct = t, t_total = t,
       cluster = tr$cluster, connectivity = tr$connectivity, radii = tr$radii,
       unassigned_rows = tr$unassigned_rows, singleton_lost_rows = tr$singleton_lost_rows,
       quant_used = tab$quant, quant_label = tab$quant_label, quant_source = quant_source,
       min_cls_used = min.cls)
}

# ---------------------------------------------------------------------------
# Registry append (never reassign; double-sourcing this file must not
# duplicate entries -- guarded by name presence, not by a sourced-once flag,
# so it is correct even if harness.R's METHOD_REGISTRY were rebuilt between
# two sources of this file).
# ---------------------------------------------------------------------------

WP0_METHODS <- list(
  "U-MCCD"   = umccd_method,
  "SU-MCCD"  = sumccd_method,
  "UN-MCCD"  = unmccd_method,
  "SUN-MCCD" = sunmccd_method
)
WP0_THRESHOLDS <- list(
  "U-MCCD" = 0.5, "SU-MCCD" = 0.5, "UN-MCCD" = 0.5, "SUN-MCCD" = 0.5
)

new_method_names <- setdiff(names(WP0_METHODS), names(METHOD_REGISTRY))
if (length(new_method_names) > 0) {
  METHOD_REGISTRY <- c(METHOD_REGISTRY, WP0_METHODS[new_method_names])
}

new_threshold_names <- setdiff(names(WP0_THRESHOLDS), names(REAL_DATA_THRESHOLDS))
if (length(new_threshold_names) > 0) {
  REAL_DATA_THRESHOLDS <- c(REAL_DATA_THRESHOLDS, WP0_THRESHOLDS[new_threshold_names])
}

cat(sprintf("wp0_mccd_methods.R: METHOD_REGISTRY now has %d entries (%s); REAL_DATA_THRESHOLDS has %d.\n",
            length(METHOD_REGISTRY), paste(names(METHOD_REGISTRY), collapse = ", "), length(REAL_DATA_THRESHOLDS)))
