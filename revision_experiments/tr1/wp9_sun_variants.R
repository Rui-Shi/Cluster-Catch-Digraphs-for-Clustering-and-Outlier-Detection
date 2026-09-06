#!/usr/bin/env Rscript
# revision_experiments/tr1/wp9_sun_variants.R
#
# WP9 (R3.10a): SUN-MCCD component ablation -- three toggles (NND statistic,
# centre-point removal, radius-search direction). See WP9_PROTOCOL.md for the
# full design, the dropped-Holm note, and the centre-off indexing decision;
# this file is the implementation those decisions describe.
#
# READ-ONLY SOURCES THIS FILE COPIES FROM (never edited):
#   R/ccds/UN_CCD.R            -- nnccd.radi() (L233-323), copied and
#                                  parameterised as nnccd.radi.ablate() below.
#   methods/outlier_detection/SUN-MCCD.R -- SUNMCCD_outlier(), substituted via
#                                  environment surgery, never redefined.
#
# SUBSTITUTION RULE (WP9_PROTOCOL.md "Substitution mechanism"): the ablated
# nnccd.radi is never assigned into .GlobalEnv. `nnccd_clustering_quantile`
# and `SUNMCCD_outlier` both resolve their internal calls lexically through
# their own closure environment (.GlobalEnv, since both were created by a
# top-level source()); run_sunmccd_ablation() below makes COPIES of those two
# functions whose `environment()` points into new, private child environments
# that carry just the one overridden binding each, leaving every other symbol
# lookup (nnccd.silhouette, connected.ksccd.m, ksccd.connected, dist, ...) to
# fall through to .GlobalEnv exactly as it always did. Nothing in .GlobalEnv
# is ever touched, so re-sourcing SUN-MCCD.R / UN_CCD.R (as
# wp0_mccd_methods.R's `if (!exists(...))` guards do) cannot revert it --
# there is nothing to revert.
#
# Assumes revision_experiments/shared/harness.R and
# revision_experiments/tr1/wp0_mccd_methods.R have already been sourced (for
# SUNMCCD_outlier, nnccd_clustering_quantile, nnccd.silhouette,
# connected.ksccd.m, ksccd.connected, get_simul, nn_quant_label_paper_SUN,
# evaluate, mccd_translate, METHOD_REGISTRY, REAL_DATA_THRESHOLDS).
# ---------------------------------------------------------------------------

stopifnot(
  "wp9_sun_variants.R: source harness.R and wp0_mccd_methods.R first" =
    all(sapply(c("nnccd_clustering_quantile", "SUNMCCD_outlier", "nnccd.silhouette",
                 "get_simul", "nn_quant_label_paper_SUN", "evaluate", "mccd_translate",
                 "METHOD_REGISTRY", "REAL_DATA_THRESHOLDS"), exists))
)

# ---------------------------------------------------------------------------
# 1. nnccd.radi.ablate() -- parameterised copy of nnccd.radi()
#    (R/ccds/UN_CCD.R:233-323). Copies the "quantile == 'lower'" branch only,
#    matching the original: every call site in this codebase passes
#    quantile = "lower" and the original function returns an all-zero R for
#    any other value (no other branch is implemented there either), so
#    reproducing exactly that scope is faithful, not a simplification.
#
# NEW ARGUMENTS (the three toggles):
#   stat           "both" (current: reject CSR if EITHER mean or median NND
#                  falls below its own lower quantile) / "mean" / "median"
#                  (reject using only that one statistic).
#   remove_centre  TRUE (current: the ball's own centre point i is excluded
#                  from the observed NND statistic, o.d[2:j] / o.d[j:(n-1)])
#                  / FALSE (centre included, o.d[1:j] / o.d[j:n]).
#
# CENTRE-OFF ENVELOPE INDEX (WP9_PROTOCOL.md "centre-off indexing decision").
# The envelope's convention is "entry k models the statistic of a k-point
# CSR cloud". The on-variant's envelope index already equals its observed
# cloud size (ascend: index j-1 for a (j-1)-point cloud; descend: index
# `size-1` for a `size`-point cloud, per that branch's own pre-existing
# offset). Including the centre adds exactly one point to the observed
# cloud, so the off-variant reads the envelope ONE ENTRY HIGHER than the
# on-variant would at the same step -- ascend: index j instead of j-1;
# descend: rev(...)[j+1] instead of rev(...)[j+2] (rev()'s subscript
# DECREASES by 1 to move the underlying envelope index UP by 1). This is
# implemented via a single `centre_shift <- if (remove_centre) 0L else 1L`
# added to the on-variant's index arithmetic, so the two directions share
# one rule instead of two independently-decided ones.
# ---------------------------------------------------------------------------

nnccd.radi.ablate <- function(dx, quantile = "lower", method = c("ascend", "descend"),
                               low.num, quant, simul = NULL, niter, scores = FALSE,
                               stat = c("both", "mean", "median"), remove_centre = TRUE) {
  method <- match.arg(method)
  stat   <- match.arg(stat)

  ddx <- as.matrix(dist(dx))
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0, n)

  centre_shift <- if (isTRUE(remove_centre)) 0L else 1L

  # rejection (ascend) / acceptance (descend) predicates for the requested stat
  reject_ascend <- function(avg_obs, med_obs, avg_bound, med_bound) {
    switch(stat,
           both   = (avg_obs < avg_bound) | (med_obs < med_bound),
           mean   = (avg_obs < avg_bound),
           median = (med_obs < med_bound))
  }
  accept_descend <- function(avg_obs, med_obs, avg_bound, med_bound) {
    switch(stat,
           both   = (avg_obs > avg_bound) & (med_obs > med_bound),
           mean   = (avg_obs > avg_bound),
           median = (med_obs > med_bound))
  }

  if (quantile == "lower") {
    if (!is.null(simul)) {
      NN.envelop <- list(average = simul$average[1:n], median = simul$median[1:n])
    } else {
      NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)
    }

    for (i in 1:n) {
      if (method == "ascend") {
        o.d <- order(ddx[i, ])
        for (j in low.num:n) {
          r <- ddx[i, o.d[j]]
          cloud_idx <- if (remove_centre) o.d[2:j] else o.d[1:j]
          NN.dist.obs <- NNDest.dist.f(ddx[cloud_idx, cloud_idx], r)

          env_idx <- (j - 1L) + centre_shift
          lower.bound.ave <- NN.envelop$average[env_idx]
          lower.bound.med <- NN.envelop$median[env_idx]

          if (reject_ascend(NN.dist.obs$averge, NN.dist.obs$median, lower.bound.ave, lower.bound.med)) {
            if (j == low.num) R[i] <- 0 else R[i] <- ddx[i, o.d[j - 1]]
            break
          }
        }
      }
      if (method == "descend") {
        o.d <- order(ddx[i, ], decreasing = TRUE)
        for (j in 1:(n - low.num)) {
          r <- ddx[i, o.d[j]]
          cloud_idx <- if (remove_centre) o.d[j:(n - 1)] else o.d[j:n]
          NN.dist.obs <- NNDest.dist.f(ddx[cloud_idx, cloud_idx], r)

          # rev(v)[m] = v[n - m + 1]; decreasing m by centre_shift moves the
          # underlying envelope index UP by centre_shift (see header note).
          rev_sub <- (j + 2L) - centre_shift
          lower.bound.ave <- rev(NN.envelop$average)[rev_sub]
          lower.bound.med <- rev(NN.envelop$median)[rev_sub]

          if (accept_descend(NN.dist.obs$averge, NN.dist.obs$median, lower.bound.ave, lower.bound.med)) {
            R[i] <- r
            break
          }
        }
      }
      if (scores && R[i] == 0) R[i] <- sort(ddx[i, ])[2]  # avoid 0 radius, matches original scores branch
    }
  }
  list(R = R, KS = NULL)
}

# ---------------------------------------------------------------------------
# 2. Explicit, non-global substitution -- see file header and
#    WP9_PROTOCOL.md "Substitution mechanism".
# ---------------------------------------------------------------------------

#' Build a `nnccd_clustering_quantile` whose internal `nnccd.radi` call
#' resolves to `nnccd.radi.ablate(..., stat=, remove_centre=)` instead of the
#' stock `nnccd.radi`, WITHOUT touching .GlobalEnv.
#'
#' @return a function with the same signature/semantics as
#'   nnccd_clustering_quantile, a fresh copy with a rebound environment.
make_ablated_clustering_fn <- function(stat, remove_centre) {
  nnradi_override <- function(dx, quantile = "lower", method = "ascend", low.num, quant,
                               simul = NULL, niter, scores = FALSE) {
    nnccd.radi.ablate(dx = dx, quantile = quantile, method = method, low.num = low.num,
                       quant = quant, simul = simul, niter = niter, scores = scores,
                       stat = stat, remove_centre = remove_centre)
  }

  env_radi <- new.env(parent = environment(nnccd_clustering_quantile))
  env_radi$nnccd.radi <- nnradi_override
  env_radi$nnccd.radi.ablate <- nnccd.radi.ablate  # in case of nested lookups

  nccq_copy <- nnccd_clustering_quantile
  environment(nccq_copy) <- env_radi
  nccq_copy
}

#' Run the SUN-MCCD pipeline with the ablated radius search substituted in,
#' via a second layer of environment surgery around SUNMCCD_outlier (which
#' itself calls nnccd_clustering_quantile unqualified).
#'
#' @param datax data matrix/data.frame passed straight to SUNMCCD_outlier.
#' @param simul the NN quantile table (get_simul("NN", d, ...)$simul).
#' @param min.cls S_min proportion (see CLAUDE.md -- a proportion of n).
#' @param method "ascend" or "descend" -- the direction toggle.
#' @param low.num SUN-MCCD's own low.num (default 3, per SUN-MCCD.R).
#' @param stat "both"/"mean"/"median" -- the NND-statistic toggle.
#' @param remove_centre TRUE/FALSE -- the centre-removal toggle.
#' @return the same list SUNMCCD_outlier returns: clusters, label, radii.
run_sunmccd_ablation <- function(datax, simul = NULL, min.cls = 0, method = "ascend",
                                  low.num = 3, stat = "both", remove_centre = TRUE) {
  nccq_copy <- make_ablated_clustering_fn(stat = stat, remove_centre = remove_centre)

  env_sun <- new.env(parent = environment(SUNMCCD_outlier))
  env_sun$nnccd_clustering_quantile <- nccq_copy

  sun_copy <- SUNMCCD_outlier
  environment(sun_copy) <- env_sun

  sun_copy(datax = datax, simul = simul, min.cls = min.cls, method = method, low.num = low.num)
}

# ---------------------------------------------------------------------------
# 3. Harness-compatible method wrapper, mirroring sunmccd_method()
#    (wp0_mccd_methods.R) exactly except for the call into
#    run_sunmccd_ablation() instead of straight SUNMCCD_outlier(), and the
#    three extra toggle arguments.
# ---------------------------------------------------------------------------

sunmccd_ablate_method <- function(X, d, Y = NULL, method = "ascend", min.cls = 0,
                                   low.num = 3, quant = NULL,
                                   stat = "both", remove_centre = TRUE, ...) {
  X <- as.matrix(X)
  rownames(X) <- as.character(seq_len(nrow(X)))
  q_label <- if (is.null(quant)) nn_quant_label_paper_SUN(d) else quant
  tab <- get_simul("NN", d, quant = q_label, n = nrow(X))
  t0 <- Sys.time()
  res <- run_sunmccd_ablation(datax = X, simul = tab$simul, min.cls = min.cls,
                               method = method, low.num = low.num,
                               stat = stat, remove_centre = remove_centre)
  t <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  tr <- mccd_translate(res, nrow(X))
  list(score = tr$score, t_construct = t, t_total = t,
       cluster = tr$cluster, connectivity = tr$connectivity, radii = tr$radii,
       unassigned_rows = tr$unassigned_rows, singleton_lost_rows = tr$singleton_lost_rows,
       quant_used = tab$quant, quant_label = tab$quant_label,
       min_cls_used = min.cls, stat = stat, remove_centre = remove_centre, method_used = method)
}

# ---------------------------------------------------------------------------
# 4. The 5-variant grid (one-at-a-time ablation; see WP9_PROTOCOL.md)
# ---------------------------------------------------------------------------

WP9_VARIANTS <- list(
  stock              = list(stat = "both",   remove_centre = TRUE,  method = "ascend"),
  stat_mean          = list(stat = "mean",   remove_centre = TRUE,  method = "ascend"),
  stat_median        = list(stat = "median", remove_centre = TRUE,  method = "ascend"),
  centre_off         = list(stat = "both",   remove_centre = FALSE, method = "ascend"),
  direction_descend  = list(stat = "both",   remove_centre = TRUE,  method = "descend")
)

# ---------------------------------------------------------------------------
# 5. Synthetic generators -- verbatim copy of 55_wp2c_simulation_arm.R's
#    gen_uniform()/gen_gaussian() (themselves lifted from
#    09_wp3_synthetic.R:124-195, tracing to the original
#    RKCCD_OOS_IOS/Simulation/{Uniform,Gaussian}/10d/10d_2cls_n500_cont5%.R
#    drivers). n, d, cont are arguments; every other constant is unchanged.
# ---------------------------------------------------------------------------

gen_uniform <- function(seed, n, d, cont) {
  cls_dis = 3; otl_dis = 2; r_min = 0.7; r_max = 1.3
  mu1 = rep(3, d); mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = apply(rbind(mu1, mu2), 2, mean)
  n1 = round(n * (1 - cont) * 0.5); n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)
  set.seed(seed)
  data1 = rpoisball.unit(n1, d) * runif(1, r_min, r_max) + matrix(rep(mu1, n1), ncol = d, byrow = TRUE)
  data2 = rpoisball.unit(n2, d) * runif(1, r_min, r_max) + matrix(rep(mu2, n2), ncol = d, byrow = TRUE)
  i = 0; outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier = rbind(outlier, temp); i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

gen_gaussian <- function(seed, n, d, cont) {
  cls_dis = 3; otl_dis = 2; r_min = 0.7; r_max = 1.3
  mu1 = rep(3, d); mu2 = c(3 + cls_dis, rep(3, d - 1))
  mu = apply(rbind(mu1, mu2), 2, mean)
  n1 = round(n * (1 - cont) * 0.5); n2 = round(n * (1 - cont) * 0.5) - 1
  n0 = round(n * cont)
  noise_level = 0.01
  sigma = 1 / sqrt(qchisq(1 - noise_level, d))
  set.seed(seed)
  data1 = mvrnorm(n1, mu1, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  data2 = mvrnorm(n2, mu2, diag(d) * (sigma * runif(1, r_min, r_max))^2)
  i = 0; outlier = NULL
  while (i < n0) {
    temp = rpoisball.unit(1, d) * 5 + mu
    if (sqrt(sum((temp - mu1)^2)) > otl_dis & sqrt(sum((temp - mu2)^2)) > otl_dis) {
      outlier = rbind(outlier, temp); i = i + 1
    }
  }
  rownames(outlier) = NULL
  list(X = rbind(data1, data2, outlier), n = n1 + n2 + n0, n0 = n0)
}

WP9_SETTINGS <- do.call(rbind, lapply(c("uniform", "gaussian"), function(gen)
  do.call(rbind, lapply(c(3L, 10L), function(dd)
    data.frame(generator = gen, d = dd, n_nominal = 200L, contam = 0.05,
               stringsAsFactors = FALSE)))))
WP9_SETTINGS$setting_id <- sprintf("%s_d%d", WP9_SETTINGS$generator, WP9_SETTINGS$d)

WP9_BASE_SEED <- 9100L

cat("wp9_sun_variants.R: nnccd.radi.ablate, run_sunmccd_ablation, sunmccd_ablate_method, ",
    "WP9_VARIANTS (", length(WP9_VARIANTS), "), WP9_SETTINGS (", nrow(WP9_SETTINGS), " rows) ready.\n", sep = "")
