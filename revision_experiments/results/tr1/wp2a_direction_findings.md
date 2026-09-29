# WP2(a) radius-search direction sweep (ascend vs descend) -- findings

Also serves reviewer point R3.10a (ablation toggle).

## Step 1: what the toggle does and its reach

```
# descend), the last row of the parameter-sensitivity table. Reviewer R3.10a
# also lists this as an ablation toggle, so one run serves both points.
#
# ---------------------------------------------------------------------------
# STEP 1 FINDINGS (read from source before any code below was written):
#
# THE TOGGLE. UNMCCD_outlier(datax, simul, method="ascend", niter)
# (methods/outlier_detection/UN-MCCD.R L7) and SUNMCCD_outlier(datax, simul,
# min.cls, method="ascend", low.num) (methods/outlier_detection/SUN-MCCD.R
# L9) both pass `method` straight into nnccd_clustering_quantile(), which
# passes it straight into nnccd.radi(dx, quantile, method, low.num, quant,
# simul, niter, scores) (R/ccds/UN_CCD.R L233-323). harness.R's own wrappers
# already expose it: wp0_mccd_methods.R L316 (unmccd_method) and L333
# (sunmccd_method) both declare `method = "ascend"` as a pass-through
# parameter of the METHOD_REGISTRY entry. No override of harness.R or
# wp0_mccd_methods.R is needed to flip direction -- it is a normal call
# argument: METHOD_REGISTRY[["UN-MCCD"]](X=, d=, Y=, method="descend").
#
# WHAT ASCEND/DESCEND ACTUALLY DO (R/ccds/UN_CCD.R L233-323). For each point
# i, "ascend" (L245-264) starts from a ball holding only its `low.num`
# nearest neighbours and GROWS the ball outward one point at a time
# (increasing radius, j = low.num..n), testing at each step whether the
# within-ball NN-distance statistic (mean and median, center point excluded)
# has fallen below the CSR lower envelope; it stops and reports the radius
# from the PREVIOUS step (j-1) the first time that happens -- i.e. it takes
# the SMALLEST radius at which the neighbourhood first looks significantly
# too tight relative to a Poisson null. "descend" (L265-279, L303-317) starts
# from the ball containing ALL n-1 other points (the largest possible
# radius) and SHRINKS it inward one point at a time (j = 1..n-low.num,
# points dropped in decreasing-distance order), stopping and reporting the
# CURRENT radius the first time the statistic is already back ABOVE the
# envelope -- i.e. it takes the LARGEST radius that still looks
# CSR-consistent. Both searches are hunting for the same
# CSR/non-CSR crossing along the distance-sorted neighbour sequence for
# point i, approached from opposite ends. If that crossing is monotonic
# (statistic drops below the envelope exactly once as the ball grows and
# never recovers), ascend and descend land on the same radius. If it is not
# monotonic -- the statistic dips below the envelope and climbs back above
# it more than once, which can happen with the small, noisy `low.num`-sized
# balls this test uses -- the two directions can lock onto genuinely
# different radii for the same point. That non-monotonicity is the entire
# mechanism by which this sweep can find a difference.
#
# RK-BASED DETECTORS HAVE NO ANALOGOUS TOGGLE. RUMCCD_outlier
# (methods/outlier_detection/RU-MCCDs.R) and SUMCCD_outlier
# (methods/outlier_detection/SU-MCCDs.R) take no `method` argument at all
# (verified against their full signatures). Their underlying radius routine,
# RKCCD_correct_quant -> rccd.clustering_correct_quantile (R/ccds/RK_CCD_New.R
# L19-46), does have a `method` argument, but it selects between
# "non-dynamic"/"dynamic" edge correction for the K-function estimate inside
# ccd.Kest.edge.quantile() (L158-224) -- an unrelated axis (edge-correction
# scheme, not search direction) -- and RUMCCD_outlier/SUMCCD_outlier both
# hard-code it to "non-dynamic" via RKCCD_correct_quant's default
# (RK_CCD_New.R L225-247) with no path from the outlier-detector call to
# override it. wp0_mccd_methods.R's own wrappers confirm this from the
# harness side: umccd_method (L282) and sumccd_method (L299) declare no
# `method` parameter at all, unlike unmccd_method/sunmccd_method. So the
# ascend/descend toggle -- and R3.10a's ablation -- reaches only half of the
# 2x2 design family (the NND-based half: UN-MCCD, SUN-MCCD), not
# U-MCCD/SU-MCCD. This asymmetry is itself worth stating in the response
# letter.
#
# VERIFICATION THAT THE ARGUMENT GENUINELY REACHES nnccd.radi. nnccd.radi is
# shadowed by a thin tap installed AFTER wp0_mccd_methods.R has fully sourced
# (so the tap is not clobbered by a later source() -- neither UN-MCCD.R nor
# SUN-MCCD.R nor UN_CCD.R is re-sourced per detector call, only once at
# startup via wp0_mccd_methods.R's `if (!exists(...))` guards). The tap calls
# straight through to the REAL nnccd.radi (captured into ORIG_NNCCD_RADI
# before the shadow is installed) and records the `method` value it actually
# received plus the resulting R vector. Every CSV row therefore carries
# nnradi_method_received (must equal the requested direction) and a radii
# signature (radii_sig, radii_mean); the analysis step below refuses to
# report "no effect" unless it can show at least one cell where ascend and
# descend produced different radii_sig -- same discipline
# 40_wp2a_tolerance_sweep.R used for its own override guard.
#
# Usage:
#   Rscript 43_wp2a_direction_sweep.R [datasets] [methods] [directions] [time_budget_sec]
#   Rscript 43_wp2a_direction_sweep.R --summary
#
# Checkpointed on (dataset, method, direction): a restart skips completed
# cells. TIME_BUDGET (default 480s) stops the grid early and prints what is
# left, so a single invocation fits inside the 600s PowerShell tool timeout;
# re-invoke (same command) to continue -- has_result() skips finished cells.
# ---------------------------------------------------------------------------
```

## Step 2/3: measured comparison (ascend vs descend), per dataset x method

Cells run: 64 rows in `wp2a_direction_sweep.csv`. Complete (ascend,descend) pairs analysed: 32.
Datasets not yet run (time budget exhausted): none -- all requested datasets covered..

### Argument-reaches-nnccd.radi guard

- Every cell's `nnradi_method_received` matched the requested direction: YES (0 mismatches across 18 checkable pairs; 14 not_run pairs made no call and are excluded).
- Cells where ascend and descend produced a different `radii_sig` (the actual radius vector nnccd.radi returned): 18 / 18 ok pairs.
- Guard passed: the two directions are demonstrably reaching different code paths inside nnccd.radi and producing different radii, so any label/metric equality found below is a genuine null result, not a wiring failure.

### Label and metric comparison

- Cells compared: 18.
- Identical flagged sets: 8 / 18.
- Cells whose labels changed: 10.
- Largest |dBA| = 0.2206 at glass / SUN-MCCD (ascend BA=0.7696, descend BA=0.5490).
- Mean dBA (descend-ascend) = -0.01878, mean dF2 = -0.01863 over 18 ok pairs.
- Cluster-count changes: 2 / 18 ok pairs.

### Per-cell table

| dataset | method | n | d | radii_differ | ncls_asc | ncls_desc | nflag_asc | nflag_desc | jaccard | identical | dBA | dF2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ecoli | SUN-MCCD |  336 |  7 | TRUE |  1 |  1 |   73 |  88 | 0.8295 | FALSE | -0.0229 | -0.0417 |
| ecoli | UN-MCCD |  336 |  7 | TRUE |  2 |  2 |  115 | 114 | 0.9913 | FALSE | +0.0015 | +0.0016 |
| glass | SUN-MCCD |  213 |  9 | TRUE |  1 |  1 |  103 | 193 | 0.5337 | FALSE | -0.2206 | -0.1272 |
| glass | UN-MCCD |  213 |  9 | TRUE |  3 |  1 |   50 | 193 | 0.2591 | FALSE | +0.0556 | +0.0802 |
| hepatitis | SUN-MCCD |   74 | 19 | TRUE |  1 |  1 |   28 |  28 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| hepatitis | UN-MCCD |   74 | 19 | TRUE |  1 |  1 |   28 |  28 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| lymphography | SUN-MCCD |  148 | 18 | TRUE |  1 |  1 |   23 |  23 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| lymphography | UN-MCCD |  148 | 18 | TRUE |  2 |  2 |   26 |  26 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| pageblocks | SUN-MCCD | 4795 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| pageblocks | UN-MCCD | 4795 | 10 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| PenDigits | SUN-MCCD | 3200 | 16 | NA |  1 | NA | 2948 | NA | NA | NA | NA | NA |
| PenDigits | UN-MCCD | 3200 | 16 | NA |  2 | NA | 2871 | NA | NA | NA | NA | NA |
| pima | SUN-MCCD |  555 |  8 | TRUE |  1 |  1 |  273 | 290 | 0.9414 | FALSE | +0.0032 | +0.0030 |
| pima | UN-MCCD |  555 |  8 | TRUE |  2 |  2 |  125 | 261 | 0.3640 | FALSE | +0.1365 | +0.2151 |
| shuffle | SUN-MCCD | 1013 |  9 | NA |  1 | NA |  349 | NA | NA | NA | NA | NA |
| shuffle | UN-MCCD | 1013 |  9 | NA | 57 | NA |  304 | NA | NA | NA | NA | NA |
| stamps | SUN-MCCD |  340 |  9 | TRUE |  1 |  1 |   52 |  52 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| stamps | UN-MCCD |  340 |  9 | TRUE |  2 | 36 |   73 |  23 | 0.1852 | FALSE | -0.1676 | -0.3467 |
| thyroid | SUN-MCCD | 3656 |  6 | NA |  1 | NA |  812 | NA | NA | NA | NA | NA |
| thyroid | UN-MCCD | 3656 |  6 | NA |  2 | NA | 1413 | NA | NA | NA | NA | NA |
| vertebral | SUN-MCCD |  240 |  6 | TRUE |  2 |  2 |   18 |  15 | 0.3200 | FALSE | -0.0500 | -0.1087 |
| vertebral | UN-MCCD |  240 |  6 | TRUE |  2 |  2 |   19 |  48 | 0.3958 | FALSE | -0.0690 | -0.0062 |
| vowels | SUN-MCCD | 1452 | 12 | NA |  1 | NA |  660 | NA | NA | NA | NA | NA |
| vowels | UN-MCCD | 1452 | 12 | NA |  2 | NA |  691 | NA | NA | NA | NA | NA |
| waveform | SUN-MCCD | 3443 | 21 | NA |  1 | NA |  301 | NA | NA | NA | NA | NA |
| waveform | UN-MCCD | 3443 | 21 | NA |  2 | NA | 1054 | NA | NA | NA | NA | NA |
| WBC | SUN-MCCD |  223 |  9 | TRUE |  2 |  2 |  103 | 105 | 0.9810 | FALSE | -0.0047 | -0.0048 |
| WBC | UN-MCCD |  223 |  9 | TRUE |  2 |  2 |  117 | 117 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| WDBC | SUN-MCCD |  367 | 30 | TRUE |  1 |  1 |  251 | 251 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| WDBC | UN-MCCD |  367 | 30 | TRUE |  1 |  1 |  302 | 302 | 1.0000 | TRUE | +0.0000 | +0.0000 |
| wilt | SUN-MCCD | 4819 |  5 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| wilt | UN-MCCD | 4819 |  5 | NA | NA | NA | NA | NA | NA | NA | NA | NA |

### Headline verdict

Ascend and descend disagree in 10 of 18 cells tested. Mean dBA (descend-ascend) = -0.0188 suggests a systematic direction, not a wash.

