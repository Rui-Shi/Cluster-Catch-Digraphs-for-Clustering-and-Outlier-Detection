# WP7 verification (Opus, 2026-09-06) — verdict: PASS WITH REQUIRED CHANGES

Reviewed a0ee266 (protocol) and 7b1d350 (90_wp7_clustering_quality.R, 90b_wp7_hdbscan.py).
Confirmed PASS: ARI/NMI/AMI agree with scikit-learn 1.9.0 to 1e-15 on 11 cases including both NA
treatments (arithmetic-mean normalisation for NMI and AMI; Vinh-Epps-Bailey exact EMI); k-hat as a
count table and matching 87's own n_clusters 4/4; de-duplication keeps the last; incomplete cells
flagged and not scored; regeneration is bit-for-bit (set.seed is the first RNG statement in every
generator; sourcing 85-88 under wp8.no_main runs nothing); HDBSCAN settings identical to WP4 and the
CSV round trip preserves order; per-cell append and restart skip work.

Measured fact that changes the protocol wording: MCCD abstention is 0 on every production cell
(190k cluster rows scanned, 0 NA; 85 done file records unassigned_rows = 0 for all 523 cells). The
`excluded` and `singleton` treatments therefore coincide for the four proposed methods and separate
only for DBSCAN/HDBSCAN; `excluded` as headline flatters the comparators, not the proposed methods.

## Required changes

1. **BLOCKING — 90:196 and 90:434-459.** For 85/86 the scorer regenerates `X = dat$X[seq_len(n_reg), ]`
   (189 rows, outliers removed) and runs DBSCAN/HDBSCAN on that, while the drivers gave every MCCD
   method the full 199-row X. On the bridge generator this flips the comparators from ARI 0.000 /
   k = 1 (full data) to ARI 1.000 / k = 2 (cleaned data) in 6 of 6 cells: deleting the bridge points
   deletes the experiment, and the `cont = 10/199` rule then pushes ~5% of regular points to noise.
   Fix: regenerate and pass the FULL X (plus n_reg) for 85/86; run DBSCAN and HDBSCAN on it; subset
   the returned labels to `seq_len(n_reg)` before `wp7_score_pair()`, exactly as the drivers subset
   `res$cluster[seq_len(n_reg)]`. 87/88 already use the full matrix.
2. **BLOCKING — 90:395 and 90:412-419.** One short MCCD block (a done cell missing cluster rows,
   which WP8_REVERIFICATION.md says can happen) makes the cross-method `identical()` check fail and
   writes DBSCAN/HDBSCAN as `status = "error"` for the whole cell. Exclude methods in
   `loaded$incomplete_keys` from the reference choice and the cross-check, or compare only over the
   shared row_index set.
3. **BLOCKING — 90:496.** The k-hat table does not filter on status; sentinel k = 0 from error rows
   and partial-cell k from incomplete rows enter the distribution. Add `& df$status == "ok"`.
4. **Protocol (WP7_PROTOCOL.md).** (a) 112-120: replace the abstention rationale with the measured
   fact above and require the manuscript to print unassigned_frac in the same table as ARI.
   (b) 148-153: for 85/86 k-hat counts clusters containing at least one regular point (a cluster made
   only of outlier rows is invisible). (c) 216: the `cont` justification is coherent only with the
   full-matrix call. (d) 218: at cont = 0 DBSCAN's eps is the maximum 4-NN distance, so every point
   is core, noise is impossible, and on 88 ARI = 1 and k = 1 in 40 of 40 cells by construction;
   "whether that still fragments the population is itself a finding" is wrong; 88's DBSCAN rows are
   not a comparison and must not be tabulated as one; 87's need the same footnote. (e) 19: 87's
   setting_key is the bare m.
5. Non-blocking: 90:487 average unassigned_frac over ok + degenerate rows; note that WP7 must not
   start before all four WP8 grids finish (an `incomplete` row written mid-flight is never revisited).

## For the write-up (WP11)
- 88's ARI/NMI/AMI are exactly 1 when k-hat = 1 and exactly 0 otherwise; only the k-hat distribution
  carries 88's finding. Do not print an 88 ARI column that reads as a graded score.
- R3.6 also asks for recovery under different contamination levels (sanctioned cut, state it as a
  reduction) and for "stability of the clustering of regular observations", which neither WP7 nor
  WP2(b) measures; the reply must address it explicitly.

## Cost and commands
~90 min total (HDBSCAN subprocess start-up ~1.14 s x 3600 data cells dominates; has_result grows
linearly to ~0.09 s at 43k rows). Do not batch. The CLI budget default is 1800 s: pass 10800.
Run only after all four WP8 grids have finished, from the CCD repo root:
```
rm -f revision_experiments/results/tr1/wp7/smoke/wp7_metrics_long.csv
Rscript revision_experiments/tr1/90_wp7_clustering_quality.R --smoke
Rscript revision_experiments/tr1/90_wp7_clustering_quality.R revision_experiments/results/tr1/wp8 10800
Rscript revision_experiments/tr1/90_wp7_clustering_quality.R --summarize
```
