# WP8 verification (Opus, 2026-09-05) — verdicts and required changes

Reviewed commit `5ed690c`: `WP8_PROTOCOL.md`, `85_wp8_boundary_fp.R`, `86_wp8_outlier_types.R`,
`87_wp8_small_cluster.R`, `88_wp8_csr_violation.R`, smoke outputs. Configuration (S_min = 0.05,
alpha via paper resolvers: RK 99/999, UN 90/99, SUN 90/999 at d = 3/10; `get_simul(n = nrow(X))`;
baselines at study settings; generators byte-faithful to `55`) all PASS. The items below are what
must change before the real run.

## Verdicts

| Script | Verdict |
|---|---|
| 85 boundary FP | PASS WITH REQUIRED CHANGES |
| 86 outlier types | **FAIL** (`local` and `bridge` generators do not produce the named phenomenon) |
| 87 small cluster | PASS WITH REQUIRED CHANGES (cluster 3 confounds size with density; ~6.7 h cost) |
| 88 CSR violation | PASS WITH REQUIRED CHANGES |

## Cross-cutting (all four scripts)

- **Seed bug (FAIL, all four).** `s_i` is the row index of a settings table built from the `dims`
  CLI argument, so the same setting gets a different seed depending on how the run is chunked
  (e.g. 85 gaussian_d3 rep 7: DIMS="3,10" -> seed 308508; DIMS="3" -> 208508). `has_result` keys
  do not include `seed`, so chunked runs silently mix two data sets under one `setting_id`.
  Fix: build `CANON <- build_settings(<full declared grid>)` once; `SETTINGS <- CANON[CANON$d %in%
  DIMS, , drop = FALSE]`; `s_i <- match(s$setting_id, CANON$setting_id)` before the seed line.
  Locations: 85:224,227,234; 86:231,234,241; 87:132,135-136,143; 88:162,165,172.
- **Add a settings selector and a `budget` (seconds) argument** mirroring
  `55_wp2c_simulation_arm.R` (`run [settings] [reps] [budget] [cores]`) so <= 10-minute chunks
  are reachable without varying `dims`.
- **Block writes.** `append_result` accepts a list of equal-length vectors (multi-row); verified
  6.4 ms per 189-row block vs 773 ms row by row. Write each cell's cluster rows (and 85's bin
  rows) as ONE block. ~1.27 M single-row appends otherwise, on the G: enclosure.
- **`--summarize` mode in every script (FAIL, none has one).** Means AND standard errors over
  reps; filter `status == "ok"`; de-duplicate on the key keeping the last occurrence (retried
  cells leave an error row and an ok row).
- Document the `detected_cluster = NA` (unassigned) convention in the protocol's output schema.
- Data are regenerated once per method (not per rep); harmless, leave it.

## 85_wp8_boundary_fp.R

1. Seed canonicalisation (above).
2. After line ~168 add `stopifnot(length(res$score) == dat$n, !anyNA(res$score))` (85 is the only
   script without the guard, and its keys-only done file can never trip the NA check).
3. **Equal-width `r/r_max` bins fail at d = 10**: 66% of regular points land in bin 10 and bins
   1-5 hold ~0.1% combined, so there is no interior to contrast. Add an equal-mass binning
   alongside: new column `bin_type`, 20 rows per cell: `"width"` = current `ceiling(10*rn)`,
   `"mass"` = `ceiling(10 * rank(r, ties.method = "first") / length(r))`. Update the protocol
   schema. The rank version also fixes the unbounded `r_max` in the Gaussian arm.
4. Lines ~173-189: write the bin rows as one block and the 189 cluster rows as one block (also
   removes the partial-write duplication window; the done marker currently gates only the marker).
5. Add `unassigned_rows = res$unassigned_rows` to the cluster block or done row.
6. `--summarize`: per (setting_id, d, method, bin_type, bin): `sum(n_flagged)/sum(n_regular)` as
   the point estimate (never a mean of per-rep ratios), plus mean and SE of per-rep ratios over
   non-empty reps and the rep count; de-duplicate; drop bin rows whose cell has no done marker.

## 86_wp8_outlier_types.R (FAIL)

Measured: `local` outlier NN distance / median regular NN distance = 1.28 (d=3), **0.61 (d=10)**;
the 0.3 clearance never bites (0 forced keeps per rep); smoke shows TPR = 0 for all nine methods
including LOF. `bridge`: 6.2/10 (d=3) and 5.2/10 (d=10) labelled outliers lie INSIDE a cluster's
realised ball, capping TPR near 0.4-0.5 for every method. `collective` is fair (gap 0.67-0.69).

1. `local` (lines ~96-97): replace the fixed clearance with a shell construction — host cluster
   drawn as a dense core of radius `0.5*scale`; each local outlier at a uniform direction with
   radius `runif(1, 0.75, 1.0) * scale` (inside nominal support, outside the core). **Mandatory
   acceptance check before launch:** `mean(outlier NN dist) / median(host regular NN dist)` over
   20 reps at d = 3 and d = 10 must be >= 2; record the values in the protocol; abort otherwise.
2. `bridge` (lines ~116-119): span only the gap between the realised cluster surfaces. Expose the
   jitter (`attr(pts, "scale") <- scale` in `draw_cluster`, line ~78); with `s1`, `s2`:
   `t_lo <- (s1+0.25)/CLS_DIS`, `t_hi <- 1-(s2+0.25)/CLS_DIS`, `t_i <- t_lo + (i-0.5)/n0*(t_hi-t_lo)`.
   Assert no bridge point is within `s_k` of `mu_k`; record `n_inside` per row.
3. `collective` (line ~133): `center <- mu1 + u * (s1 + 0.5)` (realised radius, fixed 0.5 stand-off);
   correct the protocol's "just past the cluster's typical edge" to the measured stand-off.
4. Seed canonicalisation; 5. cluster rows as one block; 6. `--summarize` (mean and SE of
   TPR/TNR/BA/F2 per type x d x method).
7. Re-run `--smoke` at BOTH d = 3 and d = 10 for all three types; require non-zero TPR for at
   least LOF on `local`.

## 87_wp8_small_cluster.R

1. Cluster 3 gets the same radius as clusters 1-2, so lowering `m` lowers density, not just size
   (NN spacing ratio vs background at d=3: 4.42 at m=2, 2.26 at m=10, 1.71 at m=20). Density-match:
   `data3 <- rpoisball.unit(m, d) * (runif(1, R_MIN, R_MAX) * (m/n1)^(1/d)) + mu3`. Record the
   deviation as a dated appended note in the protocol.
2. Record `n_unassigned_third = sum(is.na(res$cluster[third_idx]))` and
   `singleton_lost = res$singleton_lost_rows` per row (frac_flagged = 0 otherwise conflates "own
   cluster recognised" with "never covered"; `mccd_translate` scores unclaimed rows 0).
3. Seed canonicalisation against a canonical SIZES x c(3,10) table.
4. Add a WP7 cluster file `87_small_cluster_clusters.csv` (same schema, `true_cluster` in
   {1,2,3}, one block per cell).
5. Add an `m` selector argument (4th positional, comma list) for chunking.
6. `--summarize`: mean and SE of `frac_flagged` and `n_clusters` per (m, d, method), plus the
   fraction of reps with `n_clusters == 3`.
7. **Cost ~6.7 h**, 96% of it U-/SU-MCCD at d = 3 (11 s each). Decision taken: run the two RK
   methods at 50 reps and the two NND methods at 100 reps; state it in the protocol.

## 88_wp8_csr_violation.R

1. `gauss_unequal` confounds anisotropy with the Gaussian's own radial gradient. Decision taken:
   add a fourth generator `gauss_equal` (`sigma <- rep(1, d)`) so anisotropy is identified by the
   `gauss_unequal - gauss_equal` contrast; fix the protocol text (lines ~219-223).
2. DBSCAN is fed `Y = rep(1, n)` so its oracle contamination is 0 and its flag rate is 0 by
   construction. Keep it but write `note = "oracle contamination = 0; flag rate 0 by construction"`
   and footnote it in any table.
3. MST at fixed `thresh = 1.2` flags 39.5% of pure CSR; the paper's real-data script tunes thresh
   per data set (1.05-1.6). Decision taken: sweep `thresh in {1.05, 1.2, 1.4, 1.6, 2.0}` for MST
   in this experiment (cheap; 0.07 s per fit) and report the best per generator plus the fixed
   1.2 value, labelled.
4. Seed canonicalisation; 5. 200 cluster rows as one block; 6. `--summarize` (mean and SE of
   flag_rate per generator x d x method, plus the paired `generator - uniform_control` delta at the
   same seed with its own SE).

## Cost (measured, one core, n = 199)

9-method cell: 23.3 s at d = 3 (U/SU-MCCD 11 s each), 7.3 s at d = 10. As written: 85 2.1 h,
86 3.1 h, 87 6.7 h, 88 3.1 h (15 h); with block writes and the 87 decision above, roughly 11 h.
Long grids are launched from the main session as background jobs (per-cell append, resumable),
not by a subagent.
