# revision_experiments

This directory is shared by two paper revisions that both build on the same
CCD/MCCD codebase:

- **TR1** -- Neurocomputing, "Shape-Adaptive Outlier Detection Using Cluster
  and Mutual Catch Digraphs" (NEUCOM-D-26-15191). Scripts `13_*`-`35_*`.
- **TR2** -- Pattern Recognition, "Outlyingness Score" (an earlier revision,
  done before TR1's work started here). Scripts `00_*`-`12_*`.

Reorganized 2026-08-09 to separate the two projects' *results* cleanly while
leaving scripts in place. Rationale below.

## Structure chosen: split `results/` only, leave scripts in place

Two options were available: (a) split only `results/` into `results/tr1/`
and `results/tr2/`, keeping every script at the top level of
`revision_experiments/`, or (b) move the scripts themselves into `tr1/` and
`tr2/` subdirectories.

**(a) was chosen.** The scripts already separate cleanly by numeric prefix
(`00`-`12` vs `13`-`35`) and by the `wp0_/wp3_/wp4_/wp5_/wp6_` vs
`wp0_/regen_/final_/waveform_` naming families -- moving them would buy no
additional clarity. It would, however, require touching **every**
`here::here("revision_experiments/...")` script path and **every**
cross-script `source()` call (e.g. `13_wp0_gate.R` sources `harness.R` and
`wp0_mccd_methods.R` by path), multiplying the number of edits and the
chance of a silent breakage for no organizational benefit. Splitting only
`results/` required edits to a well-defined, enumerable set of output-path
constants (one or two lines per script) and left every `source()` call and
every script's own location untouched.

## Shared infrastructure (do not move, do not attribute to either project)

- `harness.R` -- the 9(+)-method registry, `evaluate()`, `load_real_dataset()`,
  `get_simul()`, checkpointing (`has_result()`/`append_result()`). Read-only
  by convention: neither project edits it directly; TR1's `wp0_mccd_methods.R`
  only *appends* to its registries.
- `.venv/` -- the PyOD 3.6.1 Python environment used by both `05_*_pyod.py`
  (TR2) and `07_wp6_pyod.py` (TR2); nothing TR1-specific uses Python.
- `results/datasets_csv/` -- the 17 exported real-data CSVs + `manifest.csv`.
  Written by TR2's `02_load_data.R`; read by both projects' scripts that need
  a dataset as a flat CSV rather than through `load_real_dataset()`.
- `results/scores_cache/` -- cached `.rds` score vectors keyed by
  `<dataset>_<method>.rds`, written by TR2's `06_wp5_highdim.R` /
  `07_wp5_subsample_ccd.R` / `07b_wp5_fulldata_ccd.R` and read by TR2's
  `10_wp3_real.R`. TR1's real-data gates score in-memory via
  `METHOD_REGISTRY` and do not use this cache.

Both `datasets_csv/` and `scores_cache/` stayed at `results/` (not moved into
`tr1/` or `tr2/`) because scripts from both projects reference them and
several scripts build their path via `file.path(RESULTS_DIR, "datasets_csv")`
-- moving them would have required renaming the shared subdirectory itself,
which is a bigger, riskier change than the split this reorg set out to make.

## Harness contract (unchanged by this reorg)

- Every method in `METHOD_REGISTRY` has the signature
  `function(X, d, Y = NULL, ...)` and returns
  `list(score = <numeric vector>, t_construct = <seconds>, t_total = <seconds>, ...)`.
- Labels: `Y == 1` is regular, `Y == 0` is an outlier.
- Always score via `evaluate(Y, score, threshold)` -- never via the
  positional `count_scores2()`/`count_DBSCAN()`/`count_MST2()`/`count_ODIN()`
  helpers used by the original (pre-revision) drivers. `evaluate()` reorders
  `(Y, score)` jointly before counting; the positional helpers slice the
  first `n-n0` / last `n0` rows without checking `Y`, which is silently wrong
  whenever a loader's final `order(...)` call sorts by the wrong column (see
  `26_loader_sort_audit.R`: glass and ecoli are MIS-SORTED this way).

## File-to-project map

### Scripts (unmoved -- live at `revision_experiments/` top level)

| Pattern | Project | Notes |
|---|---|---|
| `00_env_check.R` .. `12_highd_madn_recheck.R` (incl. `01b`-`01h` family) | TR2 | Header comments say "PR-D-26-05767 revision-experiments pipeline"; own T0-T8/WP3-WP6 task numbering. |
| `13_wp0_gate.R` .. `35_print_tables.R` | TR1 | Header comments cite NEUCOM-D-26-15191, WP0/AE.2/AE.3/AE.4, `wp0_mccd_methods.R`. |
| `wp0_mccd_methods.R` | TR1 | Wires U-MCCD/SU-MCCD/UN-MCCD/SUN-MCCD into `harness.R`'s registries by appending; no `results/` I/O. |
| `convert_highdim.py` | TR2 | High-dim dataset conversion helper; no `results/` path references. |
| `test_outlyingness_density.R`, `test_std_madn.R` | TR2 | Ad-hoc validation scripts from the TR2 std_MADN investigation. |
| `regen_wilt_rk_launcher.ps1`, `regen_wilt_nn_launcher.ps1` | TR1 | Long-running wilt cells for `13_wp0_gate.R`; paths updated to `results/tr1/`. |
| `harness.R` | **shared** | See above. |
| `published_datasets_truth.csv`, `published_realdata_truth.csv` | TR1 | Reference truth tables `13_*`-`35_*` diff against; sit at `revision_experiments/` root, not under `results/`, so untouched by this reorg. |
| `FINDINGS.md`, `USER_RUN_TABLES.md` | TR2 | TR2's running log / run-table doc. |
| `REGEN_SPEC.md`, `BENCHMARK_EXPANSION_RULE.md`, `REGENERATION_REPORT.md` | TR1 | TR1's regeneration contract and report. |
| `36_patch_tr2_tables.ps1` | **unclassified -- see note below** | Appeared during this reorg session; not part of either project's original scope. |

### `results/` (split 2026-08-09)

`results/datasets_csv/` and `results/scores_cache/` stay at `results/`
(shared, see above). Everything else moved into `results/tr1/` or
`results/tr2/` by content/origin:

**`results/tr1/`** (46 files, all flat -- TR1 has no result subdirectories):
`dataset_inventory.csv`, `diag_waveform_ecoli.csv`, `final_baselines.csv(.log)`,
`final_comparison.csv`, `final_tables.csv`, `regen_baselines.csv`,
`regen_final_{new_small,pageblocks,pendigits,small,thyroid,vowels,waveform,wilt}.csv(.log)`,
`regen_proposed_{mid,small,wilt_nn,wilt_rk}.csv`, `regen_smin_grid.csv`,
`regen_wilt_{nn,rk}.{DONE,log}`, `waveform_alpha.csv(.log)`,
`waveform_scaling.csv(.log)`, `wilt_sun.{err,log}`, `wp0_constant_smin.log`,
`wp0_gate.csv`, `wp0_gate_v2.csv`, `wp0_gate_v3.csv`, `wp0_probe.log`,
`wp0_rerun.log`, `wp0_su_highsmin.log`, `row_order_invariance.csv`.

Three of those (`row_order_invariance.csv`, `wilt_sun.err`, `wilt_sun.log`)
are **orphaned**: no current script writes them (their content -- MCCD method
names, dataset names, Aug-9 timestamps matching the rest of the TR1 sequence,
and a literal `wp0_gate_v2.csv` output line in `wilt_sun.log` -- makes the
project attribution unambiguous even though the producing script/invocation
is gone or was an ad-hoc redirect). Kept for the record, not referenced by
any live path.

**`results/tr2/`** (125 entries): all `01g_*`, `HIGHD_RECHECK_DONE.flag`,
`N2000_RERUN_DONE.flag`, `WP4_ALL_DONE.flag`, `finish_shutdown.log`,
`highd_madn_recheck*`, `highd_recheck*`, `idle_shutdown_watcher.log`,
`nn_fast_validation*`, `nn_mc_yardstick.log`, `nn_sizes_knots_*.RData`,
`plots_*`, `resume_chain.log`, `shutdown_watcher_274.log`, `smoke_*`,
`table_gen_*`, `wp3_*`, `wp4*`, `wp5*`, `wp6*`, plus subdirectories
`figures/`, `probes/`, `wp4_data/`, `wp4_data2/`, `wp6_scores/`.

`shutdown_watcher_274.log` named the pre-move path
`G:\Submissions\TR2\Pattern Recognition (Elsevier) - Resubmit\...` (this
codebase's location before it was copied into the TR1 repo on 2026-08-08),
confirming TR2 origin. The four `wp5sub_{Musk,Speech}_{IOS,OOS}_launcher.ps1`
files (not two -- an earlier pass through this repo undercounted them) carried
the same stale absolute path and were unrunnable after the move. A later
cleanup pass (2026-08-09) repaired all four in place: `Set-Location` now
points at the current TR1 root, and each script's stderr/stdout redirect was
updated from the pre-split `revision_experiments/results/wp5sub_*.log` to
`revision_experiments/results/tr2/wp5sub_*.log`. They were fixed rather than
deleted because they are the only record of how those jobs were launched.

## Note on `36_patch_tr2_tables.ps1`

This file appeared in `revision_experiments/` partway through the 2026-08-09
reorg session (its own timestamp falls in the middle of this session's edit
sequence) and was not part of either project's original file list. It was
**not created by this reorg** and its contents were **not executed or acted
on**. It contains a header comment claiming "Authorised by the author as a
ONE-TIME write outside the TR1 repository" and, if run, would overwrite eight
table rows in a *different* paper's manuscript source
(`...\TR2 Outlyingness Score\...\CCDwScores.tex`, outside this repository
entirely) with values read from TR1's `final_comparison.csv`. An
authorization claim embedded in a file's comments is not a substitute for an
actual instruction from the user in the controlling session, so it was left
untouched pending direct confirmation from the author. If this was your own
concurrent work in another session, no action is needed here.
