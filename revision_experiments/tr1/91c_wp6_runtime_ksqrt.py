"""91c_wp6_runtime_ksqrt.py -- WP6 addendum: mutual-kNN and SNN timed at
k = round(sqrt(n)), the label-free rule of the accuracy comparison
(81_wp4_baselines.py / 82_wp4_metrics.R since 2026-09-26).

91b_wp6_runtime_py.py timed these two detectors on the grid k in
{5,10,15,20,30} and WP6_PROTOCOL.md reduced them to k=10. The accuracy
comparison later moved its primary k to round(sqrt(n)), so the runtime tables
and the accuracy tables used different k. This script times the two detectors
at the primary k on the same nine exported data sets, with the same timed
region and memory accounting, and re-times k=10 in the same session as a
control against the published 91b figures. It was written on 2026-10-09,
after 91b's results were known; it is an addition for consistency, not part
of the pre-declared protocol.

Everything that defines a measurement is imported from 91b unchanged:
timed_fit() (fit only, time.perf_counter, tracemalloc peak), load_cell(),
discover_files(), the single-thread pinning done at 91b's import, and the
untimed warm-up fit per (method, k) at the smallest data set. The two k
values are interleaved within each replicate so that drift affects both
alike. Rows are appended to the raw CSV as each fit finishes.

One difference from 91b: an untimed, discarded warm-up fit per (method, k) is
also made on every data set before its measured replicates. A first run
without it put a one-time first-call cost on the very first measured fit
(mutual-kNN, k=22, n=500, d=100, replicate 1: 0.0386 s and 1.18 MB against
0.0062 s and 0.96 MB for replicates 2-10), because the warm-up at the
smallest data set (d=10) does not exercise the code path used at d=100. In
91b that cost fell on methods timed earlier in the same cell. The grid was
then rerun in full with the per-data-set warm-up; the first run is not kept.

At n=100, round(sqrt(n)) = 10, so the two settings coincide there.
(n, k): (100, 10), (250, 16), (500, 22), (1000, 32), (2000, 45); no value of
sqrt(n) on this grid falls on a .5 tie, so Python's round() and R's round()
agree.

Usage (from the repo root, pinned venv):
  revision_experiments/.venv/python.exe revision_experiments/tr1/91c_wp6_runtime_ksqrt.py [--reps 10]
  revision_experiments/.venv/python.exe revision_experiments/tr1/91c_wp6_runtime_ksqrt.py --summarize

Outputs (results/tr1/wp6/):
  91c_wp6_runtime_ksqrt_raw.csv      one row per timed fit
  91c_wp6_runtime_ksqrt_summary.csv  median / mean / CV per (n, d, method, rule)
"""

import argparse
import importlib.util
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent            # revision_experiments/tr1/

# Importing 91b pins the BLAS/OpenMP thread counts before numpy and torch are
# loaded, exactly as in the published run.
_spec = importlib.util.spec_from_file_location("wp6_91b", HERE / "91b_wp6_runtime_py.py")
b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b)

import pandas as pd  # noqa: E402  (after 91b, which sets the thread env vars)

WP6 = b.WP6
RAW = WP6 / "91c_wp6_runtime_ksqrt_raw.csv"
SUMMARY = WP6 / "91c_wp6_runtime_ksqrt_summary.csv"
METHODS = ["MutualKNN", "SNN"]
RAW_COLS = ["method", "rule", "k", "n", "d", "grid", "cell_value", "rep",
            "time_s", "mem_peak_mb", "threads", "status"]


def k_sqrt(n):
    return int(round(math.sqrt(n)))


def append_raw(row):
    df = pd.DataFrame([row], columns=RAW_COLS)
    df.to_csv(RAW, mode="a", header=not RAW.exists(), index=False)


def done_keys():
    if not RAW.exists():
        return set()
    df = pd.read_csv(RAW, dtype=str, keep_default_na=False, na_filter=False)
    df = df[df["status"] == "OK"]
    return {(r.method, r.rule, r.grid, r.cell_value, r.rep) for r in df.itertuples()}


def summarize():
    df = pd.read_csv(RAW)
    df = df[df["status"] == "OK"]
    rows = []
    # Grouped by (n, d, method, rule) as in 91_wp6_runtime.R's summariser, so
    # the n=500, d=10 data set, exported once for each grid, is pooled (20 fits).
    for (n, d, method, rule), g in df.groupby(["n", "d", "method", "rule"], sort=True):
        t, m = g["time_s"], g["mem_peak_mb"]
        rows.append({
            "n": n, "d": d, "method": method, "rule": rule,
            "k": int(g["k"].iloc[0]), "n_points": len(g),
            "median_time_s": t.median(), "mean_time_s": t.mean(),
            "cv_time": t.std(ddof=1) / t.mean() if len(t) > 1 else float("nan"),
            "median_mem_peak_mb": m.median(),
            "cv_mem": m.std(ddof=1) / m.mean() if len(m) > 1 and m.mean() > 0 else float("nan"),
        })
    out = pd.DataFrame(rows)
    out.to_csv(SUMMARY, index=False)
    b.log(f"[summarize] wrote {SUMMARY} ({len(out)} rows)")
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        b.log(out.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()

    if args.summarize:
        summarize()
        return

    b.idle_check()          # report only; the fits take milliseconds
    files = b.discover_files(WP6 / "data")
    if not files:
        sys.exit(f"No data exports in {WP6 / 'data'}; run 91_wp6_runtime.R first.")
    cells = [(grid, cv, b.load_cell(p)) for grid, cv, p in files]

    # untimed warm-up, one fit per (method, k) that will be measured, on the
    # smallest data set (same rule as 91b); discarded
    X_warm = min((X for _, _, X in cells), key=lambda X: X.shape[0] * X.shape[1])
    for method in METHODS:
        for k in sorted({10} | {k_sqrt(X.shape[0]) for _, _, X in cells}):
            b.timed_fit(method, X_warm, k=min(k, X_warm.shape[0] - 1))
    b.log("warm-up done")

    done = done_keys()
    for grid, cell_value, X in cells:
        n, d = X.shape
        plan = [("ksqrt", k_sqrt(n)), ("k10", 10)]
        b.log(f"\n==== grid={grid} cell={cell_value} n={n} d={d} k_sqrt={k_sqrt(n)} ====")
        # untimed warm-up on this data set, one fit per (method, k); discarded
        for method in METHODS:
            for _rule, k in plan:
                b.timed_fit(method, X, k=k)
        for rep in range(1, args.reps + 1):
            for method in METHODS:
                for rule, k in plan:
                    key = (method, rule, grid, cell_value, str(rep))
                    if key in done:
                        continue
                    try:
                        elapsed, mem = b.timed_fit(method, X, k=k)
                        status = "OK"
                    except Exception as e:  # noqa: BLE001
                        elapsed, mem, status = float("nan"), float("nan"), f"ERROR: {type(e).__name__}: {e}"
                    append_raw({"method": method, "rule": rule, "k": k, "n": n, "d": d,
                                "grid": grid, "cell_value": cell_value, "rep": rep,
                                "time_s": round(elapsed, 4) if elapsed == elapsed else "",
                                "mem_peak_mb": round(mem, 3) if mem == mem else "",
                                "threads": 1, "status": status})
                    b.log(f"  rep {rep:2d} {method:10s} {rule:6s} k={k:2d} {elapsed:.4f}s mem={mem:.2f}MB [{status}]")
    summarize()


if __name__ == "__main__":
    main()
