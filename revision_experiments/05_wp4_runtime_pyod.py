"""05_wp4_runtime_pyod.py -- WP4: PyOD runtime grid (task T4).

Times ECOD, LUNAR, and AutoEncoder *fits* (fit only, features only -- labels
are never passed to .fit()) on the datasets exported by 04_wp4_runtime.R to
results/wp4_data2/<grid>_<cell_value>_rep1.csv, 10 reps each by default.

REDESIGN (2026-07-19, "runtime2"): 04_wp4_runtime.R's grids changed --
Grid 1 ("n"): n in {100,250,500,1000,2000} at d=10 (was d=20, n to 4000);
Grid 2 ("d"): d in {5,10,50,100,500} at n=500 (was n=1000, d to 1000). This
script has NO hardcoded grid of its own -- it only discovers whatever
*_rep1.csv files exist under DATA_DIR and times them -- so "updating its
grids" means pointing DATA_DIR at the NEW export directory (wp4_data2, not
the original wp4_data) and writing to NEW raw/aggregate files
(wp4_runtime2_pyod_raw.csv / wp4_runtime2_pyod.csv). wp4_data2 exists
because the old and new grids' cell_value labels collide (e.g. both have a
"n_250" cell, at d=20 vs d=10) -- reusing the old export directory would
silently time the WRONG dataset for several new-grid cells. The ORIGINAL
wp4_data/, wp4_runtime_pyod_raw.csv, and wp4_runtime_pyod.csv are untouched
by this file.

Comparability notes (documented per task brief):
  - Single-threaded: OMP/MKL/OPENBLAS/NUMEXPR thread env vars are pinned to 1
    BEFORE importing numpy/torch, and torch.set_num_threads(1) /
    torch.set_num_interop_threads(1) are applied. The R methods in
    04_wp4_runtime.R are likewise single-threaded (no doParallel in the timed
    path; isotree pinned to nthreads=1), so wall-clock numbers are
    single-worker comparable across the two scripts.
  - Data asymmetry: R times each rep on a FRESH dataset draw
    (seed = 1000*rep + cell_index), while this script times all reps on the
    single exported rep-1 file per cell (identical data across its reps; reps
    vary only model seeds and system noise). Documented in
    04_wp4_runtime.R's header as well.
  - Timing covers model.fit(X) only: model construction, data loading, and
    score post-processing are excluded.

Seeding: per rep, numpy/random/torch global RNGs are seeded with the rep
index, and rep is also passed as random_state to LUNAR/AutoEncoder (ECOD is
deterministic; its reps only measure timing noise).

Checkpointing: each (grid, cell_value, method, rep) appends immediately to
results/wp4_runtime_pyod_raw.csv and is skipped when already present, so
the script resumes cleanly and never duplicates rows. Failed fits are
recorded with status FAILED:<error> and do not abort the sweep.

Aggregation (end of every invocation): mean/sd fit seconds per
(grid, cell_value, method) over OK rows -> results/wp4_runtime_pyod.csv
(all grids in one file; the `grid` column separates them, micro included
if present in the raw file).

Usage (from the CLONE repo root, using the pinned venv):
  & "revision_experiments\\.venv\\Scripts\\python.exe" `
      "revision_experiments/05_wp4_runtime_pyod.py" [--reps 10] [--include-micro]

Files whose names start with "micro_" (written by 04's hidden micro mode)
are excluded unless --include-micro is given.
"""

import os

# Must precede numpy/torch imports for single-threaded comparability.
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import random
import re
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import torch

torch.set_num_threads(1)
try:
    torch.set_num_interop_threads(1)
except RuntimeError:
    # Can only be set once per process, before parallel work begins; if the
    # import graph already initialized interop threading, proceed (intra-op
    # threading above is the setting that matters for these small models).
    pass

from pyod.models.auto_encoder import AutoEncoder
from pyod.models.ecod import ECOD
from pyod.models.lunar import LUNAR

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "results" / "tr2" / "wp4_data2"
RAW_PATH = HERE / "results" / "tr2" / "wp4_runtime2_pyod_raw.csv"
AGG_PATH = HERE / "results" / "tr2" / "wp4_runtime2_pyod.csv"
ERROR_LOG_PATH = HERE / "results" / "tr2" / "wp4_pyod_errors.log"

METHODS = ["ECOD", "LUNAR", "AutoEncoder"]
FILE_RE = re.compile(r"^(n|d|micro)_(.+)_rep1\.csv$")

RAW_COLS = ["grid", "cell_value", "n", "d", "method", "rep", "seed",
            "fit_seconds", "status", "timestamp"]


def log(msg):
    print(msg, flush=True)


def seed_everything(seed):
    s = int(seed)
    np.random.seed(s)
    random.seed(s)
    torch.manual_seed(s)


def build_model(method, seed):
    # pyod defaults throughout (contamination only affects labels_, not the
    # fit or decision_scores_ -- verified for these models in 07_wp6_pyod.py).
    if method == "ECOD":
        return ECOD()
    if method == "LUNAR":
        return LUNAR(random_state=seed, verbose=0)
    if method == "AutoEncoder":
        return AutoEncoder(random_state=seed, verbose=0)
    raise ValueError(f"unknown method {method}")


def discover_files(include_micro):
    out = []
    for p in sorted(DATA_DIR.glob("*_rep1.csv")):
        m = FILE_RE.match(p.name)
        if not m:
            log(f"[warn] unrecognized file name skipped: {p.name}")
            continue
        grid, cell_value = m.group(1), m.group(2)
        if grid == "micro" and not include_micro:
            continue
        out.append((grid, cell_value, p))
    return out


def load_existing_keys():
    if not RAW_PATH.exists():
        return set()
    try:
        df = pd.read_csv(RAW_PATH)
    except Exception:
        return set()
    return {
        (str(r["grid"]), str(r["cell_value"]), str(r["method"]), int(r["rep"]))
        for _, r in df.iterrows()
    }


def append_row(row):
    df = pd.DataFrame([row], columns=RAW_COLS)
    header = not RAW_PATH.exists()
    df.to_csv(RAW_PATH, mode="a", header=header, index=False)


def aggregate():
    if not RAW_PATH.exists():
        log("[aggregate] no raw file; nothing to aggregate")
        return
    df = pd.read_csv(RAW_PATH)
    if df.empty:
        return
    df["cell_value"] = df["cell_value"].astype(str)
    ok = df[df["status"] == "OK"]
    g_ok = (
        ok.groupby(["grid", "cell_value", "n", "d", "method"], as_index=False)
        .agg(n_reps_ok=("fit_seconds", "size"),
             mean_fit_seconds=("fit_seconds", "mean"),
             sd_fit_seconds=("fit_seconds", "std"))
    )
    g_status = (
        df.groupby(["grid", "cell_value", "method"], as_index=False)["status"]
        .agg(lambda s: ";".join(sorted(set(s))))
        .rename(columns={"status": "statuses"})
    )
    out = g_ok.merge(g_status, on=["grid", "cell_value", "method"], how="outer")
    out["_ord"] = pd.to_numeric(out["cell_value"], errors="coerce")
    out = out.sort_values(["grid", "_ord", "cell_value", "method"]).drop(columns="_ord")
    out["mean_fit_seconds"] = out["mean_fit_seconds"].round(4)
    out["sd_fit_seconds"] = out["sd_fit_seconds"].round(4)
    out.to_csv(AGG_PATH, index=False)
    log(f"[aggregate] wrote {AGG_PATH} ({len(out)} rows)")


def main():
    ap = argparse.ArgumentParser(description="WP4 PyOD runtime grid (fit-only timing)")
    ap.add_argument("--reps", type=int, default=10, help="timing reps per (cell, method); default 10")
    ap.add_argument("--include-micro", action="store_true",
                    help="also time the micro_* validation exports (excluded by default)")
    args = ap.parse_args()

    log(f"numpy {np.__version__}, pandas {pd.__version__}, torch {torch.__version__}")
    import pyod
    log(f"pyod {pyod.__version__}; torch threads = {torch.get_num_threads()}")

    files = discover_files(args.include_micro)
    if not files:
        log(f"No wp4_data exports found in {DATA_DIR}. Run 04_wp4_runtime.R first.")
        return
    log(f"{len(files)} dataset files: {[p.name for _, _, p in files]}")

    existing = load_existing_keys()
    n_done = n_skip = n_fail = 0

    for grid, cell_value, path in files:
        df = pd.read_csv(path)
        X = df.drop(columns=["label"]).values.astype(np.float64)
        n, d = X.shape
        log(f"\n==== {path.name} (grid={grid}, cell={cell_value}, n={n}, d={d}) ====")

        for method in METHODS:
            for rep in range(1, args.reps + 1):
                key = (grid, cell_value, method, rep)
                if key in existing:
                    log(f"  {method:12s} rep {rep}/{args.reps} [checkpoint skip]")
                    n_skip += 1
                    continue
                seed_everything(rep)
                try:
                    model = build_model(method, rep)   # construction untimed
                    t0 = time.perf_counter()
                    model.fit(X)                        # fit only, features only
                    elapsed = time.perf_counter() - t0
                    status = "OK"
                    n_done += 1
                except Exception as e:  # noqa: BLE001 -- never abort the sweep
                    elapsed = float("nan")
                    status = f"FAILED: {type(e).__name__}: {e}"
                    with open(ERROR_LOG_PATH, "a", encoding="utf-8") as f:
                        f.write(f"\n=== {grid} {cell_value} {method} rep{rep} ===\n"
                                f"{traceback.format_exc()}\n")
                    n_fail += 1
                append_row({
                    "grid": grid, "cell_value": cell_value, "n": n, "d": d,
                    "method": method, "rep": rep, "seed": rep,
                    "fit_seconds": round(elapsed, 4) if elapsed == elapsed else "",
                    "status": status,
                    "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
                existing.add(key)
                shown = f"{elapsed:.3f}s" if elapsed == elapsed else "NA"
                log(f"  {method:12s} rep {rep}/{args.reps} fit={shown} [{status if status == 'OK' else status[:80]}]")

    log(f"\nDone: {n_done} timed, {n_skip} checkpoint-skipped, {n_fail} failed.")
    if n_fail:
        log(f"Tracebacks in {ERROR_LOG_PATH}")
    aggregate()


if __name__ == "__main__":
    main()
