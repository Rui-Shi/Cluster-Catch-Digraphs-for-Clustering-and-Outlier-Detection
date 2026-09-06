"""
81_wp4_baselines.py -- WP4 (AE.3, R1.3, R1.5, R3.3, R6.5): eight modern /
task-aligned competitors on the sixteen Section 6 real-data sets.

Everything this script does is declared in advance in
revision_experiments/tr1/WP4_PROTOCOL.md. Read that first; this file is the
implementation, not the specification.

Methods and cells (24 per data set, 384 in total):

  ECOD        pyod defaults, deterministic, 1 fit
  COPOD       pyod defaults, deterministic, 1 fit
  HDBSCAN     hdbscan defaults; score = outlier_scores_ (GLOSH); + native noise labels
  OPTICS      sklearn defaults; score = reachability_ with non-finite entries
              replaced by 1.01 x max finite; + native noise labels
  MutualKNN   k in {5,10,15,20,30}; score = k - m_i; + native labels (m_i == 0)
  SNN         k in {5,10,15,20,30}; score = -(SNN density), Ertoz-Steinbach-Kumar
  DIF         pyod defaults, seeds 1..5
  LUNAR       pyod defaults, seeds 1..5

Order within a data set is cheapest-first, and data sets run smallest-n first,
so partial results are useful early.

LABEL POLARITY -- two opposite conventions, never mixed:
  * results/tr1/wp4/data/<dataset>.csv column `label`: 1 = REGULAR, 0 = OUTLIER
    (this repo's / RealData_Collection.R's convention). This script only reads
    it, to compute n and the true contamination rate for the log.
  * every `is_outlier` column this script writes: 1 = OUTLIER, 0 = regular
    (PyOD / sklearn noise-label convention). Native threshold-free variants only.

SCORE POLARITY: higher = more outlying, for every method, without exception.
Densities are negated.

THRESHOLDS: none are applied here. This script emits raw scores; T1
(contamination = 0.1) and T2 (contamination = true rate) labels, and all
TPR/TNR/BA/F2 metrics, are derived downstream in R by the same evaluate()
that scores the nine methods already in the study.

Checkpointing: one line is APPENDED to results/tr1/wp4/fit_log.csv the moment
a cell finishes -- nothing is buffered until the end, because the volume this
repo lives on drops off the bus under sustained write load. A cell already
logged "ok" is skipped on restart, so the driver is resumable and is meant to
be called repeatedly under a short wall-clock budget.

Usage (from anywhere; paths are resolved from this file's location):
    .venv/python.exe revision_experiments/tr1/81_wp4_baselines.py [options]

    --budget SECONDS   stop starting new cells once this much wall clock has
                       elapsed (default 400), so a single invocation fits
                       comfortably inside a 10-minute cap
    --dataset NAME     restrict to one data set (repeatable)
    --method NAME      restrict to one method (repeatable)
    --status           print remaining/completed counts and exit
    --data-dir PATH    override the data directory this script reads
                       <dataset>.csv / manifest.csv from. Default (omit this
                       flag): results/tr1/wp4/data, i.e. CURRENT BEHAVIOUR IS
                       UNCHANGED. Added for WP5 (tr1/84b_wp5_metrics.R),
                       which points this at results/tr1/wp5/data so the same
                       eight competitors, under the identical settings and
                       thresholding rule declared in WP4_PROTOCOL.md, run
                       against the four real data sets above d = 21.
    --out-dir PATH     override the output directory this script writes
                       scores/, fit_log.csv, fit_errors.log and versions.txt
                       into. Default (omit this flag): results/tr1/wp4, i.e.
                       CURRENT BEHAVIOUR IS UNCHANGED. WP5 points this at
                       results/tr1/wp5.
"""

import argparse
import csv
import importlib.metadata as importlib_metadata
import random
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WP4 = HERE.parent / "results" / "tr1" / "wp4"
DATA_DIR = WP4 / "data"
SCORES_DIR = WP4 / "scores"
FIT_LOG_PATH = WP4 / "fit_log.csv"
ERROR_LOG_PATH = WP4 / "fit_errors.log"
VERSIONS_PATH = WP4 / "versions.txt"

K_GRID = [5, 10, 15, 20, 30]
SEEDS = [1, 2, 3, 4, 5]
NA = "NA"
MAX_RETRIES = 2

FIT_LOG_HEADER = [
    "dataset", "method", "k", "seed", "n", "d",
    "fit_seconds", "status", "error", "note",
]


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# cells
# ---------------------------------------------------------------------------

def cells_for_dataset(dataset):
    """Cheapest-first within a data set."""
    out = [
        (dataset, "ECOD", NA, NA),
        (dataset, "COPOD", NA, NA),
        (dataset, "HDBSCAN", NA, NA),
        (dataset, "OPTICS", NA, NA),
    ]
    for k in K_GRID:
        out.append((dataset, "MutualKNN", str(k), NA))
    for k in K_GRID:
        out.append((dataset, "SNN", str(k), NA))
    for s in SEEDS:
        out.append((dataset, "DIF", NA, str(s)))
    for s in SEEDS:
        out.append((dataset, "LUNAR", NA, str(s)))
    return out


def cell_stem(dataset, method, k, seed):
    stem = f"{dataset}_{method}"
    if k != NA:
        stem += f"_k{k}"
    if seed != NA:
        stem += f"_seed{seed}"
    return stem


# ---------------------------------------------------------------------------
# fit log -- append-only, one line per finished cell
# ---------------------------------------------------------------------------

def read_log_state():
    """(cells logged ok, {cell: number of logged failures})."""
    done = set()
    fails = {}
    if not FIT_LOG_PATH.exists():
        return done, fails
    with open(FIT_LOG_PATH, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (row["dataset"], row["method"], row["k"], row["seed"])
            if row.get("status") == "ok":
                done.add(key)
            else:
                fails[key] = fails.get(key, 0) + 1
    return done, fails


def read_done_cells():
    return read_log_state()[0]


def append_fit_log(dataset, method, k, seed, n, d, fit_seconds, status, error, note):
    new = not FIT_LOG_PATH.exists()
    with open(FIT_LOG_PATH, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if new:
            w.writerow(FIT_LOG_HEADER)
        w.writerow([
            dataset, method, k, seed, n, d,
            "" if fit_seconds is None else f"{fit_seconds:.3f}",
            status, error, note,
        ])
        f.flush()


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def load_dataset(name):
    df = pd.read_csv(DATA_DIR / f"{name}.csv")
    y = df["label"].to_numpy(dtype=int)          # 1 = regular, 0 = outlier
    X = df.drop(columns=["label"]).to_numpy(dtype=np.float64)
    return X, y


def save_scores(stem, scores):
    pd.DataFrame({"score": np.asarray(scores, dtype=float).ravel()}).to_csv(
        SCORES_DIR / f"{stem}.csv", index=False
    )


def save_native_labels(stem, is_outlier):
    """PyOD polarity: 1 = outlier, 0 = regular."""
    pd.DataFrame({"is_outlier": np.asarray(is_outlier, dtype=int).ravel()}).to_csv(
        SCORES_DIR / f"{stem}_labels.csv", index=False
    )


def seed_everything(seed):
    import torch
    s = int(seed)
    np.random.seed(s)
    random.seed(s)
    torch.manual_seed(s)


# ---------------------------------------------------------------------------
# neighbour-based detectors (own code)
# ---------------------------------------------------------------------------

def knn_index(X, k):
    """Exact Euclidean kNN indices, self excluded, shape (n, k).

    Ties are broken by index order, which is sklearn's own behaviour. The
    query asks for k+1 neighbours and drops self; if self is somehow absent
    (should not happen -- the loader deduplicates rows) the farthest of the
    k+1 is dropped instead, which keeps the row length at exactly k.
    """
    from sklearn.neighbors import NearestNeighbors

    n = X.shape[0]
    nn = NearestNeighbors(n_neighbors=min(k + 1, n)).fit(X)
    _, idx = nn.kneighbors(X)
    out = np.empty((n, k), dtype=np.int64)
    for i in range(n):
        row = idx[i]
        row = row[row != i]
        out[i] = row[:k]
    return out


def knn_adjacency(idx, n):
    """Sparse boolean A with A[i, j] = 1 iff j in kNN(i)."""
    from scipy import sparse

    k = idx.shape[1]
    rows = np.repeat(np.arange(n, dtype=np.int64), k)
    cols = idx.ravel()
    data = np.ones(rows.shape[0], dtype=np.int32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def mutual_knn_scores(X, k):
    """score = k - m_i, m_i = #{j : j in kNN(i) and i in kNN(j)}."""
    n = X.shape[0]
    A = knn_adjacency(knn_index(X, k), n)
    m = np.asarray(A.multiply(A.T).sum(axis=1)).ravel()
    return (k - m).astype(float), m.astype(int)


def snn_scores(X, k):
    """Ertoz-Steinbach-Kumar SNN density; score = -density."""
    n = X.shape[0]
    A = knn_adjacency(knn_index(X, k), n)
    S = (A @ A.T).tocsr()               # S[i, j] = |kNN(i) cap kNN(j)|
    density = np.asarray(A.multiply(S).sum(axis=1)).ravel()
    return -density.astype(float)


# ---------------------------------------------------------------------------
# one cell
# ---------------------------------------------------------------------------

def run_cell(dataset, method, k, seed, X, y):
    """Fit one cell, write its outputs, return (fit_seconds, note)."""
    stem = cell_stem(dataset, method, k, seed)
    t0 = time.time()
    note = ""

    if method in ("ECOD", "COPOD", "DIF", "LUNAR"):
        if method == "ECOD":
            from pyod.models.ecod import ECOD
            model = ECOD()
        elif method == "COPOD":
            from pyod.models.copod import COPOD
            model = COPOD()
        elif method == "DIF":
            from pyod.models.dif import DIF
            seed_everything(seed)
            model = DIF(random_state=int(seed))
        else:
            from pyod.models.lunar import LUNAR
            seed_everything(seed)
            model = LUNAR(random_state=int(seed), verbose=0)
        model.fit(X)
        scores = np.asarray(model.decision_scores_, dtype=float).ravel()
        elapsed = time.time() - t0
        if scores.shape[0] != X.shape[0]:
            raise RuntimeError(f"decision_scores_ length {scores.shape[0]} != n {X.shape[0]}")
        save_scores(stem, scores)

    elif method == "HDBSCAN":
        import hdbscan
        model = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=None)
        model.fit(X)
        scores = np.asarray(model.outlier_scores_, dtype=float).ravel()
        elapsed = time.time() - t0
        n_nonfinite = int((~np.isfinite(scores)).sum())
        n_noise = int((model.labels_ == -1).sum())
        n_clusters = int(len(set(model.labels_.tolist()) - {-1}))
        note = f"n_clusters={n_clusters};n_noise={n_noise};n_nonfinite_glosh={n_nonfinite}"
        save_scores(stem, scores)
        save_native_labels(stem, (model.labels_ == -1).astype(int))

    elif method == "OPTICS":
        from sklearn.cluster import OPTICS
        model = OPTICS(min_samples=5, xi=0.05)
        model.fit(X)
        reach = np.asarray(model.reachability_, dtype=float).ravel()
        elapsed = time.time() - t0
        finite = np.isfinite(reach)
        n_inf = int((~finite).sum())
        if finite.any():
            fill = float(np.max(reach[finite])) * 1.01
        else:
            fill = 1.0
        scores = reach.copy()
        scores[~finite] = fill
        n_noise = int((model.labels_ == -1).sum())
        n_clusters = int(len(set(model.labels_.tolist()) - {-1}))
        note = (f"n_clusters={n_clusters};n_noise={n_noise};"
                f"n_inf_reachability={n_inf};inf_fill={fill:.6g}")
        save_scores(stem, scores)
        save_native_labels(stem, (model.labels_ == -1).astype(int))

    elif method == "MutualKNN":
        scores, m = mutual_knn_scores(X, int(k))
        elapsed = time.time() - t0
        note = f"n_zero_mutual_degree={int((m == 0).sum())}"
        save_scores(stem, scores)
        save_native_labels(stem, (m == 0).astype(int))

    elif method == "SNN":
        scores = snn_scores(X, int(k))
        elapsed = time.time() - t0
        save_scores(stem, scores)

    else:
        raise ValueError(f"unknown method {method}")

    n_nonfinite = int((~np.isfinite(np.asarray(scores, dtype=float))).sum())
    if n_nonfinite:
        note = (note + ";" if note else "") + f"n_nonfinite_score={n_nonfinite}"
    return elapsed, note


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def record_versions():
    import sklearn
    import scipy
    import torch
    import pyod

    lines = [
        f"python      {sys.version}",
        f"numpy       {np.__version__}",
        f"pandas      {pd.__version__}",
        f"scipy       {scipy.__version__}",
        f"scikit-learn {sklearn.__version__}",
        f"torch       {torch.__version__}",
        f"pyod        {pyod.__version__}",
        f"hdbscan     {importlib_metadata.version('hdbscan')}",
    ]
    VERSIONS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for ln in lines:
        log("  " + ln)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=float, default=400.0)
    ap.add_argument("--dataset", action="append", default=None)
    ap.add_argument("--method", action="append", default=None)
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--data-dir", type=str, default=None,
                     help="override data dir (default: results/tr1/wp4/data)")
    ap.add_argument("--out-dir", type=str, default=None,
                     help="override output dir for scores/fit_log/etc (default: results/tr1/wp4)")
    args = ap.parse_args()

    # Overriding these module-level constants here, before anything below
    # reads them, is what lets WP5 point this driver at its own data/output
    # folders without touching any method, hyperparameter, seed, threshold
    # rule, or file-naming convention inside the script. Omitting both flags
    # leaves every path exactly as it was (results/tr1/wp4/...).
    global DATA_DIR, SCORES_DIR, FIT_LOG_PATH, ERROR_LOG_PATH, VERSIONS_PATH
    if args.data_dir is not None:
        DATA_DIR = Path(args.data_dir)
    if args.out_dir is not None:
        out_dir = Path(args.out_dir)
        SCORES_DIR = out_dir / "scores"
        FIT_LOG_PATH = out_dir / "fit_log.csv"
        ERROR_LOG_PATH = out_dir / "fit_errors.log"
        VERSIONS_PATH = out_dir / "versions.txt"

    SCORES_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(DATA_DIR / "manifest.csv").sort_values("n")
    datasets = list(manifest["dataset"])
    if args.dataset:
        datasets = [d for d in datasets if d in set(args.dataset)]

    all_cells = []
    for ds in datasets:
        for c in cells_for_dataset(ds):
            if args.method and c[1] not in set(args.method):
                continue
            all_cells.append(c)

    done, fails = read_log_state()
    # A cell that has already failed MAX_RETRIES times is not retried again:
    # it stays in the log as a failure and is reported as one.
    todo = [c for c in all_cells if c not in done and fails.get(c, 0) < MAX_RETRIES]
    dead = [c for c in all_cells if c not in done and fails.get(c, 0) >= MAX_RETRIES]

    if args.status:
        log(f"{len(all_cells)} cells in scope, {len(done & set(all_cells))} done, "
            f"{len(todo)} remaining, {len(dead)} given up on after {MAX_RETRIES} failures")
        if todo:
            by_method = {}
            for c in todo:
                by_method[c[1]] = by_method.get(c[1], 0) + 1
            log("  remaining by method: " + ", ".join(f"{m}={v}" for m, v in sorted(by_method.items())))
            log(f"  next: {todo[0]}")
        return 0

    log(f"WP4 baselines: {len(all_cells)} cells in scope, {len(todo)} remaining, budget {args.budget:.0f}s")
    if not todo:
        log("nothing to do")
        return 0
    record_versions()

    t_start = time.time()
    cache_name, cache_X, cache_y = None, None, None
    n_ok = n_fail = 0

    for (ds, method, k, seed) in todo:
        if time.time() - t_start > args.budget:
            log(f"budget reached; stopping cleanly with {len(todo) - n_ok - n_fail} cells left in this scope")
            break

        if cache_name != ds:
            cache_X, cache_y = load_dataset(ds)
            cache_name = ds
        X, y = cache_X, cache_y
        n, d = X.shape

        try:
            elapsed, note = run_cell(ds, method, k, seed, X, y)
            append_fit_log(ds, method, k, seed, n, d, elapsed, "ok", "", note)
            n_ok += 1
            log(f"ok    {ds:13s} {method:10s} k={k:3s} seed={seed:3s} n={n:5d} d={d:3d} t={elapsed:8.2f}s {note}")
        except Exception as e:  # noqa: BLE001 -- a failed cell must not abort the sweep
            tb = traceback.format_exc()
            with open(ERROR_LOG_PATH, "a", encoding="utf-8") as f:
                f.write(f"\n=== {ds} {method} k={k} seed={seed} ===\n{tb}\n")
            msg = f"{type(e).__name__}: {e}".replace("\n", " ")[:300]
            append_fit_log(ds, method, k, seed, n, d, None, "failed", msg, "")
            n_fail += 1
            log(f"FAIL  {ds:13s} {method:10s} k={k:3s} seed={seed:3s}: {msg}")

    log(f"\nthis invocation: {n_ok} ok, {n_fail} failed, {time.time() - t_start:.0f}s wall clock")
    remaining = len([c for c in all_cells if c not in read_done_cells()])
    log(f"remaining in scope: {remaining}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
