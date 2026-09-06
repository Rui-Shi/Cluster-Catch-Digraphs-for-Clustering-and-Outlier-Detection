"""91b_wp6_runtime_py.py -- WP6: Python-side runtime/memory grid for the 8
WP4 competitors (AE.5, R3.10 runtime/memory/scalability half, R5.5, R5.6
scalability angle). See tr1/WP6_PROTOCOL.md, declared before this file was
written -- that document is the specification; this is the implementation.

Times ECOD, COPOD, HDBSCAN, OPTICS, MutualKNN, SNN, DIF, LUNAR -- the same
eight competitors and settings as tr1/81_wp4_baselines.py (read-only; this
script does not import or edit it, it reuses its documented settings) -- on
the per-cell datasets 91_wp6_runtime.R exports to results/tr1/wp6/data/
(or results/tr1/wp6/smoke/data/ for --smoke).

Grid (fixed by 91_wp6_runtime.R's CELLS table, discovered from file names,
not hardcoded here): n in {100,250,500,1000,2000} at d=10 ("n" files),
d in {5,10,50,100} at n=500 ("d" files). MutualKNN/SNN: one timed fit per
k in {5,10,15,20,30}. DIF/LUNAR: one timed fit per seed in {1,...,5}.
ECOD/COPOD/HDBSCAN/OPTICS: one fit each (deterministic; reps vary timing
noise only). All timed `--reps` times per (cell, method, k-or-seed).

Comparability with the R side (91_wp6_runtime.R):
  - Single-threaded: OMP/MKL/OPENBLAS/NUMEXPR/VECLIB env vars pinned to 1
    BEFORE numpy/torch are imported (thread pools are fixed at import time),
    plus torch.set_num_threads(1) / set_num_interop_threads(1). hdbscan and
    sklearn's OPTICS have no separate thread knob beyond BLAS/OpenMP, which
    the env vars cover.
  - Timed region: construction (model = Method(...)) is NOT timed; only the
    fit call (`model.fit(X)`) is, via time.perf_counter(), matching
    81_wp4_baselines.py's "timing covers model.fit(X) only" convention and
    05_wp4_runtime_pyod.py's "fit only, features only" rule.
  - Data asymmetry (declared, not hidden): the R side draws a fresh dataset
    per rep; this script times all its reps against the SAME rep-1 export
    per cell, so its reps vary only model-construction seed / system noise,
    identical to 04_wp4_runtime.R / 05_wp4_runtime_pyod.py's own asymmetry.

Memory: tracemalloc peak (not psutil RSS delta -- see WP6_PROTOCOL.md #5 for
why: RSS delta is contaminated by the interpreter's own one-time import
footprint and by OS allocator behaviour). tracemalloc.start() immediately
before model.fit(X), tracemalloc.get_traced_memory() -> (current, peak)
immediately after, tracemalloc.stop(). Known limitation, stated in advance:
PyTorch's own C++-side allocations inside DIF/LUNAR are invisible to
tracemalloc, so those two methods' mem_peak_mb under-counts their true
footprint.

Checkpointing: keys-only done file (results/tr1/wp6/91b_wp6_runtime_py_done.csv,
or the smoke-directory equivalent for --smoke), same reason as the R side's
done file: table_load_s-equivalent columns are legitimately absent/NA for
every Python method, and checkpointing against a payload-bearing raw file
would make has_result()-style "any NA is incomplete" logic misfire. Rows are
appended to the raw CSV as soon as each (cell, method, variant, rep) finishes.

Usage (from the CLONE repo root, using the pinned venv):
  revision_experiments/.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --smoke
  revision_experiments/.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py --idle-check
  revision_experiments/.venv/python.exe revision_experiments/tr1/91b_wp6_runtime_py.py [--reps 10] [--force]

NOT RUN FOR REAL by the session that wrote this file -- needs an idle
machine; only --smoke was exercised. See tr1/WP6_PROTOCOL.md #10.
"""

import os

# Must precede numpy/torch imports for single-threaded comparability
# (matches 05_wp4_runtime_pyod.py's ordering).
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import random
import re
import subprocess
import sys
import time
import traceback
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd
import torch

torch.set_num_threads(1)
try:
    torch.set_num_interop_threads(1)
except RuntimeError:
    pass   # can only be set once per process; harmless if already fixed

HERE = Path(__file__).resolve().parent.parent   # revision_experiments/
WP6 = HERE / "results" / "tr1" / "wp6"

K_GRID = [5, 10, 15, 20, 30]
SEEDS = [1, 2, 3, 4, 5]
METHODS_SIMPLE = ["ECOD", "COPOD", "HDBSCAN", "OPTICS"]   # one fit each
METHODS_K = ["MutualKNN", "SNN"]                          # one fit per k
METHODS_SEED = ["DIF", "LUNAR"]                           # one fit per seed

RAW_COLS = ["method", "variant", "n", "d", "cell_value", "rep", "seed",
            "time_s", "mem_peak_mb", "table_load_s", "threads", "status"]
DONE_COLS = ["method", "variant", "n", "d", "cell_value", "rep"]

FILE_RE = re.compile(r"^(n|d)_(.+)_rep1\.csv$")


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description="WP6 Python-side runtime/memory grid")
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--smoke", action="store_true",
                    help="one cell (n=100,d=10), 1 rep, all methods, "
                         "results/tr1/wp6/smoke/ -- run twice to see the "
                         "checkpoint skip fire")
    ap.add_argument("--idle-check", action="store_true",
                    help="list other Rscript/python processes and exit; runs nothing")
    ap.add_argument("--force", action="store_true",
                    help="skip the idle-machine check before a real run")
    return ap.parse_args()


# ---------------------------------------------------------------------------
# Idle-machine check (WP6_PROTOCOL.md #6). psutil is not in this venv
# (verified: ModuleNotFoundError), so this uses the documented tasklist
# fallback, matching 91_wp6_runtime.R's own approach.
# ---------------------------------------------------------------------------
WATCH_IMAGES = {"Rscript.exe", "Rterm.exe", "R.exe", "python.exe", "pythonw.exe"}


def idle_check():
    self_pid = os.getpid()
    try:
        out = subprocess.run(["tasklist", "/FO", "CSV", "/NH"],
                             capture_output=True, text=True, check=False).stdout
    except Exception as e:  # noqa: BLE001
        log(f"[idle-check] tasklist unavailable ({e}) -- cannot verify; proceeding as if idle.")
        return True
    hits = []
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip('"') for p in line.split('","')]
        if len(parts) < 2:
            continue
        image = parts[0].strip('"')
        try:
            pid = int(parts[1])
        except ValueError:
            continue
        if pid == self_pid:
            continue
        if image in WATCH_IMAGES:
            hits.append((image, pid))
    if not hits:
        log("[idle-check] no other Rscript/Rterm/R/python/pythonw processes found. Machine looks idle.")
        return True
    log("[idle-check] OTHER PROCESSES FOUND (machine is NOT idle):")
    for image, pid in hits:
        log(f"    {image}  PID={pid}")
    return False


# ---------------------------------------------------------------------------
# data discovery
# ---------------------------------------------------------------------------
def discover_files(data_dir):
    out = []
    for p in sorted(data_dir.glob("*_rep1.csv")):
        m = FILE_RE.match(p.name)
        if not m:
            log(f"[warn] unrecognized file name skipped: {p.name}")
            continue
        grid, cell_value = m.group(1), m.group(2)
        out.append((grid, cell_value, p))
    return out


def load_cell(path):
    df = pd.read_csv(path)
    X = df.drop(columns=["label"]).to_numpy(dtype=np.float64)
    return X


# ---------------------------------------------------------------------------
# checkpoint files
# ---------------------------------------------------------------------------
def read_done_keys(done_path):
    if not done_path.exists():
        return set()
    try:
        # keep_default_na=False/na_filter=False: an empty `variant` field
        # (the 4 one-fit methods use "") must round-trip as "", not NaN --
        # pandas' default NA-sniffing otherwise turns "" into NaN even under
        # dtype=str, which silently breaks the key match against a freshly
        # built key of ("", ...) and made every ECOD/COPOD/HDBSCAN/OPTICS
        # smoke rep look undone forever (caught by the second --smoke run).
        df = pd.read_csv(done_path, dtype=str, keep_default_na=False, na_filter=False)
    except Exception:
        return set()
    return {tuple(r[c] for c in DONE_COLS) for _, r in df.iterrows()}


def append_done(done_path, key):
    row = dict(zip(DONE_COLS, key))
    df = pd.DataFrame([row], columns=DONE_COLS)
    header = not done_path.exists()
    done_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(done_path, mode="a", header=header, index=False)


def append_raw(raw_path, row):
    df = pd.DataFrame([row], columns=RAW_COLS)
    header = not raw_path.exists()
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(raw_path, mode="a", header=header, index=False)


# ---------------------------------------------------------------------------
# neighbour-based detectors (own code -- verbatim settings from
# 81_wp4_baselines.py, reproduced here since that file is off-limits to
# import from / edit for this work package's scripts)
# ---------------------------------------------------------------------------
def knn_index(X, k):
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
    from scipy import sparse

    k = idx.shape[1]
    rows = np.repeat(np.arange(n, dtype=np.int64), k)
    cols = idx.ravel()
    data = np.ones(rows.shape[0], dtype=np.int32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def fit_mutual_knn(X, k):
    n = X.shape[0]
    A = knn_adjacency(knn_index(X, k), n)
    m = np.asarray(A.multiply(A.T).sum(axis=1)).ravel()
    return (k - m).astype(float)


def fit_snn(X, k):
    n = X.shape[0]
    A = knn_adjacency(knn_index(X, k), n)
    S = (A @ A.T).tocsr()
    density = np.asarray(A.multiply(S).sum(axis=1)).ravel()
    return -density.astype(float)


def seed_everything(seed):
    s = int(seed)
    np.random.seed(s)
    random.seed(s)
    torch.manual_seed(s)


# ---------------------------------------------------------------------------
# one timed fit
# ---------------------------------------------------------------------------
def timed_fit(method, X, k=None, seed=None):
    """Returns (elapsed_s, mem_peak_mb, status)."""
    if method == "ECOD":
        from pyod.models.ecod import ECOD
        model = ECOD()
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    elif method == "COPOD":
        from pyod.models.copod import COPOD
        model = COPOD()
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    elif method == "HDBSCAN":
        import hdbscan
        model = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=None)
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    elif method == "OPTICS":
        from sklearn.cluster import OPTICS
        model = OPTICS(min_samples=5, xi=0.05)
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    elif method == "MutualKNN":
        tracemalloc.start()
        t0 = time.perf_counter()
        fit_mutual_knn(X, int(k))
        elapsed = time.perf_counter() - t0
    elif method == "SNN":
        tracemalloc.start()
        t0 = time.perf_counter()
        fit_snn(X, int(k))
        elapsed = time.perf_counter() - t0
    elif method == "DIF":
        from pyod.models.dif import DIF
        seed_everything(seed)
        model = DIF(random_state=int(seed))
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    elif method == "LUNAR":
        from pyod.models.lunar import LUNAR
        seed_everything(seed)
        model = LUNAR(random_state=int(seed), verbose=0)
        tracemalloc.start()
        t0 = time.perf_counter()
        model.fit(X)
        elapsed = time.perf_counter() - t0
    else:
        raise ValueError(f"unknown method {method}")

    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    mem_peak_mb = peak / (1024.0 * 1024.0)
    return elapsed, mem_peak_mb


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def cells_and_methods():
    """Yields (method, variant_label, k_or_None, seed_or_None)."""
    out = []
    for m in METHODS_SIMPLE:
        out.append((m, "", None, None))
    for m in METHODS_K:
        for k in K_GRID:
            out.append((m, f"k{k}", k, None))
    for m in METHODS_SEED:
        for s in SEEDS:
            out.append((m, f"seed{s}", None, s))
    return out


def main():
    args = parse_args()

    if args.idle_check:
        ok = idle_check()
        sys.exit(0 if ok else 1)

    if args.smoke:
        data_dir = WP6 / "smoke" / "data"
        raw_path = WP6 / "smoke" / "91b_wp6_runtime_py_raw.csv"
        done_path = WP6 / "smoke" / "91b_wp6_runtime_py_done.csv"
        reps = 1
        # smoke needs the n=100,d=10 export; 91_wp6_runtime.R --smoke writes
        # it to results/tr1/wp6/smoke/data/n_100_rep1.csv. If it is not
        # there yet (e.g. this script run standalone), export it here so
        # --smoke is self-contained.
        data_dir.mkdir(parents=True, exist_ok=True)
        smoke_file = data_dir / "n_100_rep1.csv"
        if not smoke_file.exists():
            log(f"[smoke] {smoke_file} not found -- run 91_wp6_runtime.R --smoke first, "
                "or this script will generate an equivalent file itself.")
            rng = np.random.RandomState(60101)
            n, d = 100, 10
            X = rng.normal(loc=3.0, scale=1.0, size=(n, d))
            y = np.ones(n, dtype=int)
            y[-5:] = 0
            pd.DataFrame(np.hstack([X, y.reshape(-1, 1)]),
                        columns=[f"V{i+1}" for i in range(d)] + ["label"]).to_csv(
                smoke_file, index=False)
    else:
        if not args.force:
            ok = idle_check()
            if not ok:
                log("Refusing to start the real WP6 Python grid: other processes are "
                    "present (see list above). Wait for them to finish, or pass --force.")
                sys.exit(1)
        data_dir = WP6 / "data"
        raw_path = WP6 / "91b_wp6_runtime_py_raw.csv"
        done_path = WP6 / "91b_wp6_runtime_py_done.csv"
        reps = args.reps

    log(f"numpy {np.__version__}, pandas {pd.__version__}, torch {torch.__version__}")
    import pyod
    log(f"pyod {pyod.__version__}; torch threads = {torch.get_num_threads()}")

    files = discover_files(data_dir)
    if not files:
        log(f"No cell exports found in {data_dir}. Run 91_wp6_runtime.R "
            f"{'--smoke' if args.smoke else ''} first.")
        sys.exit(1)
    log(f"{len(files)} dataset files: {[p.name for _, _, p in files]}")

    done = read_done_keys(done_path)
    plan = cells_and_methods()
    n_done = n_skip = n_fail = 0

    for grid, cell_value, path in files:
        X = load_cell(path)
        n, d = X.shape
        log(f"\n==== {path.name} (grid={grid}, cell={cell_value}, n={n}, d={d}) ====")

        for method, variant, k, seed_base in plan:
            for rep in range(1, reps + 1):
                key = (method, variant, str(n), str(d), cell_value, str(rep))
                if key in done:
                    log(f"  {method:10s} {variant:6s} rep {rep}/{reps} [checkpoint skip]")
                    n_skip += 1
                    continue
                seed = seed_base if seed_base is not None else rep
                try:
                    elapsed, mem_mb = timed_fit(method, X, k=k, seed=seed)
                    status = "OK"
                    n_done += 1
                except Exception as e:  # noqa: BLE001 -- never abort the sweep
                    elapsed, mem_mb = float("nan"), float("nan")
                    status = f"ERROR: {type(e).__name__}: {e}"
                    err_log = raw_path.parent / "91b_wp6_pyod_errors.log"
                    err_log.parent.mkdir(parents=True, exist_ok=True)
                    with open(err_log, "a", encoding="utf-8") as f:
                        f.write(f"\n=== {grid} {cell_value} {method} {variant} rep{rep} ===\n"
                                f"{traceback.format_exc()}\n")
                    n_fail += 1

                append_raw(raw_path, {
                    "method": method, "variant": variant, "n": n, "d": d,
                    "cell_value": cell_value, "rep": rep, "seed": seed,
                    "time_s": round(elapsed, 4) if elapsed == elapsed else "",
                    "mem_peak_mb": round(mem_mb, 3) if mem_mb == mem_mb else "",
                    "table_load_s": "",     # not applicable -- Python methods load no quantile table
                    "threads": 1,
                    "status": status,
                })
                if status == "OK":
                    append_done(done_path, key)
                    done.add(key)
                shown = f"{elapsed:.4f}s mem={mem_mb:.2f}MB" if elapsed == elapsed else "NA"
                log(f"  {method:10s} {variant:6s} rep {rep}/{reps} {shown} [{status if status == 'OK' else status[:80]}]")

    log(f"\nDone: {n_done} timed, {n_skip} checkpoint-skipped, {n_fail} failed.")
    if args.smoke:
        log("Smoke run complete. Re-run with --smoke to confirm the checkpoint skip fires.")


if __name__ == "__main__":
    main()
