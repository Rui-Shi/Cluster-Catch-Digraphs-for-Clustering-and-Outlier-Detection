"""
08_wp6_pyod_baselines.py -- WP6: three PyOD baselines (ECOD, LUNAR, AutoEncoder) on the
14 real datasets, for the PR-D-26-05767 revision-experiments pipeline (task T6).

What this does
---------------
For each of the 14 datasets in results/datasets_csv/ (excluding manifest.csv):

  - ECOD (deterministic): one fit, saved under seed slot "seed0". A separate,
    unconditional determinism check fits ECOD twice on one dataset (Hepatitis)
    and asserts the two decision_scores_ arrays are identical.
  - LUNAR and AutoEncoder: 5 seeds each (seed in 1..5). Both numpy's and
    torch's global RNGs are seeded per fit (seed_everything), and random_state
    is also passed to the constructor, since PyOD's own randomness (e.g.
    LUNAR's negative sampling, AutoEncoder's weight init / batch shuffling)
    is not fully covered by random_state alone in all pyod versions -- see
    the empirical reproducibility check performed during development.

All three models are constructed with pyod's own defaults except
contamination=0.1 (also pyod's default) and verbose=0 (silenced trainer
chatter; does not affect the fitted model or its scores). Models are fit on
the feature matrix only -- labels are never passed to .fit().

decision_scores_ polarity: verified against pyod's own BaseDetector
docstring ("The higher, the more abnormal") and empirically confirmed here
(ECOD AUC on Musk with y=1 for true outliers computed as ~0.96 during
development), so scores are saved as-is: higher = more outlying.

contamination and label derivation: verified (BaseDetector._process_decision_
scores, and inspection of LUNAR.fit / AutoEncoder.fit source) that the
`contamination` constructor argument does NOT influence model fitting or
decision_scores_ for any of these three models -- it is only consumed
downstream to set threshold_/labels_ via
    threshold_ = percentile(decision_scores_, 100 * (1 - contamination))
    labels_    = (decision_scores_ > threshold_).astype(int)
So a single fit's raw decision_scores_ is reused to derive BOTH the
contamination=0.1 label vector and the true-outlier-rate label vector,
following pyod's own thresholding convention exactly. These companion
label vectors use PyOD's own polarity (1 = outlier, 0 = inlier) -- the
OPPOSITE of this repo's dataset `label` column convention (1 = regular,
0 = outlier) -- to avoid silently relabeling pyod's native output. Columns
are named accordingly (is_outlier_contam01_pyod / is_outlier_truerate_pyod).

Outputs
-------
  results/wp6_scores/<dataset>_<method>_seed<k>.csv
      one column `score` (raw decision_scores_, higher = more outlying),
      row order = dataset row order.
  results/wp6_scores/<dataset>_<method>_seed<k>_labels.csv
      companion file, columns is_outlier_contam01_pyod, is_outlier_truerate_pyod
      (PyOD polarity: 1 = outlier).
  results/wp6_fit_log.csv
      dataset, method, seed, n, d, fit_seconds, status (one row per attempted
      (dataset, method, seed) combo; rewritten incrementally so progress
      survives interruption).
  results/wp6_fit_errors.log
      full tracebacks for any failed fit, appended.
  results/wp6_determinism_check.txt
      result of the ECOD double-fit determinism check.

Checkpointing: a (dataset, method, seed) combo is skipped if both its score
file and label file already exist on disk, so the script resumes cleanly
after interruption. Re-running does not duplicate rows in wp6_fit_log.csv.
"""

import random
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from pyod.models.auto_encoder import AutoEncoder
from pyod.models.ecod import ECOD
from pyod.models.lunar import LUNAR

HERE = Path(__file__).resolve().parent.parent
DATA_DIR = HERE / "results" / "datasets_csv"
SCORES_DIR = HERE / "results" / "tr2" / "wp6_scores"
FIT_LOG_PATH = HERE / "results" / "tr2" / "wp6_fit_log.csv"
ERROR_LOG_PATH = HERE / "results" / "tr2" / "wp6_fit_errors.log"
DETERMINISM_LOG_PATH = HERE / "results" / "tr2" / "wp6_determinism_check.txt"

LUNAR_AE_SEEDS = [1, 2, 3, 4, 5]
DETERMINISM_DATASET = "Hepatitis"  # smallest dataset -> fast double-fit check


def log(msg):
    print(msg, flush=True)


def seed_everything(seed):
    s = int(seed) if seed is not None else 0
    np.random.seed(s)
    random.seed(s)
    torch.manual_seed(s)


def load_manifest():
    mf = pd.read_csv(DATA_DIR / "manifest.csv")
    return mf.set_index("dataset")


def dataset_csv_names():
    names = []
    for p in sorted(DATA_DIR.glob("*.csv")):
        if p.stem == "manifest":
            continue
        names.append(p.stem)
    return names


def load_dataset(name):
    df = pd.read_csv(DATA_DIR / f"{name}.csv")
    y = df["label"].values.astype(int)  # repo convention: 1 = regular, 0 = outlier
    X = df.drop(columns=["label"]).values.astype(np.float64)
    return X, y


def score_path(dataset, method, seed):
    return SCORES_DIR / f"{dataset}_{method}_seed{seed}.csv"


def labels_path(dataset, method, seed):
    return SCORES_DIR / f"{dataset}_{method}_seed{seed}_labels.csv"


def save_scores(path, scores):
    pd.DataFrame({"score": scores}).to_csv(path, index=False)


def save_labels(path, scores, contamination_default, contamination_true):
    # Replicate pyod's own BaseDetector._process_decision_scores convention.
    thr01 = np.percentile(scores, 100 * (1 - contamination_default))
    thr_true = np.percentile(scores, 100 * (1 - contamination_true))
    lab01 = (scores > thr01).astype(int)
    lab_true = (scores > thr_true).astype(int)
    pd.DataFrame(
        {
            "is_outlier_contam01_pyod": lab01,
            "is_outlier_truerate_pyod": lab_true,
        }
    ).to_csv(path, index=False)


def build_model(method, seed):
    if method == "ECOD":
        return ECOD(contamination=0.1)
    if method == "LUNAR":
        return LUNAR(contamination=0.1, random_state=seed, verbose=0)
    if method == "AutoEncoder":
        return AutoEncoder(contamination=0.1, random_state=seed, verbose=0)
    raise ValueError(f"unknown method {method}")


class FitLog:
    """Incrementally-persisted fit log, keyed by (dataset, method, seed) so
    resumed runs neither duplicate rows nor lose prior entries."""

    def __init__(self, path):
        self.path = path
        self.rows = {}
        if path.exists():
            existing = pd.read_csv(path)
            for _, r in existing.iterrows():
                key = (r["dataset"], r["method"], int(r["seed"]))
                self.rows[key] = r.to_dict()

    def has(self, key):
        return key in self.rows

    def add(self, dataset, method, seed, n, d, fit_seconds, status):
        key = (dataset, method, seed)
        self.rows[key] = dict(
            dataset=dataset,
            method=method,
            seed=seed,
            n=n,
            d=d,
            fit_seconds=fit_seconds,
            status=status,
        )
        self.flush()

    def flush(self):
        df = pd.DataFrame(list(self.rows.values()))
        if len(df):
            df = df.sort_values(["dataset", "method", "seed"]).reset_index(drop=True)
        df.to_csv(self.path, index=False)


def run_determinism_check():
    log(f"=== ECOD determinism check on {DETERMINISM_DATASET} ===")
    X, _ = load_dataset(DETERMINISM_DATASET)
    m1 = ECOD(contamination=0.1)
    m1.fit(X)
    m2 = ECOD(contamination=0.1)
    m2.fit(X)
    s1 = np.asarray(m1.decision_scores_, dtype=float)
    s2 = np.asarray(m2.decision_scores_, dtype=float)
    identical = bool(np.array_equal(s1, s2))
    max_abs_diff = float(np.max(np.abs(s1 - s2))) if s1.shape == s2.shape else float("nan")
    msg = (
        f"dataset={DETERMINISM_DATASET} n={X.shape[0]} d={X.shape[1]}\n"
        f"identical (np.array_equal): {identical}\n"
        f"max_abs_diff: {max_abs_diff}\n"
    )
    log(msg)
    DETERMINISM_LOG_PATH.write_text(msg, encoding="utf-8")
    assert identical, "ECOD is not deterministic across repeated fits -- investigate before proceeding"


def run_one(fit_log, dataset, method, seed, X, y, contamination_true):
    n, d = X.shape
    sp = score_path(dataset, method, seed)
    lp = labels_path(dataset, method, seed)
    key = (dataset, method, seed)

    if sp.exists() and lp.exists():
        log(f"SKIP (exists)  {dataset:15s} {method:12s} seed{seed}")
        if not fit_log.has(key):
            fit_log.add(dataset, method, seed, n, d, float("nan"), "OK (skipped, pre-existing output)")
        return

    seed_everything(seed)
    t0 = time.time()
    try:
        model = build_model(method, seed)
        model.fit(X)
        elapsed = time.time() - t0
        scores = np.asarray(model.decision_scores_, dtype=float).ravel()
        if scores.shape[0] != n:
            raise RuntimeError(f"decision_scores_ length {scores.shape[0]} != n {n}")
        save_scores(sp, scores)
        save_labels(lp, scores, 0.1, contamination_true)
        fit_log.add(dataset, method, seed, n, d, elapsed, "OK")
        log(f"OK             {dataset:15s} {method:12s} seed{seed} n={n} d={d} t={elapsed:.2f}s")
    except Exception as e:  # noqa: BLE001 -- deliberately broad, must not abort the sweep
        elapsed = time.time() - t0
        tb = traceback.format_exc()
        with open(ERROR_LOG_PATH, "a", encoding="utf-8") as f:
            f.write(f"\n=== {dataset} {method} seed{seed} ===\n{tb}\n")
        fit_log.add(dataset, method, seed, n, d, elapsed, f"FAILED: {type(e).__name__}: {e}")
        log(f"FAIL           {dataset:15s} {method:12s} seed{seed}: {type(e).__name__}: {e}")


def main():
    SCORES_DIR.mkdir(parents=True, exist_ok=True)

    log(f"numpy {np.__version__}, pandas {pd.__version__}, torch {torch.__version__}")
    import pyod

    log(f"pyod {pyod.__version__}")
    log(f"ECOD defaults:        {ECOD().get_params()}")
    log(f"LUNAR defaults:       {LUNAR().get_params()}")
    log(f"AutoEncoder defaults: {AutoEncoder().get_params()}")

    run_determinism_check()

    manifest = load_manifest()
    names = dataset_csv_names()
    log(f"\n{len(names)} datasets: {names}\n")

    fit_log = FitLog(FIT_LOG_PATH)

    for dataset in names:
        X, y = load_dataset(dataset)
        n = len(y)
        n0_data = int((y == 0).sum())
        n0_manifest = int(manifest.loc[dataset, "n_outliers"]) if dataset in manifest.index else None
        if n0_manifest is not None and n0_manifest != n0_data:
            log(
                f"WARNING: {dataset} manifest n_outliers={n0_manifest} != "
                f"data label==0 count={n0_data}; using data-derived value"
            )
        contamination_true = n0_data / n

        # ECOD: single deterministic fit, seed slot 0
        run_one(fit_log, dataset, "ECOD", 0, X, y, contamination_true)

        # LUNAR, AutoEncoder: 5 seeds each
        for method in ["LUNAR", "AutoEncoder"]:
            for seed in LUNAR_AE_SEEDS:
                run_one(fit_log, dataset, method, seed, X, y, contamination_true)

    n_ok = sum(1 for r in fit_log.rows.values() if str(r["status"]).startswith("OK"))
    n_fail = sum(1 for r in fit_log.rows.values() if str(r["status"]).startswith("FAILED"))
    log(f"\nDone. {n_ok} OK, {n_fail} FAILED, {len(fit_log.rows)} total logged combos.")
    if n_fail:
        log(f"See {ERROR_LOG_PATH} for tracebacks.")


if __name__ == "__main__":
    main()
