"""
84a_wp5_fetch_convert.py -- WP5 (AE.4, R1.8, R3.8, R5.6): data acquisition and
conversion for the four real data sets above d = 21.

Everything this script does is declared in advance in
revision_experiments/tr1/WP5_PROTOCOL.md. Read that first; this file is the
implementation, not the specification.

What it does, per WP5_PROTOCOL.md:

  1. Downloads 20_letter.npz and 24_mnist.npz from the ADBench GitHub mirror
     (the same mirror already used for Musk/Speech/InternetAds) into
     results/tr1/wp5/data/raw/.
  2. Verifies letter's (n, d, n_outliers) against the inventory (1600, 32,
     100) -- fatal on mismatch. For mnist, verifies d = 100 (fatal) and
     prints the observed n / n_outliers against the literature figures
     (n=7603, n_outliers=700) as an informational check only (protocol S6:
     no file in this repo pre-declares an authoritative mnist outlier count).
  3. Flips label polarity source (1=outlier) -> repo convention (1=regular,
     0=outlier).
  4. Applies scale_R_safe() column-wise to letter (full n) and to mnist
     (after subsampling) -- the numpy reimplementation of
     revision_experiments/tr2/02_load_data.R:142-153, matching R's mad()
     default constant (1.4826) exactly.
  5. Draws mnist's n=1000 contamination-preserving subsample, seed 20260905
     (protocol S2), before scaling (scaling is computed on the subsample's
     own columns, since that subsample is the data set WP5 actually studies
     -- not on the full n=7603 set it was drawn from).
  6. Copies musk (Musk_sub1000.csv) and arrhythmia (Arrhythmia.csv) from
     results/datasets_csv/ into results/tr1/wp5/data/, unchanged, after
     verifying their (n, d, n_outliers) against the inventory (1000, 166, 32
     and 452, 274, 66) -- fatal on mismatch. No re-scaling.
  7. Writes results/tr1/wp5/data/manifest.csv.

This is data preparation, not the experiment -- it is run directly (not left
for a human to launch), per the WP5 task brief.

Usage (from anywhere; paths resolved from this file's location):
    .venv/python.exe revision_experiments/tr1/84a_wp5_fetch_convert.py [--skip-download]

    --skip-download   reuse the raw .npz files already in data/raw/ instead
                       of re-fetching (idempotent re-runs during development;
                       still re-verifies and re-converts).
"""

import argparse
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
WP5_DIR = HERE.parent / "results" / "tr1" / "wp5"
RAW_DIR = WP5_DIR / "data" / "raw"
DATA_DIR = WP5_DIR / "data"
EXISTING_DATASETS_CSV = HERE.parent / "results" / "datasets_csv"

ADBENCH_BASE = "https://raw.githubusercontent.com/Minqi824/ADBench/main/adbench/datasets/Classical"
LETTER_URL = f"{ADBENCH_BASE}/20_letter.npz"
MNIST_URL = f"{ADBENCH_BASE}/24_mnist.npz"

# WP5_PROTOCOL.md S1 -- letter's expected (n, d, n_outliers) is a hard gate.
LETTER_EXPECTED = dict(n=1600, d=32, n_outliers=100)
# mnist: only d is gated (S6); n / n_outliers are printed against the
# literature figure for the record, not enforced.
MNIST_EXPECTED_D = 100
MNIST_LITERATURE_N = 7603
MNIST_LITERATURE_N_OUTLIERS = 700

# WP5_PROTOCOL.md S1 -- musk/arrhythmia are reused as-is; still gated.
MUSK_EXPECTED = dict(n=1000, d=166, n_outliers=32)
ARRHYTHMIA_EXPECTED = dict(n=452, d=274, n_outliers=66)

MNIST_SUBSAMPLE_SEED = 20260905
MNIST_SUBSAMPLE_N = 1000


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------

def download(url, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"  fetching {url}")
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    if len(data) == 0:
        raise RuntimeError(f"downloaded 0 bytes from {url}")
    dest.write_bytes(data)
    log(f"  wrote {dest} ({len(data)} bytes)")


# ---------------------------------------------------------------------------
# scale_R_safe -- numpy reimplementation of
# revision_experiments/tr2/02_load_data.R:142-153. R's mad(x) default is
# 1.4826 * median(|x - median(x)|); replicated exactly here so this script's
# scaling is numerically identical to the R function it cites, not merely
# similar in spirit.
# ---------------------------------------------------------------------------

MAD_CONSTANT = 1.4826

def r_mad(x):
    med = np.median(x)
    return MAD_CONSTANT * np.median(np.abs(x - med))


def scale_R_safe_col(x):
    med = np.median(x)
    madn = r_mad(x)
    if madn > 0:
        return (x - med) / madn, "madn"
    sdx = np.std(x, ddof=1) if len(x) > 1 else 0.0
    if sdx > 0:
        return (x - med) / sdx, "sd_fallback"
    return x - med, "constant_fallback"


def scale_R_safe_matrix(X, name):
    n, d = X.shape
    out = np.empty_like(X, dtype=float)
    branch_counts = {"madn": 0, "sd_fallback": 0, "constant_fallback": 0}
    for j in range(d):
        out[:, j], branch = scale_R_safe_col(X[:, j])
        branch_counts[branch] += 1
    log(f"  [{name}] scale_R_safe: {branch_counts['madn']}/{d} columns via MADN, "
        f"{branch_counts['sd_fallback']}/{d} via SD fallback, "
        f"{branch_counts['constant_fallback']}/{d} constant fallback")
    if np.isnan(out).any():
        raise RuntimeError(f"[{name}] NaN produced by scale_R_safe_matrix -- should be impossible")
    return out


# ---------------------------------------------------------------------------
# label polarity: source (1=outlier, 0=regular) -> repo (1=regular, 0=outlier)
# ---------------------------------------------------------------------------

def flip_to_repo_polarity(y_source):
    y_source = np.asarray(y_source).astype(int)
    if not set(np.unique(y_source).tolist()) <= {0, 1}:
        raise RuntimeError(f"unexpected label values: {np.unique(y_source)}")
    return 1 - y_source  # 1=outlier -> 0; 0=regular -> 1


# ---------------------------------------------------------------------------
# contamination-preserving subsample, WP5_PROTOCOL.md S2
# ---------------------------------------------------------------------------

def contamination_preserving_subsample(X, y_repo, n_target, seed):
    """y_repo: repo polarity (1=regular, 0=outlier)."""
    rng = np.random.default_rng(seed)
    idx_reg = np.where(y_repo == 1)[0]
    idx_out = np.where(y_repo == 0)[0]
    n_out_full, n_reg_full = len(idx_out), len(idx_reg)
    n_full = n_out_full + n_reg_full
    cont = n_out_full / n_full
    n_out_sub = int(round(cont * n_target))
    n_reg_sub = n_target - n_out_sub
    if n_out_sub > n_out_full or n_reg_sub > n_reg_full:
        raise RuntimeError(
            f"subsample request exceeds pool size: need {n_out_sub} of {n_out_full} "
            f"outliers, {n_reg_sub} of {n_reg_full} regulars"
        )
    pick_out = rng.choice(idx_out, size=n_out_sub, replace=False)
    pick_reg = rng.choice(idx_reg, size=n_reg_sub, replace=False)
    pick = np.concatenate([pick_out, pick_reg])
    pick.sort()  # stable, reproducible row order; not semantically required
    log(f"  contamination-preserving subsample: full n={n_full} (cont={cont:.4%}) "
        f"-> n={n_target} ({n_out_sub} outliers, {n_reg_sub} regulars), seed={seed}")
    return X[pick, :], y_repo[pick]


# ---------------------------------------------------------------------------
# gate: verify (n, d, n_outliers) against the protocol's expected values
# ---------------------------------------------------------------------------

def gate(name, n, d, n_outliers, expected, fatal=True):
    ok = (n == expected["n"] and d == expected["d"] and n_outliers == expected["n_outliers"])
    tag = "MATCH" if ok else ("MISMATCH -- FATAL" if fatal else "MISMATCH -- informational only")
    log(f"  [{name}] n={n} d={d} n_outliers={n_outliers} | expected "
        f"n={expected['n']} d={expected['d']} n_outliers={expected['n_outliers']} | {tag}")
    if not ok and fatal:
        raise RuntimeError(f"[{name}] gate failed: does not match WP5_PROTOCOL.md's expected values")
    return ok


def write_csv(X_repo, y_repo, path: Path):
    d = X_repo.shape[1]
    df = pd.DataFrame(X_repo, columns=[f"V{j+1}" for j in range(d)])
    df["label"] = y_repo.astype(int)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    log(f"  wrote {path} ({df.shape[0]} rows x {df.shape[1]} cols incl. label)")


# ---------------------------------------------------------------------------
# per-dataset handlers
# ---------------------------------------------------------------------------

def handle_letter(skip_download):
    log("\n=== letter ===")
    raw_path = RAW_DIR / "letter.npz"
    if not skip_download or not raw_path.exists():
        download(LETTER_URL, raw_path)
    else:
        log(f"  --skip-download: reusing {raw_path}")

    d = np.load(raw_path)
    X, y_source = d["X"], d["y"]
    n, dd = X.shape
    n_outliers = int((y_source == 1).sum())
    gate("letter", n, dd, n_outliers, LETTER_EXPECTED, fatal=True)

    if np.isnan(X).any():
        raise RuntimeError("letter: NaN in raw feature matrix")

    y_repo = flip_to_repo_polarity(y_source)
    X_scaled = scale_R_safe_matrix(X, "letter")
    write_csv(X_scaled, y_repo, DATA_DIR / "letter.csv")

    return dict(
        dataset="letter", n=n, d=dd, n_outliers=n_outliers,
        contamination=n_outliers / n,
        preprocessing="robust median/MADN per column (scale_R_safe, SD/constant fallback)",
        source=LETTER_URL, subsample_seed="",
    )


def handle_mnist(skip_download):
    log("\n=== mnist ===")
    raw_path = RAW_DIR / "mnist.npz"
    if not skip_download or not raw_path.exists():
        download(MNIST_URL, raw_path)
    else:
        log(f"  --skip-download: reusing {raw_path}")

    d = np.load(raw_path)
    X_full, y_source_full = d["X"], d["y"]
    n_full, dd = X_full.shape
    n_outliers_full = int((y_source_full == 1).sum())
    log(f"  [mnist, full source set] n={n_full} d={dd} n_outliers={n_outliers_full} "
        f"({n_outliers_full / n_full:.4%})")
    log(f"  [mnist] literature reference: n={MNIST_LITERATURE_N} "
        f"n_outliers={MNIST_LITERATURE_N_OUTLIERS} "
        f"({MNIST_LITERATURE_N_OUTLIERS / MNIST_LITERATURE_N:.4%}) -- informational only, not gated")
    if dd != MNIST_EXPECTED_D:
        raise RuntimeError(f"mnist: d={dd}, expected {MNIST_EXPECTED_D} -- this is the one hard gate")
    if np.isnan(X_full).any():
        raise RuntimeError("mnist: NaN in raw feature matrix")

    y_repo_full = flip_to_repo_polarity(y_source_full)
    X_sub, y_sub = contamination_preserving_subsample(
        X_full, y_repo_full, MNIST_SUBSAMPLE_N, MNIST_SUBSAMPLE_SEED
    )
    n_outliers_sub = int((y_sub == 0).sum())

    X_scaled = scale_R_safe_matrix(X_sub, "mnist")
    write_csv(X_scaled, y_sub, DATA_DIR / "mnist.csv")

    return dict(
        dataset="mnist", n=MNIST_SUBSAMPLE_N, d=dd, n_outliers=n_outliers_sub,
        contamination=n_outliers_sub / MNIST_SUBSAMPLE_N,
        preprocessing="robust median/MADN per column (scale_R_safe, SD/constant fallback), "
                      "computed on the n=1000 subsample",
        source=MNIST_URL, subsample_seed=MNIST_SUBSAMPLE_SEED,
    )


def handle_reused(name, src_csv, out_name, expected):
    log(f"\n=== {name} (reused, no re-scaling) ===")
    src_path = EXISTING_DATASETS_CSV / src_csv
    if not src_path.exists():
        raise RuntimeError(f"[{name}] expected existing CSV not found: {src_path}")
    df = pd.read_csv(src_path)
    if df.columns[-1] != "label":
        raise RuntimeError(f"[{name}] {src_path} does not end in a 'label' column")
    n, d = df.shape[0], df.shape[1] - 1
    n_outliers = int((df["label"] == 0).sum())
    gate(name, n, d, n_outliers, expected, fatal=True)

    out_path = DATA_DIR / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    log(f"  copied {src_path} -> {out_path} unchanged ({n} rows x {d+1} cols incl. label)")

    return dict(
        dataset=name, n=n, d=d, n_outliers=n_outliers, contamination=n_outliers / n,
        preprocessing=f"unchanged from results/datasets_csv/{src_csv} "
                      "(robust median/MADN per column, existing SD/constant fallback)",
        source=f"results/datasets_csv/{src_csv} (see that folder's own manifest.csv for provenance)",
        subsample_seed=20260716 if name == "musk" else "",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    log("=== WP5 data acquisition and conversion ===")
    log(f"Output dir: {DATA_DIR}")

    rows = []
    rows.append(handle_letter(args.skip_download))
    rows.append(handle_mnist(args.skip_download))
    rows.append(handle_reused("musk", "Musk_sub1000.csv", "musk.csv", MUSK_EXPECTED))
    rows.append(handle_reused("arrhythmia", "Arrhythmia.csv", "arrhythmia.csv", ARRHYTHMIA_EXPECTED))

    manifest = pd.DataFrame(rows, columns=[
        "dataset", "n", "d", "n_outliers", "contamination", "preprocessing",
        "source", "subsample_seed",
    ])
    manifest_path = DATA_DIR / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    log(f"\nwrote {manifest_path}")
    log("\n=== all four data sets converted ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
