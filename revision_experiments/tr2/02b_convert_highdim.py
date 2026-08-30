"""
02b_convert_highdim.py

Task T1 (high-dimensional dataset acquisition) for the PR-D-26-05767 revision
experiments pipeline.

Converts the four raw ODDS/ADBench-format outlier-detection benchmark files
(Musk, Speech, InternetAds as .npz with arrays X/y; Arrhythmia as .mat with
arrays X/y) into plain CSV files that R can read with base read.csv(), with
no dependency on R.matlab's npz support (there is none) or Python<->R bridge
packages.

Output convention (deliberately mirrors the *source* ODDS/ADBench polarity,
NOT the RealData_Collection.R polarity): each row is the raw feature vector
followed by a final column named "label" where 1 = outlier, 0 = inlier. This
matches the y=1-means-outlier convention used by both ODDS .mat files and
ADBench .npz files. The R loader (02_load_data.R) is responsible for flipping
this to the repo's own convention (1 = regular, 0 = outlier) exactly as
RealData_Collection.R already does for glass/vertebral/vowels/thyroid.mat via
`ifelse(y == "1", 0, 1)`.

Run from anywhere; paths are resolved relative to this script's location.
"""

import os
import sys

import numpy as np

try:
    import scipy.io as sio
except ImportError:
    sio = None

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DATA_DIR = os.path.join(REPO_ROOT, "data", "outlier_detection")

# (name, source filename, loader) -- loader returns (X, y) with y in {0,1}, 1=outlier
DATASETS = [
    ("musk", "musk.npz", "npz"),
    ("speech", "speech.npz", "npz"),
    ("internetads", "internetads.npz", "npz"),
    ("arrhythmia", "arrhythmia.mat", "mat"),
]

EXPECTED = {
    "musk": dict(n=3062, d=166, n_outliers=97),
    "arrhythmia": dict(n=452, d=274, n_outliers=66),
    "speech": dict(n=3686, d=400, n_outliers=61),
    "internetads": dict(n=1966, d=1555, n_outliers=368),
}


def load_npz(path):
    d = np.load(path)
    X = np.asarray(d["X"], dtype=float)
    y = np.asarray(d["y"]).ravel().astype(int)
    return X, y


def load_mat(path):
    if sio is None:
        raise RuntimeError(
            "scipy is not available in this venv, and file '%s' is a .mat file. "
            "Install scipy or provide a .npz alternative." % path
        )
    d = sio.loadmat(path)
    X = np.asarray(d["X"], dtype=float)
    y = np.asarray(d["y"]).ravel().astype(int)
    return X, y


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    print("=== 02b_convert_highdim.py ===")
    print("DATA_DIR:", DATA_DIR)

    all_ok = True
    for name, fname, kind in DATASETS:
        src_path = os.path.join(DATA_DIR, fname)
        if not os.path.exists(src_path):
            print("MISSING SOURCE FILE: %s (dataset %s)" % (src_path, name))
            all_ok = False
            continue

        if kind == "npz":
            X, y = load_npz(src_path)
        elif kind == "mat":
            X, y = load_mat(src_path)
        else:
            raise ValueError(kind)

        n, d = X.shape
        n_out = int((y == 1).sum())

        nan_in_X = bool(np.isnan(X).any())
        uniq_y = sorted(np.unique(y).tolist())

        print("\n--- %s ---" % name)
        print("  source: %s" % fname)
        print("  n=%d  d=%d  n_outliers(y==1)=%d" % (n, d, n_out))
        print("  y unique values: %s" % uniq_y)
        print("  any NaN in X: %s" % nan_in_X)

        exp = EXPECTED.get(name)
        if exp is not None:
            match = (n == exp["n"]) and (d == exp["d"]) and (n_out == exp["n_outliers"])
            print(
                "  EXPECTED n=%d d=%d n_outliers=%d -> MATCH=%s"
                % (exp["n"], exp["d"], exp["n_outliers"], match)
            )
            if not match:
                all_ok = False

        if nan_in_X:
            print("  WARNING: NaN present in feature matrix for %s" % name)
            all_ok = False

        # Write CSV: feature columns V1..Vd + final column 'label' (1=outlier, matches source polarity)
        colnames = ["V%d" % (j + 1) for j in range(d)]
        out_path = os.path.join(DATA_DIR, "%s_raw.csv" % name)

        # Use numpy for fast CSV writing; header manually prepended.
        header = ",".join(colnames + ["label"])
        out_arr = np.hstack([X, y.reshape(-1, 1).astype(float)])
        np.savetxt(out_path, out_arr, delimiter=",", header=header, comments="", fmt="%.10g")

        print("  wrote: %s (%d rows x %d cols incl. label)" % (out_path, out_arr.shape[0], out_arr.shape[1]))

    print("\n=== SUMMARY ===")
    if all_ok:
        print("All datasets converted and match expected gate values. OK.")
    else:
        print("One or more datasets FAILED expected-value checks or had missing files. See above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
