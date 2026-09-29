#!/usr/bin/env python
"""revision_experiments/tr1/90b_wp7_hdbscan.py

WP7 helper (R3.6, clustering quality) -- fits HDBSCAN on one regenerated
WP8 cell and writes back native cluster labels. Called per cell, via a CSV
round trip, from 90_wp7_clustering_quality.R's system2() call. See
WP7_PROTOCOL.md section 6.

Defaults are IDENTICAL to WP4's own HDBSCAN use (tr1/81_wp4_baselines.py,
WP4_PROTOCOL.md section 2): min_cluster_size=5, min_samples=None. This is
the same hdbscan 0.8.44 install WP4 already verified in this .venv.

Usage:
    python.exe 90b_wp7_hdbscan.py <input_csv> <output_csv>

<input_csv>  -- feature matrix only, one row per point, header row present
                (column names are irrelevant, only values are read).
<output_csv> -- one column, "label" -- hdbscan's native cluster labels,
                -1 = noise. Not remapped in any way: ARI/NMI/AMI/k-hat are
                all invariant to a relabelling, so no unassigned -> NA
                conversion happens here; that is done in R after reading
                this file back (label == -1 -> NA).

Exits non-zero with a one-line message on stderr on any failure (bad input,
fit exception), so the R caller's system2(..., ) exit-status check can
distinguish "wrote a labels file" from "did not."
"""
import sys


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: 90b_wp7_hdbscan.py <input_csv> <output_csv>\n")
        sys.exit(2)
    in_csv, out_csv = sys.argv[1], sys.argv[2]

    try:
        import pandas as pd
        import hdbscan
    except Exception as e:  # pragma: no cover - environment problem, not data
        sys.stderr.write(f"import failure: {e}\n")
        sys.exit(3)

    try:
        X = pd.read_csv(in_csv).values
    except Exception as e:
        sys.stderr.write(f"failed to read {in_csv}: {e}\n")
        sys.exit(4)

    if X.shape[0] < 2:
        sys.stderr.write(f"input has {X.shape[0]} rows, need >= 2\n")
        sys.exit(5)

    try:
        model = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=None)
        labels = model.fit_predict(X)
    except Exception as e:
        sys.stderr.write(f"HDBSCAN fit failed: {e}\n")
        sys.exit(6)

    try:
        pd.DataFrame({"label": labels}).to_csv(out_csv, index=False)
    except Exception as e:
        sys.stderr.write(f"failed to write {out_csv}: {e}\n")
        sys.exit(7)


if __name__ == "__main__":
    main()
