#!/usr/bin/env python3
"""Line plots of the Supplement (Lineplots_Clustering_<method>_<dist>.png), drawn from the printed tables.

The six figures plot the mean ARI of Supplement Tables S2-S4 (uniform) and S6-S8 (Gaussian):
  panel 1: 3 clusters, x = n, one line per d;  panel 2: 3 clusters, x = d, one line per n;
  panel 3: n = 200, x = number of clusters, one line per d.
The values are read from Supplementary_Material.tex (each table found by its \\label), so the figures always
show the printed table values. The plotting code is that of the authors' notebook
figure_unused/Line_Plots_Clustering.ipynb (cells 0-5; they differ only in title, figure size and file name),
whose ARI values were typed in by hand and had gone out of date (2026-11-03); only the data source changed.
After drawing, every plotted point is read back from the figure and compared with its table value; any
difference stops the script.

USAGE
  python 39_line_plots.py <Supplementary_Material.tex> <out_dir>
Run 2026-11-03 with Python 3.13.13, seaborn 0.13.2, matplotlib 3.11.2, pandas 3.1.0, numpy 2.5.4.
PNGs are saved at 100 dpi with a tight bounding box, which gives the about 2000-pixel-wide files the
Supplement includes.

OUTPUT (<out_dir>)
  Lineplots_Clustering_{KS,RK,UN}-CCD_{Uniform,Gaussian}.png
  lineplot_data.csv     method, dist, n, clusters, d, ari_printed (360 rows: 6 tables x 60 ARI values)
  plotted_points.json   every line of every panel as drawn (read back from matplotlib)
"""
import csv
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

D_VALUES = [2, 3, 5, 10, 20]
N_VALUES = [50, 100, 200, 500]
K_VALUES = [2, 3, 5]

# (table label, method, distribution, notebook figsize, notebook suptitle)
FIGURES = [
    ("tab:uni_KS-CCDs_Apdx", "KS-CCD", "Uniform", (20, 5), "The Line Plots of ARIs for KS-CCDs Under the Uniform Settings"),
    ("tab:uni_RK-CCDs_Apdx", "RK-CCD", "Uniform", (20, 6), "The Line Plots of ARIs for RK-CCDs Under the Uniform Settings"),
    ("tab:uni_UN-CCDs_Apdx", "UN-CCD", "Uniform", (20, 6), "The Line Plots of ARIs for UN-CCDs Under the Uniform Settings"),
    ("tab:Gau_KS-CCDs_Apdx", "KS-CCD", "Gaussian", (20, 6), "The Line Plots of ARIs for KS-CCDs Under the Gaussian Settings"),
    ("tab:Gau_RK-CCDs_Apdx", "RK-CCD", "Gaussian", (20, 6), "The Line Plots of ARIs for RK-CCDs Under the Gaussian Settings"),
    ("tab:Gau_UN-CCDs_Apdx", "UN-CCD", "Gaussian", (20, 6), "The Line Plots of ARIs for UN-CCDs Under the Gaussian Settings"),
]


def fail(msg):
    sys.exit("FAIL: " + msg)


def strip_markup(text):
    """Remove \\textcolor{red}{...} and \\textbf{...} wrappers, keeping their content."""
    prev = None
    while prev != text:
        prev = text
        text = re.sub(r"\\textcolor\{[^{}]*\}\{([^{}]*)\}", r"\1", text)
        text = re.sub(r"\\textbf\{([^{}]*)\}", r"\1", text)
    return text.strip()


def parse_table(tex, label, method, dist):
    """The 60 ARI values (4 n x 3 clusters x 5 d) of one table, as printed."""
    key = "\\label{%s}" % label
    if tex.count(key) != 1:
        fail("%s occurs %d times in the .tex" % (key, tex.count(key)))
    lpos = tex.index(key)
    tstart = tex.rfind("\\begin{tabular}", 0, lpos)
    tend = tex.find("\\end{tabular}", tstart)
    tabenv = tex.rfind("\\begin{table}", 0, lpos)
    if tstart < 0 or tend < 0 or tend > lpos or tabenv > tstart:
        fail("no tabular in the table environment before %s" % key)
    body = tex[tstart:tend]
    body = re.sub(r"^\\begin\{tabular\}\{[^}]*\}", "", body)
    body = re.sub(r"\\cline\{[^}]*\}", "", body)
    body = body.replace("\\hline", "")
    dvals = [int(x) for x in re.findall(r"\$d=(\d+)\$", body)]
    if dvals != D_VALUES:
        fail("%s: header d values %r, expected %r" % (label, dvals, D_VALUES))
    rows, cur_n = [], None
    for chunk in re.split(r"\\\\", body):
        chunk = chunk.strip()
        if not re.search(r"\d+\s+clusters", chunk):
            continue
        cells = [strip_markup(c) for c in chunk.split("&")]
        if len(cells) != 17:
            fail("%s: row has %d cells, expected 17: %r" % (label, len(cells), chunk[:120]))
        nm = re.search(r"n\s*=\s*(\d+)", cells[0])
        if nm:
            cur_n = int(nm.group(1))
        km = re.fullmatch(r"(\d+) clusters", cells[1])
        if not km:
            fail("%s: bad cluster label %r" % (label, cells[1]))
        for di, d in enumerate(D_VALUES):
            ari = cells[2 + 3 * di]  # columns per d: ARI, Sil, SR
            if not re.fullmatch(r"\d\.\d{3}", ari):
                fail("%s: ARI cell %r is not a 3-decimal number" % (label, ari))
            rows.append({"method": method, "dist": dist, "n": cur_n, "clusters": int(km.group(1)), "d": d,
                         "ari_printed": ari})
    expected = {(n, k, d) for n in N_VALUES for k in K_VALUES for d in D_VALUES}
    if len(rows) != 60 or {(r["n"], r["clusters"], r["d"]) for r in rows} != expected:
        fail("%s: %d ARI values, expected exactly 60 (4 n x 3 clusters x 5 d)" % (label, len(rows)))
    return rows


def notebook_data(rows):
    """The notebook's data dict, built from the table, in the notebook's row order (d, then n, then clusters)."""
    v = {(r["d"], r["n"], r["clusters"]): float(r["ari_printed"]) for r in rows}
    data = {"n": [], "d": [], "clusters": [], "ARI": []}
    for d in D_VALUES:
        for n in N_VALUES:
            for k in K_VALUES:
                data["n"].append(n)
                data["d"].append(d)
                data["clusters"].append(k)
                data["ARI"].append(v[(d, n, k)])
    return data


def draw(data, figsize, suptitle, out_png):
    """The plotting code of notebook cells 0-5, unchanged except figsize/suptitle/file name and savefig settings."""
    df = pd.DataFrame(data)
    df['clusters'] = df['clusters'].astype(int)
    df['d'] = df['d'].astype(int)
    df_c3 = df[df['clusters'] == 3].copy()
    df_n200 = df[df['n'] == 200].copy()
    font_sizes = {
        'suptitle': 24,
        'title': 24,
        'label': 22,
        'ticks': 20,
        'legend_title': 18,
        'legend_items': 16
    }
    plt.style.use('seaborn-v0_8-white')
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    fig.suptitle(suptitle, fontsize=font_sizes['suptitle'])
    sns.lineplot(data=df_c3, x='n', y='ARI', style='d', ax=axes[0], color='black', markers=True, dashes=True)
    axes[0].set_title('ARI vs. Sample Size $n$ (fix 3 clusters)', fontsize=font_sizes['title'])
    axes[0].set_xlabel('Data Size ($n$)', fontsize=font_sizes['label'])
    axes[0].set_ylabel('ARI', fontsize=font_sizes['label'])
    axes[0].tick_params(axis='both', which='major', labelsize=font_sizes['ticks'])
    axes[0].set_xticks([50, 100, 200, 500])
    axes[0].legend(title=r'Dims ($d$)', loc='best', title_fontsize=font_sizes['legend_title'], fontsize=font_sizes['legend_items'])
    sns.lineplot(data=df_c3, x='d', y='ARI', style='n', ax=axes[1], color='black', markers=True, dashes=True)
    axes[1].set_title('ARI vs. Dimensions $d$ (fix 3 clusters)', fontsize=font_sizes['title'])
    axes[1].set_xlabel('Dimensions ($d$)', fontsize=font_sizes['label'])
    axes[1].set_ylabel('ARI', fontsize=font_sizes['label'])
    axes[1].tick_params(axis='both', which='major', labelsize=font_sizes['ticks'])
    axes[1].set_xticks([2, 3, 5, 10, 20])
    axes[1].legend(title=r'Samples ($n$)', loc='best', title_fontsize=font_sizes['legend_title'], fontsize=font_sizes['legend_items'])
    sns.lineplot(data=df_n200, x='clusters', y='ARI', style='d', ax=axes[2], color='black', markers=True, dashes=True)
    axes[2].set_title('ARI vs. Number of Clusters (fix $n=200$)', fontsize=font_sizes['title'])
    axes[2].set_xlabel('Number of Clusters', fontsize=font_sizes['label'])
    axes[2].set_ylabel('ARI', fontsize=font_sizes['label'])
    axes[2].tick_params(axis='both', which='major', labelsize=font_sizes['ticks'])
    axes[2].set_xticks([2, 3, 5])
    axes[2].legend(title=r'Dims ($d$)', loc='best', title_fontsize=font_sizes['legend_title'], fontsize=font_sizes['legend_items'])
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(out_png, dpi=100, bbox_inches='tight')
    return fig


def read_back(fig):
    """Every data line of every panel, labelled by the legend entry with the same marker and line style."""
    panels = []
    for pi, ax in enumerate(fig.axes):
        leg = ax.get_legend()
        handles, texts = list(leg.legend_handles), [t.get_text() for t in leg.get_texts()]
        lines = []
        for ln in ax.get_lines():
            xy = np.asarray(ln.get_xydata(), dtype=float)
            if xy.size == 0:
                continue  # seaborn's empty legend proxy
            hits = [texts[i] for i, h in enumerate(handles)
                    if (h.get_marker(), h.get_linestyle()) == (ln.get_marker(), ln.get_linestyle())]
            if len(hits) != 1:
                fail("panel %d: a data line matches %d legend entries" % (pi, len(hits)))
            lines.append({"legend": hits[0], "x": xy[:, 0].tolist(), "y": xy[:, 1].tolist()})
        panels.append({"panel": pi + 1, "title": ax.get_title(), "legend_title": leg.get_title().get_text(),
                       "lines": lines})
    return panels


def check(panels, rows, name):
    """Each plotted table value must be drawn at its place: panel 1 (3 clusters; x = n, line = d),
    panel 2 (3 clusters; x = d, line = n), panel 3 (n = 200; x = clusters, line = d)."""
    expected = []
    for r in rows:
        d, n, k, v = r["d"], r["n"], r["clusters"], float(r["ari_printed"])
        if k == 3:
            expected += [(0, str(d), n, v), (1, str(n), d, v)]
        if n == 200:
            expected.append((2, str(d), k, v))
    drawn = {(pi, ln["legend"], x): y for pi, p in enumerate(panels) for ln in p["lines"]
             for x, y in zip(ln["x"], ln["y"])}
    if len(drawn) != len(expected):
        fail("%s: %d points drawn, %d expected" % (name, len(drawn), len(expected)))
    for pi, lab, x, v in expected:
        y = drawn.get((pi, lab, float(x)))
        if y is None or abs(y - v) > 1e-9:
            fail("%s panel %d line %s x=%s: table %s, drawn %s" % (name, pi + 1, lab, x, v, y))
    return len(expected)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    tex_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)
    tex = open(tex_path, encoding="utf-8").read()
    all_rows, dump = [], {}
    for label, method, dist, figsize, suptitle in FIGURES:
        rows = parse_table(tex, label, method, dist)
        name = "Lineplots_Clustering_%s_%s" % (method, dist)
        plt.close("all")
        fig = draw(notebook_data(rows), figsize, suptitle, os.path.join(out_dir, name + ".png"))
        panels = read_back(fig)
        npts = check(panels, rows, name)
        plt.close(fig)
        dump[name] = {"table": label, "suptitle": fig.get_suptitle(), "panels": panels}
        all_rows += rows
        print("%s: %s, 60 ARI values, %d plotted points equal the table" % (name, label, npts))
    with open(os.path.join(out_dir, "lineplot_data.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["method", "dist", "n", "clusters", "d", "ari_printed"])
        w.writeheader()
        w.writerows(all_rows)
    with open(os.path.join(out_dir, "plotted_points.json"), "w", encoding="utf-8") as fh:
        json.dump(dump, fh, indent=1)
    print("wrote 6 figures, lineplot_data.csv (%d rows), plotted_points.json to %s" % (len(all_rows), out_dir))


if __name__ == "__main__":
    main()
