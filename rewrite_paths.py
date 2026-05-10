"""Comprehensive path rewriter (re-run May 2026).

Rewrites every absolute-path call across all .R files to a portable
here::here("...") form anchored at the repo root. Handles:

    source(), load(), setwd(), save.image(), save()

Run once after directory reorganization. Maps machine-specific prefixes
(/mmfs1/, /media/rui/exNVME/, /media/shirui001/..., G:/code_working_folder/,
OneDrive variants) to the new repo layout, including historical method
aliases (M-FCCDs <- SU-MCCDs, M-FNNCCDs <- SUN-MCCD, M-CCDs <- RU-MCCDs).
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Machine-specific prefixes observed in absolute path calls.
# Order matters: longer prefixes first to avoid partial matches.
PREFIXES = [
    "/mmfs1/home/rzs0112/code_working_folder/",
    "/mmfs1/home/rzs0112/exNVME/code_working_folder/",
    "/media/rui/exNVME/code_working_folder/",
    "/media/shirui001/Documents/share/code_working folder/",
    "/media/shirui001/WD SN350/code_working folder/",
    "G:/code_working_folder/",
    "G:/OneDrive - Auburn University/Research Outliers Detection/code/",
    "G:/OneDrive - Auburn University/Research_Outliers Detection/code/",
    "C:/Users/shiru/OneDrive - Auburn University/Research Outliers Detection/code/",
]

# Map "logical path under code_working_folder/" -> repo-relative path.
# Order: longest/most-specific first.
RULES = [
    # Specific top-level driver scripts
    ("CCDs_Clustering/RK-CCDs.R",            "methods/clustering/RK-CCDs.R"),
    ("CCDs_Clustering/UN-CCDs.R",            "methods/clustering/UN-CCDs.R"),
    ("CCDs_Clustering/KS-CCDs.R",            "methods/clustering/KS-CCDs.R"),
    ("Outlyingness_Score/RKCCD_OOS_IOS.R",   "methods/outlyingness_scores/RKCCD_OOS_IOS.R"),
    ("Outlyingness_Score/UNCCD_OOS_IOS.R",   "methods/outlyingness_scores/UNCCD_OOS_IOS.R"),
    ("Outlyingness_Score/Outlyingness_Score.R", "methods/outlyingness_scores/Outlyingness_Score.R"),
    # Real-data collection scripts
    ("Algo_Compare_Clustering/Real_Datasets/Real_Data_Collection.R", "data/clustering/Real_Data_Collection.R"),
    ("Algo_Compare/Real Datasets/RealData_Collection.R",             "data/outlier_detection/RealData_Collection.R"),
    ("Algo_Compare_OD/Real_Datasets",                                "data/outlier_detection"),
    # Quantile-table subfolders moved to R/ root (NOT under R/general_functions/)
    ("general functions/NN-test_quantile/",  "R/NN-test_quantile/"),
    ("general functions/RK-test_quantile/",  "R/RK-test_quantile/"),
    ("general_functions/NN-test_quantile/",  "R/NN-test_quantile/"),
    ("general_functions/RK-test_quantile/",  "R/RK-test_quantile/"),
    # Library code
    ("ccds/",                       "R/ccds/"),
    ("general functions/",          "R/general_functions/"),
    ("general_functions/",          "R/general_functions/"),
    # Driver-script directories (catch-all for OS/Outlyingness_Score that didn't match above)
    ("Outlyingness_Score/",         "methods/outlyingness_scores/"),
    # Historical aliases for outlier-detection methods (from arXiv:2409.11596 era)
    ("M-FCCDs/M-FCCDs.R",          "methods/outlier_detection/SU-MCCDs.R"),
    ("M-FCCDs/",                   "simulations/outlier_detection/SU-MCCDs/"),
    ("M-FNNCCDs/M-FNNCCDs.R",      "methods/outlier_detection/SUN-MCCD.R"),
    ("M-FNNCCDs/",                 "simulations/outlier_detection/SUN-MCCDs/"),
    ("M-CCDs/M-CCDs.R",            "methods/outlier_detection/RU-MCCDs.R"),
    ("M-CCDs/",                    "simulations/outlier_detection/RU-MCCDs/"),
    # M-NNCCDs is the historical name of UN-MCCDs (NN-based, non-shape-adaptive,
    # mirroring the M-FNNCCDs <- SUN-MCCD pairing).
    ("M-NNCCDs/M-NNCCDs.R",        "methods/outlier_detection/UN-MCCD.R"),
    ("M-NNCCDs/",                  "simulations/outlier_detection/UN-MCCDs/"),
    # Lowercase variant of UNCCD_OOS_IOS — observed once
    ("UnCCD_OOS_IOS/",             "simulations/outlyingness_scores/UNCCD_OOS_IOS/"),
    # Per-method directories that may appear in source()/setwd() calls
    ("RU-MCCDs/RU-MCCDs.R",        "methods/outlier_detection/RU-MCCDs.R"),
    ("SU-MCCDs/SU-MCCDs.R",        "methods/outlier_detection/SU-MCCDs.R"),
    ("SUN-MCCDs/SUN-MCCD.R",       "methods/outlier_detection/SUN-MCCD.R"),
    ("UN-MCCDs/UN-MCCD.R",         "methods/outlier_detection/UN-MCCD.R"),
    ("RU-MCCDs/",                  "simulations/outlier_detection/RU-MCCDs/"),
    ("SU-MCCDs/",                  "simulations/outlier_detection/SU-MCCDs/"),
    ("SUN-MCCDs/",                 "simulations/outlier_detection/SUN-MCCDs/"),
    ("UN-MCCDs/",                  "simulations/outlier_detection/UN-MCCDs/"),
    # Simulation-experiment trees (used by setwd()/load()/save.image())
    ("CCDs_Clustering/KS_CCDs/",                    "simulations/clustering/KS_CCDs/"),
    ("CCDs_Clustering/RK_CCDs/",                    "simulations/clustering/RK_CCDs/"),
    ("CCDs_Clustering/UN_CCDs/",                    "simulations/clustering/UN_CCDs/"),
    ("CCDs_Clustering/Algo_Compare_Clustering/",    "simulations/clustering/Algo_Compare_Clustering/"),
    ("CCDs_Outlier_Detection/",                     "simulations/outlier_detection/"),
    ("CCDs_Outlyingness_Score/RKCCD_OOS_IOS/",      "simulations/outlyingness_scores/RKCCD_OOS_IOS/"),
    ("CCDs_Outlyingness_Score/UNCCD_OOS_IOS/",      "simulations/outlyingness_scores/UNCCD_OOS_IOS/"),
    ("CCDs_Outlyingness_Score/Gaussian_cutoffs/",   "simulations/outlyingness_scores/Gaussian_cutoffs/"),
    ("CCDs_Outlyingness_Score/Uniform_cutoffs/",    "simulations/outlyingness_scores/Uniform_cutoffs/"),
    # Baseline detector trees
    ("Algo_Compare_Clustering/",    "simulations/clustering/Algo_Compare_Clustering/"),
    ("Algo_Compare/",               "simulations/outlier_detection/Algo_Compare_OutlierDetection/"),
    # Threshold-related references
    ("RKCCD_OOS_IOS/",              "simulations/outlyingness_scores/RKCCD_OOS_IOS/"),
    ("UNCCD_OOS_IOS/",              "simulations/outlyingness_scores/UNCCD_OOS_IOS/"),
    # OneDrive layout aliases for library code
    ("KS-CCDs & mKNN/functions/",   "R/ccds/"),
    ("KS-CCDs mKNN/functions/",     "R/ccds/"),
    ("Algorithm 3/functions/",      "R/ccds/"),
    ("Heuristic 1/functions/",      "R/ccds/"),
    # KS-MCGs (predecessor name of D-MCGs)
    ("KS-MCGs/simulations/",        "simulations/outlier_detection/D-MCGs/Simulations/"),
    ("KS-MCGs/",                    "simulations/outlier_detection/D-MCGs/"),
    # Unbundled companion files (still mapped to plausible repo locations)
    ("NNCCD_OOS_IOS/",              "simulations/outlyingness_scores/NNCCD_OOS_IOS/"),
]

# Match function calls with a string-literal first argument.
PATH_CALL_RX = re.compile(
    r'(?P<lead>^\s*#?\s*)'
    r'(?P<fn>source|load|setwd|save\.image|save)\s*\(\s*'
    r'(?P<q>["\'])(?P<path>[^"\']+)(?P=q)'
    r'(?P<rest>(?:\s*,[^)]*)?)\s*\)',
    re.MULTILINE,
)


def strip_prefix(path: str) -> str | None:
    for p in PREFIXES:
        if path.startswith(p):
            return path[len(p):]
    return None


def map_logical(logical: str) -> str | None:
    for old, new in RULES:
        if logical.startswith(old):
            return new + logical[len(old):]
    return None


def is_absolute(path: str) -> bool:
    return path.startswith(("/", "G:/", "C:/", "D:/", "E:/", "F:/")) or (
        len(path) > 2 and path[1] == ":" and path[2] in "\\/"
    )


def process_file(rfile: Path):
    try:
        text = rfile.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = rfile.read_text(encoding="latin-1")

    rewrites = 0
    orphans: list[tuple[int, str]] = []

    def repl(m: re.Match) -> str:
        nonlocal rewrites
        lead = m.group("lead")
        fn = m.group("fn")
        path = m.group("path")
        rest = m.group("rest") or ""
        if not is_absolute(path):
            return m.group(0)
        logical = strip_prefix(path)
        if logical is None:
            line_no = text[: m.start()].count("\n") + 1
            orphans.append((line_no, path))
            return m.group(0)
        new = map_logical(logical)
        if new is None:
            line_no = text[: m.start()].count("\n") + 1
            orphans.append((line_no, path))
            return m.group(0)
        rewrites += 1
        return f'{lead}{fn}(here::here("{new}"){rest})'

    new_text = PATH_CALL_RX.sub(repl, text)
    if new_text != text:
        rfile.write_text(new_text, encoding="utf-8")
    return rewrites, orphans


def main() -> None:
    rfiles = [p for p in ROOT.rglob("*.R") if ".git" not in p.parts]
    print(f"Scanning {len(rfiles)} R files")
    total_rewrites = 0
    all_orphans: list[tuple[Path, int, str]] = []
    files_modified = 0
    for f in rfiles:
        before = total_rewrites
        r, o = process_file(f)
        total_rewrites += r
        if r:
            files_modified += 1
        for line_no, p in o:
            all_orphans.append((f, line_no, p))

    print(f"Rewrites:           {total_rewrites}")
    print(f"Files modified:     {files_modified}")
    print(f"Orphans:            {len(all_orphans)}")

    if all_orphans:
        by_logical: dict[str, int] = {}
        for f, line, path in all_orphans:
            logical = strip_prefix(path) or path
            by_logical[logical] = by_logical.get(logical, 0) + 1
        print("\n--- distinct orphan paths (no rule matched) ---")
        for logical, count in sorted(by_logical.items(), key=lambda x: -x[1]):
            print(f"  {count:5d}x  {logical}")


if __name__ == "__main__":
    main()
