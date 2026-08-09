# 39_verify_manuscript_tables.R -- does the manuscript actually agree with the
# result files?
#
# The two supplement tables and the main-text aggregate table were written by
# hand from 35_print_tables.R output. Hand transcription is exactly the failure
# mode this revision is correcting, so it should not be trusted on its own.
# This parses the LaTeX back out and diffs it against final_comparison.csv.
#
# Read-only. Reads the manuscript repo; writes nothing.

suppressMessages(library(here))

ROOT <- dirname(here::here())          # manuscript repo root
SUP  <- file.path(ROOT, "SupplementaryMaterial.tex")
MAIN <- file.path(ROOT, "CCD_OutlierDetection_Neurocomputing.tex")

cmp <- read.csv(here::here("revision_experiments/results/tr1/final_comparison.csv"),
                stringsAsFactors = FALSE)

M  <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD", "LOF", "DBSCAN", "MST", "ODIN", "iForest")
DS <- c("hepatitis", "lymphography", "glass", "WBC", "vertebral", "ecoli",
        "stamps", "WDBC", "pima", "Shuttle", "vowels", "PenDigits",
        "waveform", "thyroid", "pageblocks", "wilt")
# the manuscript calls it Shuttle; the pipeline calls it shuffle
key <- function(ds) if (ds == "Shuttle") "shuffle" else ds

get <- function(ds, m, metric) {
  r <- cmp[cmp$dataset == key(ds) & cmp$method == m, ]
  if (!nrow(r)) return(NA_real_)
  as.numeric(r[[metric]][1])
}

lines_sup <- readLines(SUP, warn = FALSE)

# Pull the numeric cells out of a "  <dataset> & a & b & ... \\ \hline" row.
row_nums <- function(lns, ds) {
  pat <- paste0("^\\s*", ds, "\\s+&")
  hit <- grep(pat, lns)
  lapply(hit, function(i) {
    body <- sub("\\\\\\\\.*$", "", lns[i])
    body <- sub(paste0("^\\s*", ds, "\\s*&"), "", body)
    v <- suppressWarnings(as.numeric(trimws(strsplit(body, "&", fixed = TRUE)[[1]])))
    v[!is.na(v)]
  })
}

bad <- 0; checked <- 0

for (tbl in 1:2) {
  metrics <- if (tbl == 1) c("TPR", "TNR") else c("BA", "F2")
  cat(sprintf("\n=== supplement table %d (%s / %s) ===\n", tbl, metrics[1], metrics[2]))
  for (ds in DS) {
    cand <- Filter(function(v) length(v) == 18, row_nums(lines_sup, ds))
    if (length(cand) < tbl) { cat(sprintf("  %-13s NO ROW FOUND\n", ds)); bad <- bad + 1; next }
    v <- cand[[tbl]]
    for (j in seq_along(M)) {
      for (k in 1:2) {
        want <- get(ds, M[j], metrics[k])
        got  <- v[(j - 1) * 2 + k]
        checked <- checked + 1
        if (is.na(want) || abs(want - got) > 5e-4) {
          bad <- bad + 1
          cat(sprintf("  MISMATCH %-13s %-9s %-3s  tex=%.3f  csv=%s\n",
                      ds, M[j], metrics[k], got,
                      if (is.na(want)) "NA" else sprintf("%.3f", want)))
        }
      }
    }
  }
}
cat(sprintf("\n  supplement: %d cells checked, %d mismatches\n", checked, bad))

# ---- main-text aggregate table -------------------------------------------
cat("\n=== main text: aggregate table ===\n")
lines_main <- readLines(MAIN, warn = FALSE)
abad <- 0
for (m in M) {
  want <- c(mean(sapply(DS, function(d) get(d, m, "F2"))),
            median(sapply(DS, function(d) get(d, m, "F2"))),
            mean(sapply(DS, function(d) get(d, m, "BA"))),
            median(sapply(DS, function(d) get(d, m, "BA"))))
  pat <- paste0("^\\s*(\\\\textbf\\{)?", m, "(\\})?\\s+&")
  hit <- grep(pat, lines_main)
  hit <- hit[sapply(hit, function(i) length(gregexpr("&", lines_main[i])[[1]]) == 4)]
  if (!length(hit)) { cat(sprintf("  %-9s NO ROW FOUND\n", m)); abad <- abad + 1; next }
  body <- sub("\\\\\\\\.*$", "", lines_main[hit[1]])
  got <- suppressWarnings(as.numeric(gsub("[^0-9.]", "",
           trimws(strsplit(body, "&", fixed = TRUE)[[1]])[-1])))
  # Tolerance is a hair over half an ulp of the last printed digit. Several
  # of these aggregates land exactly on a rounding tie -- SU-MCCD's median BA
  # is 0.69250 and SUN-MCCD's median F2 is 0.32850 -- so a flat 5e-4 flags
  # correctly-rounded values.
  d <- abs(got - want)
  cat(sprintf("  %-9s tex %s | csv %s | %s\n", m,
              paste(sprintf("%.3f", got),  collapse = " "),
              paste(sprintf("%.3f", want), collapse = " "),
              if (all(d <= 5.001e-4)) "ok" else "MISMATCH"))
  if (!all(d <= 5.001e-4)) abad <- abad + 1
}
cat(sprintf("\n  aggregate table: %d mismatched rows\n", abad))

# ---- derived claims -------------------------------------------------------
cat("\n=== derived claims ===\n")
wins <- sum(sapply(DS, function(ds) {
  v <- sapply(M, function(m) get(ds, m, "F2"))
  names(which.max(v)) %in% M[1:4]
}))
dims <- c(hepatitis=19, lymphography=18, glass=9, WBC=9, vertebral=6, ecoli=7,
          stamps=9, WDBC=30, pima=8, Shuttle=9, vowels=12, PenDigits=16,
          waveform=21, thyroid=6, pageblocks=10, wilt=5)
low <- sum(sapply(names(dims)[dims <= 10], function(ds) {
  v <- sapply(M, function(m) get(ds, m, "F2")); names(which.max(v)) %in% M[1:4] }))
cat(sprintf("  proposed win outright        : %d of 16   (text says seven)\n", wins))
cat(sprintf("  proposed win among d <= 10   : %d of 10   (text says five)\n", low))
cat(sprintf("  proposed win among d >= 12   : %d of 6    (text says two)\n", wins - low))

cat(sprintf("\nTOTAL PROBLEMS: %d\n", bad + abad))
