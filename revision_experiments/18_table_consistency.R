# 18_table_consistency.R -- audit the two supplement result tables for internal
# inconsistency. BA is defined as (TPR + TNR)/2, so every cell across the two
# tables is over-determined and any typo in one of the three shows up as a
# mismatch. Three TNR typos were already found by hand in the SUN-MCCD row;
# this checks all 72 (method, dataset) pairs so none are missed.
#
# Reads the .tex directly. Writes nothing.

SM <- "G:/Submissions/TR1/TR1_Neurocomputing_resubmit/SupplementaryMaterial.tex"
tex <- readLines(SM, warn = FALSE)

DATASETS <- c("hepatitis", "glass", "vertebral", "ecoli",
              "stamps", "vowels", "waveform", "wilt")
METHODS  <- c("U-MCCD", "SU-MCCD", "UN-MCCD", "SUN-MCCD",
              "LOF", "DBSCAN", "MST", "ODIN", "iForest")

# grab the numeric cells of a given method's row inside a given line range
row_vals <- function(lines, meth) {
  trimmed <- trimws(lines)
  hit <- lines[startsWith(trimmed, paste0(meth, " ")) |
               startsWith(trimmed, paste0(meth, "\t"))]
  if (!length(hit)) return(NULL)
  nums <- regmatches(hit[1], gregexpr("[0-9]+\\.[0-9]{3}", hit[1]))[[1]]
  as.numeric(nums)
}

t1 <- tex[554:572]   # TPR / TNR
t2 <- tex[574:591]   # BA / F2

cat(sprintf("\n%-9s %-10s %7s %7s %7s %8s %8s\n",
            "method", "dataset", "TPR", "TNR", "BA", "(T+N)/2", "delta"))
bad <- 0
for (m in METHODS) {
  a <- row_vals(t1, m); b <- row_vals(t2, m)
  if (is.null(a) || is.null(b) || length(a) < 16 || length(b) < 16) {
    cat(sprintf("%-9s  PARSE FAILED (got %d and %d values)\n", m,
                length(a), length(b))); next
  }
  for (i in seq_along(DATASETS)) {
    tpr <- a[2 * i - 1]; tnr <- a[2 * i]; ba <- b[2 * i - 1]
    calc <- (tpr + tnr) / 2
    d <- abs(round(calc, 3) - ba)
    if (d > 0.0005) {
      bad <- bad + 1
      cat(sprintf("%-9s %-10s %7.3f %7.3f %7.3f %8.4f %8.4f  <<< INCONSISTENT\n",
                  m, DATASETS[i], tpr, tnr, ba, calc, d))
    }
  }
}
cat(sprintf("\n%d inconsistent (method, dataset) pairs out of %d\n",
            bad, length(METHODS) * length(DATASETS)))
