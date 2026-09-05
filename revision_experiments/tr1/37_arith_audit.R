# 37_arith_audit.R -- internal-consistency audit of the final 16x9 table.
#
# Prompted by the TR2 session, which noticed that LOF/vertebral reports
# BA = 0.488 while its own TPR = 0.033 and TNR = 0.938 imply 0.486.
#
# BA is defined as (TPR + TNR)/2 everywhere in this project, so every row must
# satisfy that identity. The stored values are rounded to 3 d.p., which admits
# an error of at most 0.001 in the reconstruction; anything larger is a defect,
# not rounding.
#
# F2 cannot be checked this way -- it needs the confusion matrix, not just the
# two rates -- but it can be reconstructed from TPR, TNR and the known n, n0:
#     TP = TPR*n0,  FP = (1-TNR)*(n-n0),  F2 = 5*TP / (5*TP + 4*FN + FP)
# with FN = n0 - TP. Same 3 d.p. tolerance, loosened to 0.002 because two
# rounded inputs propagate.
#
# Read-only. Writes nothing.

suppressMessages(library(here))

cmp <- read.csv(here::here("revision_experiments", "results", "tr1",
                           "final_comparison.csv"), stringsAsFactors = FALSE)

# n and n0 as the loader returns them (BENCHMARK_EXPANSION_RULE.md).
meta <- data.frame(
  dataset = c("hepatitis","lymphography","glass","WBC","vertebral","ecoli",
              "stamps","WDBC","pima","shuffle","vowels","PenDigits",
              "waveform","thyroid","pageblocks","wilt"),
  n  = c(74,148,213,223,240,336,340,367,555,1013,1452,3200,3443,3656,4795,4819),
  n0 = c(7,6,9,10,30,8,31,10,55,13,46,20,100,93,510,257),
  stringsAsFactors = FALSE
)

cmp <- merge(cmp, meta, by = "dataset", all.x = TRUE)
if (any(is.na(cmp$n))) {
  stop("unmatched dataset(s): ",
       paste(unique(cmp$dataset[is.na(cmp$n)]), collapse = ", "))
}

cmp$BA_implied <- (cmp$TPR + cmp$TNR) / 2
cmp$BA_err     <- abs(cmp$BA - cmp$BA_implied)

TP <- cmp$TPR * cmp$n0
FN <- cmp$n0 - TP
FP <- (1 - cmp$TNR) * (cmp$n - cmp$n0)
cmp$F2_implied <- ifelse(TP == 0, 0, 5 * TP / (5 * TP + 4 * FN + FP))
cmp$F2_err     <- abs(cmp$F2 - cmp$F2_implied)

cat("rows checked:", nrow(cmp), "\n\n")

cat("=== BA inconsistent with its own TPR and TNR (tol 0.001) ===\n")
bad_ba <- cmp[cmp$BA_err > 0.0011, ]
if (nrow(bad_ba) == 0) {
  cat("  none\n")
} else {
  bad_ba <- bad_ba[order(-bad_ba$BA_err), ]
  cat(sprintf("%-13s %-9s %7s %7s %8s %8s %8s\n",
              "dataset", "method", "TPR", "TNR", "BA", "implied", "err"))
  for (i in seq_len(nrow(bad_ba))) {
    r <- bad_ba[i, ]
    cat(sprintf("%-13s %-9s %7.3f %7.3f %8.3f %8.3f %8.3f\n",
                r$dataset, r$method, r$TPR, r$TNR, r$BA, r$BA_implied, r$BA_err))
  }
}

cat("\n=== F2 inconsistent with TPR, TNR, n, n0 (tol 0.002) ===\n")
bad_f2 <- cmp[cmp$F2_err > 0.0021, ]
if (nrow(bad_f2) == 0) {
  cat("  none\n")
} else {
  bad_f2 <- bad_f2[order(-bad_f2$F2_err), ]
  cat(sprintf("%-13s %-9s %7s %7s %8s %8s %8s\n",
              "dataset", "method", "TPR", "TNR", "F2", "implied", "err"))
  for (i in seq_len(nrow(bad_f2))) {
    r <- bad_f2[i, ]
    cat(sprintf("%-13s %-9s %7.3f %7.3f %8.3f %8.3f %8.3f\n",
                r$dataset, r$method, r$TPR, r$TNR, r$F2, r$F2_implied, r$F2_err))
  }
}

cat("\n=== worst 8 BA residuals, whether or not flagged ===\n")
top <- cmp[order(-cmp$BA_err), ][1:8, ]
cat(sprintf("%-13s %-9s %8s %8s %8s\n",
            "dataset", "method", "BA", "implied", "err"))
for (i in seq_len(nrow(top))) {
  r <- top[i, ]
  cat(sprintf("%-13s %-9s %8.3f %8.3f %8.3f\n",
              r$dataset, r$method, r$BA, r$BA_implied, r$BA_err))
}

cat("\ndone\n")
