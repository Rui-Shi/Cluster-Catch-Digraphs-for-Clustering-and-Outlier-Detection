# Verify Ceyhan's new manuscript sentence against the implementation:
#   "If a tie block reaches either endpoint of the ranked sample, we retain
#    the tied score and use the standard midrank convention."
#
# Sources the real file so we test shipped code, not a transcription.
SRC <- here::here("methods", "outlyingness_scores", "Outlyingness_Score.R")
source(SRC)

report <- function(label, scores, vd) {
  out <- break_ties_eq10(scores, vd)
  changed <- !isTRUE(all.equal(scores, out))
  cat(sprintf("\n--- %s ---\n", label))
  cat("  in  :", paste(sprintf("%.4f", scores), collapse = "  "), "\n")
  cat("  out :", paste(sprintf("%.4f", out),    collapse = "  "), "\n")
  cat("  tie retained? ", if (changed) "NO  <-- tie was BROKEN" else "YES", "\n")
  invisible(out)
}

cat("========================================================\n")
cat(" Does the code retain tied scores at the endpoints?\n")
cat("========================================================\n")

# 1. block at the BOTTOM of the ranked sample (i == 1)
report("BOTTOM endpoint block: three points tied at 5, then 10, 20",
       c(5, 5, 5, 10, 20), c(1, 2, 3, 1, 1))

# 2. block at the TOP of the ranked sample (j == n)
report("TOP endpoint block: 1, 2, then three points tied at 7",
       c(1, 2, 7, 7, 7), c(1, 1, 1, 2, 3))

# 3. interior block (control - manuscript's stated case)
report("INTERIOR block (control): 1, then three tied at 5, then 10",
       c(1, 5, 5, 5, 10), c(1, 2, 3, 1, 1))

# 4. whole vector tied (the genuinely degenerate case)
report("WHOLE vector tied at 4",
       c(4, 4, 4, 4), c(1, 2, 3, 4))

# ---- consequence for a thresholded decision -------------------------------
cat("\n\n========================================================\n")
cat(" Consequence: can an endpoint block cross a cutoff?\n")
cat("========================================================\n")

# top block tied at 7, cutoff 6.5: all three exceed it before tie-breaking
s <- c(1, 2, 7, 7, 7); v <- c(1, 1, 1, 2, 3); cutoff <- 6.5
o <- break_ties_eq10(s, v)
cat(sprintf("\n  cutoff = %.2f\n", cutoff))
cat("  flagged BEFORE tie-break :", sum(s > cutoff), "of 3 tied points\n")
cat("  flagged AFTER  tie-break :", sum(o[3:5] > cutoff), "of 3 tied points\n")
cat("  scores after:", paste(sprintf("%.4f", o[3:5]), collapse = "  "), "\n")

# bottom block tied at 5, cutoff 5.5
s2 <- c(5, 5, 5, 10, 20); v2 <- c(1, 2, 3, 1, 1); cut2 <- 5.5
o2 <- break_ties_eq10(s2, v2)
cat(sprintf("\n  cutoff = %.2f\n", cut2))
cat("  flagged BEFORE tie-break :", sum(s2[1:3] > cut2), "of 3 tied points\n")
cat("  flagged AFTER  tie-break :", sum(o2[1:3] > cut2), "of 3 tied points\n")
cat("  scores after:", paste(sprintf("%.4f", o2[1:3]), collapse = "  "), "\n")

cat("\n\n--- what the code actually does at an endpoint (L208-209) ---\n")
cat("  lo = if(i>1) s[i-1] else s[i]   # falls back to the block's OWN tied value\n")
cat("  hi = if(j<n) s[j+1] else s[j]   # falls back to the block's OWN tied value\n")
cat("  guard: rewrite iff is.finite(lo) && is.finite(hi) && hi>lo && den_sum>0\n")
cat("  At an endpoint the bracket is half-anchored at the tied value t,\n")
cat("  so hi>lo still holds and the block IS rewritten.\n")
cat("  Only a FULLY tied vector gives lo==hi and is left alone.\n")
