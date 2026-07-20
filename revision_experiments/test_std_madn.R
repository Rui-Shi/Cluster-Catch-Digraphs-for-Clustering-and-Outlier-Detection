# test_std_madn.R -- regression test for the std_MADN degenerate-case fallback
#
# Manuscript convention (CCDwScores.tex line 625, Ceyhan's first-cycle item
# B3): standardize by (x - Med)/MADN; "In the rare degenerate case where
# MADN = 0 ... we replace the denominator with the sample standard deviation
# SD; if both MADN and SD vanish ... we set IOS^std = 0 for every x in C."
# So the fallback chain is MADN -> SD -> 0, and the numerator (centering on
# the median) is unchanged in the SD branch (only the DENOMINATOR is replaced).
#
# Bug (found 2026-07-19, WP5 full-data Musk): the code
# (Outlyingness_Score.R:56-63) implemented only the two ENDS of the chain --
# MADN != 0 -> (x-Med)/MADN, else -> all zeros -- with no SD step. The two
# rules agree unless MADN = 0 but SD > 0 (>=50% of values tie at the median
# while the rest vary). That regime shows up on real high-d data: full-data
# Musk raw OOS had 65% of points tied at the median (MADN = 0) yet SD = 17.25,
# so the code zeroed all 3062 scores instead of SD-standardizing them.
#
# This test pins all four branches. It FAILS on the pre-fix code at the
# MADN=0 & SD>0 case (returns zeros; manuscript says SD-scale).
#
# Run: Rscript "revision_experiments/test_std_madn.R"

source(here::here("methods/outlyingness_scores/Outlyingness_Score.R"))

fail <- 0
check <- function(ok, label) {
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) fail <<- fail + 1
}

# --- (1) MADN = 0 but SD > 0: replace denominator with SD, keep median center
# x: six values tied at the median (1), four spread above -> mad(x) = 0,
# sd(x) > 0. Manuscript: result = (x - median(x)) / sd(x).
x1 <- c(1, 1, 1, 1, 1, 1, 2, 3, 4, 5)
stopifnot(mad(x1) == 0, sd(x1) > 0)                 # precondition: the divergence regime
got1 <- std_MADN(x1)
exp1 <- (x1 - median(x1)) / sd(x1)                  # median-centered, SD-scaled
check(isTRUE(all.equal(got1, exp1)),
      "MADN=0 & SD>0 -> (x - median)/SD (manuscript B3 SD fallback)")
check(any(got1 != 0),
      "MADN=0 & SD>0 -> NOT all zeros (the pre-fix bug)")

# --- (2) MADN > 0: unchanged behavior, exactly (x - median)/MADN -----------
x2 <- as.numeric(1:11)                              # mad > 0
got2 <- std_MADN(x2)
exp2 <- (x2 - median(x2)) / mad(x2)
check(isTRUE(all.equal(got2, exp2)),
      "MADN>0 -> (x - median)/MADN (backward-compatible)")

# --- (3) both MADN = 0 and SD = 0 (constant): all zeros --------------------
x3 <- rep(3, 5)
got3 <- std_MADN(x3)
check(all(got3 == 0) && length(got3) == 5,
      "constant vector -> all zeros")

# --- (4) singleton cluster: mad = 0, sd = NA -> all zeros (length 1) --------
x4 <- 7
got4 <- std_MADN(x4)
check(length(got4) == 1 && got4 == 0,
      "singleton -> 0 (sd = NA treated as degenerate)")

cat(if (fail == 0) "\nALL PASS\n" else sprintf("\n%d FAILURE(S)\n", fail))
if (fail > 0) quit(status = 1)
