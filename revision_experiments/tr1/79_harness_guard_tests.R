#!/usr/bin/env Rscript
# revision_experiments/tr1/79_harness_guard_tests.R
#
# Guard tests for the 2026-09-05 harness audit fixes. Every assertion here
# corresponds to a defect that was found in the harness and silently changed
# or corrupted a number rather than failing:
#
#   1. Two competing alpha schedules in one session (harness rk_quant_for_d /
#      nn_quant_for_d vs. wp0's paper resolvers), disagreeing at d = 10.
#   2. A quantile table shorter than the data set it serves, padded with NA by
#      nnccd.radi() rather than refused.
#   3. NA or short score vectors, which count_scores2() counts as DETECTED
#      OUTLIERS because its predictions initialise to 0 = outlier.
#   4. append_result() writing positionally against a header it never checks.
#   5. has_result() reading a truncated line as a complete row of NAs, so a
#      restart skips a cell that was never actually computed.
#
# Read-only apart from one throwaway CSV in the session scratchpad. Runs in
# seconds; no detector is executed here (see 03_smoke_test.R / 13_wp0_gate.R
# for that).

suppressMessages(library(here))
source(here::here("revision_experiments/shared/harness.R"))
source(here::here("revision_experiments/tr1/wp0_mccd_methods.R"))

SCRATCH <- "C:/Users/shiru/AppData/Local/Temp/claude/G--Submissions-TR1-TR1-Neurocomputing-resubmit/672ac85e-3720-4f2e-bf2c-19ec764b119f/scratchpad"

PASS <- 0L; FAIL <- 0L
ok <- function(label, cond) {
  cond <- isTRUE(cond)
  if (cond) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  [%s] %s\n", if (cond) "PASS" else "FAIL", label))
}
#' TRUE if expr stops, and (when `pattern` is given) the message matches it.
stops_with <- function(expr, pattern = NULL) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  if (is.null(e)) return(FALSE)
  if (is.null(pattern)) return(TRUE)
  grepl(pattern, conditionMessage(e))
}

cat("=== 1. the deleted buckets are really gone ===\n")
ok("rk_quant_for_d() does not exist", !exists("rk_quant_for_d"))
ok("nn_quant_for_d() does not exist", !exists("nn_quant_for_d"))
ok("the three paper resolvers exist",
   all(sapply(c("rk_quant_label_paper", "nn_quant_label_paper_UN",
                "nn_quant_label_paper_SUN"), exists)))
ok("check_simul_extent() exists", exists("check_simul_extent"))

cat("\n=== 2. the resolvers reproduce the manuscript schedule ===\n")
# Authority: CCD_OutlierDetection_Neurocomputing.tex, Section "Uniform Cluster
# Settings" (RK: alpha = 1% for d < 10, 0.1% for d >= 10; UN-MCCD: alpha =
# 15/10/5/1/0.1% at d = 2,3,5,10,{20,50,100}; SUN-MCCD: as UN-MCCD but 0.1%
# already at d = 10) and Table tab:alpha_real, which tabulates exactly this
# for the sixteen real data sets:
#   d = 5-9   -> RK 1%,   UN 5%,   SUN 5%
#   d = 10-19 -> RK 0.1%, UN 1%,   SUN 0.1%
#   d >= 20   -> RK 0.1%, UN 0.1%, SUN 0.1%
# Filename tokens are the QUANTILE, i.e. alpha = 1% <-> "99", 0.1% <-> "999",
# 5% <-> "95", 10% <-> "90", 15% <-> "85".
D    <- c(2,     5,    9,    10,    12,    19,    20,    21,    30,    50,    100)
E_RK <- c("99",  "99", "99", "999", "999", "999", "999", "999", "999", "999", "999")
E_UN <- c("85",  "95", "95", "99",  "99",  "99",  "999", "999", "999", "999", "999")
E_SN <- c("85",  "95", "95", "999", "999", "999", "999", "999", "999", "999", "999")
got_rk <- vapply(D, rk_quant_label_paper,     character(1))
got_un <- vapply(D, nn_quant_label_paper_UN,  character(1))
got_sn <- vapply(D, nn_quant_label_paper_SUN, character(1))
print(data.frame(d = D, RK = got_rk, UN = got_un, SUN = got_sn,
                 exp_RK = E_RK, exp_UN = E_UN, exp_SUN = E_SN), row.names = FALSE)
ok("RK schedule matches the manuscript",  identical(got_rk, E_RK))
ok("UN-MCCD schedule matches the manuscript",  identical(got_un, E_UN))
ok("SUN-MCCD schedule matches the manuscript", identical(got_sn, E_SN))
ok("UN and SUN diverge exactly on d = 10..19",
   identical(D[got_un != got_sn], c(10, 12, 19)))

cat("\n=== 3. get_simul() refuses a table shorter than the data ===\n")
# NN-test-simul_19d_999%.RData was regenerated at exactly hepatitis's size
# (74 entries) under a generic filename. It serves hepatitis and nothing
# larger.
ok("quant is required (no silent default)",
   stops_with(get_simul("NN", 19), "`quant` is required"))
ok("NN d=19 999% loads for n = 74",
   { t <- get_simul("NN", 19, "999", n = 74); length(t$simul$average) == 74 })
ok("NN d=19 999% is refused for n = 75",
   stops_with(get_simul("NN", 19, "999", n = 75), "too short for the data"))
ok("the refusal names the extent and the requirement",
   stops_with(get_simul("NN", 19, "999", n = 75), "average = 74.*required n = 75|required n = 75(.|\n)*average = 74"))
ok("a full-length table is accepted at the same d",
   { t <- get_simul("NN", 19, "99", n = 5000); length(t$simul$average) == 5000 })

cat("\n=== 4. evaluate() refuses scores it cannot count honestly ===\n")
# count_scores2() initialises label_pred to 0 (= OUTLIER) and fills by
# which(), which drops NA. An NA score or a short score vector therefore
# lands in the outlier class with no warning at all.
Y_t <- c(1, 1, 1, 1, 0, 0)          # 4 regular, 2 outliers, outliers last
s_t <- c(0, 0, 0, 3, 3, 0)          # threshold 2 -> flags positions 4 and 5
m   <- evaluate(Y_t, s_t, 2)
cat("  hand example: "); print(round(m, 6))
# By hand: TNR = 3/4 (position 4 is a false alarm), TPR = 1/2 (position 6 is
# missed), BA = 0.625; precision = 2*0.5 / (2*0.5 + 0.25*4) = 0.5, recall =
# 0.5, F2 = 5*0.25/(4*0.5 + 0.5) = 0.5.
ok("TPR/TNR/BA/F2 on the hand-computable example",
   isTRUE(all.equal(unname(m), c(0.5, 0.75, 0.625, 0.5))))
ok("names are TPR/TNR/BA/F2", identical(names(m), c("TPR", "TNR", "BA", "F2")))
ok("stops on an NA score",
   stops_with(evaluate(Y_t, replace(s_t, 3, NA), 2), "not finite"))
ok("stops on a NaN score",
   stops_with(evaluate(Y_t, replace(s_t, 3, NaN), 2), "not finite"))
ok("the NA message names the count and position",
   stops_with(evaluate(Y_t, replace(s_t, 3, NA), 2), "1 of 6 score.*position\\(s\\) 3"))
ok("stops on a short score vector",
   stops_with(evaluate(Y_t, s_t[1:4], 2), "length mismatch"))
ok("the length message names both lengths",
   stops_with(evaluate(Y_t, s_t[1:4], 2), "length\\(score\\) = 4, length\\(Y\\) = 6"))

cat("\n=== 5. append_result() / has_result() survive a truncated file ===\n")
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(SCRATCH, "79_guard_test.csv")
if (file.exists(csv)) unlink(csv)
append_result(csv, list(dataset = "WBC", method = "LOF", BA = 0.9, F2 = 0.8))
append_result(csv, list(dataset = "ecoli", method = "LOF", BA = 0.7, F2 = 0.6))
# Simulate an append cut off mid-line by a drive drop: 3 fields, not 4, and
# no trailing newline.
cat("glass,MST,0.7", file = csv, append = TRUE)
cat("  file on disk:\n"); cat(paste0("    ", readLines(csv, warn = FALSE), collapse = "\n"), "\n")

warn <- NULL
r_good <- withCallingHandlers(
  has_result(csv, list(dataset = "WBC", method = "LOF")),
  warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
r_trunc <- suppressWarnings(has_result(csv, list(dataset = "glass", method = "MST")))
r_absent <- suppressWarnings(has_result(csv, list(dataset = "pima", method = "LOF")))
ok("complete row is found",            isTRUE(r_good))
ok("truncated row is treated as absent", identical(r_trunc, FALSE))
ok("absent key is absent",             identical(r_absent, FALSE))
ok("all three answers are strict TRUE/FALSE, never NA",
   all(vapply(list(r_good, r_trunc, r_absent), function(x) is.logical(x) && !is.na(x), logical(1))))
ok("a warning names the truncated line number",
   length(warn) > 0 && any(grepl("line\\(s\\) 4", warn)))
cat("  warning: ", if (length(warn)) warn[1] else "(none)", "\n", sep = "")

# A row whose keys match but whose payload is NA is a partial write, not a
# result: NA_real_ makes as.data.frame produce a real NA in the BA column.
csv2 <- file.path(SCRATCH, "79_guard_test_partial.csv")
if (file.exists(csv2)) unlink(csv2)
append_result(csv2, list(dataset = "WBC", method = "LOF", BA = 0.9, F2 = 0.8))
append_result(csv2, list(dataset = "pima", method = "LOF", BA = NA_real_, F2 = NA_real_))
ok("key-matching row with an NA payload is treated as absent",
   identical(suppressWarnings(has_result(csv2, list(dataset = "pima", method = "LOF"))), FALSE))

cat("\n=== 6. append_result() refuses a row whose names do not match ===\n")
ok("reordered names are refused",
   stops_with(append_result(csv2, list(method = "MST", dataset = "glass", BA = 0.5, F2 = 0.4)),
              "do not match the existing header"))
ok("a dropped column is refused",
   stops_with(append_result(csv2, list(dataset = "glass", method = "MST", BA = 0.5)),
              "do not match the existing header"))
ok("an added column is refused",
   stops_with(append_result(csv2, list(dataset = "glass", method = "MST", BA = 0.5, F2 = 0.4, seed = 1)),
              "do not match the existing header"))
ok("a correctly-ordered row still appends",
   { append_result(csv2, list(dataset = "glass", method = "MST", BA = 0.5, F2 = 0.4))
     isTRUE(suppressWarnings(has_result(csv2, list(dataset = "glass", method = "MST")))) })

unlink(c(csv, csv2))

cat(sprintf("\n=== %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
