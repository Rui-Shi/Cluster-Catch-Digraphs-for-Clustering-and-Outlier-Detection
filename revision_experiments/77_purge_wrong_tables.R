#!/usr/bin/env Rscript
# 77_purge_wrong_tables.R -- remove every NN quantile table that is a duplicate
# or carries a level name its contents do not support.
#
# WHAT IS LEFT TO REMOVE, AND WHY IT IS WRONG
# -------------------------------------------
# The "99999" generator family (R/ccds/NN-test_quantile/<d>d_99999%.R, d = 10..35)
# all share one shape: a 0.99 block that is COMMENTED OUT, then a live 0.999
# block. The commented-out form is the repair; the files on disk predate it.
# In the unrepaired form one Monte Carlo run was saved under both names, which
# is why so many pairs are byte-identical -- and the surviving content is the
# FIRST computation, at 0.99.
#
# That is not read off the scripts, it is measured. At d = 12, 16, 18, 19 and
# 21 the shipped "999" file was read directly against independent replicate
# draws and returned 0.0100-0.0105 every time (WP2_RESULTS section 9, and
# again as 72's T3). Five measured cases, all identical pairs, all 1%. d = 30
# does NOT contradict this: it is a single file, never a duplicate, and it
# measured at its nominal 0.1% -- consistent with the repair having been re-run
# there and not for the pairs.
#
#   => an identical pair holds the 1% run, and its "999" member is mislabelled.
#
# TWO GROUPS, BOTH DELETED
#   A. Identical pairs, d = 11, 13, 14, 15, 17. Delete the "999" member. The
#      content survives byte-for-byte in the "99" twin, which is correctly
#      named, so nothing is lost.
#   B. Mislabelled singletons, d = 22-28. These were identical pairs until 70
#      removed the "99" member, leaving the "999" name on 1% content. Deleting
#      loses the table outright -- which is the right call here: no data set
#      and no simulation dimension touches d = 22-28, and the alternative,
#      renaming them to "99", would stamp a level that has been inferred at
#      these dimensions rather than measured. Asserting an unmeasured level is
#      the exact fault being corrected.
#
# NOTHING REACHABLE IS TOUCHED. Every deletion is checked against what the
# resolvers can actually ask for, over every data set dimension and the
# simulation grid. A dimension that any method could reach is refused.
#
# get_simul() stop()s on a missing file, so after this the affected dimensions
# fail loudly instead of silently returning a 1% table under a 0.1% name.
#
# Usage:
#   Rscript revision_experiments/77_purge_wrong_tables.R           # dry run
#   Rscript revision_experiments/77_purge_wrong_tables.R --apply

suppressMessages(library(here))
source(here::here("revision_experiments/harness.R"))
source(here::here("revision_experiments/wp0_mccd_methods.R"))

APPLY <- "--apply" %in% commandArgs(trailingOnly = TRUE)
NNDIR <- here::here("R/NN-test_quantile")
MANI  <- here::here("revision_experiments/DELETED_DUPLICATE_TABLES.md")
load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }

inv   <- read.csv(here::here("revision_experiments/results/tr1/dataset_inventory.csv"),
                  stringsAsFactors = FALSE)
SIM_D <- c(2, 3, 5, 10, 20, 50, 100)

# every (d, token) any method can request, over real data and simulations
reachable <- unique(do.call(rbind, lapply(sort(unique(c(inv$d, SIM_D))), function(d)
  data.frame(d = d, token = c(nn_quant_label_paper_UN(d), nn_quant_label_paper_SUN(d)),
             stringsAsFactors = FALSE))))
is_reachable <- function(d, tok) any(reachable$d == d & reachable$token == tok)

# ---- full sweep of what is on disk -----------------------------------------
fs <- list.files(NNDIR, pattern = "^NN-test-simul_[0-9]+d_[0-9]+%\\.RData$")
ds <- as.integer(sub("^NN-test-simul_([0-9]+)d_.*$", "\\1", fs))
cat(sprintf("%d NN quantile tables on disk, %d distinct dimensions\n\n", length(fs), length(unique(ds))))

cand <- list()
for (d in sort(unique(ds))) {
  g  <- fs[ds == d]
  tk <- sub("^.*_([0-9]+)%\\.RData$", "\\1", g)
  o  <- order(as.numeric(paste0("0.", tk))); g <- g[o]; tk <- tk[o]

  if (length(g) == 2) {
    a <- load1(file.path(NNDIR, g[1])); b <- load1(file.path(NNDIR, g[2]))
    same <- identical(a$average, b$average) && identical(a$median, b$median)
    rm(a, b)
    if (!same) next                                   # genuinely distinct: keep both
    hi <- which.max(as.numeric(paste0("0.", tk)))     # the "999"-ish member is the wrong one
    cand[[length(cand) + 1]] <- data.frame(
      d = d, file = g[hi], token = tk[hi], group = "A: identical pair",
      reason = sprintf("byte-identical to %s, which is correctly named", g[-hi]),
      stringsAsFactors = FALSE)
  } else if (length(g) == 1 && tk == "999" && d >= 22 && d <= 28) {
    cand[[length(cand) + 1]] <- data.frame(
      d = d, file = g, token = tk, group = "B: mislabelled singleton",
      reason = "was an identical pair until 70 removed its twin; holds the 1% run",
      stringsAsFactors = FALSE)
  }
}
if (!length(cand)) { cat("nothing to purge\n"); quit(status = 0) }
p <- do.call(rbind, cand)

p$reachable <- mapply(is_reachable, p$d, p$token)
cat("=== candidates ===\n")
print(p[, c("d", "file", "token", "group", "reachable")], row.names = FALSE)

if (any(p$reachable)) {
  cat("\n*** ABORT: a candidate is reachable by a resolver -- deleting it would break a method ***\n")
  print(p[p$reachable, ], row.names = FALSE)
  quit(status = 1)
}
cat(sprintf("\n  %d candidates, none reachable by any resolver at any study dimension\n", nrow(p)))

if (!APPLY) { cat("\nDRY RUN -- nothing deleted. Re-run with --apply.\n"); quit(status = 0) }

w <- function(...) cat(sprintf(...), file = MANI, append = TRUE)
w("\n## Purged %s by `77_purge_wrong_tables.R`\n\n", format(Sys.Date()))
w("Tables whose level name their contents do not support. An identical pair\n")
w("holds the 1%% run -- measured at d = 12, 16, 18, 19, 21, which all read\n")
w("0.0100-0.0105 against independent draws -- so the \"999\" member of such a\n")
w("pair is mislabelled. None of these is reachable by either resolver at any\n")
w("data set dimension or simulation grid point; that was re-checked here.\n\n")
w("| d | deleted | token | bytes | MD5 | group | reason |\n|---|---|---|---|---|---|---|\n")

cat("\n=== deleting ===\n")
for (i in seq_len(nrow(p))) {
  f <- file.path(NNDIR, p$file[i])
  sz <- file.size(f); md5 <- unname(tools::md5sum(f))
  ok <- file.remove(f)
  cat(sprintf("  d=%-3d %-34s removed=%s\n", p$d[i], p$file[i], ok))
  w("| %d | `%s` | %s | %s | `%s` | %s | %s |\n", p$d[i], p$file[i], p$token[i],
    format(sz, big.mark = ","), substr(md5, 1, 16), p$group[i], p$reason[i])
}
w("\n")
cat(sprintf("\nmanifest appended: %s\nnow re-run 69_verify_table_inventory.R\n", MANI))
