#!/usr/bin/env Rscript
# 75_install_regen999.R -- move verified alpha=0.1% tables from staging into the
# live table directory R/NN-test_quantile.
#
# Installing is the irreversible step. No quantile table is tracked by git
# (.gitignore line 1 is *.RData), so an overwritten file is gone; the only
# reason the last accidental deletion was recoverable is that it went through
# the Recycle Bin. Three guards, in order:
#
#   1. REFUSE ANYTHING NOT VERIFIED. A dimension is installed only if
#      72_verify_regen999.R recorded pass = TRUE for it. Installing a table on
#      the strength of its filename is precisely the defect being repaired.
#   2. NEVER DESTROY THE LAST COPY OF A SHIPPED TABLE. The file about to be
#      overwritten is the one that produced the published numbers. Where an
#      identical sibling exists -- at d = 12, 16, 18, 19 the "999" file is
#      byte-identical to its "99" twin, one MC run saved twice -- the content
#      survives in the twin and overwriting costs nothing; identity is
#      re-checked HERE, at install time, not taken from an earlier run. Where
#      no twin exists (d=21, whose "99" was removed as a duplicate, leaving
#      "999" as the sole carrier of the shipped 1% content) the file is first
#      renamed to the name that honestly describes it, "99", which is also the
#      baseline 72's T1 needs.
#   3. RECORD WHAT MOVED. A tracked manifest, since the files themselves are
#      invisible to git.
#
# Usage:
#   Rscript revision_experiments/75_install_regen999.R           # dry run
#   Rscript revision_experiments/75_install_regen999.R --apply

suppressMessages(library(here))
APPLY <- "--apply" %in% commandArgs(trailingOnly = TRUE)

STAGE <- here::here("R/NN-test_quantile_regen999")
LIVE  <- here::here("R/NN-test_quantile")
VCSV  <- here::here("revision_experiments/results/tr1/wp2c_regen999_verify.csv")
MANI  <- here::here("revision_experiments/INSTALLED_REGEN999_TABLES.md")

load1 <- function(p) { e <- new.env(); load(p, envir = e); get("simul", envir = e) }
livef <- function(d, tok) file.path(LIVE, sprintf("NN-test-simul_%dd_%s%%.RData", d, tok))

if (!file.exists(VCSV)) stop("no verification record at ", VCSV, " -- run 72 first")
v <- read.csv(VCSV, stringsAsFactors = FALSE)
ok_d <- v$d[isTRUE(v$pass) | v$pass %in% c(TRUE, "TRUE")]
cat(sprintf("dimensions verified by 72: %s\n", paste(ok_d, collapse = ", ")))

#' Newest staged table for a dimension: a completed full run writes directly
#' into <d>d/, a partial pool into <d>d/pooled_n<niter>/. Prefer the deepest.
staged <- function(d) {
  direct <- file.path(STAGE, sprintf("%dd", d), sprintf("NN-test-simul_%dd_999%%.RData", d))
  pools  <- list.files(file.path(STAGE, sprintf("%dd", d)), pattern = "^pooled_n[0-9]+$",
                       full.names = TRUE)
  cands <- character(0)
  if (file.exists(direct)) cands <- direct
  if (length(pools)) {
    depth <- as.integer(sub(".*pooled_n", "", pools))
    pf <- file.path(pools[order(-depth)], sprintf("NN-test-simul_%dd_999%%.RData", d))
    cands <- c(cands, pf[file.exists(pf)])
  }
  if (!length(cands)) return(NA_character_)
  cands[1]
}

plan <- list()
for (d in sort(ok_d)) {
  src <- staged(d)
  if (is.na(src)) { cat(sprintf("  d=%-3d no staged table found, skipped\n", d)); next }
  tgt  <- livef(d, "999"); twin <- livef(d, "99")
  have_twin <- file.exists(twin)
  identical_to_twin <- NA
  if (file.exists(tgt) && have_twin) {
    a <- load1(tgt); b <- load1(twin)
    identical_to_twin <- identical(a$average, b$average) && identical(a$median, b$median)
    rm(a, b)
  }
  action <- if (!file.exists(tgt)) "install (no existing file)"
            else if (isTRUE(identical_to_twin)) "overwrite (content survives in the 99 twin)"
            else if (!have_twin) "rename existing 999 -> 99, then install"
            else "*** HALT: existing 999 differs from 99 and would be lost ***"
  plan[[length(plan) + 1]] <- data.frame(
    d = d, src = src, tgt = tgt, twin_exists = have_twin,
    identical_to_twin = identical_to_twin, action = action,
    new_len = length(load1(src)$average),
    new_niter = { m <- attr(load1(src), "nn_multiquant_meta"); if (is.null(m$niter)) NA else m$niter },
    stringsAsFactors = FALSE)
  cat(sprintf("  d=%-3d %s\n        src   %s\n        niter %s, length %d\n",
              d, action, sub(".*regen999.", "", src),
              plan[[length(plan)]]$new_niter, plan[[length(plan)]]$new_len))
}
p <- do.call(rbind, plan)
if (is.null(p)) { cat("nothing to install\n"); quit(status = 0) }
if (any(grepl("HALT", p$action))) { cat("\n*** ABORT: an existing table would be destroyed ***\n"); quit(status = 1) }

if (!APPLY) { cat("\nDRY RUN -- nothing written. Re-run with --apply.\n"); quit(status = 0) }

w <- function(...) cat(sprintf(...), file = MANI, append = TRUE)
if (!file.exists(MANI)) {
  cat("", file = MANI)
  w("# Installed alpha=0.1%% NND quantile tables\n\n")
  w("Replacements for the \"999\" tables that were byte-identical to their \"99\"\n")
  w("siblings and held 1%%, not the 0.1%% the manuscript claims for SUN-MCCD at\n")
  w("d >= 10 (WP2_RESULTS section 9). Generated by `71_regen_nnd_alpha001.R`,\n")
  w("accepted by `72_verify_regen999.R`, installed by this script. `*.RData` is\n")
  w("gitignored, so this file is the only record in version control.\n\n")
  w("| date | d | niter | length | replaced | disposition of the old file |\n")
  w("|---|---|---|---|---|---|\n")
}

cat("\n=== installing ===\n")
for (i in seq_len(nrow(p))) {
  d <- p$d[i]; tgt <- p$tgt[i]; twin <- livef(d, "99")
  disp <- "no previous file"
  if (file.exists(tgt)) {
    if (isTRUE(p$identical_to_twin[i])) {
      disp <- sprintf("overwritten; identical content retained as `%s`", basename(twin))
    } else {
      file.rename(tgt, twin)
      disp <- sprintf("renamed to `%s` (its true level)", basename(twin))
    }
  }
  file.copy(p$src[i], tgt, overwrite = TRUE)
  chk <- load1(tgt)
  cat(sprintf("  d=%-3d installed, length %d, min nonzero %.6g   [%s]\n",
              d, length(chk$average), min(chk$average[chk$average > 0]), disp))
  w("| %s | %d | %s | %d | `%s` | %s |\n", format(Sys.Date()), d, p$new_niter[i],
    length(chk$average), basename(tgt), disp)
}
cat(sprintf("\nmanifest: %s\n", MANI))
cat("now re-run 69_verify_table_inventory.R\n")
