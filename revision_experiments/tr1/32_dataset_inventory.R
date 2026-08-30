#!/usr/bin/env Rscript
# 32_dataset_inventory.R -- structural inventory of every data set the loader
# builds, for declaring the benchmark-expansion rule BEFORE any results exist.
#
# The revision expands the real-data study from 8 sets to answer AE.3 (limited
# baselines / narrow evaluation) and AE.4 (weak high-dimensional validation).
# To keep the expansion defensible, the inclusion rule must be fixed on
# structural properties only -- n, d, contamination, whether the loader sorts
# it correctly -- and never on how any method scores. This script prints ONLY
# those structural properties. No detector is run.
#
# Read-only. Writes results/dataset_inventory.csv.

suppressMessages(library(here))
source(here::here("revision_experiments", "shared", "harness.R"))

# Every object RealData_Collection.R builds, in file order.
CANDIDATES <- c("hepatitis", "glass", "vertebral", "ecoli", "stamps",
                "vowels", "waveform", "wilt",              # the current eight
                "WBC", "WDBC", "lymphography", "pima",
                "thyroid", "shuffle", "PenDigits", "pageblocks")

IN_STUDY <- c("hepatitis", "glass", "vertebral", "ecoli",
              "stamps", "vowels", "waveform", "wilt")

rows <- list()
cat(sprintf("\n%-13s %7s %5s %6s %8s %12s %-8s\n",
            "dataset", "n", "d", "n0", "contam%", "sorted?", "status"))
cat(strrep("-", 66), "\n")

for (nm in CANDIDATES) {
  out <- tryCatch({
    dat <- load_real_dataset(nm)
    Y <- dat$Y; n <- length(Y); n0 <- sum(Y == 0)
    sorted_ok <- identical(as.numeric(Y),
                           as.numeric(c(rep(1, n - n0), rep(0, n0))))
    list(n = n, d = dat$d, n0 = n0, contam = 100 * n0 / n,
         sorted = sorted_ok, status = "ok")
  }, error = function(e) {
    list(n = NA, d = NA, n0 = NA, contam = NA, sorted = NA,
         status = substr(conditionMessage(e), 1, 30))
  })

  cat(sprintf("%-13s %7s %5s %6s %8s %12s %-8s\n",
              nm,
              ifelse(is.na(out$n), "-", format(out$n)),
              ifelse(is.na(out$d), "-", format(out$d)),
              ifelse(is.na(out$n0), "-", format(out$n0)),
              ifelse(is.na(out$contam), "-", sprintf("%.1f", out$contam)),
              ifelse(is.na(out$sorted), "-", ifelse(out$sorted, "OK", "MIS-SORTED")),
              if (nm %in% IN_STUDY) "in study" else out$status))

  rows[[length(rows) + 1]] <- data.frame(
    dataset = nm, n = out$n, d = out$d, n0 = out$n0,
    contam_pct = out$contam, sorted_ok = out$sorted,
    in_study = nm %in% IN_STUDY, load_status = out$status,
    stringsAsFactors = FALSE)
}

inv <- do.call(rbind, rows)
write.csv(inv, here::here("revision_experiments/results/tr1/dataset_inventory.csv"),
          row.names = FALSE)

cat("\nNOTE: 'MIS-SORTED' means the loader's final order(...) call uses a\n",
    "feature column rather than the label column. It does NOT affect our\n",
    "numbers -- evaluate() reorders (Y, score) jointly before counting -- but\n",
    "it invalidates any published figure that was counted positionally.\n", sep = "")
cat("\nWritten: results/tr1/dataset_inventory.csv\n")
