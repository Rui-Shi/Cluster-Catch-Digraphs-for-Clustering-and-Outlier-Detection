#!/usr/bin/env Rscript
# revision_experiments/tr3/12_real_reference.R
#
# Full-precision reference values for the real-data gate, from the ORIGINAL
# code: each requested "##### <name> data #####" block of
# simulations/clustering/<M>_CCDs/Simulations/Real_Data/Real_Data_<M>CCD.R is
# evaluated verbatim, in .GlobalEnv as the script ran it. The paper prints
# these values to 3 decimals and only the KS-CCD knit output
# (Real_Data_KSCCD.html) keeps all digits, so 20_gate.R checks 11_real.R
# against this file to every digit and against the paper to 3 decimals.
#
# What differs from running the script itself
#   - the script's header (sourcing the method, the loader, the packages) is
#     replaced by tr3_source_method() and tr3_load_real(): the loader's
#     setwd() names a folder that is gone, and readxl is not installed;
#   - KS-CCDs: the cluster has tr3_cores() workers instead of 20, and
#     delta_range is taken from the script's header line. Neither changes a
#     result.
#
# USAGE
#   Rscript revision_experiments/tr3/12_real_reference.R <method> [datasets]
#   (default datasets: iris,seeds,R15)
#
# OUTPUT  revision_experiments/results/tr3/real_reference.csv
#         method, dataset, ari, sil, k_hat, delta_best (KS-CCDs), script

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript revision_experiments/tr3/12_real_reference.R <method> [datasets]")
source(here::here("revision_experiments/tr3/tr3_lib.R"))

METHOD   <- args[1]
DATASETS <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else c("iris", "seeds", "R15")
stopifnot(METHOD %in% TR3_METHODS)
OUT_CSV <- tr3_path("real_reference.csv")

SCRIPT <- unique(vapply(DATASETS, function(ds) tr3_real_setting(METHOD, ds)$script, ""))
stopifnot(length(SCRIPT) == 1)
src <- readLines(here::here(SCRIPT), warn = FALSE)
hdr <- grep("^##### \\S+ data #####\\s*$", src)
blk <- sub("^##### (\\S+) data #####\\s*$", "\\1", src[hdr])
tail_line <- grep("^(t2 = Sys\\.time\\(\\)|stopCluster\\(cl\\))", src)
tail_line <- tail_line[tail_line > max(hdr)][1]
ends <- c(hdr[-1] - 1, if (is.na(tail_line)) length(src) else tail_line - 1)
missing <- setdiff(DATASETS, blk)
if (length(missing)) stop("12_real_reference.R: no block in ", SCRIPT, " for ", paste(missing, collapse = ", "))

tr3_source_method(METHOD, "realdata")
real <- tr3_load_real(DATASETS)
if (METHOD == "KS-CCDs") {
  eval(parse(text = grep("^delta_range\\s*=", src[seq_len(hdr[1] - 1)], value = TRUE)))
  cl <- makeCluster(tr3_cores())
  registerDoParallel(cl)
}

for (ds in DATASETS) {
  i <- match(ds, blk)
  assign(ds, real[[ds]], envir = globalenv())
  for (v in c("ari", "sil", "n_clusters", "delta_best")) if (exists(v, inherits = FALSE)) rm(list = v)
  cat(sprintf("12_real_reference.R: %s, block lines %d-%d of %s\n", ds, hdr[i], ends[i], SCRIPT))
  eval(parse(text = src[hdr[i]:ends[i]]), envir = globalenv())
  tr3_append(OUT_CSV, data.frame(method = METHOD, dataset = ds, ari = fmt17(ari), sil = fmt17(sil),
                                 k_hat = n_clusters,
                                 delta_best = if (exists("delta_best", inherits = FALSE)) fmt17(delta_best) else "NA",
                                 script = SCRIPT, stringsAsFactors = FALSE))
}
if (METHOD == "KS-CCDs") stopCluster(cl)
