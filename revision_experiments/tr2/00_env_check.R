# 00_env_check.R
# Environment check for the PR-D-26-05767 revision-experiments pipeline (Task T0).
# (a) loads all required R packages
# (b) loads the two sample quantile lookup tables and asserts key components are
#     non-empty and NA-free
# (c) prints sessionInfo()
# (d) ends with "ENV CHECK OK"

pkgs <- c(
  "here", "MASS", "igraph", "foreach", "doParallel", "cluster",
  "mclust", "dbscan", "isotree", "solitude", "mvtnorm", "spatstat",
  "farff", "foreign", "R.matlab", "dplyr", "ggplot2"
)

cat("=== (a) Loading packages ===\n")
for (p in pkgs) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat("Loaded:", p, "\n")
}

cat("\n=== (b) Loading sample quantile lookup tables ===\n")

# Resolve paths relative to this script's directory, robust to how the script
# is invoked (Rscript path/to/00_env_check.R from any working directory).
this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg) == 1) normalizePath(file_arg) else NA_character_
}, error = function(e) NA_character_)

if (!is.na(this_file)) {
  script_dir <- dirname(this_file)
} else {
  script_dir <- getwd()
}
repo_root <- normalizePath(file.path(script_dir, ".."))

rk_path <- file.path(repo_root, "R", "RK-test_quantile", "RK-test-simul_10d_999%.RData")
nn_path <- file.path(repo_root, "R", "NN-test_quantile", "NN-test-simul_10d_99%.RData")

stopifnot("RK sample quantile table not found" = file.exists(rk_path))
stopifnot("NN sample quantile table not found" = file.exists(nn_path))

e_rk <- new.env()
load(rk_path, envir = e_rk)
stopifnot("simul object missing from RK table" = exists("simul", envir = e_rk))
simul_rk <- get("simul", envir = e_rk)

stopifnot(
  "RK simul$Kest.m is empty" = length(simul_rk$Kest.m) > 0,
  "RK simul$Kest.m has NAs" = !any(is.na(simul_rk$Kest.m)),
  "RK simul$quan is empty" = length(simul_rk$quan) > 0,
  "RK simul$quan[[1]] has NAs" = !any(is.na(simul_rk$quan[[1]])),
  "RK simul$r is empty" = length(simul_rk$r) > 0,
  "RK simul$r has NAs" = !any(is.na(simul_rk$r))
)
cat("RK-test-simul_10d_999%.RData: OK (Kest.m ",
    paste(dim(simul_rk$Kest.m), collapse = "x"),
    ", quan[[1]] ", paste(dim(simul_rk$quan[[1]]), collapse = "x"),
    ", r length ", length(simul_rk$r), ")\n", sep = "")

e_nn <- new.env()
load(nn_path, envir = e_nn)
stopifnot("simul object missing from NN table" = exists("simul", envir = e_nn))
simul_nn <- get("simul", envir = e_nn)

stopifnot(
  "NN simul$average is empty" = length(simul_nn$average) > 0,
  "NN simul$average has NAs" = !any(is.na(simul_nn$average)),
  "NN simul$median is empty" = length(simul_nn$median) > 0,
  "NN simul$median has NAs" = !any(is.na(simul_nn$median))
)
cat("NN-test-simul_10d_99%.RData: OK (average length ",
    length(simul_nn$average), ", median length ", length(simul_nn$median),
    ")\n", sep = "")

cat("\n=== (c) sessionInfo() ===\n")
print(sessionInfo())

cat("\nENV CHECK OK\n")
