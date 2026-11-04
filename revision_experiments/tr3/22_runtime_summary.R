#!/usr/bin/env Rscript
# revision_experiments/tr3/22_runtime_summary.R
#
# Summarises results/tr3/runtime_reps.csv (22_runtime.R): per method and n, the
# valid repetitions, the median, minimum and maximum in minutes, the published
# value (the "Time difference" of the SLURM log) and the ratio, and whether the
# repetition's means of ARI, Sil and SR (and KS-CCDs' best delta) equal the
# log's printed values. Writes results/tr3/runtime_summary.csv.
#   Rscript revision_experiments/tr3/22_runtime_summary.R [csv name] [summary name]

root <- here::here()
csv  <- file.path(root, "revision_experiments/results/tr3", if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "runtime_reps.csv")
d <- read.csv(csv, stringsAsFactors = FALSE, colClasses = "character")
for (k in c("n", "rep", "iteN", "cores")) d[[k]] <- as.integer(d[[k]])
d$minutes <- as.numeric(d$minutes)
d$valid <- d$valid == "TRUE"
d <- d[d$iteN == 500, ]

logs <- list(
  "UN-CCDs" = c("50" = "UN_CCDs/Simulations/Uniform/5d/slurm-1443072.out", "100" = "UN_CCDs/Simulations/Uniform/5d/slurm-1443073.out",
                "200" = "UN_CCDs/Simulations/Uniform/5d/slurm-1443074.out", "500" = "UN_CCDs/Simulations/Uniform/5d/slurm-1443075.out"),
  "RK-CCDs" = c("50" = "RK_CCDs/Simulations/Uniform/5d/slurm-1442240.out", "100" = "RK_CCDs/Simulations/Uniform/5d/slurm-1442241.out",
                "200" = "RK_CCDs/Simulations/Uniform/5d/slurm-1442242.out", "500" = "RK_CCDs/Simulations/Uniform/5d/slurm-1443104.out"),
  "KS-CCDs" = c("50" = "KS_CCDs/Simulations/Uniform/5d/slurm-1450264.out", "100" = "KS_CCDs/Simulations/Uniform/5d/slurm-1450265.out",
                "200" = "KS_CCDs/Simulations/Uniform/5d/slurm-1450266.out", "500" = "KS_CCDs/Simulations/Uniform/5d/slurm-1450267.out"))

read_log <- function(f) {
  x <- readLines(file.path(root, "simulations/clustering", f), warn = FALSE)
  td <- sub("Time difference of ", "", x[grepl("^Time difference of", x)])
  val <- as.numeric(sub(" .*", "", td)); unit <- sub("^[0-9.e+-]+ ", "", td)
  mins <- val * c(secs = 1/60, mins = 1, hours = 60)[unit]
  r <- x[grepl("mean ARI is", x)]
  num <- function(p) as.numeric(sub(p, "\\1", r))
  list(minutes = mins,
       ari = num(".*mean ARI is ([0-9.e-]+) ,.*"), sil = num(".*silhouette index is ([0-9.e-]+) ,.*"),
       sr = num(".*clusters is ([0-9.e-]+) %.*")/100,
       delta = suppressWarnings(as.numeric(sub(".*best delta: ([0-9.e-]+).*", "\\1", r[grepl("best delta", r)]))))
}

rows <- list()
for (m in names(logs)) for (n in names(logs[[m]])) {
  L <- read_log(logs[[m]][[n]])
  x <- d[d$method == m & d$n == as.integer(n), ]
  v <- x[x$valid, ]
  same <- NA
  if (nrow(v)) {
    same <- all(vapply(seq_len(nrow(v)), function(i) {
      a <- as.numeric(unlist(v[i, c("ari_2", "ari_3", "ari_5")])); s <- as.numeric(unlist(v[i, c("sil_2", "sil_3", "sil_5")]))
      r <- as.numeric(unlist(v[i, c("sr_2", "sr_3", "sr_5")]))
      ok <- isTRUE(all.equal(a, L$ari, tolerance = 1e-12)) && isTRUE(all.equal(s, L$sil, tolerance = 1e-12)) &&
        isTRUE(all.equal(r, L$sr, tolerance = 1e-12))
      if (m == "KS-CCDs") ok <- ok && isTRUE(all.equal(as.numeric(strsplit(v$best_delta[i], ";")[[1]]), L$delta, tolerance = 1e-12))
      ok }, logical(1)))
  }
  rows[[length(rows) + 1]] <- data.frame(
    method = m, n = as.integer(n), reps_total = nrow(x), reps_valid = nrow(v),
    cores = paste(unique(v$cores), collapse = "/"),
    median_min = if (nrow(v)) median(v$minutes) else NA, min_min = if (nrow(v)) min(v$minutes) else NA,
    max_min = if (nrow(v)) max(v$minutes) else NA, published_min = L$minutes,
    ratio_new_over_published = if (nrow(v)) median(v$minutes) / L$minutes else NA,
    results_equal_log = same, stringsAsFactors = FALSE)
}
out <- do.call(rbind, rows)
out[sapply(out, is.numeric)] <- lapply(out[sapply(out, is.numeric)], function(z) round(z, 3))
print(out, row.names = FALSE)
write.csv(out, file.path(root, "revision_experiments/results/tr3", if (length(commandArgs(TRUE)) > 1) commandArgs(TRUE)[2] else "runtime_summary.csv"), row.names = FALSE)
cat("\nMachine:", paste(unique(d$cpu), collapse = "; "), "|", paste(unique(d$r_version), collapse = "; "), "|",
    paste(unique(d$ram_gb), collapse = "/"), "GB |", paste(unique(d$os), collapse = "; "), "| power plan:", paste(unique(d$power_plan), collapse = "; "), "\n")
