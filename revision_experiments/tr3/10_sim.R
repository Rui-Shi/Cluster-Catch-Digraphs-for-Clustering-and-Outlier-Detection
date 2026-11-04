#!/usr/bin/env Rscript
# revision_experiments/tr3/10_sim.R
#
# One simulation cell group of the UN-CCD clustering paper (Tables S2-S4 and
# S6-S8), run the way the original driver
# simulations/clustering/<M>_CCDs/Simulations/<Dist>/<d>d/<d>d_235clusters_n<n>.R
# ran it, with every replicate written to a CSV.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/10_sim.R <method> <dist> <d> <n> [clusters] [variant]
#     <method>   UN-CCDs | RK-CCDs | KS-CCDs
#     <dist>     uniform | gaussian
#     [clusters] comma-separated subset of 2,3,5 (default all three)
#     [variant]  published (default), or one of VARIANTS below; each variant
#                has a protocol file in this folder, written before its runs
#   TR3_CORES overrides the number of workers (default detectCores()). The
#   results do not depend on it: data are generated on the main process and
#   the clustering is deterministic.
#
# What is reproduced from the script
#   - the settings of its row in settings_simulation.csv (00_check_settings.R
#     checks them against the script);
#   - set.seed(seed) once, then the 2-, 3- and 5-cluster lists, always all
#     three (tr3_sim_data());
#   - foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% over a
#     PSOCK cluster, with one method loaded, as the script did;
#   - per replicate: ARI = adjustedRandIndex(label, label_T), Sil = si,
#     k_hat = si.ind, and the summary line it printed;
#   - KS-CCDs: for every m in the 200-value grid, the mean si over the iteN
#     data sets, then m = grid[which.max(...)] and a final run at that m.
# The one change: UNMCCD_outlier(), which no sourced file defines, is
# UNCCD_clustering() with the same arguments (guide section 6.1).
#
# OUTPUT (revision_experiments/results/tr3/)
#   sim_replicates.csv  one row per replicate: variant, method, dist, d, n,
#                       n_clusters, iteN, rep, ari, sil, k_hat, seconds,
#                       param (KS-CCDs: the chosen m), data_fp
#   sim_done.csv        one row per finished cell; reruns skip these cells
#   sim_ks_search.csv   KS-CCDs: one row per grid value (delta_index, delta,
#                       mean_sil, mean_ari, sr, seconds); a rerun resumes
#                       the search at the first missing grid value

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript revision_experiments/tr3/10_sim.R <method> <dist> <d> <n> [clusters] [variant]")
}
source(here::here("revision_experiments/tr3/tr3_lib.R"))

METHOD   <- args[1]
DIST     <- args[2]
D        <- args[3]
NN       <- args[4]
CLUSTERS <- if (length(args) >= 5 && nzchar(args[5])) as.numeric(strsplit(args[5], ",")[[1]]) else c(2, 3, 5)
VARIANT  <- if (length(args) >= 6 && nzchar(args[6])) args[6] else "published"
stopifnot(METHOD %in% TR3_METHODS, all(CLUSTERS %in% c(2, 3, 5)))

# Each variant replaces named settings of the published cell, for the cells
# its protocol covers and nowhere else.
VARIANTS <- list(
  iteN500 = list(protocol = "ITEN500_PROTOCOL.md", set = list(iteN = "500"),
                 cells = c("RK-CCDs gaussian 3 500", "RK-CCDs uniform 3 500")))
if (VARIANT != "published") {
  v <- VARIANTS[[VARIANT]]
  if (is.null(v)) stop("10_sim.R: variant '", VARIANT, "' has no protocol yet")
  if (!file.exists(file.path(TR3_DIR, v$protocol))) stop("10_sim.R: protocol ", v$protocol, " is missing")
  if (!paste(METHOD, DIST, D, NN) %in% v$cells) {
    stop("10_sim.R: variant ", VARIANT, " covers only ", paste(v$cells, collapse = "; "))
  }
}

REP_CSV    <- tr3_path("sim_replicates.csv")
DONE_CSV   <- tr3_path("sim_done.csv")
SEARCH_CSV <- tr3_path("sim_ks_search.csv")

s <- tr3_sim_setting(METHOD, DIST, D, NN)
if (VARIANT != "published") for (k in names(v$set)) s[[k]] <- v$set[[k]]
tr3_source_method(METHOD, "simulation")

par <- switch(METHOD,
  "UN-CCDs" = list(simul = tr3_load_table(s$table), direction = s$direction),
  "RK-CCDs" = list(simul = tr3_load_table(s$table), quant = as.numeric(s$quant)),
  "KS-CCDs" = list(grid = eval(parse(text = s$delta_grid))))

t_start <- Sys.time()
data <- tr3_sim_data(s)
cores <- tr3_cores()
cat(sprintf("10_sim.R: %s %s d=%s n=%s iteN=%s clusters=%s variant=%s cores=%d table=%s\n",
            METHOD, DIST, D, NN, s$iteN, paste(CLUSTERS, collapse = ","), VARIANT, cores,
            if (nzchar(s$table)) s$table else "-"))

cl <- makeCluster(cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

cell_key <- function(N) list(variant = VARIANT, method = METHOD, dist = DIST, d = D, n = NN,
                             n_clusters = as.character(N))

# KS search rows of one cell, one per grid index (the last copy written).
ks_rows <- function(key) {
  h <- tr3_read(SEARCH_CSV)
  if (is.null(h)) return(h)
  h <- h[Reduce(`&`, Map(function(k, v) h[[k]] == v, names(key), key)), , drop = FALSE]
  h <- h[!duplicated(h$delta_index, fromLast = TRUE), , drop = FALSE]
  h[order(as.integer(h$delta_index)), , drop = FALSE]
}

# Every foreach() below is called at top level, as in the scripts: foreach
# exports to the workers only what the calling environment holds, following
# global functions into their bodies (tr3_one -> ccd_cluster -> the method).
# Called inside a function it would not export tr3_one at all.
for (N in CLUSTERS) {
  key <- cell_key(N)
  if (has_result(DONE_CSV, key)) {
    cat(sprintf("  %d clusters: done, skipped\n", N))
    next
  }
  g <- data[[as.character(N)]]
  data_list <- g$data
  label_T <- g$label_T
  t0 <- Sys.time()

  p <- par
  delta_best <- NULL
  if (METHOD == "KS-CCDs") {
    grid <- par$grid
    todo <- setdiff(seq_along(grid), as.integer(ks_rows(key)$delta_index))
    if (length(todo) < length(grid)) {
      cat(sprintf("  KS search: %d of %d grid values already done\n",
                  length(grid) - length(todo), length(grid)))
    }
    for (i in todo) {
      ti <- proc.time()[["elapsed"]]
      p <- list(m = grid[i])
      r <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% tr3_one(x, label_T, METHOD, p)
      r <- do.call(rbind, r)
      tr3_append(SEARCH_CSV, data.frame(key, delta_index = i, delta = fmt17(grid[i]),
                                        mean_sil = fmt17(mean(r[, "sil"])),
                                        mean_ari = fmt17(mean(r[, "ari"])),
                                        sr = fmt17(mean(r[, "k_hat"] == N)),
                                        seconds = fmt17(proc.time()[["elapsed"]] - ti),
                                        stringsAsFactors = FALSE))
    }
    h <- ks_rows(key)
    if (!identical(as.integer(h$delta_index), seq_along(grid))) {
      stop("10_sim.R: the KS search for this cell is incomplete after the run")
    }
    # the script: delta_range[which.max(Sil_seq)], Sil_seq in grid order
    delta_best <- grid[which.max(as.numeric(h$mean_sil))]
    p <- list(m = delta_best)
  }
  r <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% tr3_one(x, label_T, METHOD, p)
  r <- do.call(rbind, r)

  rows <- data.frame(key, iteN = s$iteN, rep = seq_len(nrow(r)),
                     ari = fmt17(r[, "ari"]), sil = fmt17(r[, "sil"]),
                     k_hat = fmt17(r[, "k_hat"]), seconds = fmt17(r[, "seconds"]),
                     param = if (is.null(delta_best)) "NA" else fmt17(delta_best),
                     data_fp = fmt17(vapply(data_list, tr3_fp, numeric(1))),
                     stringsAsFactors = FALSE)
  tr3_append(REP_CSV, rows)
  tr3_append(DONE_CSV, data.frame(key, iteN = s$iteN, n_rows = nrow(rows), cores = cores,
                                  wall_secs = fmt17(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
                                  finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                  stringsAsFactors = FALSE))
  tr3_print_summary(METHOD, r[, "ari"], r[, "sil"], r[, "k_hat"], N, delta_best)
}

print(Sys.time() - t_start)
