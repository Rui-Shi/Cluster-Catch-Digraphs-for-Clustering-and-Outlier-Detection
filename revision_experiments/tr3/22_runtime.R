#!/usr/bin/env Rscript
# revision_experiments/tr3/22_runtime.R
#
# Re-times one cell of tab:runtime_comparison of the UN-CCD clustering paper:
# the total run time in minutes of one simulation driver, uniform clusters,
# d = 5, 500 repetitions, number of clusters 2, 3 and 5, at one n.
#
# WHAT THE PUBLISHED NUMBERS ARE (no timing script existed)
#   They are the "Time difference" lines of the SLURM logs of the d = 5 uniform
#   drivers simulations/clustering/<M>_CCDs/Simulations/Uniform/5d/
#   5d_235clusters_n<n>.R (UN 1443072-75, RK 1442240-42 and 1443104, KS
#   1450264-67; all in the repository since 2025-09-12, i.e. before the
#   incremental radius search of 059c798 on 2026-09-23). In each driver
#     t1 = Sys.time()     comes after the sourcing of the method, the loading
#                         of the quantile table (offline) and library()
#     ... for each of 2, 3 and 5 clusters: generate iteN data sets,
#         makeCluster(detectCores()), foreach %dopar% over the data sets
#         (KS-CCDs: first one pass per value of the 200-value grid, then a
#         final pass at the best value), stopCluster, ARI, one printed line
#     t2 = Sys.time(); t2 - t1
#   So the time covers data generation, three PSOCK cluster start-ups and the
#   clustering of 3 x 500 data sets, and excludes table loading.
#   (The paper says "number of clusters varying from 2 to 5"; the drivers use
#   2, 3 and 5.)
#
# WHAT THIS SCRIPT DOES
#   The same, once, in a fresh R process, with the settings of the driver's row
#   in settings_simulation.csv (as 10_sim.R) and the methods called as the
#   driver called them (UN-CCDs: UNCCD_clustering(), the driver's
#   UNMCCD_outlier() being undefined; RK-CCDs: RKCCD_clustering(); KS-CCDs:
#   KSCCD_clustering()). It splits the total into data generation and the
#   three clustering blocks, and appends one row to
#   revision_experiments/results/tr3/runtime_reps.csv.
#   Workers: TR3_CORES, default detectCores() (as the drivers). The RK-CCD
#   table is loaded as the driver loaded it, with its raw draws (Kest.m), so
#   the export of the table to the workers is inside the timed region, as in
#   the published runs; TR3_RUNTIME_DROP_KEST=1 drops them (as 10_sim.R does).
#
# LOAD GUARD
#   Before the clock starts, the script starts 22_runtime_monitor.ps1, which
#   samples every 3 s the R processes that are not this driver's and the CPU
#   load. If one is running at the start, the script stops (override for tests:
#   TR3_RUNTIME_ALLOW_LOAD=1). If one shows up during the run, the row is
#   written with valid = FALSE and is left out of the medians.
#
# USAGE (from the clone root, in PowerShell, one timing process at a time)
#   $env:CCD_QUANTILE_GENERATE = "FALSE"
#   & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" revision_experiments/tr3/22_runtime.R <method> <n> <rep> [iteN]
#     <method>  UN-CCDs | RK-CCDs | KS-CCDs
#     <n>       50 | 100 | 200 | 500
#     <rep>     repetition number (the cell is skipped if (method, n, rep, iteN) is already in the CSV)
#     [iteN]    default 500 (the published design)
#   TR3_RUNTIME_OUT      result file name in results/tr3 (default runtime_reps.csv)
#   TR3_RUNTIME_FORCE=1  run even if the row exists
#   22_runtime_run.ps1 runs all cells and repetitions in order.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript revision_experiments/tr3/22_runtime.R <method> <n> <rep> [iteN]")
source(here::here("revision_experiments/tr3/tr3_lib.R"))

METHOD <- args[1]
NN     <- args[2]
REP    <- as.integer(args[3])
ITEN   <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "500"
stopifnot(METHOD %in% TR3_METHODS, NN %in% c("50", "100", "200", "500"), !is.na(REP))
OUT_CSV  <- tr3_path(Sys.getenv("TR3_RUNTIME_OUT", "runtime_reps.csv"))
KEY      <- list(method = METHOD, n = NN, rep = as.character(REP), iteN = ITEN)
if (Sys.getenv("TR3_RUNTIME_FORCE", "") != "1" && has_result(OUT_CSV, KEY)) {
  cat("22_runtime.R: already done:", paste(names(KEY), unlist(KEY), sep = "=", collapse = " "), "\n")
  quit(save = "no", status = 0)
}

s <- tr3_sim_setting(METHOD, "uniform", "5", NN)
s$iteN <- ITEN
cores <- tr3_cores()

# ---- before the clock: method, packages, quantile table (offline) -----------
tr3_source_method(METHOD, "simulation")
load_table <- function(file) {
  variant <- if (startsWith(file, "RK-")) "RK" else "NN"
  e <- new.env()
  load_quantile_table(here::here(sprintf("R/%s-test_quantile/%s", variant, file)), envir = e)
  tab <- get("simul", envir = e)
  if (variant == "RK" && Sys.getenv("TR3_RUNTIME_DROP_KEST", "") == "1") tab$Kest.m <- NULL
  tab
}
simul <- if (METHOD == "KS-CCDs") NULL else load_table(s$table)
direction <- s$direction                                  # UN-CCDs
quant <- if (nzchar(s$quant)) as.numeric(s$quant) else NA # RK-CCDs
delta_range <- if (METHOD == "KS-CCDs") eval(parse(text = s$delta_grid)) else NULL

# ---- machine record and load guard -------------------------------------------
ps_out <- function(cmd) {
  tryCatch(paste(system2("powershell", c("-NoProfile", "-Command", shQuote(cmd)), stdout = TRUE, stderr = FALSE),
                 collapse = " "), error = function(e) "NA")
}
cpu_name <- trimws(ps_out("(Get-CimInstance Win32_Processor | Select-Object -First 1).Name"))
ram_gb   <- trimws(ps_out("[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)"))
plan     <- trimws(ps_out("(powercfg /getactivescheme) -replace '.*\\((.*)\\).*','$1'"))
os_name  <- trimws(ps_out("(Get-CimInstance Win32_OperatingSystem).Caption + ' ' + (Get-CimInstance Win32_OperatingSystem).Version"))

mon_file <- tempfile("runtime_monitor_", fileext = ".csv")
system2("powershell", c("-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
                        shQuote(file.path(TR3_DIR, "22_runtime_monitor.ps1")),
                        "-DriverPid", Sys.getpid(), "-Out", shQuote(mon_file), "-Every", 3),
        wait = FALSE, stdout = FALSE, stderr = FALSE)
for (i in 1:60) { if (file.exists(mon_file) && length(readLines(mon_file, warn = FALSE)) >= 2) break; Sys.sleep(0.5) }
mon_read <- function() {
  m <- tryCatch(utils::read.csv(mon_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(m) || !nrow(m)) data.frame(time = NA, foreign_r = NA, foreign_ids = "", sys_cpu = NA) else m
}
m0 <- mon_read()
allow_load <- Sys.getenv("TR3_RUNTIME_ALLOW_LOAD", "") == "1"
if (is.na(m0$foreign_r[1])) stop("22_runtime.R: the load monitor did not start")
if (m0$foreign_r[1] > 0 && !allow_load) {
  stop("22_runtime.R: other R processes are running (", m0$foreign_ids[1], "); not timing under load")
}
cat(sprintf("22_runtime.R: %s uniform d=5 n=%s iteN=%s rep=%d cores=%d | %s | %s\n",
            METHOD, NN, ITEN, REP, cores, R.version.string, cpu_name))

# ---- the timed region (the driver's t1 ... t2) --------------------------------
mean_ari <- mean_sil <- sr <- t_blk <- c("2" = NA_real_, "3" = NA_real_, "5" = NA_real_)
delta_best <- c("2" = NA_real_, "3" = NA_real_, "5" = NA_real_)
t1 <- Sys.time()
data <- tr3_sim_data(s)                      # set.seed(2024), then the 2-, 3-, 5-cluster lists
t_data <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

# Every foreach() below is called at top level, as in the drivers.
for (N in c(2, 3, 5)) {
  tb <- Sys.time()
  data_list <- data[[as.character(N)]]$data
  label_T <- data[[as.character(N)]]$label_T
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  if (METHOD == "KS-CCDs") {
    Sil_seq <- c()
    for (delta in delta_range) {
      clustering_result <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% KSCCD_clustering(datax = x, m = delta)
      Sil_seq <- c(Sil_seq, mean(sapply(clustering_result, function(x) x$si)))
    }
    delta_best[as.character(N)] <- delta_range[which.max(Sil_seq)]
    delta <- delta_best[[as.character(N)]]
    clustering_result <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% KSCCD_clustering(datax = x, m = delta)
  } else if (METHOD == "RK-CCDs") {
    clustering_result <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% RKCCD_clustering(datax = x, simul = simul, quant = quant)
  } else {
    clustering_result <- foreach(x = data_list, .packages = c("MASS", "cluster")) %dopar% UNCCD_clustering(datax = x, simul = simul, method = direction)
  }
  stopCluster(cl)
  ari <- sapply(clustering_result, function(x) adjustedRandIndex(x$label, label_T))
  det <- sapply(clustering_result, function(x) x$si.ind)
  sil <- sapply(clustering_result, function(x) x$si)
  tr3_print_summary(METHOD, ari, sil, det, N, if (METHOD == "KS-CCDs") delta_best[[as.character(N)]] else NULL)
  mean_ari[as.character(N)] <- mean(ari); mean_sil[as.character(N)] <- mean(sil); sr[as.character(N)] <- mean(det == N)
  t_blk[as.character(N)] <- as.numeric(difftime(Sys.time(), tb, units = "secs"))
}
t2 <- Sys.time()
t_total <- as.numeric(difftime(t2, t1, units = "secs"))
print(t2 - t1)

# ---- after the clock ------------------------------------------------------------
Sys.sleep(1)
m <- mon_read()
m <- m[!is.na(m$foreign_r), , drop = FALSE]
foreign_max <- if (nrow(m)) max(m$foreign_r) else NA
valid <- !is.na(foreign_max) && foreign_max == 0 && nrow(m) >= 1
tr3_append(OUT_CSV, data.frame(
  method = METHOD, dist = "uniform", d = "5", n = NN, rep = REP, iteN = ITEN, cores = cores,
  seconds = fmt17(t_total), minutes = fmt17(t_total / 60),
  secs_data = fmt17(t_data), secs_2cls = fmt17(t_blk[["2"]]), secs_3cls = fmt17(t_blk[["3"]]), secs_5cls = fmt17(t_blk[["5"]]),
  ari_2 = fmt17(mean_ari[["2"]]), ari_3 = fmt17(mean_ari[["3"]]), ari_5 = fmt17(mean_ari[["5"]]),
  sil_2 = fmt17(mean_sil[["2"]]), sil_3 = fmt17(mean_sil[["3"]]), sil_5 = fmt17(mean_sil[["5"]]),
  sr_2 = fmt17(sr[["2"]]), sr_3 = fmt17(sr[["3"]]), sr_5 = fmt17(sr[["5"]]),
  best_delta = if (METHOD == "KS-CCDs") paste(fmt17(delta_best), collapse = ";") else "NA",
  valid = valid, n_samples = nrow(m), foreign_r_max = foreign_max,
  sys_cpu_mean = fmt17(mean(m$sys_cpu, na.rm = TRUE)), load_override = allow_load,
  rk_kest_dropped = (METHOD == "RK-CCDs" && Sys.getenv("TR3_RUNTIME_DROP_KEST", "") == "1"),
  r_version = R.version.string, cpu = cpu_name, ram_gb = ram_gb, power_plan = plan, os = os_name,
  started = format(t1, "%Y-%m-%d %H:%M:%S"), stringsAsFactors = FALSE))
cat(sprintf("22_runtime.R: %.2f minutes, valid=%s (monitor: %d samples, max foreign R %s)\n",
            t_total / 60, valid, nrow(m), foreign_max))
