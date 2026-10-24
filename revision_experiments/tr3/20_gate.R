#!/usr/bin/env Rscript
# revision_experiments/tr3/20_gate.R
#
# Reproduction gate for the TR3 drivers, in the style of tr1's WP0 gate. No
# new experiment runs through 10_sim.R or 11_real.R until every check here
# passes. The original scripts reproduce these numbers on this machine
# (CLUSTERING_REPRODUCTION_GUIDE.md section 9), so a mismatch is a driver bug
# until shown otherwise; nothing is tuned to force a match.
#
# Checks
#   Simulation cells, uniform clusters, d = 2, n = 50 (10_sim.R), against the
#   original SLURM logs, to all 15 printed significant digits; for each of 2,
#   3 and 5 clusters: mean ARI, mean Sil, SR, and for KS-CCDs the chosen m:
#     UN-CCDs  simulations/clustering/UN_CCDs/Simulations/Uniform/2d/slurm-1442992.out
#     RK-CCDs  simulations/clustering/RK_CCDs/Simulations/Uniform/2d/slurm-1432233.out
#     KS-CCDs  simulations/clustering/KS_CCDs/Simulations/Uniform/2d/slurm-1450256.out
#              (the 200-value grid; slurm-1442358.out is from an older grid)
#   Real data, iris, seeds and R15 (11_real.R), ARI, Sil and k_hat:
#     - against the manuscript's table Real_Data_Result_Clustering, to the 3
#       decimals printed there (values copied below);
#     - UN-CCDs and RK-CCDs against results/tr3/real_reference.csv (the
#       original script blocks, run by 12_real_reference.R), to 15 digits;
#     - KS-CCDs against the knit output Real_Data_KSCCD.html, to 15 digits,
#       and the chosen m.
#
# USAGE   Rscript revision_experiments/tr3/20_gate.R
# OUTPUT  revision_experiments/results/tr3/gate.csv (rewritten on every run):
#         one row per comparison, with pass = TRUE, FALSE or NA (not run yet)

source(here::here("revision_experiments/tr3/tr3_lib.R"))
OUT_CSV <- tr3_path("gate.csv")

GATE_SIM <- data.frame(
  method = c("UN-CCDs", "RK-CCDs", "KS-CCDs"), dist = "uniform", d = "2", n = "50",
  log = c("simulations/clustering/UN_CCDs/Simulations/Uniform/2d/slurm-1442992.out",
          "simulations/clustering/RK_CCDs/Simulations/Uniform/2d/slurm-1432233.out",
          "simulations/clustering/KS_CCDs/Simulations/Uniform/2d/slurm-1450256.out"),
  stringsAsFactors = FALSE)

# UNCCD4Clustering.tex, table Real_Data_Result_Clustering, columns iris, seeds, R15
PAPER <- utils::read.csv(colClasses = "character", text = "
method,dataset,ari,sil,k_hat
KS-CCDs,iris,0.561,0.578,2
KS-CCDs,seeds,0.497,0.468,2
KS-CCDs,R15,0.989,0.753,15
RK-CCDs,iris,0.557,0.393,3
RK-CCDs,seeds,0.500,0.164,3
RK-CCDs,R15,0.993,0.753,15
UN-CCDs,iris,0.599,0.445,3
UN-CCDs,seeds,0.239,0.286,2
UN-CCDs,R15,0.975,0.742,15")
KS_HTML <- "simulations/clustering/KS_CCDs/Simulations/Real_Data/Real_Data_KSCCD.html"

out <- list()
add <- function(check, method, cell, metric, source, expected, got, pass) {
  out[[length(out) + 1]] <<- data.frame(check = check, method = method, cell = cell, metric = metric,
                                        source = source, expected = as.character(expected),
                                        got = as.character(got), pass = pass, stringsAsFactors = FALSE)
}
# equal to the 15 significant digits that paste()/print() show
same15 <- function(got, expected) {
  if (is.na(got)) return(NA)
  # both sides rounded: logs carry 15 digits, real_reference.csv carries 17
  isTRUE(as.numeric(sprintf("%.15g", got)) == as.numeric(sprintf("%.15g", as.numeric(expected))))
}
show15 <- function(x) if (is.na(x)) "NA" else sprintf("%.15g", x)

# ---------------------------------------------------------------------------
# simulation cells
# ---------------------------------------------------------------------------
parse_log <- function(path) {
  txt <- readLines(here::here(path), warn = FALSE)
  m <- regmatches(txt, regexec(paste0(
    "The mean ARI is ([^ ]+) , The mean average silhouette index is ([^ ]+) , ",
    "and, the percentatage of times identifying the correct number of clusters is ([^ ]+) %",
    "(?:, best delta: ([^ \"]+))?"), txt, perl = TRUE))
  m <- Filter(length, m)
  if (length(m) != 3) stop("20_gate.R: expected 3 result lines in ", path, ", found ", length(m))
  data.frame(N = c(2, 3, 5), ari = vapply(m, `[`, "", 2), sil = vapply(m, `[`, "", 3),
             sr = vapply(m, `[`, "", 4), delta = vapply(m, `[`, "", 5), stringsAsFactors = FALSE)
}

keys <- c("variant", "method", "dist", "d", "n", "n_clusters")
reps <- tr3_read_done(tr3_path("sim_replicates.csv"), tr3_path("sim_done.csv"), keys, "rep")
for (g in seq_len(nrow(GATE_SIM))) {
  G <- GATE_SIM[g, ]
  lg <- parse_log(G$log)
  for (j in seq_len(nrow(lg))) {
    N <- lg$N[j]
    cell <- sprintf("%s d=%s n=%s, %d clusters", G$dist, G$d, G$n, N)
    r <- if (is.null(reps)) NULL else
      reps[reps$variant == "published" & reps$method == G$method & reps$dist == G$dist &
             reps$d == G$d & reps$n == G$n & reps$n_clusters == as.character(N), ]
    ran <- !is.null(r) && nrow(r) > 0
    if (ran) {
      r <- r[order(as.integer(r$rep)), ]
      if (!identical(as.integer(r$rep), seq_len(as.integer(r$iteN[1])))) stop("20_gate.R: incomplete cell ", cell)
      got <- c(ari = mean(as.numeric(r$ari)), sil = mean(as.numeric(r$sil)),
               sr = mean(as.numeric(r$k_hat) == N) * 100, delta = as.numeric(r$param[1]))
    } else {
      got <- c(ari = NA, sil = NA, sr = NA, delta = NA)
    }
    for (mt in c("ari", "sil", "sr")) {
      add("simulation", G$method, cell, mt, basename(G$log), lg[[mt]][j], show15(got[[mt]]),
          same15(got[[mt]], lg[[mt]][j]))
    }
    if (nzchar(lg$delta[j])) {
      add("simulation", G$method, cell, "delta", basename(G$log), lg$delta[j], show15(got[["delta"]]),
          same15(got[["delta"]], lg$delta[j]))
    }
  }
}

# ---------------------------------------------------------------------------
# real data
# ---------------------------------------------------------------------------
parse_ks_html <- function(path) {
  txt <- readLines(here::here(path), warn = FALSE)
  txt <- unlist(strsplit(txt, "\n", fixed = TRUE))
  # printed output only, not the echoed code
  hits <- grep("## \\[1\\] &quot;\\S+ x KS-CCDs: the ARI is", txt)
  rows <- lapply(hits, function(h) {
    m <- regmatches(txt[h], regexec(paste0("## \\[1\\] &quot;(\\S+) x KS-CCDs: the ARI is (\\S+) the average ",
                                           "silhouette index is (\\S+) number of clusters detected (\\d+)&quot;"),
                                    txt[h]))[[1]]
    # the chosen m is printed on the last "## [1] <number>" line before the result
    prev <- regmatches(txt[seq_len(h - 1)],
                       regexpr("## \\[1\\] [-0-9.e]+(?=\\s*($|<))", txt[seq_len(h - 1)], perl = TRUE))
    data.frame(dataset = m[2], ari = m[3], sil = m[4], k_hat = m[5],
               delta = sub("^## \\[1\\] ", "", prev[length(prev)]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

real <- tr3_read_done(tr3_path("real.csv"), tr3_path("real_done.csv"), c("variant", "method", "dataset"))
ref <- tr3_read(tr3_path("real_reference.csv"))
if (!is.null(ref)) ref <- ref[!duplicated(ref[c("method", "dataset")], fromLast = TRUE), ]
ks_html <- parse_ks_html(KS_HTML)

for (i in seq_len(nrow(PAPER))) {
  P <- PAPER[i, ]
  r <- if (is.null(real)) NULL else
    real[real$variant == "published" & real$method == P$method & real$dataset == P$dataset, ]
  ran <- !is.null(r) && nrow(r) == 1
  got <- if (ran) c(ari = as.numeric(r$ari), sil = as.numeric(r$sil), k_hat = as.numeric(r$k_hat),
                    delta = as.numeric(r$param)) else c(ari = NA, sil = NA, k_hat = NA, delta = NA)

  for (mt in c("ari", "sil")) {
    g3 <- if (is.na(got[[mt]])) "NA" else sprintf("%.3f", got[[mt]])
    add("real", P$method, P$dataset, mt, "manuscript, 3 dp", P[[mt]], g3,
        if (is.na(got[[mt]])) NA else g3 == P[[mt]])
  }
  add("real", P$method, P$dataset, "k_hat", "manuscript", P$k_hat, got[["k_hat"]],
      if (is.na(got[["k_hat"]])) NA else got[["k_hat"]] == as.numeric(P$k_hat))

  if (P$method == "KS-CCDs") {
    h <- ks_html[ks_html$dataset == P$dataset, ]
    if (nrow(h) != 1) stop("20_gate.R: no KS-CCD result for ", P$dataset, " in ", KS_HTML)
    for (mt in c("ari", "sil", "delta")) {
      add("real", P$method, P$dataset, mt, basename(KS_HTML), h[[mt]], show15(got[[mt]]),
          same15(got[[mt]], h[[mt]]))
    }
    add("real", P$method, P$dataset, "k_hat", basename(KS_HTML), h$k_hat, got[["k_hat"]],
        if (is.na(got[["k_hat"]])) NA else got[["k_hat"]] == as.numeric(h$k_hat))
  } else {
    f <- if (is.null(ref)) NULL else ref[ref$method == P$method & ref$dataset == P$dataset, ]
    have_ref <- !is.null(f) && nrow(f) == 1
    for (mt in c("ari", "sil", "k_hat")) {
      e <- if (have_ref) f[[mt]] else "not run"
      add("real", P$method, P$dataset, mt, "real_reference.csv (original script block)",
          e, show15(got[[mt]]), if (have_ref) same15(got[[mt]], e) else NA)
    }
  }
}

res <- do.call(rbind, out)
dir.create(TR3_RES, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, OUT_CSV, row.names = FALSE)

np <- sum(res$pass %in% TRUE); nf <- sum(res$pass %in% FALSE); nn <- sum(is.na(res$pass))
by_cell <- aggregate(pass ~ check + method + cell, transform(res, pass = ifelse(is.na(pass), "not run",
                     ifelse(pass, "pass", "FAIL"))), function(p) paste(sort(unique(p)), collapse = "/"))
print(by_cell[order(by_cell$check, by_cell$method), ], row.names = FALSE)
if (nf > 0) print(res[res$pass %in% FALSE, c("method", "cell", "metric", "source", "expected", "got")], row.names = FALSE)
cat(sprintf("\n20_gate.R: %d comparisons: %d pass, %d FAIL, %d not run -> %s\n", nrow(res), np, nf, nn,
            if (nf > 0) "GATE FAILS" else if (nn > 0) "GATE INCOMPLETE" else "GATE PASSES"))
cat("written:", OUT_CSV, "\n")
