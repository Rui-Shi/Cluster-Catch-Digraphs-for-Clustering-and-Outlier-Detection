#!/usr/bin/env Rscript
# revision_experiments/tr3/38_descend_pointer_check.R
#
# The descending radius search of UN-CCDs (d >= 10) in R/ccds/UN_CCD.R uses
# sorted neighbour lists since 2026-11-02 (nnccd.sorted.nbrs(), pointer per
# ball member), O(n^3) per data set instead of O(n^4) in the worst case. This
# script checks that nnccd.radi(method = "descend") returns EXACTLY the radii
# of the 2026-10-23 version, which is copied below verbatim (renamed with the
# suffix .20261023), and reports the time of both.
#
# Cases (low.num = 3, scores = FALSE unless stated; production tables, never generated)
#   sim   the first 2 data sets of each of 2, 3, 5 clusters of every UN-CCD
#         simulation cell with d = 10 or 20 (uniform, gaussian; n = 50, 100,
#         200, 500), drawn by tr3_sim_data() as 10_sim.R draws them
#   real  the six data sets of 31_hd_real.R with their tables
#   edge  d = 10 table: duplicated rows, heavy triplicates, lattice data with
#         many tied distances, n = 4 to 8, low.num 2 and 4, scores = TRUE
# A case passes if the radii are identical(), or if both versions stop with
# the same error message (degenerate inputs such as heavy duplicates).
#
# USAGE   Rscript revision_experiments/tr3/38_descend_pointer_check.R [groups]
#           [groups] comma-separated subset of sim,real,edge (default all)
# OUTPUT  revision_experiments/results/tr3/descend_pointer_check.csv
#         (rewritten on every run); stops with an error if any case fails

source(here::here("revision_experiments/tr3/tr3_lib.R"))
source(here::here("revision_experiments/tr3/30_hd_load.R"))
args <- commandArgs(trailingOnly = TRUE)
GROUPS <- if (length(args) >= 1 && nzchar(args[1])) strsplit(args[1], ",")[[1]] else c("sim", "real", "edge")
stopifnot(all(GROUPS %in% c("sim", "real", "edge")))
OUT_CSV <- tr3_path("descend_pointer_check.csv")

tr3_source_method("UN-CCDs", "simulation")
stopifnot(exists("nnccd.sorted.nbrs"), "S" %in% names(formals(nnccd.radi.descend.one)))

# ---- reference: R/ccds/UN_CCD.R at commit fddba5d (search code of ff343d6, 2026-10-23), verbatim ----
nnccd.radi.descend.one.20261023 <- function(i, ddx, o.d, n, low.num, env.ave.rev, env.med.rev,
                                   nn1.val, nn1.idx, nn2.val) {
  mem <- o.d[1:(n - 1)]                       # members at j = 1
  e <- o.d[n]                                 # the one point outside the ball
  nn <- ifelse(nn1.idx[mem] == e, nn2.val[mem], nn1.val[mem])   # position k <-> point mem[k]
  for (j in 1:(n - low.num)) {
    if (j > 1) {                              # point o.d[j-1] leaves; members are mem[j:(n-1)]
      x <- o.d[j - 1]
      pos <- j:(n - 1)
      hit <- pos[nn[pos] == ddx[mem[pos], x]]
      for (k in hit) {
        others <- mem[pos[pos != k]]
        nn[k] <- if (length(others)) min(ddx[mem[k], others]) else Inf
      }
    }
    r <- ddx[i, o.d[j]]
    cur <- nn[j:(n - 1)]
    if (mean(cur) / r > env.ave.rev[j + 1] & median(cur) / r > env.med.rev[j + 1]) {
      return(r)
    }
  }
  0
}

nnccd.radi.20261023 <- function(dx, quantile="lower", method="ascend", low.num, quant, simul=NULL, niter, scores=F){

  ddx <- as.matrix(dist(dx)) # the distance matrix
  n <- nrow(dx)
  d <- ncol(dx)
  R <- rep(0,n)

  if(quantile=="lower"){
    # a table not on disk (placeholder) or generated too short is generated at n (quantile_table.R)
    if(!is.null(simul) && !is.null(attr(simul, "ccd_table")) && .ccd_qt_extent(simul, "NN") < n) simul <- .ccd_qt_grow(simul, "NN", d, n)
    if(!is.null(simul)) {NN.envelop <- list(average=simul$average[1:n],median=simul$median[1:n])}
    else {NN.envelop <- NNDest.simpois.lower.quant(n, d, quant, niter)}
    if(method=="descend"){
      env.ave.rev <- rev(NN.envelop$average); env.med.rev <- rev(NN.envelop$median)
      dd0 <- ddx; diag(dd0) <- Inf
      nn1.idx <- apply(dd0, 1, which.min)
      nn1.val <- dd0[cbind(1:n, nn1.idx)]
      dd0[cbind(1:n, nn1.idx)] <- Inf
      nn2.val <- apply(dd0, 1, min)
      rm(dd0)
    }
    for(i in 1:n){
      if(method == "ascend"){
        o.d <- order(ddx[i,]) # the ascending distance order for i_th object
        R[i] <- nnccd.radi.ascend.one(i, ddx, o.d, n, low.num, NN.envelop$average, NN.envelop$median)
      }
      if(method=="descend"){
        o.d <- order(ddx[i,], decreasing=T) # the descending distance order for i_th object
        R[i] <- nnccd.radi.descend.one.20261023(i, ddx, o.d, n, low.num, env.ave.rev, env.med.rev,
                                       nn1.val, nn1.idx, nn2.val)
      }
      if(scores && R[i]==0){R[i]=sort(ddx[i,])[2]} # avoid 0 radius (necessary for outlyingness scores!)
    }
  }
  return(list(R=R,KS=NULL))
}
# ---- end of reference ----

rows <- list()
run_case <- function(group, case, X, simul, low.num = 3, scores = FALSE) {
  one <- function(f) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(f(X, quantile = "lower", method = "descend", low.num = low.num, quant = NULL,
                    simul = simul, niter = 1000, scores = scores)$R,
                  error = function(e) structure(conditionMessage(e), class = "err"))
    list(r = r, t = proc.time()[["elapsed"]] - t0)
  }
  a <- one(nnccd.radi.20261023)
  b <- one(nnccd.radi)
  ea <- inherits(a$r, "err"); eb <- inherits(b$r, "err")
  both_error <- ea && eb
  same <- if (both_error) identical(unclass(a$r), unclass(b$r)) else if (ea || eb) FALSE else identical(a$r, b$r)
  row <- data.frame(group = group, case = case, n = nrow(X), d = ncol(X), low_num = low.num, scores = scores,
                    pass = same, both_error = both_error,
                    error = if (ea || eb) paste(if (ea) unclass(a$r) else "-", "|", if (eb) unclass(b$r) else "-") else "",
                    zero_radii = if (!eb) sum(b$r == 0) else NA,
                    sec_20261023 = fmt17(a$t), sec_new = fmt17(b$t), stringsAsFactors = FALSE)
  rows[[length(rows) + 1]] <<- row
  cat(sprintf("%-5s %-28s n=%4d d=%2d low=%d scores=%-5s pass=%-5s%s  %.2fs -> %.2fs\n", group, case, nrow(X),
              ncol(X), low.num, scores, same, if (both_error) " (both error)" else "", a$t, b$t))
  flush.console()
}

if ("sim" %in% GROUPS) {
  for (dist in c("uniform", "gaussian")) for (d in c(10, 20)) for (n in c(50, 100, 200, 500)) {
    s <- tr3_sim_setting("UN-CCDs", dist, d, n)
    stopifnot(s$direction == "descend")
    simul <- tr3_load_table(s$table)
    data <- tr3_sim_data(s)
    for (N in c(2, 3, 5)) for (r in 1:2) {
      run_case("sim", sprintf("%s_d%d_n%d_k%d_rep%d", dist, d, n, N, r), data[[as.character(N)]]$data[[r]], simul)
    }
  }
}

if ("real" %in% GROUPS) {
  hd <- utils::read.csv(file.path(TR3_DIR, "settings_realdata_hd.csv"), colClasses = "character")
  real <- tr3_load_hd(TR3_HD_DATASETS)
  for (ds in TR3_HD_DATASETS) {
    s <- hd[hd$method == "UN-CCDs" & hd$dataset == ds, ]
    stopifnot(nrow(s) == 1, s$direction == "descend")
    D <- real[[ds]]
    run_case("real", ds, D[, -ncol(D)], tr3_load_table(s$table))
  }
}

if ("edge" %in% GROUPS) {
  s <- tr3_sim_setting("UN-CCDs", "uniform", 10, 200)
  simul <- tr3_load_table(s$table)
  X <- tr3_sim_data(s)[["3"]]$data[[1]]
  run_case("edge", "dup10_d10", rbind(X, X[1:10, ]), simul)
  run_case("edge", "triplicates_d10", rbind(X[1:50, ], X[1:50, ], X[1:50, ]), simul)
  set.seed(20261102)
  run_case("edge", "lattice012_d10", matrix(sample(0:2, 150 * 10, replace = TRUE), 150), simul)
  run_case("edge", "lattice01_d10", matrix(sample(0:1, 120 * 10, replace = TRUE), 120), simul)
  for (m in 4:8) run_case("edge", sprintf("small_n%d_d10", m), matrix(runif(m * 10), m), simul)
  for (low in c(2, 4)) run_case("edge", "sim_lownum", X, simul, low.num = low)
  run_case("edge", "sim_scores", X, simul, scores = TRUE)
  run_case("edge", "dup10_scores", rbind(X, X[1:10, ]), simul, scores = TRUE)
}

out <- do.call(rbind, rows)
utils::write.csv(out, OUT_CSV, row.names = FALSE)
cat(sprintf("\n%d cases: %d pass (%d of them both error with the same message), %d fail; written %s\n",
            nrow(out), sum(out$pass), sum(out$pass & out$both_error), sum(!out$pass), OUT_CSV))
for (g in unique(out$group)) {
  o <- out[out$group == g & !out$both_error, ]
  cat(sprintf("  %-5s total seconds: 2026-10-23 %.1f, new %.1f\n", g,
              sum(as.numeric(o$sec_20261023)), sum(as.numeric(o$sec_new))))
}
if (any(!out$pass)) stop("38_descend_pointer_check.R: radii differ in ", sum(!out$pass), " case(s)")
