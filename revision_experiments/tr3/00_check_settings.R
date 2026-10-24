#!/usr/bin/env Rscript
# revision_experiments/tr3/00_check_settings.R
#
# Checks the declared settings tables against the scripts they describe, and
# stops at the first disagreement:
#   settings_simulation.csv  vs the 120 drivers
#                            simulations/clustering/{UN,RK,KS}_CCDs/Simulations/
#                            {Uniform,Gaussian}/<d>d/<d>d_235clusters_n<n>.R
#   settings_realdata.csv    vs the blocks of the three Real_Data_<M>CCD.R scripts
#
# Two kinds of check:
#   1. Values. Each setting is read from the script text: n, d, iteN,
#      cls_dis, r_min, r_max, sig, the seed, the quantile-table file, the RK
#      level (the scripts' "if(d<10) 0.99 else 0.999", evaluated at d), the
#      UN search direction, the KS grid expression, and the size of the first
#      of three clusters (round(n/N) - 1, or + 1 at n = 100).
#   2. Structure. With those values blanked and whitespace removed, the 20
#      drivers of each (method, distribution) must be the same text, and so
#      must the blocks of each real-data script. So the table is all that
#      differs between them, and 10_sim.R / 11_real.R reproduce every driver
#      once they reproduce one.
#
# USAGE  Rscript revision_experiments/tr3/00_check_settings.R

source(here::here("revision_experiments/tr3/tr3_lib.R"))

fail <- function(...) stop(sprintf(...), call. = FALSE)

# value of the first "<var> = <value>" line, comments dropped; NA if none
assigned <- function(txt, var) {
  m <- regmatches(txt, regexec(sprintf("^\\s*%s\\s*=\\s*(.*?)\\s*(#.*)?$", var), txt, perl = TRUE))
  m <- Filter(length, m)
  if (length(m)) m[[1]][2] else NA_character_
}
count <- function(txt, pattern) sum(lengths(regmatches(txt, gregexpr(pattern, txt, fixed = TRUE))))
squash <- function(txt) gsub("\\s+", "", txt)

# ---------------------------------------------------------------------------
# simulation drivers
# ---------------------------------------------------------------------------
s <- tr3_settings("simulation")
folder <- c("UN-CCDs" = "UN_CCDs", "RK-CCDs" = "RK_CCDs", "KS-CCDs" = "KS_CCDs")
expect <- expand.grid(n = c(50, 100, 200, 500), d = c(2, 3, 5, 10, 20),
                      dist = c("uniform", "gaussian"), method = TR3_METHODS, stringsAsFactors = FALSE)
expect$driver <- sprintf("simulations/clustering/%s/Simulations/%s/%dd/%dd_235clusters_n%d.R",
                         folder[expect$method], tools::toTitleCase(expect$dist), expect$d, expect$d, expect$n)
if (nrow(s) != 120 || !setequal(s$driver, expect$driver) || anyDuplicated(s$driver)) {
  fail("settings_simulation.csv must list each of the 120 drivers exactly once")
}
on_disk <- list.files(here::here("simulations/clustering"), pattern = "_235clusters_n[0-9]+\\.R$",
                      recursive = TRUE)
if (!setequal(file.path("simulations/clustering", on_disk), expect$driver)) {
  fail("the drivers on disk are not the 120 expected ones")
}

templates <- list()
for (i in seq_len(nrow(s))) {
  r <- s[i, ]
  txt <- readLines(here::here(r$driver), warn = FALSE)
  chk <- function(ok, what) if (!isTRUE(ok)) fail("%s: %s", r$driver, what)

  for (v in c("n", "d", "iteN", "cls_dis", "r_min", "r_max")) {
    chk(identical(assigned(txt, v), r[[v]]), sprintf("%s = %s, table says %s", v, assigned(txt, v), r[[v]]))
  }
  chk(identical(assigned(txt, "sig"), if (r$dist == "gaussian") r$sig else NA_character_), "sig")
  chk(r$dist == "uniform" || r$sig != "NA", "sig missing for a Gaussian driver")
  chk(count(txt, sprintf("set.seed(%s)", r$seed)) == 1 && count(txt, "set.seed(") == 1, "seed")

  tab <- regmatches(txt, regexec('load_quantile_table\\(here::here\\("R/(NN|RK)-test_quantile/([^"]+)"\\)\\)', txt))
  tab <- vapply(Filter(length, tab), `[`, "", 3)
  if (r$method == "KS-CCDs") {
    chk(length(tab) == 0 && r$table == "", "KS-CCDs load no table")
  } else {
    chk(length(tab) == 1 && tab == r$table, sprintf("table %s, settings say %s", paste(tab, collapse = ","), r$table))
    vdir <- if (r$method == "UN-CCDs") "NN" else "RK"
    chk(startsWith(r$table, paste0(vdir, "-")) &&
          file.exists(here::here(sprintf("R/%s-test_quantile/%s", vdir, r$table))), "table file missing")
  }

  if (r$method == "UN-CCDs") {
    chk(identical(assigned(txt, "method"), sprintf('"%s"', r$direction)), "search direction")
    chk(count(txt, "UNMCCD_outlier(datax=x, simul=simul, method=method)") == 3, "UN call")
  } else {
    chk(is.na(assigned(txt, "method")) && r$direction == "", "direction set for a non-UN driver")
  }
  if (r$method == "RK-CCDs") {
    chk(count(squash(paste(txt, collapse = "")),
              "if(d<10){quant=0.99#thelevelofK-test}else{quant=0.999}") == 1, "RK level rule")
    chk(r$quant == if (as.numeric(r$d) < 10) "0.99" else "0.999", "RK level")
    chk(count(txt, "RKCCD_clustering(datax=x, simul=simul, quant=quant)") == 3, "RK call")
  } else {
    chk(r$quant == "", "level set for a non-RK driver")
  }
  if (r$method == "KS-CCDs") {
    chk(identical(assigned(txt, "delta_range"), r$delta_grid), "KS grid")
    chk(count(txt, "KSCCD_clustering(datax=x, m=delta)") == 3, "KS search call")
    chk(all(vapply(c(2, 3, 5), function(N) count(txt, sprintf("KSCCD_clustering(datax=x, m=delta_best_%dcls)", N)) == 1, TRUE)),
        "KS final call")
    chk(all(vapply(c(2, 3, 5), function(N) count(txt, sprintf("delta_best_%dcls = delta_range[which.max(Sil_seq_%dcls)]", N, N)) == 1, TRUE)),
        "KS selection rule")
  } else {
    chk(r$delta_grid == "", "grid set for a non-KS driver")
  }

  n1 <- regmatches(txt, regexec("^\\s*n1\\s*=\\s*round\\(n/N\\)\\s*([+-]\\s*[0-9]+)?\\s*$", txt))
  n1 <- vapply(Filter(length, n1), `[`, "", 2)
  chk(length(n1) == 3 && n1[1] == "" && n1[3] == "" &&
        as.numeric(gsub("\\s", "", n1[2])) == as.numeric(r$n1_adj_3cls), "size of the first of three clusters")
  for (N in c(2, 3, 5)) chk(sum(grepl(sprintf("^\\s*N\\s*=\\s*%d\\s*$", N), txt)) == 1, sprintf("N = %d block", N))

  gen <- if (r$dist == "uniform") "lapply(1:iteN,uniform_cls,n_cls=n_cls,d=d,centers=MU,noise=F,rad_range=c(r_min,r_max))" else
    "lapply(1:iteN,Gaussian_cls,n_cls=n_cls,d=d,centers=MU,sigma=sig,noise=F,rad_range=c(r_min,r_max))"
  chk(count(squash(paste(txt, collapse = "")), gen) == 3, "data-generation call")

  # structure: blank the values checked above, drop whitespace
  t2 <- sub("^\\s*(n|d|iteN|cls_dis|r_min|r_max|sig|method|delta_range)\\s*=.*$", "\\1=<V>", txt, perl = TRUE)
  t2 <- sub("^(\\s*n1\\s*=\\s*round\\(n/N\\))\\s*[+-]\\s*[0-9]+\\s*$", "\\1<ADJ>", t2)
  t2 <- gsub("load_quantile_table\\(.*\\)", "load_quantile_table(<T>)", t2)
  t2 <- gsub("setwd\\(.*\\)", "setwd(<W>)", t2)
  t2 <- gsub("save\\.image\\(.*\\)", "save.image(<S>)", t2)
  key <- paste(r$method, r$dist)
  templates[[key]] <- unique(c(templates[[key]], paste(squash(t2), collapse = "")))
}
for (k in names(templates)) {
  if (length(templates[[k]]) != 1) fail("%s: the drivers differ in more than the declared settings (%d texts)",
                                        k, length(templates[[k]]))
}
cat(sprintf("00_check_settings.R: settings_simulation.csv matches all %d drivers; one text per (method, dist): %s\n",
            nrow(s), paste(names(templates), collapse = "; ")))

# ---------------------------------------------------------------------------
# real-data scripts
# ---------------------------------------------------------------------------
rs <- tr3_settings("realdata")
for (sc in unique(rs$script)) {
  src <- readLines(here::here(sc), warn = FALSE)
  hdr <- grep("^##### \\S+ data #####\\s*$", src)
  blk <- sub("^##### (\\S+) data #####\\s*$", "\\1", src[hdr])
  sub_rs <- rs[rs$script == sc, ]
  if (!setequal(blk, sub_rs$dataset) || anyDuplicated(sub_rs$dataset)) {
    fail("%s: blocks (%s) and settings rows (%s) differ", sc, paste(blk, collapse = ","),
         paste(sub_rs$dataset, collapse = ","))
  }
  method <- unique(sub_rs$method)
  stopifnot(length(method) == 1)
  # the last block ends before the script's closing t2 = Sys.time() / stopCluster(cl)
  tail_line <- grep("^(t2 = Sys\\.time\\(\\)|stopCluster\\(cl\\))", src)
  tail_line <- tail_line[tail_line > max(hdr)][1]
  ends <- c(hdr[-1] - 1, if (is.na(tail_line)) length(src) else tail_line - 1)
  grid <- assigned(src[seq_len(hdr[1] - 1)], "delta_range")
  texts <- character(0)
  for (j in seq_along(hdr)) {
    b <- src[hdr[j]:ends[j]]
    r <- sub_rs[sub_rs$dataset == blk[j], ]
    chk <- function(ok, what) if (!isTRUE(ok)) fail("%s, %s: %s", sc, blk[j], what)
    tab <- regmatches(b, regexec('load_quantile_table\\(here::here\\("R/(NN|RK)-test_quantile/([^"]+)"\\)\\)', b))
    tab <- vapply(Filter(length, tab), `[`, "", 3)
    if (method == "KS-CCDs") {
      chk(length(tab) == 0 && r$table == "" && identical(grid, r$delta_grid), "KS grid / no table")
      chk(count(b, "Sil_seq = foreach(delta=delta_range,.packages = c(\"MASS\",\"cluster\")) %dopar% KSCCD_clustering(X, m=delta)$si") == 1,
          "KS search call")
      chk(count(b, "delta_best = delta_range[which.max(Sil_seq)]") == 1, "KS selection rule")
    } else {
      chk(length(tab) == 1 && tab == r$table, "table")
    }
    if (method == "UN-CCDs") {
      chk(identical(assigned(b, "method"), sprintf('"%s"', r$direction)), "search direction")
      chk(count(b, "UNCCD_clustering(datax=X, simul=simul, method=method)$label") == 1, "UN call")
    }
    if (method == "RK-CCDs") {
      chk(identical(assigned(b, "quant"), r$quant), "RK level")
      chk(count(b, "RKCCD_clustering(datax=X, simul=simul, quant=quant)$label") == 1, "RK call")
    }
    chk(count(b, "sil = mean(silhouette(") == 1 && count(b, "n_clusters = length(unique(") == 1, "metrics")
    t2 <- gsub(blk[j], "<DS>", b, fixed = TRUE)
    t2 <- sub("^\\s*(method|quant)\\s*=.*$", "\\1=<V>", t2)
    t2 <- gsub("load_quantile_table\\(.*\\)", "load_quantile_table(<T>)", t2)
    texts <- unique(c(texts, paste(squash(t2), collapse = "")))
  }
  if (length(texts) != 1) fail("%s: the blocks differ in more than the declared settings", sc)
}
cat(sprintf("00_check_settings.R: settings_realdata.csv matches all %d blocks of the 3 real-data scripts\n", nrow(rs)))
cat("00_check_settings.R: OK\n")
