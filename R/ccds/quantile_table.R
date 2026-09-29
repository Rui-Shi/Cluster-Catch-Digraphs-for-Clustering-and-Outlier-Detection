# quantile_table.R -- load a Monte Carlo quantile table, or generate it on the
# spot when the file is not there.
#
# The RK-based and NND-based spatial-randomness tests read their critical
# values from precomputed tables:
#   R/RK-test_quantile/RK-test-simul_<d>d_<level>%.RData   (Ripley's K envelope)
#   R/NN-test_quantile/NN-test-simul_<d>d_<level>%.RData   (mean/median NN distance)
# These tables are not in the public repository (the RK ones are 100-700 MB
# each), so the code that runs the detectors -- the simulation and real-data
# drivers, shared/harness.R and the revision scripts -- loads them through
# this file. (Audit scripts that inspect the shipped files themselves, e.g.
# the table-provenance checks in revision_experiments/tr1/, still need them.)
#
#   * load_quantile_table(path)  replaces load(path) in the simulation and
#     real-data drivers. If the file exists it is load()-ed exactly as before.
#     If not, `simul` is set to a placeholder that records the variant, d and
#     level; the construction fills it on first use at the data's own size
#     (the size is not known where the driver loads the table).
#   * ccd_quantile_table(variant, d, quant, n) returns a table covering at
#     least n points: the shipped file if long enough; a shipped file that is
#     too short, completed with generated entries past its end; a previously
#     generated table from the cache; or a new one. shared/harness.R
#     get_simul() uses it for missing or short tables.
#   * nnccd.radi() (UN_CCD.R) and ccd.Kest.edge.quantile() (RK_CCD_New.R)
#     call .ccd_qt_grow() only for placeholders and generated tables that are
#     shorter than the data, and the RK loop also when it actually reaches a
#     row past the end of a shipped table (where it used to stop with
#     "subscript out of bounds"). A shipped table is otherwise used exactly as
#     before, so published results do not change when the tables are present.
#
# The estimator is the original one. The per-iteration bodies below are
# verbatim transcriptions of the Monte Carlo loops in
# NNDestP.simpois.lower.quant (R/ccds/NN_Dist_Est.R) and
# KestP.simpois.edge.quantile (R/ccds/Kest.R); the reductions are the same
# quantile() calls. For d >= 342 the RK weight overflows in gamma(); the
# log-space-stable body of revision_experiments/tr2/01_gen_quantile_table.R is
# used there. Entry j of a table depends only on j points, so a table of
# extent n is statistically the first n entries of a larger one, and a short
# shipped table can be completed with generated entries.
#
# A generated table reproduces a published result up to Monte Carlo error,
# not bit for bit. RK tables drop the raw draw matrix Kest.m to stay small
# (the detectors read only $quan and $r).
#
# Settings (R option, else environment variable, else default):
#   ccd.quantile.niter  / CCD_QUANTILE_NITER  Monte Carlo iterations (as the
#                          shipped tables: NN 10000; RK 2000 for d < 50,
#                          10000 for d >= 50)
#   ccd.quantile.n      / CCD_QUANTILE_N      extent when the data size is not
#                          known (1000); an RK placeholder is first filled to
#                          min(n, this) and grown only if the search needs more
#   ccd.quantile.cores  / CCD_QUANTILE_CORES  worker processes (all cores - 1;
#                          1 inside a parallel worker)
#   ccd.quantile.seed   / CCD_QUANTILE_SEED   base seed (20260927; + d)
#   ccd.quantile.cache  / CCD_QUANTILE_CACHE  cache folder (default
#                          R/<RK|NN>-test_quantile/generated/)
#   ccd.quantile.generate / CCD_QUANTILE_GENERATE  FALSE stops instead of
#                          generating (TRUE)
# Cached files: <RK|NN>-test-simul_<d>d_<level>%_n<n>.RData (gitignored),
# reused by any later call that needs no more than n points.

if (!exists(".ccd_qt_memo", inherits = TRUE)) .ccd_qt_memo <- new.env(parent = emptyenv())

.ccd_qt_setting <- function(opt, env, default, as = as.numeric) {
  v <- getOption(opt)
  if (is.null(v)) {
    e <- Sys.getenv(env, "")
    v <- if (nzchar(e)) e else default
  }
  as(v)
}

.ccd_qt_generate_allowed <- function() {
  v <- getOption("ccd.quantile.generate")
  if (is.null(v)) {
    e <- toupper(Sys.getenv("CCD_QUANTILE_GENERATE", "TRUE"))
    v <- !(e %in% c("FALSE", "0", "NO"))
  }
  isTRUE(v)
}

.ccd_qt_default_niter <- function(variant, d) {
  def <- if (variant == "NN" || d >= 50) 10000 else 2000
  as.integer(.ccd_qt_setting("ccd.quantile.niter", "CCD_QUANTILE_NITER", def))
}

.ccd_qt_default_n <- function() as.integer(.ccd_qt_setting("ccd.quantile.n", "CCD_QUANTILE_N", 1000))

# extent to generate when the RK search needs row j of n: the default extent
# first (the search usually stops early), all n points only if it goes further
.ccd_qt_rk_target <- function(j, n) {
  dn <- .ccd_qt_default_n()
  if (j <= dn) min(n, dn) else n
}

# TRUE inside a forked (mclapply/doParallel on Unix) or PSOCK worker:
# generation there runs serially instead of opening a cluster of its own.
.ccd_qt_in_worker <- function() {
  child <- tryCatch(utils::getFromNamespace("isChild", "parallel")(), error = function(e) FALSE)
  isTRUE(child) || any(grepl("^MASTER=", commandArgs()))
}

# 0.99 -> "99", 0.999 -> "999", 0.9 -> "90", 0.85 -> "85" (the percentage with
# the dot removed, as in the shipped file names)
.ccd_qt_label <- function(q) {
  gsub("\\.", "", formatC(q * 100, format = "f", digits = 6, drop0trailing = TRUE))
}
# a level given as a file label ("99", "999") or a probability (0.99, "0.99")
.ccd_qt_prob <- function(quant) {
  if (is.character(quant)) {
    quant <- if (grepl(".", quant, fixed = TRUE)) as.numeric(quant) else as.numeric(paste0("0.", quant))
  }
  stopifnot(length(quant) == 1L, is.finite(quant), quant > 0, quant < 1)
  quant
}

# points a table covers at this level (0 for NULL, a placeholder, or an RK
# table without this level)
.ccd_qt_extent <- function(simul, variant, quant = NULL) {
  if (is.null(simul)) return(0L)
  if (variant == "NN") {
    return(as.integer(min(length(simul$average), length(simul$median))))
  }
  tab <- if (is.null(quant)) {
    if (length(simul$quan)) simul$quan[[1]] else NULL
  } else simul$quan[[as.character(quant)]]
  if (is.null(tab)) 0L else as.integer(nrow(tab))
}

.ccd_qt_tag <- function(simul, variant, d, quant, pending = FALSE) {
  attr(simul, "ccd_table") <- list(variant = variant, d = as.integer(d), quant = quant, pending = pending)
  simul
}

# placeholder for a table that is not on disk: extent 0, level recorded
.ccd_qt_placeholder <- function(variant, d, quant) {
  s <- if (variant == "NN") list(average = numeric(0), median = numeric(0)) else list(quan = list(), r = seq(0.1, 1, 0.1))
  .ccd_qt_tag(s, variant, d, quant, pending = TRUE)
}

# --- per-iteration bodies (verbatim; see header) ---------------------------

# NNDestP.simpois.lower.quant, SimuOnce (R/ccds/NN_Dist_Est.R)
.ccd_qt_nn_once <- function(n, d) {
  rpoisball.unit <- function(n,d){
    # inner ball and outer ball values
    r1 <- runif(n,0,1)^(1/d)
    norm.data <- matrix(MASS::mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
    data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
    data1 <- apply(data1,1,function(x) x*r1)
    return(data1)
  }
  data.simu.list = lapply(1:n, rpoisball.unit, d=d)
  NN.dist.temp = sapply(2:n, function(x){
    data.temp = data.simu.list[[x]]
    data.dist = as.matrix(dist(data.temp))
    diag(data.dist) = Inf
    NN.dist.ttemp = apply(data.dist, 1, min) # the Nearest Neighbor distance for each points
    NN.dist.ttemp.ave = mean(NN.dist.ttemp)
    NN.dist.ttemp.med = median(NN.dist.ttemp)
    return(c(NN.dist.ttemp.ave, NN.dist.ttemp.med))
  })
  NN.dist.temp.ave = c(0, NN.dist.temp[1,])
  NN.dist.temp.med = c(0, NN.dist.temp[2,])
  return(list(ave = NN.dist.temp.ave, med = NN.dist.temp.med))
}

# KestP.simpois.edge.quantile, SimuOnce (R/ccds/Kest.R)
.ccd_qt_rk_once <- function(m, d, rn) {
  r <- seq(1/rn,1,1/rn)
  rpoisball.unit <- function(n,d){
    # inner ball and outer ball values
    r1 <- runif(n,0,1)^(1/d)
    norm.data <- matrix(MASS::mvrnorm(n,rep(0,d),diag(d)),ncol=d,byrow=T)
    data1 <- apply(norm.data,1,function(x) x/sqrt(sum(x^2)))
    data1 <- apply(data1,1,function(x) x*r1)
    return(data1)
  }
  temp <- rpoisball.unit(m,d)

  # distances
  temp.dist <- as.matrix(dist(temp))

  # calculate weights for correction
  cons <- (sqrt(pi)*gamma((d+1)/2))/(2*gamma(d/2+1))
  integrand <- function(t) sin(t)^d
  ftemp <- sapply(temp.dist,function(t){
    return(integrate(integrand,0,acos(t/2))$value)
  },simplify = TRUE)
  ftemp <- matrix(ftemp,nrow=nrow(temp.dist),byrow = FALSE)
  ftemp <- cons*(1/ftemp)

  # analyze
  diag(temp.dist) <- Inf
  result <- sapply(r,function(x){
    Mtemp <- (temp.dist < x)
    Mtemp[lower.tri(Mtemp)] <- 0
    ftemp[lower.tri(ftemp)] <- 0
    Mtemp <- Mtemp*ftemp
    sumM <- cumsum(2*colSums(Mtemp))
    return(sumM/((1:m)*(1:m)))
  },simplify=TRUE)

  return(as.vector(result))
}

# log-space-stable RK body for d >= 342 (revision_experiments/tr2/
# 01_gen_quantile_table.R, KestP.simpois.edge.quantile.stable; validated there
# against the original at lower d by 01c_validate_rk_stable.R)
.ccd_qt_rk_once_stable <- function(m, d, rn) {
  r <- seq(1 / rn, 1, 1 / rn)
  rpoisball.unit <- function(n, d) {
    r1 <- runif(n, 0, 1)^(1 / d)
    norm.data <- matrix(MASS::mvrnorm(n, rep(0, d), diag(d)), ncol = d, byrow = T)
    data1 <- apply(norm.data, 1, function(x) x / sqrt(sum(x^2)))
    data1 <- apply(data1, 1, function(x) x * r1)
    return(data1)
  }
  temp <- rpoisball.unit(m, d)
  temp.dist <- as.matrix(dist(temp))
  log_cons <- log(sqrt(pi)) + lgamma((d + 1) / 2) - log(2) - lgamma(d / 2 + 1)
  log_ftemp_one <- function(t) {
    b <- acos(min(max(t / 2, -1), 1))
    if (b <= 0) return(Inf)
    lsb <- log(sin(b))
    I_scaled <- integrate(function(u) exp(d * (log(sin(u)) - lsb)), 0, b)$value
    logI <- d * lsb + log(I_scaled)
    return(log_cons - logI)
  }
  ftemp <- sapply(temp.dist, log_ftemp_one, simplify = TRUE)
  ftemp <- matrix(exp(ftemp), nrow = nrow(temp.dist), byrow = FALSE)
  diag(temp.dist) <- Inf
  result <- sapply(r, function(x) {
    Mtemp <- (temp.dist < x)
    Mtemp[lower.tri(Mtemp)] <- 0
    W <- ftemp
    W[Mtemp == 0] <- 0
    sumM <- cumsum(2 * colSums(W))
    return(sumM / ((1:m) * (1:m)))
  }, simplify = TRUE)
  return(as.vector(result))
}

# --- the Monte Carlo run -----------------------------------------------------

# One draw per iteration, serially under set.seed(seed) (cores = 1) or over a
# PSOCK cluster with clusterSetRNGStream(cl, seed) (cores > 1); reproducible
# for a fixed (niter, cores, seed). The caller's random-number state is
# restored afterwards.
.ccd_qt_draws <- function(once, args, niter, cores, seed) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  if (cores <= 1L) {
    set.seed(seed)
    return(lapply(seq_len(niter), function(i) do.call(once, args)))
  }
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, seed)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS)))
  fenv <- new.env(parent = globalenv())
  fenv$.once <- once; fenv$.args <- args
  fn <- function(i) do.call(.once, .args)
  environment(fn) <- fenv
  parallel::parLapply(cl, seq_len(niter), fn)
}

#' Generate one quantile table.
#' @param variant "RK" or "NN"; @param d dimension; @param quant level as a
#'   probability (0.99, 0.999, ...); @param n extent (points covered)
#' @return the `simul` object in the shape the detectors read:
#'   NN list(average, median); RK list(quan = list("<quant>" = n x rn), r)
ccd_generate_quantile_table <- function(variant, d, quant, n,
                                        niter = NULL, cores = NULL, seed = NULL) {
  variant <- match.arg(variant, c("RK", "NN"))
  n <- as.integer(n); d <- as.integer(d); quant <- .ccd_qt_prob(quant)
  stopifnot(n >= 2L, d >= 1L)
  niter <- as.integer(if (is.null(niter)) .ccd_qt_default_niter(variant, d) else niter)
  if (is.null(cores)) {
    cores <- if (.ccd_qt_in_worker()) 1L else
      .ccd_qt_setting("ccd.quantile.cores", "CCD_QUANTILE_CORES",
                      max(1L, parallel::detectCores() - 1L))
  }
  cores <- as.integer(max(1L, min(cores, niter)))
  if (is.null(seed)) seed <- .ccd_qt_setting("ccd.quantile.seed", "CCD_QUANTILE_SEED", 20260927) + d
  seed <- as.integer(seed)
  rn <- 10L   # every shipped RK table: r = seq(0.1, 1, by = 0.1)

  message(sprintf(paste0("[quantile table] generating %s d=%d level=%s%% for n=%d ",
                         "(niter=%d, cores=%d, seed=%d); this can take a while"),
                  variant, d, .ccd_qt_label(quant), n, niter, cores, seed))
  t0 <- Sys.time()
  if (variant == "NN") {
    out <- .ccd_qt_draws(.ccd_qt_nn_once, list(n = n, d = d), niter, cores, seed)
    ave <- do.call(rbind, lapply(out, `[[`, "ave"))
    med <- do.call(rbind, lapply(out, `[[`, "med"))
    # reduction as in NNDestP.simpois.lower.quant (lower tail, 1 - quant)
    quant.ave.lower <- sapply(1:n, function(x) quantile(ave[, x], 1 - quant))
    names(quant.ave.lower) <- NULL
    quant.med.lower <- sapply(1:n, function(x) quantile(med[, x], 1 - quant))
    names(quant.med.lower) <- NULL
    simul <- list(average = quant.ave.lower, median = quant.med.lower)
  } else {
    once <- if (d >= 342L) .ccd_qt_rk_once_stable else .ccd_qt_rk_once
    out <- .ccd_qt_draws(once, list(m = n, d = d, rn = rn), niter, cores, seed)
    Kest.m <- do.call(rbind, out)
    # reduction as in KestP.simpois.edge.quantile
    temp <- apply(Kest.m, 2, quantile, probs = as.numeric(quant))
    quan <- list()
    quan[[as.character(quant)]] <- matrix(temp, nrow = n)
    simul <- list(quan = quan, r = seq(1/rn, 1, 1/rn))
  }
  simul <- .ccd_qt_tag(simul, variant, d, quant)
  attr(simul, "generated") <- list(variant = variant, d = d, quant = quant, n = n,
                                   niter = niter, cores = cores, seed = seed,
                                   seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                                   date = format(Sys.time()), R = R.version.string)
  message(sprintf("[quantile table] done in %.0f s", attr(simul, "generated")$seconds))
  simul
}

.ccd_qt_dir <- function(variant) here::here(sprintf("R/%s-test_quantile", variant))

.ccd_qt_cache_dir <- function(dir) {
  v <- getOption("ccd.quantile.cache")
  if (is.null(v)) v <- Sys.getenv("CCD_QUANTILE_CACHE", "")
  if (nzchar(v)) v else file.path(dir, "generated")
}

# entries 1..L of `short` followed by entries L+1..n of `long` (both at the
# same level); rows are independent, so the result is a valid table
.ccd_qt_splice <- function(short, long, variant, quant, n) {
  L <- .ccd_qt_extent(short, variant, quant)
  if (variant == "NN") {
    out <- list(average = c(short$average[seq_len(L)], long$average[(L + 1):n]),
                median  = c(short$median[seq_len(L)],  long$median[(L + 1):n]))
  } else {
    key <- as.character(quant)
    out <- short
    out$Kest.m <- NULL
    out$quan[[key]] <- rbind(short$quan[[key]][seq_len(L), , drop = FALSE],
                             long$quan[[key]][(L + 1):n, , drop = FALSE])
  }
  out
}

#' Return a table covering at least n points (see header).
#' @param quant level as a probability or a file label ("99", "999")
#' @param n points to cover; NULL = the shipped file as is, or the default
#'   extent (ccd.quantile.n) when it has to be generated
#' @param dir directory of the shipped tables (default R/<variant>-test_quantile)
#' @return list(simul, file, generated)
ccd_quantile_table <- function(variant, d, quant, n = NULL, dir = NULL,
                               niter = NULL, cores = NULL, seed = NULL) {
  variant <- match.arg(variant, c("RK", "NN"))
  d <- as.integer(d); quant <- .ccd_qt_prob(quant)
  label <- .ccd_qt_label(quant)
  if (is.null(dir)) dir <- .ccd_qt_dir(variant)
  if (!is.null(n) && is.na(n)) n <- NULL
  need <- if (is.null(n)) 0L else as.integer(n)
  fname <- sprintf("%s-test-simul_%dd_%s%%.RData", variant, d, label)
  cache <- .ccd_qt_cache_dir(dir)

  key <- paste(variant, d, label, normalizePath(dir, mustWork = FALSE), cache, sep = "|")
  hit <- .ccd_qt_memo[[key]]
  if (!is.null(hit) && .ccd_qt_extent(hit$simul, variant, quant) >= need) return(hit)

  read_simul <- function(path) { e <- new.env(); load(path, envir = e); get("simul", envir = e) }
  keep <- function(res) { assign(key, res, envir = .ccd_qt_memo); res }

  find_cached <- function(size) {
    pat <- sprintf("^%s-test-simul_%dd_%s%%_n([0-9]+)\\.RData$", variant, d, label)
    cached <- list.files(cache, pattern = pat)
    if (!length(cached)) return(NULL)
    ext <- as.integer(sub(pat, "\\1", cached))
    ok <- which(ext >= size)
    if (!length(ok)) return(NULL)
    f <- file.path(cache, cached[ok][which.min(ext[ok])])
    list(simul = .ccd_qt_tag(read_simul(f), variant, d, quant), file = f)
  }
  # a generated table covering `size` points: from the cache, else new
  generated <- function(size) {
    hit <- find_cached(size)
    if (!is.null(hit)) return(hit)
    if (!.ccd_qt_generate_allowed()) {
      stop(sprintf(paste0("quantile table %s is missing or covers fewer than n = %d points, ",
                          "and generation is switched off (ccd.quantile.generate / CCD_QUANTILE_GENERATE)."),
                   file.path(dir, fname), size))
    }
    # one process generates a given (variant, d, level) at a time; parallel
    # workers and concurrent jobs wait for it and then reuse its file
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    lock <- file.path(cache, sprintf(".lock_%s-test-simul_%dd_%s", variant, d, label))
    told <- FALSE
    while (!dir.create(lock, showWarnings = FALSE)) {
      age <- as.numeric(difftime(Sys.time(), file.mtime(lock), units = "hours"))
      if (is.na(age) || age > 48) { unlink(lock, recursive = TRUE); next }   # left by a crashed run
      if (!told) { message("[quantile table] waiting for another process generating ", basename(lock),
                           " (delete that folder if no run is active)"); told <- TRUE }
      Sys.sleep(10)
      hit <- find_cached(size)
      if (!is.null(hit)) return(hit)
    }
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    hit <- find_cached(size)
    if (!is.null(hit)) return(hit)
    simul <- ccd_generate_quantile_table(variant, d, quant, size, niter, cores, seed)
    f <- file.path(cache, sprintf("%s-test-simul_%dd_%s%%_n%d.RData", variant, d, label, size))
    tmp <- tempfile(pattern = ".partial_", tmpdir = cache, fileext = ".RData")
    save(simul, file = tmp)
    if (!file.rename(tmp, f)) unlink(tmp)
    list(simul = simul, file = f)
  }

  shipped <- file.path(dir, fname)
  if (file.exists(shipped)) {
    s <- read_simul(shipped)
    L <- .ccd_qt_extent(s, variant, quant)
    if (L >= need) return(keep(list(simul = s, file = shipped, generated = FALSE)))
    if (L > 0L) {
      # too short: keep the shipped entries, generate the ones past its end
      g <- generated(need)
      sp <- .ccd_qt_tag(.ccd_qt_splice(s, g$simul, variant, quant, need), variant, d, quant)
      return(keep(list(simul = sp, file = g$file, generated = TRUE)))
    }
  }
  g <- generated(if (need >= 2L) need else .ccd_qt_default_n())
  keep(list(simul = g$simul, file = g$file, generated = TRUE))
}

#' Called by nnccd.radi() and ccd.Kest.edge.quantile(): a table covering
#' `need` points at the level of `simul` (NN) or `quant` (RK). Placeholders
#' and generated tables are replaced; a shipped table keeps its entries and
#' is completed past its end.
.ccd_qt_grow <- function(simul, variant, d, need, quant = NULL) {
  info <- attr(simul, "ccd_table")
  if (is.null(quant)) quant <- info$quant
  res <- ccd_quantile_table(variant, d, quant, n = need)$simul
  if (!is.null(info) || .ccd_qt_extent(simul, variant, quant) == 0L) return(res)
  .ccd_qt_tag(.ccd_qt_splice(simul, res, variant, quant, need), variant, d, quant)
}

#' Drop-in replacement for load(path) on a quantile-table file.
#' If the file exists it is load()-ed into envir exactly as load() would.
#' Otherwise `simul` in envir becomes a placeholder that the construction
#' fills at the data's own size. The file name gives variant, d and level
#' (<RK|NN>-test-simul_<d>d_<level>%.RData).
load_quantile_table <- function(path, envir = parent.frame()) {
  if (file.exists(path)) return(invisible(load(path, envir = envir)))
  pat <- "^(RK|NN)-test-simul_([0-9]+)d_([0-9]+)%\\.RData$"
  b <- basename(path)
  if (!grepl(pat, b)) stop("load_quantile_table(): not a quantile-table file name: ", path)
  variant <- sub(pat, "\\1", b)
  d <- as.integer(sub(pat, "\\2", b))
  quant <- .ccd_qt_prob(sub(pat, "\\3", b))
  # a driver written for another machine: this repository's copy, if any
  local_copy <- file.path(.ccd_qt_dir(variant), b)
  if (file.exists(local_copy)) return(invisible(load(local_copy, envir = envir)))
  if (!.ccd_qt_generate_allowed()) stop("quantile table not found and generation is switched off: ", path)
  message(sprintf("[quantile table] %s not found; it will be generated when the data size is known", b))
  assign("simul", .ccd_qt_placeholder(variant, d, quant), envir = envir)
  invisible("simul")
}
