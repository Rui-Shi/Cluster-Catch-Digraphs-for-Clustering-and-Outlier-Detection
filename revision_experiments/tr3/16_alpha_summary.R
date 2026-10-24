#!/usr/bin/env Rscript
# revision_experiments/tr3/16_alpha_summary.R
#
# Step 3 of ALPHA_PROTOCOL.md: reproduction check first, then the criteria.
#
# USAGE (from the clone root, in PowerShell)
#   Rscript revision_experiments/tr3/16_alpha_summary.R
#
# 1. Reproduction check. At the script's own level, config P (user_knowledge and
#    wholesale) and config A (the other seven, which P shares with A) must equal
#    the `published` UN-CCDs rows of results/tr3/real.csv to 15 significant
#    digits (ARI, Sil, k_hat). Since the rule was adopted (2026-10-23), real.csv
#    holds the new rows for user_knowledge, wholesale and asymmetric, so the
#    rows from before are read from results/tr3/_before_alpha_adopt_20261023/
#    when that folder exists. Written to results/tr3/alpha_repro_check.csv. If
#    a cell differs the script stops with exit code 1 and does not analyse.
#    The same file records a second check: wherever config A used a generated
#    table, config G used the same file and must give the identical cell.
# 2. Criteria per config and level, and the reference (the script's own levels),
#    written to results/tr3/alpha_summary.csv and printed.
#
# Composition. A config that covers only some data sets (P, S) is scored over
# all data sets, taking its own rows where it has them and config A's otherwise.
# Failure rule: a (config, data set) with an errored cell at any of the five
# levels is dropped from that config's rankings (all five levels), and the
# failed cells are listed. The reference is scored over the same data sets.
# Mean rank of the five levels: per data set, rank by ARI, 1 = best, ties
# averaged, then the mean over the data sets. `mean_rank_with_ref` ranks the
# five levels and the reference together (6 candidates).
#
# Switch for tests: TR3_ALPHA_RES (outputs and real_alpha.csv; real.csv is read
# from results/tr3 either way).

source(here::here("revision_experiments/tr3/tr3_lib.R"))
RES <- Sys.getenv("TR3_ALPHA_RES", TR3_RES)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

DS <- data.frame(
  dataset = c("iris", "seeds", "user_knowledge", "wholesale", "breast_cancer",
              "aggregation", "asymmetric", "R15", "D31"),
  true_k  = c(3, 3, 4, 3, 2, 7, 5, 15, 31),
  script_level = c(10, 5, 5, 5, 5, 15, 20, 15, 15),
  stringsAsFactors = FALSE)
PAPER7 <- c("iris", "seeds", "user_knowledge", "wholesale", "asymmetric", "R15", "D31")
ALL9 <- DS$dataset
LEVELS <- c(20, 15, 10, 5, 1)
CONFIGS <- c("A", "P", "G", "S")

ra <- tr3_read(file.path(RES, "real_alpha.csv"))
if (is.null(ra)) stop("16_alpha_summary.R: no real_alpha.csv in ", RES)
ra$error[is.na(ra$error)] <- ""
cell <- function(cfg, lev, ds) {
  r <- ra[ra$config == cfg & ra$level == as.character(lev) & ra$dataset == ds, , drop = FALSE]
  if (nrow(r) != 1) return(NULL)
  r
}
num <- function(x) suppressWarnings(as.numeric(x))
g15 <- function(x) ifelse(is.na(x), "NA", sprintf("%.15g", x))

# ---- 1. reproduction check ---------------------------------------------------
PUB_DIR <- tr3_path("_before_alpha_adopt_20261023")
if (!dir.exists(PUB_DIR)) PUB_DIR <- TR3_RES
pub <- tr3_read_done(file.path(PUB_DIR, "real.csv"), file.path(PUB_DIR, "real_done.csv"),
                     c("variant", "method", "dataset"))
pub <- pub[pub$variant == "published" & pub$method == "UN-CCDs", ]
chk <- list()
add <- function(check, ds, cfg, lev, metric, a, b) {
  chk[[length(chk) + 1]] <<- data.frame(check = check, dataset = ds, config = cfg, level = lev,
                                        metric = metric, reference = a, rerun = b,
                                        equal = identical(a, b), stringsAsFactors = FALSE)
}
for (i in seq_len(nrow(DS))) {
  ds <- DS$dataset[i]; lev <- DS$script_level[i]
  cfg <- if (ds %in% c("user_knowledge", "wholesale")) "P" else "A"
  p <- pub[pub$dataset == ds, ]
  r <- cell(cfg, lev, ds)
  for (m in c("ari", "sil", "k_hat")) {
    a <- if (nrow(p) == 1) g15(num(p[[m]])) else "MISSING_PUBLISHED"
    b <- if (is.null(r)) "MISSING_CELL" else if (nzchar(r$error)) paste("ERROR:", r$error) else g15(num(r[[m]]))
    add("published_vs_rerun", ds, cfg, lev, m, a, b)
  }
}
for (lev in LEVELS) for (ds in ALL9) {
  a <- cell("A", lev, ds); g <- cell("G", lev, ds)
  if (is.null(a) || is.null(g) || a$source != "generated") next
  for (m in c("ari", "sil", "k_hat", "error")) add("G_equals_A", ds, "G_vs_A", lev, m, a[[m]], g[[m]])
}
chk <- do.call(rbind, chk)
utils::write.csv(chk, file.path(RES, "alpha_repro_check.csv"), row.names = FALSE)
bad <- chk[!chk$equal, ]
cat(sprintf("reproduction check: %d comparisons, %d differ\n", nrow(chk), nrow(bad)))
if (nrow(bad)) {
  print(bad)
  cat("STOP: a check failed; no analysis. See alpha_repro_check.csv\n")
  quit(save = "no", status = 1)
}

# ---- 2. criteria -------------------------------------------------------------
compose <- function(cfg) {
  out <- list()
  for (lev in LEVELS) for (ds in ALL9) {
    r <- cell(cfg, lev, ds)
    if (is.null(r)) r <- cell("A", lev, ds)
    failed <- is.null(r) || nzchar(r$error)
    out[[length(out) + 1]] <- data.frame(
      level = lev, dataset = ds, failed = failed,
      error = if (is.null(r)) "missing cell" else r$error,
      ari = if (failed) NA_real_ else num(r$ari), sil = if (failed) NA_real_ else num(r$sil),
      k_hat = if (failed) NA_real_ else num(r$k_hat), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
reference <- function() {
  out <- list()
  for (i in seq_len(nrow(DS))) {
    ds <- DS$dataset[i]; lev <- DS$script_level[i]
    cfg <- if (ds %in% c("user_knowledge", "wholesale")) "P" else "A"
    r <- cell(cfg, lev, ds)
    out[[i]] <- data.frame(dataset = ds, ari = num(r$ari), sil = num(r$sil), k_hat = num(r$k_hat),
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}
REF <- reference()

summ <- list()
perds <- list()
for (cfg in CONFIGS) {
  C <- compose(cfg)
  C$true_k <- DS$true_k[match(C$dataset, DS$dataset)]
  bad_ds <- unique(C$dataset[C$failed])
  valid <- setdiff(ALL9, bad_ds)
  d7 <- intersect(PAPER7, valid); d9 <- intersect(ALL9, valid)
  failed_txt <- if (any(C$failed)) paste(sprintf("%s@%s%%: %s", C$dataset[C$failed], C$level[C$failed],
                                                  C$error[C$failed]), collapse = " | ") else ""
  cat(sprintf("\n== config %s: data sets scored (paper 7): %s | all 9: %d | failed cells: %s\n", cfg,
              paste(d7, collapse = ","), length(d9), if (nzchar(failed_txt)) failed_txt else "none"))
  # ARI matrix over the paper-7 data sets: rows = 5 levels (+ reference), cols = data sets
  A7 <- sapply(d7, function(ds) sapply(LEVELS, function(l) C$ari[C$level == l & C$dataset == ds]))
  rownames(A7) <- LEVELS
  R7 <- sapply(d7, function(ds) REF$ari[REF$dataset == ds])
  A7r <- rbind(A7, script = R7)
  rank5 <- apply(A7, 2, function(v) rank(-v, ties.method = "average"))
  rank6 <- apply(A7r, 2, function(v) rank(-v, ties.method = "average"))
  one <- function(lev, ari, sil, khat, k_true, r5, r6) {
    # k_hat = 1 counts as Sil 0 and is flagged
    s0 <- ifelse(is.na(sil) & khat == 1, 0, sil)
    data.frame(config = cfg, level = lev,
               n_ds7 = length(d7), ds7 = paste(d7, collapse = " "),
               mean_ari7 = mean(ari), mean_sil7 = mean(s0), k_correct7 = sum(khat == k_true),
               mean_rank7 = r5, mean_rank_with_ref = r6, k1_count7 = sum(khat == 1),
               stringsAsFactors = FALSE)
  }
  rows <- list()
  for (lev in LEVELS) {
    sub <- C[C$level == lev & C$dataset %in% d7, ]
    sub <- sub[match(d7, sub$dataset), ]
    rows[[length(rows) + 1]] <- one(as.character(lev), sub$ari, sub$sil, sub$k_hat, sub$true_k,
                                    mean(rank5[as.character(lev), ]), mean(rank6[as.character(lev), ]))
    sub9 <- C[C$level == lev & C$dataset %in% d9, ]
    rows[[length(rows)]]$mean_ari9 <- mean(sub9$ari)
    rows[[length(rows)]]$n_ds9 <- length(d9)
    rows[[length(rows)]]$failed_cells <- failed_txt
  }
  rf <- REF[match(d7, REF$dataset), ]
  rows[[length(rows) + 1]] <- one("script", rf$ari, rf$sil, rf$k_hat, DS$true_k[match(d7, DS$dataset)],
                                  NA_real_, mean(rank6["script", ]))
  r9 <- REF[REF$dataset %in% d9, ]
  rows[[length(rows)]]$mean_ari9 <- mean(r9$ari)
  rows[[length(rows)]]$n_ds9 <- length(d9)
  rows[[length(rows)]]$failed_cells <- failed_txt
  S <- do.call(rbind, rows)
  lv <- S$level != "script"
  best <- max(S$mean_ari7[lv])
  S$best_level <- ifelse(lv & S$mean_ari7 == best, "BEST", "")
  S$tied_within_0.001 <- ifelse(lv & S$mean_ari7 >= best - 0.001, "yes", "")
  summ[[cfg]] <- S
  # per-data-set cells of this composite
  P <- C[, c("level", "dataset", "ari", "sil", "k_hat", "failed", "error")]
  P$config <- cfg
  perds[[cfg]] <- P
  show <- S[, c("level", "mean_ari7", "mean_sil7", "k_correct7", "mean_rank7", "mean_rank_with_ref",
                "mean_ari9", "k1_count7", "best_level", "tied_within_0.001")]
  num_cols <- c("mean_ari7", "mean_sil7", "mean_rank7", "mean_rank_with_ref", "mean_ari9")
  show[num_cols] <- lapply(show[num_cols], function(x) round(x, 4))
  print(show, row.names = FALSE)
  tied <- S$level[lv & S$mean_ari7 >= best - 0.001]
  cat(sprintf("config %s: highest mean ARI over 7 data sets at alpha = %s%% (%.4f); within 0.001 of it: %s\n",
              cfg, S$level[lv & S$mean_ari7 == best][1], best, paste0(tied, "%", collapse = ", ")))
}
summary_df <- do.call(rbind, summ)
for (col in c("mean_ari7", "mean_sil7", "mean_rank7", "mean_rank_with_ref", "mean_ari9"))
  summary_df[[col]] <- fmt17(summary_df[[col]])
utils::write.csv(summary_df, file.path(RES, "alpha_summary.csv"), row.names = FALSE)

pd <- do.call(rbind, perds)
pd$ari <- fmt17(pd$ari); pd$sil <- fmt17(pd$sil)
utils::write.csv(pd, file.path(RES, "alpha_cells_composed.csv"), row.names = FALSE)

cat("\nConfig A per data set (ARI / Sil / k_hat):\n")
PA <- perds[["A"]]
for (ds in ALL9) {
  x <- PA[PA$dataset == ds, ]
  cat(sprintf("%-15s %s\n", ds, paste(sprintf("%s%%: %.3f/%.3f/%s", x$level, x$ari, x$sil, x$k_hat), collapse = "  ")))
}
cat("\n16_alpha_summary.R done\n")
