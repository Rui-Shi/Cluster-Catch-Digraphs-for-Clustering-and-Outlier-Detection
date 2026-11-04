#!/usr/bin/env Rscript
# revision_experiments/tr3/24_alpha_supplement.R
#
# Two supplement tables from runs that already exist (nothing is rerun):
#   1. UN-CCDs on the seven real data sets of the paper at alpha = 20, 15, 10,
#      5 and 1%, with the settings the paper now uses (tables of the data's own
#      d, ascending search for every d < 10). Source: real_alpha.csv, config A,
#      except seeds, which comes from config S (config A ran seeds descending).
#   2. The same cells with a shipped table against an independently generated
#      table of the same d and level (niter 10000, seed 20260927 + d; config G).
#      seeds is left out: its config G run used the descending search.
#      Plus the largest relative difference of the table entries
#      (alpha_table_check.csv) for the d of these data sets.
# Checks: at the levels the paper's rule sets (15% for d = 2, 10% for d = 3, 4,
# 5% for 5 <= d <= 9), every cell of table 1 must equal the UN-CCD row of
# real.csv (ARI and Sil within 1e-12, same k_hat).
#
# USAGE (from the clone root): Rscript revision_experiments/tr3/24_alpha_supplement.R
# OUTPUT: results/tr3/supp_alpha.tex (two tabular bodies) and a summary on stdout.

source(here::here("revision_experiments/tr3/tr3_lib.R"))

PAPER  <- c("iris", "seeds", "user_knowledge", "wholesale", "asymmetric", "R15", "D31")
LABEL  <- c(iris = "iris", seeds = "seeds", user_knowledge = "knowledge", wholesale = "wholesale",
            asymmetric = "asymmetric", R15 = "R15", D31 = "D31")
TRUE_K <- c(iris = 3, seeds = 3, user_knowledge = 4, wholesale = 3, asymmetric = 5, R15 = 15, D31 = 31)
LEVELS <- c("20", "15", "10", "5", "1")
rule_level <- function(d) if (d == 2) "15" else if (d <= 4) "10" else "5"

ra   <- utils::read.csv(tr3_path("real_alpha.csv"), colClasses = "character")
real <- utils::read.csv(tr3_path("real.csv"), colClasses = "character")
tchk <- utils::read.csv(tr3_path("alpha_table_check.csv"), colClasses = "character")

cell <- function(config, ds, lv) {
  r <- ra[ra$config == config & ra$dataset == ds & ra$level == lv, ]
  if (nrow(r) != 1) stop("expected one row: ", config, " ", ds, " ", lv)
  if (nzchar(r$error)) stop("failed cell: ", config, " ", ds, " ", lv)
  r
}
current <- function(ds, lv) cell(if (ds == "seeds") "S" else "A", ds, lv)

f3 <- function(x) { x <- as.numeric(x); s <- sprintf("%.3f", x); s[s == "-0.000"] <- "0.000"; s }
red <- function(s) sprintf("\\textcolor{red}{%s}", s)

# --- check against the published UN-CCD rows --------------------------------
for (ds in PAPER) {
  r <- current(ds, rule_level(as.integer(current(ds, "5")$d)))
  p <- real[real$variant == "published" & real$method == "UN-CCDs" & real$dataset == ds, ]
  stopifnot(nrow(p) == 1,
            abs(as.numeric(r$ari) - as.numeric(p$ari)) < 1e-12,
            abs(as.numeric(r$sil) - as.numeric(p$sil)) < 1e-12,
            r$k_hat == p$k_hat)
}
cat("check passed: the rule's cells equal real.csv for all seven data sets\n")

# --- table 1: alpha sweep ----------------------------------------------------
lines1 <- character(0)
summ <- data.frame(level = LEVELS, mean_ari = NA_real_, k_correct = NA_integer_)
for (ds in PAPER) {
  d <- as.integer(current(ds, "5")$d)
  cells <- vapply(LEVELS, function(lv) {
    r <- current(ds, lv)
    s <- c(f3(r$ari), f3(r$sil), r$k_hat)
    if (lv == rule_level(d)) s <- sprintf("\\textbf{%s}", s)
    paste(red(s), collapse = " & ")
  }, "")
  lines1 <- c(lines1, sprintf("  %s & %s & %s & %s \\\\ \\hline",
                              red(LABEL[[ds]]), red(d), red(TRUE_K[[ds]]), paste(cells, collapse = " & ")))
}
for (i in seq_along(LEVELS)) {
  rs <- lapply(PAPER, function(ds) current(ds, LEVELS[i]))
  summ$mean_ari[i]  <- mean(vapply(rs, function(r) as.numeric(r$ari), 0))
  summ$k_correct[i] <- sum(vapply(seq_along(PAPER), function(j) as.integer(rs[[j]]$k_hat) == TRUE_K[[PAPER[j]]], TRUE))
}
rule_rows <- lapply(PAPER, function(ds) current(ds, rule_level(as.integer(current(ds, "5")$d))))
rule_mean <- mean(vapply(rule_rows, function(r) as.numeric(r$ari), 0))
rule_k    <- sum(vapply(seq_along(PAPER), function(j) as.integer(rule_rows[[j]]$k_hat) == TRUE_K[[PAPER[j]]], TRUE))
lines1 <- c(lines1, sprintf("  \\multicolumn{3}{|c|}{%s} & %s \\\\ \\hline", red("Mean \\textit{ARI}; $\\hat{k}=k$"),
  paste(vapply(seq_along(LEVELS), function(i)
    sprintf("\\multicolumn{3}{c|}{%s}", red(sprintf("%s; %d of 7", f3(summ$mean_ari[i]), summ$k_correct[i]))), ""),
    collapse = " & ")))

# --- table 2: shipped against regenerated tables -----------------------------
lines2 <- character(0)
n_cells <- 0; n_same <- 0
for (ds in setdiff(PAPER, "seeds")) {
  for (lv in LEVELS) {
    a <- cell("A", ds, lv); g <- cell("G", ds, lv)
    if (a$source != "shipped") next
    if (g$source != "generated") stop("config G did not use a generated table: ", ds, " ", lv)
    n_cells <- n_cells + 1
    same <- a$ari == g$ari && a$sil == g$sil && a$k_hat == g$k_hat
    n_same <- n_same + same
    s <- c(LABEL[[ds]], paste0(lv, "\\%"), f3(a$ari), f3(g$ari), f3(a$sil), f3(g$sil), a$k_hat, g$k_hat)
    if (lv == rule_level(as.integer(a$d))) s[2] <- sprintf("\\textbf{%s}", s[2])
    lines2 <- c(lines2, sprintf("  %s \\\\ \\hline", paste(red(s), collapse = " & ")))
  }
}
ds_d <- unique(vapply(PAPER, function(ds) current(ds, "5")$d, ""))
tk <- tchk[tchk$d %in% ds_d, ]
tk$worst <- pmax(as.numeric(tk$max_rel_diff_average), as.numeric(tk$max_rel_diff_median))

out <- tr3_path("supp_alpha.tex")
writeLines(c("% generated by revision_experiments/tr3/24_alpha_supplement.R; do not edit",
             "% table 1: alpha sweep", lines1, "", "% table 2: shipped vs regenerated tables", lines2), out)

cat("\nalpha sweep, seven data sets (mean ARI; k_hat = true k):\n")
print(transform(summ, mean_ari = round(mean_ari, 4)), row.names = FALSE)
cat(sprintf("rule by d: mean ARI %.4f; k_hat = true k on %d of 7\n", rule_mean, rule_k))
cat(sprintf("\nregenerated tables: %d cells with a shipped table (seeds left out); %d identical in ARI, Sil and k_hat\n",
            n_cells, n_same))
cat(sprintf("largest relative difference of table entries for d in {%s}: %.4f (d = %s, %s); next largest %.4f\n",
            paste(sort(as.integer(ds_d)), collapse = ", "), max(tk$worst),
            tk$d[which.max(tk$worst)], tk$level[which.max(tk$worst)], sort(tk$worst, decreasing = TRUE)[2]))
cat("written:", out, "\n")
