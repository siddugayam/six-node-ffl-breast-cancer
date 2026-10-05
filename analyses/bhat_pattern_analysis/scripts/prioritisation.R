## Section P of SETTINGS.md: the paper's prioritisation (scripts/06_node_prioritisation/01_prioritisation.R) with Bhat FFL participation.
## The paper's script is read and evaluated with two changes only (the F1 precedent of analyses/six_node_pattern):
##   its output folder, and one inserted line after its FFL count 'fps' is computed, which replaces that count by the
##   node's number of instances in one Bhat network (the n_ffl column of the network's node file; 0 if absent).
## A run with nothing inserted must reproduce the stored results/v5/node_prioritisation_full.csv and top10_per_class.csv
## byte for byte. Ranks are recomputed among the rankable nodes (priority not NA), overall and within class, ties 'min',
## as scripts/06_node_prioritisation/40_node_compendium_assemble.py does. Rank change = paper rank - version rank (positive = moved up).
## usage: Rscript prioritisation.R <analysis root> <folder of node files> <node-file suffix> <output folder> <label>
##   node files: <folder>/<network><suffix> with columns node, type, n_ffl (the network deposit: suffix _nodes.tsv)
suppressPackageStartupMessages({ library(data.table); library(tools) })
a <- commandArgs(TRUE); REV <- a[1]; NDIR <- a[2]; SUF <- a[3]; OUT <- a[4]; LABEL <- a[5]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(OUT, "prioritisation.log"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG); flush(LOG) }
SRC <- file.path(REV, "scripts", "06_node_prioritisation", "01_prioritisation.R")
code <- readLines(SRC)
r5 <- grep('R5 <- file.path(REV,"results/v5")', code, fixed = TRUE); stopifnot(length(r5) == 1)
fl <- grep("fps <- allroles[, .(ffl_cores=.N), by=name]", code, fixed = TRUE); stopifnot(length(fl) == 1)
run_version <- function(dir, fps_new) {
  cc <- code
  cc[r5] <- sub('R5 <- file.path(REV,"results/v5")', 'R5 <- OUTDIR', cc[r5], fixed = TRUE)
  if (!is.null(fps_new)) cc <- append(cc, "fps <- FPS_NEW   # inserted: Bhat FFL participation", after = fl)
  env <- new.env(); env$OUTDIR <- dir; env$FPS_NEW <- fps_new
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  capture.output(eval(parse(text = cc), envir = env))
  fread(file.path(dir, "node_prioritisation_full.csv"))
}
NETS <- as.vector(outer(c("miRNA_FFL", "TF_FFL", "composite_FFL"), 3:6, function(k, n) paste0(n, "node_", k)))
## the paper's run
P0 <- run_version(file.path(OUT, "runs", "paper"), NULL)
same <- c(md5sum(file.path(OUT, "runs", "paper", "node_prioritisation_full.csv")) == md5sum(file.path(REV, "results/v5/node_prioritisation_full.csv")),
          md5sum(file.path(OUT, "runs", "paper", "top10_per_class.csv")) == md5sum(file.path(REV, "results/v5/top10_per_class.csv")))
chk <- sprintf("CHECK P: the unmodified run reproduces results/v5/node_prioritisation_full.csv byte for byte: %s; top10_per_class.csv: %s", same[1], same[2])
say(chk); stopifnot(all(same))
DOM <- c("T_topology", "E_evidence", "D_expression", "R_replication", "M_multiomic", "C_clinical", "F_functional")
rk <- function(P) {                                      # ranks among rankable nodes, ties 'min'
  X <- P[!is.na(priority)]
  X[, rank_all := frank(-priority, ties.method = "min")]
  X[, rank_class := frank(-priority, ties.method = "min"), by = type]
  X
}
R0 <- rk(P0); stopifnot(nrow(R0) == 583)
top0 <- fread(file.path(REV, "results/v5/top10_per_class.csv"))[, .(type, name)]
## Table S4 within-class rank against the ranks among rankable nodes
S4 <- fread(file.path(REV, "results/v5/tables/TableS4_node_prioritisation_full.csv"))
s4 <- merge(S4[, .(name, type, rank_within_type)], R0[, .(name, rank_class)], by = "name")
s4[, off := suppressWarnings(as.integer(rank_within_type)) - rank_class]
s4chk <- s4[, .(nodes = .N, offsets = paste(names(table(off)), table(off), sep = " x", collapse = "; ")), by = type]
chk <- c(chk, sprintf("CHECK Table S4 rank_within_type minus the rank among rankable nodes, %s (%d rankable): %s", s4chk$type, s4chk$nodes, s4chk$offsets))
TAB <- list(); SUMM <- list(); MOV <- list(); T4 <- list()
for (net in NETS) {
  nd <- fread(file.path(NDIR, paste0(net, SUF)))
  fps_new <- nd[, .(name = node, ffl_cores = as.integer(n_ffl))]
  P1 <- run_version(file.path(OUT, "runs", net), fps_new)
  R1 <- rk(P1); stopifnot(nrow(R1) == 583, setequal(R1$name, R0$name))
  M <- merge(R0[, .(name, type, paper_score = priority, paper_rank_overall = rank_all, paper_rank_class = rank_class)],
             R1[, c("name", "ffl_cores", "T_topology", "priority", "rank_all", "rank_class", DOM[-1]), with = FALSE], by = "name")
  setnames(M, c("ffl_cores", "priority", "rank_all", "rank_class"), c("ffl_count_used", "score", "rank_overall", "rank_class"))
  M[, in_this_network := ifelse(name %in% nd$node, "yes", "no")]
  M[, rank_change_overall := paper_rank_overall - rank_overall]; M[, rank_change_class := paper_rank_class - rank_class]
  setcolorder(M, c("name", "type", "in_this_network", "ffl_count_used", "T_topology", "paper_score", "score", "paper_rank_overall",
                   "rank_overall", "rank_change_overall", "paper_rank_class", "rank_class", "rank_change_class", DOM[-1]))
  setorder(M, rank_overall, name)
  fwrite(M, file.path(OUT, paste0("P_full_ranking_", net, ".csv")))
  sp <- function(x) suppressWarnings(cor(x$paper_score, x$score, method = "spearman"))
  SUMM[[net]] <- data.table(network = net, spearman_all = sp(M), spearman_TF = sp(M[type == "TF"]), spearman_Gene = sp(M[type == "Gene"]),
                            spearman_miRNA = sp(M[type == "miRNA"]), median_abs_change_overall = median(abs(M$rank_change_overall)),
                            p90_abs_change_overall = unname(quantile(abs(M$rank_change_overall), 0.9)),
                            median_abs_change_class = median(abs(M$rank_change_class)), p90_abs_change_class = unname(quantile(abs(M$rank_change_class), 0.9)),
                            n_move_50_or_more = sum(abs(M$rank_change_overall) >= 50))
  MOV[[net]] <- M[abs(rank_change_overall) >= 50, .(network = net, name, type, paper_rank_overall, rank_overall, rank_change_overall,
                                                    direction = ifelse(rank_change_overall > 0, "up", "down"))][order(-abs(rank_change_overall))]
  top1 <- fread(file.path(OUT, "runs", net, "top10_per_class.csv"))[, .(type, name)]      # written by the paper's code
  for (ty in c("TF", "Gene", "miRNA")) {
    a0 <- top0[type == ty, name]; a1 <- top1[type == ty, name]
    T4[[paste(net, ty)]] <- data.table(network = net, set = paste0("top 10 ", ty), n_kept = length(intersect(a0, a1)),
                                       entering = paste(setdiff(a1, a0), collapse = "; "), leaving = paste(setdiff(a0, a1), collapse = "; "))
  }
  for (k in c(10L, 20L)) {
    a0 <- R0[order(rank_all, name)][seq_len(k), name]; a1 <- R1[order(rank_all, name)][seq_len(k), name]
    T4[[paste(net, k)]] <- data.table(network = net, set = paste0("overall top ", k, " (extra)"), n_kept = length(intersect(a0, a1)),
                                      entering = paste(setdiff(a1, a0), collapse = "; "), leaving = paste(setdiff(a0, a1), collapse = "; "))
  }
  say(net, ": Spearman ", sprintf("%.4f", SUMM[[net]]$spearman_all), "; thirty kept ", sum(rbindlist(T4)[network == net & startsWith(set, "top 10")]$n_kept), " of 30")
}
fwrite(rbindlist(SUMM), file.path(OUT, "P_summary.csv"))
fwrite(rbindlist(MOV), file.path(OUT, "P_movers_50.csv"))
fwrite(rbindlist(T4), file.path(OUT, "P_table4_comparison.csv"))
writeLines(chk, file.path(OUT, "P_checks.txt")); for (x in chk[-1]) say(x)
say("label: ", LABEL)
close(LOG)
