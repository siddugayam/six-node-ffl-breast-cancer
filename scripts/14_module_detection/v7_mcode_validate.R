## Validation of the R port against the reference Java implementation's own
## unit tests (MCODEAlgorithmTest.java) and against ProNet::mcode.
source("scripts/14_module_detection/mcode.R")
suppressPackageStartupMessages(library(igraph))
out <- list()

## --- Java test 1: testCompleteGraphWithDefaultParameters -------------------
## expects 1 cluster, cluster score 16, every node score 15.0
g <- make_full_graph(16)
info <- mcode_score_graph(g)
cl <- mcode_find_clusters(g, info)
t1 <- c(n_clusters = length(cl), score = cl[[1]]$score, size = cl[[1]]$size,
        edges = cl[[1]]$edges, node_score = unique(info$score))
cat("TEST 1 complete graph K16, defaults:\n"); print(t1)
stopifnot(length(cl) == 1, abs(cl[[1]]$score - 16) < 1e-12,
          cl[[1]]$size == 16, cl[[1]]$edges == 120,
          all(abs(info$score - 15) < 1e-12))
out$test1 <- "PASS (1 cluster, score 16, 16 nodes, 120 edges, node score 15)"

## --- Java test 2: testCompleteGraphIncludingLoops --------------------------
## expects 1 cluster, score 14.118 (+/- 0.0009)
info2 <- mcode_score_graph(g, include_loops = TRUE)
cl2 <- mcode_find_clusters(g, info2, include_loops = TRUE)
cat("TEST 2 complete graph K16, includeLoops=TRUE: score =", cl2[[1]]$score, "\n")
stopifnot(length(cl2) == 1, abs(cl2[[1]]$score - 14.118) < 0.0009)
out$test2 <- sprintf("PASS (score %.4f, Java expects 14.118 +/- 0.0009)", cl2[[1]]$score)

## --- hand-checked small case A: two K4s joined by a bridge ---------------
## Hand trace: every one of the 8 vertices scores 3 (closed neighbourhood of a
## bridge endpoint has highest 3-core = its own K4, density 1).  With a uniform
## score landscape the seed threshold 3*(1-0.2)=2.4 admits every vertex, so
## MCODE walks straight across the bridge and returns ONE 8-node cluster,
## density 13/28, score 3.714.  This is real MCODE behaviour, not a bug, and it
## is the failure mode that matters for a network with many tied scores.
e <- rbind(t(combn(1:4,2)), t(combn(5:8,2)), c(4,5))
g3 <- graph_from_edgelist(e, directed = FALSE)
i3 <- mcode_score_graph(g3); cl3 <- mcode_find_clusters(g3, i3)
cat("TEST 3 two bridged K4s: n =", length(cl3),
    " size =", cl3[[1]]$size, " score =", round(cl3[[1]]$score,4),
    " node scores =", paste(unique(i3$score), collapse=","), "\n")
stopifnot(length(cl3) == 1, cl3[[1]]$size == 8,
          abs(cl3[[1]]$score - 13/28*8) < 1e-9, all(i3$score == 3))
out$test3 <- "PASS (bridged K4s merge into one 8-node cluster, score 3.7143 - hand traced)"

## --- hand-checked small case B: K5 with a pendant path -------------------
## Hand trace: nodes 1-5 score 4; the pendant path nodes score 2/3 and 0.
## Seed 1, threshold 3.2 -> cluster {1..5}, score 5.  The {6,7} attempt has no
## 2-core and is filtered.  Exactly one cluster.
e <- rbind(t(combn(1:5,2)), c(1,6), c(6,7), c(7,8))
g4 <- graph_from_edgelist(e, directed = FALSE)
i4 <- mcode_score_graph(g4); cl4 <- mcode_find_clusters(g4, i4)
cat("TEST 4 K5 + pendant path: n =", length(cl4), " size =", cl4[[1]]$size,
    " score =", cl4[[1]]$score, " scores =", paste(round(i4$score,3), collapse=","), "\n")
stopifnot(length(cl4) == 1, cl4[[1]]$size == 5, abs(cl4[[1]]$score - 5) < 1e-12,
          setequal(cl4[[1]]$nodes, 1:5),
          all(abs(i4$score[1:5] - 4) < 1e-12), abs(i4$score[6] - 2/3) < 1e-12,
          i4$score[8] == 0)
out$test4 <- "PASS (K5 + pendant path -> one 5-node cluster, score 5 - hand traced)"

## --- haircut / k-core helper ---------------------------------------------
e <- rbind(t(combn(1:4,2)), c(1,5))
g5 <- graph_from_edgelist(e, directed = FALSE)
stopifnot(setequal(mcode_kcore_vids(g5, 2L), 1:4), is.null(mcode_kcore(g5, 4L)))
out$test5 <- "PASS (2-core haircut drops the singly-connected vertex)"

## --- ProNet cross-check on a random graph ----------------------------------
## ProNet implements the paper pseudocode (threshold re-referenced at each
## recursion step, no degree cutoff / k-core / max depth), so exact agreement
## is not expected; this checks the port is in the same family.
set.seed(1)
gr <- sample_gnp(120, 0.06)
V(gr)$name <- paste0("n", seq_len(vcount(gr)))
mine <- mcode_find_clusters(gr)
ok <- requireNamespace("ProNet", quietly = TRUE)
if (ok) {
  pn <- ProNet::mcode(gr, vwp = 0.2, haircut = TRUE, fluff = FALSE, loops = FALSE)
  jac <- function(a, b) length(intersect(a,b)) / length(union(a,b))
  best <- sapply(mine, function(m) max(c(0, sapply(pn$COMPLEX, function(p) jac(m$nodes, p)))))
  cat("TEST 6 vs ProNet: mine n =", length(mine), " ProNet n =", length(pn$COMPLEX),
      " best Jaccard per cluster:", paste(round(best,2), collapse=","), "\n")
  out$test6 <- sprintf("ProNet 1.0.0 run: %d vs %d clusters, median best Jaccard %.2f",
                       length(mine), length(pn$COMPLEX), median(best))
} else out$test6 <- "ProNet not available"

writeLines(paste(names(out), unlist(out), sep = ": "),
           "results/v7/mcode_00_validation.txt")
cat("\n--- validation written ---\n")
