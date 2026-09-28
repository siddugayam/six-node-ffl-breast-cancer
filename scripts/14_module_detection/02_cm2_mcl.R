#!/usr/bin/env Rscript
# v7 / 01 : Markov clustering (van Dongen 2000), the algorithm clusterMaker2
#           implements as its "MCL Cluster" command.
#
# Implemented here directly so that the inflation sweep (1.8 - 4.0) can be run
# reproducibly from the command line. Cross-checked against the CRAN `MCL`
# package (v1.0) at inflation 2.0 in 02_cm2_battery.R, and against the
# clusterMaker2 2.3.4 Java implementation driven through CyREST in
# 03_clustermaker2_cyrest.R.
#
# Parameters follow clusterMaker2 / mcl defaults unless stated:
#   expansion power  = 2          (fixed in both mcl and clusterMaker2)
#   self-loops       = max edge weight per node (= 1 here; clusterMaker2 adds
#                      loops the same way so that singleton nodes converge)
#   pruning          = 1e-6       (clusterMaker2 "edge weight cutoff" 1e-15;
#                      1e-6 is mcl's default -P equivalent and is tested for
#                      sensitivity)
#   convergence      = chaos < 1e-6 or 200 iterations

mcl_vandongen <- function(A, inflation, expansion = 2, add_loops = TRUE,
                          prune = 1e-6, tol = 1e-6, max_iter = 200) {
  A <- as.matrix(A); storage.mode(A) <- "double"
  n <- nrow(A)
  if (add_loops) diag(A) <- pmax(apply(A, 1, max), 1)
  ncol_norm <- function(M) {
    cs <- colSums(M); cs[cs <= 0] <- 1
    sweep(M, 2, cs, "/")
  }
  M <- ncol_norm(A)
  it <- 0L; chaos <- Inf
  repeat {
    it <- it + 1L
    E <- M
    if (expansion > 1) for (k in seq_len(expansion - 1)) E <- E %*% M
    E <- E ^ inflation
    E[E < prune] <- 0
    E <- ncol_norm(E)
    chaos <- max(apply(E, 2, function(x) max(x) - sum(x^2)))
    M <- E
    if (chaos < tol || it >= max_iter) break
  }
  B <- (M > prune)
  B <- B | t(B)
  g <- igraph::graph_from_adjacency_matrix(B, mode = "undirected", diag = FALSE)
  memb <- igraph::components(g)$membership
  names(memb) <- rownames(A)
  list(membership = memb, iterations = it, chaos = chaos,
       converged = chaos < tol)
}
