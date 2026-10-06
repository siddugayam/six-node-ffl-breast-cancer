## ---------------------------------------------------------------------------
## v7 -- characterisation of the MCODE result: agreement with the paper's
## prioritisation, placement of the TF arm, edge-evidence provenance of the top
## module, and a degree-preserving null for the top module's MCODE score.
## ---------------------------------------------------------------------------
source("scripts/14_module_detection/mcode.R")
suppressPackageStartupMessages(library(igraph))
set.seed(20260912)
OUT <- "results/v7"; p <- function(x) file.path(OUT, x)
load(p("mcode_workspace.RData"))
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors = FALSE)
edges <- read.delim("data/canonical_edges.tsv", stringsAsFactors = FALSE)
nm <- V(g)$name; ty <- V(g)$type; N <- vcount(g)
pri30 <- refsets$prioritised_30; ms_ex <- refsets$MS_exemplar_module
mir29 <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
hyper <- function(hit, K, n, N) if (hit == 0) 1 else phyper(hit-1, K, N-K, n, lower.tail = FALSE)

## ---- A. agreement with the prioritised 30 --------------------------------
memb <- strsplit(res_def$members, ";")
rows <- list()
for (k in seq_along(memb)) {
  u <- unique(unlist(memb[seq_len(k)]))
  rows[[k]] <- data.frame(top_k_clusters = k, n_nodes = length(u),
    n_prioritised = sum(pri30 %in% u), expected = length(pri30)*length(u)/N,
    fold = (sum(pri30 %in% u)/length(u))/(length(pri30)/N),
    p_hyper = hyper(sum(pri30 %in% u), length(pri30), length(u), N),
    jaccard = length(intersect(pri30,u))/length(union(pri30,u)),
    hits = paste(intersect(pri30, u), collapse=";"), stringsAsFactors = FALSE)
}
agree <- do.call(rbind, rows)
write.csv(agree, p("mcode_10_agreement_prioritised30.csv"), row.names = FALSE)
cat("== AGREEMENT WITH THE 30 PRIORITISED NODES (cumulative over ranked clusters) ==\n")
print(agree[, 1:7])

## per-cluster, and the reverse direction: where do the 30 sit?
clust_of <- setNames(rep(NA_integer_, N), nm)
for (k in seq_along(memb)) clust_of[memb[[k]]] <- k
w <- data.frame(node = pri30, type = nodes$type[match(pri30, nodes$name)],
                mcode_cluster = unname(clust_of[pri30]))
w$in_any_cluster <- !is.na(w$mcode_cluster)
write.csv(w, p("mcode_11_prioritised30_placement.csv"), row.names = FALSE)
cat("\nprioritised 30 falling in ANY MCODE cluster:", sum(w$in_any_cluster), "/30\n")
cat("clusters they land in:", paste(sort(unique(na.omit(w$mcode_cluster))), collapse=","), "\n")
cat("nodes covered by the 13 clusters:", sum(!is.na(clust_of)), "of", N, "\n")

## ---- B. placement of the named axis members ------------------------------
named <- c("COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1","EZH2",
           "hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-101",
           "hsa-let-7b","hsa-let-7e")
plc <- data.frame(node = named, type = nodes$type[match(named, nodes$name)],
                  degree = degree(g)[match(named, nm)],
                  mcode_node_score = info$score[match(named, nm)],
                  mcode_cluster = unname(clust_of[named]),
                  in_MS_exemplar = named %in% ms_ex)
write.csv(plc, p("mcode_12_named_axis_placement.csv"), row.names = FALSE)
cat("\n== PLACEMENT OF THE NAMED AXIS MEMBERS ==\n"); print(plc)

## ---- C. what the top module actually is ----------------------------------
c1 <- memb[[1]]
sub1 <- induced_subgraph(g, match(c1, nm))
el <- as_edgelist(sub1)
ev <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors = FALSE)
key <- function(a,b) paste(pmin(a,b), pmax(a,b))
evk <- setNames(ev$tier, key(ev$source, ev$target))
de <- edges[key(edges$source, edges$target) %in% key(el[,1], el[,2]), ]
de$tier <- evk[key(de$source, de$target)]
de$tier[is.na(de$tier)] <- "not_a_miRNA_target_edge"
write.csv(de, p("mcode_13_top_cluster_edges.csv"), row.names = FALSE)
cat("\n== TOP MODULE (C1) EDGE PROVENANCE ==\n")
print(table(de$edge_type, de$tier))
cat("C1 undirected edges:", ecount(sub1), " density:", edge_density(sub1), "\n")
cat("C1 members:", paste(c1, collapse=", "), "\n")

## ---- D. degree-preserving null for the top module's score ----------------
B <- 200
nullres <- do.call(rbind, lapply(seq_len(B), function(b) {
  gr <- rewire(g, keeping_degseq(loops = FALSE, niter = 20 * ecount(g)))
  ir <- mcode_score_graph(gr); cr <- mcode_find_clusters(gr, ir)
  if (length(cr) == 0) return(data.frame(b=b, n_clusters=0, top_score=NA, top_size=NA,
                                         max_size=NA, n_score_ge_6.667=0))
  sc <- sapply(cr, function(x) x$score); sz <- sapply(cr, function(x) x$size)
  data.frame(b = b, n_clusters = length(cr), top_score = max(sc),
             top_size = sz[which.max(sc)], max_size = max(sz),
             n_score_ge_6.667 = sum(sc >= res_def$score[1] - 1e-9))
}))
write.csv(nullres, p("mcode_14_degree_preserving_null.csv"), row.names = FALSE)
obs <- res_def$score[1]
pnull <- (sum(nullres$top_score >= obs - 1e-9) + 1) / (B + 1)
cat("\n== DEGREE-PRESERVING NULL (", B, "rewirings, degree sequence preserved ) ==\n")
cat(sprintf("observed top MCODE score %.4f; null top score mean %.3f sd %.3f max %.3f; p = %.4f\n",
            obs, mean(nullres$top_score), sd(nullres$top_score), max(nullres$top_score), pnull))
cat(sprintf("observed n clusters %d; null mean %.1f (sd %.1f)\n",
            nrow(res_def), mean(nullres$n_clusters), sd(nullres$n_clusters)))

## ---- E. rank of the collagen/miR-29 module across the cutoff sweep -------
ax <- sens[sens$collagens_and_mir29_together, ]
axtab <- do.call(rbind, lapply(split(ax, ax$run), function(d) {
  d <- d[order(d$rank), ][1, ]
  data.frame(run = d$run, axis_cluster_rank = d$rank, size = d$size,
             score = d$score, n_miR29 = d$n_miR29, n_MS_exemplar = d$n_MS_exemplar,
             members = d$members, stringsAsFactors = FALSE) }))
axtab <- axtab[order(axtab$run), ]
write.csv(axtab, p("mcode_15_axis_cluster_across_cutoffs.csv"), row.names = FALSE)
cat("\n== COLLAGEN + miR-29 CLUSTER ACROSS NODE-SCORE CUTOFFS ==\n")
print(axtab[, c("run","axis_cluster_rank","size","score","n_miR29","n_MS_exemplar")])
cat("\nmembers at each cutoff:\n"); for(i in 1:nrow(axtab)) cat(" ", axtab$run[i], ":", axtab$members[i], "\n")
cat("\nDONE\n")
