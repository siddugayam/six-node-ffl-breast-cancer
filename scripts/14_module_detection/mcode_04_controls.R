## ---------------------------------------------------------------------------
## v7 -- controls for the MCODE result.
## The manuscript's canonical network contains a 30-edge miRNA-miRNA layer that
## exists ONLY among the nine miRNAs of the manuscript's own exemplar module,
## and which the project's own audit (results/00_DATA_AUDIT.md, item E3) shows
## is GeneMANIA co-expression for 28 of the 30 pairs, with no co-transcriptional
## basis.  Half the edges of MCODE's top-scoring cluster are those edges.  These
## controls ask whether MCODE still finds the module without them.
## ---------------------------------------------------------------------------
source("scripts/14_module_detection/mcode.R")
suppressPackageStartupMessages(library(igraph))
set.seed(20260912)
OUT <- "results/v7"; p <- function(x) file.path(OUT, x)
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors = FALSE)
edges <- read.delim("data/canonical_edges.tsv", stringsAsFactors = FALSE)
L     <- read.delim("data/layer_miRNA_miRNA.tsv", stringsAsFactors = FALSE)
ev    <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors = FALSE)
mm    <- read.csv("results/ffl_module_membership.csv", stringsAsFactors = FALSE)
pri30 <- read.csv("results/v5/tables/Table4_prioritised_30.csv", stringsAsFactors = FALSE)$name
ms_ex <- mm$node[mm$in_MS_exemplar_module]
mir29 <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"); collagen <- c("COL1A1","COL3A1")
k <- function(a,b) paste(pmin(a,b), pmax(a,b))

build <- function(el) {
  gd <- graph_from_data_frame(el[, c("source","target")], directed = TRUE,
                              vertices = nodes[, c("name","type")])
  g <- simplify(as_undirected(gd, mode = "collapse"), remove.multiple = TRUE, remove.loops = TRUE)
  V(g)$type <- nodes$type[match(V(g)$name, nodes$name)]; g }

hyper <- function(hit,K,n,N) if (hit==0) 1 else phyper(hit-1,K,N-K,n,lower.tail=FALSE)
run <- function(g, tag, note) {
  info <- mcode_score_graph(g); cl <- mcode_find_clusters(g, info)
  nm <- V(g)$name; ty <- V(g)$type; N <- vcount(g)
  if (length(cl) == 0) return(data.frame(network = tag, note = note, rank = NA,
      size = 0, score = NA, stringsAsFactors = FALSE))
  do.call(rbind, lapply(cl, function(c) {
    v <- nm[c$nodes]
    data.frame(network = tag, note = note, rank = c$rank, seed = nm[c$seed],
      size = c$size, edges = c$edges, score = c$score,
      n_TF = sum(ty[c$nodes]=="TF"), n_miRNA = sum(ty[c$nodes]=="miRNA"),
      n_Gene = sum(ty[c$nodes]=="Gene"),
      COL1A1 = "COL1A1" %in% v, COL3A1 = "COL3A1" %in% v, n_miR29 = sum(mir29 %in% v),
      collagens_and_mir29_together = all(collagen %in% v) && any(mir29 %in% v),
      n_MS_exemplar = sum(ms_ex %in% v),
      p_MS_exemplar = hyper(sum(ms_ex %in% v), length(ms_ex), c$size, N),
      n_prioritised30 = sum(pri30 %in% v),
      members = paste(v, collapse=";"), stringsAsFactors = FALSE) })) }

## --- network variants ------------------------------------------------------
is_mm <- edges$edge_type == "miRNA_miRNA"
genomic <- k(L$miRNA_1, L$miRNA_2)[L$threshold_10kb]
supported_mm <- is_mm & (k(edges$source, edges$target) %in% genomic)
strongk <- k(ev$source, ev$target)[ev$tier == "strong"]

variants <- list()
variants[["as_published"]] <- list(el = edges, note = "canonical network as used by the paper")
variants[["drop_all_miRNA_miRNA"]] <- list(el = edges[!is_mm, ],
  note = "all 30 miRNA-miRNA edges removed")
variants[["keep_only_coclustered_miRNA_miRNA"]] <- list(el = edges[!is_mm | supported_mm, ],
  note = "only the 2 genomically co-clustered miRNA-miRNA pairs kept (audit E3)")
## replace the hand-drawn layer with a consistently built one: all 10kb
## polycistronic pairs among canonical miRNA nodes
gen_el <- L[L$threshold_10kb & L$miRNA_1 %in% nodes$name & L$miRNA_2 %in% nodes$name,
            c("miRNA_1","miRNA_2")]
names(gen_el) <- c("source","target"); gen_el$edge_type <- "miRNA_miRNA"
variants[["consistent_genomic_miRNA_layer"]] <- list(
  el = rbind(edges[!is_mm, c("source","target","edge_type")], gen_el),
  note = sprintf("hand-drawn layer replaced by all %d 10-kb polycistronic pairs", nrow(gen_el)))
variants[["drop_miRNA_miRNA_and_gene_gene"]] <- list(
  el = edges[!is_mm & edges$edge_type != "gene_gene", ],
  note = "miRNA-miRNA layer and the single COL1A1-COL3A1 gene_gene edge removed")
keepstrong <- (k(edges$source, edges$target) %in% strongk) | (edges$edge_type != "miRNA_target")
variants[["strong_miRNA_target_only"]] <- list(el = edges[keepstrong, ],
  note = "miRNA-target edges restricted to the strong evidence tier")
variants[["strong_and_no_miRNA_miRNA"]] <- list(el = edges[keepstrong & !is_mm, ],
  note = "strong miRNA-target tier AND no miRNA-miRNA edges")

res <- do.call(rbind, lapply(names(variants), function(v) {
  g <- build(variants[[v]]$el)
  cat(sprintf("%-36s %d nodes %d undirected edges\n", v, vcount(g), ecount(g)))
  run(g, v, variants[[v]]$note) }))
write.csv(res, p("mcode_16_controls.csv"), row.names = FALSE)

cat("\n== CONTROLS: does the collagen + miR-29 module survive? ==\n")
smry <- do.call(rbind, lapply(split(res, factor(res$network, levels=names(variants))), function(d) {
  ax <- d[which(d$collagens_and_mir29_together), ]
  ex <- d[which.max(d$n_MS_exemplar), ]
  data.frame(network = d$network[1], n_clusters = sum(d$size > 0),
             top_size = d$size[1], top_score = round(d$score[1], 3), top_seed = d$seed[1],
             axis_found = nrow(ax) > 0,
             axis_rank = if (nrow(ax)) ax$rank[1] else NA,
             axis_size = if (nrow(ax)) ax$size[1] else NA,
             axis_score = if (nrow(ax)) round(ax$score[1],3) else NA,
             best_exemplar_overlap = ex$n_MS_exemplar,
             best_exemplar_p = signif(ex$p_MS_exemplar, 3),
             stringsAsFactors = FALSE) }))
print(smry, row.names = FALSE)
write.csv(smry, p("mcode_17_controls_summary.csv"), row.names = FALSE)

cat("\n-- top 3 clusters of each control --\n")
for (v in names(variants)) { d <- res[res$network == v, ]
  cat("\n##", v, "--", d$note[1], "\n")
  for (i in seq_len(min(3, nrow(d)))) if (d$size[i] > 0)
    cat(sprintf("  C%d size=%d score=%.3f : %s\n", d$rank[i], d$size[i], d$score[i], d$members[i])) }
cat("\nDONE\n")
