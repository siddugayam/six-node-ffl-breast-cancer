#!/usr/bin/env Rscript
# v7 / 09 : sensitivity to network composition. 4,819 of the 6,859 directed edges
#           are miRNA->target, so the community structure of the whole graph is
#           dominated by shared-miRNA-target bipartite structure. If the reactive
#           stroma programme is a protein-level module, it should be easier to see
#           in the TF/gene-only subgraph. Also: TF->target layer only.
suppressMessages({library(igraph); library(data.table)})
setwd("/path/to/revision")
set.seed(20260912)

nodes <- fread("data/canonical_nodes.tsv"); ntype <- setNames(nodes$type, nodes$name)
g <- readRDS("results/v7/cm2_graph_undirected.rds")
ECM <- c("COL1A1","COL3A1","COL7A1","COL18A1","FN1","SPARC","POSTN","LOX","MMP2","MMP14",
         "THBS1","THBS2","TGFB2","TGFBI","SERPINE1","ACTA2","TAGLN","VIM","CDH11","FBN1",
         "CTGF","PDGFRB","ITGA5","ITGB1")
EMTall <- unique(fread("results/v7/cm2_genesets_msigdb.csv")[
  gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", gene_symbol])

report <- function(gg, tag) {
  gg <- induced_subgraph(gg, which(degree(gg) > 0))
  V(gg)$name -> vn
  if (!("COL1A1" %in% vn)) { cat(tag, ": COL1A1 absent\n"); return(NULL) }
  prot <- vn[ntype[vn] %in% c("Gene", "TF")]
  ecm <- intersect(ECM, vn); emt <- intersect(EMTall, vn)
  res <- rbindlist(lapply(list(
      list("Louvain",     function() membership(cluster_louvain(gg))),
      list("FastGreedy",  function() membership(cluster_fast_greedy(gg))),
      list("Walktrap",    function() membership(cluster_walktrap(gg, steps = 4))),
      list("Infomap",     function() membership(cluster_infomap(gg, nb.trials = 10))),
      list("LabelProp",   function() membership(cluster_label_prop(gg))),
      list("Leiden",      function() membership(cluster_leiden(gg, objective_function = "modularity",
                                                               resolution = 1, n_iterations = 10)))),
    function(a) {
      m <- a[[2]](); names(m) <- vn
      cl <- m["COL1A1"]; mem <- names(m)[m == cl]; pm <- intersect(mem, prot)
      obs <- max(table(m[ecm]))
      null <- replicate(3000, max(table(m[sample(prot, length(ecm))])))
      data.table(network = tag, algorithm = a[[1]], n = length(vn), k = max(m),
                 modularity = round(modularity(gg, m), 4),
                 COL1A1_module = length(mem),
                 COL3A1_same = unname(m["COL3A1"] == cl),
                 miR29a_same = if ("hsa-miR-29a" %in% vn) unname(m["hsa-miR-29a"] == cl) else NA,
                 ECM_in_module = length(intersect(mem, ecm)),
                 ECM_max_any_module = obs, ECM_exp = round(mean(null), 1),
                 ECM_p = (sum(null >= obs) + 1) / 3001,
                 EMT_in_module = length(intersect(mem, emt)),
                 EMT_exp = round(length(pm) * length(emt) / length(prot), 1))
    }))
  res
}

e <- fread("data/canonical_edges.tsv")
prot_nodes <- nodes[type %in% c("Gene", "TF"), name]
g_prot <- induced_subgraph(g, which(V(g)$name %in% prot_nodes))
g_tft  <- simplify(as_undirected(graph_from_data_frame(
            e[edge_type %in% c("TF_target", "gene_gene"), .(source, target)], directed = TRUE),
          mode = "collapse"))
g_mir  <- simplify(as_undirected(graph_from_data_frame(
            e[edge_type %in% c("miRNA_target", "miRNA_miRNA"), .(source, target)], directed = TRUE),
          mode = "collapse"))

out <- rbindlist(list(report(g, "full (587 nodes)"),
                      report(g_prot, "TF+gene only"),
                      report(g_tft,  "TF->target layer only"),
                      report(g_mir,  "miRNA->target layer only")), fill = TRUE)
fwrite(out, "results/v7/cm2_subnetwork_sensitivity.csv")
print(out, nrows = 40)
