## F2: does the 70-hub list of the TCGA-BRCA Cox screen change without the legacy miRNA-miRNA edges?
## The hub list is built by scripts/10_survival_and_clinical/cox_models/10_survival_cox_hubs.R (lines 27-70) on a RECONSTRUCTED network:
## miRNA_target from data/canonical_edges.tsv, TF_target, TF_miRNA, gene_gene and miRNA_miRNA from the
## layer files.  The 30 legacy GeneMANIA edges are canonical_edges.tsv edge_type miRNA_miRNA, which that
## script never reads.  Here: (1) the stored hub list is rebuilt with the same code and compared with
## results/network_topology_hubs.csv; (2) the legacy pairs present in the reconstructed network are listed;
## (3) sensitivity: those pairs are removed and the hub list rebuilt.  Read-only on the project.
suppressPackageStartupMessages({ library(data.table); library(igraph) })
BASE <- "/path/to/revision"
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv")); cedge <- fread(file.path(BASE, "data/canonical_edges.tsv"))
mir_t <- cedge[edge_type == "miRNA_target", .(source, target, edge_type = "miRNA_target")]
tf_t  <- fread(file.path(BASE, "data/layer_TF_target.tsv"))[, .(source, target, edge_type = "TF_target")]
tf_m  <- fread(file.path(BASE, "data/layer_TF_miRNA.tsv"))[, .(source, target, edge_type = "TF_miRNA")]
gg    <- fread(file.path(BASE, "data/layer_gene_gene.tsv"))[, .(source, target, edge_type = "gene_gene")]
mm    <- fread(file.path(BASE, "data/layer_miRNA_miRNA.tsv"))[, .(source = miRNA_1, target = miRNA_2, edge_type = "miRNA_miRNA")]
recon <- unique(rbindlist(list(mir_t, tf_t, tf_m, gg, mm))); recon <- recon[source %in% nodes$name & target %in% nodes$name]
hubs <- function(re) {
  g <- graph_from_data_frame(re[, .(source, target)], directed = TRUE, vertices = data.frame(name = nodes$name))
  topo <- data.table(name = V(g)$name, degree_tot = degree(g, mode = "all"), betweenness = betweenness(g, directed = TRUE))
  hub_deg <- topo[order(-degree_tot)][1:50, name]; hub_btw <- topo[order(-betweenness)][1:50, name]
  named <- c("NFKB1", "RELA", "SP1", "ETS1", "COL1A1", "COL3A1", "VEGFA", "CCND2", "MYC", "E2F1", "TP53", "hsa-miR-130a",
             "hsa-miR-124", "hsa-miR-101", "hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-let-7b", "hsa-let-7e")
  list(topo = topo, hubs = sort(union(union(hub_deg, hub_btw), named)))
}
st <- fread(file.path(BASE, "results/network_topology_hubs.csv"))
h0 <- hubs(recon)
cat("rebuilt hub list =", length(h0$hubs), "; identical to the stored is_hub list:", setequal(h0$hubs, st[is_hub == TRUE, name]), "\n")
m <- merge(h0$topo, st[, .(name, d = degree_tot, b = betweenness)], by = "name")
cat("degree and betweenness identical to stored:", all(m$degree_tot == m$d), isTRUE(all.equal(m$betweenness, m$b)), "\n")
leg <- cedge[edge_type == "miRNA_miRNA", .(a = pmin(source, target), b = pmax(source, target))]
rk <- recon[, .(a = pmin(source, target), b = pmax(source, target), edge_type)]
ov <- merge(rk, leg, by = c("a", "b"))
cat("legacy GeneMANIA pairs present in the reconstructed network:", nrow(ov), "\n"); print(ov)
cat("edge types of the reconstructed network:\n"); print(recon[, .N, by = edge_type])
re2 <- recon[!paste(pmin(source, target), pmax(source, target)) %in% paste(ov$a, ov$b)]
h1 <- hubs(re2)
cat("sensitivity, those pairs removed: hub list", length(h1$hubs), "identical:", setequal(h1$hubs, h0$hubs),
    "; added:", paste(setdiff(h1$hubs, h0$hubs), collapse = ","), "; lost:", paste(setdiff(h0$hubs, h1$hubs), collapse = ","), "\n")
fwrite(data.table(hub = h0$hubs), file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "f2_hub_list_70.csv"))
