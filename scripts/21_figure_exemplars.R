# Redraw the three published higher-order circuits exactly as deposited, but DIRECTED and
# SIGNED, with the figure/legend swap
# corrected.
suppressMessages({library(data.table); library(igraph); library(ggplot2); library(ggraph)
                  library(tidygraph); library(patchwork); library(ggrepel)})
REV  <- "/path/to/revision"
REPO <- "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR"
FIG  <- file.path(REV,"figures")
nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
ntype <- setNames(nodes$type, nodes$name)
canon <- fromJSON <- NULL
mm <- jsonlite::fromJSON(file.path(REV,"data/name_map.json"))$merge_map

read_sif <- function(pfx) {
  x <- fread(file.path(REPO,"SIF_files",paste0(pfx,".sif")), header=FALSE, sep="\t")
  d <- data.table(source=unlist(mm[x$V1]), target=unlist(mm[x$V3]))
  d <- unique(d[!is.na(source) & !is.na(target) & source!=target])
  d[, src_type := ntype[source]][, tgt_type := ntype[target]]
  d[, edge_type := fifelse(src_type=="miRNA" & tgt_type=="miRNA", "miRNA-miRNA co-cluster",
                    fifelse(src_type=="miRNA", "miRNA ⊣ target",
                     fifelse(tgt_type=="miRNA", "TF → miRNA",
                      fifelse(src_type=="Gene",  "gene → gene", "TF → target"))))]
  d[, effect := fifelse(src_type=="miRNA" & tgt_type!="miRNA", "repress", "activate")]
  d[]
}

EPAL <- c(`TF → target`="#EA580C", `TF → miRNA`="#F59E0B",
          `miRNA ⊣ target`="#2563EB", `miRNA-miRNA co-cluster`="#7C3AED",
          `gene → gene`="#059669")
PAL  <- c(TF="#C2410C", miRNA="#1D4ED8", Gene="#047857")

panel <- function(pfx, title, sub) {
  d  <- read_sif(pfx)
  vs <- data.table(name=unique(c(d$source,d$target)))
  vs[, type := ntype[name]]
  vs[, short := sub("^hsa-", "", name)]
  g  <- as_tbl_graph(graph_from_data_frame(d[,.(source,target,edge_type,effect)],
                                           directed=TRUE, vertices=vs))
  set.seed(7)
  ggraph(g, layout="sugiyama") +
    geom_edge_fan(aes(colour=edge_type, linetype=effect),
                  arrow=arrow(length=unit(2.1,"mm"), type="closed"),
                  end_cap=circle(3.6,"mm"), start_cap=circle(3.6,"mm"),
                  width=0.45, alpha=0.9) +
    geom_node_point(aes(fill=type, shape=type), size=4.6, colour="grey15", stroke=0.3) +
    geom_node_text(aes(label=short), size=2.5, repel=TRUE, max.overlaps=Inf,
                   point.padding=unit(3,"pt"), min.segment.length=unit(4,"pt"),
                   segment.size=0.16, segment.colour="grey55", fontface="bold") +
    scale_shape_manual(values=c(TF=21, miRNA=24, Gene=22), name="Node class") +
    scale_fill_manual(values=PAL, name="Node class") +
    scale_edge_colour_manual(values=EPAL, name="Interaction") +
    scale_edge_linetype_manual(values=c(activate="solid", repress="dashed"), name="Effect") +
    labs(title=title, subtitle=sub) +
    theme_void(base_size=9) +
    theme(plot.title=element_text(size=9.6, face="bold"),
          plot.subtitle=element_text(size=7.6, colour="grey35"),
          legend.position="right", legend.key.size=unit(8,"pt"),
          legend.text=element_text(size=7), legend.title=element_text(size=7.6),
          plot.margin=margin(6,6,6,6))
}

p1 <- panel("4-TF", "a  4-node circuit  (published as Figure 5)",
            "11 nodes / 18 edges — 4 TFs, 5 miRNAs, 2 genes")
p2 <- panel("5-TF", "b  5-node circuit  (not shown in the published manuscript)",
            "14 nodes / 41 edges — 3 TFs, 9 miRNAs, 2 genes")
p3 <- panel("6-TF", "c  6-node circuit  (published as Figure 4)",
            "10 nodes / 30 edges — 3 TFs, 5 miRNAs, 2 genes")

out <- (p1 / p2 / p3) + plot_annotation(
  title = "The three deposited higher-order circuits, redrawn as directed signed graphs",
  caption = paste0(
    "Every transcription factor acts on COL1A1 while every miRNA acts on COL3A1, so no target is co-regulated by both a TF and a miRNA:\n",
    "none of the three circuits contains a classical miRNA-TF-gene feed-forward loop core. They connect only through COL1A1 → COL3A1, an edge for which\n",
    "TRRUST holds no record in either direction and STRING supports only an undirected association. Legends of the published Figures 4\n",
    "and 5 were transposed: the legend labelled 5-node describes panel c, and the legend labelled 6-node describes panel a."),
  theme = theme(plot.title=element_text(size=11.5, face="bold"),
                plot.caption=element_text(size=7.4, colour="grey30", hjust=0)))
ggsave(file.path(FIG,"Figure3_exemplars_directed.png"), out, width=8.6, height=11.2, dpi=400, bg="white")
ggsave(file.path(FIG,"Figure3_exemplars_directed.pdf"), out, width=8.6, height=11.2, bg="white")
cat("wrote Figure3_exemplars_directed\n")
