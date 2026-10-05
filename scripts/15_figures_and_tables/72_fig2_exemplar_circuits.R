# Fig. 2 of the research-article version (1 Oct 2026): the exemplar four-, five- and six-node circuits as directed signed graphs
# (formerly Fig. S2; same data, palette and layout as 21_figure_exemplars.R, redrawn at print width with 8 pt labels, italic
# gene symbols, one shared key). Circuits are read from the SIF files of the first submission; node classes from
# data/canonical_nodes.tsv. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(igraph); library(ggplot2); library(ggraph); library(tidygraph); library(patchwork); library(ggrepel)})
REV <- Sys.getenv("FFL_REV", unset = "/path/to/revision")
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
SIFDIR <- Sys.getenv("FFL_SIF", unset = file.path(REV,"analyses/original_submission_code/SIF_files"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
TXT <- 8 / .pt
nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
ntype <- setNames(nodes$type, nodes$name)
mm <- jsonlite::fromJSON(file.path(REV,"data/name_map.json"))$merge_map

read_sif <- function(pfx) {
  x <- fread(file.path(SIFDIR, paste0(pfx,".sif")), header=FALSE, sep="\t")
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
ECOL <- c(`TF → target`=EPAL[["TF_target"]], `TF → miRNA`=EPAL[["TF_miRNA"]], `miRNA ⊣ target`=EPAL[["miRNA_target"]],
          `miRNA-miRNA co-cluster`=EPAL[["miRNA_miRNA"]], `gene → gene`=EPAL[["gene_gene"]])

panel <- function(pfx, title, sub) {
  d  <- read_sif(pfx)
  vs <- data.table(name=unique(c(d$source,d$target)))
  vs[, type := ntype[name]]
  vs[, short := sub("^hsa-", "", name)]
  g  <- as_tbl_graph(graph_from_data_frame(d[,.(source,target,edge_type,effect)], directed=TRUE, vertices=vs))
  set.seed(7)
  ggraph(g, layout="sugiyama") +
    geom_edge_fan(aes(colour=edge_type, linetype=effect), arrow=arrow(length=unit(2.0,"mm"), type="closed"),
                  end_cap=circle(3.4,"mm"), start_cap=circle(3.4,"mm"), width=0.4, alpha=0.9) +
    geom_node_point(aes(fill=type, shape=type), size=4.4, colour="grey15", stroke=0.3) +
    geom_node_text(aes(label=short, fontface=ifelse(type=="miRNA","plain","italic")), size=TXT, repel=TRUE, max.overlaps=Inf,
                   bg.colour="white", bg.r=0.16, point.padding=unit(5,"pt"), box.padding=0.35, force=2, max.iter=30000, min.segment.length=unit(4,"pt"),
                   segment.size=0.16, segment.colour="grey55", seed=7) +
    scale_shape_manual(values=c(TF=21, miRNA=24, Gene=22), limits=c("TF","miRNA","Gene"), name="Node class") +
    scale_fill_manual(values=PAL, limits=c("TF","miRNA","Gene"), name="Node class") +
    scale_edge_colour_manual(values=ECOL, limits=names(ECOL), drop=FALSE, name="Interaction (dashed lines: repression by a miRNA)") +
    scale_edge_linetype_manual(values=c(activate="solid", repress="dashed"), guide="none") +
    guides(edge_colour=guide_legend(order=2, nrow=3, byrow=FALSE, title.position="top", override.aes=list(linewidth=0.8)),
           shape=guide_legend(order=1, title.position="top"), fill=guide_legend(order=1, title.position="top")) +
    labs(title=title, subtitle=sub) +
    theme_void(base_size=9.5) +
    theme(plot.title=element_text(size=9.5, face="bold", colour=INK), plot.subtitle=element_text(size=8, colour=MUTED),
          legend.position="bottom", legend.box="horizontal", legend.key.size=unit(9,"pt"), legend.spacing.x=unit(6,"pt"),
          legend.text=element_text(size=8, colour=INK), legend.title=element_text(size=8.5, colour=INK),
          plot.tag=element_text(size=11.5, face="bold", colour=INK), plot.margin=margin(4,6,4,6))
}
p1 <- panel("4-TF", "four-node circuit", "11 nodes, 18 edges: 4 TFs, 5 miRNAs, 2 genes") + theme(legend.position="none")
p2 <- panel("5-TF", "five-node circuit", "14 nodes, 41 edges: 3 TFs, 9 miRNAs, 2 genes") + theme(legend.position="none")
p3 <- panel("6-TF", "six-node circuit", "10 nodes, 30 edges: 3 TFs, 5 miRNAs, 2 genes")
out <- (p1 / p2 / p3) + plot_layout(heights=c(1, 1.35, 1.05), guides="collect") +
  plot_annotation(tag_levels="a", theme=theme(legend.position="bottom"))
ggsave(file.path(FIG,"Fig2_exemplar_circuits.png"), out, width=6.85, height=9.0, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"Fig2_exemplar_circuits.pdf"), out, width=6.85, height=9.0, bg="white", device=cairo_pdf)
cat("wrote Fig2_exemplar_circuits\n")
