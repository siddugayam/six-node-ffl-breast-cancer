# Fig. S1 (the global regulatory network, FigureS1_global_network), with a census panel and directed diagrams
# of the exemplar modules; Figs S2 and 2 of the paper are drawn by scripts 70 and 72.
suppressMessages({library(data.table); library(igraph); library(ggplot2); library(ggraph)
                  library(tidygraph); library(patchwork); library(jsonlite); library(scales)})
REV <- "/path/to/revision"
FIG <- file.path(REV,"figures"); dir.create(FIG, showWarnings=FALSE)
lg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
edges <- fread(file.path(REV,"data/canonical_edges.tsv"))
lg("nodes",nrow(nodes),"edges",nrow(edges))

# augment with rebuilt layers, as used for the census
add <- function(f, s, t, et) {
  p <- file.path(REV,"data",f); if (!file.exists(p)) return(NULL)
  d <- fread(p); if (!all(c(s,t) %in% names(d))) return(NULL)
  data.table(source=d[[s]], target=d[[t]], edge_type=et)
}
extra <- rbindlist(list(
  add("layer_TF_target.tsv","source","target","TF_target"),
  add("layer_gene_gene.tsv","source","target","gene_gene"),
  add("layer_miRNA_miRNA.tsv","miRNA_1","miRNA_2","miRNA_miRNA")), use.names=TRUE, fill=TRUE)
extra <- extra[source %in% nodes$name & target %in% nodes$name & source!=target]
E <- unique(rbind(edges[,.(source,target,edge_type)], extra))
lg("augmented edges", nrow(E), " by type:", paste(names(table(E$edge_type)), table(E$edge_type), collapse=" "))

PAL <- c(TF="#C2410C", miRNA="#1D4ED8", Gene="#047857")
EPAL <- c(TF_target="#EA580C", TF_miRNA="#F59E0B", miRNA_target="#2563EB",
          miRNA_miRNA="#7C3AED", gene_gene="#059669")

# ---------------------------------------------------------------- Fig A: global network
g <- graph_from_data_frame(E[,.(source,target,edge_type)], directed=TRUE,
                           vertices=nodes[,.(name,type)])
g <- delete_vertices(g, V(g)[degree(g)==0])
V(g)$deg <- degree(g)
lab <- head(V(g)$name[order(-V(g)$deg)], 28)
V(g)$lab <- ifelse(V(g)$name %in% lab, V(g)$name, NA)
set.seed(42)
pA <- ggraph(as_tbl_graph(g), layout="fr") +
  geom_edge_fan(aes(colour=edge_type), width=0.16, alpha=0.22, show.legend=TRUE) +
  geom_node_point(aes(fill=type, size=deg), shape=21, colour="grey20", stroke=0.18) +
  geom_node_text(aes(label=lab), size=2.3, repel=TRUE, max.overlaps=40,
                 family="sans", segment.size=0.15, segment.colour="grey55") +
  scale_edge_colour_manual(values=EPAL, name="Interaction") +
  scale_fill_manual(values=PAL, name="Node class") +
  scale_size_continuous(range=c(0.7,6.5), name="Degree") +
  theme_void(base_size=9) +
  theme(legend.position="right", legend.key.size=unit(8,"pt"),
        plot.margin=margin(4,4,4,4)) +
  guides(fill=guide_legend(override.aes=list(size=3.4)),
         edge_colour=guide_legend(override.aes=list(width=1.4, alpha=1)))
ggsave(file.path(FIG,"FigureS1_global_network.png"), pA, width=10, height=7.2, dpi=400, bg="white")
ggsave(file.path(FIG,"FigureS1_global_network.pdf"), pA, width=10, height=7.2, bg="white")
lg("wrote FigureS1_global_network  (nodes",vcount(g),"edges",ecount(g),")")

# ---------------------------------------------------------------- Fig B: FFL census
cen <- rbindlist(lapply(3:6, function(k){
  f <- file.path(REV, sprintf("results/ffl_census_k%d.json", k))
  if (!file.exists(f)) return(NULL)
  d <- fromJSON(f)
  data.table(n=k, subgraphs=d$est_subgraphs, ffl=d$est_ffl,
             exhaustive = all(unlist(d$probs)==1))
}))
if (nrow(cen)) {
  fwrite(cen, file.path(REV,"results/ffl_census_summary.csv"))
  pB <- ggplot(cen, aes(factor(n), ffl)) +
    geom_col(aes(fill=exhaustive), width=0.62) +
    geom_text(aes(label=label_number(scale_cut=cut_short_scale())(ffl)),
              vjust=-0.5, size=3.1) +
    scale_y_log10(labels=label_number(scale_cut=cut_short_scale()),
                  expand=expansion(mult=c(0,0.16))) +
    scale_fill_manual(values=c(`TRUE`="#0F172A",`FALSE`="#64748B"),
                      labels=c(`TRUE`="exhaustive",`FALSE`="RAND-ESU estimate"),
                      name=NULL) +
    labs(x="Module size (nodes)", y="n-node FFLs (log scale)",
         title="Systematic census of n-node feed-forward loops",
         subtitle="Modules satisfying D1-D4; counts grow ~10x per added node") +
    theme_minimal(base_size=10) +
    theme(panel.grid.major.x=element_blank(), legend.position="top")
  ggsave(file.path(FIG,"Figure2B_ffl_census.png"), pB, width=5.4, height=4.2, dpi=400, bg="white")
  lg("wrote Figure2B_ffl_census"); print(cen)
}

# ---------------------------------------------------------------- Fig C: directed module diagrams
draw_module <- function(members, title) {
  sub <- E[source %in% members & target %in% members]
  if (!nrow(sub)) return(NULL)
  sub[, sign := ifelse(edge_type=="miRNA_target", "repress", "activate")]
  vs <- nodes[name %in% members, .(name,type)]
  gg <- as_tbl_graph(graph_from_data_frame(sub[,.(source,target,edge_type,sign)],
                                           directed=TRUE, vertices=vs))
  ggraph(gg, layout="sugiyama") +
    geom_edge_link(aes(colour=edge_type, linetype=sign),
                   arrow=arrow(length=unit(2.4,"mm"), type="closed"),
                   end_cap=circle(6,"mm"), start_cap=circle(6,"mm"),
                   width=0.5, alpha=0.85) +
    geom_node_point(aes(fill=type, shape=type), size=8.5, colour="grey20", stroke=0.3) +
    geom_node_text(aes(label=name), size=2.15, colour="white", fontface="bold") +
    scale_shape_manual(values=c(TF=21, miRNA=24, Gene=22), name="Class") +
    scale_fill_manual(values=PAL, name="Class") +
    scale_edge_colour_manual(values=EPAL, name="Interaction") +
    scale_edge_linetype_manual(values=c(activate="solid", repress="dashed"),
                               name="Effect") +
    labs(title=title) + theme_void(base_size=9) +
    theme(plot.title=element_text(size=9, face="bold"), legend.position="right")
}

# the authors' own published exemplars, redrawn correctly and directed
ex <- list(
  `4-node exemplar (published Fig. 5)` = c("ETS1","NFKB1","RELA","SP1","COL1A1","COL3A1",
      "hsa-let-7b","hsa-let-7e","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"),
  `6-node exemplar (published Fig. 4)` = c("NFKB1","RELA","SP1","COL1A1","COL3A1",
      "hsa-let-7b","hsa-let-7e","hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"))
ps <- Filter(Negate(is.null), Map(draw_module, ex, names(ex)))
if (length(ps)) {
  ggsave(file.path(FIG,"Figure3_published_exemplars_directed.png"),
         wrap_plots(ps, ncol=1) + plot_annotation(
           title="Published higher-order circuits, redrawn as directed signed graphs",
           subtitle="Dashed = repression. Note the TFs act on COL1A1 while the miRNAs act on COL3A1,\nso no target is co-regulated by both: these circuits contain no 3-node FFL core.",
           theme=theme(plot.title=element_text(size=11,face="bold"),
                       plot.subtitle=element_text(size=8.4, colour="grey30"))),
         width=8.4, height=9.6, dpi=400, bg="white")
  lg("wrote Figure3_published_exemplars_directed")
}
lg("DONE")
