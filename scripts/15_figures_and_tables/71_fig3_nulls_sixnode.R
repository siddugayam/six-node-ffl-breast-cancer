# Fig. 3: (a) over-representation of three-node FFLs and of reciprocal TF-miRNA pairs against the three null models;
# (b) the six-node composite pattern of Bhat et al. (2024) under the bipartite-layer and core-conditional nulls, as tested and in
# three variants (physical STRING links only, each co-transcribed cluster counted once, six-node modules under the topological
# definition; every value is read from
# results/tissue_and_cohorts/fig3b_six_node_pattern.csv, whose last column names the source file and line). Drawn at 174 mm, 8 pt text, 600 dpi.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)})
REV <- Sys.getenv("FFL_REV", unset = "/path/to/six-node-ffl-breast-cancer")
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
BASE <- 9.5; TXT <- 8 / .pt

# ---- a: three null models ---------------------------------------------------------------
RR <- file.path(REV,"results/expression_validation/census_tables")
m <- fread(file.path(RR,"motif_nulls_all_graphs.csv"))
m <- m[graph=="dep_noexemplar" & metric %in% c("Composite-FFL (once per reciprocal pair)","TF-FFL",
                                    "miRNA-FFL","All 3-node FFL modules","reciprocal TF<->miRNA pairs")]
m[, null := fcase(grepl("NULL-A", null_desc), "Degree-preserving null",
                  grepl("NULL-B", null_desc), "Reciprocity-preserving null",
                  grepl("NULL-C", null_desc), "Bipartite-layer null")]
m[, null := factor(null, levels=c("Degree-preserving null","Reciprocity-preserving null","Bipartite-layer null"))]
m[, lab := fcase(metric=="Composite-FFL (once per reciprocal pair)","Composite-FFL",
                 metric=="All 3-node FFL modules","all 3-node FFLs",
                 metric=="reciprocal TF<->miRNA pairs","reciprocal TF–miRNA pairs",
                 default=metric)]
m[, lab := factor(lab, levels=c("reciprocal TF–miRNA pairs","Composite-FFL","all 3-node FFLs","miRNA-FFL","TF-FFL"))]
m[, dir := fifelse(fold>1.02, "enriched", fifelse(fold<0.98, "depleted", "not different"))]
pA <- ggplot(m, aes(x=log2(fold), y=lab, fill=dir)) +
  geom_vline(xintercept=0, colour=RULE, linewidth=0.5) + geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.2f×", fold), hjust=ifelse(log2(fold)>0, -0.15, 1.15)), size=TXT, colour=INK) +
  facet_wrap(~null, nrow=1, labeller=labeller(null=label_wrap_gen(22))) +
  scale_fill_manual(values=c(enriched=unname(SEM["pos"]), depleted=unname(SEM["neg"]), `not different`=unname(SEM["null"])), name=NULL) +
  scale_x_continuous(breaks=c(-5,0,5), labels=label_number(accuracy=1, style_negative="minus"), expand=expansion(mult=c(0.5,0.42))) +
  labs(x=expression(log[2]*" fold change vs randomised networks"), y=NULL, tag="a") +
  theme_pub(BASE) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())

# ---- b: six-node composite pattern ----------------------------------------------------
b <- fread(file.path(REV,"results/tissue_and_cohorts/fig3b_six_node_pattern.csv"))
setnames(b, "ratio_observed_over_null_mean", "ratio")
b[, sig := p_upper < 0.05]
b[, variant := sub("topological-definition six-node modules", "six-node modules of the topological definition", variant)]
b[, variant := gsub("gene-gene", "gene–gene", variant, fixed=TRUE)]   # en dash, as in the text
b[, null := factor(fcase(null=="NULL-C","Bipartite-layer null", null=="NULL-L","Core-conditional null"),
                      levels=c("Bipartite-layer null","Core-conditional null"))]
lv <- b[, .(row_order=unique(row_order)), by=variant][order(-row_order)]$variant
b[, variant := factor(variant, levels=lv)]
b[, vlab := stringr::str_wrap(as.character(variant), 30)]
lvl <- sapply(lv, function(v) stringr::str_wrap(v, 30))   # first row of the table at the top
b[, vlab := factor(vlab, levels=lvl)]
pB <- ggplot(b, aes(ratio, vlab)) +
  geom_vline(xintercept=1, colour=MUTED, linewidth=0.4, linetype="dashed") +
  geom_line(aes(group=vlab), colour=RULE, linewidth=1.6) +
  geom_point(aes(shape=null, fill=sig, colour=sig), size=2.8, stroke=0.6) +
  geom_text(aes(label=sprintf("%.2f×", ratio), vjust=ifelse(null=="Bipartite-layer null", 1.9, -1.15)), size=TXT, colour=INK) +
  scale_shape_manual(values=c(`Bipartite-layer null`=21, `Core-conditional null`=23), name="null model") +
  scale_fill_manual(values=c(`TRUE`=unname(SEM["pos"]), `FALSE`="white"), labels=c(`TRUE`="p < 0.05", `FALSE`="p ≥ 0.05"), name="upper-tail p") +
  scale_colour_manual(values=c(`TRUE`=unname(SEM["pos"]), `FALSE`=unname(SEM["null"])), guide="none") +
  scale_x_continuous(trans="log2", breaks=c(0.5,1,2,4), limits=c(0.7,5), expand=expansion(mult=c(0.02,0.04))) +
  guides(fill=guide_legend(order=1, override.aes=list(shape=21, colour=unname(c(SEM["null"], SEM["pos"])), fill=c("white", unname(SEM["pos"])))),
         shape=guide_legend(order=2, override.aes=list(fill="grey60", colour="grey30"))) +
  labs(x="observed count / mean count in randomised networks (log scale)", y=NULL, tag="b") +
  theme_pub(BASE) + theme(legend.position="bottom", legend.box="horizontal", legend.justification="left", panel.grid.major.y=element_blank())
out <- pA / pB + plot_layout(heights=c(1.15, 1))
ggsave(file.path(FIG,"Fig3_nulls_sixnode.png"), out, width=6.85, height=5.6, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(FIG,"Fig3_nulls_sixnode.pdf"), out, width=6.85, height=5.6, bg="white", device=cairo_pdf)
cat("wrote Fig3_nulls_sixnode\n")
