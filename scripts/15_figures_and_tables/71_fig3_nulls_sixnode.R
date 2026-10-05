# Fig. 3 of the research-article version (1 Oct 2026): (a) over-representation of three-node FFLs and of reciprocal TF-miRNA pairs
# against the three null models (formerly Fig. 1c, same code and data); (b) the six-node composite pattern of Bhat et al. (2024)
# under NULL-C and NULL-L, as tested and in the three variants that limit the result (new panel; every value is read from
# results/v6/fig3b_six_node_pattern.csv, whose last column names the source file and line). Drawn at 174 mm, 8 pt text, 600 dpi.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)})
REV <- Sys.getenv("FFL_REV", unset = "/path/to/revision")
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
BASE <- 9.5; TXT <- 8 / .pt

# ---- a: three null models ---------------------------------------------------------------
RR <- file.path(REV,"results/v2/census_rerun")
m <- fread(file.path(RR,"motif_nulls_all_graphs.csv"))
m <- m[graph=="dep_nolegacy" & metric %in% c("Composite-FFL (once per reciprocal pair)","TF-FFL",
                                    "miRNA-FFL","All 3-node FFL modules","reciprocal TF<->miRNA pairs")]
m[, null := fcase(grepl("NULL-A", null_desc), "A: degree-preserving",
                  grepl("NULL-B", null_desc), "B: reciprocity preserved",
                  grepl("NULL-C", null_desc), "C: B + bipartite curveball")]
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
b <- fread(file.path(REV,"results/v6/fig3b_six_node_pattern.csv"))
setnames(b, "ratio_observed_over_null_mean", "ratio")
b[, sig := p_upper < 0.05]
b[, null := factor(null, levels=c("NULL-C","NULL-L"))]
lv <- b[, .(row_order=unique(row_order)), by=variant][order(-row_order)]$variant
b[, variant := factor(variant, levels=lv)]
b[, vlab := stringr::str_wrap(as.character(variant), 30)]
lvl <- sapply(lv, function(v) stringr::str_wrap(v, 30))   # first row of the table at the top
b[, vlab := factor(vlab, levels=lvl)]
pB <- ggplot(b, aes(ratio, vlab)) +
  geom_vline(xintercept=1, colour=MUTED, linewidth=0.4, linetype="dashed") +
  geom_line(aes(group=vlab), colour=RULE, linewidth=1.6) +
  geom_point(aes(shape=null, fill=sig, colour=sig), size=2.8, stroke=0.6) +
  geom_text(aes(label=sprintf("%.2f×", ratio), vjust=ifelse(null=="NULL-C", 1.9, -1.15)), size=TXT, colour=INK) +
  scale_shape_manual(values=c(`NULL-C`=21, `NULL-L`=23), name="null model") +
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
