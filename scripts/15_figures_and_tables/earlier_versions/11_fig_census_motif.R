# Figure: FFL census across module sizes, edge-class saturation, and motif significance
# against three null models.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/15_figures_and_tables/earlier_versions/10_theme.R"))
FIG <- file.path(REV,"figures/v5"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

# ---- panel A: census -----------------------------------------------------------------
cen <- data.table(n=3:7,
                  ffl=c(6037, 1.20e5, 1.54e6, 1.94e7, 2.70e8),
                  exact=c(TRUE,FALSE,FALSE,FALSE,FALSE),
                  maxclass=c(3,5,6,6,6))
pA <- ggplot(cen, aes(factor(n), ffl)) +
  geom_col(aes(fill=exact), width=0.62) +
  geom_text(aes(label=label_number(scale_cut=cut_short_scale(), accuracy=0.1)(ffl)),
            vjust=-0.6, size=2.7, colour=INK) +
  scale_y_log10(labels=label_number(scale_cut=cut_short_scale()),
                expand=expansion(mult=c(0,0.18))) +
  scale_fill_manual(values=c(`TRUE`=INK, `FALSE`=MUTED),
                    labels=c(`TRUE`="exhaustive", `FALSE`="RAND-ESU estimate"), name=NULL) +
  labs(x="module size (nodes)", y=expression(italic(n)*"-node FFLs  (log scale)"),
       title="a   Census of n-node feed-forward loops",
       subtitle="counts rise ~10-fold per added node") +
  theme_pub() + theme(panel.grid.major.x=element_blank(), legend.position="top", legend.justification="left",
                      legend.background=element_rect(fill="white", colour=NA))

# ---- panel B: edge-class saturation --------------------------------------------------
pB <- ggplot(cen, aes(factor(n), maxclass)) +
  geom_col(fill=SEM["accent"], width=0.62, alpha=0.85) +
  geom_hline(yintercept=6, linetype="dashed", colour=SEM["neg"], linewidth=0.4) +
  annotate("text", x=0.6, y=6.32, label="saturation: 6 of 7 classes",
           size=2.4, colour=SEM["neg"], hjust=0) +
  annotate("segment", x=2.35, xend=2.9, y=7.15, yend=6.2,
           arrow=arrow(length=unit(1.5,"mm"), type="closed"), colour=INK, linewidth=0.35) +
  annotate("text", x=2.3, y=7.35, label="first reached at n = 5", size=2.5, colour=INK, hjust=1) +
  scale_y_continuous(limits=c(0,7.6), breaks=0:7, expand=expansion(mult=c(0,0.02))) +
  labs(x="module size (nodes)", y="max. distinct edge classes per module",
       title="b   Class diversity saturates at five nodes",
       subtitle="above n = 5, modules gain multiplicity\nnot new interaction types") +
  theme_pub() + theme(panel.grid.major.x=element_blank())

# ---- panel C: motif significance -----------------------------------------------------
m <- fread(file.path(REV,"results/v2/verify_motif.csv"))
m <- m[graph=="dep" & metric %in% c("Composite-FFL (once per reciprocal pair)","TF-FFL",
                                    "miRNA-FFL","All 3-node FFL modules","reciprocal TF<->miRNA pairs")]
m[, null := fcase(grepl("NULL-A", null_desc), "A: full randomisation",
                  grepl("NULL-B", null_desc), "B: reciprocity preserved",
                  grepl("NULL-C", null_desc), "C: B + bipartite curveball")]
m[, lab := fcase(metric=="Composite-FFL (once per reciprocal pair)","Composite-FFL",
                 metric=="All 3-node FFL modules","all 3-node FFLs",
                 metric=="reciprocal TF<->miRNA pairs","reciprocal TF-miRNA pairs",
                 default=metric)]
m[, lab := factor(lab, levels=c("reciprocal TF-miRNA pairs","Composite-FFL","all 3-node FFLs",
                                "miRNA-FFL","TF-FFL"))]
m[, dir := fifelse(fold>1.02, "enriched", fifelse(fold<0.98, "depleted", "not different"))]
pC <- ggplot(m, aes(x=log2(fold), y=lab, fill=dir)) +
  geom_vline(xintercept=0, colour=RULE, linewidth=0.5) +
  geom_col(width=0.6) +
  geom_text(aes(label=sprintf("%.2f×", fold),
                hjust=ifelse(log2(fold)>0, -0.15, 1.15)), size=2.4, colour=INK) +
  facet_wrap(~null, nrow=1) +
  scale_fill_manual(values=c(enriched=unname(SEM["pos"]), depleted=unname(SEM["neg"]),
                             `not different`=unname(SEM["null"])), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0.22,0.22))) +
  labs(x=expression(log[2]*" fold change vs randomised networks"), y=NULL,
       title="c   Motif over-representation depends entirely on the null model",
       subtitle="1,000 randomisations per null. Null A destroys the reciprocal TF-miRNA pairing that the study's own miRNA-TF filter imposes, and so\nmanufactures a 3.3-fold composite enrichment while simultaneously making the total FFL count come out depleted - an incoherent result\nthat is diagnostic of an inappropriate null.",
       caption="Under the two reciprocity-preserving nulls the genuine excess is 6-16 % (all 3-node FFLs: +13 % and +16 %, both p = 0.001). The network carries\n4.9-fold more reciprocal TF-miRNA pairs than chance (Z = 88), which is the hypergeometric miRNA-TF filter rather than biology.") +
  theme_pub() + theme(legend.position="top", panel.grid.major.y=element_blank())

out <- (pA | pB) / pC + plot_layout(heights=c(1, 1.15))
ggsave(file.path(FIG,"Fig3_census_and_motif_significance.png"), out, width=11.2, height=8.2, dpi=400, bg="white")
ggsave(file.path(FIG,"Fig3_census_and_motif_significance.pdf"), out, width=11.2, height=8.2, bg="white")
cat("wrote Fig3_census_and_motif_significance\n")
