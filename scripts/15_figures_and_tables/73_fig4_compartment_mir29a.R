# Fig. 4: (a) mediation of the TF-collagen and miR-29-collagen associations by stromal content across 42 estimates;
# (b) spatial compartment occupancy, (c) patient-to-xenograft change and (d) bulk versus within-stroma TF-COL1A1
# correlation; (e) miR-29a survival meta-analysis. One theme (theme_pub), 8 pt text, 174 mm, 600 dpi.
# Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales); library(ragg)})
REV <- Sys.getenv("FFL_REV", unset = "/path/to/revision")
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
FIG <- file.path(REV,"figures/final"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
BASE <- 9.5; TXT <- 8 / .pt
R6 <- file.path(REV,"results/v6")
ital <- theme(axis.text.y = element_text(face = "italic"))   # gene symbols on the y axis
HI <- "#B03A2E"; BL <- "#2E5E8A"; OR <- "#D98C5F"; GR <- "#7A7A7A"

# ---- a: mediation across 42 estimates --------------------------------------------------
mm <- fread(file.path(REV,"results/v3/deconv_mediation_extended.csv")); mm <- mm[!is.na(ACME) & !is.na(ADE)]
mm <- mm[!(mediator %in% c("MCPcounter_Fibroblasts_recomputed","NNLS_Wu641_full_CAFs","QPROG_Wu641_full_CAFs"))]
want <- data.table(x=c("ETS1","NFKB1","SP1","RELA","hsa-miR-29a","hsa-miR-29a","hsa-miR-29b"),
                   y=c("COL1A1","COL1A1","COL1A1","COL1A1","COL3A1","COL1A1","COL3A1"),
                   path=c("ETS1>COL1A1","NFKB1>COL1A1","SP1>COL1A1","RELA>COL1A1","miR-29a|COL3A1","miR-29a|COL1A1","miR-29b|COL3A1"),
                   expr=c("italic('ETS1')~' → '~italic('COL1A1')","italic('NFKB1')~' → '~italic('COL1A1')","italic('SP1')~' → '~italic('COL1A1')","italic('RELA')~' → '~italic('COL1A1')",
                          "'miR-29a'~' ⊣ '~italic('COL3A1')","'miR-29a'~' ⊣ '~italic('COL1A1')","'miR-29b'~' ⊣ '~italic('COL3A1')"))
med <- merge(want, mm, by=c("x","y"), all.x=TRUE)[, .(acme=median(ACME), lo=min(ACME), hi=max(ACME), ade=median(ADE),
      n_inc=sum(toupper(as.character(inconsistent)) %in% c("TRUE","1")), n=.N), by=path]
med <- med[match(want$path, path)]
med[, kind := ifelse(grepl("miR", path), "miRNA arm", "TF arm")]; med[, lab := paste0(n_inc, "/", n)]
med[, path := factor(path, levels=rev(want$path))]
labmap <- setNames(want$expr, want$path)
pa <- ggplot(med, aes(y=path)) + geom_vline(xintercept=0, colour=RULE, linewidth=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.16, colour=MUTED, linewidth=0.5) +
  geom_point(aes(x=acme, colour=kind), size=2.8) + geom_point(aes(x=ade), shape=4, size=2.2, colour=INK, stroke=0.7) +
  geom_text(aes(x=0.80, label=lab), size=TXT, colour=MUTED, hjust=0) +
  scale_colour_manual(values=c(`TF arm`=unname(SEM["neg"]), `miRNA arm`=unname(SEM["pos"])), name=NULL) +
  scale_y_discrete(labels=function(v) parse(text=unname(labmap[v]))) +
  scale_x_continuous(breaks=seq(-0.25,0.75,0.25), labels=label_number(accuracy=0.01, style_negative="minus"),
                     limits=c(-0.35,0.90), expand=expansion(0)) + labs(x="effect size", y=NULL, tag="a") +
  theme_pub(BASE) + theme(legend.position="top", legend.justification="left", panel.grid.major.y=element_blank())

# ---- b: spatial compartment AUC ---------------------------------------------------------
a <- fread(file.path(R6, "spatial_compartment_marker_enrichment_meta.csv"))
a <- a[scheme == "marker_tertile"]
keep <- c("COL1A1","DCN","COL3A1","POSTN","CXCL12","FN1","ACTA2","PDGFRB","ETS1","EGR2","NFKB1","BRCA1","E2F1","EZH2","SREBF1","ESR1","GATA3","EPCAM","KRT8")
a <- a[gene %in% keep]
a[, cls := fifelse(gene %in% c("COL1A1","COL3A1"), "collagen",
            fifelse(gene %in% c("DCN","POSTN","CXCL12","FN1","ACTA2","PDGFRB"), "stromal",
            fifelse(gene %in% c("EPCAM","KRT8"), "epithelial marker", "TF")))]
a[, gene := factor(gene, levels = a[order(mean_auc)]$gene)]
pb <- ggplot(a, aes(mean_auc, gene, colour = cls)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey55", linewidth = .3) +
  geom_segment(aes(x = .5, xend = mean_auc, yend = gene), linewidth = .5) + geom_point(size = 1.9) +
  scale_colour_manual(values = c(collagen = HI, stromal = OR, TF = BL, `epithelial marker` = GR), name = NULL) +
  scale_x_continuous(limits = c(.15, .92)) + labs(x = "AUC (stromal vs epithelial spots)", y = NULL, tag = "b") +
  theme_pub(BASE) + ital + theme(panel.grid.major.y = element_blank(), legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6)))

# ---- c: patient-to-xenograft percentile change ------------------------------------------
b <- fread(file.path(R6, "pdx_matched_pairs_rank_change.csv"))
sel <- c("COL3A1","POSTN","PDGFRB","CXCL12","COL1A1","FN1","DCN","LUM","ACTA2","PTPRC","E2F1","EZH2","MYBL2","BRCA1","EPCAM","KRT8","ERBB2","MKI67")
b <- b[gene %in% sel]
b[, cls := fifelse(gene %in% c("COL3A1","POSTN","PDGFRB","CXCL12","COL1A1","FN1"), "stromal arm",
            fifelse(gene %in% c("E2F1","EZH2","MYBL2","BRCA1"), "proliferative arm",
            fifelse(gene %in% c("EPCAM","KRT8","ERBB2","MKI67"), "carcinoma control", "stromal/immune control")))]
b[, gene := factor(gene, levels = b[order(median_paired_delta)]$gene)]
pc <- ggplot(b, aes(y = gene, colour = cls)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = .3) +
  geom_segment(aes(x = median_pctrank_originator, xend = median_pctrank_pdx, yend = gene),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .45) +
  geom_point(aes(x = median_pctrank_originator), size = 1.4) +
  scale_colour_manual(values = c(`stromal arm` = HI, `stromal/immune control` = OR, `proliferative arm` = BL, `carcinoma control` = GR), name = NULL) +
  scale_x_continuous(limits = c(0, 102)) + labs(x = "expression percentile within sample", y = NULL, tag = "c") +
  theme_pub(BASE) + ital + theme(panel.grid.major.y = element_blank(), legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6)))

# ---- d: bulk versus within-stroma TF-COL1A1 correlation ---------------------------------
cc <- fread(file.path(R6, "spatial_tf_collagen_correlation_meta.csv"))
cc <- cc[scheme == "marker_tertile" & outcome == "COL1A1"]
cc[, gene := factor(tf, levels = cc[order(rho_tcga_bulk)]$tf)]
pd <- ggplot(cc, aes(y = gene)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = .3) +
  geom_segment(aes(x = rho_tcga_bulk, xend = rho_stromal_meta, yend = gene), colour = "grey65", linewidth = .5) +
  geom_point(aes(x = rho_tcga_bulk), colour = BL, size = 1.8) + geom_point(aes(x = rho_stromal_meta), colour = HI, size = 1.8) +
  annotate("text", x = -.33, y = 12.9, label = "TCGA bulk", colour = BL, size = TXT, hjust = 0) +
  annotate("text", x =  .33, y = 12.9, label = "within stromal spots", colour = HI, size = TXT, hjust = 1) +
  scale_x_continuous(labels = label_number(accuracy = 0.1, style_negative = "minus")) +
  coord_cartesian(xlim = c(-.33, .33), ylim = c(1, 13.2), clip = "off") +
  labs(x = expression(Spearman~rho~with~italic("COL1A1")), y = NULL, tag = "d") +
  theme_pub(BASE) + ital + theme(panel.grid.major.y = element_blank())

# ---- e: miR-29a survival meta-analysis ------------------------------------------------
d <- fread(file.path(R6, "mirna_meta_forest_table.csv"))
d <- d[miRNA == "miR-29a" & analysis == "primary_endpoint_univariate"]
d[, label := sub("^TCGA_BRCA", "TCGA-BRCA", label)]; d[row_type != "study", label := "pooled"]
d[, label := factor(label, levels = rev(c(d[row_type == "study"]$label, d[row_type != "study"]$label)))]
d[, isp := row_type != "study"]
pe <- ggplot(d, aes(HR, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .45, colour = "grey35") +
  geom_point(aes(size = ifelse(isp, 3.1, 1.9), shape = isp, colour = ifelse(isp, HI, BL))) +
  scale_colour_identity() + scale_size_identity() + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_x_continuous(trans = "log", breaks = c(.6, .8, 1, 1.25)) +
  labs(x = "hazard ratio per SD (95% CI)", y = NULL, tag = "e") + theme_pub(BASE) + theme(panel.grid.major.y = element_blank())

out <- pa / (pb | pc) / (pd | pe) + plot_layout(heights = c(0.85, 1.25, 1))
ggsave(file.path(FIG,"Fig4_compartment_mir29a.png"), out, width = 6.85, height = 9.0, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG,"Fig4_compartment_mir29a.pdf"), out, width = 6.85, height = 9.0, bg = "white", device = cairo_pdf)
cat("wrote Fig4_compartment_mir29a\n")
