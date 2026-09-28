# Figure 4: spatial compartment occupancy, patient-to-xenograft change, bulk versus within-stroma
# TF-COL1A1 correlation and the miR-29a survival meta-analysis.
# Revised 2026-09-25: drawn at print width (174 mm) with 8 pt text (9 pt axis titles); gene symbols
# italic; cohort labels as in the text. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
R6 <- "results/v6"; OUT <- "figures/final"
TXT <- 8 / .pt   # in-panel text 8 pt
th <- theme_bw(base_size = 9) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.tag = element_text(face = "bold", size = 11.5),
  plot.tag.position = "topleft",
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  axis.text = element_text(size = 8), axis.title = element_text(size = 9), legend.position = "none")
ital <- theme(axis.text.y = element_text(face = "italic"))   # gene symbols on the y axis

## ---- (a) spatial compartment AUC -------------------------------------------
a <- fread(file.path(R6, "spatial_compartment_marker_enrichment_meta.csv"))
a <- a[scheme == "marker_tertile"]
keep <- c("COL1A1","DCN","COL3A1","POSTN","CXCL12","FN1","ACTA2","PDGFRB","ETS1",
          "EGR2","NFKB1","BRCA1","E2F1","EZH2","SREBF1","ESR1","GATA3","EPCAM","KRT8")
a <- a[gene %in% keep]
a[, cls := fifelse(gene %in% c("COL1A1","COL3A1"), "collagen",
            fifelse(gene %in% c("DCN","POSTN","CXCL12","FN1","ACTA2","PDGFRB"), "stromal",
            fifelse(gene %in% c("EPCAM","KRT8"), "epithelial marker", "TF")))]
a[, gene := factor(gene, levels = a[order(mean_auc)]$gene)]
pal <- c(collagen = "#B03A2E", stromal = "#D98C5F", TF = "#2E5E8A", `epithelial marker` = "#7A7A7A")
pa <- ggplot(a, aes(mean_auc, gene, colour = cls)) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey55", linewidth = .3) +
  geom_segment(aes(x = .5, xend = mean_auc, yend = gene), linewidth = .5) +
  geom_point(size = 1.9) + scale_colour_manual(values = pal) +
  scale_x_continuous(limits = c(.15, .92)) +
  labs( tag = "a",
       x = "AUC", y = NULL) + th + ital

## ---- (b) PDX paired percentile change --------------------------------------
b <- fread(file.path(R6, "pdx_matched_pairs_rank_change.csv"))
sel <- c("COL3A1","POSTN","PDGFRB","CXCL12","COL1A1","FN1","DCN","LUM","ACTA2","PTPRC",
         "E2F1","EZH2","MYBL2","BRCA1","EPCAM","KRT8","ERBB2","MKI67")
b <- b[gene %in% sel]
b[, cls := fifelse(gene %in% c("COL3A1","POSTN","PDGFRB","CXCL12","COL1A1","FN1"), "stromal arm",
            fifelse(gene %in% c("E2F1","EZH2","MYBL2","BRCA1"), "proliferative arm",
            fifelse(gene %in% c("EPCAM","KRT8","ERBB2","MKI67"), "carcinoma control", "stromal/immune control")))]
b[, gene := factor(gene, levels = b[order(median_paired_delta)]$gene)]
pal2 <- c(`stromal arm` = "#B03A2E", `stromal/immune control` = "#D98C5F",
          `proliferative arm` = "#2E5E8A", `carcinoma control` = "#7A7A7A")
pb <- ggplot(b, aes(y = gene, colour = cls)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = .3) +
  geom_segment(aes(x = median_pctrank_originator, xend = median_pctrank_pdx, yend = gene),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .45) +
  geom_point(aes(x = median_pctrank_originator), size = 1.4) +
  scale_colour_manual(values = pal2) + scale_x_continuous(limits = c(0, 102)) +
  labs( tag = "b",
       x = "expression percentile within sample", y = NULL) + th + ital +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8), legend.key.size = unit(3, "mm"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.6)))

## ---- (c) bulk vs within-stroma TF-COL1A1 correlation -----------------------
cc <- fread(file.path(R6, "spatial_tf_collagen_correlation_meta.csv"))
cc <- cc[scheme == "marker_tertile" & outcome == "COL1A1"]
cc[, gene := factor(tf, levels = cc[order(rho_tcga_bulk)]$tf)]
pc <- ggplot(cc, aes(y = gene)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = .3) +
  geom_segment(aes(x = rho_tcga_bulk, xend = rho_stromal_meta, yend = gene),
               colour = "grey65", linewidth = .5) +
  geom_point(aes(x = rho_tcga_bulk), colour = "#2E5E8A", size = 1.8) +
  geom_point(aes(x = rho_stromal_meta), colour = "#B03A2E", size = 1.8) +
  annotate("text", x = -.33, y = 12.75, label = "TCGA bulk", colour = "#2E5E8A", size = TXT, hjust = 0) +
  annotate("text", x =  .33, y = 12.75, label = "within stromal spots", colour = "#B03A2E", size = TXT, hjust = 1) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1, style_negative = "minus")) +
  coord_cartesian(xlim = c(-.33, .33), ylim = c(1, 13.1), clip = "off") +
  labs( tag = "c",
       x = expression(Spearman~rho~with~italic("COL1A1")), y = NULL) + th + ital

## ---- (d) miR-29a forest ----------------------------------------------------
d <- fread(file.path(R6, "mirna_meta_forest_table.csv"))
d <- d[miRNA == "miR-29a" & analysis == "primary_endpoint_univariate"]
d[, label := sub("^TCGA_BRCA", "TCGA-BRCA", label)]; d[row_type != "study", label := "pooled"]
d[, label := factor(label, levels = rev(c(d[row_type == "study"]$label, d[row_type != "study"]$label)))]
d[, isp := row_type != "study"]
pd <- ggplot(d, aes(HR, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .45, colour = "grey35") +
  geom_point(aes(size = ifelse(isp, 3.1, 1.9), shape = isp,
                 colour = ifelse(isp, "#B03A2E", "#2E5E8A"))) +
  scale_colour_identity() + scale_size_identity() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
  scale_x_continuous(trans = "log", breaks = c(.6, .8, 1, 1.25)) +
  labs(tag = "d",
       x = "hazard ratio per SD (95% CI)", y = NULL) + th

fig <- (pa | pb) / (pc | pd) + plot_annotation(theme = theme(plot.margin = margin(2,2,2,2)))
ggsave(file.path(OUT, "Fig8_independent_validation.png"), fig, width = 6.85, height = 6.6, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(OUT, "Fig8_independent_validation.pdf"), fig, width = 6.85, height = 6.6, bg = "white", device = cairo_pdf)
cat("Fig8 written\n")
