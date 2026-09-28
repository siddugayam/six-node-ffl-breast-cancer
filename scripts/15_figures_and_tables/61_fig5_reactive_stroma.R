# Figure 5: FFL module scores against the Farmer stroma-related signature, and the pCR meta-analysis.
# Revised 2026-09-25: drawn at print width (174 mm) with 8 pt text (9 pt axis titles); gene symbols
# italic; signature names as in the legend. Run in a UTF-8 locale (LC_ALL=en_US.UTF-8).
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
REV <- "/path/to/revision"
source(file.path(REV,"scripts/15_figures_and_tables/10_theme.R"))
OUT <- file.path(REV,"figures/final")
BASE <- 9.5; TXT <- 8 / .pt   # theme_pub(9.5): tick labels 8 pt; in-panel text 8 pt
# y-axis labels: gene symbols italic (plotmath), everything else as plain text
ylab <- function(v) parse(text=ifelse(grepl("^COL1A1 \\+ COL3A1", v),
  paste0("italic('COL1A1')*' + '*italic('COL3A1')*'", sub("^COL1A1 \\+ COL3A1", "", v), "'"),
  paste0("'", v, "'")))

## ---- (a) module scores vs the Farmer stroma signature -----------------------
a <- fread(file.path(REV,"results/v5/farmer_module_score_correlation.csv"))[method=="ssGSEA"]
## add the two reference points the legend promises: the collagen pair and an epithelial panel
ref <- fread(file.path(REV,"results/v5/farmer_correlations_tcga.csv"))[
  block=="signature_score" & farmer_score=="ssGSEA" &
  partner %in% c("COLLAGEN_PAIR","FINAK_REF_EPITHELIUM_E")]
if (!nrow(ref)) ref <- fread(file.path(REV,"results/v5/farmer_correlations_tcga.csv"))[
  block=="signature_score" & farmer_score=="meanZ" &
  partner %in% c("COLLAGEN_PAIR","FINAK_REF_EPITHELIUM_E")]
a <- rbind(a, data.table(module=ref$partner, method="ssGSEA", spearman_rho=ref$rho), fill=TRUE)
a[, lab := fifelse(module=="FFL_3node","3-node FFL",
            fifelse(module=="FFL_higher_order_only","higher-order only",
            fifelse(module=="MS_exemplar_module","exemplar module",
            fifelse(module=="COLLAGEN_PAIR","COL1A1 + COL3A1 (reference)",
            fifelse(module=="FINAK_REF_EPITHELIUM_E","epithelial panel (reference)",
            fifelse(module=="FFL_4node","4-node FFL",
            fifelse(module=="FFL_5node","5-node FFL",
            fifelse(module=="FFL_6node_all","6-node FFL", gsub("_"," ",module)))))))))]
a <- a[order(spearman_rho)]
a[, lab := factor(lab, levels=lab)]
a[, hl := grepl("COL1A1|epithelial", lab)]
pa <- ggplot(a, aes(spearman_rho, lab)) +
  geom_vline(xintercept=0, colour=RULE, linewidth=0.4) +
  geom_segment(aes(x=0, xend=spearman_rho, yend=lab, colour=hl), linewidth=0.6) +
  geom_point(aes(colour=hl), size=2.2) +
  geom_text(aes(label=sub("-", "\u2212", sprintf("%.2f", spearman_rho)), hjust=ifelse(spearman_rho < 0, 1.35, -0.35)),
            size=TXT, colour=MUTED) +
  scale_colour_manual(values=c(`TRUE`=unname(SEM["neg"]), `FALSE`=unname(SEM["null"])), guide="none") +
  scale_x_continuous(limits=c(-0.15,1.0), expand=expansion(mult=c(0.01,0.12))) +
  scale_y_discrete(labels=ylab) +
  labs(x=expression(Spearman~rho~with~the~Farmer~stroma~signature), y=NULL, tag="a") +
  theme_pub(BASE) + theme(panel.grid.major.y=element_blank())

## ---- (b) pCR meta-analysis --------------------------------------------------
b <- fread(file.path(REV,"results/v5/farmer_neoadjuvant_meta.csv"))
keep <- c("COLLAGEN_PAIR","MS_exemplar_module","FFL_3node","FFL_higher_order_only","FFL_6node_all",
          "FARMER_STROMAL","CHANG_WOUND_UP_VANTVEER","HALLMARK_EMT","FINAK_SDPP","CAF_scRNA_50")
b <- b[score %in% keep]
b[, lab := fifelse(score=="COLLAGEN_PAIR","COL1A1 + COL3A1",
            fifelse(score=="MS_exemplar_module","exemplar module",
            fifelse(score=="FFL_3node","3-node FFL",
            fifelse(score=="FFL_higher_order_only","higher-order only",
            fifelse(score=="FFL_6node_all","6-node FFL",
            fifelse(score=="CHANG_WOUND_UP_VANTVEER","wound response (control)",
            fifelse(score=="FARMER_STROMAL","Farmer stroma signature",
            fifelse(score=="HALLMARK_EMT","Hallmark EMT",
            fifelse(score=="FINAK_SDPP","SDPP",
                    gsub("_"," ",score))))))))))]
b <- b[order(-OR_random)]
b[, lab := factor(lab, levels=lab)]
b[, hl := grepl("COL1A1", lab)]
pb <- ggplot(b, aes(OR_random, lab)) +
  geom_vline(xintercept=1, linetype=2, colour=MUTED, linewidth=0.4) +
  geom_errorbarh(aes(xmin=random_lo, xmax=random_hi), height=0, linewidth=0.55, colour=MUTED) +
  geom_point(aes(colour=hl), size=2.2) +
  scale_colour_manual(values=c(`TRUE`=unname(SEM["neg"]), `FALSE`=unname(SEM["null"])), guide="none") +
  scale_x_continuous(trans="log", breaks=c(0.7,0.85,1,1.2,1.5,2,2.5)) +
  scale_y_discrete(labels=ylab) +
  labs(x="odds ratio per SD for pathological complete response (95% CI)", y=NULL, tag="b") +
  theme_pub(BASE) + theme(panel.grid.major.y=element_blank())

fig <- pa / pb + plot_layout(heights=c(1,1))
ggsave(file.path(OUT,"Fig5_reactive_stroma.png"), fig, width=6.85, height=6.4, dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(OUT,"Fig5_reactive_stroma.pdf"), fig, width=6.85, height=6.4, bg="white", device=cairo_pdf)
cat("wrote Fig5_reactive_stroma\n")
