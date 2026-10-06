#!/usr/bin/env Rscript
# v4/06 -- Figures: dot plots, enrichment curves, ridge plots, FFL-class NES
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(ggplot2); library(ggridges)
  library(patchwork); library(cowplot)
})
ROOT <- "/path/to/revision"
RES  <- file.path(ROOT, "results", "v4")
FIG  <- file.path(ROOT, "figures", "v4"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(ROOT, "cache", "v4")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
theme_set(theme_bw(base_size = 9))

sets <- readRDS(file.path(CACHE, "genesets.rds"))
A    <- fread(file.path(RES, "gsea_all_results.csv"))
LISTS <- c("tumour_vs_normal", "mir29a_corr",
           "ffl_signed_influence")
NICE <- c(tumour_vs_normal = "Tumour vs normal (limma t)",
          mir29a_corr = "Gene ~ miR-29a-3p (Spearman rho)",
          ffl_signed_influence = "FFL signed influence")
PRIMARY <- c("H", "C2_CP_KEGG", "C2_CP_REACTOME", "C6")

read_rank <- function(tag) {
  d <- fread(file.path(RES, paste0("rank_", tag, ".csv")))
  d <- d[!is.na(stat)][!duplicated(feature)]
  sort(setNames(d$stat, d$feature), decreasing = TRUE)
}
short <- function(x, n = 52) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)

for (L in LISTS) {
  st <- read_rank(L)
  sub <- A[ranked_list == L & collection %in% PRIMARY & padj < 0.05]
  if (nrow(sub) == 0) { msg("no significant primary-collection terms for ", L); next }

  ## ---------------- dot plot ------------------------------------------------
  up <- head(sub[NES > 0][order(pval)], 15)
  dn <- head(sub[NES < 0][order(pval)], 15)
  dd <- rbind(up, dn)
  dd[, lab := short(paste0(pathway, "  [", collection, "]"), 62)]
  dd[, lab := factor(lab, levels = dd[order(NES)]$lab)]
  p <- ggplot(dd, aes(NES, lab, size = leadingEdge_size, colour = -log10(padj))) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point() +
    scale_colour_viridis_c(option = "C", name = expression(-log[10]~FDR)) +
    scale_size_continuous(name = "leading\nedge n", range = c(1.5, 5)) +
    labs(x = "Normalised enrichment score (NES)", y = NULL,
         title = paste0("GSEA: ", NICE[[L]]),
         subtitle = "Top 15 positive and top 15 negative terms (Hallmark / KEGG / Reactome / C6), FDR < 0.05") +
    theme(axis.text.y = element_text(size = 6.5), plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, paste0("fig_gsea_dotplot_", L, ".png")), p,
         width = 9.2, height = 6.6, dpi = 300)
  ggsave(file.path(FIG, paste0("fig_gsea_dotplot_", L, ".pdf")), p, width = 9.2, height = 6.6)

  ## ---------------- enrichment curves for the top terms ---------------------
  top6 <- rbind(head(sub[NES > 0][order(pval)], 3), head(sub[NES < 0][order(pval)], 3))
  pl <- lapply(seq_len(nrow(top6)), function(i) {
    pw <- top6$pathway[i]; cl <- top6$collection[i]
    g  <- sets[[cl]][[pw]]
    plotEnrichment(g, st) +
      labs(title = short(pw, 46),
           subtitle = sprintf("NES %.2f | FDR %.2g | size %d",
                              top6$NES[i], top6$padj[i], top6$size[i]),
           x = "rank in ordered gene list", y = "enrichment score") +
      theme(plot.title = element_text(size = 7.5, face = "bold"),
            plot.subtitle = element_text(size = 6.5))
  })
  pp <- wrap_plots(pl, ncol = 3) +
    plot_annotation(title = paste0("Enrichment plots — ", NICE[[L]]),
                    theme = theme(plot.title = element_text(face = "bold", size = 11)))
  ggsave(file.path(FIG, paste0("fig_gsea_enrichment_", L, ".png")), pp,
         width = 10.5, height = 5.6, dpi = 300)
  ggsave(file.path(FIG, paste0("fig_gsea_enrichment_", L, ".pdf")), pp, width = 10.5, height = 5.6)

  ## ---------------- ridge plot ---------------------------------------------
  rid <- rbind(head(sub[NES > 0][order(pval)], 10), head(sub[NES < 0][order(pval)], 10))
  rl <- rbindlist(lapply(seq_len(nrow(rid)), function(i) {
    le <- strsplit(rid$leadingEdge[i], ";", fixed = TRUE)[[1]]
    data.table(pathway = short(rid$pathway[i], 46), NES = rid$NES[i],
               value = as.numeric(st[le]))
  }))
  rl[, pathway := factor(pathway, levels = unique(rl[order(NES)]$pathway))]
  pr <- ggplot(rl, aes(x = value, y = pathway, fill = NES)) +
    ggridges::geom_density_ridges(scale = 2.1, alpha = 0.85, linewidth = 0.25) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0) +
    labs(x = paste0("ranking statistic (", NICE[[L]], ")"), y = NULL,
         title = paste0("Leading-edge distributions — ", NICE[[L]]),
         subtitle = "Top 10 positive and top 10 negative terms; density of the ranking statistic over leading-edge genes") +
    theme(axis.text.y = element_text(size = 6.5), plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, paste0("fig_gsea_ridge_", L, ".png")), pr,
         width = 8.6, height = 6.4, dpi = 300)
  ggsave(file.path(FIG, paste0("fig_gsea_ridge_", L, ".pdf")), pr, width = 8.6, height = 6.4)
  msg("figures written for ", L)
}

## ================= FFL class NES ============================================
fc <- A[collection == "FFL_CLASS" & !grepl("__|exemplar|ALL_NETWORK", pathway)]
fc[, pathway := factor(pathway, levels = c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node"))]
fc[, ranked_list := factor(ranked_list, levels = LISTS, labels = NICE[LISTS])]
fc[, sig := ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "ns")))]
pf <- ggplot(fc, aes(pathway, NES, fill = NES)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sig, vjust = ifelse(NES > 0, -0.3, 1.2)), size = 2.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~ ranked_list, scales = "free_y", ncol = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0) +
  labs(x = "FFL motif class (protein-coding node set)", y = "NES",
       title = "Are the FFL motif-class node sets enriched in each ranking?",
       subtitle = "fgsea, 10,000 nPermSimple, BH within collection. * FDR<0.05  ** <0.01  *** <0.001") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        strip.text = element_text(size = 7), legend.position = "none",
        plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_gsea_ffl_class_NES.png"), pf, width = 9.5, height = 5.4, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_ffl_class_NES.pdf"), pf, width = 9.5, height = 5.4)

## FFL class enrichment curves against tumour-vs-normal
st1 <- read_rank("tumour_vs_normal")
fc1 <- A[collection == "FFL_CLASS" & ranked_list == "tumour_vs_normal" &
           pathway %in% c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node")]
pl <- lapply(seq_len(nrow(fc1)), function(i) {
  plotEnrichment(sets$FFL_CLASS[[fc1$pathway[i]]], st1) +
    labs(title = fc1$pathway[i],
         subtitle = sprintf("NES %.2f | FDR %.2g | n=%d", fc1$NES[i], fc1$padj[i], fc1$size[i]),
         x = "rank (tumour vs normal t)", y = "ES") +
    theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5))
})
pp <- wrap_plots(pl, ncol = 3) +
  plot_annotation(title = "FFL motif-class node sets against the tumour-vs-normal ranking",
                  theme = theme(plot.title = element_text(face = "bold", size = 11)))
ggsave(file.path(FIG, "fig_gsea_ffl_class_enrichment.png"), pp, width = 10.5, height = 5.6, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_ffl_class_enrichment.pdf"), pp, width = 10.5, height = 5.6)
msg("FFL class figures written")

## ================= Hallmark NES heatmap across lists ========================
hh <- A[collection == "H"]
hm <- dcast(hh, pathway ~ ranked_list, value.var = "NES")
hp <- dcast(hh, pathway ~ ranked_list, value.var = "padj")
mh <- melt(hm, id.vars = "pathway", variable.name = "ranked_list", value.name = "NES")
mp <- melt(hp, id.vars = "pathway", variable.name = "ranked_list", value.name = "padj")
mh <- merge(mh, mp, by = c("pathway", "ranked_list"))
ord <- hm[order(-hm[["tumour_vs_normal"]])]$pathway
mh[, pathway := factor(sub("^HALLMARK_", "", pathway), levels = sub("^HALLMARK_", "", ord))]
mh[, ranked_list := factor(ranked_list, levels = LISTS, labels = NICE[LISTS])]
ph <- ggplot(mh, aes(ranked_list, pathway, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(!is.na(padj) & padj < 0.05, "*", "")), size = 2.4, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey93") +
  labs(x = NULL, y = NULL, title = "Hallmark NES across all three ranked lists",
       subtitle = "* FDR < 0.05 (BH within collection)") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7),
        axis.text.y = element_text(size = 6), plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_gsea_hallmark_NES_heatmap.png"), ph, width = 7.6, height = 8.2, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_hallmark_NES_heatmap.pdf"), ph, width = 7.6, height = 8.2)

## ================= leading-edge recurrence ==================================
lf <- fread(file.path(RES, "gsea_leadingedge_gene_frequency.csv"))
tp <- head(lf[order(-n_ranked_lists, -n_sets_total)], 40)
tp[, gene := factor(gene, levels = rev(tp$gene))]
pg <- ggplot(tp, aes(n_sets_total, gene, fill = factor(n_ranked_lists))) +
  geom_col(width = 0.75) +
  scale_fill_brewer(palette = "YlOrRd", name = "ranked lists\ndriven") +
  labs(x = "number of significant gene sets whose leading edge contains this gene",
       y = NULL, title = "Genes driving the most signatures across all three ranked lists",
       subtitle = "Leading-edge recurrence (FDR<0.05 sets, MSigDB collections only)") +
  theme(axis.text.y = element_text(size = 6.5), plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_gsea_leadingedge_recurrence.png"), pg, width = 8, height = 7.2, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_leadingedge_recurrence.pdf"), pg, width = 8, height = 7.2)
msg("DONE 06")
