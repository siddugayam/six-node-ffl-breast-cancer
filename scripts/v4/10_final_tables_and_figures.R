#!/usr/bin/env Rscript
# v4/10 -- (D) one consolidated tidy results table across every analysis in v4,
#              plus the two figures that carry the negative/positive controls.
suppressPackageStartupMessages({
  library(data.table); library(fgsea); library(ggplot2); library(patchwork)
})
ROOT  <- "/path/to/revision"
RES   <- file.path(ROOT, "results", "v4")
FIG   <- file.path(ROOT, "figures", "v4"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(ROOT, "cache", "v4")
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
theme_set(theme_bw(base_size = 9))

sets <- readRDS(file.path(CACHE, "genesets.rds"))

## ============================ consolidated tidy table =======================
mt  <- fread(file.path(RES, "gsea_master_table.csv"))       # 5 primary lists x 14 collections + FFL derived
dg  <- fread(file.path(RES, "gsea_mir130a_control_diagnostic_all.csv"))
prov<- fread(file.path(RES, "ranklist_provenance.csv"))

DIAG_DETAIL <- c(
  mir130a_corr_adj_epi = "1066 paired primary tumours; partial Spearman | epithelial score",
  mir130a_corr_adj_all = "1066 paired primary tumours; partial Spearman | epithelial + CAF + immune scores",
  mir29a_corr_adj_all  = "1066 paired primary tumours; partial Spearman | epithelial + CAF + immune scores",
  mir130a_corr_luma    = "414 PAM50 LumA tumours; Spearman rho vs hsa-miR-130a-3p",
  mir130a_corr_lumb    = "182 PAM50 LumB tumours; Spearman rho vs hsa-miR-130a-3p",
  mir130a_corr_basal   = "132 PAM50 Basal tumours; Spearman rho vs hsa-miR-130a-3p")
dg2 <- dg[!ranked_list %in% c("mir130a_corr", "mir29a_corr")]
dg2[, `:=`(analysis = "control diagnostic (sensitivity rankings)",
           collection_label = c(C3_MIR_LEGACY = "miRNA targets (legacy seed)",
                                C3_MIR_MIRDB = "miRNA targets (miRDB)",
                                H = "Hallmark")[collection],
           ranking_statistic = "Spearman / partial-Spearman rho vs the miRNA",
           ranking_detail = DIAG_DETAIL[ranked_list],
           padj_global = NA_real_, n_features_ranked = NA_integer_, n_samples = NA_integer_,
           minSize_used = 10L, maxSize_used = 500L, nPermSimple = 10000L,
           method = "fgseaMultilevel(eps=0)",
           significant_FDR05 = padj < 0.05,
           direction = ifelse(NES > 0, "enriched at the positive/high end",
                              "enriched at the negative/low end"))]
mt[, analysis := ifelse(collection == "FFL_DERIVED", "FFL derived node sets",
                 ifelse(collection == "FFL_CLASS", "FFL motif-class node sets (analysis C)",
                        "primary GSEA (analyses 1-5)"))]
ALL <- rbind(mt, dg2[, names(mt), with = FALSE], use.names = TRUE)
setcolorder(ALL, c("analysis", "ranked_list", "ranking_statistic", "ranking_detail",
                   "n_features_ranked", "n_samples", "collection", "collection_label",
                   "pathway", "size", "ES", "NES", "pval", "padj", "padj_global",
                   "log2err", "direction", "significant_FDR05",
                   "leadingEdge_size", "leadingEdge"))
setorder(ALL, analysis, ranked_list, collection, pval)
fwrite(ALL, file.path(RES, "gsea_MASTER_TABLE_ALL.csv"))
msg("gsea_MASTER_TABLE_ALL.csv: ", nrow(ALL), " rows | FDR<0.05: ", sum(ALL$significant_FDR05))

## a compact, human-readable headline table (no leading-edge strings)
HEAD <- ALL[significant_FDR05 == TRUE,
            .(analysis, ranked_list, collection, collection_label, pathway, size,
              NES, pval, padj, leadingEdge_size)]
setorder(HEAD, ranked_list, collection, pval)
top <- rbind(HEAD[NES > 0, head(.SD, 15), by = .(ranked_list, collection)],
             HEAD[NES < 0, head(.SD, 15), by = .(ranked_list, collection)])
setorder(top, ranked_list, collection, -NES)
fwrite(top, file.path(RES, "gsea_MASTER_TABLE_TOP15.csv"))
msg("gsea_MASTER_TABLE_TOP15.csv: ", nrow(top), " rows")

## per ranked-list x collection counts, including the diagnostics
cnt <- ALL[, .(n_tested = .N, n_sig_FDR05 = sum(significant_FDR05),
               n_sig_pos = sum(significant_FDR05 & NES > 0),
               n_sig_neg = sum(significant_FDR05 & NES < 0),
               top_pos = pathway[which.max(ifelse(significant_FDR05, NES, NA))],
               top_neg = pathway[which.min(ifelse(significant_FDR05, NES, NA))]),
           by = .(analysis, ranked_list, collection, collection_label)]
setorder(cnt, analysis, ranked_list, collection)
fwrite(cnt, file.path(RES, "gsea_MASTER_COUNTS.csv"))
print(cnt[analysis == "primary GSEA (analyses 1-5)"], nrows = 80)

## ============================ figure: the two controls side by side =========
read_rank <- function(tag) {
  d <- fread(file.path(RES, paste0("rank_", tag, ".csv")))
  d <- d[!is.na(stat)][!duplicated(feature)]
  sort(setNames(d$stat, d$feature), decreasing = TRUE)
}
panel <- function(cl, pw, rl, title, src) {
  st  <- read_rank(rl)
  row <- src[ranked_list == rl & collection == cl & pathway == pw]
  plotEnrichment(sets[[cl]][[pw]], st) +
    labs(title = title,
         subtitle = if (nrow(row)) sprintf("NES %+.2f | p = %.3g | FDR = %.3g | n = %d",
                                           row$NES[1], row$pval[1], row$padj[1], row$size[1])
                    else "not tested",
         x = "rank in ordered list", y = "enrichment score") +
    theme(plot.title = element_text(size = 7.6, face = "bold"),
          plot.subtitle = element_text(size = 6.6))
}
A  <- fread(file.path(RES, "gsea_all_results.csv"))
D  <- fread(file.path(RES, "gsea_mir130a_control_diagnostic_all.csv"))
pl <- list(
  panel("C3_MIR_MIRDB", "MIR130A_3P", "mir130a_corr",
        "FAILS  •  miR-130a-3p targets (miRDB)\nvs gene ~ miR-130a-3p rho", A),
  panel("C3_MIR_LEGACY", "TTGCACT_MIR130A_MIR301_MIR130B", "mir130a_corr",
        "FAILS  •  miR-130a seed family (legacy)\nvs gene ~ miR-130a-3p rho", A),
  panel("C3_MIR_MIRDB", "MIR130A_3P", "mir130a_high_vs_low",
        "FAILS  •  miR-130a-3p targets (miRDB)\nvs miR-130a high-vs-low tertile t", A),
  panel("C3_MIR_MIRDB", "MIR29A_3P", "mir29a_corr",
        "PASSES  •  miR-29a-3p targets (miRDB)\nvs gene ~ miR-29a-3p rho", A),
  panel("C3_MIR_LEGACY", "TGGTGCT_MIR29A_MIR29B_MIR29C", "mir29a_corr",
        "PASSES  •  miR-29a seed family (legacy)\nvs gene ~ miR-29a-3p rho", A),
  panel("C3_MIR_MIRDB", "MIR130A_3P", "mir130a_corr_adj_all",
        "STILL FAILS  •  miR-130a-3p targets (miRDB)\nvs composition-adjusted miR-130a rho", D))
pp <- wrap_plots(pl, ncol = 3) +
  plot_annotation(
    title = "Internal positive control: do a miRNA's own predicted targets sit at the low end of its expression axis?",
    subtitle = "A working axis gives a NEGATIVE NES. miR-29a-3p passes on the identical construction; miR-130a-3p does not, and does not recover after adjusting for epithelial/CAF/immune composition.",
    theme = theme(plot.title = element_text(face = "bold", size = 10.5),
                  plot.subtitle = element_text(size = 7.6)))
ggsave(file.path(FIG, "fig_gsea_mir_control_130a_vs_29a.png"), pp, width = 11.5, height = 6.4, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_mir_control_130a_vs_29a.pdf"), pp, width = 11.5, height = 6.4)
msg("control comparison figure written")

## ============================ figure: FFL class vs both backgrounds =========
nb <- fread(file.path(RES, "gsea_ffl_class_network_background_test.csv"))
tw <- fread(file.path(RES, "gsea_ffl_class_transcriptome_background_test.csv"))
KEEP <- c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node",
          "4node_not_3node", "5node_not_3node", "6node_not_3node",
          "higherorder_not_3node", "ALL_NETWORK_NODES")
b <- merge(nb[set %in% KEEP, .(set, n_in_bg, network_bg = bg_mean_absT, p_network = p_vs_network_bg)],
           tw[set %in% KEEP, .(set, n, obs_mean_absT, transcriptome_bg = transcriptome_mean_absT,
                               p_transcriptome = p_vs_transcriptome)],
           by = "set")
## ALL_NETWORK_NODES *is* the within-network background, so it has no network-background test
b[set == "ALL_NETWORK_NODES", `:=`(network_bg = NA_real_, p_network = NA_real_)]
b[, set := factor(set, levels = rev(KEEP))]
bl <- melt(b, id.vars = c("set", "obs_mean_absT"),
           measure.vars = c("network_bg", "transcriptome_bg"),
           variable.name = "background", value.name = "bg_mean")
bl <- bl[!is.na(bg_mean)]
bl[, background := factor(background, levels = c("transcriptome_bg", "network_bg"),
                          labels = c("size-matched draw from the transcriptome (20,250 genes)",
                                     "size-matched draw from the 364 network nodes"))]
pb <- ggplot(bl, aes(y = set)) +
  geom_segment(aes(x = bg_mean, xend = obs_mean_absT, yend = set), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(x = bg_mean), colour = "grey45", size = 1.7, shape = 1) +
  geom_point(aes(x = obs_mean_absT), colour = "#B2182B", size = 2.3) +
  facet_wrap(~ background, ncol = 2) +
  labs(x = expression("mean |t| (tumour vs normal)"), y = NULL,
       title = "Are the FFL motif-class node sets more disease-associated than chance?",
       subtitle = "Red = observed set mean |t|; open circle = mean of 10,000 size-matched random draws.\nEvery class beats the transcriptome; only the 3-node classes beat a within-network background.") +
  theme(axis.text.y = element_text(size = 7), strip.text = element_text(size = 7.2),
        plot.title = element_text(face = "bold"))
ggsave(file.path(FIG, "fig_gsea_ffl_class_backgrounds.png"), pb, width = 9.4, height = 4.6, dpi = 300)
ggsave(file.path(FIG, "fig_gsea_ffl_class_backgrounds.pdf"), pb, width = 9.4, height = 4.6)
msg("FFL background figure written")
fwrite(b, file.path(RES, "gsea_ffl_class_background_summary.csv"))

## ============================ file inventory ================================
fs <- data.table(file = list.files(RES, pattern = "\\.(csv|rds)$", full.names = FALSE))
fs[, size_kb := round(file.size(file.path(RES, file)) / 1024, 1)]
fs[, n_rows := vapply(file, function(f) if (grepl("\\.csv$", f))
    tryCatch(nrow(fread(file.path(RES, f), select = 1L)), error = function(e) NA_integer_)
    else NA_integer_, integer(1))]
setorder(fs, file)
fwrite(fs, file.path(RES, "gsea_OUTPUT_INVENTORY.csv"))
print(fs, nrows = 100)
msg("DONE 10")
