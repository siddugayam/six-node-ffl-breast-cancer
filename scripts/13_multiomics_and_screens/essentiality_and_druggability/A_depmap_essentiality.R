#!/usr/bin/env Rscript
## A) DepMap CRISPR essentiality of miRNA-TF-gene FFL network nodes
suppressPackageStartupMessages({library(data.table)})
set.seed(1)
REV  <- "/path/to/revision"
OUT  <- file.path(REV,"results/multiomics")
RAW  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
MODEL<- "/path/to/home/Desktop/DD/R_GPR/ESIA_BRCA/brca_eisa_pilot/depmap/Model.csv"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"), "|", ..., "\n")

## ---- annotation ----------------------------------------------------------
mod <- fread(MODEL)
msg("Model.csv rows:", nrow(mod))
breast_ids <- mod[OncotreeLineage=="Breast", ModelID]
msg("Breast models in annotation:", length(breast_ids))

## ---- CRISPR gene effect --------------------------------------------------
msg("reading CRISPRGeneEffect.csv ...")
ce <- fread(file.path(RAW,"CRISPRGeneEffect.csv"))
setnames(ce, 1, "ModelID")
msg("CRISPR matrix:", nrow(ce), "cell lines x", ncol(ce)-1, "genes")

gcols <- setdiff(names(ce),"ModelID")
gsym  <- sub(" \\(\\d+\\)$","", gcols)
stopifnot(!any(duplicated(gsym)))

ce_breast_ids <- intersect(ce$ModelID, breast_ids)
msg("Breast lines with CRISPR data:", length(ce_breast_ids))
is_breast <- ce$ModelID %in% ce_breast_ids
msg("Non-breast lines with CRISPR data:", sum(!is_breast))

## lineage for every screened line (for reporting)
lin <- mod[match(ce$ModelID, ModelID), OncotreeLineage]
msg("Screened lines with no lineage annotation:", sum(is.na(lin)))

## ---- network nodes -------------------------------------------------------
nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
topo  <- fread(file.path(REV,"results/network_topology_hubs.csv"))
mods  <- readRDS(file.path(REV,"data/ffl_module_sets.rds"))
prot  <- nodes[type %in% c("TF","Gene")]
msg("protein-coding network nodes (TF+Gene):", nrow(prot))

in_mat <- prot$name %in% gsym
msg("network protein nodes present in DepMap matrix:", sum(in_mat))
missing <- prot$name[!in_mat]
writeLines(missing, file.path(OUT,"depmap_missing_network_genes.txt"))
msg("missing (written to file):", paste(missing, collapse=", "))

sel <- prot$name[in_mat]
idx <- match(sel, gsym)
M   <- as.matrix(ce[, gcols[idx], with=FALSE])   # lines x selected genes
colnames(M) <- sel

## ---- per-gene statistics -------------------------------------------------
Mb <- M[is_breast, , drop=FALSE]
Mo <- M[!is_breast, , drop=FALSE]

res <- data.table(gene = sel)
res[, n_breast_screened   := colSums(!is.na(Mb))]
res[, n_other_screened    := colSums(!is.na(Mo))]
res[, mean_chronos_breast := colMeans(Mb, na.rm=TRUE)]
res[, median_chronos_breast := apply(Mb,2,median,na.rm=TRUE)]
res[, sd_chronos_breast   := apply(Mb,2,sd,na.rm=TRUE)]
res[, frac_essential_breast := colMeans(Mb < -0.5, na.rm=TRUE)]
res[, n_essential_breast  := colSums(Mb < -0.5, na.rm=TRUE)]
res[, mean_chronos_other  := colMeans(Mo, na.rm=TRUE)]
res[, frac_essential_other:= colMeans(Mo < -0.5, na.rm=TRUE)]
res[, delta_breast_minus_other := mean_chronos_breast - mean_chronos_other]

## Wilcoxon rank-sum: breast vs all other lineages (selective essentiality)
wp <- vapply(seq_len(ncol(M)), function(j){
  b <- Mb[,j]; o <- Mo[,j]
  b <- b[!is.na(b)]; o <- o[!is.na(o)]
  if(length(b)<3 || length(o)<3) return(NA_real_)
  suppressWarnings(wilcox.test(b,o)$p.value)
}, numeric(1))
res[, wilcox_p_breast_vs_other := wp]
res[, wilcox_q_breast_vs_other := p.adjust(wp, method="BH")]
res[, selective_breast := !is.na(wilcox_q_breast_vs_other) &
      wilcox_q_breast_vs_other < 0.05 & delta_breast_minus_other < 0]

## common-essential flag (pan-lineage)
res[, mean_chronos_all := colMeans(M, na.rm=TRUE)]
res[, frac_essential_all := colMeans(M < -0.5, na.rm=TRUE)]
res[, essential_class := fifelse(frac_essential_all >= 0.9, "common_essential",
                          fifelse(frac_essential_breast >= 0.5, "essential_in_most_breast",
                          fifelse(frac_essential_breast > 0, "essential_in_some_breast",
                                  "non_essential")))]

## ---- annotate with network attributes ------------------------------------
res <- merge(res, topo[, .(gene=name, node_type=type, degree_tot, betweenness,
                           hub_by_degree, hub_by_betweenness, hub_named_in_MS, is_hub)],
             by="gene", all.x=TRUE)
res[, in_N3          := gene %in% mods$N3]
res[, in_N4          := gene %in% mods$N4]
res[, in_N5          := gene %in% mods$N5]
res[, in_N6          := gene %in% mods$N6]
res[, higher_order_only := gene %in% mods$HIGHER_ONLY]
res[, in_MS_module   := gene %in% mods$MS_MODULE]
res[, module_class   := fifelse(higher_order_only, "higher_order_only",
                         fifelse(in_N3, "in_3node_FFL", "not_in_any_FFL"))]
## DE status
de <- fread(file.path(REV,"results/BRCA_DEX_genes.csv"))
denm <- names(de); fcol <- denm[1]
res <- merge(res, de[, .(gene=get(fcol), logFC_TvsN=logFC, adjP_TvsN=adj.P.Val)],
             by="gene", all.x=TRUE)

setorder(res, mean_chronos_breast)
fwrite(res, file.path(OUT,"depmap_essentiality_nodes.csv"))
msg("wrote depmap_essentiality_nodes.csv rows:", nrow(res))

## ---- named hubs ----------------------------------------------------------
named <- c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","VEGFA","CCND2","MYC",
           "E2F1","TP53","EZH2","STAT3","HIF1A","TGFBR2")
nm <- res[gene %in% named]
setorder(nm, mean_chronos_breast)
fwrite(nm, file.path(OUT,"depmap_named_hubs.csv"))
msg("named hubs found:", nrow(nm), "of", length(named),
    "| missing:", paste(setdiff(named, res$gene), collapse=","))
print(nm[, .(gene, node_type, n_breast_screened, mean_chronos_breast,
             frac_essential_breast, mean_chronos_other, delta_breast_minus_other,
             wilcox_q_breast_vs_other, selective_breast, essential_class, is_hub)])

## ---- comparison 1: hubs vs non-hubs --------------------------------------
cmp <- list()
h  <- res[is_hub==TRUE,  mean_chronos_breast]
nh <- res[is_hub==FALSE, mean_chronos_breast]
w1 <- wilcox.test(h, nh)
cmp[[1]] <- data.table(comparison="hub vs non-hub (network protein nodes)",
  group1="hub", n1=length(h), median1=median(h), mean1=mean(h),
  group2="non-hub", n2=length(nh), median2=median(nh), mean2=mean(nh),
  statistic=unname(w1$statistic), p_value=w1$p.value, test="Wilcoxon rank-sum (two-sided)")

## also fraction-essential comparison
h2  <- res[is_hub==TRUE,  frac_essential_breast]
nh2 <- res[is_hub==FALSE, frac_essential_breast]
w1b <- wilcox.test(h2, nh2)
cmp[[2]] <- data.table(comparison="hub vs non-hub (fraction breast lines essential)",
  group1="hub", n1=length(h2), median1=median(h2), mean1=mean(h2),
  group2="non-hub", n2=length(nh2), median2=median(nh2), mean2=mean(nh2),
  statistic=unname(w1b$statistic), p_value=w1b$p.value, test="Wilcoxon rank-sum (two-sided)")

## ---- comparison 2: higher-order vs 3-node only ---------------------------
ho <- res[higher_order_only==TRUE, mean_chronos_breast]
n3 <- res[in_N3==TRUE & higher_order_only==FALSE, mean_chronos_breast]
w2 <- wilcox.test(ho, n3)
cmp[[3]] <- data.table(comparison="higher-order-only nodes vs 3-node-FFL nodes (mean Chronos, breast)",
  group1="higher_order_only", n1=length(ho), median1=median(ho), mean1=mean(ho),
  group2="in_3node_FFL", n2=length(n3), median2=median(n3), mean2=mean(n3),
  statistic=unname(w2$statistic), p_value=w2$p.value, test="Wilcoxon rank-sum (two-sided)")

ho2 <- res[higher_order_only==TRUE, frac_essential_breast]
n32 <- res[in_N3==TRUE & higher_order_only==FALSE, frac_essential_breast]
w2b <- wilcox.test(ho2, n32)
cmp[[4]] <- data.table(comparison="higher-order-only vs 3-node-FFL (fraction breast lines essential)",
  group1="higher_order_only", n1=length(ho2), median1=median(ho2), mean1=mean(ho2),
  group2="in_3node_FFL", n2=length(n32), median2=median(n32), mean2=mean(n32),
  statistic=unname(w2b$statistic), p_value=w2b$p.value, test="Wilcoxon rank-sum (two-sided)")

## N6-member vs N3-member (overlapping sets, stated as such)
a <- res[in_N6==TRUE, mean_chronos_breast]; b <- res[in_N3==TRUE, mean_chronos_breast]
w3 <- wilcox.test(a,b)
cmp[[5]] <- data.table(comparison="members of N6 modules vs members of N3 modules (overlapping sets)",
  group1="N6", n1=length(a), median1=median(a), mean1=mean(a),
  group2="N3", n2=length(b), median2=median(b), mean2=mean(b),
  statistic=unname(w3$statistic), p_value=w3$p.value, test="Wilcoxon rank-sum (two-sided)")

## MS 10-node module vs rest of network
a <- res[in_MS_module==TRUE, mean_chronos_breast]; b <- res[in_MS_module==FALSE, mean_chronos_breast]
w4 <- wilcox.test(a,b)
cmp[[6]] <- data.table(comparison="manuscript 10-node module vs all other network protein nodes",
  group1="MS_module", n1=length(a), median1=median(a), mean1=mean(a),
  group2="other", n2=length(b), median2=median(b), mean2=mean(b),
  statistic=unname(w4$statistic), p_value=w4$p.value, test="Wilcoxon rank-sum (two-sided)")

## network nodes vs genome-wide background (are network genes essential at all?)
gw <- colMeans(M_all <- as.matrix(ce[is_breast, gcols, with=FALSE]), na.rm=TRUE)
netmean <- res$mean_chronos_breast
bg <- gw[!(names(gw) %in% paste0(res$gene," "))]
bgv <- gw[ !(sub(" \\(\\d+\\)$","",names(gw)) %in% res$gene) ]
w5 <- wilcox.test(netmean, bgv)
cmp[[7]] <- data.table(comparison="network protein nodes vs all other genes genome-wide (breast lines)",
  group1="network", n1=length(netmean), median1=median(netmean), mean1=mean(netmean),
  group2="genome_background", n2=length(bgv), median2=median(bgv,na.rm=TRUE), mean2=mean(bgv,na.rm=TRUE),
  statistic=unname(w5$statistic), p_value=w5$p.value, test="Wilcoxon rank-sum (two-sided)")

cmpdt <- rbindlist(cmp)
fwrite(cmpdt, file.path(OUT,"depmap_group_comparisons.csv"))
msg("wrote depmap_group_comparisons.csv rows:", nrow(cmpdt))
print(cmpdt)

## ---- lineage composition summary ----------------------------------------
lsum <- data.table(lineage=lin)[, .N, by=lineage][order(-N)]
fwrite(lsum, file.path(OUT,"depmap_screened_lineages.csv"))
msg("wrote depmap_screened_lineages.csv rows:", nrow(lsum))

## breast line list
bl <- mod[ModelID %in% ce_breast_ids, .(ModelID, StrippedCellLineName, OncotreeSubtype,
                                        OncotreeCode, PrimaryOrMetastasis)]
fwrite(bl, file.path(OUT,"depmap_breast_lines_used.csv"))
msg("wrote depmap_breast_lines_used.csv rows:", nrow(bl))

## top essential network genes in breast
setorder(res, mean_chronos_breast)
msg("=== 20 most essential network nodes in breast lines ===")
print(res[1:20, .(gene,node_type,mean_chronos_breast,frac_essential_breast,essential_class,is_hub)])
msg("DONE")
