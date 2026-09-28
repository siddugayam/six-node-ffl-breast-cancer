#!/usr/bin/env Rscript
# Consolidated independent-cohort replication table.
set.seed(42)
RES <- "/path/to/revision/results/multiomics/"
R0  <- "/path/to/revision/results/"
DATA<- "/path/to/revision/data/"

tcga_de <- read.csv(paste0(R0,"BRCA_DEX_genes.csv"), stringsAsFactors=FALSE)
hubs <- read.csv(paste0(R0,"network_topology_hubs.csv"), stringsAsFactors=FALSE)
hub_genes <- hubs$name[hubs$is_hub & hubs$type!="miRNA"]
cat("protein-coding hubs:", length(hub_genes), "\n")
rows <- list()

## ---- (i) hub DE direction in the three GEO tumour-vs-normal cohorts ----
geo <- read.csv(paste0(R0,"external_GEO_DE.csv"), stringsAsFactors=FALSE)
plat <- c(GSE42568="Affymetrix HG-U133 Plus 2.0 (GPL570)",
          GSE45827="Affymetrix HG-U133 Plus 2.0 (GPL570)",
          GSE10780="Affymetrix HG-U133 Plus 2.0 (GPL570)")
for (g in unique(geo$gse)) {
  d <- geo[geo$gse==g, ]
  m <- merge(d, tcga_de, by.x="feature", by.y="feature", suffixes=c("_ext","_tcga"))
  m <- m[m$feature %in% hub_genes, ]
  sig <- m[m$adj.P.Val_tcga < 0.05, ]
  same <- sign(sig$logFC_ext)==sign(sig$logFC_tcga)
  samesig <- same & sig$adj.P.Val_ext < 0.05
  bp <- if (nrow(sig)>0) binom.test(sum(same), nrow(sig), 0.5, alternative="greater")$p.value else NA
  rows[[length(rows)+1]] <- data.frame(cohort=g, accession=g, platform=unname(plat[g]),
    n=d$n_normal[1]+d$n_tumour[1], n_detail=paste0(d$n_normal[1]," normal / ",d$n_tumour[1]," tumour"),
    analysis="hub_DE_tumour_vs_normal", n_features_tested=nrow(m), n_TCGA_significant=nrow(sig),
    n_replicate_same_direction=sum(same), pct_replicate_same_direction=round(100*mean(same),1),
    n_replicate_same_direction_and_sig=sum(samesig),
    pct_replicate_same_direction_and_sig=round(100*mean(samesig),1),
    binom_p=signif(bp,3), stringsAsFactors=FALSE)
}
## ---- MET500 metastasis vs TCGA primary (hub direction) ----
met <- read.csv(paste0(R0,"external_MET500_breast.csv"), stringsAsFactors=FALSE)
mh <- met[met$gene %in% hub_genes, ]
rows[[length(rows)+1]] <- data.frame(cohort="MET500", accession="MET500 breast metastases",
  platform="RNA-seq (poly-A / exome capture), log2", n=mh$n_met[1],
  n_detail=paste0(mh$n_met[1]," metastases vs ",mh$n_prim[1]," TCGA primaries"),
  analysis="hub_rank_shift_metastasis_vs_primary", n_features_tested=nrow(mh),
  n_TCGA_significant=NA, n_replicate_same_direction=sum(mh$wilcox_q_BH<0.05),
  pct_replicate_same_direction=round(100*mean(mh$wilcox_q_BH<0.05),1),
  n_replicate_same_direction_and_sig=NA, pct_replicate_same_direction_and_sig=NA,
  binom_p=NA, stringsAsFactors=FALSE)

## ---- (ii) named-axis correlation replication ----
cc <- read.csv(paste0(RES,"external_cohort_correlations.csv"), stringsAsFactors=FALSE)
ref <- cc[cc$cohort=="TCGA-BRCA (reference)", ]
ref$key <- paste(ref$source, ref$target)
info <- read.csv(paste0(RES,"external_cohort_info.csv"), stringsAsFactors=FALSE)
for (co in setdiff(unique(cc$cohort), "TCGA-BRCA (reference)")) {
  d <- cc[cc$cohort==co & !is.na(cc$rho), ]; d$key <- paste(d$source, d$target)
  m <- merge(d, ref[,c("key","rho","p")], by="key", suffixes=c("_ext","_tcga"))
  sig <- m[m$p_tcga < 0.05, ]
  same <- sign(sig$rho_ext)==sign(sig$rho_tcga)
  samesig <- same & sig$p_ext < 0.05
  ii <- info[info$cohort==co, ]
  rows[[length(rows)+1]] <- data.frame(cohort=co, accession=ii$accession, platform=ii$platform,
    n=ii$n_tumours, n_detail=paste0(ii$n_tumours," tumours; ",ii$n_with_miRNA," with miRNA"),
    analysis="named_axis_correlation_direction", n_features_tested=nrow(m),
    n_TCGA_significant=nrow(sig), n_replicate_same_direction=sum(same),
    pct_replicate_same_direction=round(100*mean(same),1),
    n_replicate_same_direction_and_sig=sum(samesig),
    pct_replicate_same_direction_and_sig=round(100*mean(samesig),1),
    binom_p=signif(binom.test(sum(same), nrow(sig), 0.5, alternative="greater")$p.value,3),
    stringsAsFactors=FALSE)
}
## ---- (iii) METABRIC survival ----
sv <- read.csv(paste0(RES,"metabric_survival_replication.csv"), stringsAsFactors=FALSE)
tc <- read.csv(paste0(R0,"survival_cox_hubs.csv"), stringsAsFactors=FALSE)
for (endp in c("OS","RFS")) {
  d <- sv[sv$endpoint==endp & sv$model=="univariate" & !is.na(sv$q_BH), ]
  tce <- tc[tc$endpoint==ifelse(endp=="OS","OS","PFI") & tc$is_hub, ]
  m <- merge(d, tce, by="feature")
  same <- sign(log(m$HR_per_SD.x)) == sign(log(m$HR_per_SD.y))
  rows[[length(rows)+1]] <- data.frame(cohort="METABRIC", accession="cBioPortal brca_metabric",
    platform="Illumina HT-12 v3 microarray", n=max(d$n),
    n_detail=paste0(max(d$n)," tumours, ",max(d$n_event)," ",endp," events"),
    analysis=paste0("hub_Cox_",endp), n_features_tested=nrow(d),
    n_TCGA_significant=sum(m$q_value<0.05),
    n_replicate_same_direction=sum(same), pct_replicate_same_direction=round(100*mean(same),1),
    n_replicate_same_direction_and_sig=sum(d$q_BH<0.05),
    pct_replicate_same_direction_and_sig=round(100*mean(d$q_BH<0.05),1),
    binom_p=signif(binom.test(sum(same), length(same), 0.5, alternative="greater")$p.value,3),
    stringsAsFactors=FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, paste0(RES,"external_cohort_replication.csv"), row.names=FALSE)
cat("\nexternal_cohort_replication.csv rows:", nrow(res), "\n")
print(res[,c("cohort","platform","n","analysis","n_features_tested","n_TCGA_significant",
             "pct_replicate_same_direction","pct_replicate_same_direction_and_sig","binom_p")],
      row.names=FALSE)
