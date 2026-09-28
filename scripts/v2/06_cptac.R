## Task E: CPTAC-BRCA replication of the miR-29 -> collagen axis at PROTEIN level,
## and global miRNA-target concordance mRNA vs protein.  (Verifies N12.)
setwd("/path/to/revision")
BASE <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
PRO <- file.path(BASE,"Proteome_BCM_GENCODE_v34_harmonized_v1/Proteome_BCM_GENCODE_v34_harmonized_v1",
                 "BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
RNA <- file.path(BASE,"RNA_BCM_v1/RNA_BCM_v1","BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt")
MIR <- file.path(BASE,"miRNA_BCM_v1","BRCA_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt")
ANN <- file.path(BASE,"Proteome_BCM_GENCODE_v34_harmonized_v1/Proteome_BCM_GENCODE_v34_harmonized_v1",
                 "README/Gene_annotation_and_representable_isoform_mapping_table.txt")
rd <- function(f){ x <- read.delim(f, row.names=1, check.names=FALSE); as.matrix(x) }
P <- rd(PRO); R <- rd(RNA); MI <- rd(MIR)
cat("protein:", dim(P), " rna:", dim(R), " mirna:", dim(MI), "\n")
ann <- read.delim(ANN, stringsAsFactors=FALSE)
map <- setNames(ann$gene_name, ann$gene)                       # versioned ENSG -> symbol
sym <- function(M){ s <- map[rownames(M)]
  bad <- is.na(s); if(any(bad)){ s2 <- setNames(ann$gene_name, sub("\\..*","",ann$gene))
    s[bad] <- s2[sub("\\..*","",rownames(M)[bad])] }
  M <- M[!is.na(s) & !duplicated(s), , drop=FALSE]
  rownames(M) <- s[!is.na(s) & !duplicated(s)]; M }
P <- sym(P); R <- sym(R)
cat("after symbol mapping -> protein:", dim(P), " rna:", dim(R), "\n")
cs_pm <- intersect(colnames(P), colnames(MI)); cs_rm <- intersect(colnames(R), colnames(MI))
cs_all <- Reduce(intersect, list(colnames(P), colnames(R), colnames(MI)))
cat("samples: protein&miRNA =", length(cs_pm), " rna&miRNA =", length(cs_rm),
    " all three =", length(cs_all), "\n")

## ---- map network miRNA labels onto CPTAC mature names (arm chosen by abundance)
armmap <- function(nm){
  if (nm %in% rownames(MI)) return(nm)
  cand <- c(paste0(nm,"-3p"), paste0(nm,"-5p")); cand <- cand[cand %in% rownames(MI)]
  if (!length(cand)) return(NA_character_)
  cand[which.max(rowMeans(MI[cand,,drop=FALSE]))]
}
mirs <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
for (m in mirs) cat(m, "->", armmap(m), "\n")

sp <- function(a,b) suppressWarnings(cor(a,b,method="spearman"))
sp_p <- function(a,b){ k <- complete.cases(a,b); if(sum(k)<10) return(c(NA,NA,NA))
  ct <- suppressWarnings(cor.test(a[k],b[k],method="spearman", exact=FALSE))
  c(unname(ct$estimate), ct$p.value, sum(k)) }
## Steiger test for two dependent correlations sharing one variable
steiger <- function(r12,r13,r23,n){
  R <- (1-r12^2-r13^2-r23^2)+2*r12*r13*r23
  t2 <- (r12-r13)*sqrt(((n-1)*(1+r23))/((2*((n-1)/(n-3))*R)+(((r12+r13)/2)^2)*(1-r23)^3))
  c(t=t2, p=2*pt(-abs(t2), n-3))
}

targets <- c("COL1A1","COL3A1","COL1A2","COL4A1","COL5A1","COL5A2","COL6A3","COL11A1","FN1","SPARC","LOX")
targets <- intersect(targets, intersect(rownames(P), rownames(R)))
out <- NULL
for (m in mirs){ mm <- armmap(m); if(is.na(mm)) next
  for (g in targets){
    a <- sp_p(MI[mm, cs_rm], R[g, cs_rm])      # mRNA level
    b <- sp_p(MI[mm, cs_pm], P[g, cs_pm])      # protein level
    ss <- cs_all
    r12 <- sp(MI[mm,ss], R[g,ss]); r13 <- sp(MI[mm,ss], P[g,ss]); r23 <- sp(R[g,ss], P[g,ss])
    st <- steiger(r12,r13,r23,length(ss))
    out <- rbind(out, data.frame(miRNA=m, arm=mm, target=g,
      rho_mRNA=a[1], p_mRNA=a[2], n_mRNA=a[3],
      rho_protein=b[1], p_protein=b[2], n_protein=b[3],
      rho_mRNA_common=r12, rho_prot_common=r13, rho_mRNA_prot=r23, n_common=length(ss),
      steiger_t=unname(st[1]), steiger_p=unname(st[2]),
      protein_stronger = abs(b[1]) > abs(a[1]))) }
}
print(out[,c("miRNA","arm","target","rho_mRNA","p_mRNA","n_mRNA","rho_protein","p_protein",
             "n_protein","steiger_p","protein_stronger")], digits=3, row.names=FALSE)
cat("\nmiR-29 x ECM pairs where |rho_protein| > |rho_mRNA|:",
    sum(out$protein_stronger), "of", nrow(out),
    " binom p =", signif(binom.test(sum(out$protein_stronger), nrow(out), 0.5)$p.value,3), "\n")
core <- out[out$target %in% c("COL1A1","COL3A1","COL1A2","COL5A1","COL5A2"),]
cat("core fibrillar collagens only:", sum(core$protein_stronger), "of", nrow(core),
    " binom p =", signif(binom.test(sum(core$protein_stronger), nrow(core), 0.5)$p.value,3), "\n")

## ---- variant: collapse a network miRNA label by AVERAGING all mature arms ----
armavg <- function(nm, cols){ a <- grep(paste0("^", nm, "($|-)"), rownames(MI), value=TRUE)
  a <- a[grepl(paste0("^", nm, "(-[0-9]+)?-(3p|5p)$|^", nm, "$"), a)]
  if(!length(a)) return(NULL); list(v=colMeans(MI[a, cols, drop=FALSE]), arms=a) }
v2 <- NULL
for (m in mirs) for (g in c("COL1A1","COL3A1")){
  aR <- armavg(m, cs_rm); aP <- armavg(m, cs_pm); if(is.null(aR)) next
  x <- sp_p(aR$v, R[g, cs_rm]); y <- sp_p(aP$v, P[g, cs_pm])
  v2 <- rbind(v2, data.frame(miRNA=m, arms=paste(aR$arms, collapse="+"), target=g,
    rho_mRNA=x[1], p_mRNA=x[2], rho_protein=y[1], p_protein=y[2], n=y[3])) }
cat("\n-- VARIANT: miRNA label collapsed as the MEAN of all mature arms --\n")
print(v2, digits=4, row.names=FALSE)
write.csv(v2, "results/v2/v2_cptac_mir29_armaveraged.csv", row.names=FALSE)
write.csv(out, "results/v2/v2_cptac_mir29_axis.csv", row.names=FALSE)

## ---------------- global: all miRNA_target edges in the network ------------
E <- read.delim("data/canonical_edges.tsv", stringsAsFactors=FALSE)
tier <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors=FALSE)
E <- E[E$edge_type=="miRNA_target",]
E$tier <- tier$tier[match(paste(E$source,E$target), paste(tier$source,tier$target))]
E$arm <- vapply(E$source, armmap, character(1))
E <- E[!is.na(E$arm) & E$target %in% rownames(P) & E$target %in% rownames(R), ]
cat("\nnetwork miRNA_target edges testable in CPTAC:", nrow(E),
    " (", length(unique(E$arm)), "miRNAs,", length(unique(E$target)), "targets )\n")
E$rho_mRNA <- mapply(function(a,b) sp(MI[a,cs_rm], R[b,cs_rm]), E$arm, E$target)
E$rho_prot <- mapply(function(a,b) sp(MI[a,cs_pm], P[b,cs_pm]), E$arm, E$target)
g <- function(v) c(n=sum(!is.na(v)), frac_neg=mean(v<0, na.rm=TRUE), mean=mean(v, na.rm=TRUE))
cat("\n-- global miRNA_target concordance (expected sign = negative) --\n")
sm <- rbind(mRNA=g(E$rho_mRNA), protein=g(E$rho_prot))
print(sm, digits=4)
cat("paired Wilcoxon rho_protein vs rho_mRNA p =",
    signif(wilcox.test(E$rho_prot, E$rho_mRNA, paired=TRUE)$p.value,3),
    "; |rho| p =", signif(wilcox.test(abs(E$rho_prot), abs(E$rho_mRNA), paired=TRUE)$p.value,3), "\n")
cat("prop.test frac_neg protein vs mRNA p =",
    signif(prop.test(c(sum(E$rho_prot<0,na.rm=TRUE), sum(E$rho_mRNA<0,na.rm=TRUE)),
                     c(sum(!is.na(E$rho_prot)), sum(!is.na(E$rho_mRNA))))$p.value,3), "\n")
cat("fraction of edges where |rho_protein| > |rho_mRNA| :",
    round(mean(abs(E$rho_prot) > abs(E$rho_mRNA), na.rm=TRUE),4), "\n")
by_t <- do.call(rbind, lapply(split(E, E$tier), function(d)
  data.frame(tier=d$tier[1], n=nrow(d), fneg_mRNA=mean(d$rho_mRNA<0,na.rm=TRUE),
             fneg_prot=mean(d$rho_prot<0,na.rm=TRUE),
             mean_mRNA=mean(d$rho_mRNA,na.rm=TRUE), mean_prot=mean(d$rho_prot,na.rm=TRUE))))
print(by_t, digits=4, row.names=FALSE)
write.csv(E[,c("source","arm","target","tier","rho_mRNA","rho_prot")],
          "results/v2/v2_cptac_global_edges.csv", row.names=FALSE)
cat("DONE 06\n")
