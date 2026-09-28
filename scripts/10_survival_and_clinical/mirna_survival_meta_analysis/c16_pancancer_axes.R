## ==========================================================================
## c16_pancancer_axes.R -- Part D
## miR-29a -| COL1A1/COL3A1 and miR-101 -| EZH2 in every TCGA cohort with
## enough samples; per-cohort CAF score; does axis strength track CAF content?
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(ppcor)})
setwd("/path/to/revision")
OUT<-"results/v6"; C<-"cache/v6/pancancer"; ML<-"/path/to/home/Desktop/DD/R_GPR/ML"
MINN <- 50

G <- fread(file.path(C,"pancan_gene_subset.tsv"))
gs <- G[[1]]; G <- as.matrix(G[,-1]); rownames(G) <- gs
cat("gene matrix:", dim(G), "\n")
MI <- fread(file.path(ML,"pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz"))
mr <- MI[[1]]; MI <- as.matrix(MI[,-1]); rownames(MI) <- mr
cat("miRNA matrix:", dim(MI), "\n")

SV <- fread(file.path(ML,"Survival_SupplementalTable_S1_20171025_xena_sp.txt"))
PH <- fread(file.path(ML,"TCGA_phenotype_denseDataOnlyDownload.tsv.gz"))
setnames(SV, "cancer type abbreviation", "cohort")

common <- intersect(colnames(G), colnames(MI))
cat("samples with both assays:", length(common), "\n")
ph <- PH[match(common, PH$sample)]
sv <- SV[match(common, SV$sample)]
keep <- !is.na(sv$cohort) & ph$sample_type %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Metastatic")
## restrict to primary tumours; SKCM is predominantly Metastatic, LAML is blood
keep <- keep & (ph$sample_type=="Primary Tumor" |
                (sv$cohort=="SKCM" & ph$sample_type=="Metastatic") |
                (sv$cohort=="LAML" & grepl("Primary Blood", ph$sample_type)))
common <- common[keep]; ph<-ph[keep]; sv<-sv[keep]
cat("primary-tumour samples used:", length(common), "cohorts:", length(unique(sv$cohort)), "\n")
G <- G[, common, drop=FALSE]; MI <- MI[, common, drop=FALSE]

CAFG <- readLines("results/v2/cafA_gene_list_v3.txt"); CAFG <- CAFG[CAFG %in% rownames(G)]
cat("CAF-A genes measurable:", length(CAFG), "\n")
## pan-cancer comparable CAF score: z each gene across ALL primary tumours, mean
zc <- function(x){ s<-sd(x,na.rm=TRUE); if(!is.finite(s)||s==0) return(rep(NA,length(x))); (x-mean(x,na.rm=TRUE))/s }
Z <- t(apply(G[CAFG,,drop=FALSE], 1, zc))
cafscore <- colMeans(Z, na.rm=TRUE)

MIR <- c(`miR-29a`="hsa-miR-29a-3p", `miR-101`="hsa-miR-101-3p",
         `miR-29b`="hsa-miR-29b-3p", `miR-29c`="hsa-miR-29c-3p",
         `miR-21`="hsa-miR-21-5p", `miR-141`="hsa-miR-141-3p")
stopifnot(all(MIR %in% rownames(MI)))
AXES <- list(c("miR-29a","COL1A1"), c("miR-29a","COL3A1"), c("miR-101","EZH2"),
             c("miR-29b","COL1A1"), c("miR-29c","COL1A1"),
             c("miR-21","COL1A1"), c("miR-141","EZH2"))

rows <- list()
for(ck in sort(unique(sv$cohort))){
  i <- which(sv$cohort==ck); if(length(i)<MINN) next
  cs <- cafscore[i]
  for(a in AXES){
    m <- as.numeric(MI[MIR[[a[1]]], i]); g <- as.numeric(G[a[2], i])
    ok <- is.finite(m)&is.finite(g)&is.finite(cs)
    if(sum(ok)<MINN) next
    r <- suppressWarnings(cor.test(m[ok], g[ok], method="spearman", exact=FALSE))
    pc <- try(pcor.test(rank(m[ok]), rank(g[ok]), rank(cs[ok]), method="pearson"), silent=TRUE)
    rows[[length(rows)+1]] <- data.frame(cohort=ck, n=sum(ok),
      miRNA=a[1], gene=a[2], axis=paste0(a[1]," -| ",a[2]),
      rho=unname(r$estimate), p=r$p.value,
      rho_partialCAF=if(inherits(pc,"try-error")) NA else pc$estimate,
      p_partialCAF=if(inherits(pc,"try-error")) NA else pc$p.value,
      caf_mean=mean(cs[ok]), caf_sd=sd(cs[ok]),
      caf_vs_gene_rho=suppressWarnings(cor(cs[ok], g[ok], method="spearman")),
      stringsAsFactors=FALSE)
  }
}
PC <- rbindlist(rows)
PC[, q := p.adjust(p,"BH"), by=axis]
fwrite(PC, file.path(OUT,"pancancer_axes.csv"))

## cohort-level CAF summary
CS <- data.table(sample=common, cohort=sv$cohort, caf=cafscore)
CSU <- CS[, .(n=.N, caf_mean=mean(caf,na.rm=TRUE), caf_median=median(caf,na.rm=TRUE),
              caf_sd=sd(caf,na.rm=TRUE)), by=cohort][order(-caf_mean)]
fwrite(CSU, file.path(OUT,"pancancer_caf_score_by_cohort.csv"))

## does axis strength track CAF content across cohorts?
tr <- list()
for(ax in unique(PC$axis)){
  s <- PC[axis==ax]
  if(nrow(s)<8) next
  a <- cor.test(s$caf_mean, s$rho, method="spearman", exact=FALSE)
  b <- cor.test(s$caf_sd,  s$rho, method="spearman", exact=FALSE)
  tr[[length(tr)+1]] <- data.frame(axis=ax, k_cohorts=nrow(s),
    median_rho=median(s$rho), n_negative=sum(s$rho<0), n_sig_negative=sum(s$rho<0 & s$q<0.05),
    rho_vs_cafmean=unname(a$estimate), p_vs_cafmean=a$p.value,
    rho_vs_cafsd=unname(b$estimate), p_vs_cafsd=b$p.value,
    median_rho_partialCAF=median(s$rho_partialCAF, na.rm=TRUE), stringsAsFactors=FALSE)
}
TR <- rbindlist(tr); fwrite(TR, file.path(OUT,"pancancer_axes_vs_caf.csv"))
cat("\n=== axis summary across TCGA cohorts (n>=50) ===\n"); print(as.data.frame(TR))
cat("\n=== BRCA vs the rest, per axis ===\n")
print(as.data.frame(PC[axis %in% c("miR-29a -| COL1A1","miR-29a -| COL3A1","miR-101 -| EZH2")][
  , .(cohort,axis,n,rho=round(rho,3),q=signif(q,3),rho_pCAF=round(rho_partialCAF,3),caf=round(caf_mean,3))][
  order(axis,rho)]))
