## ===========================================================================
## PART C. CPTAC pan-cancer: miR-29 -> collagen at PROTEIN vs mRNA level.
## Independent implementation, written fresh 2026-09-09.
## ===========================================================================
setwd("/path/to/revision")
OUT <- "results/v2"
CP  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
MIR <- file.path(CP, "miRNA_BCM_v1")
PRO <- file.path(CP, "Proteome_BCM_GENCODE_v34_harmonized_v1/Proteome_BCM_GENCODE_v34_harmonized_v1")
RNA <- file.path(CP, "RNA_BCM_v1/RNA_BCM_v1")
options(width=200)

## cohort key: miRNA-file name  ->  proteome/RNA-file names
COH <- data.frame(
  cohort   = c("BRCA","ccRCC","COAD","GBM","HGSC","HNSCC","LSCC","LUAD","PDAC","UCEC"),
  mir_tag  = c("BRCA","ccRCC","COAD","GBM","HGSC","HNSCC","LSCC","LUAD","PDAC","UCEC"),
  prot_tag = c("BRCA","CCRCC","COAD","GBM","OV",  "HNSCC","LSCC","LUAD","PDAC","UCEC"),
  rna_tag  = c("BRCA","ccRCC","COAD","GBM","HGSC","HNSCC","LSCC","LUAD","PDAC","UCEC"),
  stringsAsFactors = FALSE)

ann <- read.delim(file.path(PRO,"README","Gene_annotation_and_representable_isoform_mapping_table.txt"),
                  stringsAsFactors=FALSE, check.names=FALSE)
ens2sym <- setNames(ann$gene_name, sub("\\..*$","",ann$gene))
cat("annotation rows:", nrow(ann), " unique ENSG:", length(unique(names(ens2sym))), "\n")

readmat <- function(path){
  d <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
  as.matrix(d)
}
tosym <- function(M){
  s <- ens2sym[sub("\\..*$","",rownames(M))]
  k <- !is.na(s) & s != ""
  M <- M[k,,drop=FALSE]; s <- s[k]
  ## collapse duplicate symbols by the row with the most non-missing values
  o <- order(s, -rowSums(is.finite(M)))
  M <- M[o,,drop=FALSE]; s <- s[o]
  keep <- !duplicated(s)
  M <- M[keep,,drop=FALSE]; rownames(M) <- s[keep]
  M
}
sp <- function(x,y){
  k <- is.finite(x)&is.finite(y); n <- sum(k)
  if (n < 12) return(c(rho=NA,p=NA,n=n))
  ct <- suppressWarnings(cor.test(x[k],y[k],method="spearman",exact=FALSE))
  c(rho=unname(ct$estimate), p=unname(ct$p.value), n=n)
}
## Steiger's test for two dependent correlations sharing variable x
steiger <- function(r12, r13, r23, n){
  if (any(!is.finite(c(r12,r13,r23))) || n < 15) return(NA)
  R <- (1-r12^2-r13^2-r23^2) + 2*r12*r13*r23
  t <- (r12-r13)*sqrt(((n-1)*(1+r23)) / (2*R*(n-1)/(n-3) + ((r12+r13)^2/4)*(1-r23)^3))
  2*pt(-abs(t), df=n-3)
}

MIRSETS <- list(
  "hsa-miR-29a"    = c("hsa-miR-29a-5p","hsa-miR-29a-3p"),
  "hsa-miR-29b"    = c("hsa-miR-29b-1-5p","hsa-miR-29b-2-5p","hsa-miR-29b-3p"),
  "hsa-miR-29c"    = c("hsa-miR-29c-5p","hsa-miR-29c-3p"),
  "hsa-miR-29a-3p" = "hsa-miR-29a-3p",
  "hsa-miR-29b-3p" = "hsa-miR-29b-3p",
  "hsa-miR-29c-3p" = "hsa-miR-29c-3p")
TARGETS <- c("COL1A1","COL3A1","COL1A2","COL5A1","COL5A2","COL6A3","COL11A1","FN1","SPARC","LOX")
CAFB    <- c("DCN","LUM","FAP","THY1")
cafA    <- readLines(file.path(OUT,"cafA_gene_list_v3.txt"))

inv <- list(); res <- list()
for (i in seq_len(nrow(COH))){
  ch <- COH$cohort[i]
  fm <- file.path(MIR, sprintf("%s_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt", COH$mir_tag[i]))
  fp <- file.path(PRO, sprintf("%s_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt", COH$prot_tag[i]))
  fr <- file.path(RNA, sprintf("%s_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt", COH$rna_tag[i]))
  hasm <- file.exists(fm); hasp <- file.exists(fp); hasr <- file.exists(fr)
  if (!(hasm && hasp)) { cat("SKIP", ch, "mir:",hasm," prot:",hasp,"\n"); next }
  Mi <- readmat(fm); Pr <- tosym(readmat(fp))
  Rn <- if (hasr) tosym(readmat(fr)) else NULL
  smp <- intersect(colnames(Mi), colnames(Pr))
  smr <- if (!is.null(Rn)) Reduce(intersect, list(colnames(Mi), colnames(Pr), colnames(Rn))) else character(0)
  cat(sprintf("%-6s miRNA %d x %d | prot %d x %d | rna %s | overlap mir&prot = %d | all three = %d\n",
      ch, nrow(Mi), ncol(Mi), nrow(Pr), ncol(Pr),
      if(is.null(Rn)) "-" else paste(dim(Rn),collapse="x"), length(smp), length(smr)))
  inv[[length(inv)+1]] <- data.frame(cohort=ch, n_mirna_samples=ncol(Mi), n_protein_samples=ncol(Pr),
      n_rna_samples=if(is.null(Rn)) 0L else ncol(Rn), n_overlap_mirna_protein=length(smp),
      n_overlap_all3=length(smr), n_mirnas=nrow(Mi), n_proteins=nrow(Pr), stringsAsFactors=FALSE)
  if (length(smp) < 25) { cat("  too few overlapping samples, skipping stats\n"); next }

  Mi2 <- Mi[, smp, drop=FALSE]; Pr2 <- Pr[, smp, drop=FALSE]
  Rn2 <- if (!is.null(Rn) && length(smr) >= 25) Rn[, smr, drop=FALSE] else NULL

  ## CAF score from the PROTEOME (protein-level stromal content), and from RNA
  zc <- function(X, gs){
    gs <- intersect(gs, rownames(X)); if (!length(gs)) return(NULL)
    S <- X[gs,,drop=FALSE]
    mu <- apply(S,1,mean,na.rm=TRUE); sd <- apply(S,1,sd,na.rm=TRUE)
    ok <- is.finite(sd) & sd>0
    if (!any(ok)) return(NULL)
    list(score=colMeans(((S-mu)/sd)[ok,,drop=FALSE], na.rm=TRUE), n=sum(ok))
  }
  cafP <- zc(Pr2, cafA); cafPB <- zc(Pr2, CAFB)
  cafR <- if (!is.null(Rn2)) zc(Rn2, cafA) else NULL

  for (ml in names(MIRSETS)){
    arms <- intersect(MIRSETS[[ml]], rownames(Mi2))
    if (!length(arms)) next
    xm <- colMeans(Mi2[arms,,drop=FALSE], na.rm=TRUE)
    xr <- if (!is.null(Rn2)) colMeans(Mi[arms, smr, drop=FALSE], na.rm=TRUE) else NULL
    for (tg in TARGETS){
      yp <- if (tg %in% rownames(Pr2)) Pr2[tg,] else NULL
      yr <- if (!is.null(Rn2) && tg %in% rownames(Rn2)) Rn2[tg,] else NULL
      if (is.null(yp) && is.null(yr)) next
      ap <- if (!is.null(yp)) sp(xm, yp) else c(rho=NA,p=NA,n=0)
      ar <- if (!is.null(yr)) sp(xr, yr) else c(rho=NA,p=NA,n=0)
      ## partial (protein-level CAF adjustment)
      pp <- c(rho=NA,p=NA)
      if (!is.null(yp) && !is.null(cafP)){
        k <- is.finite(xm)&is.finite(yp)&is.finite(cafP$score); nn <- sum(k)
        if (nn>=15){
          rx <- rank(xm[k]); ry <- rank(yp[k]); rz <- rank(cafP$score[k])
          rxy<-cor(rx,ry); rxz<-cor(rx,rz); ryz<-cor(ry,rz)
          rr <- (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
          tt <- rr*sqrt((nn-3)/(1-rr^2))
          pp <- c(rho=rr, p=2*pt(-abs(tt), df=nn-3))
        }
      }
      ## Steiger test protein vs mRNA (same miRNA x, dependent samples smr)
      st <- NA
      if (!is.null(yp) && !is.null(yr) && length(smr)>=25){
        yp2 <- Pr2[tg, smr]; yr2 <- Rn2[tg, smr]; xm2 <- colMeans(Mi[arms, smr, drop=FALSE], na.rm=TRUE)
        k <- is.finite(xm2)&is.finite(yp2)&is.finite(yr2)
        if (sum(k) >= 20)
          st <- steiger(cor(xm2[k],yp2[k],method="spearman"),
                        cor(xm2[k],yr2[k],method="spearman"),
                        cor(yp2[k],yr2[k],method="spearman"), sum(k))
      }
      res[[length(res)+1]] <- data.frame(cohort=ch, miRNA=ml, arms=paste(arms,collapse="+"),
        target=tg, n_protein=ap["n"], rho_protein=ap["rho"], p_protein=ap["p"],
        n_mrna=ar["n"], rho_mrna=ar["rho"], p_mrna=ar["p"],
        rho_protein_partialCAF=pp["rho"], p_protein_partialCAF=pp["p"],
        steiger_p_protein_vs_mrna=st,
        cafA_proteins_used=if(is.null(cafP)) 0L else cafP$n, stringsAsFactors=FALSE)
    }
  }
  ## record the CAF/collagen relationship at protein level
  if (!is.null(cafP)){
    for (tg in c("COL1A1","COL3A1")){
      if (tg %in% rownames(Pr2)){
        a <- sp(cafP$score, Pr2[tg,])
        res[[length(res)+1]] <- data.frame(cohort=ch, miRNA="CAF_A_protein", arms="score",
          target=tg, n_protein=a["n"], rho_protein=a["rho"], p_protein=a["p"],
          n_mrna=0, rho_mrna=NA, p_mrna=NA, rho_protein_partialCAF=NA, p_protein_partialCAF=NA,
          steiger_p_protein_vs_mrna=NA, cafA_proteins_used=cafP$n, stringsAsFactors=FALSE)
      }
    }
  }
}
INV <- do.call(rbind, inv); rownames(INV) <- NULL
RES <- do.call(rbind, res); rownames(RES) <- NULL
RES$q_BH_protein <- p.adjust(RES$p_protein, method="BH")
write.csv(INV, file.path(OUT,"cptac_pancancer_inventory_v3.csv"), row.names=FALSE)
write.csv(RES, file.path(OUT,"cptac_pancancer_mir29.csv"), row.names=FALSE)
cat("\nwrote cptac_pancancer_mir29.csv rows =", nrow(RES), "\n\n")

cat("=== INVENTORY ===\n"); print(INV)
cat("\n=== miR-29 (arm-averaged) -> COL1A1 / COL3A1 at PROTEIN level, per cohort ===\n")
k <- RES$miRNA %in% c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c") & RES$target %in% c("COL1A1","COL3A1")
print(RES[k, c("cohort","miRNA","target","n_protein","rho_protein","p_protein",
               "n_mrna","rho_mrna","p_mrna","rho_protein_partialCAF","steiger_p_protein_vs_mrna")], digits=3)
cat("\n=== dominant-arm (-3p) version ===\n")
k3 <- RES$miRNA %in% c("hsa-miR-29a-3p","hsa-miR-29b-3p","hsa-miR-29c-3p") & RES$target %in% c("COL1A1","COL3A1")
print(RES[k3, c("cohort","miRNA","target","n_protein","rho_protein","p_protein","rho_mrna",
                "rho_protein_partialCAF")], digits=3)

## per-cohort summary: protein vs mRNA |rho| across all miR-29 x ECM pairs
cat("\n=== per-cohort protein-vs-mRNA comparison (miR-29 arm-averaged x 10 ECM targets) ===\n")
sm <- do.call(rbind, lapply(split(RES[RES$miRNA %in% c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"),],
                                  RES$cohort[RES$miRNA %in% c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")]), function(d){
  d <- d[is.finite(d$rho_protein) & is.finite(d$rho_mrna), ]
  if (!nrow(d)) return(NULL)
  bt <- binom.test(sum(abs(d$rho_protein) > abs(d$rho_mrna)), nrow(d))
  neg <- binom.test(sum(d$rho_protein < 0), nrow(d))
  data.frame(cohort=d$cohort[1], n_pairs=nrow(d),
    mean_rho_protein=mean(d$rho_protein), mean_rho_mrna=mean(d$rho_mrna),
    n_protein_stronger=sum(abs(d$rho_protein)>abs(d$rho_mrna)),
    binom_p_stronger=bt$p.value,
    n_negative_protein=sum(d$rho_protein<0), binom_p_negative=neg$p.value,
    wilcox_p_rho=suppressWarnings(wilcox.test(d$rho_protein, d$rho_mrna, paired=TRUE)$p.value),
    stringsAsFactors=FALSE)
}))
rownames(sm) <- NULL
print(sm, digits=3)
write.csv(sm, file.path(OUT,"cptac_pancancer_protein_vs_mrna_v3.csv"), row.names=FALSE)
cat("\nDONE\n")
