## ==========================================================================
## N1_build_cohorts.R
## Standardise every transcriptomic cohort into
##   list(name, accession, platform, X = symbols x samples (log2), pheno, notes)
## and save cache/newcohorts/cohorts.rds
## Every number printed here is computed in this run.
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
setwd("/path/to/revision")
source("scripts/07_expression_validation/external_cohorts/N0_geo_utils.R")
CA <- "cache/newcohorts"
cohorts <- list()
add <- function(o){ cohorts[[o$name]] <<- o
  cat(sprintf("[COHORT] %-12s acc=%-28s platform=%-38s genes=%6d samples=%5d\n",
              o$name, o$accession, o$platform, nrow(o$X), ncol(o$X))) }

## ------------------------------------------------------------ TCGA-BRCA ---
g <- readRDS("data/brca_gene_expr.rds")
ph <- readRDS("data/brca_pheno.rds")
sv <- as.data.frame(readRDS("data/brca_survival_clean.rds"))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]
nor <- ph$sample[ph$sample_type=="Solid Tissue Normal"]
cat("TCGA primary tumours:", length(tum), " normals:", length(nor), "\n")
Xt <- g[, intersect(colnames(g), tum), drop=FALSE]
pt <- data.frame(sample=colnames(Xt), stringsAsFactors=FALSE)
m <- match(pt$sample, sv$sample_id)
pt$OS_time <- sv$OS.time[m]/30.4375; pt$OS_event <- sv$OS[m]
pt$age <- sv$age_at_initial_pathologic_diagnosis[m]
st <- sv$ajcc_pathologic_tumor_stage[m]
pt$stage <- ifelse(grepl("Stage I[^IV]*$|Stage I$",st),"I",
            ifelse(grepl("Stage II",st),"II", ifelse(grepl("Stage III",st),"III",
            ifelse(grepl("Stage IV",st),"IV", NA))))
add(list(name="TCGA_BRCA", accession="TCGA-BRCA (UCSC Xena)",
         platform="Illumina HiSeq RNA-seq (log2 RSEM+1)", X=Xt, pheno=pt,
         normals=g[, intersect(colnames(g), nor), drop=FALSE],
         notes="reference cohort; tumour-vs-normal DE available"))

## ------------------------------------------------------------- METABRIC ---
mb <- readRDS("data/metabric.rds")
Xm <- mb$M; cl <- as.data.frame(mb$cl)
stopifnot(all(colnames(Xm) %in% cl$SAMPLE_ID) || all(colnames(Xm) %in% cl$PATIENT_ID))
key <- if(all(colnames(Xm) %in% cl$PATIENT_ID)) "PATIENT_ID" else "SAMPLE_ID"
mm <- match(colnames(Xm), cl[[key]])
pm <- data.frame(sample=colnames(Xm),
                 OS_time=cl$OS_time[mm], OS_event=cl$OS_event[mm],
                 RFS_time=cl$RFS_time[mm], RFS_event=cl$RFS_event[mm],
                 age=cl$age[mm], grade=cl$grade[mm], stage=cl$stage_group[mm],
                 stringsAsFactors=FALSE)
add(list(name="METABRIC", accession="cBioPortal brca_metabric",
         platform="Illumina HT-12 v3 microarray", X=Xm, pheno=pm,
         notes="OS and RFS available"))

## ----------------------------------------------------------------- GTEx ---
f <- file.path(CA,"gtex_breast_tpm.gct.gz")
if(file.exists(f) && file.size(f) > 4e7){
  dt <- fread(cmd=paste("zcat", shQuote(f), "| tail -n +3"), sep="\t", header=TRUE)
  sym <- dt$Description
  V <- as.matrix(dt[, -(1:3), with=FALSE])
  rownames(V) <- sym
  Xg <- log2(V + 1)
  mu <- rowMeans(Xg); ord <- order(rownames(Xg), -mu)
  Xg <- Xg[ord, , drop=FALSE]; Xg <- Xg[!duplicated(rownames(Xg)), , drop=FALSE]
  keep <- rowMeans(V[match(rownames(Xg), rownames(V)), , drop=FALSE] > 0.5) > 0.2
  cat("GTEx genes before/after expression filter:", nrow(Xg), sum(keep), "\n")
  Xg <- Xg[keep, , drop=FALSE]
  sa <- fread(file.path(CA,"GTEx_pheno.txt"), sep="\t", quote="")
  sa <- sa[match(colnames(Xg), sa$SAMPID), ]
  subj <- sub("^(GTEX-[^-]+).*$","\\1", colnames(Xg))
  sj <- fread(file.path(CA,"GTEx_subject.txt"), sep="\t")
  pg <- data.frame(sample=colnames(Xg), subject=subj,
                   tissue=sa$SMTSD, rin=sa$SMRIN, ischemia=sa$SMTSISCH,
                   sex=sj$SEX[match(subj, sj$SUBJID)],
                   age_bin=sj$AGE[match(subj, sj$SUBJID)], stringsAsFactors=FALSE)
  add(list(name="GTEx_breast", accession="GTEx v8 breast mammary tissue",
           platform="Illumina TrueSeq RNA-seq (log2 TPM+1)", X=Xg, pheno=pg,
           notes="NON-TUMOUR baseline; no survival"))
} else cat("!! GTEx file missing/short\n")

## ------------------------------------------------------- GEO array sets ---
gpl570 <- NULL
gplf <- file.path(CA,"GPL570.annot.gz")
if(file.exists(gplf) && file.size(gplf) > 8e6) gpl570 <- read_gpl_annot(gplf)
if(!is.null(gpl570)) cat("GPL570 annotation rows:", nrow(gpl570), "\n")

## (GEO array series are parsed in N1b/N1c, not here)

saveRDS(cohorts, file.path(CA,"cohorts_base.rds"), compress=FALSE)
cat("\nSaved base cohorts:", paste(names(cohorts), collapse=", "), "\n")
