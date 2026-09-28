## ==========================================================================
## M3_build_mirna_cohorts.R
## Assemble every breast cohort on disk that carries a miRNA assay, with its
## phenotype/outcome table. Saves cache/v4/multicohort/mirna_cohorts.rds and
## results/v4/multicohort_mirna_inventory.csv
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
source("scripts/v2/N0_geo_utils.R")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
num <- function(v){ v[v %in% c(".","NA","na","not available","Unknown","unknown","")] <- NA
                    suppressWarnings(as.numeric(v)) }

## ---- canonical mature-miRNA name resolution ------------------------------
## Old arrays annotate the dominant (3p) arm as "hsa-miR-130a"; the passenger
## arm is "hsa-miR-130a*". Never accept a "*", "-5p" or "-pre" row as the 3p arm.
resolve <- function(rn, base){
  cand <- c(base, paste0(base,"-3p"), tolower(base))
  hit <- rn[rn %in% cand]
  if(!length(hit)){
    ## permissive: exact base ignoring case, still excluding star/5p/pre
    hit <- rn[tolower(rn) %in% tolower(cand)]
  }
  hit <- hit[!grepl("\\*|-5p$|-pre$|pre$", hit)]
  unique(hit)
}


## ---- platform annotation reader that prefers a real mature-miRNA column ---
read_mirna_annot <- function(path){
  con <- if(grepl("\\.gz$",path)) gzfile(path,"rt") else file(path,"rt")
  L <- readLines(con); close(con)
  b <- grep("^!platform_table_begin", L)[1]; e <- grep("^!platform_table_end", L)[1]
  dt <- data.table::fread(text=paste(L[(b+1):(e-1)], collapse="\n"), sep="\t",
                          quote="", header=TRUE, colClasses="character")
  nm <- colnames(dt)
  ord <- c(grep("^miRNA_ID$", nm), grep("^miRNA_ID_LIST$", nm),
           grep("^ILMN_Gene$", nm), grep("^Gene Symbol$", nm), grep("^SYMBOL$", nm))
  ord <- unique(ord); stopifnot(length(ord)>0)
  sym <- rep(NA_character_, nrow(dt))
  for(i in ord){ v <- trimws(dt[[i]]); take <- is.na(sym) & v!="" & grepl("^hsa-", v); sym[take] <- v[take] }
  cat(sprintf("  annot %s: %d rows, columns tried [%s], hsa- symbols resolved %d\n",
      basename(path), nrow(dt), paste(nm[ord], collapse=","), sum(!is.na(sym))))
  data.frame(ID=as.character(dt[[1]]), Symbol=sym, stringsAsFactors=FALSE)
}

## ---- collapse a probe matrix to unique mature-miRNA rows ------------------
collapse_mirna <- function(X, ann, label){
  sy <- ann$Symbol[match(rownames(X), ann$ID)]
  ok <- !is.na(sy) & grepl("^hsa-", sy) & !grepl(",", sy)
  cat(sprintf("  %s: probes %d -> annotated single hsa- miRNA %d\n", label, nrow(X), sum(ok)))
  X <- X[ok,,drop=FALSE]; sy <- sy[ok]
  mu <- rowMeans(X, na.rm=TRUE); o <- order(sy, -mu)
  X <- X[o,,drop=FALSE]; sy <- sy[o]
  keep <- !duplicated(sy); X <- X[keep,,drop=FALSE]; rownames(X) <- sy[keep]
  X
}

CO <- list()
addc <- function(o){ CO[[o$name]] <<- o
  cat(sprintf("[miRNA COHORT] %-14s acc=%-12s plat=%-46s miRNAs=%5d n=%4d  endpoints=%s\n",
    o$name,o$accession,substr(o$platform,1,46),nrow(o$M),ncol(o$M),
    paste(o$endpoints, collapse="/"))) }

## ============================================================== TCGA-BRCA ==
mi <- readRDS("data/brca_mirna_expr.rds")
ph <- readRDS("data/brca_pheno.rds")
sv <- as.data.frame(readRDS("data/brca_survival_clean.rds"))
clin <- fread("data/brca_clinicalMatrix.tsv", sep="\t", quote="")
cat("TCGA miRNA matrix:", dim(mi), "\n")
tum <- ph$sample[ph$sample_type=="Primary Tumor"]; nor <- ph$sample[ph$sample_type=="Solid Tissue Normal"]
Mt <- mi[, intersect(colnames(mi), tum), drop=FALSE]
Mn <- mi[, intersect(colnames(mi), nor), drop=FALSE]
cat("TCGA miRNA tumours:", ncol(Mt), " normals:", ncol(Mn), "\n")
m <- match(colnames(Mt), sv$sample_id)
sidcol <- names(clin)[1]
mc <- match(colnames(Mt), clin[[sidcol]])
pt <- data.frame(sample=colnames(Mt),
  OS_time=sv$OS.time[m]/30.4375, OS_event=sv$OS[m],
  PFI_time=sv$PFI.time[m]/30.4375, PFI_event=sv$PFI[m],
  DSS_time=sv$DSS.time[m]/30.4375, DSS_event=sv$DSS[m],
  age=sv$age_at_initial_pathologic_diagnosis[m],
  stage=sv$ajcc_pathologic_tumor_stage[m],
  pam50=clin$PAM50Call_RNAseq[mc], ER=clin$ER_Status_nature2012[mc],
  stringsAsFactors=FALSE)
cat("TCGA PAM50 available:", sum(!is.na(pt$pam50) & pt$pam50!=""), "\n")
print(table(pt$pam50, useNA="ifany")); print(table(pt$ER, useNA="ifany"))
addc(list(name="TCGA_BRCA", accession="TCGA-BRCA", platform="Illumina HiSeq miRNA-seq (log2 RPM)",
  M=Mt, normals=Mn, pheno=pt, endpoints=c("OS","PFI","DSS"),
  notes="reference miRNA cohort; matched mRNA available"))

## =============================================================== GSE19783 ==
g <- readRDS("cache/external2/gse19783.rds")
M <- g$GSE19783_miRNA
p <- parse_series_matrix("cache/external2/GSE19783-GPL8227_series_matrix.txt.gz")
mm <- p$meta
d <- data.frame(title=trimws(sub("\\s*\\(miRNA\\)\\s*$","", mm$Sample_title)), gsm=mm$geo_accession,
  subtype=pull_char(mm,"breast cancer subtype"),
  death  =pull_char(mm,"death status"),
  dfs    =num(pull_char(mm,"disease free survival time \\(months\\)")),
  ER     =pull_char(mm,"estrogen receptor status"),
  HER2   =pull_char(mm,"her2 \\(fish\\) status"),
  TP53   =pull_char(mm,"tp53 mutation status"), stringsAsFactors=FALSE)
cat("GSE19783 series-matrix samples:", nrow(d), " matched to miRNA columns:",
    sum(colnames(M) %in% d$title), "\n")
stopifnot(sum(colnames(M) %in% d$title) == ncol(M))
d <- d[match(colnames(M), d$title), ]
d$sample <- colnames(M)
## endpoint: any-cause death within the DFS follow-up (the deposited fields)
d$OS_event <- ifelse(is.na(d$death), NA, as.integer(grepl("^Dead", d$death)))
d$BCSS_event <- ifelse(is.na(d$death), NA, as.integer(d$death=="Dead of BC"))
d$OS_time <- d$dfs; d$BCSS_time <- d$dfs
print(table(d$death, useNA="ifany"))
cat("GSE19783 OS events:", sum(d$OS_event, na.rm=TRUE),
    " BCSS events:", sum(d$BCSS_event, na.rm=TRUE),
    " with time:", sum(is.finite(d$OS_time)), "\n")
addc(list(name="GSE19783", accession="GSE19783", platform="Agilent-019118 human miRNA 2.0 (GPL8227)",
  M=M, normals=NULL, pheno=d, endpoints=c("OS","BCSS"),
  notes="101 tumours, matched mRNA on GPL6480; time variable is the deposited disease-free-survival time"))

## =============================================================== GSE22216 ==
p <- parse_series_matrix(file.path(CA,"GSE22216_series_matrix.txt.gz"))
ann <- read_mirna_annot(file.path(CA,"GPL8178_family.soft.gz"))
X <- collapse_mirna(p$X, ann, "GSE22216")
mm <- p$meta
d <- data.frame(sample=mm$geo_accession, gsm=mm$geo_accession,
  age  =num(pull_char(mm,"patient age")),
  size =num(pull_char(mm,"tumour size")),
  node =num(pull_char(mm,"nodes involved")),
  ER   =pull_char(mm,"er status"),
  grade=num(pull_char(mm,"tumour grade")),
  DRFS_event=num(pull_char(mm,"distant-relapse event")),
  DRFS_time =num(pull_char(mm,"distant-relapse free survival")),
  stringsAsFactors=FALSE)
d$DRFS_time <- d$DRFS_time*12   ## deposited in years
d <- d[match(colnames(X), d$sample), ]
stopifnot(identical(d$sample, colnames(X)))
cat("GSE22216 n:", nrow(d), " DRFS events:", sum(d$DRFS_event, na.rm=TRUE),
    " median follow-up (months):", round(median(d$DRFS_time, na.rm=TRUE),1), "\n")
addc(list(name="GSE22216", accession="GSE22216 (= GSE22220 miRNA arm)",
  platform="Illumina Human v1 MicroRNA array (GPL8178)",
  M=X, normals=NULL, pheno=d, endpoints=c("DRFS"),
  notes="210 early primary breast cancers (Buffa et al.); matched mRNA is GSE22219"))

## =============================================================== GSE37405 ==
PL <- c(GPL13703="GSE37405-GPL13703_series_matrix.txt.gz",
        GPL14149="GSE37405-GPL14149_series_matrix.txt.gz",
        GPL15462="GSE37405-GPL15462_series_matrix.txt.gz")
Xs <- list(); ds <- list()
for(pl in names(PL)){
  pp <- parse_series_matrix(file.path(CA, PL[pl]))
  an <- read_mirna_annot(file.path(CA, paste0(pl,"_self.txt")))
  x <- collapse_mirna(pp$X, an, paste0("GSE37405 ",pl))
  mm <- pp$meta
  dd <- data.frame(sample=mm$geo_accession, platform=pl,
    slide=pull_char(mm,"slide id"),
    age  =num(pull_char(mm,"age at op")),
    size =num(pull_char(mm,"size \\(mm\\)")),
    node =num(pull_char(mm,"nodal-status, positive")),
    t_ok =num(pull_char(mm,"time.to.ok")),
    t_rec=num(pull_char(mm,"time.to.recurrence \\(years\\)")),
    t_dth=num(pull_char(mm,"time.to.death \\(years\\)")),
    stringsAsFactors=FALSE)
  cat(sprintf("GSE37405 %s: probes %d -> miRNAs %d ; samples %d ; miR-130a present %s\n",
      pl, nrow(pp$X), nrow(x), ncol(x), length(resolve(rownames(x),"hsa-miR-130a"))>0))
  Xs[[pl]] <- x; ds[[pl]] <- dd
}
common <- Reduce(intersect, lapply(Xs, rownames))
cat("GSE37405 miRNAs common to all three platforms:", length(common), "\n")
Xc <- do.call(cbind, lapply(Xs, function(x) x[common,,drop=FALSE]))
dc <- rbindlist(ds); dc <- as.data.frame(dc)
stopifnot(identical(dc$sample, colnames(Xc)))
## RFS: recurrence if time.to.recurrence present, else censored at time.to.ok
dc$RFS_event <- as.integer(is.finite(dc$t_rec))
dc$RFS_time  <- ifelse(is.finite(dc$t_rec), dc$t_rec, dc$t_ok)*12
dc$OS_event  <- as.integer(is.finite(dc$t_dth))
dc$OS_time   <- ifelse(is.finite(dc$t_dth), dc$t_dth, dc$t_ok)*12
cat("GSE37405 total n:", nrow(dc), " RFS events:", sum(dc$RFS_event),
    " OS events:", sum(dc$OS_event),
    " usable RFS time:", sum(is.finite(dc$RFS_time)),
    " usable OS time:", sum(is.finite(dc$OS_time)), "\n")
print(table(dc$platform))
addc(list(name="GSE37405", accession="GSE37405", platform="Exiqon miRCURY LNA (GPL13703/14149/15462)",
  M=Xc, normals=NULL, pheno=dc, endpoints=c("RFS","OS"),
  notes="high-risk ER+ breast cancers, adjuvant tamoxifen mono-therapy (DBCG); three array generations, analysed with platform as a stratum"))

## =============================================================== GSE26659 ==
p <- parse_series_matrix(file.path(CA,"GSE26659_series_matrix.txt.gz"))
X <- p$X; mm <- p$meta
tt <- mm$Sample_title
istum <- grepl("^Breast cancer patient", tt)
d <- data.frame(sample=mm$geo_accession, title=tt, tumour=istum,
  age=num(pull_char(mm,"age \\(years\\)")),
  grade=num(pull_char(mm,"grade")),
  node=pull_char(mm,"lymph node \\(ln\\) status"),
  ER=pull_char(mm,"estrogen receptor \\(er\\)"),
  relapse=pull_char(mm,"relapse"), stringsAsFactors=FALSE)
d$relapse_bin <- ifelse(d$relapse=="yes",1L, ifelse(d$relapse=="no",0L,NA))
cat("GSE26659: tumours", sum(istum), " non-tumour", sum(!istum),
    " relapse yes/no:", sum(d$relapse_bin==1,na.rm=TRUE),"/",sum(d$relapse_bin==0,na.rm=TRUE), "\n")
addc(list(name="GSE26659", accession="GSE26659", platform="Agilent-019118 human miRNA 2.0 (GPL8227)",
  M=X[,istum,drop=FALSE], normals=X[,!istum,drop=FALSE], pheno=d[istum,],
  endpoints=c("relapse72m_binary"),
  notes="77 ductal carcinomas + 17 mammoplasties; only a binary 72-month relapse flag, no time -> logistic, not Cox"))

## =============================================================== GSE40525 ==
p <- parse_series_matrix(file.path(CA,"GSE40525_series_matrix.txt.gz"))
X <- p$X; mm <- p$meta
tis <- pull_char(mm,"tissue")
cat("GSE40525 tissue values:\n"); print(table(tis, useNA="ifany"))
istum <- grepl("^Breast primary tumor", tis)   ## "Peritumor breast tissue" is the normal arm
cat("GSE40525 tumour/normal split:", sum(istum), "/", sum(!istum), "\n")
d <- data.frame(sample=mm$geo_accession, tissue=tis, tumour=istum,
  ER=pull_char(mm,"estrogen receptor status"),
  HER2=pull_char(mm,"her2 status"),
  subtype=pull_char(mm,"subtype \\(of tumor or matched tumor\\)"), stringsAsFactors=FALSE)
addc(list(name="GSE40525", accession="GSE40525", platform="Agilent-019118 human miRNA 2.0 (GPL8227)",
  M=X[,istum,drop=FALSE], normals=X[,!istum,drop=FALSE], pheno=d[istum,],
  endpoints=character(0), notes="tumour + matched normal breast; no outcome data"))

## =============================================================== GSE45666 ==
p <- parse_series_matrix(file.path(CA,"GSE45666_series_matrix.txt.gz"))
an <- read_mirna_annot(file.path(CA,"GPL14767_self.txt"))
X <- collapse_mirna(p$X, an, "GSE45666")
## deposited on a raw (background-subtracted) intensity scale -> log2, floored at 1
cat("GSE45666 raw value range:", signif(range(X, na.rm=TRUE),4), "-> log2(pmax(x,1))\n")
X <- log2(pmax(X, 1))
cat("GSE45666 log2 value range:", signif(range(X, na.rm=TRUE),4), "\n")
mm <- p$meta
tis <- pull_char(mm,"tissue type")
cat("GSE45666 tissue types:\n"); print(table(tis, useNA="ifany"))
cat("GSE45666 probes ->", nrow(X), "unique hsa- miRNAs\n")
istum <- grepl("^Tumor", tis)
d <- data.frame(sample=mm$geo_accession, tissue=tis, tumour=istum,
  patient=pull_char(mm,"patient id"),
  ER=pull_char(mm,"er"), grade=num(pull_char(mm,"grade")),
  stage=pull_char(mm,"Stage"), stringsAsFactors=FALSE)
addc(list(name="GSE45666", accession="GSE45666", platform="Agilent-021827 human miRNA G4470C (GPL14767)",
  M=X[,istum,drop=FALSE], normals=X[,!istum,drop=FALSE], pheno=d[istum,],
  endpoints=character(0), notes="tumour + adjacent normal; no outcome data"))

## ---- inventory ------------------------------------------------------------
inv <- rbindlist(lapply(CO, function(o){
  r130 <- resolve(rownames(o$M),"hsa-miR-130a")
  r29  <- resolve(rownames(o$M),"hsa-miR-29a")
  data.table(cohort=o$name, accession=o$accession, platform=o$platform,
    n_tumour=ncol(o$M), n_normal=if(is.null(o$normals)) 0L else ncol(o$normals),
    n_miRNA=nrow(o$M), endpoints=paste(o$endpoints, collapse=";"),
    mir130a_row=paste(r130, collapse=";"), mir29a_row=paste(r29, collapse=";"),
    notes=o$notes)
}))
print(inv[, .(cohort,n_tumour,n_normal,n_miRNA,endpoints,mir130a_row)])
fwrite(inv, file.path(OUT,"multicohort_mirna_inventory.csv"))
saveRDS(CO, file.path(CA,"mirna_cohorts.rds"))
cat("\nSaved", length(CO), "miRNA cohorts\n")
