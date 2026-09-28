## ==========================================================================
## M4_build_mrna_cohorts.R
## Take the 8 cohorts assembled by the parallel v2 workflow
## (cache/newcohorts/cohorts_base.rds), replace the 173-gene SCAN-B slice with
## the 1,897-gene slice streamed in M2, and attach ER / subtype / covariates.
## Saves cache/v4/multicohort/mrna_cohorts.rds
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
CO <- readRDS("cache/newcohorts/cohorts_base.rds")
cat("cohorts loaded:", paste(names(CO), collapse=", "), "\n")

## ---------------- SCAN-B: swap in the wider gene slice ---------------------
dt <- fread(file.path(CA,"GSE96058_subset_v4.csv"), header=TRUE)
sym <- as.character(dt[[1]]); X <- as.matrix(dt[,-1,with=FALSE]); rownames(X) <- sym
cat("SCAN-B wide slice:", nrow(X), "genes x", ncol(X), "columns\n")
isrep <- grepl("repl$", colnames(X))
cat("technical replicate columns dropped:", sum(isrep), "\n")
X <- X[, !isrep, drop=FALSE]
stopifnot(!any(duplicated(rownames(X))), !any(duplicated(colnames(X))))
old <- CO$GSE96058$X
cat("previous SCAN-B slice:", nrow(old), "genes;", "columns identical:",
    identical(colnames(old), colnames(X)), "\n")
sh <- intersect(rownames(old), rownames(X))
cat("genes shared with the previous slice:", length(sh),
    " max |difference| on shared values:",
    max(abs(old[sh, colnames(X), drop=FALSE] - X[sh,,drop=FALSE]), na.rm=TRUE), "\n")
CO$GSE96058$X <- X
stopifnot(identical(CO$GSE96058$pheno$sample, colnames(X)))

## ---------------- covariate / subtype harmonisation ------------------------
norm_er <- function(v){
  v <- tolower(trimws(as.character(v)))
  out <- rep(NA_character_, length(v))
  out[v %in% c("1","pos","positive","er+","p")] <- "Positive"
  out[v %in% c("0","neg","negative","er-","n")] <- "Negative"
  out
}
harm <- function(nm){
  ph <- CO[[nm]]$pheno
  for(k in c("ER","grade","size","node","age","subtype")) if(is.null(ph[[k]])) ph[[k]] <- NA
  ph$cohort <- nm
  ph
}

## TCGA -- PAM50 + ER from the Xena clinical matrix
clin <- fread("data/brca_clinicalMatrix.tsv", sep="\t", quote="")
sidcol <- names(clin)[1]
ph <- CO$TCGA_BRCA$pheno
m <- match(ph$sample, clin[[sidcol]])
ph$subtype <- clin$PAM50Call_RNAseq[m]
ph$ER <- norm_er(clin$ER_Status_nature2012[m])
ph$grade <- NA_real_; ph$size <- NA_real_; ph$node <- NA_character_
sv <- as.data.frame(readRDS("data/brca_survival_clean.rds"))
ms <- match(ph$sample, sv$sample_id)
ph$PFI_time <- sv$PFI.time[ms]/30.4375; ph$PFI_event <- sv$PFI[ms]
ph$DSS_time <- sv$DSS.time[ms]/30.4375; ph$DSS_event <- sv$DSS[ms]
CO$TCGA_BRCA$pheno <- ph
cat("\nTCGA subtype:\n"); print(table(ph$subtype, useNA="ifany"))
cat("TCGA ER:\n"); print(table(ph$ER, useNA="ifany"))

## METABRIC -- CLAUDIN_SUBTYPE + ER_STATUS
mb <- readRDS("data/metabric.rds"); cl <- as.data.frame(mb$cl)
key <- if(all(CO$METABRIC$pheno$sample %in% cl$PATIENT_ID)) "PATIENT_ID" else "SAMPLE_ID"
m <- match(CO$METABRIC$pheno$sample, cl[[key]])
CO$METABRIC$pheno$subtype <- cl$CLAUDIN_SUBTYPE[m]
CO$METABRIC$pheno$ER <- norm_er(cl$ER_STATUS[m])
CO$METABRIC$pheno$size <- if(!is.null(cl$TUMOR_SIZE)) cl$TUMOR_SIZE[m] else NA
CO$METABRIC$pheno$node <- if(!is.null(cl$LYMPH_NODES_EXAMINED_POSITIVE)) cl$LYMPH_NODES_EXAMINED_POSITIVE[m] else NA
cat("\nMETABRIC subtype:\n"); print(table(CO$METABRIC$pheno$subtype, useNA="ifany"))
cat("METABRIC ER:\n"); print(table(CO$METABRIC$pheno$ER, useNA="ifany"))

## SCAN-B
CO$GSE96058$pheno$ER <- norm_er(CO$GSE96058$pheno$ER)
CO$GSE96058$pheno$subtype <- CO$GSE96058$pheno$pam50
cat("\nSCAN-B subtype:\n"); print(table(CO$GSE96058$pheno$subtype, useNA="ifany"))
cat("SCAN-B ER:\n"); print(table(CO$GSE96058$pheno$ER, useNA="ifany"))

## GEO arrays
for(nm in c("GSE21653","GSE22219")) CO[[nm]]$pheno$ER <- norm_er(CO[[nm]]$pheno$ER)
CO$GSE58812$pheno$ER <- "Negative"      ## series is 100% triple-negative
CO$GSE58812$pheno$subtype <- "TNBC"
for(nm in names(CO)) CO[[nm]]$pheno <- harm(nm)
for(nm in c("GSE21653","GSE22219","GSE20685","GSE58812"))
  cat(sprintf("%s ER: %s\n", nm, paste(names(table(CO[[nm]]$pheno$ER, useNA="ifany")),
      table(CO[[nm]]$pheno$ER, useNA="ifany"), sep="=", collapse=" ")))

## ---------------- inventory ------------------------------------------------
EPMAP <- list(TCGA_BRCA=c("OS","PFI","DSS"), METABRIC=c("OS","RFS"),
  GSE96058="OS", GSE20685="OS", GSE58812=c("OS","DMFS"),
  GSE21653="DFS", GSE22219="DRFS", GTEx_breast=character(0))
inv <- rbindlist(lapply(names(CO), function(nm){
  o <- CO[[nm]]; ph <- o$pheno
  eps <- EPMAP[[nm]]
  ev <- sapply(eps, function(e) sum(ph[[paste0(e,"_event")]], na.rm=TRUE))
  data.table(cohort=nm, accession=o$accession, platform=o$platform,
    n=ncol(o$X), n_genes=nrow(o$X),
    n_normal=if(is.null(o$normals)) 0L else ncol(o$normals),
    endpoints=paste(eps, collapse=";"),
    events=paste(sprintf("%s=%d", eps, ev), collapse=";"),
    has_ER=sum(!is.na(ph$ER)), has_subtype=sum(!is.na(ph$subtype) & ph$subtype!=""),
    has_age=sum(is.finite(suppressWarnings(as.numeric(ph$age)))),
    has_grade=sum(is.finite(suppressWarnings(as.numeric(ph$grade)))),
    has_size=sum(is.finite(suppressWarnings(as.numeric(ph$size)))),
    notes=o$notes)
}))
print(inv[, .(cohort,n,n_genes,n_normal,endpoints,events,has_ER,has_subtype)])
fwrite(inv, file.path(OUT,"multicohort_mrna_inventory.csv"))
saveRDS(CO, file.path(CA,"mrna_cohorts.rds"))
cat("\nSaved", length(CO), "mRNA cohorts\n")
