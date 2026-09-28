## ==========================================================================
## M12 -- one consolidated cohort inventory for the multi-cohort arm:
## accession, n, platform, survival availability, miRNA availability.
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
CA <- "cache/v4/multicohort"; OUT <- "results/v4"
source("scripts/07_expression_validation/external_cohorts/N0_geo_utils.R")
CO <- readRDS(file.path(CA,"mrna_cohorts.rds"))
MI <- readRDS(file.path(CA,"mirna_cohorts.rds"))

## verify GSE22216 is the miRNA arm of GSE22220
a <- system("zcat cache/v4/multicohort/GSE22216_series_matrix.txt.gz | grep '^!Sample_geo_accession'", intern=TRUE)
b <- system("zcat cache/external2/GSE22220-GPL8178_series_matrix.txt.gz | grep '^!Sample_geo_accession'", intern=TRUE)
ga <- gsub('"','',strsplit(a,"\t")[[1]][-1]); gb <- gsub('"','',strsplit(b,"\t")[[1]][-1])
cat("GSE22216 GSMs:", length(ga), " GSE22220-GPL8178 GSMs:", length(gb),
    " identical sets:", setequal(ga,gb), "\n")

EPMAP <- list(TCGA_BRCA=c("OS","PFI","DSS"), METABRIC=c("OS","RFS"), GSE96058="OS",
  GSE20685="OS", GSE58812=c("OS","DMFS"), GSE21653="DFS", GSE22219="DRFS",
  GTEx_breast=character(0))
rows <- list()
for(nm in names(CO)){
  o <- CO[[nm]]; ph <- o$pheno; eps <- EPMAP[[nm]]
  ev <- if(length(eps)) sapply(eps, function(e) sum(ph[[paste0(e,"_event")]], na.rm=TRUE)) else integer(0)
  rows[[length(rows)+1]] <- data.table(cohort=nm, accession=o$accession, layer="mRNA",
    platform=o$platform, n_tumour=ncol(o$X),
    n_normal=if(is.null(o$normals)) 0L else ncol(o$normals),
    n_features=nrow(o$X),
    survival=if(length(eps)) paste(sprintf("%s (%d events)", eps, ev), collapse="; ") else "none",
    miRNA_assay=ifelse(nm %in% c("TCGA_BRCA","GSE22219"),
      ifelse(nm=="TCGA_BRCA","yes - matched TCGA miRNA-seq",
             "yes - matched GSE22216/GSE22220 miRNA array (205 of 210 patients linked)"), "no"),
    used_for=paste(c(if(length(eps)) "C/D/E survival meta" else NULL,
      if(nm=="GTEx_breast") "non-tumour baseline; A3 cross-study hub DE" else NULL,
      if(nm=="TCGA_BRCA") "reference cohort; target set derived here" else NULL), collapse="; "),
    notes=o$notes)
}
for(nm in names(MI)){
  o <- MI[[nm]]; ph <- o$pheno; eps <- o$endpoints
  ev <- sapply(eps, function(e){
    if(e=="relapse72m_binary") return(sum(ph$relapse_bin, na.rm=TRUE))
    sum(ph[[paste0(e,"_event")]], na.rm=TRUE)})
  rows[[length(rows)+1]] <- data.table(cohort=nm, accession=o$accession, layer="miRNA",
    platform=o$platform, n_tumour=ncol(o$M),
    n_normal=if(is.null(o$normals)) 0L else ncol(o$normals),
    n_features=nrow(o$M),
    survival=if(length(eps)) paste(sprintf("%s (%d events)", eps, ev), collapse="; ") else "none",
    miRNA_assay="yes",
    used_for=paste(c(if(length(eps) && !identical(eps,"relapse72m_binary")) "B Cox + HR meta" else NULL,
      if(identical(eps,"relapse72m_binary")) "B logistic only (no time variable)" else NULL,
      if(!is.null(o$normals) && ncol(o$normals)>0) "A2 + B tumour-vs-normal" else NULL,
      if(nm %in% c("GSE19783","TCGA_BRCA")) "C0 target-score validation (matched mRNA)" else NULL),
      collapse="; "),
    notes=o$notes)
}
INV <- rbindlist(rows)
INV <- INV[order(layer, -n_tumour)]
print(INV[, .(cohort,layer,n_tumour,n_normal,n_features,survival,miRNA_assay)])
fwrite(INV, file.path(OUT,"multicohort_COHORT_INVENTORY.csv"))
cat("\ntotal tumours, mRNA layer:", sum(INV[layer=="mRNA" & cohort!="GTEx_breast"]$n_tumour), "\n")
cat("total tumours, miRNA layer:", sum(INV[layer=="miRNA"]$n_tumour), "\n")
cat("miRNA cohorts with a time-to-event outcome:",
    paste(INV[layer=="miRNA" & survival!="none" & !grepl("relapse72m", survival)]$cohort, collapse=", "), "\n")
