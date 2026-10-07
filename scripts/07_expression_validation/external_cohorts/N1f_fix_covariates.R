## ==========================================================================
## N1f_fix_covariates.R -- recodes three clinical covariates of the cohorts used
## by N2 (identified from the Cox convergence warnings of N2).
##   (1) TCGA-BRCA stage: "Stage III" is matched before "Stage II", so that
##       Stage IIIA/IIIB/IIIC patients are coded III
##       (274 patients).
##   (2) GSE21653: the literal string "NA" in node / ER / size is read as
##       missing rather than as its own factor level (`nodeNA`).
##   (3) GSE20685: stage level "1c" has n=2 and no events, giving an infinite
##       Cox coefficient; sub-stage letters are collapsed to the main stage.
## Every count printed below is computed in this run.
## ==========================================================================
setwd("/path/to/revision")
CA <- "cache/newcohorts"
co <- readRDS(file.path(CA,"cohorts_base.rds"))

## ---------------------------------------------------------------- (1) -----
sv <- as.data.frame(readRDS("data/brca_survival_clean.rds"))
ph <- co$TCGA_BRCA$pheno
st <- sv$ajcc_pathologic_tumor_stage[match(ph$sample, sv$sample_id)]
old <- ph$stage
new <- ifelse(grepl("Stage IV", st), "IV",
       ifelse(grepl("Stage III", st), "III",
       ifelse(grepl("Stage II", st),  "II",
       ifelse(grepl("Stage I",  st),  "I", NA))))
cat("(1) TCGA stage re-derived. old vs new cross-tab:\n")
print(table(old=old, new=new, useNA="ifany"))
cat("    patients moved II -> III:", sum(old=="II" & new=="III", na.rm=TRUE), "\n")
co$TCGA_BRCA$pheno$stage <- new

## ---------------------------------------------------------------- (2) -----
nastr <- function(v){
  if(!is.character(v)) return(v)
  v[trimws(v) %in% c("NA","na","N/A",".","","null","NULL","unknown","Unknown")] <- NA
  v
}
for(nm in names(co)){
  ph <- co[[nm]]$pheno
  changed <- character(0)
  for(cv in c("grade","stage","ER","node","size")){
    if(!(cv %in% names(ph))) next
    before <- sum(is.na(ph[[cv]]))
    ph[[cv]] <- nastr(ph[[cv]])
    after <- sum(is.na(ph[[cv]]))
    if(after > before) changed <- c(changed, sprintf("%s(+%d)", cv, after-before))
  }
  co[[nm]]$pheno <- ph
  if(length(changed)) cat("(2)", nm, ": literal-NA strings converted ->",
                          paste(changed, collapse=" "), "\n")
}

## ---------------------------------------------------------------- (3) -----
s <- co$GSE20685$pheno$stage
cat("(3) GSE20685 stage before:", paste(paste0(names(table(s)),"=",as.integer(table(s))), collapse=" "), "\n")
co$GSE20685$pheno$stage <- sub("^([0-9]+).*$", "\\1", s)
s2 <- co$GSE20685$pheno$stage
cat("    GSE20685 stage after :", paste(paste0(names(table(s2)),"=",as.integer(table(s2))), collapse=" "), "\n")

## ---------------------------------------------------------------- (4) -----
## METABRIC codes 12 patients as stage "0" (DCIS). In the complete-case OS set
## only 2 of them remain, and because "0" sorts first it becomes the factor
## reference level -- every stage coefficient then carries SE ~1.0 and coxph
## warns that the coefficients may be infinite. Merge stage 0 into stage I.
s0 <- co$METABRIC$pheno$stage
cat("(4) METABRIC stage before:", paste(paste0(names(table(s0)),"=",as.integer(table(s0))), collapse=" "), "\n")
s0[!is.na(s0) & s0=="0"] <- "I"
co$METABRIC$pheno$stage <- s0
cat("    METABRIC stage after :", paste(paste0(names(table(s0)),"=",as.integer(table(s0))), collapse=" "), "\n")

saveRDS(co, file.path(CA,"cohorts_base.rds"), compress=FALSE)
cat("\nre-saved cohorts_base.rds ; cohorts:", paste(names(co), collapse=", "), "\n")
