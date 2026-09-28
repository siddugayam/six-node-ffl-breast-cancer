## ==========================================================================
## N1e_build_scanb.R -- add SCAN-B / GSE96058 (3,273 tumours, RNA-seq) to the
## cohort list, from the streamed gene subset + the two series matrices.
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
setwd("/path/to/revision")
source("scripts/07_expression_validation/external_cohorts/N0_geo_utils.R")
CA <- "cache/newcohorts"
cohorts <- readRDS(file.path(CA,"cohorts_base.rds"))
num <- function(v) suppressWarnings(as.numeric(v))

## ---- expression subset ----------------------------------------------------
dt <- fread(file.path(CA,"GSE96058_subset.csv"), header=TRUE)
sym <- as.character(dt[[1]])
X <- as.matrix(dt[, -1, with=FALSE]); rownames(X) <- sym
cat("SCAN-B subset:", nrow(X), "genes x", ncol(X), "columns\n")
isrep <- grepl("repl$", colnames(X))
cat("technical replicate columns dropped:", sum(isrep), "\n")
X <- X[, !isrep, drop=FALSE]
stopifnot(!any(duplicated(colnames(X))), !any(duplicated(rownames(X))))
cat("unique tumours:", ncol(X), "\n")
cat("value range:", round(range(X, na.rm=TRUE),3), "(FPKM, log2(x+0.1) transformed by the depositors)\n")

## ---- metadata from both series matrices -----------------------------------
meta <- NULL
for(f in c("GSE96058_sm_GPL11154.txt.gz","GSE96058_sm_GPL18573.txt.gz")){
  p <- parse_series_matrix(file.path(CA,f))
  m <- p$meta
  m$title <- m$Sample_title
  d <- data.frame(
    title = m$Sample_title, gsm = m$geo_accession,
    instrument = pull_char(m,"instrument model"),
    age   = num(pull_char(m,"age at diagnosis")),
    size  = num(pull_char(m,"tumor size")),
    node  = pull_char(m,"lymph node status"),
    ER    = pull_char(m,"er status"),
    PGR   = pull_char(m,"pgr status"),
    HER2  = pull_char(m,"her2 status"),
    grade = pull_char(m,"nhg"),
    pam50 = pull_char(m,"pam50 subtype"),
    OS_time  = num(pull_char(m,"overall survival days"))/30.4375,
    OS_event = num(pull_char(m,"overall survival event")),
    endocrine = pull_char(m,"endocrine treated"),
    chemo     = pull_char(m,"chemo treated"),
    stringsAsFactors=FALSE)
  cat(f, ": ", nrow(d), " samples\n", sep="")
  meta <- rbind(meta, d)
}
cat("total series-matrix samples:", nrow(meta), "\n")
meta$grade[meta$grade %in% c("NA","")] <- NA
meta$grade <- num(sub("^G","", meta$grade))
meta$ER[meta$ER=="NA"] <- NA; meta$node[meta$node=="NA"] <- NA
mm <- match(colnames(X), meta$title)
cat("expression columns matched to series-matrix titles:", sum(!is.na(mm)), "of", ncol(X), "\n")
stopifnot(all(!is.na(mm)))
ph <- meta[mm, ]; ph$sample <- colnames(X); rownames(ph) <- NULL

cat("OS available:", sum(is.finite(ph$OS_time) & !is.na(ph$OS_event)),
    " events:", sum(ph$OS_event, na.rm=TRUE), "\n")
print(summary(ph[, c("age","size","grade","OS_time")]))
print(table(ph$ER, useNA="ifany"))

cohorts$GSE96058 <- list(name="GSE96058", accession="GSE96058 (SCAN-B)",
  platform="Illumina HiSeq 2000 / NextSeq 500 RNA-seq (log2 FPKM)",
  X=X, pheno=ph,
  notes="3,273 primary breast tumours, population-based Swedish SCAN-B cohort; OS")
saveRDS(cohorts, file.path(CA,"cohorts_base.rds"), compress=FALSE)
cat("\ncohorts now:", paste(names(cohorts), collapse=", "), "\n")
