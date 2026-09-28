## N1c: add GSE20685 (GPL570) and GSE22219 (GPL6098) to the cohort list
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
setwd("/path/to/revision")
source("scripts/v2/N0_geo_utils.R")
CA <- "cache/newcohorts"
cohorts <- readRDS(file.path(CA,"cohorts_base.rds"))
num <- function(v) suppressWarnings(as.numeric(v))

gpl570 <- read_gpl_annot(file.path(CA,"GPL570.annot.gz"))
gpl6098 <- read_gpl_soft(file.path(CA,"GPL6098.soft"), sym_regex="^Symbol$")
cat("GPL6098 annotation rows:", nrow(gpl6098),
    " unique symbols:", length(unique(gpl6098$Symbol)), "\n")

PRIORITY <- local({
  h <- read.csv("results/network_topology_hubs.csv", stringsAsFactors=FALSE)
  unique(c(h$name[h$is_hub & h$type!="miRNA"],
           readLines(file.path(CA,"cafA_signature_genes.txt")),
           "DCN","LUM","FAP","THY1","COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1"))
})
cat("priority symbol set for probe collapsing:", length(PRIORITY), "genes\n")

build <- function(gse, plat_label, annot, pheno_fun, notes){
  fp <- file.path(CA, paste0(gse,"_series_matrix.txt.gz"))
  p <- parse_series_matrix(fp)
  cat("\n###", gse, ": probes =", nrow(p$X), " samples =", ncol(p$X),
      " platform(s) =", paste(unique(p$platform), collapse=","), "\n")
  X <- p$X
  mx <- max(X, na.rm=TRUE)
  if(mx > 100){ X[X < 1] <- 1; X <- log2(X)
    cat("  linear scale (max=", round(mx,1), ") -> log2\n") } else
    cat("  already log-scale (max=", round(mx,2), ")\n")
  X <- X[rowSums(is.na(X)) < 0.2*ncol(X), , drop=FALSE]
  Y <- collapse_to_symbol(X, annot, priority=PRIORITY)
  cat("  probes ->", nrow(Y), "unique gene symbols\n")
  ph <- pheno_fun(p$meta); ph$sample <- colnames(Y)
  stopifnot(nrow(ph)==ncol(Y))
  list(name=gse, accession=gse, platform=plat_label, X=Y, pheno=ph, notes=notes)
}

ph20685 <- function(meta){
  data.frame(
    age = num(pull_char(meta,"age at diagnosis")),
    OS_event = num(pull_char(meta,"event_death")),
    OS_time  = num(pull_char(meta,"follow_up_duration \\(years\\)"))*12,
    metastasis_event = num(pull_char(meta,"event_metastasis")),
    stage = pull_char(meta,"t_stage"),
    node  = pull_char(meta,"n_stage"),
    stringsAsFactors=FALSE)
}
ph22219 <- function(meta){
  g <- pull_char(meta,"tumour grade"); g[g=="."] <- NA
  data.frame(
    age  = num(pull_char(meta,"patient age")),
    size = num(pull_char(meta,"tumour size")),
    node = num(pull_char(meta,"nodes involved")),
    ER   = pull_char(meta,"er status"),
    grade= num(g),
    DRFS_event = num(pull_char(meta,"distant-relapse event")),
    DRFS_time  = num(pull_char(meta,"distant-relapse free survival"))*12,
    stringsAsFactors=FALSE)
}

for(s in list(
  list(g="GSE20685", plat="Affymetrix HG-U133 Plus 2.0 (GPL570)", a=gpl570, f=ph20685,
       notes="327 primary breast tumours (Taiwan); OS from follow-up duration"),
  list(g="GSE22219", plat="Illumina humanRef-8 v1.0 beadchip (GPL6098)", a=gpl6098, f=ph22219,
       notes="216 primary breast tumours (Uppsala/UK); distant-relapse-free survival"))){
  o <- tryCatch(build(s$g, s$plat, s$a, s$f, s$notes),
                error=function(e){cat("!! FAILED", s$g, conditionMessage(e), "\n"); NULL})
  if(!is.null(o)){ cohorts[[s$g]] <- o
    print(summary(o$pheno[sapply(o$pheno, is.numeric)])) }
}
saveRDS(cohorts, file.path(CA,"cohorts_base.rds"), compress=FALSE)
cat("\ncohorts now:", paste(names(cohorts), collapse=", "), "\n")
