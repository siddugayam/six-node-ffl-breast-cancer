## ==========================================================================
## N1b_build_geo.R -- parse the GEO array series into standardised cohorts,
## collapse probes to gene symbols, extract survival + clinical covariates.
## ==========================================================================
suppressPackageStartupMessages({library(data.table); library(matrixStats)})
setwd("/path/to/revision")
source("scripts/v2/N0_geo_utils.R")
CA <- "cache/newcohorts"
cohorts <- readRDS(file.path(CA,"cohorts_base.rds"))

gpl570 <- read_gpl_annot(file.path(CA,"GPL570.annot.gz"))
cat("GPL570 annotation rows:", nrow(gpl570),
    " with a symbol:", sum(gpl570$Symbol != "" & !is.na(gpl570$Symbol)), "\n")

num <- function(v) suppressWarnings(as.numeric(v))

PRIORITY <- local({
  h <- read.csv("results/network_topology_hubs.csv", stringsAsFactors=FALSE)
  unique(c(h$name[h$is_hub & h$type!="miRNA"],
           readLines(file.path(CA,"cafA_signature_genes.txt")),
           "DCN","LUM","FAP","THY1","COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1"))
})
cat("priority symbol set for probe collapsing:", length(PRIORITY), "genes\n")

build <- function(gse, plat_label, annot, pheno_fun, notes){
  fp <- file.path(CA, paste0(gse,"_series_matrix.txt.gz"))
  if(!file.exists(fp)){ cat("!! missing", fp, "\n"); return(NULL) }
  p <- tryCatch(parse_series_matrix(fp), error=function(e){cat("!! parse fail",gse,conditionMessage(e),"\n"); NULL})
  if(is.null(p)) return(NULL)
  cat("\n###", gse, ": probes =", nrow(p$X), " samples =", ncol(p$X),
      " platform(s) =", paste(unique(p$platform), collapse=","), "\n")
  X <- p$X
  ## log2 if the matrix is clearly on a linear scale
  mx <- max(X, na.rm=TRUE)
  if(mx > 100){ X[X < 1] <- 1; X <- log2(X); cat("  linear scale detected (max=",
      round(mx,1), ") -> log2 transformed\n") } else cat("  already log-scale (max=", round(mx,2), ")\n")
  X <- X[rowSums(is.na(X)) < 0.2*ncol(X), , drop=FALSE]
  Y <- collapse_to_symbol(X, annot, priority=PRIORITY)
  cat("  probes ->", nrow(Y), "unique gene symbols\n")
  ph <- pheno_fun(p$meta)
  ph$sample <- colnames(Y)
  stopifnot(nrow(ph)==ncol(Y))
  list(name=gse, accession=gse, platform=plat_label, X=Y, pheno=ph, notes=notes)
}

## --------------------------------------------------------------- GSE21653
ph21653 <- function(meta){
  d <- data.frame(
    age   = num(pull_char(meta,"age at diagnosis")),
    grade = num(pull_char(meta,"sbr grade")),
    ER    = pull_char(meta,"er ihc"),
    node  = pull_char(meta,"pn"),
    size  = pull_char(meta,"pt"),
    DFS_event = num(pull_char(meta,"dfs evt")),
    DFS_time  = num(pull_char(meta,"dfs time \\(months\\)")),
    stringsAsFactors=FALSE)
  d
}
## --------------------------------------------------------------- GSE58812
ph58812 <- function(meta){
  data.frame(
    age = num(pull_char(meta,"age at diag")),
    OS_event  = num(pull_char(meta,"death")),
    OS_time   = num(pull_char(meta,"os \\(days\\)"))/30.4375,
    DMFS_event= num(pull_char(meta,"meta")),
    DMFS_time = num(pull_char(meta,"mfs \\(days\\)"))/30.4375,
    stringsAsFactors=FALSE)
}

specs <- list(
  GSE21653 = list(plat="Affymetrix HG-U133 Plus 2.0 (GPL570)", annot=gpl570,
                  f=ph21653, notes="266 primary breast tumours (Marseille); DFS"),
  GSE58812 = list(plat="Affymetrix HG-U133 Plus 2.0 (GPL570)", annot=gpl570,
                  f=ph58812, notes="107 triple-negative breast tumours (Nantes); OS + DMFS")
)
for(g in names(specs)){
  s <- specs[[g]]
  o <- build(g, s$plat, s$annot, s$f, s$notes)
  if(!is.null(o)){
    cohorts[[g]] <- o
    cat("  pheno summary:\n"); print(summary(o$pheno[sapply(o$pheno,is.numeric)]))
  }
}
saveRDS(cohorts, file.path(CA,"cohorts_base.rds"), compress=FALSE)
cat("\ncohorts now:", paste(names(cohorts), collapse=", "), "\n")
