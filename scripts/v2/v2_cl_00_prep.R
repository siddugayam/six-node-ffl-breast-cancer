#!/usr/bin/env Rscript
## Cell-line "stroma-free test": data preparation
## DepMap Public 24Q4 (figshare article 27993248), md5-verified.
suppressPackageStartupMessages({library(data.table)})
options(warn=1)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4")
OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

## ---------- Model annotation ----------
mod <- fread(file.path(DM,"Model.csv"))
msg("Model.csv rows:", nrow(mod))
breast <- mod[OncotreeLineage=="Breast"]
msg("Breast models:", nrow(breast),
    "| cell lines:", sum(breast$ModelType=="Cell Line"),
    "| organoids:", sum(breast$ModelType=="Organoid"),
    "| non-cancerous:", sum(breast$OncotreePrimaryDisease=="Non-Cancerous"))

## ---------- Expression ----------
msg("reading OmicsExpressionProteinCodingGenesTPMLogp1.csv ...")
ex <- fread(file.path(DM,"OmicsExpressionProteinCodingGenesTPMLogp1.csv"))
setnames(ex, 1, "ModelID")
msg("expression matrix:", nrow(ex), "models x", ncol(ex)-1, "protein-coding genes")
gcols <- setdiff(names(ex),"ModelID")
gsym  <- sub(" \\(\\d+\\)$","", gcols)
stopifnot(!any(duplicated(gsym)))
## HARD ASSERTION: feature names must be gene symbols, not indices
stopifnot(all(grepl("^[A-Za-z]", gsym)))
stopifnot(all(c("COL1A1","COL3A1","ETS1","NFKB1","RELA","SP1","EPCAM","ACTB") %in% gsym))
msg("gene-symbol assertion passed; e.g.", paste(head(gsym,5),collapse=", "))

E <- as.matrix(ex[, gcols, with=FALSE]); rownames(E) <- ex$ModelID; colnames(E) <- gsym
saveRDS(list(E=E, models=ex$ModelID, genes=gsym), file.path(DM,"expr_24Q4.rds"))
msg("saved expr_24Q4.rds ; range", paste(round(range(E,na.rm=TRUE),3),collapse=" .. "))

## ---------- CCLE miRNA (Nanostring, 2018-11-03) ----------
gct <- readLines(file.path(DM,"CCLE_miRNA_20181103.gct"))
dims <- strsplit(gct[2],"\t")[[1]]
msg("CCLE miRNA gct declares", dims[1], "features x", dims[2], "cell lines")
hdr <- strsplit(gct[3],"\t")[[1]]
body <- do.call(rbind, strsplit(gct[-(1:3)], "\t"))
mirn <- body[,1]; mdesc <- body[,2]
M <- matrix(as.numeric(body[,-(1:2)]), nrow=nrow(body))
rownames(M) <- mirn; colnames(M) <- hdr[-(1:2)]
msg("miRNA matrix:", nrow(M), "x", ncol(M), "| value range",
    paste(round(range(M,na.rm=TRUE),3),collapse=" .. "))
saveRDS(list(M=M, desc=mdesc), file.path(DM,"ccle_mirna.rds"))
msg("miRNA rownames example:", paste(head(mirn,4),collapse=" | "))
msg("done")
