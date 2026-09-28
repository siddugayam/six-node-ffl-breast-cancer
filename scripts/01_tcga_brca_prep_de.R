#!/usr/bin/env Rscript
## TCGA-BRCA expression matrices + differential expression (tumour vs normal)
## Shared substrate for downstream FFL validation analyses.
suppressPackageStartupMessages({library(data.table); library(limma)})

ML   <- "/path/to/home/Desktop/DD/R_GPR/ML"
REV  <- "/path/to/revision"
DATA <- file.path(REV,"data"); RES <- file.path(REV,"results"); LOG <- file.path(REV,"logs")
dir.create(DATA,showWarnings=FALSE,recursive=TRUE); dir.create(RES,showWarnings=FALSE,recursive=TRUE)
dir.create(LOG,showWarnings=FALSE,recursive=TRUE)
logf <- file.path(LOG,"tcga_brca_prep.log")
say <- function(...){ msg <- paste0(format(Sys.time(),"%H:%M:%S"),"  ",paste0(...,collapse=""))
  cat(msg,"\n"); cat(msg,"\n",file=logf,append=TRUE) }
cat("=== run",format(Sys.time()),"===\n",file=logf,append=TRUE)

GENE_F <- file.path(ML,"EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
MIR_F  <- file.path(ML,"pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz")
PHE_F  <- file.path(ML,"TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
SUR_F  <- file.path(ML,"Survival_SupplementalTable_S1_20171025_xena_sp.txt")

## ---- 1. phenotype -----------------------------------------------------------
phe <- fread(cmd=paste("zcat",shQuote(PHE_F)), sep="\t", header=TRUE)
say("phenotype rows: ",nrow(phe)," cols: ",paste(names(phe),collapse=","))
brca <- phe[`_primary_disease`=="breast invasive carcinoma"]
say("BRCA phenotype rows: ",nrow(brca))
print(table(brca$sample_type))
brca_tn <- brca[sample_type %in% c("Primary Tumor","Solid Tissue Normal")]
say("BRCA Primary Tumor+Solid Tissue Normal in phenotype: ",nrow(brca_tn))

## ---- 2. gene expression, subset columns EARLY -------------------------------
hdr_g <- names(fread(cmd=paste("zcat",shQuote(GENE_F)), sep="\t", nrows=0))
say("gene matrix header cols: ",length(hdr_g)," first: ",hdr_g[1])
keep_g <- which(hdr_g %in% brca_tn$sample)
keep_g <- keep_g[!duplicated(hdr_g[keep_g])]        # drop duplicated barcodes if any
say("gene columns matched to BRCA T/N: ",length(keep_g),
    " (duplicated barcodes among them: ",sum(duplicated(hdr_g[which(hdr_g %in% brca_tn$sample)])),")")
g <- fread(cmd=paste("zcat",shQuote(GENE_F)), sep="\t", select=c(1L,keep_g), header=TRUE)
gid <- g[[1]]; g[,1:=NULL]
gm <- as.matrix(g); rownames(gm) <- gid; rm(g); gc()
say("gene matrix dims (raw): ",nrow(gm)," x ",ncol(gm),
    " ; duplicated gene ids: ",sum(duplicated(rownames(gm))),
    " -> ",paste(unique(rownames(gm)[duplicated(rownames(gm))]),collapse=","))
## collapse duplicated gene symbols by mean (limma/topTable silently drops ALL
## rownames if any are duplicated, so this must be done before DE)
if(any(duplicated(rownames(gm)))){
  dup <- unique(rownames(gm)[duplicated(rownames(gm))])
  keep <- !(rownames(gm) %in% dup)
  extra <- t(sapply(dup, function(s) colMeans(gm[rownames(gm)==s,,drop=FALSE],na.rm=TRUE)))
  rownames(extra) <- dup
  gm <- rbind(gm[keep,,drop=FALSE], extra)
  gm <- gm[order(rownames(gm)),,drop=FALSE]
  say("collapsed ",length(dup)," duplicated gene symbol(s) by mean; new dims: ",
      nrow(gm)," x ",ncol(gm))
}
stopifnot(!any(duplicated(rownames(gm))))

## ---- 3. miRNA expression ----------------------------------------------------
hdr_m <- names(fread(cmd=paste("zcat",shQuote(MIR_F)), sep="\t", nrows=0))
keep_m <- which(hdr_m %in% brca_tn$sample)
keep_m <- keep_m[!duplicated(hdr_m[keep_m])]
say("miRNA columns matched to BRCA T/N: ",length(keep_m))
m <- fread(cmd=paste("zcat",shQuote(MIR_F)), sep="\t", select=c(1L,keep_m), header=TRUE)
mid <- m[[1]]; m[,1:=NULL]
mm <- as.matrix(m); rownames(mm) <- mid; rm(m); gc()
say("miRNA matrix dims: ",nrow(mm)," x ",ncol(mm),
    " ; duplicated miRNA ids: ",sum(duplicated(rownames(mm))))
stopifnot(!any(duplicated(rownames(mm))))

## ---- 4. pheno / survival tables --------------------------------------------
allsamp <- union(colnames(gm), colnames(mm))
brca_pheno <- as.data.frame(brca_tn[sample %in% allsamp, .(sample, sample_type)])
sur <- fread(SUR_F, sep="\t", header=TRUE)
brca_surv <- as.data.frame(sur[sample %in% allsamp])
say("brca_pheno rows: ",nrow(brca_pheno),"  brca_survival rows: ",nrow(brca_surv))

ng_t <- sum(colnames(gm) %in% brca_tn[sample_type=="Primary Tumor",sample])
ng_n <- sum(colnames(gm) %in% brca_tn[sample_type=="Solid Tissue Normal",sample])
nm_t <- sum(colnames(mm) %in% brca_tn[sample_type=="Primary Tumor",sample])
nm_n <- sum(colnames(mm) %in% brca_tn[sample_type=="Solid Tissue Normal",sample])
say("GENE  n tumour=",ng_t,"  n normal=",ng_n)
say("miRNA n tumour=",nm_t,"  n normal=",nm_n)

saveRDS(gm, file.path(DATA,"brca_gene_expr.rds"))
saveRDS(mm, file.path(DATA,"brca_mirna_expr.rds"))
saveRDS(brca_pheno, file.path(DATA,"brca_pheno.rds"))
saveRDS(brca_surv, file.path(DATA,"brca_survival.rds"))

## ---- 5. miRNA canonical id map ---------------------------------------------
canon_mir <- function(x){
  y <- sub("-(3p|5p)$","",x)                       # strip arm suffix
  y <- sub("^(hsa-(miR|mir|let)-[0-9]+[a-zA-Z]*)-[0-9]+$","\\1",y)  # strip paralogue -1/-2/-3
  y <- sub("^hsa-mir-","hsa-miR-",y)               # upper-case miR
  y
}
map <- data.table(raw_id=rownames(mm), canonical_id=canon_mir(rownames(mm)))
fwrite(map, file.path(DATA,"mirna_id_map.tsv"), sep="\t")
say("mirna_id_map.tsv rows: ",nrow(map),"  unique canonical ids: ",uniqueN(map$canonical_id),
    "  canonical ids from >1 row: ",sum(table(map$canonical_id)>1))

## collapse by mean
idx  <- split(seq_len(nrow(mm)), map$canonical_id)
mmc  <- do.call(rbind, lapply(idx, function(i) if(length(i)==1) mm[i,] else colMeans(mm[i,,drop=FALSE],na.rm=TRUE)))
rownames(mmc) <- names(idx)
say("canonical-collapsed miRNA matrix dims: ",nrow(mmc)," x ",ncol(mmc))
saveRDS(mmc, file.path(DATA,"brca_mirna_expr_canonical.rds"))

## ---- 6. limma DE ------------------------------------------------------------
run_limma <- function(X, tumour, normal, label){
  X <- X[, c(tumour, normal), drop=FALSE]
  ok <- rowSums(!is.na(X)) >= 3 & apply(X,1,function(v) {v<-v[!is.na(v)]; length(v)>2 && sd(v)>0})
  say(label,": rows in matrix=",nrow(X),"  rows tested (non-constant, >=3 obs)=",sum(ok))
  X <- X[ok,,drop=FALSE]
  grp <- factor(c(rep("Tumor",length(tumour)), rep("Normal",length(normal))), levels=c("Normal","Tumor"))
  des <- model.matrix(~grp)
  fit <- eBayes(lmFit(X, des))
  tt  <- topTable(fit, coef=2, number=Inf, sort.by="P")
  data.table(feature=rownames(tt), logFC=tt$logFC, AveExpr=tt$AveExpr, t=tt$t,
             P.Value=tt$P.Value, adj.P.Val=tt$adj.P.Val)
}
gt <- intersect(colnames(gm), brca_tn[sample_type=="Primary Tumor",sample])
gn <- intersect(colnames(gm), brca_tn[sample_type=="Solid Tissue Normal",sample])
mt <- intersect(colnames(mmc), brca_tn[sample_type=="Primary Tumor",sample])
mn <- intersect(colnames(mmc), brca_tn[sample_type=="Solid Tissue Normal",sample])

de_g <- run_limma(gm, gt, gn, "GENES")
de_m <- run_limma(mmc, mt, mn, "miRNAs")
fwrite(de_g, file.path(RES,"BRCA_DEX_genes.csv"))
fwrite(de_m, file.path(RES,"BRCA_DEX_mirnas.csv"))

sig <- function(d) sum(abs(d$logFC)>1 & d$adj.P.Val<0.05)
say("GENES tested=",nrow(de_g),"  significant |logFC|>1 & adj.P<0.05 = ",sig(de_g))
say("miRNAs tested=",nrow(de_m)," significant |logFC|>1 & adj.P<0.05 = ",sig(de_m))

allnodes <- rbind(
  data.table(Gene=de_g$feature, logFC=de_g$logFC, adj.P.Val=de_g$adj.P.Val, class="gene"),
  data.table(Gene=de_m$feature, logFC=de_m$logFC, adj.P.Val=de_m$adj.P.Val, class="miRNA"))
fwrite(allnodes, file.path(RES,"BRCA_DEX_ALL_nodes.csv"))
say("BRCA_DEX_ALL_nodes.csv rows: ",nrow(allnodes))

## ---- 7. network node matching ----------------------------------------------
nodes <- fread(file.path(DATA,"canonical_nodes.tsv"), sep="\t", header=TRUE)
say("canonical nodes: ",nrow(nodes))
gset <- de_g$feature; mset <- de_m$feature
matched <- logical(nrow(nodes)); how <- character(nrow(nodes))
for(i in seq_len(nrow(nodes))){
  nm <- nodes$name[i]; ty <- nodes$type[i]
  if(ty=="miRNA"){
    if(nm %in% mset){ matched[i] <- TRUE; how[i] <- "exact" }
    else {
      al <- unlist(strsplit(as.character(nodes$aliases[i]),"[;,|]"))
      al <- canon_mir(trimws(al))
      hit <- al[al %in% mset]
      if(length(hit)) { matched[i] <- TRUE; how[i] <- "alias" }
    }
  } else {
    if(nm %in% gset){ matched[i] <- TRUE; how[i] <- "exact" }
    else {
      al <- trimws(unlist(strsplit(as.character(nodes$aliases[i]),"[;,|]")))
      hit <- al[al %in% gset]
      if(length(hit)) { matched[i] <- TRUE; how[i] <- "alias" }
    }
  }
}
## relaxed family match for miRNAs only (ignore trailing letter of the family,
## e.g. network hsa-miR-103 vs data hsa-miR-103a) - reported separately, NOT the headline rate
fam <- function(x) sub("^(hsa-(?:miR|let)-[0-9]+)[a-z]?$","\\1", x, perl=TRUE)
mset_fam <- fam(mset)
relaxed <- matched
for(i in which(!matched & nodes$type=="miRNA"))
  if(fam(nodes$name[i]) %in% mset_fam) relaxed[i] <- TRUE
nodes[, matched := matched][, match_mode := how][, matched_relaxed := relaxed]
fwrite(nodes, file.path(RES,"node_DE_match.tsv"), sep="\t")
tab <- nodes[, .(n=.N, matched=sum(matched)), by=type]
tab[, rate := round(100*matched/n,1)]
say("NODE MATCHING:"); print(tab); capture.output(print(tab), file=logf, append=TRUE)
say("overall matched ",sum(matched),"/",nrow(nodes)," = ",round(100*mean(matched),1),"%")
say("miRNA nodes matched with relaxed family rule: ",
    sum(relaxed & nodes$type=="miRNA"),"/",sum(nodes$type=="miRNA"))
unm <- nodes[type=="miRNA" & !matched, name]
say("unmatched network miRNAs (",length(unm),"): ",paste(head(unm,20),collapse=", "))
writeLines(unm, file.path(RES,"unmatched_network_mirnas.txt"))
unmg <- nodes[type!="miRNA" & !matched, name]
say("unmatched network genes/TFs (",length(unmg),"): ",paste(head(unmg,30),collapse=", "))
say("DONE")
