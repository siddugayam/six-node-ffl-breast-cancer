#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- part B (miRNA arm) + sensitivity analyses
## CCLE miRNA = Nanostring nCounter panel, CCLE_miRNA_20181103.gct
##   (734 miRNAs x 954 lines), from data.broadinstitute.org/ccle/
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); n<-sum(ok)
  if(n<8) return(list(rho=NA_real_,p=NA_real_,n=n))
  ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  list(rho=unname(ct$estimate),p=ct$p.value,n=n) }
pcor <- function(x,y,z){ ok<-is.finite(x)&is.finite(y)&apply(as.matrix(z),1,function(r) all(is.finite(r)))
  x<-rank(x[ok]); y<-rank(y[ok]); Z<-apply(as.matrix(z)[ok,,drop=FALSE],2,rank)
  rx<-resid(lm(x~Z)); ry<-resid(lm(y~Z)); n<-sum(ok)
  r<-cor(rx,ry); df<-n-2-ncol(as.matrix(z))
  t<-r*sqrt(df/(1-r^2)); list(rho=r,p=2*pt(-abs(t),df),n=n) }

W   <- readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod <- fread(file.path(DM,"Model.csv"))
mi  <- readRDS(file.path(DM,"ccle_mirna.rds"))
Mm  <- mi$M; rownames(Mm) <- mi$desc                    # miRNA names as rownames
msg("CCLE miRNA:", nrow(Mm),"miRNAs x",ncol(Mm),"lines; raw Nanostring counts")
Mm  <- log2(Mm + 1)                                     # log2 transform
stopifnot(all(c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c") %in% rownames(Mm)))

## map CCLE_Name -> ModelID
key <- mod[CCLEName!="" , .(CCLEName, ModelID, OncotreeLineage, ModelType,
                            OncotreePrimaryDisease, StrippedCellLineName)]
setkey(key, CCLEName)
cn  <- colnames(Mm); mid <- key[cn, ModelID]
msg("miRNA columns mapped to a DepMap ModelID:", sum(!is.na(mid)), "of", length(cn))
brs <- key[cn, OncotreeLineage]=="Breast" & !is.na(mid)
msg("breast-lineage lines in the miRNA panel:", sum(brs, na.rm=TRUE))
use <- which(!is.na(mid) & mid %in% W$breast_ca & key[cn,OncotreeLineage]=="Breast")
ids <- mid[use]
msg("MALIGNANT breast cell lines with BOTH miRNA and 24Q4 expression:", length(ids))
writeLines(paste(key[cn[use],StrippedCellLineName], ids), file.path(OUT,"celllines_mirna_lines_used.txt"))

MIR <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
TG  <- c("COL1A1","COL3A1")

## EMT / mesenchymal score (tumour-cell-intrinsic confounder, NOT a stromal score)
emt_genes <- intersect(c("VIM","CDH2","ZEB1","ZEB2","SNAI2","TWIST1","FN1","SERPINE1"), colnames(E))
zsc <- function(m) scale(m)
emt <- rowMeans(scale(E[ids, emt_genes, drop=FALSE]))
msg("EMT score genes:", paste(emt_genes, collapse=", "))

res <- rbindlist(lapply(MIR, function(m) rbindlist(lapply(TG, function(g){
  x <- Mm[m, use]; y <- E[ids, g]
  s  <- sp(x,y)
  pc <- tryCatch(pcor(x,y,cbind(emt)), error=function(e) list(rho=NA,p=NA,n=NA))
  det <- (2^y - 1) >= 1
  s2 <- if(sum(det)>=8) sp(x[det], y[det]) else list(rho=NA_real_,p=NA_real_,n=sum(det))
  data.table(miRNA=m, target=g, n=s$n, rho=s$rho, p=s$p,
             rho_partial_EMT=pc$rho, p_partial_EMT=pc$p,
             n_detectable=sum(det), rho_detectable_only=s2$rho, p_detectable_only=s2$p)
}))))
res[, q := p.adjust(p, method="BH")]
fwrite(res, file.path(OUT,"celllines_miR29_collagen_correlations.csv"))
print(res)

## -------- all network miRNA_target edges testable in cell lines -------------
ed <- fread(file.path(REV,"data/canonical_edges.tsv"))
mt <- ed[edge_type=="miRNA_target" & source %in% rownames(Mm) & target %in% colnames(E)]
msg("canonical miRNA_target edges testable in CCLE lines:", nrow(mt),
    "| distinct miRNAs:", uniqueN(mt$source), "| distinct targets:", uniqueN(mt$target))
tier <- fread(file.path(REV,"data/edge_evidence_tier.tsv"))
mt <- merge(mt, unique(tier[,.(source,target,tier)]), by=c("source","target"), all.x=TRUE)
ee <- rbindlist(lapply(seq_len(nrow(mt)), function(i){
  s <- sp(Mm[mt$source[i], use], E[ids, mt$target[i]])
  data.table(miRNA=mt$source[i], target=mt$target[i], tier=mt$tier[i],
             n=s$n, rho=s$rho, p=s$p)}))
ee <- ee[is.finite(rho)]
ee[, q := p.adjust(p, method="BH")]
fwrite(ee, file.path(OUT,"celllines_mirna_target_edge_rho.csv"))
summ <- ee[, .(n_edges=.N, frac_negative=mean(rho<0), mean_rho=mean(rho),
               median_rho=median(rho),
               binom_p=binom.test(sum(rho<0), .N, 0.5)$p.value), by=tier]
summ <- rbind(ee[, .(tier="ALL", n_edges=.N, frac_negative=mean(rho<0), mean_rho=mean(rho),
                     median_rho=median(rho),
                     binom_p=binom.test(sum(rho<0), .N, 0.5)$p.value)], summ)
fwrite(summ, file.path(OUT,"celllines_mirna_target_concordance.csv"))
print(summ)

## -------- TF_target Activation concordance in cell lines (N8 analogue) -----
tt <- ed[edge_type=="TF_target" & source %in% colnames(E) & target %in% colnames(E) &
         source!=target]
trr <- fread(file.path(REV,"data/db/trrust_human.tsv"), header=FALSE,
             col.names=c("TF","target","mode","pmid"))
md <- trr[, .(modes=paste(sort(unique(mode)),collapse=";")), by=.(TF,target)]
md[, cls := fifelse(grepl("Activation",modes) & grepl("Repression",modes), "Ambiguous",
             fifelse(grepl("Activation",modes), "Activation",
             fifelse(grepl("Repression",modes), "Repression", "Unknown")))]
tt <- merge(tt, md[,.(source=TF,target,cls)], by=c("source","target"), all.x=TRUE)
msg("TF_target edges testable in cell lines:", nrow(tt), "| Activation:", sum(tt$cls=="Activation",na.rm=TRUE))
te <- rbindlist(lapply(seq_len(nrow(tt)), function(i){
  s <- sp(E[W$breast_ca, tt$source[i]], E[W$breast_ca, tt$target[i]])
  data.table(TF=tt$source[i], target=tt$target[i], cls=tt$cls[i], n=s$n, rho=s$rho, p=s$p)}))
te <- te[is.finite(rho)]
fwrite(te, file.path(OUT,"celllines_TF_target_edge_rho.csv"))
tsum <- te[, .(n_edges=.N, frac_expected_sign=mean(if(.BY$cls=="Repression") rho<0 else rho>0),
               mean_rho=mean(rho)), by=cls]
print(tsum)
fwrite(tsum, file.path(OUT,"celllines_TF_target_concordance.csv"))
msg("done")
