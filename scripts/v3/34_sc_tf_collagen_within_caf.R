#!/usr/bin/env Rscript
## PART 2J) Does the TF -> COL1A1 association exist WITHIN the CAF compartment alone?
## Donor-level pseudobulk (Human Breast Cancer Single Cell Atlas, 138 donors, 8 studies),
## compared with the same correlation in bulk TCGA-BRCA.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

D <- fread(file.path(OUT,"sc_hbca_donor_compartment_pseudobulk.csv"))
msg("donor x compartment rows:", nrow(D), "| compartments:", paste(sort(unique(D$compartment)), collapse=","))
TFS <- c("NFKB1","RELA","SP1","ETS1","MYC","STAT3","HIF1A","TP53","MKL1","SRF","RUNX1","RUNX2",
         "TWIST1","SNAI2","ZEB1","JUN","FOS","EGR1","CEBPB","YAP1","WWTR1","SMAD3","SMAD4","TGFB1","TGFBR2")
TFS <- intersect(TFS, names(D))
TARG <- c("COL1A1","COL3A1","COL1A2","POSTN","FN1","SERPINE1","LOX")
TARG <- intersect(TARG, names(D))

logcpm <- function(x) log2(x + 1)
run <- function(dt, label, minc=50) {
  dt <- dt[n_cells >= minc]
  msg(label, ": donors with >=", minc, "cells:", nrow(dt))
  rbindlist(lapply(TFS, function(tf) rbindlist(lapply(TARG, function(tg) {
    x <- logcpm(dt[[tf]]); y <- logcpm(dt[[tg]])
    ok <- is.finite(x) & is.finite(y); if (sum(ok) < 20) return(NULL)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
    ## within-study partial: residualise both on study
    st <- factor(dt$study[ok])
    rx <- residuals(lm(rank(x[ok]) ~ st)); ry <- residuals(lm(rank(y[ok]) ~ st))
    cp <- suppressWarnings(cor.test(rx, ry, method="pearson"))
    data.table(stratum=label, TF=tf, target=tg, n_donors=sum(ok),
               rho=unname(ct$estimate), p=ct$p.value,
               rho_within_study=unname(cp$estimate), p_within_study=cp$p.value)
  }))))
}
CAF <- run(D[compartment=="CAF"], "within_CAF")
MAL <- run(D[compartment=="malignant"], "within_malignant")
## cross-compartment: malignant-cell TF vs CAF collagen, same donor
mm <- merge(D[compartment=="malignant", c("donor","study","n_cells", TFS), with=FALSE],
            D[compartment=="CAF", c("donor","n_cells", TARG), with=FALSE], by="donor", suffixes=c("_mal","_caf"))
mm <- mm[n_cells_mal>=50 & n_cells_caf>=50]
msg("donors with >=50 malignant AND >=50 CAF cells:", nrow(mm))
CROSS <- rbindlist(lapply(TFS, function(tf) rbindlist(lapply(TARG, function(tg) {
  x <- logcpm(mm[[tf]]); y <- logcpm(mm[[tg]]); ok <- is.finite(x)&is.finite(y)
  if (sum(ok)<20) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method="spearman", exact=FALSE))
  st <- factor(mm$study[ok]); rx <- residuals(lm(rank(x[ok])~st)); ry <- residuals(lm(rank(y[ok])~st))
  cp <- suppressWarnings(cor.test(rx, ry, method="pearson"))
  data.table(stratum="malignantTF_vs_CAFcollagen", TF=tf, target=tg, n_donors=sum(ok),
             rho=unname(ct$estimate), p=ct$p.value,
             rho_within_study=unname(cp$estimate), p_within_study=cp$p.value)
}))))

## bulk TCGA reference for the same TF-target pairs
G <- readRDS(file.path(REV,"data/brca_gene_expr.rds")); ph <- readRDS(file.path(REV,"data/brca_pheno.rds"))
stopifnot(!any(duplicated(rownames(G))))
tum <- intersect(colnames(G), ph$sample[ph$sample_type=="Primary Tumor"])
msg("TCGA primary tumours:", length(tum))
BULK <- rbindlist(lapply(TFS, function(tf) rbindlist(lapply(TARG, function(tg) {
  if (!(tf %in% rownames(G)) || !(tg %in% rownames(G))) return(NULL)
  ct <- suppressWarnings(cor.test(as.numeric(G[tf,tum]), as.numeric(G[tg,tum]), method="spearman", exact=FALSE))
  data.table(stratum="TCGA_bulk", TF=tf, target=tg, n_donors=length(tum),
             rho=unname(ct$estimate), p=ct$p.value, rho_within_study=NA_real_, p_within_study=NA_real_)
}))))

R <- rbindlist(list(BULK, CAF, MAL, CROSS))
R[, fdr := p.adjust(p, "BH"), by=stratum]
fwrite(R, file.path(OUT,"sc_tf_collagen_within_compartment.csv"))
msg("wrote sc_tf_collagen_within_compartment.csv:", nrow(R), "rows")

msg("\n=== TF vs COL1A1: bulk TCGA versus within-CAF pseudobulk ===")
W <- dcast(R[target=="COL1A1"], TF ~ stratum, value.var="rho")
Wp <- dcast(R[target=="COL1A1"], TF ~ stratum, value.var="p")
setnames(Wp, setdiff(names(Wp),"TF"), paste0("p_", setdiff(names(Wp),"TF")))
print(merge(W, Wp, by="TF"), digits=3, nrows=40)
msg("\n=== same for COL3A1 ===")
print(dcast(R[target=="COL3A1"], TF ~ stratum, value.var="rho"), digits=3, nrows=40)
msg("\n=== within-CAF, study-adjusted, COL1A1 (sorted) ===")
print(R[stratum=="within_CAF" & target=="COL1A1"][order(-abs(rho_within_study)),
        .(TF, n_donors, rho, p, fdr, rho_within_study, p_within_study)], digits=3, nrows=40)
msg("\nnumber of TFs positively and significantly (FDR<0.05) associated with COL1A1:")
for (s in unique(R$stratum)) msg("  ", s, ":", R[stratum==s & target=="COL1A1" & rho>0 & fdr<0.05, .N], "of", R[stratum==s & target=="COL1A1", .N])
msg("DONE 34")
