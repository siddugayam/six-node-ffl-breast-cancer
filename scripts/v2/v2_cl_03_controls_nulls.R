#!/usr/bin/env Rscript
## ============================================================================
## STROMA-FREE TEST -- power, positive controls, matched nulls, pan-cancer
## ============================================================================
suppressPackageStartupMessages({library(data.table)})
options(warn=1); set.seed(20260908)
REV <- "/path/to/revision"
DM  <- file.path(REV,"data/depmap24q4"); OUT <- file.path(REV,"results/v2")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
sp <- function(x,y){ ok<-is.finite(x)&is.finite(y); n<-sum(ok)
  if(n<8) return(c(rho=NA,p=NA,n=n))
  ct<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  c(rho=unname(ct$estimate),p=ct$p.value,n=n) }
fz <- function(r,n){ r<-unname(r); n<-unname(n); z<-atanh(r); se<-1/sqrt(n-3); c(lo=tanh(z-1.96*se), hi=tanh(z+1.96*se)) }
zdiff <- function(r1,n1,r2,n2){ r1<-unname(r1);r2<-unname(r2);n1<-unname(n1);n2<-unname(n2); z<-(atanh(r1)-atanh(r2))/sqrt(1/(n1-3)+1/(n2-3))
  c(z=z, p=2*pnorm(-abs(z))) }

W  <- readRDS(file.path(DM,"cl_workspace.rds")); E<-W$E
mod<- fread(file.path(DM,"Model.csv"))
mi <- readRDS(file.path(DM,"ccle_mirna.rds")); Mm<-mi$M; rownames(Mm)<-mi$desc
Mm <- log2(Mm+1)
key<- mod[CCLEName!="", .(CCLEName,ModelID,OncotreeLineage)]; setkey(key,CCLEName)
cn <- colnames(Mm); mid <- key[cn,ModelID]; lin <- key[cn,OncotreeLineage]
ok_all <- which(!is.na(mid) & mid %in% rownames(E))
msg("ALL CCLE lines with miRNA + 24Q4 expression:", length(ok_all))
ids_all <- mid[ok_all]
br <- which(!is.na(mid) & mid %in% W$breast_ca & lin=="Breast"); ids_br <- mid[br]
mes_lin <- c("Soft Tissue","Bone","Fibroblast","Pleura")
me <- which(!is.na(mid) & mid %in% rownames(E) & lin %in% mes_lin); ids_me <- mid[me]
msg("mesenchymal-lineage lines (Soft Tissue/Bone/Fibroblast/Pleura) with both:", length(me))

## ---------------- 1. positive controls for the miRNA platform --------------
ctrl <- data.table(
  miRNA = c("hsa-miR-200c","hsa-miR-200b","hsa-miR-200a","hsa-miR-141",
            "hsa-miR-200c","hsa-miR-200b","hsa-let-7a","hsa-let-7b",
            "hsa-miR-21","hsa-miR-17","hsa-miR-155"),
  target= c("ZEB1","ZEB1","ZEB1","ZEB1","ZEB2","ZEB2","HMGA2","HMGA2",
            "PDCD4","CDKN1A","CEBPB"))
ctrl <- ctrl[miRNA %in% rownames(Mm) & target %in% colnames(E)]
cc <- rbindlist(lapply(seq_len(nrow(ctrl)), function(i){
  a <- sp(Mm[ctrl$miRNA[i], br],    E[ids_br,  ctrl$target[i]])
  b <- sp(Mm[ctrl$miRNA[i], ok_all],E[ids_all, ctrl$target[i]])
  data.table(miRNA=ctrl$miRNA[i], target=ctrl$target[i],
             rho_breast50=a["rho"], p_breast50=a["p"], n_breast=a["n"],
             rho_allCCLE=b["rho"], p_allCCLE=b["p"], n_all=b["n"])}))
fwrite(cc, file.path(OUT,"celllines_mirna_positive_controls.csv")); print(cc)

## ---------------- 2. miR-29 -> collagen: breast / pan-cancer / mesenchymal --
MIR <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
TG  <- c("COL1A1","COL3A1","COL1A2","COL5A1","COL4A1","FN1","SPARC","LOX","MMP2")
TG  <- intersect(TG, colnames(E))
grid <- CJ(m=MIR, g=TG, sorted=FALSE)
mk <- function(idx, ids, lab) rbindlist(lapply(seq_len(nrow(grid)), function(i){
  s <- sp(Mm[grid$m[i], idx], E[ids, grid$g[i]])
  ci <- if(is.finite(s["rho"])) fz(s["rho"], s["n"]) else c(lo=NA,hi=NA)
  data.table(set=lab, miRNA=grid$m[i], target=grid$g[i], n=s["n"], rho=s["rho"],
             ci_lo=ci["lo"], ci_hi=ci["hi"], p=s["p"])}))
m29 <- rbind(mk(br,ids_br,"breast_cancer_lines"),
             mk(ok_all,ids_all,"all_CCLE_lines"),
             mk(me,ids_me,"mesenchymal_lineage_lines"))
m29[, q := p.adjust(p,"BH"), by=set]
fwrite(m29, file.path(OUT,"celllines_miR29_multiset.csv"))
print(m29[target %in% c("COL1A1","COL3A1","COL1A2")][order(target,miRNA,set)])

## ---------------- 3. TCGA vs cell line: formal difference ------------------
X  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
Xm <- readRDS(file.path(REV,"data/brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(REV,"data/brca_pheno.rds"))
tum<- intersect(ph$sample[ph$sample_type=="Primary Tumor"], colnames(X))
pair <- intersect(tum, colnames(Xm))
msg("TCGA paired tumours (gene+miRNA):", length(pair))
mmap <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c")
stopifnot(all(mmap %in% rownames(Xm)))
cmp <- rbindlist(lapply(c("ETS1","NFKB1","RELA","SP1"), function(tf)
        rbindlist(lapply(c("COL1A1","COL3A1"), function(g){
  a <- sp(X[tf,tum], X[g,tum]); b <- sp(E[W$breast_ca,tf], E[W$breast_ca,g])
  z <- zdiff(a["rho"],a["n"],b["rho"],b["n"])
  data.table(axis=paste0(tf,"->",g), arm="TF",
             tcga_rho=a["rho"], tcga_n=a["n"], cl_rho=b["rho"], cl_n=b["n"],
             cl_ci_lo=fz(b["rho"],b["n"])["lo"], cl_ci_hi=fz(b["rho"],b["n"])["hi"],
             z_diff=z["z"], p_diff=z["p"])}))))
cmp <- rbind(cmp, rbindlist(lapply(mmap, function(m)
        rbindlist(lapply(c("COL1A1","COL3A1"), function(g){
  a <- sp(Xm[m,pair], X[g,pair]); b <- sp(Mm[m,br], E[ids_br,g])
  z <- zdiff(a["rho"],a["n"],b["rho"],b["n"])
  data.table(axis=paste0(m,"->",g), arm="miRNA",
             tcga_rho=a["rho"], tcga_n=a["n"], cl_rho=b["rho"], cl_n=b["n"],
             cl_ci_lo=fz(b["rho"],b["n"])["lo"], cl_ci_hi=fz(b["rho"],b["n"])["hi"],
             z_diff=z["z"], p_diff=z["p"])})))))
cmp <- rbind(cmp, {a<-sp(X["COL1A1",tum],X["COL3A1",tum]); b<-sp(E[W$breast_ca,"COL1A1"],E[W$breast_ca,"COL3A1"])
  z<-zdiff(a["rho"],a["n"],b["rho"],b["n"])
  data.table(axis="COL1A1<->COL3A1", arm="gene_gene", tcga_rho=a["rho"],tcga_n=a["n"],
             cl_rho=b["rho"],cl_n=b["n"],cl_ci_lo=fz(b["rho"],b["n"])["lo"],
             cl_ci_hi=fz(b["rho"],b["n"])["hi"],z_diff=z["z"],p_diff=z["p"])})
cmp[, q_diff := p.adjust(p_diff,"BH")]
fwrite(cmp, file.path(OUT,"celllines_tcga_vs_cellline_named_axes.csv")); print(cmp)

## ---------------- 4. power ------------------------------------------------
pw <- data.table(n=c(50,69), rho_detectable_80pct=NA_real_)
for(i in 1:2){ n<-pw$n[i]
  f<-function(r) pnorm(abs(atanh(r))*sqrt(n-3)-1.959964)-0.80
  pw$rho_detectable_80pct[i] <- uniroot(f,c(0.01,0.95))$root }
print(pw); fwrite(pw, file.path(OUT,"celllines_power.csv"))

## ---------------- 5. matched nulls for the two concordance panels ---------
NDRAW <- 20L
decile <- function(v) cut(v, quantile(v, probs=seq(0,1,.1), na.rm=TRUE),
                          include.lowest=TRUE, labels=FALSE)
## (a) miRNA_target edges, cell lines
ee <- fread(file.path(OUT,"celllines_mirna_target_edge_rho.csv"))
mir_pool <- unique(ee$miRNA); gen_pool <- unique(ee$target)
mirmu <- rowMeans(Mm[mir_pool, br, drop=FALSE]); genmu <- colMeans(E[ids_br, gen_pool, drop=FALSE])
md <- decile(mirmu); gd <- decile(genmu)
real <- paste(ee$miRNA, ee$target)
nullr <- numeric(0); it <- 0L
for(k in 1:NDRAW) for(i in seq_len(nrow(ee))){
  it <- it+1L
  mi_c <- sample(mir_pool[md==md[match(ee$miRNA[i],mir_pool)]],1)
  gi_c <- sample(gen_pool[gd==gd[match(ee$target[i],gen_pool)]],1)
  if(paste(mi_c,gi_c) %in% real) next
  nullr <- c(nullr, suppressWarnings(cor(Mm[mi_c,br], E[ids_br,gi_c], method="spearman")))
}
nullr <- nullr[is.finite(nullr)]
msg("miRNA null draws attempted:", it, "| retained:", length(nullr))
## (b) TF_target edges, cell lines
te <- fread(file.path(OUT,"celllines_TF_target_edge_rho.csv"))
tf_pool <- unique(te$TF); tg_pool <- unique(te$target)
tfmu <- colMeans(E[W$breast_ca, tf_pool, drop=FALSE]); tgmu <- colMeans(E[W$breast_ca, tg_pool, drop=FALSE])
td <- decile(tfmu); gd2 <- decile(tgmu); realt <- paste(te$TF, te$target)
nullt <- numeric(0); it2 <- 0L
for(k in 1:NDRAW) for(i in seq_len(nrow(te))){
  it2 <- it2+1L
  a <- sample(tf_pool[td==td[match(te$TF[i],tf_pool)]],1)
  b <- sample(tg_pool[gd2==gd2[match(te$target[i],tg_pool)]],1)
  if(a==b || paste(a,b) %in% realt) next
  nullt <- c(nullt, suppressWarnings(cor(E[W$breast_ca,a], E[W$breast_ca,b], method="spearman")))
}
nullt <- nullt[is.finite(nullt)]
msg("TF null draws attempted:", it2, "| retained:", length(nullt))

nsum <- rbindlist(list(
 data.table(panel="miRNA_target_ALL", obs_frac=mean(ee$rho<0), obs_mean_rho=mean(ee$rho),
   n_edges=nrow(ee), null_frac=mean(nullr<0), null_mean_rho=mean(nullr), n_null=length(nullr),
   p=prop.test(c(sum(ee$rho<0), sum(nullr<0)), c(nrow(ee), length(nullr)))$p.value),
 data.table(panel="miRNA_target_strong", obs_frac=mean(ee[tier=="strong"]$rho<0),
   obs_mean_rho=mean(ee[tier=="strong"]$rho), n_edges=nrow(ee[tier=="strong"]),
   null_frac=mean(nullr<0), null_mean_rho=mean(nullr), n_null=length(nullr),
   p=prop.test(c(sum(ee[tier=="strong"]$rho<0), sum(nullr<0)), c(nrow(ee[tier=="strong"]), length(nullr)))$p.value),
 data.table(panel="miRNA_target_predicted_only", obs_frac=mean(ee[tier=="predicted_only"]$rho<0),
   obs_mean_rho=mean(ee[tier=="predicted_only"]$rho), n_edges=nrow(ee[tier=="predicted_only"]),
   null_frac=mean(nullr<0), null_mean_rho=mean(nullr), n_null=length(nullr),
   p=prop.test(c(sum(ee[tier=="predicted_only"]$rho<0), sum(nullr<0)), c(nrow(ee[tier=="predicted_only"]), length(nullr)))$p.value),
 data.table(panel="TF_target_Activation", obs_frac=mean(te[cls=="Activation"]$rho>0),
   obs_mean_rho=mean(te[cls=="Activation"]$rho), n_edges=nrow(te[cls=="Activation"]),
   null_frac=mean(nullt>0), null_mean_rho=mean(nullt), n_null=length(nullt),
   p=prop.test(c(sum(te[cls=="Activation"]$rho>0), sum(nullt>0)), c(nrow(te[cls=="Activation"]), length(nullt)))$p.value),
 data.table(panel="TF_target_Repression", obs_frac=mean(te[cls=="Repression"]$rho<0),
   obs_mean_rho=mean(te[cls=="Repression"]$rho), n_edges=nrow(te[cls=="Repression"]),
   null_frac=mean(nullt<0), null_mean_rho=mean(nullt), n_null=length(nullt),
   p=prop.test(c(sum(te[cls=="Repression"]$rho<0), sum(nullt<0)), c(nrow(te[cls=="Repression"]), length(nullt)))$p.value)))
fwrite(nsum, file.path(OUT,"celllines_concordance_vs_matched_null.csv")); print(nsum)
msg("done")
