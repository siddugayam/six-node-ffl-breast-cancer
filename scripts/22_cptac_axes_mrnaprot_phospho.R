#!/usr/bin/env Rscript
# 22_cptac_axes_mrnaprot_phospho.R
# C) named axes at protein level   D) per-gene mRNA-protein concordance
# E) NF-kB phospho-site activity vs collagen protein
# plus sensitivity analyses for A/B (TF mRNA -> target protein ; complete-protein subset)
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results","multiomics")
con <- file(file.path(REV,"logs","cptac_axes_phospho.log"), open="wt")
say <- function(...) { m<-paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse="")); cat(m,"\n"); cat(m,"\n",file=con); flush(con) }
set.seed(1234); say("=== 22 START ===")
b  <- readRDS(file.path(OUT,"cptac_bundle.rds"))
W  <- readRDS(file.path(OUT,"cptac_AB_workspace.rds"))
S3 <- b$s_all; S_RP <- b$s_rp; S_PPH <- b$s_pph
MIN_N <- 30L

ct <- function(x,y,label_a,label_b,label_set) {
  k <- !is.na(x) & !is.na(y)
  if (sum(k)<MIN_N || sd(x[k])==0 || sd(y[k])==0)
    return(data.frame(feature_a=label_a,feature_b=label_b,sample_set=label_set,rho=NA_real_,p=NA_real_,n=sum(k)))
  h <- suppressWarnings(cor.test(x[k],y[k],method="spearman",exact=FALSE))
  data.frame(feature_a=label_a,feature_b=label_b,sample_set=label_set,
             rho=unname(h$estimate),p=h$p.value,n=sum(k))
}

## ---------------- C) NAMED AXES AT PROTEIN LEVEL ----------------
say("=== C) named axes, protein level ===")
Pm3 <- b$Pg[,S3,drop=FALSE]; Rm3 <- b$Rg[,S3,drop=FALSE]; Mm3 <- b$MIc[,S3,drop=FALSE]
Pmr <- b$Pg[,S_RP,drop=FALSE]; Rmr <- b$Rg[,S_RP,drop=FALSE]
getP <- function(g,X) if (g %in% rownames(X)) X[g,] else rep(NA_real_,ncol(X))
axes <- list()
mir_gene <- list(c("hsa-miR-29a","COL1A1"),c("hsa-miR-29b","COL1A1"),c("hsa-miR-29c","COL1A1"),
                 c("hsa-miR-29a","COL3A1"),c("hsa-miR-29b","COL3A1"),c("hsa-miR-29c","COL3A1"),
                 c("hsa-let-7b","COL3A1"),c("hsa-let-7e","COL3A1"),c("hsa-miR-101","EZH2"),
                 c("hsa-miR-130a","VEGFA"),c("hsa-miR-124","STAT3"))
for (pr in mir_gene) {
  axes[[length(axes)+1]] <- cbind(ct(getP(pr[1],Mm3), getP(pr[2],Pm3), pr[1], paste0(pr[2],"_PROTEIN"),"RNA+prot+miR n<=101"), layer="miRNA_vs_protein")
  axes[[length(axes)+1]] <- cbind(ct(getP(pr[1],Mm3), getP(pr[2],Rm3), pr[1], paste0(pr[2],"_mRNA"),   "RNA+prot+miR n<=101"), layer="miRNA_vs_mRNA")
}
tf_gene <- list(c("NFKB1","COL1A1"),c("RELA","COL1A1"),c("SP1","COL1A1"),c("ETS1","COL1A1"),
                c("NFKB1","COL3A1"),c("RELA","COL3A1"),c("SP1","COL3A1"),c("ETS1","COL3A1"))
for (pr in tf_gene) {
  axes[[length(axes)+1]] <- cbind(ct(getP(pr[1],Pmr), getP(pr[2],Pmr), paste0(pr[1],"_PROTEIN"), paste0(pr[2],"_PROTEIN"),"RNA+prot n<=121"), layer="TFprotein_vs_protein")
  axes[[length(axes)+1]] <- cbind(ct(getP(pr[1],Rmr), getP(pr[2],Pmr), paste0(pr[1],"_mRNA"),    paste0(pr[2],"_PROTEIN"),"RNA+prot n<=121"), layer="TFmRNA_vs_protein")
  axes[[length(axes)+1]] <- cbind(ct(getP(pr[1],Rmr), getP(pr[2],Rmr), paste0(pr[1],"_mRNA"),    paste0(pr[2],"_mRNA"),   "RNA+prot n<=121"), layer="TFmRNA_vs_mRNA")
}
axes[[length(axes)+1]] <- cbind(ct(getP("COL1A1",Pmr), getP("COL3A1",Pmr),"COL1A1_PROTEIN","COL3A1_PROTEIN","RNA+prot n<=121"), layer="protein_protein")
axes[[length(axes)+1]] <- cbind(ct(getP("COL1A1",Rmr), getP("COL3A1",Rmr),"COL1A1_mRNA","COL3A1_mRNA","RNA+prot n<=121"), layer="mRNA_mRNA")
AX <- do.call(rbind, axes)
AX$fdr_within_table <- p.adjust(AX$p, "BH")
write.csv(AX, file.path(OUT,"cptac_named_axes_protein.csv"), row.names=FALSE)
say("WROTE cptac_named_axes_protein.csv rows=", nrow(AX))
for (i in seq_len(nrow(AX))) say(sprintf("%-22s %-18s %-20s rho=%s p=%s n=%s", AX$feature_a[i],AX$feature_b[i],AX$layer[i],
   ifelse(is.na(AX$rho[i]),"NA",sprintf("%+.4f",AX$rho[i])), ifelse(is.na(AX$p[i]),"NA",sprintf("%.3g",AX$p[i])), AX$n[i]))

## ---------------- D) mRNA-PROTEIN CONCORDANCE PER GENE ----------------
say("=== D) per-gene mRNA-protein Spearman (n<=121 RNA+protein samples) ===")
gs <- rownames(b$Pg)
D <- do.call(rbind, lapply(gs, function(g) {
  x <- Rmr[g,]; y <- Pmr[g,]; k <- !is.na(x)&!is.na(y)
  if (sum(k)<MIN_N || sd(x[k])==0 || sd(y[k])==0) return(data.frame(gene=g,rho=NA_real_,p=NA_real_,n=sum(k)))
  h <- suppressWarnings(cor.test(x[k],y[k],method="spearman",exact=FALSE))
  data.frame(gene=g,rho=unname(h$estimate),p=h$p.value,n=sum(k))
}))
D$node_type <- read.delim(file.path(REV,"data","canonical_nodes.tsv"),stringsAsFactors=FALSE)$type[
  match(D$gene, read.delim(file.path(REV,"data","canonical_nodes.tsv"),stringsAsFactors=FALSE)$name)]
D$fdr <- NA_real_; ok <- !is.na(D$p); D$fdr[ok] <- p.adjust(D$p[ok],"BH")
D <- D[order(-D$rho),]
write.csv(D, file.path(OUT,"cptac_mRNA_protein_concordance.csv"), row.names=FALSE)
say("WROTE cptac_mRNA_protein_concordance.csv rows=", nrow(D))
dv <- D$rho[!is.na(D$rho)]
say("genes tested ", length(dv), " of ", nrow(D), " ; median rho=", signif(median(dv),4),
    " mean=", signif(mean(dv),4), " IQR=", signif(quantile(dv,.25),4), "..", signif(quantile(dv,.75),4),
    " min=", signif(min(dv),4), " max=", signif(max(dv),4))
say("fraction rho>0: ", signif(mean(dv>0),4), " ; rho>0.3: ", signif(mean(dv>0.3),4),
    " ; rho>0.5: ", signif(mean(dv>0.5),4), " ; FDR<0.05: ", signif(mean(D$fdr[!is.na(D$fdr)]<0.05),4))
for (g in c("COL1A1","COL3A1","VEGFA","EZH2","NFKB1","RELA","SP1","ETS1","MYC","TP53")) {
  r <- D[D$gene==g,]
  if (nrow(r)==0) say(sprintf("%-7s NOT MEASURED at protein", g))
  else say(sprintf("%-7s rho=%s p=%s fdr=%s n=%d", g, ifelse(is.na(r$rho),"NA",sprintf("%+.4f",r$rho)),
        ifelse(is.na(r$p),"NA",sprintf("%.3g",r$p)), ifelse(is.na(r$fdr),"NA",sprintf("%.3g",r$fdr)), r$n))
}

## ---------------- E) PHOSPHOPROTEOME: NF-kB ----------------
say("=== E) NF-kB phospho-sites vs collagen protein ===")
PH <- b$PH
ph_ids <- rownames(PH)
ph_ens <- sub("\\..*$","", sapply(strsplit(ph_ids,"\\|"), `[`, 1))
ph_site<- sapply(strsplit(ph_ids,"\\|"), function(v) if(length(v)>=3) v[3] else NA)
want <- c("RELA","NFKB1","NFKB2","IKBKB","NFKBIA","REL")
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=want, keytype="SYMBOL", columns="ENSEMBL"))
m <- m[!is.na(m$ENSEMBL),]
say("target-gene ENSG: ", paste(paste0(m$SYMBOL,"=",m$ENSEMBL), collapse=" "))
sel <- which(ph_ens %in% m$ENSEMBL)
say("phospho rows matching NF-kB genes: ", length(sel))
PHs <- PH[sel,,drop=FALSE]
gsym <- m$SYMBOL[match(ph_ens[sel], m$ENSEMBL)]
sset <- S_PPH
colP <- b$Pg[,sset,drop=FALSE]
rowsE <- list()
for (i in seq_along(sel)) {
  x <- PHs[i, sset]
  nn <- sum(!is.na(x))
  for (tg in c("COL1A1","COL3A1","VIM","FN1")) {
    if (!(tg %in% rownames(colP))) next
    r <- ct(x, colP[tg,], paste0(gsym[i],"_",ph_site[sel][i]), paste0(tg,"_PROTEIN"), "protein+phospho n<=122")
    r$phospho_id <- ph_ids[sel][i]; r$gene <- gsym[i]; r$site <- ph_site[sel][i]; r$n_site_measured <- nn
    rowsE[[length(rowsE)+1]] <- r
  }
  # phospho vs its own total protein (stoichiometry check)
  if (gsym[i] %in% rownames(colP)) {
    r <- ct(x, colP[gsym[i],], paste0(gsym[i],"_",ph_site[sel][i]), paste0(gsym[i],"_TOTALPROTEIN"), "protein+phospho n<=122")
    r$phospho_id <- ph_ids[sel][i]; r$gene <- gsym[i]; r$site <- ph_site[sel][i]; r$n_site_measured <- nn
    rowsE[[length(rowsE)+1]] <- r
  }
}
PHOS <- do.call(rbind, rowsE)
PHOS$fdr_within_table <- NA_real_; ok <- !is.na(PHOS$p); PHOS$fdr_within_table[ok] <- p.adjust(PHOS$p[ok],"BH")
write.csv(PHOS, file.path(OUT,"cptac_nfkb_phospho_vs_collagen.csv"), row.names=FALSE)
say("WROTE cptac_nfkb_phospho_vs_collagen.csv rows=", nrow(PHOS))
say("--- sites detected per gene: ", paste(names(table(gsym)), table(gsym), sep="=", collapse=" "))
say("--- RELA sites: ", paste(unique(ph_site[sel][gsym=="RELA"]), collapse=","))
pp <- PHOS[!is.na(PHOS$rho) & PHOS$feature_b %in% c("COL1A1_PROTEIN","COL3A1_PROTEIN"),]
pp <- pp[order(pp$p),]
for (i in seq_len(nrow(pp))) say(sprintf("%-18s vs %-16s rho=%+.4f p=%.3g fdr=%.3g n=%d", pp$feature_a[i],pp$feature_b[i],pp$rho[i],pp$p[i],pp$fdr_within_table[i],pp$n[i]))
say("--- rows with too few measurements (rho NA): ", sum(is.na(PHOS$rho)))

## also: RELA/NFKB1 mRNA and total protein vs collagen protein already in AX (part C)
## ---------------- SENSITIVITY ----------------
say("=== SENSITIVITY: TF mRNA -> target PROTEIN, and complete-protein subset ===")
A <- W$A; B <- W$B
rho_pair <- function(X,a,Y,bn) { n<-length(a); rr<-rep(NA_real_,n); nn<-rep(NA_integer_,n)
  for (i in seq_len(n)) { if(!(a[i]%in%rownames(X))||!(bn[i]%in%rownames(Y))) next
    x<-X[a[i],]; y<-Y[bn[i],]; k<-!is.na(x)&!is.na(y); if(sum(k)<MIN_N) next
    if(sd(x[k])==0||sd(y[k])==0) next
    rr[i]<-suppressWarnings(cor(x[k],y[k],method="spearman")); nn[i]<-sum(k) }
  list(rho=rr,n=nn) }
Pmr2 <- b$Pg[,S_RP,drop=FALSE]; Rmr2 <- b$Rg[,S_RP,drop=FALSE]
rX <- rho_pair(Rmr2, B$source, Pmr2, B$target)
B$rho_TFmRNA_targetProtein <- rX$rho; B$n_TFmRNA_targetProtein <- rX$n
sg <- B$trrust_mode %in% c("Activation","Repression")
predsign <- c(Activation=1,Repression=-1)[B$trrust_mode]
concX <- sign(B$rho_TFmRNA_targetProtein)==predsign
SENS <- data.frame(
  test=c("TF mRNA -> target mRNA (n<=121)","TF protein -> target protein (n<=101 set)","TF mRNA -> target protein (n<=121)"),
  n=c(sum(sg&!is.na(B$rho_mRNA)), sum(sg&!is.na(B$rho_protein)), sum(sg&!is.na(B$rho_TFmRNA_targetProtein))),
  concordance=c(mean(B$conc_mRNA[sg],na.rm=TRUE), mean(B$conc_protein[sg],na.rm=TRUE), mean(concX[sg],na.rm=TRUE)),
  mean_rho=c(mean(B$rho_mRNA[sg],na.rm=TRUE), mean(B$rho_protein[sg],na.rm=TRUE), mean(B$rho_TFmRNA_targetProtein[sg],na.rm=TRUE)))
# complete-protein subset for A
complete_prot <- rownames(b$Pg)[rowSums(is.na(b$Pg[,S3,drop=FALSE]))==0]
say("network proteins with ZERO missing values in the 101-sample set: ", length(complete_prot))
iA <- A$target %in% complete_prot & !is.na(A$rho_protein)
SENS2 <- data.frame(
  test=c("miRNA->target mRNA (all usable edges)","miRNA->target protein (all usable edges)",
         "miRNA->target mRNA (complete-protein targets only)","miRNA->target protein (complete-protein targets only)"),
  n=c(sum(!is.na(A$rho_mRNA)), sum(!is.na(A$rho_protein)), sum(iA), sum(iA)),
  concordance=c(mean(A$conc_mRNA,na.rm=TRUE), mean(A$conc_protein,na.rm=TRUE),
                mean(A$conc_mRNA[iA],na.rm=TRUE), mean(A$conc_protein[iA],na.rm=TRUE)),
  mean_rho=c(mean(A$rho_mRNA,na.rm=TRUE), mean(A$rho_protein,na.rm=TRUE),
             mean(A$rho_mRNA[iA],na.rm=TRUE), mean(A$rho_protein[iA],na.rm=TRUE)))
SENSALL <- rbind(cbind(block="B_TF_target_signed",SENS), cbind(block="A_miRNA_target",SENS2))
write.csv(SENSALL, file.path(OUT,"cptac_sensitivity.csv"), row.names=FALSE)
say("WROTE cptac_sensitivity.csv rows=", nrow(SENSALL))
for (i in seq_len(nrow(SENSALL))) say(sprintf("%-22s %-55s n=%4d conc=%.4f meanRho=%+.4f",
  SENSALL$block[i],SENSALL$test[i],SENSALL$n[i],SENSALL$concordance[i],SENSALL$mean_rho[i]))
say("=== 22 DONE ==="); close(con)
