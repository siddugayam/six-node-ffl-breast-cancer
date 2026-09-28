#!/usr/bin/env Rscript
# 24_nfkb_activity_score.R -- NF-kB TRANSCRIPTIONAL ACTIVITY (target-gene signature) vs collagen,
# because RELA has no phospho-site in the CPTAC-BRCA phosphoproteome.
REV <- "/path/to/revision"; OUT <- file.path(REV,"results","multiomics")
con <- file(file.path(REV,"logs","nfkb_activity.log"), open="wt")
say <- function(...) { m<-paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse="")); cat(m,"\n"); cat(m,"\n",file=con); flush(con) }
say("=== 24 START ===")
b <- readRDS(file.path(OUT,"cptac_bundle.rds"))
edges <- read.delim(file.path(REV,"data","canonical_edges.tsv"), stringsAsFactors=FALSE)
S <- b$s_rp; Rm <- b$Rg[,S,drop=FALSE]; Pm <- b$Pg[,S,drop=FALSE]
tg <- unique(edges$target[edges$edge_type=="TF_target" & edges$source %in% c("NFKB1","RELA")])
tg <- setdiff(intersect(tg, rownames(Rm)), c("COL1A1","COL3A1","NFKB1","RELA"))
say("NF-kB (NFKB1/RELA) target genes in the canonical network measured in CPTAC RNA: ", length(tg))
say("  ", paste(sort(tg), collapse=","))
zall <- t(scale(t(Rm[tg,,drop=FALSE])))
score <- colMeans(zall, na.rm=TRUE)
score_loo <- function(g) { keep <- setdiff(tg, g); colMeans(zall[keep,,drop=FALSE], na.rm=TRUE) }
ct <- function(x,y){ k<-!is.na(x)&!is.na(y); if(sum(k)<30) return(c(NA,NA,sum(k)))
  h<-suppressWarnings(cor.test(x[k],y[k],method="spearman",exact=FALSE)); c(unname(h$estimate),h$p.value,sum(k)) }
rows <- list()
for (g in c("COL1A1","COL3A1","FN1","POSTN","VIM","VEGFA")) {
  sc <- if (g %in% tg) score_loo(g) else score   # leave-one-out: never include the tested gene
  if (g %in% rownames(Pm)) { v<-ct(sc, Pm[g,]); rows[[length(rows)+1]] <- data.frame(score=ifelse(g %in% tg,"NFKB_target_signature_mRNA_LOO","NFKB_target_signature_mRNA"), vs=paste0(g,"_PROTEIN"), rho=v[1],p=v[2],n=v[3]) }
  if (g %in% rownames(Rm)) { v<-ct(sc, Rm[g,]); rows[[length(rows)+1]] <- data.frame(score=ifelse(g %in% tg,"NFKB_target_signature_mRNA_LOO","NFKB_target_signature_mRNA"), vs=paste0(g,"_mRNA"),    rho=v[1],p=v[2],n=v[3]) }
}
for (g in c("NFKB1","RELA")) {
  v<-ct(score, Rm[g,]); rows[[length(rows)+1]] <- data.frame(score="NFKB_target_signature_mRNA", vs=paste0(g,"_mRNA"), rho=v[1],p=v[2],n=v[3])
  if (g %in% rownames(Pm)) { v<-ct(score, Pm[g,]); rows[[length(rows)+1]] <- data.frame(score="NFKB_target_signature_mRNA", vs=paste0(g,"_PROTEIN"), rho=v[1],p=v[2],n=v[3]) }
}
D <- do.call(rbind, rows); D$fdr <- p.adjust(D$p,"BH")
D$n_signature_genes <- length(tg)
## random-signature null: is the collagen association specific to NF-kB targets?
set.seed(1234)
allg <- rownames(Rm)[apply(Rm,1,function(v) sd(v,na.rm=TRUE))>0]
NB <- 1000; nullrho <- matrix(NA_real_, NB, 2, dimnames=list(NULL,c("COL1A1_PROTEIN","COL3A1_PROTEIN")))
for (i in seq_len(NB)) {
  rg <- sample(allg, length(tg)); zz <- t(scale(t(Rm[rg,,drop=FALSE]))); sc <- colMeans(zz, na.rm=TRUE)
  nullrho[i,1] <- ct(sc, Pm["COL1A1",])[1]; nullrho[i,2] <- ct(sc, Pm["COL3A1",])[1]
}
obs1 <- D$rho[D$vs=="COL1A1_PROTEIN"][1]; obs2 <- D$rho[D$vs=="COL3A1_PROTEIN"][1]
say("RANDOM-SIGNATURE NULL (", NB, " random ", length(tg), "-gene signatures):")
say("  COL1A1 protein: observed rho=", signif(obs1,4), " ; null mean=", signif(mean(nullrho[,1]),4),
    " sd=", signif(sd(nullrho[,1]),4), " ; empirical p(one-sided)=", signif((sum(nullrho[,1]>=obs1)+1)/(NB+1),4))
say("  COL3A1 protein: observed rho=", signif(obs2,4), " ; null mean=", signif(mean(nullrho[,2]),4),
    " sd=", signif(sd(nullrho[,2]),4), " ; empirical p(one-sided)=", signif((sum(nullrho[,2]>=obs2)+1)/(NB+1),4))
D$random_signature_null_mean <- NA_real_; D$random_signature_p <- NA_real_
D$random_signature_null_mean[D$vs=="COL1A1_PROTEIN"] <- mean(nullrho[,1])
D$random_signature_null_mean[D$vs=="COL3A1_PROTEIN"] <- mean(nullrho[,2])
D$random_signature_p[D$vs=="COL1A1_PROTEIN"] <- (sum(nullrho[,1]>=obs1)+1)/(NB+1)
D$random_signature_p[D$vs=="COL3A1_PROTEIN"] <- (sum(nullrho[,2]>=obs2)+1)/(NB+1)
write.csv(D, file.path(OUT,"cptac_nfkb_activity_score.csv"), row.names=FALSE)
say("WROTE cptac_nfkb_activity_score.csv rows=", nrow(D))
for (i in seq_len(nrow(D))) say(sprintf("  NFkB-target-signature vs %-16s rho=%+.4f p=%.3g fdr=%.3g n=%d", D$vs[i],D$rho[i],D$p[i],D$fdr[i],D$n[i]))
say("=== 24 DONE ==="); close(con)
