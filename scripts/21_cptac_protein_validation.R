#!/usr/bin/env Rscript
# 21_cptac_protein_validation.R
# A) miRNA->target: mRNA vs PROTEIN correlation, decile-matched null, by evidence tier
# B) TF->target at protein level, split by TRRUST mode
# C) named axes at protein level
# D) per-gene mRNA-protein concordance
# E) NF-kB phospho-site activity vs collagen protein
suppressMessages({library(data.table)})
REV <- "/path/to/revision"
OUT <- file.path(REV,"results","multiomics")
LOG <- file.path(REV,"logs","cptac_protein_validation.log")
con <- file(LOG, open="wt")
say <- function(...) { m<-paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse="")); cat(m,"\n"); cat(m,"\n",file=con); flush(con) }
set.seed(1234)
say("=== 21_cptac_protein_validation.R START ===")

b <- readRDS(file.path(OUT,"cptac_bundle.rds"))
edges <- read.delim(file.path(REV,"data","canonical_edges.tsv"), stringsAsFactors=FALSE)
nodes <- read.delim(file.path(REV,"data","canonical_nodes.tsv"), stringsAsFactors=FALSE)
tier  <- read.delim(file.path(REV,"data","edge_evidence_tier.tsv"), stringsAsFactors=FALSE)
cov_tf<- read.delim(file.path(REV,"results","TF_target_TRRUST_coverage.tsv"), stringsAsFactors=FALSE)

S3   <- b$s_all                      # RNA + protein + miRNA
S_RP <- b$s_rp                       # RNA + protein
MIN_N <- 30L
say("samples: RNA+protein+miRNA n=",length(S3),"  RNA+protein n=",length(S_RP))

Pm <- b$Pg[, S3, drop=FALSE]; Rm <- b$Rg[, S3, drop=FALSE]; Mm <- b$MIc[, S3, drop=FALSE]
nz <- function(X) apply(X,1,function(v) { v<-v[!is.na(v)]; length(v)>=MIN_N && sd(v)>0 })
prot_ok  <- rownames(Pm)[nz(Pm)]
mrna_ok  <- rownames(Rm)[nz(Rm)]
mir_all  <- intersect(nodes$name[nodes$type=="miRNA"], rownames(Mm))
mir_ok   <- mir_all[nz(Mm[mir_all,,drop=FALSE])]
gene_both <- sort(intersect(prot_ok, mrna_ok))
say("network genes usable at protein (>=",MIN_N," non-NA, sd>0): ",length(prot_ok))
say("network genes usable at mRNA: ",length(mrna_ok))
say("network genes usable at BOTH (analysis universe): ",length(gene_both))
say("network miRNAs usable: ",length(mir_ok)," of ",length(mir_all))

## ---- correlation helpers (pairwise complete, Spearman) ----
rho_pair <- function(X, a, Y, bnm) {
  n <- length(a); rr <- rep(NA_real_,n); nn <- rep(NA_integer_,n)
  for (i in seq_len(n)) {
    if (!(a[i] %in% rownames(X)) || !(bnm[i] %in% rownames(Y))) next
    x <- X[a[i],]; y <- Y[bnm[i],]; k <- !is.na(x) & !is.na(y)
    if (sum(k) < MIN_N) next
    if (sd(x[k])==0 || sd(y[k])==0) next
    rr[i] <- suppressWarnings(cor(x[k],y[k],method="spearman")); nn[i] <- sum(k)
  }
  list(rho=rr, n=nn)
}
p_from <- function(rho,n){ tt<-rho*sqrt((n-2)/pmax(1-rho^2,.Machine$double.xmin)); 2*pt(-abs(tt),df=n-2) }
auc_fun <- function(x,y){ n1<-length(x); n2<-length(y); if(n1==0||n2==0) return(NA_real_)
  r<-rank(c(x,y)); (sum(r[seq_len(n1)])-n1*(n1+1)/2)/(n1*n2) }

## ================= EDGE TABLE =================
key <- function(s,t) paste(s,t,sep="\r")
E <- edges
tk <- setNames(tier$tier, key(tier$source,tier$target))
E$evidence_tier <- unname(tk[key(E$source,E$target)]); E$evidence_tier[E$edge_type!="miRNA_target"] <- NA
tfm <- setNames(cov_tf$mode, key(cov_tf$source,cov_tf$target))
E$trrust_mode <- unname(tfm[key(E$source,E$target)]); E$trrust_mode[E$edge_type!="TF_target"] <- NA
E$trrust_mode[E$edge_type=="TF_target" & (is.na(E$trrust_mode)|E$trrust_mode=="")] <- "no_TRRUST_record"

## ---- A) miRNA -> target -----------------------------------------------------
say("=== A) miRNA -> target: mRNA vs PROTEIN ===")
A <- E[E$edge_type=="miRNA_target" & E$source %in% mir_ok & E$target %in% gene_both, ]
say("miRNA_target edges in canonical network: ", sum(E$edge_type=="miRNA_target"))
say("  ...usable here (miRNA measured AND target has BOTH mRNA and protein): ", nrow(A))
say("  tier breakdown: ", paste(names(table(A$evidence_tier)), table(A$evidence_tier), sep="=", collapse=" "))
rA_m <- rho_pair(Mm, A$source, Rm, A$target)
rA_p <- rho_pair(Mm, A$source, Pm, A$target)
A$rho_mRNA <- rA_m$rho; A$n_mRNA <- rA_m$n
A$rho_protein <- rA_p$rho; A$n_protein <- rA_p$n
A$p_mRNA <- p_from(A$rho_mRNA,A$n_mRNA); A$p_protein <- p_from(A$rho_protein,A$n_protein)
A$predicted_sign <- -1L
A$conc_mRNA <- sign(A$rho_mRNA)==A$predicted_sign
A$conc_protein <- sign(A$rho_protein)==A$predicted_sign

## ---- B) TF -> target --------------------------------------------------------
say("=== B) TF -> target: mRNA vs PROTEIN ===")
B <- E[E$edge_type=="TF_target" & E$source %in% gene_both & E$target %in% gene_both, ]
say("TF_target edges in network: ", sum(E$edge_type=="TF_target"), " ; usable (both partners have mRNA+protein): ", nrow(B))
say("  TRRUST mode: ", paste(names(table(B$trrust_mode)), table(B$trrust_mode), sep="=", collapse=" "))
rB_m <- rho_pair(Rm, B$source, Rm, B$target)
rB_p <- rho_pair(Pm, B$source, Pm, B$target)
B$rho_mRNA <- rB_m$rho; B$n_mRNA <- rB_m$n
B$rho_protein <- rB_p$rho; B$n_protein <- rB_p$n
B$p_mRNA <- p_from(B$rho_mRNA,B$n_mRNA); B$p_protein <- p_from(B$rho_protein,B$n_protein)
B$predicted_sign <- c(Activation=1L,Repression=-1L)[B$trrust_mode]
B$conc_mRNA <- ifelse(is.na(B$predicted_sign),NA,sign(B$rho_mRNA)==B$predicted_sign)
B$conc_protein <- ifelse(is.na(B$predicted_sign),NA,sign(B$rho_protein)==B$predicted_sign)

EDGE <- rbind(
  data.frame(A[,c("source","target","edge_type","evidence_tier","trrust_mode","predicted_sign",
                  "rho_mRNA","p_mRNA","n_mRNA","rho_protein","p_protein","n_protein","conc_mRNA","conc_protein")]),
  data.frame(B[,c("source","target","edge_type","evidence_tier","trrust_mode","predicted_sign",
                  "rho_mRNA","p_mRNA","n_mRNA","rho_protein","p_protein","n_protein","conc_mRNA","conc_protein")]))
EDGE$fdr_mRNA <- NA_real_; EDGE$fdr_protein <- NA_real_
for (et in unique(EDGE$edge_type)) {
  i <- EDGE$edge_type==et
  j <- i & !is.na(EDGE$p_mRNA);    EDGE$fdr_mRNA[j]    <- p.adjust(EDGE$p_mRNA[j],"BH")
  j <- i & !is.na(EDGE$p_protein); EDGE$fdr_protein[j] <- p.adjust(EDGE$p_protein[j],"BH")
}
write.csv(EDGE, file.path(OUT,"cptac_miRNA_mRNA_vs_protein.csv"), row.names=FALSE)
say("WROTE cptac_miRNA_mRNA_vs_protein.csv rows=", nrow(EDGE))

## ================= MATCHED RANDOM-PAIR NULL =================
say("=== matched random-pair null (expression-decile + degree matched, >=10 per edge) ===")
deg_out <- table(factor(edges$source, levels=nodes$name)); deg_out <- setNames(as.integer(deg_out),names(deg_out))
deg_in  <- table(factor(edges$target, levels=nodes$name)); deg_in  <- setNames(as.integer(deg_in), names(deg_in))
node_type <- setNames(nodes$type, nodes$name)
Eset <- new.env(hash=TRUE,parent=emptyenv())
for (r in seq_len(nrow(edges))) assign(key(edges$source[r],edges$target[r]),TRUE,envir=Eset)

pool <- c(gene_both, mir_ok)
platform <- setNames(c(rep("gene",length(gene_both)), rep("mirna",length(mir_ok))), pool)
ntype <- node_type[pool]; names(ntype) <- pool
mean_expr_mrna <- c(rowMeans(Rm[gene_both,,drop=FALSE],na.rm=TRUE), rowMeans(Mm[mir_ok,,drop=FALSE],na.rm=TRUE))
mean_expr_prot <- c(rowMeans(Pm[gene_both,,drop=FALSE],na.rm=TRUE), rowMeans(Mm[mir_ok,,drop=FALSE],na.rm=TRUE))
mk_dec <- function(v) { d <- setNames(rep(NA_integer_,length(pool)),pool)
  for (pl in c("gene","mirna")) { nn<-pool[platform==pl]; x<-v[nn]
    br<-unique(quantile(x,probs=seq(0,1,0.1))); d[nn]<-as.integer(cut(x,breaks=br,include.lowest=TRUE)) }
  d }
DEC <- list(mRNA=mk_dec(mean_expr_mrna), protein=mk_dec(mean_expr_prot))
say("null pool: ",length(gene_both)," gene/TF + ",length(mir_ok)," miRNA = ",length(pool)," nodes")

NPER <- 10L
make_null <- function(EDF, dec) {
  cache <- new.env(hash=TRUE,parent=emptyenv())
  get_cand <- function(nd, role) {
    ck <- paste(nd,role,sep="\r")
    if (exists(ck,envir=cache,inherits=FALSE)) return(get(ck,envir=cache))
    ty<-ntype[[nd]]; d0<-dec[[nd]]; dvec <- if(role=="out") deg_out else deg_in
    w<-0L; cands<-character(0)
    repeat { cands <- setdiff(pool[ntype==ty & !is.na(dec) & abs(dec-d0)<=w], nd)
             if (length(cands)>=8L || w>=9L) break; w<-w+1L }
    if (length(cands)>8L) { dd<-abs(log2(dvec[cands]+1)-log2(dvec[[nd]]+1))
      cands<-cands[order(dd)][seq_len(max(8L,ceiling(length(cands)/2)))] }
    res<-list(cands=cands,w=w); assign(ck,res,envir=cache); res
  }
  capS<-vector("list",nrow(EDF)); capT<-vector("list",nrow(EDF)); capI<-vector("list",nrow(EDF)); wS<-integer(nrow(EDF)); wT<-integer(nrow(EDF))
  for (j in seq_len(nrow(EDF))) {
    cs<-get_cand(EDF$source[j],"out"); ct<-get_cand(EDF$target[j],"in"); wS[j]<-cs$w; wT[j]<-ct$w
    got<-0L; tries<-0L; seen<-character(0); ss<-character(0); tt<-character(0)
    while (got<NPER && tries<400L) { tries<-tries+1L
      a<-cs$cands[sample.int(length(cs$cands),1L)]; bb<-ct$cands[sample.int(length(ct$cands),1L)]
      if (a==bb) next; kk<-key(a,bb); if (kk %in% seen) next
      if (exists(kk,envir=Eset,inherits=FALSE)) next
      seen<-c(seen,kk); ss<-c(ss,a); tt<-c(tt,bb); got<-got+1L }
    capS[[j]]<-ss; capT[[j]]<-tt; capI[[j]]<-rep(j,length(ss))
  }
  list(NL=data.frame(edge_row=unlist(capI), rnd_source=unlist(capS), rnd_target=unlist(capT), stringsAsFactors=FALSE),
       exactS=sum(wS==0), exactT=sum(wT==0), nrows=nrow(EDF))
}

# ---- nulls for A (miRNA->gene) : decile matched on mRNA and on protein
nA_m <- make_null(A, DEC$mRNA); nA_p <- make_null(A, DEC$protein)
say("A null (mRNA-decile):   ",nrow(nA_m$NL)," pairs, ",round(nrow(nA_m$NL)/nrow(A),2)," per edge; exact-decile source ",nA_m$exactS,"/",nrow(A)," target ",nA_m$exactT,"/",nrow(A))
say("A null (protein-decile):",nrow(nA_p$NL)," pairs, ",round(nrow(nA_p$NL)/nrow(A),2)," per edge; exact-decile source ",nA_p$exactS,"/",nrow(A)," target ",nA_p$exactT,"/",nrow(A))
nB_m <- make_null(B, DEC$mRNA); nB_p <- make_null(B, DEC$protein)
say("B null (mRNA-decile):   ",nrow(nB_m$NL)," pairs; exact source ",nB_m$exactS,"/",nrow(B)," target ",nB_m$exactT,"/",nrow(B))
say("B null (protein-decile):",nrow(nB_p$NL)," pairs; exact source ",nB_p$exactS,"/",nrow(B)," target ",nB_p$exactT,"/",nrow(B))

eval_null <- function(NL, Xsrc, Ytgt) { r <- rho_pair(Xsrc, NL$rnd_source, Ytgt, NL$rnd_target); NL$rho <- r$rho; NL$n <- r$n; NL }
nA_m$NL <- eval_null(nA_m$NL, Mm, Rm)   # miRNA vs mRNA
nA_p$NL <- eval_null(nA_p$NL, Mm, Pm)   # miRNA vs protein
nB_m$NL <- eval_null(nB_m$NL, Rm, Rm)
nB_p$NL <- eval_null(nB_p$NL, Pm, Pm)
say("null rho computable: A/mRNA ",sum(!is.na(nA_m$NL$rho)),"  A/protein ",sum(!is.na(nA_p$NL$rho)),
    "  B/mRNA ",sum(!is.na(nB_m$NL$rho)),"  B/protein ",sum(!is.na(nB_p$NL$rho)))
write.csv(rbind(cbind(layer="A_miRNA_mRNA",  nA_m$NL), cbind(layer="A_miRNA_protein",nA_p$NL),
                cbind(layer="B_TF_mRNA",     nB_m$NL), cbind(layer="B_TF_protein",   nB_p$NL)),
          file.path(OUT,"cptac_null_pairs.csv"), row.names=FALSE)
say("WROTE cptac_null_pairs.csv rows=", nrow(nA_m$NL)+nrow(nA_p$NL)+nrow(nB_m$NL)+nrow(nB_p$NL))

## ---- summary builder --------------------------------------------------------
summarise <- function(label, EDF, rows, assay, rho_col, conc_col, NLobj) {
  rows <- rows[!is.na(EDF[[rho_col]][rows])]
  if (length(rows) < 3) return(NULL)
  rr <- EDF[[rho_col]][rows]
  nl <- NLobj$NL[NLobj$NL$edge_row %in% rows & !is.na(NLobj$NL$rho), ]
  ps <- EDF$predicted_sign[rows]
  nulconc <- if (all(!is.na(ps))) {
      psmap <- setNames(EDF$predicted_sign[rows], as.character(rows))
      mean(sign(nl$rho)==psmap[as.character(nl$edge_row)], na.rm=TRUE)
    } else NA_real_
  cc <- EDF[[conc_col]][rows]; cc <- cc[!is.na(cc)]
  k <- sum(cc); n <- length(cc)
  bt <- if (n>0 && !is.na(nulconc) && nulconc>0 && nulconc<1) binom.test(k,n,p=nulconc,alternative="greater")$p.value else NA_real_
  bt5<- if (n>0) binom.test(k,n,p=0.5)$p.value else NA_real_
  wt <- suppressWarnings(wilcox.test(rr, nl$rho)$p.value)
  data.frame(stratum=label, assay=assay, n_edges=length(rows), n_null_pairs=nrow(nl),
             mean_rho_real=mean(rr), mean_rho_null=mean(nl$rho),
             median_rho_real=median(rr), median_rho_null=median(nl$rho),
             delta_mean_rho=mean(rr)-mean(nl$rho),
             AUC_signed_rho=auc_fun(ps*rr, sign(nl$rho)*0 + nl$rho*psmapv(EDF,rows,nl)),
             n_edges_conc_tested=n, n_concordant=k,
             concordance_rate=if(n>0) k/n else NA_real_, null_concordance_rate=nulconc,
             excess_over_null=if(n>0) k/n-nulconc else NA_real_,
             binom_p_vs_null=bt, binom_p_vs_0.5=bt5, wilcox_p_rho_vs_null=wt,
             n_FDR05_correct_sign=sum(EDF[[conc_col]][rows] & EDF[[sub("rho","fdr",rho_col)]][rows]<0.05, na.rm=TRUE),
             stringsAsFactors=FALSE)
}
psmapv <- function(EDF, rows, nl) { m <- setNames(EDF$predicted_sign[rows], as.character(rows)); unname(m[as.character(nl$edge_row)]) }
# attach fdr columns onto A and B for the FDR count
A$fdr_mRNA <- p.adjust(A$p_mRNA,"BH"); A$fdr_protein <- p.adjust(A$p_protein,"BH")
B$fdr_mRNA <- p.adjust(B$p_mRNA,"BH"); B$fdr_protein <- p.adjust(B$p_protein,"BH")

res <- list()
strat_A <- list(miRNA_target_ALL = seq_len(nrow(A)),
                miRNA_target_strong = which(A$evidence_tier=="strong"),
                miRNA_target_weak = which(A$evidence_tier=="weak"),
                miRNA_target_predicted_only = which(A$evidence_tier=="predicted_only"))
for (nm in names(strat_A)) {
  res[[paste0(nm,"_mRNA")]]    <- summarise(nm, A, strat_A[[nm]], "mRNA",   "rho_mRNA",   "conc_mRNA",   nA_m)
  res[[paste0(nm,"_protein")]] <- summarise(nm, A, strat_A[[nm]], "protein","rho_protein","conc_protein",nA_p)
}
strat_B <- list(TF_target_Activation = which(B$trrust_mode=="Activation"),
                TF_target_Repression = which(B$trrust_mode=="Repression"),
                TF_target_ALL_signed = which(B$trrust_mode %in% c("Activation","Repression")))
for (nm in names(strat_B)) {
  res[[paste0(nm,"_mRNA")]]    <- summarise(nm, B, strat_B[[nm]], "mRNA",   "rho_mRNA",   "conc_mRNA",   nB_m)
  res[[paste0(nm,"_protein")]] <- summarise(nm, B, strat_B[[nm]], "protein","rho_protein","conc_protein",nB_p)
}
SUM <- do.call(rbind, res)
write.csv(SUM, file.path(OUT,"cptac_concordance_summary.csv"), row.names=FALSE)
say("WROTE cptac_concordance_summary.csv rows=", nrow(SUM))
for (i in seq_len(nrow(SUM))) say(sprintf("%-30s %-7s n=%4d  meanRho=%+.4f (null %+.4f)  conc=%.4f (null %.4f) binomP=%.3g",
    SUM$stratum[i],SUM$assay[i],SUM$n_edges[i],SUM$mean_rho_real[i],SUM$mean_rho_null[i],
    SUM$concordance_rate[i],SUM$null_concordance_rate[i],SUM$binom_p_vs_null[i]))

## paired mRNA vs protein comparison on identical edges
pa <- A[!is.na(A$rho_mRNA) & !is.na(A$rho_protein),]
wt <- wilcox.test(pa$rho_mRNA, pa$rho_protein, paired=TRUE)
say("A) paired edges with BOTH rho: ",nrow(pa),"  mean rho mRNA=",signif(mean(pa$rho_mRNA),4),
    " protein=",signif(mean(pa$rho_protein),4)," Wilcoxon paired p=",signif(wt$p.value,4))
pb <- B[!is.na(B$rho_mRNA) & !is.na(B$rho_protein) & B$trrust_mode %in% c("Activation","Repression"),]
wtb <- wilcox.test(pb$rho_mRNA, pb$rho_protein, paired=TRUE)
say("B) paired signed TF edges with BOTH rho: ",nrow(pb),"  mean rho mRNA=",signif(mean(pb$rho_mRNA),4),
    " protein=",signif(mean(pb$rho_protein),4)," Wilcoxon paired p=",signif(wtb$p.value,4))
PAIRED <- data.frame(layer=c("miRNA_target","TF_target_signed"), n_edges=c(nrow(pa),nrow(pb)),
  mean_rho_mRNA=c(mean(pa$rho_mRNA),mean(pb$rho_mRNA)), mean_rho_protein=c(mean(pa$rho_protein),mean(pb$rho_protein)),
  conc_mRNA=c(mean(pa$conc_mRNA),mean(pb$conc_mRNA,na.rm=TRUE)), conc_protein=c(mean(pa$conc_protein),mean(pb$conc_protein,na.rm=TRUE)),
  wilcox_paired_p=c(wt$p.value,wtb$p.value))
write.csv(PAIRED, file.path(OUT,"cptac_paired_mRNA_vs_protein.csv"), row.names=FALSE)
say("WROTE cptac_paired_mRNA_vs_protein.csv rows=",nrow(PAIRED))

saveRDS(list(A=A,B=B,SUM=SUM,gene_both=gene_both,mir_ok=mir_ok), file.path(OUT,"cptac_AB_workspace.rds"))
say("=== 21 DONE ===")
close(con)
