## stroma_adjusted_concordance.R -- sign concordance of network edges by evidence tier, before and after adjustment for
## each of the 45 stromal estimates of the mediation analysis (ABSOLUTE purity and the ESTIMATE stromal score among
## them), against the five-draw within-network matched null, and with the targets correlated with the Farmer
## stroma-related score left out (see SETTINGS.md).
## The data preparation, the null draws (seeds 501-505, 10 pairs per edge) and the test are copied unchanged from
## scripts/07_expression_validation/edge_correlations/04_caf_adjusted_vs_matched_null.R; the script first checks that
## it reproduces that script's table.
## Root of the repository: environment variable FFL_ROOT. Outputs go to the folder of this script; no sample-level
## values are written.
## usage: FFL_ROOT=<repository checkout> Rscript stroma_adjusted_concordance.R
suppressPackageStartupMessages({library(matrixStats)})
ROOT <- Sys.getenv("FFL_ROOT"); stopifnot(nzchar(ROOT))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
setwd(ROOT)
K <- 10L; REPS <- 5L
gexp <- readRDS("data/brca_gene_expr.rds"); mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
samp <- sort(intersect(intersect(colnames(gexp),colnames(mexp)), tum)); N <- length(samp)
G <- gexp[,samp]; M <- mexp[,samp]; G <- G[rowSds(G)>0,]; M <- M[rowSds(M)>0,]
zr <- function(X){R <- t(apply(X,1,rank)); (R-rowMeans(R))/rowSds(R)}
ZG <- zr(G); ZM <- zr(M); dn <- N-1
cs <- readRDS("results/expression_validation/caf_scores.rds"); stopifnot(identical(cs$samp, samp))
zc <- function(v){r <- rank(v); (r-mean(r))/sd(r)}
mk <- function(v){ cg <- as.vector(ZG%*%zc(v))/dn; cm <- as.vector(ZM%*%zc(v))/dn
  function(r,si,sty,ti,tty){ a <- ifelse(sty=="miRNA",cm[si],cg[si]); b <- ifelse(tty=="miRNA",cm[ti],cg[ti])
    (r-a*b)/sqrt((1-a^2)*(1-b^2)) } }
pA <- mk(cs$cafA); pB <- mk(cs$cafB)
E <- read.delim("data/canonical_edges.tsv", stringsAsFactors=FALSE)
tier <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors=FALSE)
tft <- read.delim("data/layer_TF_target.tsv", stringsAsFactors=FALSE); tfm <- read.delim("data/layer_TF_miRNA.tsv", stringsAsFactors=FALSE)
E$tier <- ifelse(E$edge_type=="miRNA_target", tier$tier[match(paste(E$source,E$target),paste(tier$source,tier$target))], NA)
E$mode <- NA_character_
i <- E$edge_type=="TF_target"; E$mode[i] <- tft$mode[match(paste(E$source[i],E$target[i]),paste(tft$source,tft$target))]
i <- E$edge_type=="TF_miRNA";  E$mode[i] <- tfm$mode[match(paste(E$source[i],E$target[i]),paste(tfm$source,tfm$target))]
E$sty <- ifelse(E$source_type=="miRNA","miRNA","gene"); E$tty <- ifelse(E$target_type=="miRNA","miRNA","gene")
E$si <- ifelse(E$sty=="miRNA", match(E$source,rownames(M)), match(E$source,rownames(G)))
E$ti <- ifelse(E$tty=="miRNA", match(E$target,rownames(M)), match(E$target,rownames(G)))
E <- E[!is.na(E$si)&!is.na(E$ti)&E$source!=E$target,]
E$exp <- NA_real_; E$exp[E$edge_type=="miRNA_target"] <- -1
E$exp[E$edge_type%in%c("TF_target","TF_miRNA")&E$mode=="Activation"] <- 1
E$exp[E$edge_type%in%c("TF_target","TF_miRNA")&E$mode=="Repression"] <- -1
rho <- function(si,sty,ti,tty){o <- numeric(length(si))
  for(a in c("miRNA","gene")) for(b in c("miRNA","gene")){k <- which(sty==a&tty==b); if(!length(k))next
    Za <- if(a=="miRNA")ZM else ZG; Zb <- if(b=="miRNA")ZM else ZG
    o[k] <- rowSums(Za[si[k],,drop=FALSE]*Zb[ti[k],,drop=FALSE])/dn}; o}
E$rho <- rho(E$si,E$sty,E$ti,E$tty)
E$rA <- pA(E$rho,E$si,E$sty,E$ti,E$tty); E$rB <- pB(E$rho,E$si,E$sty,E$ti,E$tty)
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)
ng <- which(rownames(G)%in%nodes$name[nodes$type!="miRNA"]); nmi <- which(rownames(M)%in%nodes$name[nodes$type=="miRNA"])
qg <- cut(rowMeans(G)[ng], quantile(rowMeans(G)[ng],seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
qm <- cut(rowMeans(M)[nmi], quantile(rowMeans(M)[nmi],seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
pg <- split(ng,qg); pm <- split(nmi,qm); decg <- setNames(qg,ng); decm <- setNames(qm,nmi)
realset <- new.env(hash=TRUE,size=2e4); for(kk in paste(E$source,E$target)) assign(kk,TRUE,envir=realset)
gn <- rownames(G); mn <- rownames(M)
sets <- list(TF_target_Activation=E$edge_type=="TF_target"&E$mode=="Activation",
             TF_target_Repression=E$edge_type=="TF_target"&E$mode=="Repression",
             TF_miRNA_Activation=E$edge_type=="TF_miRNA"&E$mode=="Activation",
             TF_miRNA_Repression=E$edge_type=="TF_miRNA"&E$mode=="Repression",
             miRNA_target_all=E$edge_type=="miRNA_target",
             miRNA_target_strong=E$edge_type=="miRNA_target"&E$tier=="strong",
             miRNA_target_weak=E$edge_type=="miRNA_target"&E$tier=="weak",
             miRNA_target_predicted_only=E$edge_type=="miRNA_target"&E$tier=="predicted_only")
DRAWS <- list(); allout <- NULL
for(rep_i in seq_len(REPS)){
  set.seed(500+rep_i); n <- nrow(E); si <- integer(n*K); ti <- integer(n*K)
  for(r in seq_len(n)){
    ds <- if(E$sty[r]=="miRNA") decm[as.character(E$si[r])] else decg[as.character(E$si[r])]
    dt <- if(E$tty[r]=="miRNA") decm[as.character(E$ti[r])] else decg[as.character(E$ti[r])]
    ps <- if(E$sty[r]=="miRNA") pm[[as.character(ds)]] else pg[[as.character(ds)]]
    pt <- if(E$tty[r]=="miRNA") pm[[as.character(dt)]] else pg[[as.character(dt)]]
    j <- ((r-1)*K+1):(r*K); got <- 0L; tries <- 0L
    while(got<K && tries<500L){ tries <- tries+1L
      a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
      if(E$sty[r]==E$tty[r] && a==b) next
      an <- if(E$sty[r]=="miRNA") mn[a] else gn[a]; bn <- if(E$tty[r]=="miRNA") mn[b] else gn[b]
      if(exists(paste(an,bn),envir=realset,inherits=FALSE)) next
      got <- got+1L; si[j[got]] <- a; ti[j[got]] <- b }
    while(got<K){ a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
      if(E$sty[r]==E$tty[r] && a==b) next; got <- got+1L; si[j[got]] <- a; ti[j[got]] <- b }
  }
  idx <- rep(seq_len(n),each=K); sty <- rep(E$sty,each=K); tty <- rep(E$tty,each=K)
  DRAWS[[rep_i]] <- list(si=si, ti=ti, idx=idx, sty=sty, tty=tty)
  nr <- rho(si,sty,ti,tty); nA <- pA(nr,si,sty,ti,tty); nB <- pB(nr,si,sty,ti,tty); nex <- E$exp[idx]
  for(nmm in names(sets)){ s <- which(sets[[nmm]]&!is.na(E$exp)); if(!length(s))next; k <- idx%in%s
    row <- data.frame(rep=rep_i, set=nmm, n=length(s))
    for(tag in c("raw","cafA","cafB")){
      rv <- switch(tag, raw=E$rho[s], cafA=E$rA[s], cafB=E$rB[s])
      nv <- switch(tag, raw=nr[k],    cafA=nA[k],   cafB=nB[k])
      x1 <- sum(sign(rv)==E$exp[s]); p0 <- mean(sign(nv)==nex[k])
      row[[paste0("conc_",tag)]] <- x1/length(s); row[[paste0("null_",tag)]] <- p0
      row[[paste0("p_",tag)]] <- binom.test(x1,length(s),p=p0)$p.value
    }
    allout <- rbind(allout,row) }
}
agg <- do.call(rbind, lapply(split(allout,allout$set), function(d) data.frame(set=d$set[1], n=d$n[1],
  conc_raw=d$conc_raw[1], null_raw=mean(d$null_raw), p_raw=mean(d$p_raw),
  conc_cafA=d$conc_cafA[1], null_cafA=mean(d$null_cafA), p_cafA=mean(d$p_cafA), p_cafA_max=max(d$p_cafA),
  conc_cafB=d$conc_cafB[1], null_cafB=mean(d$null_cafB), p_cafB=mean(d$p_cafB), p_cafB_max=max(d$p_cafB))))

## ---- control: the published table of 04 -------------------------------------------------------------------------
pub <- read.csv("results/expression_validation/caf_adjusted_vs_network_null.csv", stringsAsFactors=FALSE)
pub <- pub[match(agg$set, pub$set),]; num <- setdiff(colnames(agg), "set")
relmax <- sapply(num, function(c) max(abs(agg[[c]]-pub[[c]])/pmax(abs(pub[[c]]),1e-300)))
ctrl <- data.frame(column=num, max_relative_difference=relmax, reproduced=relmax<1e-9)
write.csv(ctrl, file.path(HERE, "control_reproduction.csv"), row.names=FALSE)
cat("control: 04 table reproduced:", all(ctrl$reproduced), "\n"); stopifnot(all(ctrl$reproduced))

## ---- the 45 stromal estimates, each alone -----------------------------------------------------------------------
FIB <- read.csv("results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv", check.names=FALSE)
FIB <- FIB[match(samp, FIB$sample), ]; stopifnot(identical(FIB$sample, samp)); EST <- setdiff(colnames(FIB), "sample")
stopifnot(length(EST) == 45)
farmer <- read.csv("results/prioritisation_and_stroma/farmer_pergene_rho.csv", stringsAsFactors=FALSE)
FARMER_SET <- farmer$gene[farmer$rho_vs_farmer > 0.5]
cat("Farmer-correlated genes (rho > 0.5):", length(FARMER_SET), "; network miRNA targets among them:",
    length(intersect(FARMER_SET, unique(E$target[E$edge_type=="miRNA_target"]))), "\n")
keep_all <- rep(TRUE, nrow(E)); keep_lo <- !(E$target %in% FARMER_SET)
## rows of G and M that any edge or null pair uses (network nodes)
rowsG <- sort(unique(c(ng, E$si[E$sty=="gene"], E$ti[E$tty=="gene"]))); rowsM <- sort(unique(c(nmi, E$si[E$sty=="miRNA"], E$ti[E$tty=="miRNA"])))
res <- list()
for(est in EST){
  v <- FIB[[est]]; ok <- is.finite(v); ns <- sum(ok)
  if(ns == N){ ZGs <- ZG; ZMs <- ZM; dns <- dn } else {
    ZGs <- matrix(NA_real_, nrow(G), ns); ZMs <- matrix(NA_real_, nrow(M), ns)
    ZGs[rowsG,] <- zr(G[rowsG, ok, drop=FALSE]); ZMs[rowsM,] <- zr(M[rowsM, ok, drop=FALSE]); dns <- ns-1 }
  rhoS <- function(si,sty,ti,tty){o <- numeric(length(si))
    for(a in c("miRNA","gene")) for(b in c("miRNA","gene")){k <- which(sty==a&tty==b); if(!length(k))next
      Za <- if(a=="miRNA")ZMs else ZGs; Zb <- if(b=="miRNA")ZMs else ZGs
      o[k] <- rowSums(Za[si[k],,drop=FALSE]*Zb[ti[k],,drop=FALSE])/dns}; o}
  zv <- zc(v[ok]); cg <- as.vector(ZGs[rowsG,,drop=FALSE] %*% zv)/dns; cm <- as.vector(ZMs[rowsM,,drop=FALSE] %*% zv)/dns
  CG <- numeric(nrow(G)); CG[rowsG] <- cg; CM <- numeric(nrow(M)); CM[rowsM] <- cm
  padj <- function(r,si,sty,ti,tty){ a <- ifelse(sty=="miRNA",CM[si],CG[si]); b <- ifelse(tty=="miRNA",CM[ti],CG[ti])
    (r-a*b)/sqrt((1-a^2)*(1-b^2)) }
  er <- rhoS(E$si,E$sty,E$ti,E$tty); ea <- padj(er,E$si,E$sty,E$ti,E$tty)
  for(rep_i in seq_len(REPS)){
    D <- DRAWS[[rep_i]]; nr <- rhoS(D$si,D$sty,D$ti,D$tty); na <- padj(nr,D$si,D$sty,D$ti,D$tty); nex <- E$exp[D$idx]
    for(variant in c("all edges","Farmer-correlated targets left out")){
      kp <- if(variant=="all edges") keep_all else keep_lo
      for(nmm in names(sets)){ s <- which(sets[[nmm]] & !is.na(E$exp) & kp); if(!length(s)) next; k <- D$idx %in% s
        for(tag in c("unadjusted_same_tumours","adjusted")){
          rv <- if(tag=="adjusted") ea[s] else er[s]; nv <- if(tag=="adjusted") na[k] else nr[k]
          x1 <- sum(sign(rv)==E$exp[s]); p0 <- mean(sign(nv)==nex[k])
          res[[length(res)+1]] <- data.frame(estimate=est, n_tumours=ns, edges=variant, set=nmm, n_edges=length(s),
            rep=rep_i, analysis=tag, conc=x1/length(s), null=p0, p=binom.test(x1,length(s),p=p0)$p.value)
        } } } }
  cat(sprintf("%-36s tumours %d done\n", est, ns))
}
R <- do.call(rbind, res)
A <- do.call(rbind, lapply(split(R, list(R$estimate, R$edges, R$set, R$analysis), drop=TRUE), function(d)
  data.frame(estimate=d$estimate[1], n_tumours=d$n_tumours[1], edges=d$edges[1], set=d$set[1], n_edges=d$n_edges[1],
             analysis=d$analysis[1], conc=d$conc[1], null=mean(d$null), p_mean=mean(d$p), p_max=max(d$p))))
A <- A[order(A$edges, A$set, A$analysis, A$estimate), ]
write.csv(A, file.path(HERE, "concordance_adjusted_by_estimate.csv"), row.names=FALSE)

## ---- summary ------------------------------------------------------------------------------------------------------
S <- list()
for(variant in unique(A$edges)) for(nmm in names(sets)){
  a <- A[A$edges==variant & A$set==nmm, ]; if(!nrow(a)) next
  un <- a[a$analysis=="unadjusted_same_tumours" & a$estimate=="ESTIMATE_StromalScore", ]
  ad <- a[a$analysis=="adjusted", ]
  one <- function(e) ad[ad$estimate==e, ]
  ab <- one("ABSOLUTE_nontumour"); es <- one("ESTIMATE_StromalScore")
  aun <- a[a$analysis=="unadjusted_same_tumours" & a$estimate=="ABSOLUTE_nontumour", ]
  S[[length(S)+1]] <- data.frame(edges=variant, set=nmm, n_edges=un$n_edges,
    unadj_conc=un$conc, unadj_null=un$null, unadj_p=un$p_mean, unadj_p_max=un$p_max,
    ABSOLUTE_tumours=ab$n_tumours, ABSOLUTE_unadj_conc=aun$conc, ABSOLUTE_unadj_null=aun$null, ABSOLUTE_unadj_p=aun$p_mean,
    ABSOLUTE_conc=ab$conc, ABSOLUTE_null=ab$null, ABSOLUTE_p=ab$p_mean, ABSOLUTE_p_max=ab$p_max,
    ESTIMATE_conc=es$conc, ESTIMATE_null=es$null, ESTIMATE_p=es$p_mean, ESTIMATE_p_max=es$p_max,
    all45_conc_median=median(ad$conc), all45_conc_min=min(ad$conc), all45_conc_max=max(ad$conc),
    all45_excess_median=median(ad$conc-ad$null), all45_excess_min=min(ad$conc-ad$null), all45_excess_max=max(ad$conc-ad$null),
    all45_p_median=median(ad$p_mean), all45_p_min=min(ad$p_mean), all45_p_max=max(ad$p_mean),
    all45_n_p_below_0.05=sum(ad$p_mean<0.05), all45_n_above_null=sum(ad$conc>ad$null))
}
S <- do.call(rbind, S)
write.csv(S, file.path(HERE, "concordance_summary.csv"), row.names=FALSE)
print(S[grepl("miRNA", S$set), c("edges","set","n_edges","unadj_conc","unadj_null","unadj_p","ABSOLUTE_conc","ABSOLUTE_null",
        "ABSOLUTE_p","ESTIMATE_conc","ESTIMATE_null","ESTIMATE_p","all45_conc_median","all45_p_median","all45_n_p_below_0.05")], digits=3)
cat("DONE\n")
