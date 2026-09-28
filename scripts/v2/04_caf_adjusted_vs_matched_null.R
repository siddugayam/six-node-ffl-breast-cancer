## N9 under the SAME (within-network-node, decile-matched) null used by the paper.
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
K <- 10L; REPS <- 5L
cat("### within-network-node decile-matched null; REPS =",REPS,"x K =",K,"draws/edge\n")
gexp <- readRDS("data/brca_gene_expr.rds"); mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
samp <- sort(intersect(intersect(colnames(gexp),colnames(mexp)), tum)); N <- length(samp)
G <- gexp[,samp]; M <- mexp[,samp]; G <- G[rowSds(G)>0,]; M <- M[rowSds(M)>0,]
zr <- function(X){R <- t(apply(X,1,rank)); (R-rowMeans(R))/rowSds(R)}
ZG <- zr(G); ZM <- zr(M); dn <- N-1
cs <- readRDS("results/v2/v2_caf_scores.rds"); stopifnot(identical(cs$samp, samp))
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
allout <- NULL
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
print(agg, digits=4, row.names=FALSE)
write.csv(agg,"results/v2/v2_caf_adjusted_vs_network_null.csv", row.names=FALSE)
cat("DONE 04\n")
