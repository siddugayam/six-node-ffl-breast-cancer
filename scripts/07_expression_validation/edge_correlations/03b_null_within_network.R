## 4th null construction: random pairs drawn from the NETWORK'S OWN node set,
## decile-matched.  (Diagnostic: this is what reproduces the published null.)
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
set.seed(11); K <- 10L
cat("### K =", K, "null draws per edge (within-network pool)\n")
gexp <- readRDS("data/brca_gene_expr.rds"); mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
samp <- sort(intersect(intersect(colnames(gexp),colnames(mexp)), tum)); N <- length(samp)
G <- gexp[,samp]; M <- mexp[,samp]; G <- G[rowSds(G)>0,]; M <- M[rowSds(M)>0,]
zr <- function(X){R <- t(apply(X,1,rank)); (R-rowMeans(R))/rowSds(R)}
ZG <- zr(G); ZM <- zr(M); dn <- N-1
E <- read.delim("data/canonical_edges.tsv", stringsAsFactors=FALSE)
tier <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors=FALSE)
tft <- read.delim("data/layer_TF_target.tsv", stringsAsFactors=FALSE)
tfm <- read.delim("data/layer_TF_miRNA.tsv", stringsAsFactors=FALSE)
E$tier <- ifelse(E$edge_type=="miRNA_target", tier$tier[match(paste(E$source,E$target),paste(tier$source,tier$target))], NA)
E$mode <- NA_character_
i <- E$edge_type=="TF_target"; E$mode[i] <- tft$mode[match(paste(E$source[i],E$target[i]),paste(tft$source,tft$target))]
i <- E$edge_type=="TF_miRNA";  E$mode[i] <- tfm$mode[match(paste(E$source[i],E$target[i]),paste(tfm$source,tfm$target))]
E$sty <- ifelse(E$source_type=="miRNA","miRNA","gene"); E$tty <- ifelse(E$target_type=="miRNA","miRNA","gene")
E$si <- ifelse(E$sty=="miRNA", match(E$source,rownames(M)), match(E$source,rownames(G)))
E$ti <- ifelse(E$tty=="miRNA", match(E$target,rownames(M)), match(E$target,rownames(G)))
E <- E[!is.na(E$si)&!is.na(E$ti)&E$source!=E$target,]
cat("measurable edges:", nrow(E), "\n")
E$exp <- NA_real_; E$exp[E$edge_type=="miRNA_target"] <- -1
E$exp[E$edge_type %in% c("TF_target","TF_miRNA") & E$mode=="Activation"] <- 1
E$exp[E$edge_type %in% c("TF_target","TF_miRNA") & E$mode=="Repression"] <- -1
rho <- function(si,sty,ti,tty){o <- numeric(length(si))
  for(a in c("miRNA","gene")) for(b in c("miRNA","gene")){k <- which(sty==a&tty==b); if(!length(k))next
    Za <- if(a=="miRNA")ZM else ZG; Zb <- if(b=="miRNA")ZM else ZG
    o[k] <- rowSums(Za[si[k],,drop=FALSE]*Zb[ti[k],,drop=FALSE])/dn}; o}
E$rho <- rho(E$si,E$sty,E$ti,E$tty)
## restrict candidate pools to network nodes present in the matrices
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)
ng <- which(rownames(G) %in% nodes$name[nodes$type!="miRNA"])
nm <- which(rownames(M) %in% nodes$name[nodes$type=="miRNA"])
cat("network gene pool:", length(ng), " miRNA pool:", length(nm), "\n")
qg <- cut(rowMeans(G)[ng], quantile(rowMeans(G)[ng],seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
qm <- cut(rowMeans(M)[nm], quantile(rowMeans(M)[nm],seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
pg <- split(ng,qg); pm <- split(nm,qm)
decg <- setNames(qg, ng); decm <- setNames(qm, nm)
realset <- new.env(hash=TRUE, size=2e4)
for(kk in paste(E$source,E$target)) assign(kk, TRUE, envir=realset)
gname <- rownames(G); mname <- rownames(M)
REPS <- 5L
allout <- NULL
for (rep_i in seq_len(REPS)){
set.seed(100+rep_i)
n <- nrow(E); si <- integer(n*K); ti <- integer(n*K)
for(r in seq_len(n)){
  ds <- if(E$sty[r]=="miRNA") decm[as.character(E$si[r])] else decg[as.character(E$si[r])]
  dt <- if(E$tty[r]=="miRNA") decm[as.character(E$ti[r])] else decg[as.character(E$ti[r])]
  ps <- if(E$sty[r]=="miRNA") pm[[as.character(ds)]] else pg[[as.character(ds)]]
  pt <- if(E$tty[r]=="miRNA") pm[[as.character(dt)]] else pg[[as.character(dt)]]
  j <- ((r-1)*K+1):(r*K); got <- 0L; tries <- 0L
  while(got < K && tries < 500L){ tries <- tries+1L
    a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
    if(E$sty[r]==E$tty[r] && a==b) next
    an <- if(E$sty[r]=="miRNA") mname[a] else gname[a]
    bn <- if(E$tty[r]=="miRNA") mname[b] else gname[b]
    if(exists(paste(an,bn), envir=realset, inherits=FALSE)) next
    got <- got+1L; si[j[got]] <- a; ti[j[got]] <- b }
  while(got < K){ a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
    if(E$sty[r]==E$tty[r] && a==b) next
    got <- got+1L; si[j[got]] <- a; ti[j[got]] <- b }
}
idx <- rep(seq_len(n),each=K); sty <- rep(E$sty,each=K); tty <- rep(E$tty,each=K)
keep <- rep(TRUE, length(si))
nr <- rho(si,ti=ti,sty=sty,tty=tty); nex <- E$exp[idx]
sets <- list(TF_target_Activation=E$edge_type=="TF_target"&E$mode=="Activation",
             TF_target_Repression=E$edge_type=="TF_target"&E$mode=="Repression",
             TF_miRNA_Activation=E$edge_type=="TF_miRNA"&E$mode=="Activation",
             TF_miRNA_Repression=E$edge_type=="TF_miRNA"&E$mode=="Repression",
             miRNA_target_all=E$edge_type=="miRNA_target",
             miRNA_target_strong=E$edge_type=="miRNA_target"&E$tier=="strong",
             miRNA_target_weak=E$edge_type=="miRNA_target"&E$tier=="weak",
             miRNA_target_predicted_only=E$edge_type=="miRNA_target"&E$tier=="predicted_only")
out <- NULL
for(nmm in names(sets)){ s <- which(sets[[nmm]] & !is.na(E$exp)); if(!length(s)) next
  k <- (idx %in% s) & keep
  p0 <- mean(sign(nr[k])==nex[k]); x1 <- sum(sign(E$rho[s])==E$exp[s])
  out <- rbind(out, data.frame(set=nmm, n=length(s), conc=x1/length(s), conc_null=p0,
    p_prop=prop.test(c(x1,sum(sign(nr[k])==nex[k])), c(length(s),sum(k)))$p.value,
    p_binom=binom.test(x1,length(s),p=p0)$p.value,
    mean_rho=mean(E$rho[s]), mean_rho_null=mean(nr[k]),
    p_wilcox=wilcox.test(E$rho[s],nr[k])$p.value))
}
out$rep <- rep_i; allout <- rbind(allout, out)
}
agg <- do.call(rbind, lapply(split(allout, allout$set), function(d)
  data.frame(set=d$set[1], n=d$n[1], conc=d$conc[1],
             conc_null_mean=mean(d$conc_null), conc_null_sd=sd(d$conc_null),
             p_binom_mean=mean(d$p_binom), p_binom_min=min(d$p_binom), p_binom_max=max(d$p_binom),
             mean_rho=d$mean_rho[1], mean_rho_null=mean(d$mean_rho_null),
             p_wilcox_mean=mean(d$p_wilcox), p_wilcox_max=max(d$p_wilcox))))
cat("\n### within-network-node decile-matched null, REPS =", REPS, "x K =", K, "draws/edge\n")
print(agg, digits=4, row.names=FALSE)
write.csv(allout, "results/v2/v2_null_within_network_pool_reps.csv", row.names=FALSE)
write.csv(agg, "results/v2/v2_null_within_network_pool.csv", row.names=FALSE)
cat("DONE 03b\n")
