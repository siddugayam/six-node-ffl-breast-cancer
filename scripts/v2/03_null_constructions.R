## Sensitivity of the N8 null to how the decile-matched random pairs are drawn.
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
OUT <- "results/v2"; set.seed(11); K <- 20L
cat("### K =", K, "null draws per edge, per construction\n")
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
E$tier <- ifelse(E$edge_type=="miRNA_target",
                 tier$tier[match(paste(E$source,E$target), paste(tier$source,tier$target))], NA)
E$mode <- NA_character_
i <- E$edge_type=="TF_target"; E$mode[i] <- tft$mode[match(paste(E$source[i],E$target[i]),paste(tft$source,tft$target))]
i <- E$edge_type=="TF_miRNA";  E$mode[i] <- tfm$mode[match(paste(E$source[i],E$target[i]),paste(tfm$source,tfm$target))]
E$sty <- ifelse(E$source_type=="miRNA","miRNA","gene"); E$tty <- ifelse(E$target_type=="miRNA","miRNA","gene")
E$si <- ifelse(E$sty=="miRNA", match(E$source,rownames(M)), match(E$source,rownames(G)))
E$ti <- ifelse(E$tty=="miRNA", match(E$target,rownames(M)), match(E$target,rownames(G)))
E <- E[!is.na(E$si)&!is.na(E$ti)&E$source!=E$target,]
E$exp <- NA_real_; E$exp[E$edge_type=="miRNA_target"] <- -1
E$exp[E$edge_type %in% c("TF_target","TF_miRNA") & E$mode=="Activation"] <- 1
E$exp[E$edge_type %in% c("TF_target","TF_miRNA") & E$mode=="Repression"] <- -1
rho <- function(si,sty,ti,tty){o <- numeric(length(si))
  for(a in c("miRNA","gene")) for(b in c("miRNA","gene")){k <- which(sty==a&tty==b); if(!length(k))next
    Za <- if(a=="miRNA")ZM else ZG; Zb <- if(b=="miRNA")ZM else ZG
    o[k] <- rowSums(Za[si[k],,drop=FALSE]*Zb[ti[k],,drop=FALSE])/dn}; o}
E$rho <- rho(E$si,E$sty,E$ti,E$tty)
dg <- cut(rowMeans(G), quantile(rowMeans(G),seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
dm <- cut(rowMeans(M), quantile(rowMeans(M),seq(0,1,.1)), include.lowest=TRUE, labels=FALSE)
pg <- split(seq_along(dg),dg); pm <- split(seq_along(dm),dm)
pool <- function(ty,d) if(ty=="miRNA") pm[[as.character(d)]] else pg[[as.character(d)]]
dec  <- function(ty,i) if(ty=="miRNA") dm[i] else dg[i]
draw <- function(mode){
  n <- nrow(E); si <- integer(n*K); ti <- integer(n*K)
  for(r in seq_len(n)){ ps <- pool(E$sty[r],dec(E$sty[r],E$si[r])); pt <- pool(E$tty[r],dec(E$tty[r],E$ti[r]))
    j <- ((r-1)*K+1):(r*K)
    si[j] <- if(mode=="target_only") E$si[r] else ps[sample.int(length(ps),K,replace=TRUE)]
    ti[j] <- if(mode=="source_only") E$ti[r] else pt[sample.int(length(pt),K,replace=TRUE)]
  }
  list(si=si,ti=ti,sty=rep(E$sty,each=K),tty=rep(E$tty,each=K),idx=rep(seq_len(n),each=K))
}
sets <- list(TF_target_Activation=E$edge_type=="TF_target"&E$mode=="Activation",
             TF_target_Repression=E$edge_type=="TF_target"&E$mode=="Repression",
             TF_miRNA_Activation=E$edge_type=="TF_miRNA"&E$mode=="Activation",
             TF_miRNA_Repression=E$edge_type=="TF_miRNA"&E$mode=="Repression",
             miRNA_target_all=E$edge_type=="miRNA_target",
             miRNA_target_strong=E$edge_type=="miRNA_target"&E$tier=="strong",
             miRNA_target_weak=E$edge_type=="miRNA_target"&E$tier=="weak",
             miRNA_target_predicted_only=E$edge_type=="miRNA_target"&E$tier=="predicted_only")
out <- NULL
for (mode in c("both","source_only","target_only")){
  d <- draw(mode); nr <- rho(d$si,d$sty,d$ti,d$tty); nex <- E$exp[d$idx]
  for(nm in names(sets)){ s <- which(sets[[nm]] & !is.na(E$exp)); if(!length(s)) next
    k <- d$idx %in% s
    p1 <- mean(sign(E$rho[s])==E$exp[s]); p0 <- mean(sign(nr[k])==nex[k])
    out <- rbind(out, data.frame(null=mode, set=nm, n=length(s),
      conc=p1, conc_null=p0,
      p_prop=prop.test(c(sum(sign(E$rho[s])==E$exp[s]), sum(sign(nr[k])==nex[k])), c(length(s), sum(k)))$p.value,
      p_binom=binom.test(sum(sign(E$rho[s])==E$exp[s]), length(s), p=p0)$p.value,
      mean_rho=mean(E$rho[s]), mean_rho_null=mean(nr[k]),
      p_wilcox=wilcox.test(E$rho[s], nr[k])$p.value))
  }
}
print(out, digits=4)
write.csv(out, file.path(OUT,"v2_null_construction_sensitivity.csv"), row.names=FALSE)
cat("DONE 03\n")
