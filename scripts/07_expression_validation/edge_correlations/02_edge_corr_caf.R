## ==========================================================================
## Independent verification of N6 / N8 / N9  (tasks B and C)
## Edge-level Spearman correlations in TCGA-BRCA primary tumours,
## expression-decile-matched null, and non-circular CAF (stromal) adjustment.
## Written from scratch for this verification run.
## ==========================================================================
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
OUT <- "results/v2"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
set.seed(20260908)
K_NULL <- 20L                      # null pairs drawn per real edge (>=10 required)
cat("### null draws per edge K =", K_NULL, "\n")

## ---------------------------------------------------------------- data ----
gexp  <- readRDS("data/brca_gene_expr.rds")
mexp  <- readRDS("data/brca_mirna_expr_canonical.rds")
pheno <- readRDS("data/brca_pheno.rds")

tum <- pheno$sample[pheno$sample_type=="Primary Tumor"]
samp <- intersect(intersect(colnames(gexp), colnames(mexp)), tum)
samp <- sort(samp)
N <- length(samp)
cat("paired primary tumours with BOTH assays: n =", N, "\n")
G <- gexp[, samp, drop=FALSE]; M <- mexp[, samp, drop=FALSE]

## drop constant features (rho undefined)
G <- G[rowSds(G) > 0, , drop=FALSE]; M <- M[rowSds(M) > 0, , drop=FALSE]
cat("genes usable:", nrow(G), " miRNAs usable:", nrow(M), "\n")

## rank transform once; Spearman == Pearson on (tie-corrected) ranks
rankrows <- function(X) t(apply(X, 1, rank))
RG <- rankrows(G); RM <- rankrows(M)
zrows <- function(R){ mu <- rowMeans(R); sd <- rowSds(R); (R-mu)/sd }
ZG <- zrows(RG); ZM <- zrows(RM)
denom <- N - 1

## ------------------------------------------------------------ edge sets ---
E <- read.delim("data/canonical_edges.tsv", stringsAsFactors=FALSE)
tier <- read.delim("data/edge_evidence_tier.tsv", stringsAsFactors=FALSE)
tft  <- read.delim("data/layer_TF_target.tsv", stringsAsFactors=FALSE)
tfm  <- read.delim("data/layer_TF_miRNA.tsv", stringsAsFactors=FALSE)

E$key <- paste(E$source, E$target, E$edge_type, sep="|")
E$tier <- tier$tier[match(paste(E$source,E$target), paste(tier$source,tier$target))]
E$tier[E$edge_type!="miRNA_target"] <- NA
E$mode <- NA_character_
i <- E$edge_type=="TF_target"
E$mode[i] <- tft$mode[match(paste(E$source[i],E$target[i]), paste(tft$source,tft$target))]
i <- E$edge_type=="TF_miRNA"
E$mode[i] <- tfm$mode[match(paste(E$source[i],E$target[i]), paste(tfm$source,tfm$target))]
E$mode[is.na(E$mode) & E$edge_type %in% c("TF_target","TF_miRNA")] <- "Unknown"

cat("\n--- N6: evidence tiers of miRNA_target edges ---\n")
print(table(E$tier[E$edge_type=="miRNA_target"], useNA="ifany"))
cat("strong pct:", round(100*sum(E$tier=="strong",na.rm=TRUE)/sum(E$edge_type=="miRNA_target"),2), "\n")
cat("\n--- N7: TRRUST modes on the 770 TF_target edges ---\n")
print(table(E$mode[E$edge_type=="TF_target"], useNA="ifany"))
cat("repression pct:", round(100*sum(E$mode=="Repression" & E$edge_type=="TF_target",na.rm=TRUE)/770,2), "\n")
cat("--- TF_miRNA modes ---\n"); print(table(E$mode[E$edge_type=="TF_miRNA"], useNA="ifany"))

## which edges are measurable
srow <- function(nm, ty) ifelse(ty=="miRNA", match(nm, rownames(M)), match(nm, rownames(G)))
E$si <- srow(E$source, E$source_type); E$ti <- srow(E$target, E$target_type)
E$ok <- !is.na(E$si) & !is.na(E$ti) & !(E$source==E$target)
cat("\nmeasurable edges:", sum(E$ok), "of", nrow(E), "\n")

## ------------------------------------------------------- rho machinery ----
getZ <- function(ty) if (ty=="miRNA") ZM else ZG
rho_pairs <- function(si, sty, ti, tty){
  out <- numeric(length(si))
  for (a in c("miRNA","gene")) for (b in c("miRNA","gene")){
    k <- which(sty==a & tty==b); if(!length(k)) next
    Za <- getZ(a); Zb <- getZ(b)
    out[k] <- rowSums(Za[si[k],,drop=FALSE] * Zb[ti[k],,drop=FALSE]) / denom
  }
  out
}
## partial rho given a covariate z-scored rank vector zc
make_partial <- function(zc){
  RCg <- as.vector(ZG %*% zc)/denom   # rho(gene_i, C)
  RCm <- as.vector(ZM %*% zc)/denom
  names(RCg) <- rownames(ZG); names(RCm) <- rownames(ZM)
  function(rxy, si, sty, ti, tty){
    rxc <- ifelse(sty=="miRNA", RCm[si], RCg[si])
    ryc <- ifelse(tty=="miRNA", RCm[ti], RCg[ti])
    (rxy - rxc*ryc)/sqrt((1-rxc^2)*(1-ryc^2))
  }
}

## ----------------------------------------------- decile-matched null ------
## deciles on mean expression, computed separately within each assay
dec_g <- cut(rowMeans(G), breaks=quantile(rowMeans(G), probs=seq(0,1,0.1)),
             include.lowest=TRUE, labels=FALSE)
dec_m <- cut(rowMeans(M), breaks=quantile(rowMeans(M), probs=seq(0,1,0.1)),
             include.lowest=TRUE, labels=FALSE)
pool_g <- split(seq_len(nrow(G)), dec_g); pool_m <- split(seq_len(nrow(M)), dec_m)
getpool <- function(ty) if (ty=="miRNA") pool_m else pool_g
getdec  <- function(ty) if (ty=="miRNA") dec_m else dec_g

real_key <- new.env(hash=TRUE, size=2e4)
for (k in paste(E$source, E$target)) assign(k, TRUE, envir=real_key)
namef <- function(ty, i) if (ty=="miRNA") rownames(M)[i] else rownames(G)[i]

Eok <- E[E$ok, ]
cat("building decile-matched null: ", nrow(Eok), "edges x", K_NULL, "=",
    nrow(Eok)*K_NULL, "pairs\n")
nsi <- integer(nrow(Eok)*K_NULL); nti <- integer(nrow(Eok)*K_NULL)
nsty <- character(length(nsi)); ntty <- character(length(nsi))
nidx <- rep(seq_len(nrow(Eok)), each=K_NULL)
ptr <- 0L
for (r in seq_len(nrow(Eok))){
  sty <- Eok$source_type[r]; tty <- Eok$target_type[r]
  sty2 <- if (sty=="miRNA") "miRNA" else "gene"; tty2 <- if (tty=="miRNA") "miRNA" else "gene"
  ps <- getpool(sty2)[[ as.character(getdec(sty2)[Eok$si[r]]) ]]
  pt <- getpool(tty2)[[ as.character(getdec(tty2)[Eok$ti[r]]) ]]
  got <- 0L; tries <- 0L
  while (got < K_NULL && tries < 400L){
    tries <- tries + 1L
    a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
    if (sty2==tty2 && a==b) next
    if (exists(paste(namef(sty2,a), namef(tty2,b)), envir=real_key, inherits=FALSE)) next
    got <- got + 1L; ptr <- ptr + 1L
    nsi[ptr] <- a; nti[ptr] <- b; nsty[ptr] <- sty2; ntty[ptr] <- tty2
  }
  if (got < K_NULL){  # pad by repeating draws without the edge-exclusion
    while (got < K_NULL){
      a <- ps[sample.int(length(ps),1)]; b <- pt[sample.int(length(pt),1)]
      if (sty2==tty2 && a==b) next
      got <- got+1L; ptr <- ptr+1L
      nsi[ptr] <- a; nti[ptr] <- b; nsty[ptr] <- sty2; ntty[ptr] <- tty2
    }
  }
}
stopifnot(ptr == nrow(Eok)*K_NULL)
cat("null pairs built:", ptr, "\n")

## ------------------------------------------------------ CAF scores --------
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)
node_names <- unique(nodes$name)
gmt <- readLines(system.file("extdata","SI_geneset.gmt", package="estimate"))
strom <- strsplit(gmt[grep("StromalSignature", gmt)], "\t")[[1]][-c(1,2)]
strom <- strom[strom != ""]
cat("\nESTIMATE stromal signature:", length(strom), "genes\n")
s1 <- setdiff(strom, grep("^COL", strom, value=TRUE))
cat("  after removing collagens (", length(strom)-length(s1), "):", length(s1), "\n")
s2 <- setdiff(s1, node_names)
cat("  after removing network nodes (", length(s1)-length(s2), " removed:",
    paste(intersect(s1,node_names), collapse=","), "):", length(s2), "\n")
sigA <- intersect(s2, rownames(G))
cat("  present in expression matrix -> CAF-A uses", length(sigA), "genes\n")

panel <- c("DCN","LUM","FAP","PDGFRB","THY1","POSTN")
p2 <- setdiff(panel, node_names)
cat("marker panel:", paste(panel, collapse=","), "-> excluded as network nodes:",
    paste(intersect(panel,node_names), collapse=","), "\n")
sigB <- intersect(p2, rownames(G))
cat("  CAF-B uses", length(sigB), "genes:", paste(sigB, collapse=","), "\n")
stopifnot(!any(grepl("^COL", c(sigA,sigB))), !any(c(sigA,sigB) %in% node_names))

zsc <- function(X){ (X - rowMeans(X))/rowSds(X) }
cafA <- colMeans(zsc(G[sigA,,drop=FALSE]))
cafB <- colMeans(zsc(G[sigB,,drop=FALSE]))
cat("cor(CAF-A, CAF-B) Spearman =", round(cor(cafA,cafB,method="spearman"),3), "\n")
cat("CAF-A vs COL1A1 rho =", round(cor(cafA, G["COL1A1",], method="spearman"),3),
    "; CAF-B vs COL1A1 rho =", round(cor(cafB, G["COL1A1",], method="spearman"),3), "\n")
cat("CAF-A vs COL3A1 rho =", round(cor(cafA, G["COL3A1",], method="spearman"),3),
    "; CAF-B vs COL3A1 rho =", round(cor(cafB, G["COL3A1",], method="spearman"),3), "\n")
zc <- function(v){ r <- rank(v); (r-mean(r))/sd(r) }
partA <- make_partial(zc(cafA)); partB <- make_partial(zc(cafB))

## ------------------------------------------------------- compute rho ------
## normalise assay labels: source_type/target_type are TF/Gene/miRNA
Eok$sty2 <- ifelse(Eok$source_type=="miRNA","miRNA","gene")
Eok$tty2 <- ifelse(Eok$target_type=="miRNA","miRNA","gene")
stopifnot(all(Eok$sty2 %in% c("miRNA","gene")), all(Eok$tty2 %in% c("miRNA","gene")))
Eok$rho      <- rho_pairs(Eok$si, Eok$sty2, Eok$ti, Eok$tty2)
stopifnot(sum(Eok$rho==0) < 5)   # guard against the type-label bug
Eok$rho_cafA <- partA(Eok$rho, Eok$si, Eok$sty2, Eok$ti, Eok$tty2)
Eok$rho_cafB <- partB(Eok$rho, Eok$si, Eok$sty2, Eok$ti, Eok$tty2)

nrho      <- rho_pairs(nsi, nsty, nti, ntty)
nrho_cafA <- partA(nrho, nsi, nsty, nti, ntty)
nrho_cafB <- partB(nrho, nsi, nsty, nti, ntty)

## expected sign per edge class
expsign <- function(et, mode){
  s <- rep(NA_real_, length(et))
  s[et=="miRNA_target"] <- -1
  s[et %in% c("TF_target","TF_miRNA") & mode=="Activation"] <- +1
  s[et %in% c("TF_target","TF_miRNA") & mode=="Repression"] <- -1
  s
}
Eok$exp <- expsign(Eok$edge_type, Eok$mode)
nexp <- Eok$exp[nidx]

## ------------------------------------------------------ summarise ---------
summ <- function(label, sel){
  r  <- Eok$rho[sel];  ra <- Eok$rho_cafA[sel]; rb <- Eok$rho_cafB[sel]
  ex <- Eok$exp[sel]
  ns <- nidx %in% which(sel)
  nr <- nrho[ns]; nra <- nrho_cafA[ns]; nrb <- nrho_cafB[ns]; nex <- nexp[ns]
  cc <- function(x,e) mean(sign(x)==e, na.rm=TRUE)
  pt <- function(x1,n1,x2,n2) prop.test(c(x1,x2), c(n1,n2))$p.value
  o <- list(set=label, n_edge=sum(sel), n_null=sum(ns),
            conc=cc(r,ex), conc_null=cc(nr,nex),
            conc_p = pt(sum(sign(r)==ex,na.rm=TRUE), sum(!is.na(r)),
                        sum(sign(nr)==nex,na.rm=TRUE), sum(!is.na(nr))),
            conc_binom_p = binom.test(sum(sign(r)==ex,na.rm=TRUE), sum(!is.na(r)),
                                      p=cc(nr,nex))$p.value,
            mean_rho=mean(r), mean_rho_null=mean(nr),
            rho_wilcox_p = wilcox.test(r, nr)$p.value,
            rho_t_p = t.test(r, nr)$p.value,
            conc_cafA=cc(ra,ex), conc_null_cafA=cc(nra,nex),
            conc_p_cafA = pt(sum(sign(ra)==ex,na.rm=TRUE), sum(!is.na(ra)),
                             sum(sign(nra)==nex,na.rm=TRUE), sum(!is.na(nra))),
            mean_rho_cafA=mean(ra), mean_rho_null_cafA=mean(nra),
            rho_wilcox_p_cafA = wilcox.test(ra, nra)$p.value,
            conc_cafB=cc(rb,ex), conc_null_cafB=cc(nrb,nex),
            conc_p_cafB = pt(sum(sign(rb)==ex,na.rm=TRUE), sum(!is.na(rb)),
                             sum(sign(nrb)==nex,na.rm=TRUE), sum(!is.na(nrb))),
            mean_rho_cafB=mean(rb), mean_rho_null_cafB=mean(nrb),
            rho_wilcox_p_cafB = wilcox.test(rb, nrb)$p.value)
  as.data.frame(o, stringsAsFactors=FALSE)
}

sets <- list(
  "TF_target_Activation" = Eok$edge_type=="TF_target" & Eok$mode=="Activation",
  "TF_target_Repression" = Eok$edge_type=="TF_target" & Eok$mode=="Repression",
  "TF_target_Unknown"    = Eok$edge_type=="TF_target" & Eok$mode=="Unknown",
  "TF_target_Ambiguous"  = Eok$edge_type=="TF_target" & Eok$mode=="Ambiguous",
  "TF_miRNA_Activation"  = Eok$edge_type=="TF_miRNA"  & Eok$mode=="Activation",
  "TF_miRNA_Repression"  = Eok$edge_type=="TF_miRNA"  & Eok$mode=="Repression",
  "miRNA_target_all"     = Eok$edge_type=="miRNA_target",
  "miRNA_target_strong"  = Eok$edge_type=="miRNA_target" & Eok$tier=="strong",
  "miRNA_target_weak"    = Eok$edge_type=="miRNA_target" & Eok$tier=="weak",
  "miRNA_target_predicted_only" = Eok$edge_type=="miRNA_target" & Eok$tier=="predicted_only")
sets <- sets[sapply(sets, function(s) sum(s & !is.na(Eok$exp))>0)]
res <- do.call(rbind, lapply(names(sets), function(nm) summ(nm, sets[[nm]] & !is.na(Eok$exp))))
print(res[,c("set","n_edge","conc","conc_null","conc_p","conc_binom_p","mean_rho","mean_rho_null","rho_wilcox_p")], digits=4)
cat("\n--- after CAF-A (ESTIMATE stromal, collagens+nodes removed) ---\n")
print(res[,c("set","conc_cafA","conc_null_cafA","conc_p_cafA","mean_rho_cafA","mean_rho_null_cafA")], digits=4)
cat("\n--- after CAF-B (DCN/LUM/FAP/THY1 marker panel) ---\n")
print(res[,c("set","conc_cafB","conc_null_cafB","conc_p_cafB","mean_rho_cafB","mean_rho_null_cafB")], digits=4)

write.csv(res, file.path(OUT,"v2_edge_class_summary_genomewide_null.csv"), row.names=FALSE)
write.csv(Eok[,c("source","target","edge_type","mode","tier","exp","rho","rho_cafA","rho_cafB")],
          file.path(OUT,"v2_edge_rho.csv"), row.names=FALSE)

## ------------------------------------------------- named axes -------------
axes <- rbind(
  data.frame(x=c("ETS1","NFKB1","SP1","RELA"), y="COL1A1", xt="gene", stringsAsFactors=FALSE),
  data.frame(x=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), y="COL1A1", xt="miRNA"),
  data.frame(x=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), y="COL3A1", xt="miRNA"),
  data.frame(x="hsa-miR-101", y="EZH2", xt="miRNA"),
  data.frame(x=c("hsa-let-7b","hsa-let-7e"), y="COL3A1", xt="miRNA"),
  data.frame(x="hsa-miR-130a", y="VEGFA", xt="miRNA"),
  data.frame(x="COL1A1", y="COL3A1", xt="gene"))
axes$in_network <- paste(axes$x, axes$y) %in% paste(E$source, E$target)
axes$si <- ifelse(axes$xt=="miRNA", match(axes$x, rownames(M)), match(axes$x, rownames(G)))
axes$ti <- match(axes$y, rownames(G))
ok <- !is.na(axes$si) & !is.na(axes$ti)
axes$rho <- NA; axes$rho_cafA <- NA; axes$rho_cafB <- NA
axes$rho[ok] <- rho_pairs(axes$si[ok], axes$xt[ok], axes$ti[ok], rep("gene",sum(ok)))
axes$rho_cafA[ok] <- partA(axes$rho[ok], axes$si[ok], axes$xt[ok], axes$ti[ok], rep("gene",sum(ok)))
axes$rho_cafB[ok] <- partB(axes$rho[ok], axes$si[ok], axes$xt[ok], axes$ti[ok], rep("gene",sum(ok)))
axes$pct_retained_cafA <- round(100*axes$rho_cafA/axes$rho,1)

## ------------------------------------------- CAF tertile stratification ---
ter <- cut(cafA, breaks=quantile(cafA, c(0,1/3,2/3,1)), include.lowest=TRUE,
           labels=c("CAFlow","CAFmid","CAFhigh"))
cat("\nCAF-A tertile sizes:", paste(table(ter), collapse="/"), "\n")
for (lv in levels(ter)){
  k <- which(ter==lv)
  Gs <- G[,k,drop=FALSE]; Ms <- M[,k,drop=FALSE]
  v <- rep(NA_real_, nrow(axes))
  for (r in which(ok)){
    xv <- if (axes$xt[r]=="miRNA") Ms[axes$x[r],] else Gs[axes$x[r],]
    v[r] <- suppressWarnings(cor(xv, Gs[axes$y[r],], method="spearman"))
  }
  axes[[paste0("rho_",lv)]] <- v
}
print(axes[,c("x","y","in_network","rho","rho_cafA","rho_cafB","pct_retained_cafA",
              "rho_CAFlow","rho_CAFmid","rho_CAFhigh")], digits=3)
write.csv(axes, file.path(OUT,"v2_named_axes.csv"), row.names=FALSE)

## network-wide concordance within tertiles for TF_target Activation & miRNA_target
cat("\n--- tertile-stratified network-wide concordance ---\n")
tstr <- NULL
for (lv in levels(ter)){
  k <- which(ter==lv)
  ZGs <- zrows(t(apply(G[,k,drop=FALSE],1,rank))); ZMs <- zrows(t(apply(M[,k,drop=FALSE],1,rank)))
  dn <- length(k)-1
  rr <- numeric(nrow(Eok))
  for (a in c("miRNA","gene")) for (b in c("miRNA","gene")){
    kk <- which(ifelse(Eok$source_type=="miRNA","miRNA","gene")==a &
                ifelse(Eok$target_type=="miRNA","miRNA","gene")==b)
    if(!length(kk)) next
    Za <- if(a=="miRNA") ZMs else ZGs; Zb <- if(b=="miRNA") ZMs else ZGs
    rr[kk] <- rowSums(Za[Eok$si[kk],,drop=FALSE]*Zb[Eok$ti[kk],,drop=FALSE])/dn
  }
  for (nm in names(sets)){
    s <- sets[[nm]] & !is.na(Eok$exp)
    tstr <- rbind(tstr, data.frame(tertile=lv, set=nm, n=sum(s),
      conc=mean(sign(rr[s])==Eok$exp[s], na.rm=TRUE), mean_rho=mean(rr[s], na.rm=TRUE)))
  }
}
print(reshape(tstr[,c("tertile","set","conc")], idvar="set", timevar="tertile", direction="wide"), digits=4)
write.csv(tstr, file.path(OUT,"v2_tertile_concordance.csv"), row.names=FALSE)

saveRDS(list(cafA=cafA, cafB=cafB, samp=samp, sigA=sigA, sigB=sigB),
        file.path(OUT,"v2_caf_scores.rds"))
cat("\nDONE 02\n")
