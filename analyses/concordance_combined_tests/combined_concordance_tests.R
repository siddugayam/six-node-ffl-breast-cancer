## combined_concordance_tests.R -- sign concordance of network edges by evidence tier against the five-draw
## within-network matched null, combined over the five draws in several ways (arithmetic mean of the two-sided p values;
## one test against the pooled null pairs of all five draws; Cauchy and harmonic-mean combinations; one-sided versions),
## unadjusted and adjusted for CAF-A, CAF-B, ABSOLUTE purity, the ESTIMATE stromal score and each of the 42 distinct
## stromal estimates of the mediation analysis, with all edges and with the Farmer-correlated targets left out
## (see SETTINGS.md).
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

## ---- adjustments ----------------------------------------------------------------------------------------------
FIB <- read.csv("results/systems_and_dynamics/deconv_fibroblast_estimates_extended.csv", check.names=FALSE)
FIB <- FIB[match(samp, FIB$sample), ]; stopifnot(identical(FIB$sample, samp))
DUP <- c(MCPcounter_Fibroblasts_recomputed="MCPcounter_Fibroblasts", NNLS_Wu641_full_CAFs="NNLS_Wu_CAFs",
         QPROG_Wu641_full_CAFs="QPROG_Wu_CAFs")
dupcheck <- sapply(names(DUP), function(d) { same <- names(FIB)[sapply(names(FIB), function(c) c != d && c != "sample" &&
  isTRUE(all.equal(FIB[[c]], FIB[[d]], tolerance=0)))]; paste(same, collapse=";") })
cat("exact repeats of the dropped columns:\n"); print(dupcheck); stopifnot(all(nzchar(dupcheck)))
EST <- setdiff(setdiff(colnames(FIB), "sample"), names(DUP)); stopifnot(length(EST) == 42)
farmer <- read.csv("results/prioritisation_and_stroma/farmer_pergene_rho.csv", stringsAsFactors=FALSE)
FARMER_SET <- farmer$gene[farmer$rho_vs_farmer > 0.5]
keep_lo <- !(E$target %in% FARMER_SET)
rowsG <- sort(unique(c(ng, E$si[E$sty=="gene"], E$ti[E$tty=="gene"]))); rowsM <- sort(unique(c(nmi, E$si[E$sty=="miRNA"], E$ti[E$tty=="miRNA"])))
ADJ <- c(list(unadjusted=NULL, CAF_A=cs$cafA, CAF_B=cs$cafB), setNames(lapply(EST, function(e) FIB[[e]]), paste0("estimate:", EST)))

cauchy <- function(p){ p <- pmin(pmax(p, 1e-300), 1 - 1e-16); T <- mean(tan((0.5 - p) * pi)); 0.5 - atan(T) / pi }
hmp_vw <- function(p) min(1, exp(1) * log(length(p)) / mean(1 / p))
res <- list()
for(an in names(ADJ)){
  v <- ADJ[[an]]
  if(is.null(v)) { ok <- rep(TRUE, N) } else ok <- is.finite(v)
  ns <- sum(ok)
  if(ns == N){ ZGs <- ZG; ZMs <- ZM; dns <- dn } else {
    ZGs <- matrix(NA_real_, nrow(G), ns); ZMs <- matrix(NA_real_, nrow(M), ns)
    ZGs[rowsG,] <- zr(G[rowsG, ok, drop=FALSE]); ZMs[rowsM,] <- zr(M[rowsM, ok, drop=FALSE]); dns <- ns-1 }
  rhoS <- function(si,sty,ti,tty){o <- numeric(length(si))
    for(a in c("miRNA","gene")) for(b in c("miRNA","gene")){k <- which(sty==a&tty==b); if(!length(k))next
      Za <- if(a=="miRNA")ZMs else ZGs; Zb <- if(b=="miRNA")ZMs else ZGs
      o[k] <- rowSums(Za[si[k],,drop=FALSE]*Zb[ti[k],,drop=FALSE])/dns}; o}
  if(is.null(v)) padj <- function(r, ...) r else {
    zv <- zc(v[ok]); CG <- numeric(nrow(G)); CM <- numeric(nrow(M))
    CG[rowsG] <- as.vector(ZGs[rowsG,,drop=FALSE] %*% zv)/dns; CM[rowsM] <- as.vector(ZMs[rowsM,,drop=FALSE] %*% zv)/dns
    padj <- function(r,si,sty,ti,tty){ a <- ifelse(sty=="miRNA",CM[si],CG[si]); b <- ifelse(tty=="miRNA",CM[ti],CG[ti])
      (r-a*b)/sqrt((1-a^2)*(1-b^2)) } }
  ea <- padj(rhoS(E$si,E$sty,E$ti,E$tty), E$si,E$sty,E$ti,E$tty)
  NV <- lapply(DRAWS, function(D) list(na=padj(rhoS(D$si,D$sty,D$ti,D$tty),D$si,D$sty,D$ti,D$tty), idx=D$idx))
  for(variant in c("all edges","Farmer-correlated targets left out")){
    kp <- if(variant=="all edges") rep(TRUE, nrow(E)) else keep_lo
    for(nmm in names(sets)){ s <- which(sets[[nmm]] & !is.na(E$exp) & kp); if(!length(s)) next
      x1 <- sum(sign(ea[s])==E$exp[s]); n <- length(s)
      p2 <- p1 <- p0 <- numeric(REPS); conc_null_pairs <- 0; npairs <- 0
      for(d in seq_len(REPS)){ k <- NV[[d]]$idx %in% s; nv <- NV[[d]]$na[k]; nex <- E$exp[NV[[d]]$idx[k]]
        hit <- sum(sign(nv)==nex); p0[d] <- hit/length(nv); conc_null_pairs <- conc_null_pairs + hit; npairs <- npairs + length(nv)
        p2[d] <- binom.test(x1, n, p=p0[d])$p.value; p1[d] <- binom.test(x1, n, p=p0[d], alternative="greater")$p.value }
      pp <- conc_null_pairs/npairs
      res[[length(res)+1]] <- data.frame(adjustment=an, n_tumours=ns, edges=variant, set=nmm, n_edges=n,
        pairs_per_edge_per_draw=npairs/(REPS*n), conc=x1/n, null_mean_of_draws=mean(p0), null_pooled=pp,
        p_mean_two_sided=mean(p2), p_twice_mean_two_sided=min(1, 2*mean(p2)),
        p_pooled_two_sided=binom.test(x1, n, p=pp)$p.value, p_cauchy_two_sided=cauchy(p2), p_hmp_vovk_wang_two_sided=hmp_vw(p2),
        p_pooled_one_sided=binom.test(x1, n, p=pp, alternative="greater")$p.value, p_mean_one_sided=mean(p1),
        p_cauchy_one_sided=cauchy(p1), p_draw1=p2[1], p_draw2=p2[2], p_draw3=p2[3], p_draw4=p2[4], p_draw5=p2[5])
    } }
  cat(sprintf("%-40s tumours %d done\n", an, ns))
}
R <- do.call(rbind, res)
## the unadjusted and CAF-A/CAF-B rows must equal 04 (all edges)
chk <- merge(R[R$edges=="all edges" & R$adjustment %in% c("unadjusted","CAF_A","CAF_B"),], agg, by="set")
for(a in c("unadjusted","CAF_A","CAF_B")){ tag <- c(unadjusted="raw", CAF_A="cafA", CAF_B="cafB")[[a]]; z <- chk[chk$adjustment==a,]
  stopifnot(max(abs(z$conc - z[[paste0("conc_",tag)]])) < 1e-12, max(abs(z$p_mean_two_sided - z[[paste0("p_",tag)]])) < 1e-12) }
write.csv(R, file.path(HERE, "combined_tests_by_adjustment.csv"), row.names=FALSE)

METHODS <- c(mean_two_sided="p_mean_two_sided", twice_mean="p_twice_mean_two_sided", pooled_two_sided="p_pooled_two_sided",
             cauchy_two_sided="p_cauchy_two_sided", harmonic_mean_vovk_wang="p_hmp_vovk_wang_two_sided",
             pooled_one_sided="p_pooled_one_sided", mean_one_sided="p_mean_one_sided", cauchy_one_sided="p_cauchy_one_sided")
## summary over the 42 distinct estimates, in the layout of `concordance_summary_42_distinct.csv` plus the combined-test columns
S <- list()
for(variant in c("all edges","Farmer-correlated targets left out")) for(nmm in sort(names(sets))){
  a <- R[R$edges==variant & R$set==nmm & grepl("^estimate:", R$adjustment), ]; if(!nrow(a)) next
  stopifnot(nrow(a) == 42)
  row <- data.frame(edges=variant, set=nmm, n_estimates=nrow(a), p_median=median(a$p_mean_two_sided), p_min=min(a$p_mean_two_sided),
    p_max=max(a$p_mean_two_sided), n_p_below_0.05=sum(a$p_mean_two_sided<0.05), n_above_null=sum(a$conc>a$null_mean_of_draws),
    excess_median=median(a$conc-a$null_mean_of_draws), excess_min=min(a$conc-a$null_mean_of_draws), excess_max=max(a$conc-a$null_mean_of_draws),
    check.names=FALSE)
  for(m in names(METHODS)[-1]){ v <- a[[METHODS[[m]]]]
    row[[paste0(m,"_p_median")]] <- median(v); row[[paste0(m,"_p_min")]] <- min(v); row[[paste0(m,"_p_max")]] <- max(v)
    row[[paste0(m,"_n_below_0.05")]] <- sum(v<0.05) }
  S[[length(S)+1]] <- row }
S <- do.call(rbind, S)
write.csv(S, file.path(HERE, "combined_tests_summary_42_distinct.csv"), row.names=FALSE)
## conclusion changes relative to the arithmetic mean
C <- R[, c("adjustment","edges","set","n_edges","conc","null_mean_of_draws")]
for(m in names(METHODS)) C[[paste0("sig_",m)]] <- R[[METHODS[[m]]]] < 0.05
C$changes_vs_mean <- apply(C[, paste0("sig_", names(METHODS)[-1])] != C$sig_mean_two_sided, 1, function(z) paste(names(METHODS)[-1][z], collapse=";"))
write.csv(C, file.path(HERE, "conclusion_changes.csv"), row.names=FALSE)
cat("adjustment x set rows with any change vs the arithmetic mean:", sum(nzchar(C$changes_vs_mean)), "of", nrow(C), "\n")
print(R[R$adjustment %in% c("unadjusted","CAF_A","CAF_B") & R$edges=="all edges" & R$set %in% c("TF_target_Activation","miRNA_target_all","miRNA_target_strong","miRNA_target_weak","miRNA_target_predicted_only"),
        c("adjustment","set","n_edges","conc","null_pooled","p_mean_two_sided","p_pooled_two_sided","p_cauchy_two_sided","p_hmp_vovk_wang_two_sided","p_pooled_one_sided")], digits=3)
cat("DONE\n")
