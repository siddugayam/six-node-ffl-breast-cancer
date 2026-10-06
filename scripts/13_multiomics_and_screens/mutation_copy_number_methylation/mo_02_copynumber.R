#!/usr/bin/env Rscript
## Part C: GISTIC2 thresholded copy number for network nodes + CN-expression correlation
suppressPackageStartupMessages({library(data.table)})
RV  <- "/path/to/revision"
CA  <- file.path(RV,"cache/multiomics"); OUT <- file.path(RV,"results/multiomics")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
nodes <- fread(file.path(RV,"data/canonical_nodes.tsv"))

## ---------- miRNA precursor coordinates (hg19, UCSC wgRna = miRBase) ----------
wg <- fread(cmd=paste0("zcat ",CA,"/wgRna_hg19.txt.gz"), header=FALSE)
setnames(wg, c("bin","chrom","start","end","name","score","strand","thickStart","thickEnd","type"))
wg <- wg[type=="miRNA"]
mk <- function(x){ y <- sub("^hsa-","",x); y <- sub("^mir-","MIR",y); y <- sub("^let-","MIRLET",y); toupper(gsub("-","",y)) }
wg[, sym := mk(name)]
cat("wgRna hg19 miRNA precursors:", nrow(wg), "\n")
mat2pre <- function(mature){
  stem <- sub("-(3p|5p)$","", sub("^hsa-","", tolower(mature)))
  wg$name[grepl(paste0("^hsa-", stem, "(-[0-9]+)?$"), wg$name)]
}

## ---------- GISTIC2 ----------
g <- fread(cmd=paste0("zcat ",CA,"/Gistic2_CN_thresholded.gz")); setnames(g,1,"gene")
gm <- as.matrix(g[,-1]); rownames(gm) <- g$gene; storage.mode(gm) <- "numeric"
cat("GISTIC2:", nrow(gm), "genes x", ncol(gm), "samples; NA cells:", sum(is.na(gm)), "\n")

## ---------- refGene hg19 (transcripts kept separate; merge only within 1 Mb) ----------
rg <- fread(cmd=paste0("zcat ",CA,"/refGene_hg19.txt.gz"), header=FALSE, select=c(3,5,6,13))
setnames(rg, c("chrom","txStart","txEnd","sym"))
rg <- rg[chrom %in% paste0("chr",c(1:22,"X","Y"))]
setorder(rg, sym, chrom, txStart)
rg[, grp := cumsum(c(TRUE, !(sym[-1]==sym[-.N] & chrom[-1]==chrom[-.N] & txStart[-1]-txEnd[-.N] < 1e6)))]
rgg <- rg[, .(sym=sym[1], chrom=chrom[1], start=min(txStart), end=max(txEnd)), by=grp][,-1]
rggG <- rgg[sym %in% rownames(gm)]
cat("refGene loci:", nrow(rgg), "| with GISTIC data:", nrow(rggG), "\n")

proxy_for <- function(ch, st, en, win=1e6, k=5){
  cand <- rggG[chrom==ch & end > (st-win) & start < (en+win)]
  if(!nrow(cand)) return(list(v=NULL, used=NA_character_))
  cand[, d := pmax(0, pmax(start-en, st-end))]
  cand <- unique(cand, by="sym")[order(d)][1:min(k,.N)]
  list(v=colMeans(gm[cand$sym,,drop=FALSE], na.rm=TRUE), used=paste(cand$sym, collapse="|"))
}

ge <- readRDS(file.path(RV,"data/brca_gene_expr.rds"))
me <- readRDS(file.path(RV,"data/brca_mirna_expr_canonical.rds"))

## ---------- enumerate loci ----------
L <- list()
for(i in seq_len(nrow(nodes))){
  nm <- nodes$name[i]; ty <- nodes$type[i]
  if(ty %in% c("Gene","TF")){
    rr <- rggG[sym==nm]; if(!nrow(rr)) rr <- rgg[sym==nm]
    rr <- rr[order(-(end-start))][1]
    ch <- if(nrow(rr) && !is.na(rr$chrom)) rr$chrom else NA_character_
    if(nm %in% rownames(gm)){
      L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=nm,chrom=ch,
        lstart=if(nrow(rr)) rr$start else NA_integer_, lend=if(nrow(rr)) rr$end else NA_integer_,
        cn_source=nm, cn_method="direct_gene")
    } else if(nrow(rr) && !is.na(rr$chrom)){
      p <- proxy_for(rr$chrom, rr$start, rr$end)
      L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=nm,chrom=rr$chrom,lstart=rr$start,lend=rr$end,
        cn_source=p$used, cn_method=if(is.na(p$used)) "unmapped" else "nearest_gene_proxy")
    } else L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=nm,chrom=NA_character_,
        lstart=NA_integer_,lend=NA_integer_,cn_source=NA_character_,cn_method="unmapped")
  } else {
    pres <- mat2pre(nm)
    if(!length(pres)){
      L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=NA_character_,chrom=NA_character_,
        lstart=NA_integer_,lend=NA_integer_,cn_source=NA_character_,cn_method="no_hg19_precursor")
    } else for(pn in pres){
      w <- wg[name==pn][1]
      if(w$sym %in% rownames(gm)){
        L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=pn,chrom=w$chrom,lstart=w$start,lend=w$end,
          cn_source=w$sym, cn_method="direct_mirna_gene")
      } else {
        p <- proxy_for(w$chrom, w$start, w$end)
        L[[length(L)+1]] <- data.table(node=nm,type=ty,locus=pn,chrom=w$chrom,lstart=w$start,lend=w$end,
          cn_source=p$used, cn_method=if(is.na(p$used)) "unmapped" else "nearest_gene_proxy")
      }
    }
  }
}
L <- rbindlist(L)
cat("loci enumerated:", nrow(L), "for", uniqueN(L$node), "nodes\n")

## ---------- per-locus statistics ----------
getv <- function(src){ s <- strsplit(src,"\\|")[[1]]; s <- s[s %in% rownames(gm)]
                       if(!length(s)) return(NULL); colMeans(gm[s,,drop=FALSE], na.rm=TRUE) }
stats <- rbindlist(lapply(seq_len(nrow(L)), function(i){
  r <- L[i]; v <- if(is.na(r$cn_source)) NULL else getv(r$cn_source)
  if(is.null(v)) return(data.table(n_cn_samples=NA_integer_,frac_amp=NA_real_,frac_del=NA_real_,
      frac_highamp=NA_real_,frac_deepdel=NA_real_,mean_cn=NA_real_,
      cn_expr_rho=NA_real_,cn_expr_p=NA_real_,n_cn_expr=NA_integer_))
  ex <- if(r$type %in% c("Gene","TF")){ if(r$node %in% rownames(ge)) ge[r$node,] else NULL
        } else { if(r$node %in% rownames(me)) me[r$node,] else NULL }
  rho<-NA_real_; pv<-NA_real_; nn<-NA_integer_
  if(!is.null(ex)){ cs <- intersect(names(v), names(ex))
    if(length(cs)>=30){ ct <- suppressWarnings(cor.test(v[cs],ex[cs],method="spearman"))
      rho <- unname(ct$estimate); pv <- ct$p.value; nn <- length(cs) } }
  data.table(n_cn_samples=sum(!is.na(v)), frac_amp=mean(v>=1,na.rm=TRUE), frac_del=mean(v<=-1,na.rm=TRUE),
    frac_highamp=mean(v>=2,na.rm=TRUE), frac_deepdel=mean(v<=-2,na.rm=TRUE), mean_cn=mean(v,na.rm=TRUE),
    cn_expr_rho=rho, cn_expr_p=pv, n_cn_expr=nn)
}))
L <- cbind(L, stats)
L[!is.na(cn_expr_p), cn_expr_fdr := p.adjust(cn_expr_p,"BH")]
fwrite(L, file.path(OUT,"brca_copynumber_loci_detail.csv"))
cat("wrote brca_copynumber_loci_detail.csv rows:", nrow(L), "\n")

## ---------- node level: primary locus ----------
L[, n_loci := .N, by=node]
setorder(L, node, -cn_expr_rho, na.last=TRUE)
res <- L[, .SD[1], by=node]
allsum <- L[, .(all_loci=paste(sprintf("%s(%s):amp%.2f/del%.2f", locus, chrom, frac_amp, frac_del), collapse="; ")), by=node]
res <- merge(res, allsum, by="node", sort=FALSE)
setcolorder(res, c("node","type","locus","chrom","lstart","lend","cn_source","cn_method","n_loci"))
fwrite(res[order(-frac_amp)], file.path(OUT,"brca_copynumber_nodes.csv"))
cat("wrote brca_copynumber_nodes.csv rows:", nrow(res), "\n")
print(table(res$cn_method, useNA="always"))

cat("\n=== CN-expression concordance (primary locus) ===\n")
ok <- res[!is.na(cn_expr_rho)]
cat("nodes tested:", nrow(ok), " median rho:", round(median(ok$cn_expr_rho),3),
    " %positive:", round(100*mean(ok$cn_expr_rho>0),1),
    " %FDR<0.05 & rho>0:", round(100*mean(ok$cn_expr_fdr<0.05 & ok$cn_expr_rho>0),1), "\n")
for(t in c("Gene","TF","miRNA")){ s<-ok[type==t]
  cat(sprintf("  %-6s n=%3d median rho=%.3f  %%pos=%.1f  %%sig_pos=%.1f\n",
      t,nrow(s),median(s$cn_expr_rho),100*mean(s$cn_expr_rho>0),100*mean(s$cn_expr_fdr<0.05 & s$cn_expr_rho>0))) }

cat("\n=== FEATURED LOCI (all loci, not just primary) ===\n")
feat <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","COL1A1","COL3A1","EZH2",
          "hsa-miR-101","hsa-miR-124","hsa-let-7b","hsa-let-7e","hsa-miR-34a",
          "hsa-miR-200b","hsa-miR-200c","hsa-miR-145")
print(L[node %in% feat, .(node,locus,chrom,cn_source,cn_method,amp=round(frac_amp,3),
      del=round(frac_del,3),deepdel=round(frac_deepdel,3),rho=round(cn_expr_rho,3),
      fdr=signif(cn_expr_fdr,2))][order(node)], nrows=60)

cat("\n=== top 12 amplified network nodes ===\n")
print(res[order(-frac_amp)][1:12,.(node,type,chrom,amp=round(frac_amp,3),del=round(frac_del,3),rho=round(cn_expr_rho,3))])
cat("=== top 12 deleted network nodes ===\n")
print(res[order(-frac_del)][1:12,.(node,type,chrom,del=round(frac_del,3),amp=round(frac_amp,3),rho=round(cn_expr_rho,3))])

allamp <- rowMeans(gm>=1,na.rm=TRUE); alldel <- rowMeans(gm<=-1,na.rm=TRUE)
netg <- res[type %in% c("Gene","TF") & cn_method=="direct_gene"]
cat("\ngenome-wide mean frac_amp =", round(mean(allamp),4), " frac_del =", round(mean(alldel),4), "\n")
cat("network protein-coding mean frac_amp =", round(mean(netg$frac_amp),4)," frac_del =", round(mean(netg$frac_del),4),"\n")
cat("Wilcoxon amp p =", format.pval(wilcox.test(netg$frac_amp, allamp[!names(allamp)%in%netg$node])$p.value),
    "| del p =", format.pval(wilcox.test(netg$frac_del, alldel[!names(alldel)%in%netg$node])$p.value), "\n")
