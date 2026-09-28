#!/usr/bin/env Rscript
## 19) Build the miR-29 family target set (mirrors 01_mir130a_target_set.R) and its
##     TCGA-BRCA correlations, so that the miR-29 axis can be tested in the same way
##     as miR-130a in the functional-genomic screens and single-cell analyses.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260909)
REV <- "/path/to/revision"
OUT <- file.path(REV,"results/v3"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SCR <- "/path/to/scratch"
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

## ---------- 1. miRTarBase ----------
mtb <- fread(file.path(SCR,"mtb_29.csv"))
setnames(mtb, c("miRTarBase ID","miRNA","Species (miRNA)","Target Gene",
                "Target Gene (Entrez ID)","Species (Target Gene)","Experiments",
                "Support Type","References (PMID)"),
              c("MTI","miRNA","sp_mir","gene","entrez","sp_gene","exp","support","pmid"),
         skip_absent=TRUE)
mtb <- mtb[sp_mir=="hsa" & sp_gene=="hsa"]
msg("miRTarBase human-human miR-29 records:", nrow(mtb), "| unique genes:", uniqueN(mtb$gene))
print(mtb[, .N, by=.(miRNA)][order(-N)])

strong_pat <- "[Ll]uciferase|[Rr]eporter ?assay|[Ww]estern ?blot|qRT-PCR|qPCR|Real-time|Northern|ELISA|Immunoblot|Flow cytometry|Immunohisto|Immunocyto"
weak_pat   <- "CLIP|CLASH|Degradome|Microarray|Sequencing|pSILAC|Proteomics|Chimeric|iCLIP|PAR-CLIP|HITS-CLIP"
mtb[, is_strong := grepl(strong_pat, exp, ignore.case=TRUE)]
mtb[, is_clip   := grepl(weak_pat,   exp, ignore.case=TRUE)]
msg("records with a low-throughput assay:", sum(mtb$is_strong),
    "| high-throughput-only:", sum(!mtb$is_strong & mtb$is_clip))
mtb_g <- mtb[, .(mtb_n_records=.N, mtb_strong=any(is_strong), mtb_clip=any(is_clip),
                 mtb_functional=any(support=="Functional MTI"),
                 mtb_arms=paste(sort(unique(miRNA)), collapse=";"),
                 mtb_pmids=paste(sort(unique(sub("\\.0$","",as.character(pmid)))), collapse=";")),
             by=gene]
msg("miRTarBase unique miR-29 target genes:", nrow(mtb_g),
    "| with low-throughput assay:", sum(mtb_g$mtb_strong))

## ---------- 2. TargetScan 8 (conserved / default predictions, miR-29-3p family AGCACCA) ----------
ts <- fread(file.path(SCR,"ts_ctx_29.txt"))
setnames(ts, c("Gene Symbol","Site Type","context++ score","context++ score percentile",
               "weighted context++ score","weighted context++ score percentile","miRNA"),
             c("gene","site_type","ctx","ctx_pct","wctx","wctx_pct","miRNA"), skip_absent=TRUE)
sitelab <- c("1"="7mer-1a","2"="7mer-m8","3"="8mer-1a")
ts[, site_lab := sitelab[as.character(site_type)]]
ts1 <- ts[miRNA=="hsa-miR-29a-3p"]          # the three arms carry identical seeds -> identical sites
msg("TargetScan8 conserved sites for hsa-miR-29a-3p:", nrow(ts1), "| genes:", uniqueN(ts1$gene))
print(ts1[, .N, by=site_lab])
tsum <- fread(file.path(SCR,"ts_sum_29.txt"))
setnames(tsum, c("Gene Symbol","miRNA family","Aggregate PCT"), c("gene","fam","aggPCT"), skip_absent=TRUE)
tsum[, aggPCT := suppressWarnings(as.numeric(aggPCT))]
ts_g <- ts1[, .(ts_n_sites=.N, ts_best_site=site_lab[which.min(ctx)],
                ts_site_types=paste(sort(unique(site_lab)), collapse=";"),
                ts_ctx_sum=sum(ctx, na.rm=TRUE), ts_wctx_min=min(wctx, na.rm=TRUE)), by=gene]
ts_g <- merge(ts_g, tsum[, .(ts_aggPCT=max(aggPCT, na.rm=TRUE)), by=gene], by="gene", all.x=TRUE)
msg("TargetScan8 conserved predicted genes:", nrow(ts_g),
    "| aggregate PCT > 0.5:", sum(ts_g$ts_aggPCT > 0.5, na.rm=TRUE))

## ---------- 3. merge + tier ----------
A <- data.table(gene=sort(unique(c(mtb_g$gene, ts_g$gene))))
A <- merge(A, mtb_g, by="gene", all.x=TRUE); A <- merge(A, ts_g, by="gene", all.x=TRUE)
for (cc in c("mtb_strong","mtb_clip","mtb_functional")) A[[cc]][is.na(A[[cc]])] <- FALSE
A[, in_targetscan := !is.na(ts_n_sites)]
A[, strong := mtb_strong]
A[, clip := mtb_clip]
A[, tier := fifelse(strong, "STRONG_lowthroughput",
             fifelse(clip | !is.na(mtb_n_records), "WEAK_highthroughput", "PREDICTED_only"))]
print(A[, .N, by=tier])
fwrite(A, file.path(OUT,"screens_mir29_target_set.csv"))
msg("wrote screens_mir29_target_set.csv:", nrow(A), "genes")

## ---------- 4. TCGA-BRCA correlations vs the miR-29 family aggregate ----------
G  <- readRDS(file.path(REV,"data/brca_gene_expr.rds"))
Mi <- readRDS(file.path(REV,"data/brca_mirna_expr_canonical.rds"))
ph <- readRDS(file.path(REV,"data/brca_pheno.rds"))
stopifnot(!any(duplicated(rownames(G))))
tum <- ph$sample[ph$sample_type=="Primary Tumor"]
s <- intersect(intersect(colnames(G), colnames(Mi)), tum)
msg("paired primary tumours:", length(s))
mirs <- c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"); stopifnot(all(mirs %in% rownames(Mi)))
fam <- log2(colSums(2^Mi[mirs, s, drop=FALSE]))
X <- G[, s, drop=FALSE]
keep <- rowMeans(X > 0) >= 0.25            # expressed, same criterion as the DE background
X <- X[keep, , drop=FALSE]
msg("genes tested (expressed in >=25% of tumours):", nrow(X))
## RULE 4 ASSERTION: feature names must be symbols, not integer indices
nm <- rownames(X); stopifnot(!any(duplicated(nm)))
msg("ASSERT symbols: fraction of rownames starting with a letter =",
    round(mean(grepl("^[A-Za-z]", nm)),4), "| all-digit rownames:", sum(grepl("^[0-9]+$", nm)))
stopifnot(mean(grepl("^[A-Za-z]", nm)) > 0.95)

rx <- t(apply(X, 1, rank)); ry <- rank(fam)
rxc <- rx - rowMeans(rx); ryc <- ry - mean(ry)
rho <- as.numeric((rxc %*% ryc) / sqrt(rowSums(rxc^2) * sum(ryc^2)))
n <- length(s); tstat <- rho*sqrt((n-2)/(1-rho^2)); p <- 2*pt(-abs(tstat), n-2)
R <- data.table(gene=nm, rho_mir29=rho, p=p, mean_expr=rowMeans(X))
## spot-check against cor.test
g0 <- "COL1A1"; ct <- suppressWarnings(cor.test(as.numeric(X[g0,]), fam, method="spearman"))
msg("spot-check", g0, ": vectorised rho =", signif(R[gene==g0, rho_mir29],9),
    "| cor.test rho =", signif(unname(ct$estimate),9))
stopifnot(abs(R[gene==g0, rho_mir29] - unname(ct$estimate)) < 1e-6)

R[, in_target_set := gene %in% A$gene]
R[, tier := A$tier[match(gene, A$gene)]]
R[, fdr_all := p.adjust(p, "BH")]
R[in_target_set==TRUE, fdr_targets := p.adjust(p, "BH")]
R[tier=="STRONG_lowthroughput", fdr_strong := p.adjust(p, "BH")]
setorder(R, rho_mir29)
fwrite(R, file.path(OUT,"screens_mir29_tcga_correlations.csv"))
msg("wrote screens_mir29_tcga_correlations.csv")

anti <- R[in_target_set==TRUE & !is.na(fdr_targets) & fdr_targets<0.05 & rho_mir29<0, gene]
strong_all <- A[tier=="STRONG_lowthroughput", gene]
msg("miR-29 STRONG (low-throughput) targets:", length(strong_all),
    "| of which tested in TCGA:", sum(strong_all %in% R$gene))
msg("miR-29 anti-correlated validated/predicted targets (FDR<0.05, rho<0):", length(anti))
fwrite(data.table(set="mir29_anticorrelated", gene=anti),
       file.path(OUT,"screens_mir29_anticorrelated_genes.csv"))
msg("COL1A1 rho:", signif(R[gene=="COL1A1", rho_mir29],4), " COL3A1 rho:", signif(R[gene=="COL3A1", rho_mir29],4))
msg("DONE 19")
