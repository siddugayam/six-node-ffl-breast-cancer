#!/usr/bin/env Rscript
## PART 1B) Project Score (Sanger CRISPR-Cas9, Score2 release: Sanger v2 + Broad 21Q2),
## a second independent CRISPR screen with its own cell-line panel and analysis pipeline.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260909)
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3"); CA <- file.path(REV,"cache/v7")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")
ZP <- file.path(CA,"Project_Score2_fitness.zip"); DIR <- file.path(CA,"ps2")
dir.create(DIR, showWarnings=FALSE)
fl <- unzip(ZP, list=TRUE); setDT(fl)
msg("archive contents:"); print(fl[order(-Length)][1:min(25,nrow(fl)), .(Name, MB=round(Length/1e6,1))])
unzip(ZP, exdir=DIR, overwrite=TRUE)
files <- list.files(DIR, recursive=TRUE, full.names=TRUE)
files <- grep("__MACOSX", files, invert=TRUE, value=TRUE)
msg("extracted files:"); print(data.table(f=basename(files), MB=round(file.size(files)/1e6,1))[order(-MB)])

## pick the gene-level fitness/scaled Bayes-factor matrix
## The release ships four gene x model matrices; the scaled Bayesian factor is the one
## Project Score uses for cross-model comparison (negative = depleted / fitness gene).
pick <- grep("scaled_bayesian_factors", files, value=TRUE)[1]
stopifnot(!is.na(pick))
bin  <- grep("binary_matrix", files, value=TRUE)[1]
msg("gene x model matrix used:", basename(pick))
msg("binary significance matrix:", basename(bin))
## header block: row1 model_name, row2 model_id, row3 source, row4 qc_pass, row5 column names
hdr <- fread(pick, nrows=5, header=FALSE, sep="\t", colClasses="character")
model_name <- as.character(hdr[1, ]); model_id <- as.character(hdr[2, ])
src <- as.character(hdr[3, ]); qc <- as.character(hdr[4, ])
M <- fread(pick, skip=5, header=FALSE, sep="\t")
setnames(M, 1:3, c("gene_id","symbol","ensembl_gene_id"))
msg("matrix:", nrow(M), "genes x", ncol(M)-3, "models")
mid <- model_id[4:length(model_id)]; mnm <- model_name[4:length(model_name)]
srcv <- src[4:length(src)]; qcv <- qc[4:length(qc)]
stopifnot(length(mid) == ncol(M)-3)
setnames(M, 4:ncol(M), mid)
msg("sources:", paste(names(table(srcv)), table(srcv), collapse=" ; "),
    "| qc_pass TRUE:", sum(qcv=="TRUE"))
setnames(M, "symbol", "gene")
M <- M[, c("gene", mid), with=FALSE]
if (anyDuplicated(M$gene)) { msg("duplicated symbols:", sum(duplicated(M$gene))); M <- M[!duplicated(gene)] }
stopifnot(mean(grepl("^[A-Za-z]", M$gene)) > 0.9)   # rule 4: rows must be symbols
models <- setdiff(names(M), "gene")
msg("models in matrix:", length(models))

cmp <- fread(file.path(CA,"cmp_model_list.csv"))
bre_ids <- cmp[grepl("Breast", cancer_type, ignore.case=TRUE), model_id]
msg("Cell Model Passports breast models:", length(bre_ids))
bcols <- intersect(models, bre_ids)
## keep only QC-passing models
qcok <- mid[qcv == "TRUE"]
bcols <- intersect(bcols, qcok)
msg("breast models in Project Score that pass QC:", length(bcols))
stopifnot(length(bcols) >= 10)
msg("breast model names:", paste(mnm[match(bcols, mid)], collapse=", "))

X <- as.matrix(M[, ..bcols]); rownames(X) <- M$gene
ocols <- setdiff(models, bcols); Xo <- as.matrix(M[, ..ocols]); rownames(Xo) <- M$gene
sc_b <- rowMeans(X, na.rm=TRUE); sc_o <- rowMeans(Xo, na.rm=TRUE)
msg("value range:", paste(round(range(X, na.rm=TRUE),3), collapse=" .. "),
    "| median:", round(median(X, na.rm=TRUE),4))
msg("NOTE the scaled Bayesian factor is signed the OPPOSITE way to Chronos/DEMETER2:")
msg("  POSITIVE = depleted = fitness/essential gene; NEGATIVE = no fitness effect.")
## controls tell us the scale/direction
for (g in c("RPL7","PSMA1","POLR2A","EIF3B","MYC","PTEN","CDKN2A","KRAS")) if (g %in% rownames(X))
  msg("  control", g, ": mean breast score =", round(sc_b[g],4))
## essentiality threshold: with scaled Bayes factors / Fitness scores the common-essential
## controls sit far from 0; define "significant" empirically as the 5th percentile of the
## genome-wide breast distribution and also report a fixed -0.5 cut for comparability
q05 <- quantile(sc_b, 0.05, na.rm=TRUE)
msg("5th percentile of the genome-wide breast distribution:", round(q05,4))
## Project Score ships its own binary significance matrix (1 = significant fitness gene in
## that model); use it rather than an arbitrary threshold.
bh <- fread(bin, nrows=5, header=FALSE, sep="\t", colClasses="character")
bmid <- as.character(bh[2, ])[-(1:3)]
BM <- fread(bin, skip=5, header=FALSE, sep="\t")
setnames(BM, 1:3, c("gene_id","gene","ensembl_gene_id")); setnames(BM, 4:ncol(BM), bmid)
BM <- BM[!duplicated(gene)]
bcols_bin <- intersect(bcols, names(BM))
msg("binary matrix: genes", nrow(BM), "| breast models present:", length(bcols_bin))
BX <- as.matrix(BM[, ..bcols_bin]); rownames(BX) <- BM$gene
frac_sig_all <- rowMeans(BX == 1, na.rm=TRUE)
frac_sig <- setNames(rep(NA_real_, length(sc_b)), names(sc_b))
frac_sig[intersect(names(sc_b), names(frac_sig_all))] <- frac_sig_all[intersect(names(sc_b), names(frac_sig_all))]
msg("genome-wide: mean fraction of breast models in which a gene is a significant fitness gene:",
    round(mean(frac_sig, na.rm=TRUE), 4))

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
ctrl  <- fread(file.path(OUT,"systems_controllability_nodes.csv"))
hubs  <- ctrl[is_FFL_hub==TRUE, node]
alias <- c(MKL1="MRTFA", CTGF="CCN2", CYR61="CCN1")
look <- function(v, nm){ k <- nm; sw <- !(nm %in% names(v)) & (nm %in% names(alias)); k[sw] <- alias[nm[sw]]
                         o <- v[k]; names(o) <- nm; o }
ND <- data.table(node=nodes$name, type=nodes$type, is_FFL_hub=nodes$name %in% hubs)
ND[, score_breast := look(sc_b, node)][, score_other_lineages := look(sc_o, node)]
ND[, frac_significant := look(frac_sig, node)][, n_lines := length(bcols)]
fwrite(ND, file.path(OUT,"screens_projectscore_network_nodes.csv"))
msg("network nodes with a Project Score value:", sum(!is.na(ND$score_breast)), "of", nrow(ND))

FOC <- ND[node %in% c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC","TP53","STAT3","HIF1A","MKL1","ESR1","XIAP","PTEN")]
msg("FOCUS GENES in Project Score (breast models):"); print(FOC, digits=4)
fwrite(FOC, file.path(OUT,"screens_projectscore_focus_genes.csv"))

## set-level test, same decile-matched null construction as the CRISPR and RNAi analyses
A29  <- fread(file.path(OUT,"screens_mir29_target_set.csv"))
R29  <- fread(file.path(OUT,"screens_mir29_tcga_correlations.csv"))
universe <- R29$gene; gs <- rownames(X)
sets <- list(
  mir29_STRONG = A29[tier=="STRONG_lowthroughput", gene],
  mir29_anticorrelated = fread(file.path(OUT,"screens_mir29_anticorrelated_genes.csv"))$gene,
  FFL_hubs = hubs, network_TFs = nodes[type=="TF", name], network_genes = nodes[type=="Gene", name])
sets <- lapply(sets, function(g) intersect(intersect(unique(g), gs), universe))
mexp <- setNames(R29$mean_expr, R29$gene); pu <- intersect(gs, universe)
brk <- unique(quantile(mexp[pu], probs=seq(0,1,0.1)))
dec <- setNames(as.integer(cut(mexp[pu], breaks=brk, include.lowest=TRUE)), pu)
B <- 2000L
res <- rbindlist(lapply(names(sets), function(nm){
  g <- sets[[nm]]; if (length(g) < 5) return(NULL)
  excl <- if (grepl("^mir29", nm)) A29$gene else g
  pool <- setdiff(pu, excl); pbd <- split(pool, dec[pool])
  obs <- mean(sc_b[g], na.rm=TRUE); obsf <- mean(frac_sig[g] > 0.5, na.rm=TRUE)
  nmx <- matrix(unlist(lapply(dec[g], function(d) sample(pbd[[as.character(d)]], B, replace=TRUE))), nrow=length(g), byrow=TRUE)
  nl <- colMeans(matrix(sc_b[nmx], nrow=length(g)), na.rm=TRUE)
  nf <- colMeans(matrix(frac_sig[nmx] > 0.5, nrow=length(g)), na.rm=TRUE)
  ## higher scaled BF = MORE essential, so "more essential than null" is the upper tail
  data.table(set=nm, n=length(g), obs_mean_scaledBF_breast=obs, null_mean=mean(nl), null_sd=sd(nl),
             emp_p_more_essential=(1+sum(nl>=obs))/(1+B),
             emp_p_less_essential=(1+sum(nl<=obs))/(1+B),
             obs_frac_fitness_gene_in_most_models=obsf, null_frac=mean(nf),
             emp_p_frac=(1+sum(nf>=obsf))/(1+B))
}))
msg("Project Score set-level test vs expression-decile-matched null (B =", B, "):"); print(res, digits=3)
fwrite(res, file.path(OUT,"screens_projectscore_setlevel.csv"))
msg("DONE 21")
