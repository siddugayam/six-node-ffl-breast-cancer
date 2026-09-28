#!/usr/bin/env Rscript
# 23_rppa_second_platform.R -- TCGA RPPA as an independent protein platform
suppressMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results","multiomics")
con <- file(file.path(REV,"logs","rppa_second_platform.log"), open="wt")
say <- function(...) { m<-paste0(format(Sys.time(),"%H:%M:%S")," | ",paste0(...,collapse="")); cat(m,"\n"); cat(m,"\n",file=con); flush(con) }
set.seed(1234); say("=== 23 START ===")
f_rppa <- "/path/to/home/Desktop/DD/R_GPR/ML/TCGA-RPPA-pancan-clean.xena.gz"
f_pheno<- "/path/to/home/Desktop/DD/R_GPR/ML/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
stopifnot(file.exists(f_rppa), file.exists(f_pheno))
RP <- fread(f_rppa, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
rn <- RP[[1]]; RP <- as.matrix(RP[,-1,drop=FALSE]); rownames(RP) <- rn
say("RPPA matrix ", nrow(RP), " antibodies x ", ncol(RP), " samples")
PH <- fread(f_pheno, sep="\t", header=TRUE, data.table=FALSE)
say("phenotype rows ", nrow(PH), " cols: ", paste(colnames(PH), collapse=","))
brca <- PH[PH$`_primary_disease`=="breast invasive carcinoma", ]
say("BRCA samples in phenotype file: ", nrow(brca), " ; types: ", paste(names(table(brca$sample_type)), table(brca$sample_type), sep="=", collapse=" "))
inb <- intersect(colnames(RP), brca$sample)
st  <- setNames(brca$sample_type, brca$sample)[inb]
say("RPPA x BRCA overlap: ", length(inb), " ; by type: ", paste(names(table(st)), table(st), sep="=", collapse=" "))
n_normal <- sum(st=="Solid Tissue Normal")
say("*** BRCA SOLID TISSUE NORMAL samples present in this RPPA file: ", n_normal, " ***")
say("suffix census across the WHOLE RPPA file: ", paste(names(table(sub(".*-","",colnames(RP)))), table(sub(".*-","",colnames(RP))), sep="=", collapse=" "))

tum <- inb[st=="Primary Tumor"]
RPt <- RP[, tum, drop=FALSE]
say("BRCA primary-tumour RPPA matrix: ", nrow(RPt), " x ", ncol(RPt))

## ---- curated antibody -> network-node map ----
ab2gene <- c(
 AR="AR", AXL="AXL", BCLXL="BCL2L1", BETACATENIN="CTNNB1", CAVEOLIN1="CAV1",
 CMET="MET", CMET_pY1235="MET", CMYC="MYC", CJUN_pS73="JUN", E2F1="E2F1",
 ERALPHA="ESR1", ERALPHA_pS118="ESR1", ETS1="ETS1", EZH2="EZH2", FIBRONECTIN="FN1",
 FOXM1="FOXM1", FOXO3A="FOXO3", FOXO3A_pS318S321="FOXO3", GATA3="GATA3", GATA6="GATA6",
 HIF1ALPHA="HIF1A", IRF1="IRF1", LDHA="LDHA", MYH11="MYH11", MYOSINIIA="MYH9",
 MYOSINIIA_pS1943="MYH9", NFKBP65_pS536="RELA", NOTCH1="NOTCH1", P53="TP53",
 P62LCKLIGAND="SQSTM1", P63="TP63", PAI1="SERPINE1", PCADHERIN="CDH3", RB="RB1",
 RB_pS807S811="RB1", SHC_pY317="SHC1", SMAD3="SMAD3", SMAD4="SMAD4", SNAIL="SNAI1",
 STAT3_pY705="STAT3", STAT5ALPHA="STAT5A", SYNAPTOPHYSIN="SYP",
 THYMIDILATESYNTHASE="TYMS", XBP1="XBP1", YB1="YBX1", YB1_pS102="YBX1")
is_phospho <- grepl("_p[STY][0-9]", names(ab2gene)); names(is_phospho) <- names(ab2gene)
nodes <- read.delim(file.path(REV,"data","canonical_nodes.tsv"), stringsAsFactors=FALSE)
hubs  <- read.csv(file.path(REV,"results","ffl_network_hubs.csv"), stringsAsFactors=FALSE)
ab2gene <- ab2gene[names(ab2gene) %in% rownames(RPt)]
say("curated antibody->network-node assignments retained (antibody present on panel): ", length(ab2gene))
bad <- ab2gene[!(ab2gene %in% nodes$name)]
if (length(bad)) say("WARNING mapped to non-network symbol: ", paste(names(bad),bad,sep="->",collapse=","))
ab2gene <- ab2gene[ab2gene %in% nodes$name]
say("distinct network nodes covered by RPPA: ", length(unique(ab2gene)), " -> ", paste(sort(unique(unname(ab2gene))), collapse=","))

hub_any <- hubs$node[hubs$hub_degree | hubs$hub_betweenness | hubs$hub_centrality | hubs$hub_ffl]
say("network hubs (any hub flag): ", length(hub_any))
say("HUBS ON THE RPPA PANEL: ", paste(sort(intersect(hub_any, unname(ab2gene))), collapse=","))
say("HUBS NOT ON THE PANEL (protein-level not testable by RPPA): ", paste(sort(setdiff(hub_any[hub_any %in% nodes$name[nodes$type!='miRNA']], unname(ab2gene))), collapse=","))
say("NOTE: COL1A1 and COL3A1 are NOT on the RPPA panel (only COLLAGENVI/COL6A1 is).")

COV <- data.frame(antibody=names(ab2gene), network_node=unname(ab2gene),
                  is_phospho_antibody=is_phospho[names(ab2gene)],
                  node_type=nodes$type[match(unname(ab2gene), nodes$name)],
                  is_hub=unname(ab2gene) %in% hub_any,
                  n_BRCA_tumour_with_value=sapply(names(ab2gene), function(a) sum(!is.na(RPt[a,]))),
                  n_BRCA_solid_tissue_normal=n_normal,
                  tumour_vs_normal_testable = n_normal>0, stringsAsFactors=FALSE)
write.csv(COV, file.path(OUT,"rppa_panel_coverage.csv"), row.names=FALSE)
say("WROTE rppa_panel_coverage.csv rows=", nrow(COV))

## ---- TCGA mRNA / miRNA in the same patients ----
G <- readRDS(file.path(REV,"data","brca_gene_expr.rds"))
M <- readRDS(file.path(REV,"data","brca_mirna_expr_canonical.rds"))
phe <- readRDS(file.path(REV,"data","brca_pheno.rds"))
tcga_tum <- phe$sample[phe$sample_type=="Primary Tumor"]
s_rg <- sort(Reduce(intersect, list(tum, colnames(G), tcga_tum)))
s_rgm<- sort(Reduce(intersect, list(tum, colnames(G), colnames(M), tcga_tum)))
say("TCGA-BRCA primary tumours with RPPA + mRNA: ", length(s_rg))
say("TCGA-BRCA primary tumours with RPPA + mRNA + miRNA: ", length(s_rgm))

ctest <- function(x,y) { k <- !is.na(x)&!is.na(y)
  if (sum(k)<30 || sd(x[k])==0 || sd(y[k])==0) return(c(NA,NA,sum(k)))
  h <- suppressWarnings(cor.test(x[k],y[k],method="spearman",exact=FALSE))
  c(unname(h$estimate), h$p.value, sum(k)) }

## ---- F-a) RPPA protein vs its own mRNA (second-platform mRNA-protein concordance) ----
rows <- lapply(names(ab2gene), function(a) {
  g <- ab2gene[[a]]
  if (!(g %in% rownames(G))) return(NULL)
  v <- ctest(RPt[a, s_rg], G[g, s_rg])
  data.frame(antibody=a, gene=g, is_phospho=is_phospho[[a]], rho=v[1], p=v[2], n=v[3]) })
MP <- do.call(rbind, rows); MP$fdr <- p.adjust(MP$p,"BH"); MP <- MP[order(-MP$rho),]
write.csv(MP, file.path(OUT,"rppa_protein_vs_mRNA.csv"), row.names=FALSE)
say("WROTE rppa_protein_vs_mRNA.csv rows=", nrow(MP))
say("RPPA total-protein antibodies: median mRNA-protein rho = ", signif(median(MP$rho[!MP$is_phospho],na.rm=TRUE),4),
    " (n antibodies=", sum(!MP$is_phospho & !is.na(MP$rho)), ")")
for (i in seq_len(nrow(MP))) say(sprintf("  %-22s %-8s rho=%+.4f p=%.3g n=%d", MP$antibody[i],MP$gene[i],MP$rho[i],MP$p[i],MP$n[i]))

## ---- F-b) RELA pSer536 (independent platform) vs collagen / ECM mRNA ----
say("=== RELA pSer536 (RPPA NFKBP65_pS536) vs collagen mRNA in TCGA-BRCA tumours ===")
tg <- c("COL1A1","COL3A1","COL1A2","FN1","POSTN","VIM","VEGFA","SERPINE1","RELA","NFKB1","SP1","ETS1","IL6","CCL2","VCAM1")
E <- do.call(rbind, lapply(tg, function(g) {
  if (!(g %in% rownames(G))) return(data.frame(rppa_feature="NFKBP65_pS536", mRNA=g, rho=NA, p=NA, n=NA))
  v <- ctest(RPt["NFKBP65_pS536", s_rg], G[g, s_rg])
  data.frame(rppa_feature="NFKBP65_pS536", mRNA=g, rho=v[1], p=v[2], n=v[3]) }))
E$fdr <- p.adjust(E$p,"BH")
write.csv(E, file.path(OUT,"rppa_RELA_pS536_vs_targets.csv"), row.names=FALSE)
say("WROTE rppa_RELA_pS536_vs_targets.csv rows=", nrow(E))
for (i in seq_len(nrow(E))) say(sprintf("  NFKBP65_pS536 vs %-9s rho=%s p=%s fdr=%s n=%s", E$mRNA[i],
   ifelse(is.na(E$rho[i]),"NA",sprintf("%+.4f",E$rho[i])), ifelse(is.na(E$p[i]),"NA",sprintf("%.3g",E$p[i])),
   ifelse(is.na(E$fdr[i]),"NA",sprintf("%.3g",E$fdr[i])), E$n[i]))

## ---- F-c) miRNA -> target PROTEIN on the RPPA platform (independent replication of A) ----
say("=== miRNA->target edges whose target is on the RPPA panel: mRNA vs RPPA protein (n=", length(s_rgm), ") ===")
edges <- read.delim(file.path(REV,"data","canonical_edges.tsv"), stringsAsFactors=FALSE)
tier  <- read.delim(file.path(REV,"data","edge_evidence_tier.tsv"), stringsAsFactors=FALSE)
key <- function(a,b) paste(a,b,sep="\r")
tot_ab <- names(ab2gene)[!is_phospho[names(ab2gene)]]
gene2ab <- setNames(tot_ab, unname(ab2gene[tot_ab]))         # one total-protein antibody per gene
rppa_genes <- names(gene2ab)
nmeas <- sapply(rppa_genes, function(g) sum(!is.na(RPt[gene2ab[[g]], s_rgm])))
say("antibodies with ZERO measurements in BRCA tumours (dropped): ",
    paste(paste0(gene2ab[rppa_genes[nmeas < 30]], "(", rppa_genes[nmeas<30], ")"), collapse=","))
rppa_genes <- rppa_genes[nmeas >= 30]
rppa_genes <- rppa_genes[rppa_genes %in% rownames(G)]
say("RPPA network genes usable for edge test: ", length(rppa_genes), " -> ", paste(sort(rppa_genes), collapse=","))
mir_pool0 <- intersect(nodes$name[nodes$type=="miRNA"], rownames(M))
mir_pool0 <- mir_pool0[apply(M[mir_pool0, s_rgm, drop=FALSE],1,sd)>0]
say("network miRNAs measured with non-zero variance in the ", length(s_rgm), " matched TCGA tumours: ", length(mir_pool0))
mir_ok <- intersect(unique(edges$source[edges$edge_type=="miRNA_target"]), mir_pool0)
Ae <- edges[edges$edge_type=="miRNA_target" & edges$source %in% mir_ok & edges$target %in% rppa_genes, ]
Ae$evidence_tier <- setNames(tier$tier, key(tier$source,tier$target))[key(Ae$source,Ae$target)]
say("usable miRNA->target edges on RPPA panel: ", nrow(Ae), " (", length(unique(Ae$target)), " distinct targets, ", length(unique(Ae$source)), " miRNAs)")
ev <- function(a,b,Y,Ymat) { t(sapply(seq_along(a), function(i) ctest(M[a[i], s_rgm], Ymat[b[i], s_rgm]))) }
Rp <- t(sapply(seq_len(nrow(Ae)), function(i) ctest(M[Ae$source[i], s_rgm], RPt[gene2ab[[Ae$target[i]]], s_rgm])))
Rm <- t(sapply(seq_len(nrow(Ae)), function(i) ctest(M[Ae$source[i], s_rgm], G[Ae$target[i], s_rgm])))
Ae$rho_RPPAprotein <- Rp[,1]; Ae$p_RPPAprotein <- Rp[,2]; Ae$n_RPPAprotein <- Rp[,3]
Ae$rho_mRNA <- Rm[,1]; Ae$p_mRNA <- Rm[,2]; Ae$n_mRNA <- Rm[,3]
Ae$conc_protein <- Ae$rho_RPPAprotein < 0; Ae$conc_mRNA <- Ae$rho_mRNA < 0

## matched random-pair null within this small universe
mir_pool <- intersect(nodes$name[nodes$type=="miRNA"], rownames(M))
mir_pool <- mir_pool[apply(M[mir_pool, s_rgm, drop=FALSE],1,sd)>0]
dec_m <- as.integer(cut(rowMeans(M[mir_pool,s_rgm,drop=FALSE]), breaks=unique(quantile(rowMeans(M[mir_pool,s_rgm,drop=FALSE]),seq(0,1,.1))), include.lowest=TRUE)); names(dec_m) <- mir_pool
gmean <- sapply(rppa_genes, function(g) mean(RPt[gene2ab[[g]], s_rgm], na.rm=TRUE))
dec_g <- as.integer(cut(gmean, breaks=unique(quantile(gmean, seq(0,1,length.out=6))), include.lowest=TRUE)); names(dec_g) <- rppa_genes
Eset <- new.env(hash=TRUE,parent=emptyenv()); for (r in seq_len(nrow(edges))) assign(key(edges$source[r],edges$target[r]),TRUE,envir=Eset)
NPER <- 10L; nl <- list()
for (i in seq_len(nrow(Ae))) {
  cm <- setdiff(mir_pool[!is.na(dec_m) & abs(dec_m-dec_m[[Ae$source[i]]])<=1], Ae$source[i])
  cg <- setdiff(rppa_genes[!is.na(dec_g) & abs(dec_g-dec_g[[Ae$target[i]]])<=1], Ae$target[i])
  got<-0L; tries<-0L; seen<-character(0)
  while (got<NPER && tries<400L) { tries<-tries+1L
    a<-cm[sample.int(length(cm),1)]; bb<-cg[sample.int(length(cg),1)]; kk<-key(a,bb)
    if (kk %in% seen) next; if (exists(kk,envir=Eset,inherits=FALSE)) next
    seen<-c(seen,kk); nl[[length(nl)+1]] <- data.frame(edge_row=i, rnd_source=a, rnd_target=bb); got<-got+1L }
}
NL <- do.call(rbind, nl)
NLp <- t(sapply(seq_len(nrow(NL)), function(i) ctest(M[NL$rnd_source[i], s_rgm], RPt[gene2ab[[NL$rnd_target[i]]], s_rgm])))
NLm <- t(sapply(seq_len(nrow(NL)), function(i) ctest(M[NL$rnd_source[i], s_rgm], G[NL$rnd_target[i], s_rgm])))
NL$rho_protein <- NLp[,1]; NL$rho_mRNA <- NLm[,1]
say("null pairs: ", nrow(NL), " (", round(nrow(NL)/nrow(Ae),2), " per edge)")
mk <- function(lab, idx) {
  if (length(idx)<3) return(NULL)
  nlr <- NL[NL$edge_row %in% idx, ]
  data.frame(stratum=lab, n_edges=length(idx), n_null=nrow(nlr),
    conc_mRNA=mean(Ae$conc_mRNA[idx],na.rm=TRUE), null_conc_mRNA=mean(nlr$rho_mRNA<0,na.rm=TRUE),
    conc_RPPAprotein=mean(Ae$conc_protein[idx],na.rm=TRUE), null_conc_RPPAprotein=mean(nlr$rho_protein<0,na.rm=TRUE),
    mean_rho_mRNA=mean(Ae$rho_mRNA[idx],na.rm=TRUE), mean_rho_null_mRNA=mean(nlr$rho_mRNA,na.rm=TRUE),
    mean_rho_protein=mean(Ae$rho_RPPAprotein[idx],na.rm=TRUE), mean_rho_null_protein=mean(nlr$rho_protein,na.rm=TRUE),
    binom_p_protein_vs_null=binom.test(sum(Ae$conc_protein[idx],na.rm=TRUE), sum(!is.na(Ae$conc_protein[idx])),
                                        p=mean(nlr$rho_protein<0,na.rm=TRUE), alternative="greater")$p.value,
    binom_p_mRNA_vs_null=binom.test(sum(Ae$conc_mRNA[idx],na.rm=TRUE), sum(!is.na(Ae$conc_mRNA[idx])),
                                        p=mean(nlr$rho_mRNA<0,na.rm=TRUE), alternative="greater")$p.value) }
S <- rbind(mk("ALL", seq_len(nrow(Ae))),
           mk("strong", which(Ae$evidence_tier=="strong")),
           mk("weak", which(Ae$evidence_tier=="weak")),
           mk("predicted_only", which(Ae$evidence_tier=="predicted_only")))
write.csv(Ae, file.path(OUT,"rppa_miRNA_target_edges.csv"), row.names=FALSE)
write.csv(S,  file.path(OUT,"rppa_miRNA_target_summary.csv"), row.names=FALSE)
say("WROTE rppa_miRNA_target_edges.csv rows=", nrow(Ae), " ; rppa_miRNA_target_summary.csv rows=", nrow(S))
for (i in seq_len(nrow(S))) say(sprintf("  %-15s n=%3d  conc mRNA=%.3f (null %.3f)  conc RPPAprotein=%.3f (null %.3f) binomP_prot=%.3g",
   S$stratum[i],S$n_edges[i],S$conc_mRNA[i],S$null_conc_mRNA[i],S$conc_RPPAprotein[i],S$null_conc_RPPAprotein[i],S$binom_p_protein_vs_null[i]))
say("=== 23 DONE ==="); close(con)
