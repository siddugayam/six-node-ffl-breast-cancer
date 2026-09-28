#!/usr/bin/env Rscript
# =============================================================================
# E2_pathway_gsva.R
# PART B: metabolic pathway ACTIVITY INFERRED FROM mRNA EXPRESSION.
#   *** THIS IS INFERENCE FROM TRANSCRIPT ABUNDANCE, NOT METABOLOMICS. ***
#   It says nothing about metabolite flux or concentration.
#
# Direct test of manuscript claim (i): let-7b -> HK2 -> glycolysis.
# Test of claim (ii): antifolate / one-carbon score in tumour, vs module scores.
# Inputs for claim (iii) AGE-RAGE handled here; NRF2/KEAP1 in E3.
# Output: results/multiomics/metabolic_pathway_scores.csv
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(GSVA); library(msigdbr); library(matrixStats)
})
ROOT <- "/path/to/revision"
OUT  <- file.path(ROOT,"results/multiomics")
CACHE<- file.path(ROOT,"cache/metabolomics")
GMT  <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/PANCAN/c2.all.v2026.1.Hs.symbols.gmt"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
con <- file(file.path(ROOT,"logs/E2_pathway_gsva.log"), open="wt")
sink(con, split=TRUE); sink(con, type="message")
cat("=== E2_pathway_gsva.R ", format(Sys.time()), " ===\n")
set.seed(42)

## ------------------------------------------------------------------ data
ge <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
ph <- as.data.table(readRDS(file.path(ROOT,"data/brca_pheno.rds")))
mi <- readRDS(file.path(ROOT,"data/brca_mirna_expr_canonical.rds"))
cat("gene expr:", dim(ge), " miRNA:", dim(mi), "\n")
ph <- ph[sample %in% colnames(ge)]
ge <- ge[, ph$sample]
cat("tumour/normal:\n"); print(table(ph$sample_type))
is_tum <- ph$sample_type=="Primary Tumor"

## ------------------------------------------------------------------ gene sets
read_gmt <- function(f){
  ll <- strsplit(readLines(f), "\t")
  setNames(lapply(ll, function(x) unique(x[-c(1,2)])), sapply(ll,`[`,1))
}
c2 <- read_gmt(GMT); cat("C2 sets:", length(c2), "\n")
H  <- as.data.table(msigdbr(species="Homo sapiens", collection="H"))
hs <- split(H$gene_symbol, H$gs_name)
kegg <- jsonlite::fromJSON(file.path(CACHE,"kegg_sets.json"), simplifyVector=FALSE)

# split KEGG Antifolate resistance (hsa01523) into its two mechanistic arms:
# the folate/nucleotide enzymes+transporters, and the NF-kB/cytokine arm that
# the pathway contains only because antifolates suppress cytokine release.
AF_ALL  <- unlist(kegg$hsa01523$genes)
AF_NFKB <- c("CHUK","IKBKB","IKBKG","IL1B","IL6","NFKB1","RELA","TNF")
AF_FOL  <- setdiff(AF_ALL, AF_NFKB)
cat("\nKEGG hsa01523 Antifolate resistance: ", length(AF_ALL), " genes = ",
    length(AF_FOL), " folate/transport + ", length(AF_NFKB), " NF-kB/cytokine\n", sep="")
cat("  folate arm: ", paste(AF_FOL, collapse=", "), "\n")
cat("  NF-kB  arm: ", paste(AF_NFKB, collapse=", "), "\n")

# curated canonical NRF2 (NFE2L2) transcriptional targets - CURATED BY ME, stated as such
NRF2_CORE <- c("NQO1","GCLC","GCLM","TXNRD1","TXN","SLC7A11","GSR","G6PD","PGD",
               "ME1","AKR1B10","AKR1C1","AKR1C2","AKR1C3","SRXN1","PRDX1","FTL",
               "FTH1","HMOX1","OSGIN1","CBR1","GPX2","ABCC2","SQSTM1","GSTP1",
               "GSTM3","UGT1A6","EPHX1","CES1","TALDO1","IDH1","PGD")

sets <- list(
  KEGG_GLYCOLYSIS_GLUCONEOGENESIS_live      = unlist(kegg$hsa00010$genes),
  KEGG_ONE_CARBON_POOL_BY_FOLATE_live       = unlist(kegg$hsa00670$genes),
  KEGG_ANTIFOLATE_RESISTANCE_hsa01523_live  = AF_ALL,
  KEGG_ANTIFOLATE_RESISTANCE_FOLATE_ARM     = AF_FOL,
  KEGG_ANTIFOLATE_RESISTANCE_NFKB_ARM       = AF_NFKB,
  KEGG_AGE_RAGE_hsa04933_live               = unlist(kegg$hsa04933$genes),
  HALLMARK_GLYCOLYSIS                       = hs[["HALLMARK_GLYCOLYSIS"]],
  HALLMARK_OXIDATIVE_PHOSPHORYLATION        = hs[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
  HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY  = hs[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]],
  KEGG_OXIDATIVE_PHOSPHORYLATION            = c2[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
  KEGG_GLYCOLYSIS_GLUCONEOGENESIS_msigdb    = c2[["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]],
  KEGG_ONE_CARBON_POOL_BY_FOLATE_msigdb     = c2[["KEGG_ONE_CARBON_POOL_BY_FOLATE"]],
  KEGG_FOLATE_BIOSYNTHESIS                  = c2[["KEGG_FOLATE_BIOSYNTHESIS"]],
  REACTOME_METABOLISM_OF_FOLATE_AND_PTERINES= c2[["REACTOME_METABOLISM_OF_FOLATE_AND_PTERINES"]],
  WP_AGERAGE_PATHWAY                        = c2[["WP_AGERAGE_PATHWAY"]],
  REACTOME_KEAP1_NFE2L2_PATHWAY             = c2[["REACTOME_KEAP1_NFE2L2_PATHWAY"]],
  REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2= c2[["REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2"]],
  REACTOME_NFE2L2_REGULATING_ANTI_OXIDANT_DETOXIFICATION_ENZYMES =
                       c2[["REACTOME_NFE2L2_REGULATING_ANTI_OXIDANT_DETOXIFICATION_ENZYMES"]],
  WP_NRF2_PATHWAY                           = c2[["WP_NRF2_PATHWAY"]],
  SINGH_NFE2L2_TARGETS                      = c2[["SINGH_NFE2L2_TARGETS"]],
  WP_OXIDATIVE_STRESS_RESPONSE              = c2[["WP_OXIDATIVE_STRESS_RESPONSE"]],
  NFE2L2_CORE_TARGETS_curated_here          = NRF2_CORE
)
sets <- lapply(sets, function(x) unique(intersect(x, rownames(ge))))
cat("\nGene sets scored (size after intersecting the expression matrix):\n")
for (n in names(sets)) cat(sprintf("  %-62s %d\n", n, length(sets[[n]])))
stopifnot(all(lengths(sets) >= 5))

## ------------------------------------------------------------------ GSVA
cat("\nrunning GSVA on", ncol(ge), "samples ...\n")
gp <- gsvaParam(exprData=ge, geneSets=sets, kcdf="Gaussian", minSize=5, maxSize=500)
S  <- gsva(gp, verbose=FALSE)
cat("GSVA score matrix:", dim(S), "\n")
saveRDS(S, file.path(OUT,"metabolic_gsva_scores.rds"))

## ------------------------------------------------------------------ T vs N
tn <- rbindlist(lapply(rownames(S), function(p){
  x <- S[p, is_tum]; y <- S[p, !is_tum]
  w <- wilcox.test(x, y, exact=FALSE)
  data.table(pathway=p, n_tumour=length(x), n_normal=length(y),
             mean_tumour=mean(x), mean_normal=mean(y), delta=mean(x)-mean(y),
             W=unname(w$statistic), p=w$p.value)
}))
tn[, FDR := p.adjust(p,"BH")]
tn[, direction := fifelse(FDR<0.05 & delta>0,"UP_in_tumour",
                  fifelse(FDR<0.05 & delta<0,"DOWN_in_tumour","ns"))]
setorder(tn, p)
cat("\n================ TUMOUR vs NORMAL GSVA SCORE (Wilcoxon, BH) ================\n")
print(tn[, .(pathway, n_tumour, n_normal, mean_tumour=round(mean_tumour,4),
             mean_normal=round(mean_normal,4), delta=round(delta,4),
             p=signif(p,3), FDR=signif(FDR,3), direction)], nrows=50)

## ------------------------------------------------------------------ hub correlations
hubs <- fread(file.path(ROOT,"results/network_topology_hubs.csv"))[is_hub==TRUE]
cat("\nhubs:", nrow(hubs), " (", hubs[,.N,by=type][,paste(type,N,collapse="; ")], ")\n")
tum <- ph[sample_type=="Primary Tumor"]$sample
Stum <- S[, tum, drop=FALSE]

# gene/TF hubs: expression matrix
gh <- hubs[type %in% c("TF","Gene") & name %in% rownames(ge)]$name
# miRNA hubs: miRNA matrix, intersect samples
mi_s <- intersect(tum, colnames(mi))
cat("tumours with both mRNA and miRNA:", length(mi_s), "\n")
mh <- hubs[type=="miRNA" & name %in% rownames(mi)]$name
cat("hub TFs/Genes testable:", length(gh), " hub miRNAs testable:", length(mh), "\n")

cor_block <- function(featmat, feats, scoremat, samples, ftype){
  rbindlist(lapply(feats, function(f){
    v <- as.numeric(featmat[f, samples])
    rbindlist(lapply(rownames(scoremat), function(p){
      s <- as.numeric(scoremat[p, samples])
      ct <- suppressWarnings(cor.test(v, s, method="spearman", exact=FALSE))
      data.table(feature=f, feature_type=ftype, pathway=p, n=length(samples),
                 rho=unname(ct$estimate), p=ct$p.value)
    }))
  }))
}
cc <- rbind(cor_block(ge, gh, Stum, tum, "hub_TF_or_Gene"),
            cor_block(mi, mh, S[, mi_s, drop=FALSE], mi_s, "hub_miRNA"))
cc[, FDR := p.adjust(p,"BH")]
cat("\nhub-vs-pathway correlations computed:", nrow(cc), " (FDR<0.05:", sum(cc$FDR<0.05), ")\n")

## ------------------------------------------------------------------ CLAIM (i)
cat("\n\n########################################################################\n")
cat("### CLAIM (i)  \"modulate glycolysis via the let-7b-HK2 axis\"\n")
cat("########################################################################\n")
nodes <- fread(file.path(ROOT,"data/canonical_nodes.tsv"))
edges <- fread(file.path(ROOT,"data/canonical_edges.tsv"))
cat("\n[1] IS HK2 IN THE NETWORK AT ALL?\n")
cat("    HK2 in canonical_nodes.tsv : ", sum(nodes$name=="HK2"), " rows\n", sep="")
cat("    HK2 in canonical_edges.tsv : ", sum(edges$source=="HK2" | edges$target=="HK2"), " rows\n", sep="")
cat("    other hexokinases (HK1/HK3/HKDC1/GCK): ",
    sum(nodes$name %in% c("HK1","HK3","HKDC1","GCK")), " nodes\n", sep="")
cat("    glycolytic genes present in the network: ",
    paste(intersect(nodes$name, sets$KEGG_GLYCOLYSIS_GLUCONEOGENESIS_live), collapse=", "), "\n")
cat("    let-7b targets recorded in the network: ",
    paste(edges[source=="hsa-let-7b" & edge_type=="miRNA_target"]$target, collapse=", "), "\n")

cat("\n[2] IS HK2 DIFFERENTIALLY EXPRESSED (limma, TCGA-BRCA)?\n")
de <- fread(file.path(ROOT,"results/BRCA_DEX_genes.csv"))
print(de[feature %in% c("HK1","HK2","HK3","HKDC1","GCK","SLC2A1","LDHA","PKM","PFKP","GAPDH","ENO1")][
        , .(feature, logFC=round(logFC,4), P.Value=signif(P.Value,3), adj.P.Val=signif(adj.P.Val,3))])
dem <- fread(file.path(ROOT,"results/BRCA_DEX_mirnas.csv"))
cat("\n    let-7b DE: "); print(dem[feature=="hsa-let-7b", .(logFC=round(logFC,4), adj.P.Val=signif(adj.P.Val,3))])
# direct re-test of HK2 T vs N so the number is reproduced here, not just cited
hk <- wilcox.test(ge["HK2", is_tum], ge["HK2", !is_tum], exact=FALSE)
cat("    direct Wilcoxon HK2 tumour vs normal: p =", signif(hk$p.value,3),
    " mean T =", round(mean(ge["HK2",is_tum]),3), " mean N =", round(mean(ge["HK2",!is_tum]),3), "\n")

cat("\n[3] DOES let-7b CORRELATE WITH HK2 IN TUMOURS? (the axis itself)\n")
let7b <- as.numeric(mi["hsa-let-7b", mi_s])
hk2   <- as.numeric(ge["HK2", mi_s])
ct <- cor.test(let7b, hk2, method="spearman", exact=FALSE)
cat(sprintf("    hsa-let-7b vs HK2 : rho = %+.4f, p = %.3g, n = %d  --> predicted NEGATIVE\n",
            ct$estimate, ct$p.value, length(mi_s)))
# every hexokinase + core glycolytic gene, and let-7e for completeness
axis <- rbindlist(lapply(c("hsa-let-7b","hsa-let-7e","hsa-let-7a","hsa-let-7c"), function(m){
  rbindlist(lapply(c("HK1","HK2","HK3","HKDC1","GCK","SLC2A1","LDHA","PKM","PFKP","GAPDH","ENO1","PGK1"), function(g){
    if (!(g %in% rownames(ge)) || !(m %in% rownames(mi))) return(NULL)
    c2t <- suppressWarnings(cor.test(as.numeric(mi[m,mi_s]), as.numeric(ge[g,mi_s]),
                                     method="spearman", exact=FALSE))
    data.table(miRNA=m, gene=g, n=length(mi_s), rho=unname(c2t$estimate), p=c2t$p.value)
  }))
}))
axis[, FDR := p.adjust(p,"BH")]
cat("\n    let-7 family vs glycolytic genes (Spearman, tumours only):\n")
print(axis[order(miRNA,gene)][, .(miRNA, gene, n, rho=round(rho,4), p=signif(p,3), FDR=signif(FDR,3))], nrows=60)

cat("\n[4] DOES let-7b CORRELATE (NEGATIVELY) WITH THE GLYCOLYSIS SCORE?\n")
g_res <- cc[feature=="hsa-let-7b" & grepl("GLYCOL", pathway)]
if (!nrow(g_res)) {  # let-7b may not be flagged a hub; compute directly
  g_res <- rbindlist(lapply(grep("GLYCOL", rownames(S), value=TRUE), function(p){
    ct <- suppressWarnings(cor.test(let7b, as.numeric(S[p, mi_s]), method="spearman", exact=FALSE))
    data.table(feature="hsa-let-7b", feature_type="miRNA", pathway=p, n=length(mi_s),
               rho=unname(ct$estimate), p=ct$p.value, FDR=NA_real_)
  }))
}
print(g_res[, .(feature, pathway, n, rho=round(rho,4), p=signif(p,3))])
cat("\n    NOTE: manuscript predicts let-7b represses HK2 -> let-7b should be NEGATIVELY\n")
cat("    correlated with glycolysis. Positive rho contradicts the claim.\n")

## ------------------------------------------------------------------ CLAIM (ii)
cat("\n\n########################################################################\n")
cat("### CLAIM (ii)  \"Antifolate resistance\" enrichment => metabolic plasticity\n")
cat("########################################################################\n")
netg <- nodes[type %in% c("TF","Gene")]$name
cat("\nNetwork genes inside KEGG hsa01523 (30 genes):",
    paste(intersect(netg, AF_ALL), collapse=", "), "\n")
cat("  of which FOLATE arm  :", paste(intersect(netg, AF_FOL), collapse=", "), "\n")
cat("  of which NF-kB  arm  :", paste(intersect(netg, AF_NFKB), collapse=", "), "\n")
cat("\nTumour-vs-normal for the antifolate sets and their arms:\n")
print(tn[grepl("ANTIFOLATE|FOLATE|ONE_CARBON", pathway),
         .(pathway, delta=round(delta,4), p=signif(p,3), FDR=signif(FDR,3), direction)])

cat("\nAntifolate-pathway ENZYME genes, tumour vs normal (limma):\n")
print(de[feature %in% c("TYMS","DHFR","GART","ATIC","SHMT1","SHMT2","MTHFR","MTR",
                        "FPGS","GGH","SLC19A1","SLC46A1","FOLR1","ABCC1","ABCG2")][
        order(adj.P.Val)][, .(feature, logFC=round(logFC,3), adj.P.Val=signif(adj.P.Val,3))])

cat("\n--- antifolate/one-carbon score vs FFL MODULE SCORES ---\n")
ms <- readRDS(file.path(ROOT,"data/tcga_module_scores.rds"))
msg <- ms$gsva; common <- intersect(colnames(msg), colnames(Stum))
cat("samples shared with module-score matrix:", length(common), "\n")
modcor <- rbindlist(lapply(rownames(msg), function(m){
  rbindlist(lapply(grep("ANTIFOLATE|ONE_CARBON|FOLATE|GLYCOL|AGE_RAGE|AGERAGE", rownames(S), value=TRUE), function(p){
    ct <- suppressWarnings(cor.test(as.numeric(msg[m,common]), as.numeric(S[p,common]),
                                    method="spearman", exact=FALSE))
    data.table(module=m, pathway=p, n=length(common), rho=unname(ct$estimate), p=ct$p.value)
  }))
}))
modcor[, FDR := p.adjust(p,"BH")]
cat("\nHigher-order-only vs 3-node module, correlation with metabolic scores:\n")
print(modcor[module %in% c("FFL_3node","FFL_higher_order_only","FFL_6node_all","MS_exemplar_module")][
      order(pathway,module)][, .(module, pathway, n, rho=round(rho,3), FDR=signif(FDR,3))], nrows=60)

## ------------------------------------------------------------------ CLAIM (iii) AGE-RAGE
cat("\n\n########################################################################\n")
cat("### CLAIM (iii)  AGE-RAGE as a metabolic-stress sensor\n")
cat("########################################################################\n")
cat("\nAGER (the RAGE receptor itself) in network nodes:", sum(nodes$name=="AGER"), "\n")
cat("AGER DE in TCGA-BRCA:\n"); print(de[feature=="AGER", .(feature, logFC=round(logFC,3), adj.P.Val=signif(adj.P.Val,3))])
cat("\nNetwork genes inside KEGG hsa04933 (101 genes): n =",
    length(intersect(netg, unlist(kegg$hsa04933$genes))), "\n  ",
    paste(intersect(netg, unlist(kegg$hsa04933$genes)), collapse=", "), "\n")
cat("\nAGE-RAGE score tumour vs normal:\n")
print(tn[grepl("AGE_RAGE|AGERAGE", pathway), .(pathway, delta=round(delta,4), p=signif(p,3), FDR=signif(FDR,3), direction)])
cat("\nCorrelation of AGE-RAGE score with glycolysis / OXPHOS score (tumours), i.e.\n")
cat("is it behaving like a metabolic-stress readout or like a collagen/ECM readout?\n")
for (p2 in c("HALLMARK_GLYCOLYSIS","HALLMARK_OXIDATIVE_PHOSPHORYLATION")){
  ct <- cor.test(as.numeric(Stum["KEGG_AGE_RAGE_hsa04933_live",]), as.numeric(Stum[p2,]),
                 method="spearman", exact=FALSE)
  cat(sprintf("   AGE-RAGE vs %-38s rho = %+.3f  p = %.3g\n", p2, ct$estimate, ct$p.value))
}
for (g in c("COL1A1","COL3A1","FN1","SERPINE1")){
  ct <- cor.test(as.numeric(Stum["KEGG_AGE_RAGE_hsa04933_live",]), as.numeric(ge[g,tum]),
                 method="spearman", exact=FALSE)
  cat(sprintf("   AGE-RAGE vs %-38s rho = %+.3f  p = %.3g\n", g, ct$estimate, ct$p.value))
}

## ------------------------------------------------------------------ write
tn[, analysis := "tumour_vs_normal_GSVA"]
cc[, analysis := "hub_vs_pathway_spearman_tumours"]
axis[, analysis := "let7_vs_glycolytic_gene_spearman_tumours"]
modcor[, analysis := "module_vs_pathway_spearman_tumours"]
setnames(axis, c("miRNA","gene"), c("feature","pathway"))
axis[, `:=`(feature_type="miRNA", n=n)]
out <- rbindlist(list(
  tn[, .(analysis, feature=NA_character_, feature_type=NA_character_, pathway,
         n1=n_tumour, n2=n_normal, stat=delta, p, FDR,
         note="delta = mean GSVA tumour - mean GSVA normal")],
  cc[, .(analysis, feature, feature_type, pathway, n1=n, n2=NA_integer_, stat=rho, p, FDR,
         note="Spearman rho, primary tumours only")],
  axis[, .(analysis, feature, feature_type, pathway, n1=n, n2=NA_integer_, stat=rho, p, FDR,
           note="Spearman rho miRNA vs glycolytic GENE (not a pathway score)")],
  modcor[, .(analysis, feature=module, feature_type="FFL_module_score", pathway,
             n1=n, n2=NA_integer_, stat=rho, p, FDR,
             note="Spearman rho FFL module GSVA vs metabolic pathway GSVA")]
), use.names=TRUE)
out[, evidence_class := "INFERRED FROM mRNA EXPRESSION (GSVA) - NOT metabolomics"]
fwrite(out, file.path(OUT,"metabolic_pathway_scores.csv"))
cat("\n\nWROTE", file.path(OUT,"metabolic_pathway_scores.csv"), " rows:", nrow(out), "\n")
fwrite(as.data.table(as.data.frame(S), keep.rownames="pathway"),
       file.path(OUT,"metabolic_gsva_persample.csv"))
cat("WROTE metabolic_gsva_persample.csv  rows:", nrow(S), " cols:", ncol(S)+1, "\n")
sd <- data.table(pathway=names(sets), n_genes_in_set=lengths(sets),
                 genes=sapply(sets, paste, collapse="|"))
fwrite(sd, file.path(OUT,"metabolic_pathway_genesets.csv"))
cat("WROTE metabolic_pathway_genesets.csv  rows:", nrow(sd), "\n")
cat("\nDONE ", format(Sys.time()), "\n")
sink(type="message"); sink(); close(con)
