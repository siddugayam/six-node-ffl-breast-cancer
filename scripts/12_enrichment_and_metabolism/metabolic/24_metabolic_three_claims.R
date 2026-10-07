## 24_metabolic_three_claims.R --------------------------------------------
## Part D: tests of three metabolic associations (outputs metabolic_claim_*).
##  (i)   glycolysis / let-7b -> HK2
##  (ii)  antifolate resistance
##  (iii) AGE-RAGE as a metabolic-stress sensor (vs stromal confounding)
suppressPackageStartupMessages({library(data.table); library(limma); library(msigdbr)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
set.seed(1)
E   <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
ss  <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
gv  <- readRDS(file.path(RES,"metabolic_gsva_scores.rds"))
ph  <- readRDS(file.path(ROOT,"data/brca_pheno.rds"))
MI  <- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))
sets<- readRDS(file.path(RES,"metabolic_genesets.rds"))
tvn <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))
cor_all <- fread(file.path(RES,"metabolic_feature_correlations.csv"))
nodes <- fread(file.path(ROOT,"data/canonical_nodes.tsv"))
edges <- fread(file.path(ROOT,"data/canonical_edges.tsv"))

ph <- ph[match(colnames(E), ph$sample), ]; stopifnot(identical(ph$sample, colnames(E)))
grp <- factor(ifelse(ph$sample_type=="Primary Tumor","Tumor","Normal"), levels=c("Normal","Tumor"))
tum <- colnames(E)[grp=="Tumor"]; smp <- intersect(tum, colnames(MI))
cat("tumours:", length(tum), " paired with miRNA:", length(smp), "\n")

## gene-level tumour vs normal (limma) with the Rule-4 guard
des <- model.matrix(~ grp); fit <- eBayes(lmFit(E, des))
tt  <- topTable(fit, coef="grpTumor", number=Inf, sort.by="none")
stopifnot(identical(rownames(tt), rownames(E)), !all(grepl("^[0-9]+$", rownames(tt))))
GDE <- data.table(gene=rownames(tt), logFC=tt$logFC, t=tt$t, P=tt$P.Value,
                  FDR=p.adjust(tt$P.Value,"BH"),
                  mean_tumour=rowMeans(E[,grp=="Tumor"]), mean_normal=rowMeans(E[,grp=="Normal"]))
fwrite(GDE, file.path(RES,"metabolic_gene_level_tumour_vs_normal.csv"))

sp <- function(x,y){ ct <- suppressWarnings(cor.test(x,y,method="spearman", exact=FALSE))
                     c(rho=unname(ct$estimate), p=ct$p.value) }

## ================= CLAIM (i)  glycolysis + let-7b -> HK2 =================
cat("\n########## CLAIM (i) GLYCOLYSIS / let-7b -> HK2 ##########\n")
gly_ids <- grep("GLYCOLYSIS", rownames(ss), value=TRUE)
cat("\n-- glycolysis scores, tumour vs normal --\n")
print(tvn[set_id %in% gly_ids & score=="ssGSEA",
      .(set_id, n=n_genes_in_universe, d=round(cohens_d,3), t=round(t,2),
        p=signif(p_limma,3), FDR=signif(FDR_limma,3))][order(-d)])
tvn_ss <- tvn[score=="ssGSEA"][order(-cohens_d)][, rank_up := .I]
cat("rank of HALLMARK_GLYCOLYSIS among all", nrow(tvn_ss), "metabolic scores (by Cohen's d, up):",
    tvn_ss[set_id=="HALLMARK|HALLMARK_GLYCOLYSIS", rank_up], "\n")
cat("rank of KEGG glycolysis:", tvn_ss[grepl("hsa00010", set_id), rank_up], "\n")
cat("rank of CURATED_GLYCOLYSIS_CORE:", tvn_ss[set_id=="CURATED|CURATED_GLYCOLYSIS_CORE", rank_up], "\n")

cat("\n-- HK2 and glycolytic gene DE --\n")
gly_genes <- c("HK1","HK2","HK3","GCK","SLC2A1","PKM","LDHA","PFKP","ENO1","GAPDH","PGK1","ALDOA")
print(GDE[gene %in% gly_genes][order(-t), .(gene, logFC=round(logFC,3), t=round(t,2),
      FDR=signif(FDR,3), mean_T=round(mean_tumour,2), mean_N=round(mean_normal,2))])

cat("\n-- HK2 / let-7 in the published network --\n")
cat("HK2 a network node:", "HK2" %in% nodes$name, "\n")
cat("HK2 in any edge:", any(edges$source=="HK2" | edges$target=="HK2"), "\n")
let7_nodes <- grep("let-7", nodes$name, value=TRUE, ignore.case=TRUE)
cat("let-7 nodes:", paste(let7_nodes, collapse=", "), "\n")
l7e <- edges[grepl("let-7", source, ignore.case=TRUE) | grepl("let-7", target, ignore.case=TRUE)]
cat("edges involving let-7:", nrow(l7e), "; distinct let-7 targets:",
    length(unique(l7e[grepl("let-7",source,ignore.case=TRUE), target])), "\n")
cat("any let-7 -> HK2 edge:", nrow(l7e[target=="HK2"]), "\n")
fwrite(l7e, file.path(RES,"metabolic_claim_i_let7_network_edges.csv"))

cat("\n-- is HK2 a predicted let-7 target in MSigDB C3:MIR? --\n")
C3 <- as.data.table(msigdbr(species="Homo sapiens", collection="C3", subcollection="MIR:MIRDB"))
C3L<- as.data.table(msigdbr(species="Homo sapiens", collection="C3", subcollection="MIR:MIR_LEGACY"))
l7sets <- unique(c(grep("^LET7|LET_7", C3$gs_name, value=TRUE), grep("LET7|LET_7", C3L$gs_name, value=TRUE)))
c3all <- rbind(C3[,.(gs_name,gene_symbol)], C3L[,.(gs_name,gene_symbol)])
hk2_l7 <- c3all[gs_name %in% l7sets & gene_symbol=="HK2"]
cat("let-7 target sets in C3:MIR:", length(l7sets), "; those containing HK2:", nrow(hk2_l7), "\n")
if (nrow(hk2_l7)) print(hk2_l7$gs_name)
fwrite(data.table(let7_set=l7sets, contains_HK2 = l7sets %in% hk2_l7$gs_name),
       file.path(RES,"metabolic_claim_i_let7_targetset_membership.csv"))

cat("\n-- every let-7 arm vs glycolysis score and vs HK2 mRNA (", length(smp), "tumours) --\n")
let7 <- grep("^hsa-let-7", rownames(MI), value=TRUE)
hk2 <- E["HK2", smp]
res_i <- rbindlist(lapply(let7, function(m) {
  x <- MI[m, smp]
  a <- sp(x, ss["HALLMARK|HALLMARK_GLYCOLYSIS", smp])
  b <- sp(x, ss[grep("hsa00010", rownames(ss), value=TRUE)[1], smp])
  cc<- sp(x, ss["CURATED|CURATED_GLYCOLYSIS_CORE", smp])
  d <- sp(x, hk2)
  e <- sp(x, E["HK1", smp]); f <- sp(x, E["SLC2A1", smp]); g <- sp(x, E["LDHA", smp])
  data.table(miRNA=m, mean_log2RPM=mean(x),
    rho_HALLMARK_GLYCOLYSIS=a["rho"], p_HALLMARK_GLYCOLYSIS=a["p"],
    rho_KEGG_GLYCOLYSIS=b["rho"], p_KEGG_GLYCOLYSIS=b["p"],
    rho_CURATED_GLYCOLYSIS_CORE=cc["rho"], p_CURATED_GLYCOLYSIS_CORE=cc["p"],
    rho_HK2=d["rho"], p_HK2=d["p"], rho_HK1=e["rho"], p_HK1=e["p"],
    rho_SLC2A1=f["rho"], p_SLC2A1=f["p"], rho_LDHA=g["rho"], p_LDHA=g["p"])
}))
for (cl in grep("^p_", names(res_i), value=TRUE))
  res_i[[sub("^p_","FDR_",cl)]] <- p.adjust(res_i[[cl]], "BH")
setorder(res_i, rho_HK2)
print(res_i[, .(miRNA, mean=round(mean_log2RPM,2), rho_HK2=round(rho_HK2,3), FDR_HK2=signif(FDR_HK2,3),
                rho_HALLGLY=round(rho_HALLMARK_GLYCOLYSIS,3), FDR_HALLGLY=signif(FDR_HALLMARK_GLYCOLYSIS,3))])
cat("\nlet-7b-5p specifically: rho with HK2 =", round(res_i[miRNA=="hsa-let-7b-5p", rho_HK2],4),
    " p =", signif(res_i[miRNA=="hsa-let-7b-5p", p_HK2],3),
    "| rho with HALLMARK_GLYCOLYSIS =", round(res_i[miRNA=="hsa-let-7b-5p", rho_HALLMARK_GLYCOLYSIS],4),
    " p =", signif(res_i[miRNA=="hsa-let-7b-5p", p_HALLMARK_GLYCOLYSIS],3), "\n")
## let-7b's own DE
DEM <- fread(file.path(ROOT,"results/BRCA_DEX_mirnas.csv"))
print(DEM[grepl("let-7b|let-7a-5p", DEM[[1]])])
fwrite(res_i, file.path(RES,"metabolic_claim_i_let7_correlations.csv"))
fwrite(tvn[set_id %in% gly_ids], file.path(RES,"metabolic_claim_i_glycolysis_tvn.csv"))

## ================= CLAIM (ii) antifolate resistance ======================
cat("\n\n########## CLAIM (ii) ANTIFOLATE RESISTANCE / ONE-CARBON ##########\n")
af_ids <- grep("ANTIFOLATE|ONE_CARBON|FOLATE|hsa00670|hsa00790|hsa01523", rownames(ss), value=TRUE)
cat("\n-- folate / one-carbon scores, tumour vs normal --\n")
print(merge(tvn[set_id %in% af_ids & score=="ssGSEA",
  .(set_id, n=n_genes_in_universe, d=round(cohens_d,3), t=round(t,2), FDR=signif(FDR_limma,3))],
  tvn_ss[,.(set_id, rank_up)], by="set_id")[order(-d)])
af_genes <- c("DHFR","TYMS","SLC19A1","ABCC1","ABCC2","ABCC3","ABCC4","ABCC5","FPGS","GGH",
              "FOLR1","SLC46A1","ABCG2","MTHFD1","MTHFD2","SHMT1","SHMT2","ATIC","GART","MTR","MTHFR")
cat("\n-- antifolate determinant genes, tumour vs normal --\n")
print(GDE[gene %in% af_genes][order(-t), .(gene, logFC=round(logFC,3), t=round(t,2),
      FDR=signif(FDR,3), mean_T=round(mean_tumour,2), mean_N=round(mean_normal,2))])
cat("\n-- which antifolate determinants are network nodes? --\n")
cat(paste(intersect(af_genes, nodes$name), collapse=", "), "\n")
cat("(none listed = none in the 587-node network)\n")

mod_ids <- grep("^FFLCLASS\\|", rownames(ss), value=TRUE)
cat("\n-- folate/antifolate scores vs FFL module scores (", length(tum), "tumours) --\n")
res_ii <- rbindlist(lapply(af_ids, function(a) rbindlist(lapply(mod_ids, function(m) {
  r <- sp(ss[a,tum], ss[m,tum])
  ov <- length(intersect(sets[[a]], sets[[m]]))
  data.table(folate_set=a, module=m, rho=r["rho"], p=r["p"],
             n_shared_genes=ov, jaccard=ov/length(union(sets[[a]], sets[[m]])))
}))))
res_ii[, FDR := p.adjust(p,"BH")]
print(res_ii[folate_set %in% c("CURATED|CURATED_ANTIFOLATE_DETERMINANTS","CURATED|CURATED_ONE_CARBON_FOLATE",
             grep("hsa00670",af_ids,value=TRUE))][order(-abs(rho))][1:20,
      .(folate_set=substr(folate_set,1,45), module=sub("FFLCLASS\\|","",module),
        rho=round(rho,3), FDR=signif(FDR,3), shared=n_shared_genes)])
fwrite(res_ii, file.path(RES,"metabolic_claim_ii_folate_vs_modules.csv"))
fwrite(GDE[gene %in% af_genes], file.path(RES,"metabolic_claim_ii_determinant_genes_tvn.csv"))
fwrite(tvn[set_id %in% af_ids], file.path(RES,"metabolic_claim_ii_folate_tvn.csv"))

## ================= CLAIM (iii) AGE-RAGE ==================================
cat("\n\n########## CLAIM (iii) AGE-RAGE ##########\n")
age_id <- grep("hsa04933", rownames(ss), value=TRUE)
cat("AGE-RAGE set:", age_id, " n genes:", length(sets[[age_id]]), "\n")
coll <- grep("^COL[0-9]", sets[[age_id]], value=TRUE)
cat("collagen genes in the set:", length(coll), "/", length(sets[[age_id]]), ":", paste(coll, collapse=","), "\n")
cat("\n-- AGE-RAGE tumour vs normal --\n")
print(tvn[set_id==age_id, .(score, d=round(cohens_d,3), t=round(t,2), FDR=signif(FDR_limma,3),
     mean_T=round(mean_tumour,3), mean_N=round(mean_normal,3))])
cat("rank among all metabolic scores by Cohen's d (1 = most up):", tvn_ss[set_id==age_id, rank_up],
    "of", nrow(tvn_ss), "\n")

## AGE-RAGE vs context scores
ctx <- grep("^CONTROL\\|", rownames(ss), value=TRUE)
res_iii <- rbindlist(lapply(c(ctx, mod_ids, "HALLMARK|HALLMARK_GLYCOLYSIS",
                              "HALLMARK|HALLMARK_HYPOXIA"), function(f) {
  r <- sp(ss[age_id,tum], ss[f,tum])
  ov <- length(intersect(sets[[age_id]], sets[[f]]))
  data.table(partner=f, rho=r["rho"], p=r["p"], n_shared_genes=ov)
}))
r29  <- sp(ss[age_id,smp], MI["hsa-miR-29a-3p",smp])
res_iii <- rbind(res_iii, data.table(partner="hsa-miR-29a-3p",
                 rho=r29["rho"], p=r29["p"], n_shared_genes=NA_integer_))
res_iii[, FDR := p.adjust(p,"BH")]
cat("\n-- AGE-RAGE score vs context / module / miRNA --\n")
print(res_iii[order(-abs(rho)), .(partner=substr(partner,1,55), rho=round(rho,3),
      p=signif(p,3), FDR=signif(FDR,3), shared=n_shared_genes)])

## collagen-free AGE-RAGE re-score
suppressPackageStartupMessages(library(GSVA))
noncoll <- setdiff(sets[[age_id]], grep("^COL[0-9]", sets[[age_id]], value=TRUE))
sdv <- apply(E,1,sd); Ef <- E[sdv>0,]
sub_sets <- list(AGE_RAGE_FULL=sets[[age_id]], AGE_RAGE_NO_COLLAGEN=noncoll,
                 AGE_RAGE_COLLAGEN_ONLY=grep("^COL[0-9]", sets[[age_id]], value=TRUE))
sub_sets <- lapply(sub_sets, intersect, rownames(Ef))
ss2 <- gsva(ssgseaParam(Ef, sub_sets[lengths(sub_sets)>=5], alpha=0.25, normalize=TRUE, minSize=3), verbose=FALSE)
cafs <- ss["CONTROL|CONTROL_CAF_FULL", ]
tab <- rbindlist(lapply(rownames(ss2), function(v) {
  a <- sp(ss2[v,tum], cafs[tum]); b <- sp(ss2[v,tum], ss["CONTROL|CONTROL_EPITHELIAL",tum])
  cc<- sp(ss2[v,tum], ss["CONTROL|CONTROL_PROLIFERATION",tum])
  tv<- wilcox.test(ss2[v,grp=="Tumor"], ss2[v,grp=="Normal"])
  dd <- (mean(ss2[v,grp=="Tumor"])-mean(ss2[v,grp=="Normal"]))/sd(ss2[v,])
  data.table(variant=v, n_genes=length(sub_sets[[v]]),
    d_tumour_vs_normal=dd, p_tvn=tv$p.value,
    rho_CAF=a["rho"], p_CAF=a["p"], rho_EPI=b["rho"], rho_PROLIF=cc["rho"])
}))
cat("\n-- AGE-RAGE with and without its collagen genes --\n")
print(tab[, .(variant, n_genes, d_TvN=round(d_tumour_vs_normal,3), p_TvN=signif(p_tvn,3),
              rho_CAF=round(rho_CAF,3), rho_EPI=round(rho_EPI,3), rho_PROLIF=round(rho_PROLIF,3))])
fwrite(res_iii, file.path(RES,"metabolic_claim_iii_agerage_partners.csv"))
fwrite(tab, file.path(RES,"metabolic_claim_iii_agerage_collagen_decomposition.csv"))
fwrite(GDE[gene %in% sets[[age_id]]], file.path(RES,"metabolic_claim_iii_agerage_genes_tvn.csv"))
cat("\nAGE-RAGE member genes, tumour vs normal, top/bottom:\n")
print(GDE[gene %in% sets[[age_id]]][order(-t)][c(1:8, (.N-7):.N), .(gene, logFC=round(logFC,2), t=round(t,1), FDR=signif(FDR,2))])
cat("DONE 24\n")
