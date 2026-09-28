## 32_metabolic_verify.R --------------------------------------------------
## Independent re-derivation of the load-bearing numbers, without reusing the
## pipeline objects that produced them.
suppressPackageStartupMessages({library(data.table)})
ROOT <- "/path/to/revision"; RES <- file.path(ROOT,"results/v4")
E <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
ph<- readRDS(file.path(ROOT,"data/brca_pheno.rds")); ph <- ph[match(colnames(E), ph$sample),]
MI<- readRDS(file.path(ROOT,"data/brca_mirna_expr.rds"))
T_ <- ph$sample_type=="Primary Tumor"; N_ <- ph$sample_type=="Solid Tissue Normal"
cat("V1 sample counts: tumour", sum(T_), "normal", sum(N_), "\n")

## V2 HK2 differential expression, base-R t-test, vs the pipeline's limma value
h <- E["HK2",]; tt <- t.test(h[T_], h[N_])
gde <- fread(file.path(RES,"metabolic_gene_level_tumour_vs_normal.csv"))
cat("V2 HK2: base-R t =", round(tt$statistic,3), " p =", signif(tt$p.value,3),
    " mean diff =", round(mean(h[T_])-mean(h[N_]),4),
    " | pipeline limma t =", round(gde[gene=="HK2", t],3), " logFC =", round(gde[gene=="HK2", logFC],4), "\n")
## V2b cross-check against the manuscript's own DE table
dex <- fread(file.path(ROOT,"results/BRCA_DEX_genes.csv"))
cat("V2b HK2 in results/BRCA_DEX_genes.csv:", capture.output(print(dex[dex[[1]]=="HK2"])) [2:3], "\n")

## V3 let-7b-5p vs HK2 and vs glycolysis score, recomputed from scratch
ss <- readRDS(file.path(RES,"metabolic_ssgsea_scores.rds"))
smp <- intersect(colnames(E)[T_], colnames(MI))
cat("V3 paired tumours:", length(smp), "\n")
r1 <- cor(MI["hsa-let-7b-5p",smp], E["HK2",smp], method="spearman")
r2 <- cor(MI["hsa-let-7b-5p",smp], ss["HALLMARK|HALLMARK_GLYCOLYSIS",smp], method="spearman")
r3 <- cor(MI["hsa-let-7c-5p",smp], E["HK2",smp], method="spearman")
stored <- fread(file.path(RES,"metabolic_claim_i_let7_correlations.csv"))
cat("V3 let-7b-5p~HK2 recomputed", round(r1,4), "vs stored", round(stored[miRNA=="hsa-let-7b-5p", rho_HK2],4),
    "| let-7b-5p~HALLMARK_GLYCOLYSIS", round(r2,4), "vs", round(stored[miRNA=="hsa-let-7b-5p", rho_HALLMARK_GLYCOLYSIS],4),
    "| let-7c-5p~HK2", round(r3,4), "vs", round(stored[miRNA=="hsa-let-7c-5p", rho_HK2],4), "\n")

## V4 independent (mean-z) scoring of two claim-critical sets vs the ssGSEA result
sets <- readRDS(file.path(RES,"metabolic_genesets.rds"))
mz <- function(g){ Z <- t(scale(t(E[intersect(g, rownames(E)),,drop=FALSE]))); colMeans(Z, na.rm=TRUE) }
for (id in c("HALLMARK|HALLMARK_GLYCOLYSIS","CURATED|CURATED_GLYCOLYSIS_CORE",
             "CURATED|CURATED_ANTIFOLATE_DETERMINANTS",
             "KEGG_EXTRA|hsa04933_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS",
             "REACTOME|R-HSA-390918_PEROXISOMAL_LIPID_METABOLISM")) {
  z <- mz(sets[[id]])
  d_z  <- (mean(z[T_])-mean(z[N_]))/sd(z)
  d_ss <- fread(file.path(RES,"metabolic_tumour_vs_normal.csv"))[set_id==id & score=="ssGSEA", cohens_d]
  agree <- cor(z, ss[id,], method="spearman")
  cat(sprintf("V4 %-58s mean-z d=%+.3f  ssGSEA d=%+.3f  rho(mean-z,ssGSEA)=%.3f\n",
      substr(id,1,58), d_z, d_ss, agree))
}
## V5 AGE-RAGE ~ CAF recomputed with a completely independent CAF marker panel
caf2 <- intersect(c("FAP","PDGFRB","THY1","COL1A1","COL1A2","COL3A1","POSTN","DCN","LUM","FBLN1"), rownames(E))
cafz <- mz(caf2); tumall <- colnames(E)[T_]
cat("V5 AGE-RAGE(ssGSEA) ~ independent CAF mean-z, tumours: rho =",
    round(cor(ss["KEGG_EXTRA|hsa04933_AGE_RAGE_SIGNALING_PATHWAY_IN_DIABETIC_COMPLICATIONS",tumall],
              cafz[tumall], method="spearman"),3), "\n")
cat("DONE 32\n")
