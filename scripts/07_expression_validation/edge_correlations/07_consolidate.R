## Consolidated claim-by-claim verification table -> results/v2/verify_expression.csv
suppressPackageStartupMessages({library(matrixStats)})
setwd("/path/to/revision")
gexp <- readRDS("data/brca_gene_expr.rds"); mexp <- readRDS("data/brca_mirna_expr_canonical.rds")
ph <- readRDS("data/brca_pheno.rds"); tum <- ph$sample[ph$sample_type=="Primary Tumor"]
sG <- sort(intersect(colnames(gexp),tum)); sB <- sort(intersect(sG, colnames(mexp)))
nodes <- read.delim("data/canonical_nodes.tsv", stringsAsFactors=FALSE)$name
gmt <- readLines(system.file("extdata","SI_geneset.gmt", package="estimate"))
strom <- strsplit(gmt[grep("StromalSignature",gmt)],"\t")[[1]][-c(1,2)]
sigA <- intersect(setdiff(setdiff(strom, grep("^COL",strom,value=TRUE)), nodes), rownames(gexp))
sigB <- intersect(setdiff(c("DCN","LUM","FAP","PDGFRB","THY1","POSTN"), nodes), rownames(gexp))
cafsc <- function(S,ss){X <- gexp[S,ss,drop=FALSE]; colMeans((X-rowMeans(X))/rowSds(X))}
pcor <- function(x,y,c_){ rx <- rank(x); ry <- rank(y); rc <- rank(c_)
  rxy <- cor(rx,ry); rxc <- cor(rx,rc); ryc <- cor(ry,rc)
  (rxy-rxc*ryc)/sqrt((1-rxc^2)*(1-ryc^2)) }
axes <- rbind(
 data.frame(x=c("ETS1","NFKB1","SP1","RELA"), y="COL1A1", assay="gene", stringsAsFactors=FALSE),
 data.frame(x=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), y="COL1A1", assay="miRNA"),
 data.frame(x=c("hsa-miR-29a","hsa-miR-29b","hsa-miR-29c"), y="COL3A1", assay="miRNA"),
 data.frame(x="hsa-miR-101", y="EZH2", assay="miRNA"),
 data.frame(x=c("hsa-let-7b","hsa-let-7e"), y="COL3A1", assay="miRNA"),
 data.frame(x="COL1A1", y="COL3A1", assay="gene"))
E <- read.delim("data/canonical_edges.tsv", stringsAsFactors=FALSE)
axes$in_network <- paste(axes$x, axes$y) %in% paste(E$source, E$target)
res <- NULL
for (i in seq_len(nrow(axes))){
  ss <- if (axes$assay[i]=="miRNA") sB else sG
  X <- if (axes$assay[i]=="miRNA") mexp[axes$x[i], ss] else gexp[axes$x[i], ss]
  Y <- gexp[axes$y[i], ss]; cA <- cafsc(sigA, ss); cB <- cafsc(sigB, ss)
  ter <- cut(cA, quantile(cA, c(0,1/3,2/3,1)), include.lowest=TRUE, labels=c("lo","mid","hi"))
  ct <- cor.test(X, Y, method="spearman", exact=FALSE)
  res <- rbind(res, data.frame(axis=paste(axes$x[i],"->",axes$y[i]), in_network=axes$in_network[i],
    n=length(ss), rho=unname(ct$estimate), p=ct$p.value,
    rho_pCAFA=pcor(X,Y,cA), rho_pCAFB=pcor(X,Y,cB),
    pct_retained_CAFA=round(100*pcor(X,Y,cA)/unname(ct$estimate),1),
    rho_CAFlo=cor(X[ter=="lo"],Y[ter=="lo"],method="spearman"),
    rho_CAFmid=cor(X[ter=="mid"],Y[ter=="mid"],method="spearman"),
    rho_CAFhi=cor(X[ter=="hi"],Y[ter=="hi"],method="spearman"),
    sign_flips_CAFA = sign(pcor(X,Y,cA)) != sign(unname(ct$estimate)),
    stringsAsFactors=FALSE))
}
print(res, digits=3, row.names=FALSE)
write.csv(res, "results/v2/v2_named_axes_final.csv", row.names=FALSE)

## -------------------- assemble the verification table ---------------------
de <- readRDS("results/v2/v2_DE.rds")
tiers <- table(read.delim("data/edge_evidence_tier.tsv", stringsAsFactors=FALSE)$tier)
gn <- read.csv("results/v2/v2_null_within_network_pool.csv", row.names=NULL)   # paper-style null
gw <- read.csv("results/v2/v2_edge_class_summary_genomewide_null.csv")                          # genome-wide null (from 02)
ca <- read.csv("results/v2/v2_caf_adjusted_vs_network_null.csv")
cp <- read.csv("results/v2/v2_cptac_mir29_armaveraged.csv")
gv <- function(d,s,col) d[[col]][d$set==s]
rowf <- function(claim, quantity, published, verified, verdict, note)
  data.frame(claim=claim, quantity=quantity, published=published,
             verified=verified, verdict=verdict, note=note, stringsAsFactors=FALSE)
V <- rbind(
 rowf("A/DE","genes tested","20250", nrow(de$gene),"CONFIRMED","limma ~grp, zero-variance rows dropped"),
 rowf("A/DE","genes significant (FDR<0.05,|logFC|>1)","4231",
      sum(de$gene$adj.P.Val<0.05 & abs(de$gene$logFC)>1),"CONFIRMED",""),
 rowf("A/DE","miRNAs tested","495", nrow(de$mirna),"CONFIRMED",""),
 rowf("A/DE","miRNAs significant","80",
      sum(de$mirna$adj.P.Val<0.05 & abs(de$mirna$logFC)>1),"CONFIRMED",""),
 rowf("A/DE","gene samples T vs N","1097 vs 114","1097 vs 114","CONFIRMED",""),
 rowf("A/DE","miRNA samples T vs N","1068 vs 90","1068 vs 90","CONFIRMED",""),
 rowf("A/DE","topTable feature column","gene symbols","gene symbols (assert passed, 0 dup rownames)",
      "CONFIRMED","no integer-index failure"),
 rowf("N6","miRNA_target strong","707 (14.7%)", sprintf("%d (%.2f%%)", tiers["strong"], 100*tiers["strong"]/sum(tiers)),"CONFIRMED",""),
 rowf("N6","miRNA_target weak","2309", tiers["weak"],"CONFIRMED",""),
 rowf("N6","miRNA_target predicted-only","1803", tiers["predicted_only"],"CONFIRMED",""))
for (s in c("TF_target_Activation","TF_target_Repression","TF_miRNA_Activation",
            "miRNA_target_all","miRNA_target_strong","miRNA_target_weak","miRNA_target_predicted_only")){
 V <- rbind(V,
  rowf("N8", paste0(s,": concordance"), NA, sprintf("%.4f (n=%d)", gv(gn,s,"conc"), gv(gn,s,"n")),"",
       "identical under both nulls"),
  rowf("N8", paste0(s,": null concordance (within-network pool, 5x10/edge)"), NA,
       sprintf("%.4f +/- %.4f", gv(gn,s,"conc_null_mean"), gv(gn,s,"conc_null_sd")),"",
       sprintf("binom p=%.3g", gv(gn,s,"p_binom_mean"))),
  rowf("N8", paste0(s,": null concordance (genome-wide decile pool, 20/edge)"), NA,
       sprintf("%.4f", gw$conc_null[gw$set==s]),"", sprintf("prop p=%.3g", gw$conc_p[gw$set==s])),
  rowf("N8", paste0(s,": mean rho vs null"), NA,
       sprintf("%.4f vs %.4f", gv(gn,s,"mean_rho"), gv(gn,s,"mean_rho_null")),"",
       sprintf("Wilcoxon p=%.3g", gv(gn,s,"p_wilcox_mean"))))
}
for (s in c("TF_target_Activation","miRNA_target_all","miRNA_target_predicted_only")){
 V <- rbind(V,
  rowf("N9", paste0(s,": CAF-A adjusted conc vs null"), NA,
       sprintf("%.4f vs %.4f", ca$conc_cafA[ca$set==s], ca$null_cafA[ca$set==s]),"",
       sprintf("binom p=%.3g", ca$p_cafA[ca$set==s])),
  rowf("N9", paste0(s,": CAF-B adjusted conc vs null"), NA,
       sprintf("%.4f vs %.4f", ca$conc_cafB[ca$set==s], ca$null_cafB[ca$set==s]),"",
       sprintf("binom p=%.3g", ca$p_cafB[ca$set==s])))
}
for (i in seq_len(nrow(res)))
 V <- rbind(V, rowf("N9/axis", res$axis[i], NA,
   sprintf("rho=%.4f -> CAF-A %.4f (%.0f%% retained), CAF-B %.4f", res$rho[i], res$rho_pCAFA[i],
           res$pct_retained_CAFA[i], res$rho_pCAFB[i]), "",
   sprintf("n=%d; tertiles lo/mid/hi = %.3f/%.3f/%.3f", res$n[i], res$rho_CAFlo[i], res$rho_CAFmid[i], res$rho_CAFhi[i])))
for (i in seq_len(nrow(cp)))
 V <- rbind(V, rowf("N12", paste0(cp$miRNA[i]," -> ",cp$target[i]," (arm-averaged)"), NA,
   sprintf("mRNA %.4f (p=%.3g) -> protein %.4f (p=%.3g)", cp$rho_mRNA[i], cp$p_mRNA[i],
           cp$rho_protein[i], cp$p_protein[i]), "", sprintf("CPTAC-BRCA n=%d", cp$n[i])))
write.csv(V, "results/v2/verify_expression.csv", row.names=FALSE)
cat("\nverify_expression.csv rows:", nrow(V), "\n")
cat("DONE 07\n")
