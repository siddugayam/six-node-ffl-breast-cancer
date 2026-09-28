## v7 / BioNet step 6: is the maximum-scoring module a fibroblast-compartment
## module? Uses the project's own Wu et al. breast-atlas pseudobulk.
suppressMessages({library(data.table)})
ROOT <- "/path/to/revision"
OUT <- file.path(ROOT,"results/v7"); P <- function(...) file.path(OUT, paste0("bionet_", ...))
pb <- fread(file.path(ROOT,"results/v3/sc_wu_pseudobulk_celltype_minor_cpm.csv"))
caf <- grep("^CAFs", names(pb), value=TRUE); can <- grep("^Cancer_", names(pb), value=TRUE)
cat("CAF columns:", paste(caf,collapse=", "), "\nCancer columns:", paste(can,collapse=", "), "\n")
pb[, caf_cpm := rowMeans(.SD), .SDcols=caf]
pb[, can_cpm := rowMeans(.SD), .SDcols=can]
pb[, lr := log2((caf_cpm+1)/(can_cpm+1))]
setkey(pb, gene)
cat("\nsanity check (paper: COL1A1 ~146x, COL3A1 ~161x CAF vs carcinoma):\n")
print(pb[c("COL1A1","COL3A1","FN1","SPARC","POSTN","EZH2","MYBL2","CAV1","LPL","PPARG"),
         .(gene, caf_cpm=round(caf_cpm), can_cpm=round(can_cpm), ratio=round((caf_cpm+1)/(can_cpm+1),1), lr=round(lr,2))])

nod <- readRDS(P("graph_and_pvals.rds"))$nod
mem <- fread(P("module_primary_members.csv"))
sets <- readRDS(file.path(ROOT,"cache/v5/farmer_signature_sets.rds"))
bg  <- intersect(nod[is_mir==FALSE]$name, pb$gene)
mg  <- intersect(mem$name, bg)
ho  <- intersect(sets$FFL_higher_order_only, bg)
n3  <- intersect(sets$FFL_3node, bg)
pri <- intersect(fread(file.path(ROOT,"results/v5/tables/Table4_prioritised_30.csv"))$name, bg)

f <- function(lab, g) data.table(set=lab, n=length(g), median_log2_CAF_over_cancer=round(median(pb[g]$lr),3),
        frac_CAF_skewed=round(mean(pb[g]$lr>1),3),
        p_vs_network_background=if(identical(g,bg)) NA_real_ else
          wilcox.test(pb[g]$lr, pb[setdiff(bg,g)]$lr)$p.value)
tab <- rbind(f("network gene/TF background", bg), f("BioNet primary module (58-node)", mg),
             f("paper: FFL 3-node members", n3), f("paper: higher-order-only members", ho),
             f("paper: prioritised 30", pri))
fwrite(tab, P("compartment_skew.csv")); print(tab)

## which module genes are CAF-skewed at all?
md <- pb[mg][order(-lr)][, .(gene, caf_cpm=round(caf_cpm,1), can_cpm=round(can_cpm,1), log2_CAF_over_cancer=round(lr,2))]
fwrite(md, P("primary_module_compartment.csv"))
cat("\nmost CAF-skewed module members:\n"); print(head(md, 10))
cat("\nmost carcinoma-skewed module members:\n"); print(tail(md, 6))
