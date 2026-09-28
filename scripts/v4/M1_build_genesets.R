## ==========================================================================
## M1_build_genesets.R  -- define every gene set used by the multi-cohort
## expression / survival meta-analysis and write them to
## results/v4/multicohort_genesets.csv (+ an RDS for the analysis scripts).
## Every set is derived from a file on disk; nothing is hand-typed except the
## five literature tumour-suppressor anchors, which are named in the analysis plan.
## ==========================================================================
suppressPackageStartupMessages({library(data.table)})
setwd("/path/to/revision")
OUT <- "results/v4"; CA <- "cache/v4/multicohort"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(CA,  showWarnings=FALSE, recursive=TRUE)
S <- list()

## ---- 1. network hubs ------------------------------------------------------
H <- fread("results/network_topology_hubs.csv")
H <- unique(H, by="name")
S$HUB_PROTEIN <- sort(unique(H[is_hub==TRUE & type!="miRNA", name]))
S$HUB_MIRNA   <- sort(unique(H[is_hub==TRUE & type=="miRNA", name]))
cat("hubs: protein-coding", length(S$HUB_PROTEIN), " miRNA", length(S$HUB_MIRNA), "\n")

## ---- 2. all network nodes -------------------------------------------------
N <- fread("data/canonical_nodes.tsv")
cat("node table cols:", paste(names(N), collapse=","), "\n")
tc <- if("type" %in% names(N)) "type" else names(N)[2]
nc <- if("name" %in% names(N)) "name" else names(N)[1]
S$ALL_NETWORK_PROTEIN <- sort(unique(N[[nc]][N[[tc]] != "miRNA"]))
S$ALL_NETWORK_MIRNA   <- sort(unique(N[[nc]][N[[tc]] == "miRNA"]))
cat("network nodes: protein", length(S$ALL_NETWORK_PROTEIN), " miRNA", length(S$ALL_NETWORK_MIRNA), "\n")

## ---- 3. FFL motif-class node sets (protein-coding only) -------------------
M <- fread("results/motif_node_sets.tsv")
pc <- M[type != "miRNA"]
for(k in c("3-miR","3-TF","3-Comp","4-node","5-node","6-node")){
  nm <- paste0("FFL_", gsub("-","_",k))
  S[[nm]] <- sort(unique(pc[motif_set==k, node]))
}
S$FFL_3NODE_UNION      <- sort(unique(pc[motif_set %in% c("3-miR","3-TF","3-Comp"), node]))
S$FFL_HIGHERORDER_UNION<- sort(unique(pc[motif_set %in% c("4-node","5-node","6-node"), node]))
S$FFL_HIGHER_ONLY      <- setdiff(S$FFL_HIGHERORDER_UNION, S$FFL_3NODE_UNION)
cat("FFL 3-node union", length(S$FFL_3NODE_UNION),
    " higher-order union", length(S$FFL_HIGHERORDER_UNION),
    " higher-order-only", length(S$FFL_HIGHER_ONLY), "\n")

## ---- 4. miR-130a target sets ---------------------------------------------
AC <- fread("results/v3/mir130a_anticorrelated_validated_genes.csv")
S$MIR130A_ANTICORR <- sort(unique(AC$gene))
TS <- fread("results/v3/mir130a_target_set.csv")
cat("mir130a_target_set tiers:\n"); print(table(TS$tier))
S$MIR130A_STRONG   <- sort(unique(TS[tier=="STRONG_lowthroughput", gene]))
S$MIR130A_TS_ANCHOR<- c("PTEN","SMAD4","TGFBR2","DICER1","KLF4")
cat("miR-130a sets: anticorr", length(S$MIR130A_ANTICORR),
    " strong", length(S$MIR130A_STRONG),
    " ts_anchor", length(S$MIR130A_TS_ANCHOR), "\n")
## the anchors must be inside the strong set -- assert
cat("TS anchors that are in the STRONG tier:",
    paste(intersect(S$MIR130A_TS_ANCHOR, S$MIR130A_STRONG), collapse=","), "\n")

## ---- 5. miR-29 target / ECM module ---------------------------------------
E <- fread("data/canonical_edges.tsv")
m29 <- E[edge_type=="miRNA_target" & grepl("^hsa-miR-29[abc]$", source), unique(target)]
S$MIR29_TARGETS_NET <- sort(m29)
ecm <- tryCatch({
  suppressPackageStartupMessages(library(msigdbr))
  r <- as.data.table(msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"))
  gcol <- intersect(c("gene_symbol","human_gene_symbol"), names(r))[1]
  scol <- intersect(c("gs_name"), names(r))[1]
  unique(r[[gcol]][r[[scol]]=="REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"])
}, error=function(e){ cat("!! msigdbr failed:", conditionMessage(e), "\n"); character(0) })
cat("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION genes retrieved:", length(ecm), "\n")
S$MIR29_ECM <- sort(intersect(S$MIR29_TARGETS_NET, ecm))
cat("miR-29 network targets", length(S$MIR29_TARGETS_NET),
    " of which ECM-organisation", length(S$MIR29_ECM), ":",
    paste(S$MIR29_ECM, collapse=","), "\n")

## ---- 6. helper gene sets used for covariates / subtyping ------------------
S$CAF_A <- readLines("cache/newcohorts/cafA_signature_genes.txt")
S$EPITHELIAL <- c("EPCAM","KRT8","KRT18","KRT19","CDH1","KRT7","ELF3","CLDN4")
S$IMMUNE     <- c("PTPRC","CD3D","CD2","CD53","LCP1","CORO1A")
S$PROLIF     <- c("MKI67","CCNB1","BIRC5","AURKA","RRM2","TYMS","UBE2C","TOP2A")
## PAM50 gene list is needed only where the cohort has no subtype call
S$PAM50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
 "CDC6","CDCA1","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4",
 "FOXA1","FOXC1","GPR160","GRB7","KIF2C","KNTC2","KRT14","KRT17","KRT5","MAPT","MDM2",
 "MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6",
 "ORC6L","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")

saveRDS(S, file.path(CA,"genesets.rds"))
DT <- rbindlist(lapply(names(S), function(n) data.table(set=n, gene=S[[n]])))
fwrite(DT, file.path(OUT,"multicohort_genesets.csv"))
SUM <- data.table(set=names(S), n=sapply(S, length))
fwrite(SUM, file.path(OUT,"multicohort_geneset_sizes.csv"))
print(SUM)
## symbols we need streamed out of the 30,865-gene SCAN-B matrix
want <- sort(unique(unlist(S)))
want <- want[!grepl("^hsa-", want)]
writeLines(want, file.path(CA,"scanb_wanted_symbols.txt"))
cat("\nunique protein-coding symbols wanted from SCAN-B:", length(want), "\n")
