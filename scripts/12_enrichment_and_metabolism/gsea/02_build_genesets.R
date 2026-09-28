#!/usr/bin/env Rscript
# v4/02 -- Assemble MSigDB collections + the FFL-class node sets as a custom collection
suppressPackageStartupMessages({ library(msigdbr); library(data.table) })
ROOT <- "/path/to/revision"
RES  <- file.path(ROOT, "results", "v4"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(ROOT, "cache", "v4"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

COLLS <- list(
  list(tag = "H",              coll = "H",  sub = NULL,             label = "Hallmark"),
  list(tag = "C2_CP_KEGG",     coll = "C2", sub = "CP:KEGG_LEGACY", label = "KEGG (legacy)"),
  list(tag = "C2_CP_KEGGMED",  coll = "C2", sub = "CP:KEGG_MEDICUS",label = "KEGG MEDICUS"),
  list(tag = "C2_CP_REACTOME", coll = "C2", sub = "CP:REACTOME",    label = "Reactome"),
  list(tag = "C5_GO_BP",       coll = "C5", sub = "GO:BP",          label = "GO Biological Process"),
  list(tag = "C5_GO_MF",       coll = "C5", sub = "GO:MF",          label = "GO Molecular Function"),
  list(tag = "C5_GO_CC",       coll = "C5", sub = "GO:CC",          label = "GO Cellular Component"),
  list(tag = "C6",             coll = "C6", sub = NULL,             label = "Oncogenic signatures"),
  list(tag = "C7_IMMUNESIGDB", coll = "C7", sub = "IMMUNESIGDB",    label = "ImmuneSigDB"),
  list(tag = "C3_TFT_GTRD",    coll = "C3", sub = "TFT:GTRD",       label = "TF targets (GTRD)"),
  list(tag = "C3_TFT_LEGACY",  coll = "C3", sub = "TFT:TFT_LEGACY", label = "TF targets (legacy motif)"),
  list(tag = "C3_MIR_MIRDB",   coll = "C3", sub = "MIR:MIRDB",      label = "miRNA targets (miRDB)"),
  list(tag = "C3_MIR_LEGACY",  coll = "C3", sub = "MIR:MIR_LEGACY", label = "miRNA targets (legacy seed)")
)

sets <- list(); inv <- list()
for (cc in COLLS) {
  x <- if (is.null(cc$sub)) msigdbr(species = "Homo sapiens", collection = cc$coll)
       else msigdbr(species = "Homo sapiens", collection = cc$coll, subcollection = cc$sub)
  x <- as.data.table(x)[!is.na(gene_symbol) & gene_symbol != ""]
  L <- split(unique(x[, .(gs_name, gene_symbol)])$gene_symbol,
             unique(x[, .(gs_name, gene_symbol)])$gs_name)
  L <- lapply(L, unique)
  sets[[cc$tag]] <- L
  inv[[cc$tag]] <- data.table(collection = cc$tag, label = cc$label,
                              msigdb_collection = cc$coll,
                              msigdb_subcollection = ifelse(is.null(cc$sub), "", cc$sub),
                              n_sets = length(L),
                              median_size = median(lengths(L)),
                              db_version = unique(x$db_version)[1])
  msg(cc$tag, ": ", length(L), " sets, median size ", median(lengths(L)))
}

## ---------------- FFL class node sets (custom collection, analysis C) --------
mn <- fread(file.path(ROOT, "results", "motif_node_sets.tsv"))
msg("motif_node_sets.tsv: ", nrow(mn), " rows, sets: ",
    paste(unique(mn$motif_set), collapse = ", "))
## protein-coding nodes only (Gene + TF); miRNA nodes are not in gene-symbol space
mnp <- mn[type %in% c("Gene", "TF")]
ffl <- split(mnp$node, mnp$motif_set)
ffl <- lapply(ffl, unique)
## add the whole-network node set as a reference "set"
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
ffl[["ALL_NETWORK_NODES"]] <- unique(nodes$name[nodes$type %in% c("Gene", "TF")])
## split by node type within each 3-node class as a finer readout
for (s in c("3-miR", "3-TF", "3-Comp", "4-node", "5-node", "6-node")) {
  ffl[[paste0(s, "__TFonly")]]   <- unique(mn$node[mn$motif_set == s & mn$type == "TF"])
  ffl[[paste0(s, "__Geneonly")]] <- unique(mn$node[mn$motif_set == s & mn$type == "Gene"])
}
ffl <- ffl[lengths(ffl) >= 5]
sets[["FFL_CLASS"]] <- ffl
inv[["FFL_CLASS"]] <- data.table(collection = "FFL_CLASS", label = "FFL motif-class node sets",
                                 msigdb_collection = "custom", msigdb_subcollection = "",
                                 n_sets = length(ffl), median_size = median(lengths(ffl)),
                                 db_version = "canonical network v1 (587 nodes)")
msg("FFL_CLASS: ", length(ffl), " sets -> ",
    paste(names(ffl), lengths(ffl), sep = "=", collapse = "  "))

saveRDS(sets, file.path(CACHE, "genesets.rds"))
invdt <- rbindlist(inv)
fwrite(invdt, file.path(RES, "gsea_geneset_inventory.csv"))
print(invdt)

## record the FFL class set membership explicitly
fwrite(data.table(gs_name = rep(names(ffl), lengths(ffl)), gene = unlist(ffl, use.names = FALSE)),
       file.path(RES, "gsea_ffl_class_setmembers.csv"))
msg("DONE 02")
