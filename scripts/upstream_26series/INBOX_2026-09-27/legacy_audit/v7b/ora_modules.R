## GO-BP + KEGG over-representation of the module-detection outputs with and without the 30 legacy
## miRNA-miRNA edges, using the code, universes, thresholds and cached KEGG release of
## scripts/v7/v7_go_01_ora.R, and the core-ECM term list of scripts/v7/v7_go_05_answer.R.
## Modules (size >= 5): MCODE (stored controls, results/v7/mcode_16_controls.csv: 'as_published' and
## 'drop_all_miRNA_miRNA'), igraph Louvain r1.00 and Leiden r1.0 (sandboxes in ../v7), BioNet primary
## and the jActiveModules consensus (sandboxes here; consensus over converged seeds, see
## jam_converged_consensus.R).  Reproduction check: the full-network Leiden, BioNet and jAM modules are
## in the stored inventory, and their stored term tables must be matched.
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi); library(data.table)})
REV  <- "/path/to/revision"; V7 <- file.path(REV, "results/v7")
HERE <- file.path(REV, "INBOX_2026-09-27/legacy_audit/v7b"); CD <- file.path(REV, "INBOX_2026-09-27/legacy_audit/v7")
set.seed(1234)
P_ADJUST <- "BH"; MINGS <- 10; MAXGS <- 500; SIG_CUT <- 0.05; GENE_FLOOR <- 3
CORE_ECM <- c("GO:0030198","GO:0043062","GO:0045229","GO:0032963","GO:0032964","GO:0030199",
              "GO:0007160","hsa04512","hsa04974","hsa04510")
de <- fread(file.path(REV, "results/v2/v2_DE_genes.csv")); uniA_sym <- unique(de$feature)
map <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = uniA_sym, keytype = "SYMBOL", columns = "ENTREZID")))[!is.na(ENTREZID)]
map <- map[!duplicated(SYMBOL)]; UNIV_A <- unique(map$ENTREZID); SYM2EG <- setNames(map$ENTREZID, map$SYMBOL)
nodes <- fread(file.path(V7, "cm2_nodes.tsv")); setnames(nodes, c("name", "type"))
netpc <- nodes[type %in% c("Gene", "TF")]$name
map2 <- as.data.table(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = netpc, keytype = "SYMBOL", columns = "ENTREZID")))[!is.na(ENTREZID)][!duplicated(SYMBOL)]
UNIV_B <- unique(map2$ENTREZID)
SYM2EG <- c(SYM2EG, setNames(map2$ENTREZID, map2$SYMBOL)[setdiff(map2$SYMBOL, names(SYM2EG))])
EG2SYM <- setNames(names(SYM2EG), SYM2EG)
KG <- readRDS(file.path(REV, "cache/v7/kegg_hsa_download.rds"))
K2G <- as.data.frame(KG$KEGGPATHID2EXTID); names(K2G) <- c("term", "gene")
K2N <- as.data.frame(KG$KEGGPATHID2NAME);  names(K2N) <- c("term", "name")
run_one <- function(eg, univ, univ_name) {
  out <- list()
  go <- tryCatch(enrichGO(gene = eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = univ,
                          pAdjustMethod = P_ADJUST, pvalueCutoff = 1, qvalueCutoff = 1,
                          minGSSize = MINGS, maxGSSize = MAXGS, readable = TRUE), error = function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) { d <- as.data.table(as.data.frame(go)); d[, ontology := "GO:BP"]; out[[1]] <- d }
  kk <- tryCatch(enricher(gene = eg, TERM2GENE = K2G, TERM2NAME = K2N, universe = univ, pAdjustMethod = P_ADJUST,
                          pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MINGS, maxGSSize = MAXGS), error = function(e) NULL)
  if (!is.null(kk) && nrow(as.data.frame(kk))) { d <- as.data.table(as.data.frame(kk)); d[, ontology := "KEGG"]; out[[2]] <- d }
  if (!length(out)) return(NULL)
  r <- rbindlist(out, fill = TRUE); r[, universe := univ_name]; r[]
}

## ---------------------------------------------------------------- inventory
inv <- list(); add <- function(net, set, mod, mem) if (length(mem) >= 5)
  inv[[length(inv) + 1L]] <<- data.table(network = net, set = set, module = mod, size = length(mem), members = paste(mem, collapse = ";"))
mc <- fread(file.path(V7, "mcode_16_controls.csv"))
for (v in c("as_published", "drop_all_miRNA_miRNA")) { x <- mc[network == v]
  for (i in seq_len(nrow(x))) add(ifelse(v == "as_published", "full", "nolegacy"), "MCODE", paste0("C", x$rank[i]), strsplit(x$members[i], ";")[[1]]) }
for (mode in c("full", "nolegacy")) {
  g <- readRDS(file.path(CD, paste0("sandbox_", mode), "results/v7/cm2_graph_undirected.rds")); nm <- igraph::V(g)$name
  ig <- fread(file.path(CD, paste0("sandbox_", mode), "results/v7/cm2_igraph_partitions.csv"))
  for (run in c("Louvain_r1.00", "Leiden_mod_r1.0")) { sp <- split(ig$name, ig[[run]])
    for (k in names(sp)) add(mode, run, paste0("M", k), sp[[k]]) }
  bn <- fread(file.path(HERE, paste0("sandbox_", mode), "results/v7/bionet_module_primary_members.csv"))
  add(mode, "BioNet_primary", "primary", bn[[intersect(c("name", "node", "gene"), names(bn))[1]]])
  mods <- fread(file.path(HERE, paste0("sandbox_", mode), "results/v7/jam_modules_all_runs.csv"))
  top <- mods[variant == "primary" & rank == 1 & size >= 100]
  add(mode, "jAM_consensus", "core", Reduce(intersect, lapply(top$members, function(x) strsplit(x, ";", fixed = TRUE)[[1]])))
}
INV <- rbindlist(inv)

## ---------------------------------------------------------------- ORA
res <- list()
for (i in seq_len(nrow(INV))) {
  mem <- strsplit(INV$members[i], ";")[[1]]
  pc <- mem[nodes[match(mem, name)]$type %in% c("Gene", "TF")]
  eg <- unique(na.omit(SYM2EG[pc])); egA <- intersect(eg, UNIV_A); egB <- intersect(eg, UNIV_B)
  for (u in list(list(egA, UNIV_A, "A_DEtested_20250"), list(egB, UNIV_B, "B_network364"))) {
    if (length(u[[1]]) < 3) next
    r <- run_one(u[[1]], u[[2]], u[[3]])
    if (!is.null(r)) { r[, `:=`(network = INV$network[i], set = INV$set[i], module = INV$module[i], size = INV$size[i])]; res[[length(res) + 1L]] <- r }
  }
}
ALL <- rbindlist(res, fill = TRUE); ALL[, sig := p.adjust < SIG_CUT & Count >= GENE_FLOOR]
fwrite(ALL[, .(network, set, module, size, universe, ontology, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID, sig)],
       file.path(HERE, "ora_modules_full_term_table.csv.gz"))
summ <- ALL[ID %in% CORE_ECM, .(modules_with_sig_coreECM = uniqueN(module[sig]), best_coreECM_BH = min(p.adjust),
                               best_coreECM_term = Description[which.min(p.adjust)]), by = .(network, set, universe)]
cnt <- INV[, .(modules_ge5 = .N), by = .(network, set)]
summ <- merge(cnt, summ, by = c("network", "set"), all = TRUE); setorder(summ, set, universe, network)
fwrite(summ, file.path(HERE, "ora_coreECM_summary.csv")); fwrite(INV, file.path(HERE, "ora_module_inventory.csv"))
print(summ, nrows = 60)

## ---------------------------------------------------------------- reproduction check
st <- fread(file.path(V7, "goora_01_ORA_full_term_table.csv"))
chk <- rbindlist(lapply(list(c("Leiden_mod_r1.0", "cm2_Leiden_mod_r1.0"), c("BioNet_primary", "BioNet_primary"),
                             c("jAM_consensus", "jAM_consensus_core")), function(p) {
  a <- ALL[network == "full" & set == p[1]]; b <- st[partition == p[2]]
  if (!nrow(b)) return(data.table(set = p[1], stored_rows = 0L))
  a[, key := paste(universe, ontology, ID)]; b[, key := paste(universe, ontology, ID)]
  data.table(set = p[1], ours_rows = nrow(a), stored_rows = nrow(b),
             sig_terms_ours = sum(a$sig), sig_terms_stored = sum(b$sig),
             best_coreECM_BH_ours = suppressWarnings(min(a[ID %in% CORE_ECM & universe == "B_network364"]$p.adjust)),
             best_coreECM_BH_stored = suppressWarnings(min(b[ID %in% CORE_ECM & universe == "B_network364"]$p.adjust)))
}), fill = TRUE)
fwrite(chk, file.path(HERE, "ora_reproduction_check.csv")); print(chk)
