#!/usr/bin/env Rscript
## 41_build_signature_collection.R
## Assemble every stromal / ECM / wound-response / EMT gene set used in this analysis,
## update deprecated symbols to current HGNC symbols with org.Hs.eg.db, and record how
## many members are measurable in TCGA-BRCA.
suppressPackageStartupMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi); library(msigdbr)
})
BASE <- "/path/to/revision"
LST  <- file.path(BASE, "cache/v5/farmer/lists")
LOG  <- file.path(BASE, "logs/v5/41_build_signature_collection.log")
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG)
set.seed(42)

## ------------------------------------------------------------------ symbol updating ----
sym_ok  <- as.list(org.Hs.egSYMBOL2EG)
alias   <- as.list(org.Hs.egALIAS2EG)
eg2sym  <- as.list(org.Hs.egSYMBOL)
## manual overrides for symbols org.Hs.eg.db no longer carries; each checked against
## NCBI Gene: MGC3047 = Gene ID 84308 (limitrin), status "replaced", current ID 54587 = MXRA8.
MANUAL <- c(MGC3047 = "MXRA8")
update_symbols <- function(g, setname) {
  g <- unique(trimws(g)); g <- g[nzchar(g)]
  out <- character(0); rec <- list()
  for (x in g) {
    if (x %in% names(MANUAL)) { out <- c(out, MANUAL[[x]]); rec[[x]] <- MANUAL[[x]]; next }
    if (x %in% names(sym_ok)) { out <- c(out, x); next }
    if (x %in% names(alias)) {
      e <- alias[[x]]
      if (length(e) == 1L) { y <- eg2sym[[e]]; out <- c(out, y); rec[[x]] <- y; next }
    }
    rec[[x]] <- NA_character_
  }
  drop <- names(rec)[is.na(unlist(rec))]
  if (length(rec)) for (k in names(rec))
    logf(sprintf("   [%s] %s -> %s", setname, k, ifelse(is.na(rec[[k]]), "UNMAPPED (dropped)", rec[[k]])))
  list(genes = unique(out), n_in = length(g), n_remap = sum(!is.na(unlist(rec))), n_drop = length(drop))
}

sets <- list(); meta <- list()
add <- function(name, genes, source, pmid, doi, note, class) {
  u <- update_symbols(genes, name)
  sets[[name]] <<- u$genes
  meta[[length(meta) + 1L]] <<- data.table(signature = name, class = class,
      n_published = u$n_in, n_symbol_remapped = u$n_remap, n_unmapped_dropped = u$n_drop,
      n_current_symbols = length(u$genes), source = source, pmid = pmid, doi = doi, note = note)
  logf(sprintf("%-32s published=%3d  remapped=%2d  dropped=%2d  final=%3d", name,
               u$n_in, u$n_remap, u$n_drop, length(u$genes)))
}

## ---------------------------------------------------- 1. primary-source published lists ----
rdl <- function(f) readLines(file.path(LST, paste0(f, ".txt")))
add("FARMER_STROMAL", rdl("FARMER2009_STROMAL_50"),
    "Farmer et al., Nat Med 2009", "19122658", "10.1038/nm.1908",
    "Supplementary Table 4 col 1: 50-gene DCN-anchored stroma-related metagene", "primary")
add("WEST_DTF_FIBROMATOSIS", rdl("WEST_DTF_FIBROMATOSIS"),
    "West et al. PLoS Biol 2005 / Beck et al. Lab Invest 2008, via Farmer 2009 Suppl Table 5",
    "15869330;18414401", "10.1371/journal.pbio.0030187;10.1038/labinvest.2008.31",
    "desmoid-type fibromatosis (fibromatosis) stromal-response signature", "stromal")
add("WEST_SFT", rdl("WEST_SFT_SIGNATURE"),
    "West et al. PLoS Biol 2005, via Farmer 2009 Suppl Table 5", "15869330",
    "10.1371/journal.pbio.0030187", "solitary fibrous tumour stromal signature", "stromal")
add("FINAK_REF_STROMA_I", rdl("FINAK_REFERENCE_STROMA_I"),
    "Finak et al. Breast Cancer Res 2006, via Farmer 2009 Suppl Table 5", "17054791",
    "10.1186/bcr1608", "reference stromal gene list", "stromal")
add("FINAK_REF_EPITHELIUM_E", rdl("FINAK_REFERENCE_EPITHELIUM_E"),
    "Finak et al. Breast Cancer Res 2006, via Farmer 2009 Suppl Table 5", "17054791",
    "10.1186/bcr1608", "reference epithelial gene list (negative control)", "epithelial")
add("WINSLOW_STROMAL_SIG1", rdl("WINSLOW_STROMAL_SIG1"),
    "Winslow et al. Breast Cancer Res 2015", "25848820", "10.1186/s13058-015-0530-2",
    "Table 2 gene signature 1, collagen/CAF set (HR 1.79 for new tumour event)", "stromal")
add("TRIULZI_ECM", rdl("TRIULZI_ECM_GENELIST"),
    "Triulzi et al. PLoS One 2013", "23441215", "10.1371/journal.pone.0056761",
    "Table S1 a priori ECM gene list used to define the ECM3 subtype", "ECM")
add("FARMER_EMT", rdl("FARMER_EMT_SIGNATURE"),
    "Farmer et al. Nat Med 2009 Suppl Table 5", "19122658", "10.1038/nm.1908",
    "EMT gene set as used in the original Farmer analysis", "EMT")
add("FARMER_TGFBETA_TARGETS", rdl("FARMER_TGFBETA_TARGETS"),
    "Farmer et al. Nat Med 2009 Suppl Table 5", "19122658", "10.1038/nm.1908",
    "TGF-beta target gene set as used in the original Farmer analysis", "TGFbeta")
add("FARMER_WNT", rdl("FARMER_WNT_SIGNATURE"),
    "Farmer et al. Nat Med 2009 Suppl Table 5", "19122658", "10.1038/nm.1908",
    "WNT gene set as used in the original Farmer analysis", "WNT")
add("CHANG_WOUND_UP_VANTVEER", rdl("CHANG_WOUND_ACTIVATED_UP_VANTVEER"),
    "Chang et al. PLoS Biol 2004, Dataset S2 sheet 11 (SAM, van't Veer cohort)", "14737219",
    "10.1371/journal.pbio.0020007",
    "genes UP in wound-activated versus quiescent breast tumours", "wound")
add("CHANG_WOUND_DN_VANTVEER", rdl("CHANG_WOUND_ACTIVATED_DN_VANTVEER"),
    "Chang et al. PLoS Biol 2004, Dataset S2 sheet 11 (SAM, van't Veer cohort)", "14737219",
    "10.1371/journal.pbio.0020007",
    "genes DOWN in wound-activated versus quiescent breast tumours", "wound")
add("CHANG_WOUND_UP_SORLIE", rdl("CHANG_WOUND_ACTIVATED_UP_SORLIE"),
    "Chang et al. PLoS Biol 2004, Dataset S2 sheet 10 (SAM, Sorlie cohort)", "14737219",
    "10.1371/journal.pbio.0020007",
    "genes UP in wound-activated versus quiescent breast tumours", "wound")
add("MAMMOSPHERE", rdl("MAMMOSPHERE_MSP"),
    "Farmer et al. Nat Med 2009 Suppl Table 5", "19122658", "10.1038/nm.1908",
    "mammosphere signature (specificity control)", "other")

## ------------------------------------------------------------------- 2. MSigDB sets ----
md <- as.data.table(msigdbr(species = "Homo sapiens"))
gset <- function(nm) unique(md[gs_name == nm, gene_symbol])
add("CHANG_CSR_LATE_UP", gset("CSR_LATE_UP.V1_UP"),
    "Chang et al. PLoS Biol 2004, via MSigDB C6 CSR_LATE_UP.V1_UP", "14737219",
    "10.1371/journal.pbio.0020007",
    "fibroblast core serum response (wound response), late up-regulated", "wound")
add("CHANG_CSR_EARLY_UP", gset("CSR_EARLY_UP.V1_UP"),
    "Chang et al. PLoS Biol 2004, via MSigDB C6 CSR_EARLY_UP.V1_UP", "14737219",
    "10.1371/journal.pbio.0020007",
    "fibroblast core serum response (wound response), early up-regulated", "wound")
add("FINAK_SDPP", gset("FINAK_BREAST_CANCER_SDPP_SIGNATURE"),
    "Finak et al. Nat Med 2008, via MSigDB", "18438415", "10.1038/nm1764",
    "stroma-derived prognostic predictor", "stromal")
add("FARMER2005_STROMAL_CLUSTER4", gset("FARMER_BREAST_CANCER_CLUSTER_4"),
    "Farmer et al. Oncogene 2005, via MSigDB", "15897907", "10.1038/sj.onc.1208561",
    "cluster 4: stromal genes co-clustering across breast tumours", "stromal")
add("HALLMARK_EMT", gset("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
    "Liberzon et al. Cell Syst 2015 (MSigDB Hallmark)", "26771021",
    "10.1016/j.cels.2015.12.004", "hallmark EMT", "EMT")
add("HALLMARK_TGF_BETA", gset("HALLMARK_TGF_BETA_SIGNALING"),
    "Liberzon et al. Cell Syst 2015 (MSigDB Hallmark)", "26771021",
    "10.1016/j.cels.2015.12.004", "hallmark TGF-beta signalling", "TGFbeta")
add("HALLMARK_ANGIOGENESIS", gset("HALLMARK_ANGIOGENESIS"),
    "Liberzon et al. Cell Syst 2015 (MSigDB Hallmark)", "26771021",
    "10.1016/j.cels.2015.12.004", "hallmark angiogenesis", "other")
add("NABA_CORE_MATRISOME", gset("NABA_CORE_MATRISOME"),
    "Naba et al. Mol Cell Proteomics 2012 (MSigDB)", "22159717", "10.1074/mcp.M111.014647",
    "core matrisome", "ECM")
add("NABA_COLLAGENS", gset("NABA_COLLAGENS"),
    "Naba et al. Mol Cell Proteomics 2012 (MSigDB)", "22159717", "10.1074/mcp.M111.014647",
    "collagen genes", "ECM")

## ------------------------------------------------- 3. our own stromal/CAF estimators ----
gmt  <- readLines(system.file("extdata", "SI_geneset.gmt", package = "estimate"))
esets <- lapply(strsplit(gmt, "\t"), function(x) x[-c(1, 2)])
names(esets) <- sapply(strsplit(gmt, "\t"), `[`, 1)
add("ESTIMATE_STROMAL", esets$StromalSignature,
    "Yoshihara et al. Nat Commun 2013 (estimate package)", "24113773", "10.1038/ncomms3612",
    "ESTIMATE stromal signature", "ours")
add("ESTIMATE_STROMAL_noCOL", setdiff(esets$StromalSignature, grep("^COL", esets$StromalSignature, value = TRUE)),
    "Yoshihara et al. Nat Commun 2013, collagen genes removed here", "24113773",
    "10.1038/ncomms3612", "ESTIMATE stromal signature with all COL* genes removed", "ours")
add("ESTIMATE_IMMUNE", esets$ImmuneSignature,
    "Yoshihara et al. Nat Commun 2013 (estimate package)", "24113773", "10.1038/ncomms3612",
    "ESTIMATE immune signature", "ours")
caff <- file.path(BASE, "cache/celltype/caf_signature_scrna.txt")
stopifnot(file.exists(caff))
add("CAF_scRNA_50", readLines(caff),
    "this study: derived from GSE176078 breast scRNA-seq pseudobulk", "34493872",
    "10.1038/s41588-021-00911-1",
    "top-50 CAF-specific genes, collagen-free and network-node-free by construction", "ours")
add("CAF_MARKER_PANEL", c("FAP","PDGFRA","PDGFRB","THY1","DCN","LUM","FBLN1","VCAN","SULF1",
                          "ISLR","CDH11","MMP2"),
    "this study: canonical CAF markers", NA, NA, "marker panel used in ct_03_stromal_scores.R", "ours")

## ------------------------------------------------------------------- 4. FFL modules ----
mods <- readRDS(file.path(BASE, "data/ffl_module_sets.rds"))
nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
prot <- nodes[type != "miRNA", name]
gene_only <- function(v) intersect(v, prot)      # protein-coding members only
add("FFL_3node", gene_only(mods$N3), "this study", NA, NA,
    "protein-coding members of 3-node FFL cores", "FFL")
add("FFL_higher_order_only", gene_only(mods$HIGHER_ONLY), "this study", NA, NA,
    "protein-coding members entering only at n>=4 (higher-order-only module)", "FFL")
add("FFL_6node_all", gene_only(mods$N6), "this study", NA, NA,
    "protein-coding members of all FFL cores up to n=6", "FFL")
add("FFL_4node", gene_only(mods$N4), "this study", NA, NA, "protein-coding members up to n=4", "FFL")
add("FFL_5node", gene_only(mods$N5), "this study", NA, NA, "protein-coding members up to n=5", "FFL")
add("MS_exemplar_module", gene_only(mods$MS_MODULE), "this study", NA, NA,
    "the manuscript's exemplar six-node module (protein members)", "FFL")
add("COLLAGEN_PAIR", c("COL1A1", "COL3A1"), "this study", NA, NA,
    "the two collagens of the exemplar module", "FFL")

## ------------------------------------------------------------- 5. TCGA-BRCA coverage ----
gex <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
mt <- rbindlist(meta)
mt[, n_in_TCGA := sapply(sets[signature], function(s) length(intersect(s, rownames(gex))))]
mt[, frac_in_TCGA := round(n_in_TCGA / n_current_symbols, 3)]
setcolorder(mt, c("signature","class","n_published","n_symbol_remapped","n_unmapped_dropped",
                  "n_current_symbols","n_in_TCGA","frac_in_TCGA"))
fwrite(mt, file.path(BASE, "results/v5/farmer_signature_inventory.csv"))
saveRDS(sets, file.path(BASE, "cache/v5/farmer_signature_sets.rds"))
logf("\n=== inventory ===")
print(mt[, .(signature, class, n_published, n_current_symbols, n_in_TCGA, frac_in_TCGA)])
logf("WROTE results/v5/farmer_signature_inventory.csv and cache/v5/farmer_signature_sets.rds")
logf("FARMER_STROMAL final members: ", paste(sort(sets$FARMER_STROMAL), collapse = ", "))
