## 20_metabolic_genesets.R -------------------------------------------------
## Build the metabolic gene-set space for the metabolic-arms analysis.
## Sources (all fetched to cache/v4/metabolic, provenance recorded):
##   KEGG REST  (hsa metabolism maps 00010-01999) ......... live
##   Reactome   (R-HSA-1430728 Metabolism + descendants) .. live
##   MSigDB Hallmark (msigdbr 26.1.0) ..................... local cache
##   Human-GEM subsystems (SysBioChalmers Human-GEM.yml) .. live
##   Curated canonical single-pathway signatures .......... hard-coded, provenance flagged
##   FFL motif-class node sets (results/motif_node_sets.tsv)
suppressPackageStartupMessages({library(data.table); library(msigdbr)})
ROOT <- "/path/to/revision"
CACHE <- file.path(ROOT,"cache/v4/metabolic"); RES <- file.path(ROOT,"results/v4")
dir.create(CACHE, recursive=TRUE, showWarnings=FALSE); dir.create(RES, recursive=TRUE, showWarnings=FALSE)
set.seed(1)

E <- readRDS(file.path(ROOT,"data/brca_gene_expr.rds"))
stopifnot(!any(duplicated(rownames(E))))
UNIV <- rownames(E)
cat("universe:", length(UNIV), "symbols; looks-like-symbol frac:",
    round(mean(grepl("^[A-Za-z]", UNIV)),4), "\n")

sets <- list(); meta <- list()
add <- function(id, collection, genes, source, note="") {
  g <- sort(unique(intersect(genes, UNIV)))
  if (length(g) < 5) return(invisible(NULL))
  sets[[id]] <<- g
  meta[[length(meta)+1]] <<- data.table(set_id=id, collection=collection,
      n_genes_source=length(unique(genes)), n_genes_in_universe=length(g),
      source=source, note=note)
}

## ---------------- 1. KEGG (live REST) ------------------------------------
pw   <- fread(file.path(CACHE,"kegg_pathway_list_hsa.txt"), header=FALSE, sep="\t",
              col.names=c("pathway","name"))
lnk  <- fread(file.path(CACHE,"kegg_link_pathway_hsa.txt"), header=FALSE, sep="\t",
              col.names=c("pathway","gene"))
gl   <- fread(file.path(CACHE,"kegg_list_hsa.txt"), header=FALSE, sep="\t", fill=TRUE)
setnames(gl, 1:4, c("gene","type","position","description"))
gl[, symbol := sub(";.*$","", sub("^([^;]*).*","\\1", description))]
gl[, symbol := trimws(sapply(strsplit(symbol, ","), `[`, 1))]
kmap <- setNames(gl$symbol, gl$gene)
lnk[, pathway := sub("^path:","",pathway)]
lnk[, symbol := kmap[gene]]
lnk <- lnk[!is.na(symbol) & symbol != ""]
pw[, name_clean := sub(" - Homo sapiens \\(human\\)$","",name)]
pw[, num := as.integer(sub("^hsa","",pathway))]
## metabolism block = 00010-00999 ; global/overview + drug-resistance = 01000-01999
kegg_met <- pw[num >= 10 & num <= 1999]
cat("KEGG maps in hsa00xxx-hsa01xxx range:", nrow(kegg_met), "\n")
## explicit claim-relevant extras outside the metabolic numeric range
extra_ids <- c("hsa04933")   # AGE-RAGE signalling in diabetic complications
for (i in seq_len(nrow(kegg_met))) {
  p <- kegg_met$pathway[i]
  g <- lnk[pathway == p, unique(symbol)]
  nm <- toupper(gsub("[^A-Za-z0-9]+","_", kegg_met$name_clean[i]))
  add(paste0("KEGG|", p, "_", nm), "KEGG_METABOLISM", g,
      "KEGG REST rest.kegg.jp (fetched 2026-09-09)", kegg_met$name_clean[i])
}
for (p in extra_ids) {
  g <- lnk[pathway == p, unique(symbol)]
  nm <- toupper(gsub("[^A-Za-z0-9]+","_", pw[pathway==p, name_clean]))
  add(paste0("KEGG_EXTRA|", p, "_", nm), "KEGG_CLAIM_SPECIFIC", g,
      "KEGG REST rest.kegg.jp (fetched 2026-09-09)", pw[pathway==p, name_clean])
}

## ---------------- 2. Reactome metabolism ---------------------------------
rel <- fread(file.path(CACHE,"reactome_relation.txt"), header=FALSE, sep="\t",
             col.names=c("parent","child"))
rel <- rel[grepl("^R-HSA-", parent) & grepl("^R-HSA-", child)]
rp  <- fread(file.path(CACHE,"reactome_pathways.txt"), header=FALSE, sep="\t", quote="",
             col.names=c("id","name","species"))
rp  <- rp[species == "Homo sapiens"]
ROOT_MET <- "R-HSA-1430728"
desc <- ROOT_MET; frontier <- ROOT_MET
while (length(frontier)) {
  nxt <- rel[parent %in% frontier, unique(child)]
  nxt <- setdiff(nxt, desc)
  desc <- c(desc, nxt); frontier <- nxt
}
cat("Reactome Metabolism subtree (R-HSA-1430728 + descendants):", length(desc), "pathways\n")
n2r <- fread(file.path(CACHE,"reactome_ncbi_all.txt"), header=FALSE, sep="\t", quote="",
             col.names=c("entrez","reactome","url","name","evidence","species"))
n2r <- n2r[species == "Homo sapiens" & reactome %in% desc]
suppressPackageStartupMessages(library(org.Hs.eg.db))
e2s <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys=unique(as.character(n2r$entrez)), keytype="ENTREZID", columns="SYMBOL"))
e2sv <- setNames(e2s$SYMBOL, e2s$ENTREZID)
n2r[, symbol := e2sv[as.character(entrez)]]
n2r <- n2r[!is.na(symbol)]
rn <- setNames(rp$name, rp$id)
for (p in intersect(desc, unique(n2r$reactome))) {
  g <- n2r[reactome == p, unique(symbol)]
  nm <- toupper(gsub("[^A-Za-z0-9]+","_", rn[p]))
  nm <- substr(nm, 1, 90)
  add(paste0("REACTOME|", p, "_", nm), "REACTOME_METABOLISM", g,
      "Reactome NCBI2Reactome_All_Levels + ReactomePathwaysRelation (fetched 2026-09-09)", rn[p])
}

## ---------------- 3. Hallmark metabolic sets ------------------------------
H <- as.data.table(msigdbr(species="Homo sapiens", collection="H"))
hall_met <- c("HALLMARK_GLYCOLYSIS","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_ADIPOGENESIS","HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_HEME_METABOLISM","HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_MTORC1_SIGNALING","HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_HYPOXIA","HALLMARK_PEROXISOME","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
stopifnot(all(hall_met %in% H$gs_name))
for (s in hall_met)
  add(paste0("HALLMARK|", s), "HALLMARK_METABOLIC", H[gs_name == s, unique(gene_symbol)],
      paste0("msigdbr ", as.character(packageVersion("msigdbr"))), s)

## ---------------- 4. Human-GEM (Human1) subsystems ------------------------
yml <- readLines(file.path(CACHE,"Human-GEM.yml"), warn=FALSE)
gtab <- fread(file.path(CACHE,"humangem_genes.tsv"), fill=TRUE)
ens2sym <- gtab[, .(ensg=genes, sym=geneSymbols)]
ens2symv <- setNames(ens2sym$sym, ens2sym$ensg)
i_rule <- grep("^    - gene_reaction_rule: ", yml)
i_sub  <- grep("^    - subsystem:$", yml)
sub_by_rxn <- character(0); genes_by_rxn <- character(0)
## walk: each reaction block has one gene_reaction_rule line then one subsystem: line
## followed by "      - <name>" lines
rule_txt <- sub("^    - gene_reaction_rule: ", "", yml[i_rule])
rule_txt <- gsub('^"|"$', "", rule_txt)
sub_txt <- vapply(i_sub, function(i) {
  j <- i + 1; out <- character(0)
  while (j <= length(yml) && grepl("^      - ", yml[j])) { out <- c(out, sub("^      - ","",yml[j])); j <- j + 1 }
  paste(out, collapse=";")
}, character(1))
stopifnot(length(i_rule) == length(i_sub))
## pair by order: each reaction block has the rule line immediately before subsystem:
off <- i_sub - i_rule
cat("Human-GEM rule->subsystem offsets:", paste(names(table(off)), table(off), sep="x", collapse=" "), "\n")
stopifnot(all(off >= 1 & off <= 2), !is.unsorted(i_rule), !is.unsorted(i_sub))
hg <- data.table(rule=rule_txt, subsystem=sub_txt)
hg[, ensg := gsub("[()]|\\b(and|or)\\b","", rule)]
hgl <- hg[nzchar(trimws(ensg)) & nzchar(subsystem)]
hgl <- hgl[, .(ensg=trimws(unlist(strsplit(ensg,"[[:space:]]+")))), by=subsystem][nzchar(ensg)]
hgl[, symbol := ens2symv[ensg]]
hgl <- hgl[!is.na(symbol) & nzchar(symbol)]
cat("Human-GEM subsystems parsed:", length(unique(hgl$subsystem)), "; gene-subsystem pairs:", nrow(hgl), "\n")
for (s in unique(hgl$subsystem)) {
  g <- hgl[subsystem == s, unique(symbol)]
  nm <- toupper(gsub("[^A-Za-z0-9]+","_", s)); nm <- substr(nm,1,90)
  add(paste0("HUMANGEM|", nm), "HUMANGEM_SUBSYSTEM", g,
      "Human-GEM main branch model/Human-GEM.yml + genes.tsv (fetched 2026-09-09)", s)
}

## ---------------- 5. Curated canonical single-pathway signatures ----------
## Hand-curated canonical enzyme lists; flagged as such so nothing rests on
## curation alone (KEGG/Reactome/Human-GEM equivalents are scored in parallel).
cur <- list(
 CURATED_ONE_CARBON_FOLATE = c("MTHFD1","MTHFD2","MTHFD1L","MTHFD2L","SHMT1","SHMT2",
   "ALDH1L1","ALDH1L2","TYMS","DHFR","DHFR2","MTR","MTRR","MTHFR","GART","ATIC","PAICS",
   "FPGS","GGH","SLC19A1","SLC46A1","FOLR1","FOLR2","FOLR3","MTFMT","AMT","GLDC","GCSH",
   "DLD","MAT2A","MAT2B","AHCY","CBS","CTH","BHMT"),
 CURATED_ANTIFOLATE_DETERMINANTS = c("DHFR","TYMS","SLC19A1","ABCC1","ABCC2","ABCC3",
   "ABCC4","ABCC5","FPGS","GGH"),
 CURATED_GLUTAMINOLYSIS = c("SLC1A5","SLC7A5","SLC3A2","SLC7A11","SLC38A1","SLC38A2",
   "GLS","GLS2","GLUD1","GLUD2","GOT1","GOT2","GPT2","ASNS","CAD","NAGS","GFPT1","GFPT2"),
 CURATED_SERINE_GLYCINE = c("PHGDH","PSAT1","PSPH","SHMT1","SHMT2","GLDC","AMT","GCSH",
   "DLD","SDS","SDSL","SRR","GATM","GAMT","CBS","CTH"),
 CURATED_PENTOSE_PHOSPHATE = c("G6PD","PGLS","PGD","RPE","RPIA","TKT","TKTL1","TKTL2",
   "TALDO1","H6PD","PRPS1","PRPS2","PRPS1L1","RBKS","DERA","IDNK","NADK"),
 CURATED_LIPID_SYNTHESIS = c("ACLY","ACACA","ACACB","FASN","SCD","SCD5","ACSS2","ELOVL1",
   "ELOVL2","ELOVL3","ELOVL4","ELOVL5","ELOVL6","ELOVL7","ME1","GPAM","AGPAT1","AGPAT2",
   "DGAT1","DGAT2","LPIN1","SREBF1","INSIG1","THRSP","ACSL1","ACSL3","ACSL4","MVD",
   "HMGCR","HMGCS1","SQLE","FDFT1","LSS","IDI1"),
 CURATED_GLYCOLYSIS_CORE = c("HK1","HK2","HK3","GCK","GPI","PFKL","PFKM","PFKP","ALDOA",
   "ALDOB","ALDOC","TPI1","GAPDH","PGK1","PGAM1","ENO1","ENO2","ENO3","PKM","PKLR",
   "LDHA","LDHB","SLC2A1","SLC2A3","SLC16A1","SLC16A3","PDK1","PDK3","PFKFB3","PFKFB4"),
 CURATED_TCA_CYCLE = c("CS","ACO1","ACO2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","OGDH",
   "OGDHL","DLST","DLD","SUCLA2","SUCLG1","SUCLG2","SDHA","SDHB","SDHC","SDHD","FH",
   "MDH1","MDH2","PC","PDHA1","PDHB","PDHX","DLAT")
)
for (nm in names(cur))
  add(paste0("CURATED|", nm), "CURATED_SIGNATURE", cur[[nm]],
      "curated canonical enzyme list (this study; literature-standard members)", nm)

## ---------------- 6. FFL motif-class module sets --------------------------
mn <- fread(file.path(ROOT,"results/motif_node_sets.tsv"))
for (s in unique(mn$motif_set)) {
  if (grepl("^exemplar", s)) next
  g <- mn[motif_set == s & type %in% c("Gene","TF"), unique(node)]
  add(paste0("FFLCLASS|", gsub("[^A-Za-z0-9]+","_",s)), "FFL_MODULE", g,
      "results/motif_node_sets.tsv (protein-coding Gene+TF nodes)", s)
}
n3 <- mn[motif_set %in% c("3-miR","3-TF","3-Comp") & type %in% c("Gene","TF"), unique(node)]
ho <- mn[motif_set %in% c("4-node","5-node","6-node") & type %in% c("Gene","TF"), unique(node)]
add("FFLCLASS|3node_union", "FFL_MODULE", n3, "derived union of 3-miR/3-TF/3-Comp", "3-node union")
add("FFLCLASS|higherorder_union", "FFL_MODULE", ho, "derived union of 4/5/6-node", "higher-order union")
add("FFLCLASS|higherorder_not_3node", "FFL_MODULE", setdiff(ho, n3),
    "derived: higher-order nodes not in any 3-node core", "higher-order-only")
add("FFLCLASS|ALL_NETWORK_NODES", "FFL_MODULE",
    fread(file.path(ROOT,"data/canonical_nodes.tsv"))[[1]], "data/canonical_nodes.tsv", "all nodes")

## ---------------- 7. Composition / context control scores -----------------
ctrl <- list(
 CONTROL_EPITHELIAL = c("EPCAM","KRT8","KRT18","KRT19","CDH1","KRT7","ELF3","CLDN4","CLDN3","KRT5"),
 CONTROL_CAF_SMALL  = c("DCN","LUM","FAP","THY1","PDGFRB","POSTN","COL1A2","FBLN1","COL1A1","COL3A1"),
 CONTROL_IMMUNE_SMALL = c("PTPRC","CD3D","CD2","CD53","LCP1","CORO1A","CD8A","CD19","MS4A1"),
 CONTROL_PROLIFERATION = c("MKI67","PCNA","TOP2A","CCNB1","BIRC5","AURKA","BUB1","CDK1","RRM2","TYMS")
)
for (nm in names(ctrl))
  add(paste0("CONTROL|", nm), "CONTEXT_CONTROL", ctrl[[nm]], "curated context marker panel", nm)
cafA <- readLines(file.path(ROOT,"results/v2/cafA_gene_list_v3.txt"))
add("CONTROL|CONTROL_CAF_FULL", "CONTEXT_CONTROL", cafA,
    "results/v2/cafA_gene_list_v3.txt (project CAF panel)", "CAF panel v3")
imm <- readLines(file.path(ROOT,"results/v2/immune_gene_list_v3.txt"))
add("CONTROL|CONTROL_IMMUNE_FULL", "CONTEXT_CONTROL", imm,
    "results/v2/immune_gene_list_v3.txt (project immune panel)", "immune panel v3")

## ---------------- size filter + save --------------------------------------
M <- rbindlist(meta)
keep <- names(sets)[lengths(sets) >= 10 & lengths(sets) <= 500]
## never drop a claim-critical or module set on size
must <- grep("^CURATED\\||^FFLCLASS\\||^CONTROL\\||^KEGG_EXTRA\\|", names(sets), value=TRUE)
## global metabolic indices retained despite size
must <- c(must, grep("R-HSA-1430728_METABOLISM$|hsa01100_METABOLIC_PATHWAYS$|R-HSA-556833_METABOLISM_OF_LIPIDS$", names(sets), value=TRUE))
keep <- union(keep, intersect(must, names(sets)[lengths(sets) >= 5]))
sets_keep <- sets[keep]
M[, kept := set_id %in% keep]
cat("\n== sets built:", length(sets), " kept (10<=n<=500 or claim-critical):", length(sets_keep), "\n")
print(M[, .(built=.N, kept=sum(kept)), by=collection])
saveRDS(sets_keep, file.path(RES, "metabolic_genesets.rds"))
fwrite(M, file.path(RES, "metabolic_geneset_inventory.csv"))
long <- rbindlist(lapply(names(sets_keep), function(s)
  data.table(set_id=s, gene=sets_keep[[s]])))
fwrite(long, file.path(RES, "metabolic_geneset_members.csv"))
cat("saved:", file.path(RES,"metabolic_genesets.rds"), "\n")
cat("DONE 20\n")
