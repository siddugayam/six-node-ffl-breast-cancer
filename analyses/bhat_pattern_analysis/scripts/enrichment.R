## Section E of SETTINGS.md: enrichment of the twelve Bhat networks against the tumour vs normal signature.
##   E1 GSEA (fgseaMultilevel) of the network sets and derived sets, with the settings of scripts/12_enrichment_and_metabolism/gsea/03_run_gsea.R and
##      07_ffl_class_gsea.R; a machinery check reruns the paper's FFL-class collection.
##   E2 the network-background test of 07 (mean |t| against size-matched draws from the network's protein-coding nodes),
##      its signed version, and the transcriptome-wide version.
##   E3 ORA with the settings of scripts/12_enrichment_and_metabolism/over_representation/09_enrichment_ora.R (GO BP/MF/CC, KEGG, Reactome; universes A and B); a
##      machinery check reruns the paper's 3-Comp and 6-node motif sets. KEGG is downloaded by enrichKEGG (DOWNLOADS.tsv).
##      E2 tests the twelve network sets and the derived sets, as 07 tests its class and derived sets.
## usage: Rscript enrichment.R <analysis root> <network deposit folder> <original SIF folder> <output folder> <cores> [variant networks]
##   With a folder of variant networks (A2, A3), the networks are read from it and the machinery checks are not repeated.
suppressPackageStartupMessages({
  library(fgsea); library(data.table); library(BiocParallel); library(clusterProfiler); library(org.Hs.eg.db)
  library(AnnotationDbi); library(parallel)
})
a <- commandArgs(TRUE); ROOT <- a[1]; NET <- a[2]; SIF <- a[3]; OUT <- a[4]; NC <- as.integer(a[5])
NSRC <- if (length(a) >= 6) a[6] else NET; VARIANT <- length(a) >= 6
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(system2("python3", c(file.path(HERE, "bhat_common.py"), "check", ROOT, NET, SIF)) == 0)
LOG <- file(file.path(OUT, "enrichment.log"), "w")
say <- function(...) { m <- paste0(...); cat(m, "\n"); writeLines(m, LOG); flush(LOG) }
RES <- file.path(ROOT, "results", "v4")
NETS <- as.vector(outer(c("miRNA_FFL", "TF_FFL", "composite_FFL"), 3:6, function(k, n) paste0(n, "node_", k)))

## ---------------- sets
nd <- lapply(setNames(NETS, NETS), function(n) fread(file.path(NSRC, paste0(n, "_nodes.tsv")), colClasses = list(character = "node")))
nodes <- fread(file.path(ROOT, "data", "canonical_nodes.tsv"))
netsets <- list()
for (n in NETS) {
  netsets[[n]] <- unique(nd[[n]][type %in% c("Gene", "TF"), node])
  netsets[[paste0(n, "__TFonly")]] <- unique(nd[[n]][type == "TF", node])
  netsets[[paste0(n, "__Geneonly")]] <- unique(nd[[n]][type == "Gene", node])
}
netsets[["ALL_NETWORK_NODES"]] <- unique(nodes$name[nodes$type %in% c("Gene", "TF")])
netsets <- netsets[lengths(netsets) >= 5]
u3 <- unique(unlist(netsets[paste0("3node_", c("miRNA_FFL", "TF_FFL", "composite_FFL"))]))
hi <- unique(unlist(netsets[NETS[!startsWith(NETS, "3node")]]))
derived <- list(`3node_union` = u3)
for (n in NETS[!startsWith(NETS, "3node")]) derived[[paste0(n, "_not_3node")]] <- setdiff(netsets[[n]], u3)
derived[["higherorder_union"]] <- hi
derived[["higherorder_not_3node"]] <- setdiff(hi, u3)
derived[["network_not_in_any_Bhat_network"]] <- setdiff(netsets[["ALL_NETWORK_NODES"]], unique(c(u3, hi)))
derived[["ALL_NETWORK_NODES"]] <- netsets[["ALL_NETWORK_NODES"]]
derived <- derived[lengths(derived) >= 5]
fwrite(data.table(set = rep(c(names(netsets), names(derived)), c(lengths(netsets), lengths(derived))),
                  collection = rep(c("network", "derived"), c(sum(lengths(netsets)), sum(lengths(derived)))),
                  gene = c(unlist(netsets, use.names = FALSE), unlist(derived, use.names = FALSE))), file.path(OUT, "E_set_members.csv"))
say("network sets: ", paste(names(netsets), lengths(netsets), sep = "=", collapse = " "))
say("derived sets: ", paste(names(derived), lengths(derived), sep = "=", collapse = " "))

## ---------------- E1 GSEA
set.seed(20260908)
BPP <- MulticoreParam(workers = 8, RNGseed = 20260908)
read_rank <- function(tag) {                                   # as scripts/12_enrichment_and_metabolism/gsea/03_run_gsea.R
  d <- fread(file.path(RES, paste0("rank_", tag, ".csv")))
  stopifnot(all(c("feature", "stat") %in% names(d)))
  d <- d[!is.na(stat)]; d <- d[!duplicated(feature)]
  sort(setNames(d$stat, d$feature), decreasing = TRUE)
}
st <- read_rank("tumour_vs_normal")
say("ranked list tumour_vs_normal: n=", length(st))
run <- function(P, minS, maxS, coll) {
  r <- as.data.table(suppressWarnings(fgseaMultilevel(P, st, minSize = minS, maxSize = maxS, nPermSimple = 10000L,
                                                      eps = 0, nproc = 0, BPPARAM = BPP)))
  r[, padj := p.adjust(pval, "BH")]
  r[, `:=`(collection = coll, leadingEdge_size = lengths(leadingEdge), leadingEdge = vapply(leadingEdge, paste, "", collapse = ";"),
           minSize_used = minS, maxSize_used = maxS)]
  r
}
G1 <- run(netsets, 10L, 500L, "network")
G2 <- run(derived, 5L, 600L, "derived")
GS <- rbind(G1, G2); setorder(GS, collection, NES)
fwrite(GS, file.path(OUT, "E1_gsea.csv"))
say("E1 GSEA rows: ", nrow(GS))
## machinery check: the paper's FFL-class collection through the same code (main run only)
chk1 <- "CHECK E1: not repeated for variant networks"
if (!VARIANT) {
sets_v4 <- readRDS(file.path(ROOT, "cache", "v4", "genesets.rds"))
Pc <- run(sets_v4$FFL_CLASS, 10L, 500L, "paper_FFL_CLASS")
A <- fread(file.path(RES, "gsea_all_results.csv"))[ranked_list == "tumour_vs_normal" & collection == "FFL_CLASS"]
mm <- merge(Pc, A, by = "pathway", suffixes = c("", ".paper"))
chk1 <- sprintf("CHECK E1: paper FFL-class sets rerun %d, matched %d; max |ES difference| %.3g; max |NES difference| %.3g (limit 0.02); sets with padj < 0.05 here %d, in the paper %d",
                nrow(Pc), nrow(mm), max(abs(mm$ES - mm$ES.paper)), max(abs(mm$NES - mm$NES.paper)), sum(mm$padj < 0.05), sum(mm$padj.paper < 0.05))
say(chk1)
fwrite(mm[, .(pathway, size, size.paper, ES, ES.paper, NES, NES.paper, pval, pval.paper, padj, padj.paper)], file.path(OUT, "E1_check_paper_ffl_class.csv"))
}

## ---------------- E2 network-background test (07, lines 77-95)
set.seed(20260908)
tvn <- fread(file.path(RES, "rank_tumour_vs_normal.csv"))
tvec <- setNames(tvn$stat, tvn$feature)
bg <- intersect(netsets[["ALL_NETWORK_NODES"]], names(tvec))
NB <- 10000L
allsets <- c(netsets[NETS], derived[names(derived) != "ALL_NETWORK_NODES"])        # as 07: class sets and derived sets
nbres <- rbindlist(lapply(names(allsets), function(nm) {
  g <- intersect(allsets[[nm]], names(tvec)); n <- length(g)
  if (n < 5 || n >= length(bg)) return(data.table(set = nm, n_in_bg = n, obs_mean_absT = NA_real_, bg_mean_absT = NA_real_,
                                                  p_vs_network_bg = NA_real_, obs_mean_T = NA_real_, p_signed_vs_network_bg = NA_real_))
  obs <- mean(abs(tvec[g])); obss <- mean(tvec[g])
  draw <- replicate(NB, { s <- sample(bg, n); c(mean(abs(tvec[s])), mean(tvec[s])) })
  data.table(set = nm, n_in_bg = n, obs_mean_absT = obs, bg_mean_absT = mean(draw[1, ]),
             p_vs_network_bg = (1 + sum(draw[1, ] >= obs)) / (NB + 1), obs_mean_T = obss,
             p_signed_vs_network_bg = (1 + sum(abs(draw[2, ]) >= abs(obss))) / (NB + 1))
}))
allg <- names(tvec)
twres <- rbindlist(lapply(names(allsets), function(nm) {
  g <- intersect(allsets[[nm]], allg); n <- length(g)
  if (n < 5) return(NULL)
  obs <- mean(abs(tvec[g])); draw <- replicate(NB, mean(abs(tvec[sample(allg, n)])))
  data.table(set = nm, n = n, transcriptome_mean_absT = mean(draw), p_vs_transcriptome = (1 + sum(draw >= obs)) / (NB + 1))
}))
E2 <- merge(nbres, twres, by = "set", all.x = TRUE)
E2 <- merge(E2, GS[, .(set = pathway, collection, size, NES, pval, padj)], by = "set", all.x = TRUE)
E2[, ord := match(set, names(allsets))]; setorder(E2, ord); E2[, ord := NULL]
fwrite(E2, file.path(OUT, "E2_background.csv"))
say("E2 background rows: ", nrow(E2), " (background ", length(bg), " genes)")

## ---------------- E3 ORA (09_enrichment_ora.R settings)
set.seed(1234)
P_ADJUST <- "BH"; MINGS <- 10; MAXGS <- 500; SIG_CUT <- 0.05; LOW_SUPPORT <- 3
net_pc <- sort(unique(nodes[type %in% c("TF", "Gene"), name]))
de <- fread(file.path(ROOT, "results/BRCA_DEX_genes.csv"))
tcga_tested <- unique(de$feature); tcga_tested <- tcga_tested[grepl("^[A-Za-z]", tcga_tested)]
sym2ent <- function(x) {
  s <- suppressWarnings(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(x), keytype = "SYMBOL", columns = "ENTREZID")))
  s <- s[!is.na(s$ENTREZID), ]; s[!duplicated(s$SYMBOL), ]
}
UNIS <- list(TCGA_tested_plus_network = unique(sym2ent(sort(unique(c(tcga_tested, net_pc))))$ENTREZID),
             network_protein_coding = unique(sym2ent(net_pc)$ENTREZID))
say("universe A entrez ", length(UNIS[[1]]), "; universe B entrez ", length(UNIS[[2]]))
## KEGG: download once here (enrichKEGG caches it for the session); record it
kt <- Sys.time()
invisible(clusterProfiler:::download_KEGG("hsa", "KEGG", "kegg"))
ke <- get("KEGGPATHID2EXTID", envir = clusterProfiler:::get_KEGG_Env()); kn <- get("KEGGPATHID2NAME", envir = clusterProfiler:::get_KEGG_Env())
rel <- tryCatch(readLines(url("https://rest.kegg.jp/info/pathway"), warn = FALSE), error = function(e) paste("info not available:", conditionMessage(e)))
writeLines(rel, file.path(OUT, "E3_kegg_release.txt"))
dl <- data.table(time = format(kt, "%Y-%m-%d %H:%M:%S %z"), analysis = "E3 ORA",
                 what = c("KEGG pathway-gene links for hsa (clusterProfiler download.KEGG.Path)", "KEGG pathway names for hsa (same call)",
                          "KEGG pathway release information"),
                 source = c("https://rest.kegg.jp/link/hsa/pathway", "https://rest.kegg.jp/list/pathway/hsa", "https://rest.kegg.jp/info/pathway"),
                 size = c(paste(nrow(ke), "rows"), paste(nrow(kn), "rows"), paste(length(rel), "lines")),
                 kept_as = c("in memory only", "in memory only", "enrichment/E3_kegg_release.txt"))
fwrite(dl, file.path(OUT, "E3_downloads.tsv"), sep = "\t")
say("KEGG: ", nrow(ke), " pathway-gene links, ", nrow(kn), " pathways; release: ", paste(grep("Release", rel, value = TRUE), collapse = " "))
ora_one <- function(sn, syms) {
  m <- sym2ent(unique(syms)); ent <- unique(m$ENTREZID); out <- list()
  for (un in names(UNIS)) {
    U <- UNIS[[un]]; ent_u <- ent[ent %in% U]
    if (length(ent_u) < 2) next
    for (ont in c("BP", "MF", "CC")) {
      e <- tryCatch(enrichGO(gene = ent_u, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ont, universe = U, pAdjustMethod = P_ADJUST,
                             pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MINGS, maxGSSize = MAXGS, readable = TRUE), error = function(x) NULL)
      if (!is.null(e) && nrow(as.data.frame(e))) out[[length(out) + 1]] <- data.table(set = sn, universe = un, ontology = paste0("GO_", ont), as.data.frame(e))
    }
    ek <- tryCatch(enrichKEGG(gene = ent_u, organism = "hsa", keyType = "kegg", universe = U, pAdjustMethod = P_ADJUST,
                              pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MINGS, maxGSSize = MAXGS), error = function(x) NULL)
    if (!is.null(ek)) ek <- tryCatch(setReadable(ek, org.Hs.eg.db, "ENTREZID"), error = function(x) ek)
    if (!is.null(ek) && nrow(as.data.frame(ek))) out[[length(out) + 1]] <- data.table(set = sn, universe = un, ontology = "KEGG", as.data.frame(ek))
    er <- tryCatch(ReactomePA::enrichPathway(gene = ent_u, organism = "human", universe = U, pAdjustMethod = P_ADJUST, pvalueCutoff = 1,
                                             qvalueCutoff = 1, minGSSize = MINGS, maxGSSize = MAXGS, readable = TRUE), error = function(x) NULL)
    if (!is.null(er) && nrow(as.data.frame(er))) out[[length(out) + 1]] <- data.table(set = sn, universe = un, ontology = "Reactome", as.data.frame(er))
  }
  list(rows = rbindlist(out, fill = TRUE), input = data.table(set = sn, n_input_symbols = length(unique(syms)), n_mapped_entrez = length(ent)))
}
mn <- fread(file.path(ROOT, "results/motif_node_sets.tsv"))
jobs <- lapply(setNames(NETS, NETS), function(n) netsets[[n]])
jobs <- jobs[!vapply(jobs, is.null, logical(1))]
if (!VARIANT) jobs <- c(jobs, list(`paper:3-Comp` = mn[motif_set == "3-Comp" & type %in% c("TF", "Gene"), node],
                                   `paper:6-node` = mn[motif_set == "6-node" & type %in% c("TF", "Gene"), node]))
t0 <- Sys.time()
OR <- lapply(names(jobs), function(sn) ora_one(sn, jobs[[sn]]))      # sequential: the annotation databases are SQLite
ALL <- rbindlist(lapply(OR, `[[`, "rows"), fill = TRUE)
say("E3 ORA done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min; rows ", nrow(ALL))
keep <- c("set", "universe", "ontology", "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
ALL <- ALL[, intersect(keep, names(ALL)), with = FALSE]
ALL[, significant_BH_0.05 := p.adjust < SIG_CUT]; ALL[, low_support_lt3_genes := Count < LOW_SUPPORT]
fwrite(ALL[!startsWith(set, "paper:")], file.path(OUT, "E3_ora_full.csv"))
fwrite(rbindlist(lapply(OR, `[[`, "input")), file.path(OUT, "E3_input_sizes.csv"))
SUMM <- ALL[!startsWith(set, "paper:"), .(n_tested = .N, n_sig_terms = sum(significant_BH_0.05), n_sig_lt3_genes = sum(significant_BH_0.05 & Count < LOW_SUPPORT),
                                          min_padj = min(p.adjust)), by = .(set, universe, ontology)]
fwrite(SUMM, file.path(OUT, "E3_ora_summary.csv"))
## machinery check against the paper's stored ORA
P <- fread(file.path(ROOT, "results/enrichment_all_motifs_FULL.csv"))
chk <- c()
if (!VARIANT) for (ps in c("3-Comp", "6-node")) for (ont in c("GO_BP", "GO_MF", "GO_CC", "Reactome", "KEGG")) for (un in names(UNIS)) {
  me <- ALL[set == paste0("paper:", ps) & ontology == ont & universe == un]; pp <- P[motif_set == ps & ontology == ont & universe == un]
  m2 <- merge(me, pp, by = "ID", suffixes = c("", ".paper"))
  chk <- c(chk, sprintf("CHECK E3 %s %s %s: terms here %d, in the paper %d, shared %d; max |p difference| %.3g; significant here %d, in the paper %d",
                        ps, ont, un, nrow(me), nrow(pp), nrow(m2), if (nrow(m2)) max(abs(m2$pvalue - m2$pvalue.paper)) else NA,
                        sum(me$p.adjust < 0.05), sum(pp$p.adjust < 0.05)))
}
writeLines(c(chk1, chk), file.path(OUT, "E_checks.txt")); for (x in chk) say(x)
close(LOG)
