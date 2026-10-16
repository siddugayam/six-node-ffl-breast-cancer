## build_baseline_network.R -- a baseline regulatory network on the 587 nodes of the analysed network, with every edge
## among them taken directly from the primary interaction tables rather than from the exported FFL networks
## (see SETTINGS.md):
##   TF -> gene/TF     TRRUST version 2 human (data/db/trrust_human.tsv)
##   TF -> miRNA       TransmiR version 2.0 human (data/db/transmir_hsa.tsv)
##   miRNA -> target   multiMiR validated tables (cache/multimir/validated_all.rds; miRTarBase, TarBase, miRecords;
##                     non-functional miRTarBase records excluded) and, in the second network, TargetScan 8.0 default
##                     predictions (cache/direct/ts/Predicted_Targets_Info.default_predictions.txt, miR_Family_Info.txt)
## miRNA names are canonicalised with scripts/01_network_assembly/01_mirna_canon.R. Writes, for two networks
## (validated; validated plus predicted) and two typings (published; 22 co-regulators as genes), the null-model input
## files network_<name>_<typing>.txt (+ .meta.json), network_arc_counts.csv and analysed_network_overlap.csv.
## Root of the repository: environment variable FFL_ROOT.
## usage: FFL_ROOT=<repository checkout> Rscript build_baseline_network.R
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
ROOT <- Sys.getenv("FFL_ROOT"); stopifnot(nzchar(ROOT))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(ROOT, "scripts/01_network_assembly/01_mirna_canon.R"))
nodes <- fread(file.path(ROOT, "data/network/canonical_nodes.tsv"))
prot <- nodes[type %in% c("TF", "Gene"), name]; mirs <- nodes[type == "miRNA", name]
stopifnot(nrow(nodes) == 587)

## TF -> gene/TF (TRRUST)
tr <- fread(file.path(ROOT, "data/db/trrust_human.tsv"), header = FALSE, col.names = c("source", "target", "mode", "pmid"))
TT <- unique(tr[source %in% prot & target %in% prot & source != target, .(source, target)])
lay <- fread(file.path(ROOT, "data/network/layer_TF_target.tsv"))
chk_tt <- setequal(paste(TT$source, TT$target), unique(lay[source != target, paste(source, target)]))
cat("TRRUST arcs among the nodes:", nrow(TT), "; equal to the project's TRRUST layer:", chk_tt, "\n"); stopifnot(chk_tt)

## TF -> miRNA (TransmiR)
tm <- fread(file.path(ROOT, "data/db/transmir_hsa.tsv"), header = FALSE, quote = "")
TM <- unique(data.table(source = trimws(tm$V1), target = canon_mirna(trimws(tm$V2)))[source %in% prot & target %in% mirs])
lm_ <- fread(file.path(ROOT, "data/network/layer_TF_miRNA.tsv"))
chk_tm <- setequal(paste(TM$source, TM$target), paste(lm_$source, lm_$target))
cat("TransmiR arcs among the nodes:", nrow(TM), "; equal to the project's TransmiR layer:", chk_tm, "\n"); stopifnot(chk_tm)

## miRNA -> target, validated (multiMiR)
V <- as.data.table(readRDS(file.path(ROOT, "cache/multimir/validated_all.rds")))
V <- V[!grepl("non-functional", tolower(trimws(support_type)))]
MV <- unique(data.table(source = canon_mirna(V$mature_mirna_id), target = toupper(V$target_symbol))[source %in% mirs & target %in% prot])
cat("validated miRNA -> target arcs among the nodes:", nrow(MV), "\n")

## miRNA -> target, predicted (TargetScan 8.0 default predictions, human)
P <- fread(file.path(ROOT, "cache/direct/ts/Predicted_Targets_Info.default_predictions.txt"))[`Species ID` == 9606]
F <- fread(file.path(ROOT, "cache/direct/ts/miR_Family_Info.txt"))[`Species ID` == 9606]
F[, mir := canon_mirna(`MiRBase ID`)]
fam <- unique(F[, .(family = `miR family`, mir)])
PF <- merge(unique(P[, .(family = `miR Family`, target = toupper(`Gene Symbol`))]), fam, by = "family", allow.cartesian = TRUE)
cat("TargetScan human families:", uniqueN(P$`miR Family`), "; mapped to a human miRBase entry:",
    sum(unique(P$`miR Family`) %in% fam$family), "\n")
MP <- unique(PF[mir %in% mirs & target %in% prot, .(source = mir, target)])
cat("predicted miRNA -> target arcs among the nodes:", nrow(MP), "; not validated:", nrow(fsetdiff(MP, MV)), "\n")

NETS <- list(validated = unique(rbind(TT, TM, MV)), validated_plus_predicted = unique(rbind(TT, TM, MV, MP)))
T22 <- strsplit("APEX1 BMI1 BRCA1 CREBBP CTNNB1 EP300 EZH2 HDAC1 HDAC2 HDAC3 HDAC4 HDAC9 ILF3 MEN1 MKL1 MTA1 NCOR1 NF1 RB1 SIRT1 SUZ12 VHL", " ")[[1]]
stopifnot(all(T22 %in% nodes[type == "TF", name]))
TN <- c(miRNA = 0L, TF = 1L, Gene = 2L)
V0 <- sort(nodes$name); ty0 <- setNames(nodes$type, nodes$name)[V0]
counts <- list()
for (nm in names(NETS)) for (typing in c("published", "retyped")) {
  ty <- ty0; if (typing == "retyped") ty[T22] <- "Gene"
  A <- NETS[[nm]][order(source, target)]
  cls <- sort(unique(paste0(ty[A$source], "->", ty[A$target])))
  ci <- match(paste0(ty[A$source], "->", ty[A$target]), cls) - 1L
  ix <- setNames(seq_along(V0) - 1L, V0)
  fn <- file.path(HERE, sprintf("network_%s_%s.txt", nm, typing))
  writeLines(c(sprintf("%d %d %d", length(V0), nrow(A), length(cls)), paste(TN[ty], collapse = " "),
               sprintf("%d %d %d", ix[A$source], ix[A$target], ci)), fn)
  write_json(list(nodes = V0, classes = cls), paste0(fn, ".meta.json"), auto_unbox = TRUE)
  recip <- sum(sapply(seq_len(nrow(A)), function(i) ty[A$source[i]] == "TF" && ty[A$target[i]] == "miRNA")
               & paste(A$source, A$target) %in% paste(A$target, A$source))
  tab <- table(factor(paste0(ty[A$source], "->", ty[A$target]), levels = cls))
  cc <- as.list(as.integer(tab)); names(cc) <- names(tab)
  counts[[length(counts) + 1]] <- cbind(data.table(network = nm, typing = typing, arcs = nrow(A), reciprocal_TF_miRNA_pairs = recip),
                                        as.data.table(cc))
}
fwrite(rbindlist(counts, fill = TRUE), file.path(HERE, "network_arc_counts.csv"))
## overlap with the analysed network (Table S2 rows in the analysed network)
S2 <- fread(file.path(ROOT, "supplementary_tables/TableS2_all_interactions.csv"))[in_analysed_network == TRUE]
ov <- rbindlist(lapply(names(NETS), function(nm) { A <- paste(NETS[[nm]]$source, NETS[[nm]]$target)
  rbindlist(lapply(split(S2, S2$edge_type), function(s) data.table(network = nm, edge_type = s$edge_type[1], analysed_arcs = nrow(s),
    in_primary_tables = sum(paste(s$source, s$target) %in% A))))}))
ov <- rbind(ov, rbindlist(lapply(names(NETS), function(nm) data.table(network = nm, edge_type = "all", analysed_arcs = nrow(S2),
  in_primary_tables = sum(paste(S2$source, S2$target) %in% paste(NETS[[nm]]$source, NETS[[nm]]$target)),
  baseline_arcs = nrow(NETS[[nm]])))), fill = TRUE)
fwrite(ov, file.path(HERE, "analysed_network_overlap.csv"))
print(rbindlist(counts, fill = TRUE)); print(ov)
