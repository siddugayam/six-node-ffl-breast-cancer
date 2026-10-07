#!/usr/bin/env Rscript
# 11_ffl_module_membership.R
# Enumerate FFL node sets as defined in the Methods:
#   3-node : TF -> miRNA, TF -> Gene, miRNA -> Gene            (composite FFL)
#   4-node : 3-node + Gene2 joined to Gene1 by a gene-gene edge
#   5-node : 4-node + miRNA2 joined to miRNA1 by a miRNA-miRNA edge
#   6-node : 5-node + TF2 joined to TF1 by a TF-TF edge
# Produces the node sets needed for the 3-node vs higher-order head-to-head test.
# Output: results/ffl_module_membership.csv

suppressPackageStartupMessages({ library(data.table) })
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs", "11_ffl_module_membership.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                        cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
cat("", file = LOG); logf("START")

nodes <- fread(file.path(BASE, "data/canonical_nodes.tsv"))
ntype <- setNames(nodes$type, nodes$name)
cedge <- fread(file.path(BASE, "data/canonical_edges.tsv"))

mir_t <- unique(cedge[edge_type == "miRNA_target", .(source, target)])
tf_all<- fread(file.path(BASE, "data/layer_TF_target.tsv"))
tf_m  <- unique(fread(file.path(BASE, "data/layer_TF_miRNA.tsv"))[, .(source, target)])
gg    <- fread(file.path(BASE, "data/layer_gene_gene.tsv"))
mm    <- fread(file.path(BASE, "data/layer_miRNA_miRNA.tsv"))

# TF -> Gene subset (target must be a Gene node) and TF -> TF subset
tf_g  <- unique(tf_all[ntype[source] == "TF" & ntype[target] == "Gene", .(source, target)])
tf_tf <- unique(tf_all[ntype[source] == "TF" & ntype[target] == "TF" & source != target,
                       .(source, target)])
# gene-gene: treat as an undirected association (STRING is undirected; TRRUST tier is directed)
ggu <- unique(rbind(gg[, .(a = source, b = target)], gg[, .(a = target, b = source)]))
ggu <- ggu[ntype[a] == "Gene" & ntype[b] == "Gene" & a != b]
# miRNA-miRNA polycistronic pairs: undirected
mmu <- unique(rbind(mm[, .(a = miRNA_1, b = miRNA_2)], mm[, .(a = miRNA_2, b = miRNA_1)]))
mmu <- mmu[a != b]
# TF-TF: undirected for the purpose of "adding a TF-TF interaction pair"
tftfu <- unique(rbind(tf_tf[, .(a = source, b = target)], tf_tf[, .(a = target, b = source)]))

logf("edges available: miRNA->target=", nrow(mir_t), " TF->Gene=", nrow(tf_g),
     " TF->miRNA=", nrow(tf_m), " TF-TF(undir)=", nrow(tftfu),
     " gene-gene(undir)=", nrow(ggu), " miRNA-miRNA(undir)=", nrow(mmu))

## ---------------------------------------------- 3-node composite FFL -------
# TF t -> miRNA m ; TF t -> Gene g ; miRNA m -> Gene g
setnames(tf_m, c("TF", "miR")); setnames(tf_g, c("TF", "Gene")); setnames(mir_t, c("miR", "Gene"))
cand <- merge(tf_m, tf_g, by = "TF", allow.cartesian = TRUE)     # TF, miR, Gene candidates
ffl3 <- merge(cand, mir_t, by = c("miR", "Gene"))                 # keep those with miR->Gene
ffl3 <- unique(ffl3[, .(TF, miR, Gene)])
logf("3-node composite FFLs enumerated (exhaustive): ", nrow(ffl3))
logf("  distinct TFs=", uniqueN(ffl3$TF), " miRNAs=", uniqueN(ffl3$miR),
     " Genes=", uniqueN(ffl3$Gene))
fwrite(ffl3, file.path(BASE, "results/ffl3_composite_enumerated.csv"))

N3 <- sort(unique(c(ffl3$TF, ffl3$miR, ffl3$Gene)))
logf("N3 (nodes in >=1 3-node composite FFL) = ", length(N3))

## ---------------------------------------------- 4/5/6-node extensions ------
# 4-node: attach Gene2 to the Gene1 of a 3-node FFL
ext_G <- ggu[a %in% unique(ffl3$Gene)]
G2 <- sort(unique(ext_G$b))
# 5-node: attach miRNA2 to the miRNA1 of a 3-node FFL
ext_M <- mmu[a %in% unique(ffl3$miR)]
M2 <- sort(unique(ext_M$b))
# 6-node: attach TF2 to the TF1 of a 3-node FFL
ext_T <- tftfu[a %in% unique(ffl3$TF)]
T2 <- sort(unique(ext_T$b))
logf("extension partners: Gene2=", length(G2), " miRNA2=", length(M2), " TF2=", length(T2))

# exact instance counts under the nested definition
n_ext_g <- ext_G[, .N, by = a]; setnames(n_ext_g, c("Gene", "nG2"))
n_ext_m <- ext_M[, .N, by = a]; setnames(n_ext_m, c("miR",  "nM2"))
n_ext_t <- ext_T[, .N, by = a]; setnames(n_ext_t, c("TF",   "nT2"))
f <- merge(merge(merge(copy(ffl3), n_ext_g, by = "Gene", all.x = TRUE),
                 n_ext_m, by = "miR", all.x = TRUE), n_ext_t, by = "TF", all.x = TRUE)
f[is.na(nG2), nG2 := 0]; f[is.na(nM2), nM2 := 0]; f[is.na(nT2), nT2 := 0]
logf("EXACT instance counts under the manuscript's nested definition:")
logf("  3-node = ", nrow(f))
logf("  4-node = ", format(sum(as.numeric(f$nG2)), big.mark = ","))
logf("  5-node = ", format(sum(as.numeric(f$nG2) * as.numeric(f$nM2)), big.mark = ","))
logf("  6-node = ", format(sum(as.numeric(f$nG2) * as.numeric(f$nM2) * as.numeric(f$nT2)),
                           big.mark = ","))

N4 <- sort(unique(c(N3, G2)))
N5 <- sort(unique(c(N4, M2)))
N6 <- sort(unique(c(N5, T2)))
HIGHER_ONLY <- setdiff(N6, N3)
logf("N3=", length(N3), " N4=", length(N4), " N5=", length(N5), " N6=", length(N6),
     " HIGHER_ONLY (N6 \\ N3)=", length(HIGHER_ONLY))
logf("HIGHER_ONLY by type: ", paste(sprintf("%s=%d", names(table(ntype[HIGHER_ONLY])),
                                            table(ntype[HIGHER_ONLY])), collapse = " "))

## ------------------------------------------------------ the exemplar module --
MS_MODULE <- c("COL1A1", "COL3A1", "NFKB1", "RELA", "SP1",
               "hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c", "hsa-let-7b", "hsa-let-7e")

out <- data.table(node = nodes$name, type = nodes$type)
out[, in_3node        := node %in% N3]
out[, in_4node_set    := node %in% N4]
out[, in_5node_set    := node %in% N5]
out[, in_6node_set    := node %in% N6]
out[, higher_order_only := node %in% HIGHER_ONLY]
out[, in_MS_exemplar_module := node %in% MS_MODULE]
fwrite(out, file.path(BASE, "results/ffl_module_membership.csv"))
logf("WROTE results/ffl_module_membership.csv rows=", nrow(out))

saveRDS(list(N3 = N3, N4 = N4, N5 = N5, N6 = N6, HIGHER_ONLY = HIGHER_ONLY,
             MS_MODULE = MS_MODULE, ffl3 = ffl3),
        file.path(BASE, "data/ffl_module_sets.rds"))
logf("DONE")
