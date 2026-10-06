#!/usr/bin/env Rscript
# =============================================================================
# go_04_named_terms.R
# The manuscript's claim is specific, so test the specific terms, in every
# module, and report the ACTUAL p and BH-adjusted p whether or not they pass.
# Because go_02_ora.R ran with pvalueCutoff = 1, a term missing from a
# module's table was genuinely not tested (0 module genes in it, or the term
# fell outside minGSSize/maxGSSize) rather than filtered out for being weak.
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; V7 <- file.path(REV, "results/v7")
LOGF <- file.path(REV, "logs", "go_04_named_terms.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

NAMED <- c(
  "GO:0030198" = "extracellular matrix organization",
  "GO:0043062" = "extracellular structure organization",
  "GO:0045229" = "external encapsulating structure organization",
  "GO:0032963" = "collagen metabolic process",
  "GO:0032964" = "collagen biosynthetic process",
  "GO:0030199" = "collagen fibril organization",
  "GO:0001503" = "ossification",
  "GO:0042060" = "wound healing",
  "GO:0007160" = "cell-matrix adhesion",
  "hsa04512"   = "ECM-receptor interaction",
  "hsa04510"   = "Focal adhesion",
  "hsa04974"   = "Protein digestion and absorption",
  "hsa04151"   = "PI3K-Akt signaling pathway",
  # proliferation side
  "GO:0007049" = "cell cycle",
  "GO:0000278" = "mitotic cell cycle",
  "GO:1903047" = "mitotic cell cycle process",
  "GO:0051301" = "cell division",
  "GO:0006260" = "DNA replication",
  "GO:0007059" = "chromosome segregation",
  "GO:0008283" = "cell population proliferation",
  "GO:0008284" = "positive regulation of cell population proliferation",
  "hsa04110"   = "Cell cycle",
  "hsa03030"   = "DNA replication",
  "hsa04115"   = "p53 signaling pathway")

ALL <- fread(file.path(V7, "goora_01_ORA_full_term_table.csv"))
inv <- fread(file.path(V7, "goora_00_module_inventory_size5.csv"))

grid <- CJ(key = paste(inv$method, inv$partition, inv$module, sep = "|"),
           ID = names(NAMED), universe = unique(ALL$universe), unique = TRUE)
ALL[, key := paste(method, partition, module, sep = "|")]
M <- merge(grid, ALL[ID %in% names(NAMED),
      .(key, ID, universe, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID, sig)],
      by = c("key", "ID", "universe"), all.x = TRUE)
M[is.na(Description), Description := NAMED[ID]]
M[is.na(Count), `:=`(Count = 0L, GeneRatio = "0/NA", status = "not tested (0 module genes in term, or term outside minGSSize/maxGSSize)")]
M[is.na(status), status := ifelse(sig, "SIGNIFICANT (BH<0.05, >=3 genes)",
                            ifelse(p.adjust < 0.05 & Count < 3, "BH<0.05 but below 3-gene floor", "not significant"))]
ky <- unique(inv[, .(key = paste(method, partition, module, sep = "|"), method, partition, tier, module, size, n_geneTF)])
M <- merge(ky, M, by = "key")
setorder(M, method, partition, -size, universe, ID)
fwrite(M, file.path(V7, "goora_04_named_terms_every_module.csv"))
lg("rows: ", nrow(M), "  significant: ", sum(M$status == "SIGNIFICANT (BH<0.05, >=3 genes)", na.rm = TRUE))

A <- M[universe == "A_DEtested_20250"]
lg("\n== named terms reaching significance, universe A ==")
print(A[status == "SIGNIFICANT (BH<0.05, >=3 genes)",
        .(method, partition, module, size, ID, Description, GeneRatio, p.adjust, Count, geneID)], nrows = 300)
lg("\n== ECM/collagen named terms: best BH p across all modules, universe A ==")
print(A[grepl("^GO:0030198|^GO:0043062|^GO:0045229|^GO:00329|^GO:0030199|^hsa04512|^hsa04974|^GO:0007160", ID)][
        order(p.adjust)][1:25, .(method, module, size, Description, GeneRatio, pvalue, p.adjust, Count, geneID)])
