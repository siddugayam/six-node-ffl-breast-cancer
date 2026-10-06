#!/usr/bin/env Rscript
# go_06_answer_per_module.R -- consolidate the ORA into the answer table.
suppressPackageStartupMessages({library(data.table)}); options(width = 280)
REV <- "/path/to/revision"; V7 <- file.path(REV, "results/v7")
LOGF <- file.path(REV, "logs", "go_06_answer_per_module.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

inv <- fread(file.path(V7, "goora_00_module_inventory_size5.csv"))
F   <- fread(file.path(V7, "goora_01_ORA_full_term_table.csv"))
CORE_ECM  <- c("GO:0030198","GO:0043062","GO:0045229","GO:0032963","GO:0032964","GO:0030199",
               "GO:0007160","hsa04512","hsa04974","hsa04510")
STRICT_CC <- c("GO:0007049","GO:0000278","GO:1903047","GO:0051301","GO:0006260","GO:0007059",
               "GO:0044843","GO:0000082","hsa04110","hsa03030","hsa04115")

mk <- function(U) {
  A <- F[universe == U]
  s <- A[sig == TRUE, .(n_sig = .N), by = .(method, partition, module)]
  e <- A[sig == TRUE & ID %in% CORE_ECM][order(p.adjust)][, .SD[1],
        by = .(method, partition, module)][, .(method, partition, module,
        coreECM_term = Description, coreECM_ratio = GeneRatio, coreECM_padj = p.adjust,
        coreECM_n = Count, coreECM_genes = geneID)]
  c2 <- A[sig == TRUE & ID %in% STRICT_CC][order(p.adjust)][, .SD[1],
        by = .(method, partition, module)][, .(method, partition, module,
        CC_term = Description, CC_ratio = GeneRatio, CC_padj = p.adjust,
        CC_n = Count, CC_genes = geneID)]
  t1 <- A[sig == TRUE][order(p.adjust)][, .SD[1], by = .(method, partition, module)][,
        .(method, partition, module, top_term = Description, top_ont = ontology,
          top_ratio = GeneRatio, top_padj = p.adjust, top_n = Count, top_genes = geneID)]
  g <- fread(file.path(V7, "goora_01_module_gene_mapping.csv"))
  M <- Reduce(function(a, b) merge(a, b, by = c("method","partition","module"), all.x = TRUE),
       list(inv[, .(method, partition, tier, module, size, n_gene, n_TF, n_miRNA, n_geneTF)],
            g[, .(method, partition, module, n_tested = if (U == "A_DEtested_20250") n_in_universeA else n_in_universeB)],
            s, t1, e, c2))
  M[is.na(n_sig), n_sig := 0L]
  M[, testable := n_tested >= 3]
  M[, ANSWER := fifelse(!testable, "NOT TESTABLE (<3 annotatable genes)",
              fifelse(n_sig == 0, "ANNOTATES TO NOTHING (no term at BH<0.05 with >=3 genes)",
              fifelse(!is.na(coreECM_padj) & is.na(CC_padj), "ECM/collagen (no cell cycle)",
              fifelse(is.na(coreECM_padj) & !is.na(CC_padj), "cell cycle (no core ECM)",
              fifelse(!is.na(coreECM_padj) & !is.na(CC_padj), "both ECM and cell cycle",
                      "SOMETHING ELSE (neither ECM nor cell cycle)")))))]
  M[, universe := U]; setorder(M, method, partition, -size); M[]
}
OUT <- rbind(mk("A_DEtested_20250"), mk("B_network364"))
fwrite(OUT, file.path(V7, "goora_07_ANSWER_per_module.csv"))

for (U in c("A_DEtested_20250", "B_network364")) {
  lg("\n\n########## ANSWER TABLE, universe ", U, " ##########")
  M <- OUT[universe == U]
  print(M[, .(method = substr(method, 1, 13),
              part = sub("^(MCODE_|cm2_|BioNet_|jAM_)", "", partition), mod = module,
              N = size, genes = n_tested, sig = n_sig, ANSWER,
              ECM = substr(coreECM_term, 1, 34), ECMr = coreECM_ratio, ECMq = signif(coreECM_padj, 2),
              CC = substr(CC_term, 1, 26), CCr = CC_ratio, CCq = signif(CC_padj, 2))], nrows = 80)
  lg("\n-- tally --"); print(M[, .N, by = .(method, ANSWER)][order(method, -N)], nrows = 50)
  lg("\n-- modules annotating to NOTHING --")
  print(M[ANSWER %like% "NOTHING", .(method, partition, module, size, n_tested,
      best_term = top_term, note = "no term passes both BH<0.05 and >=3 genes")], nrows = 50)
  lg("\n-- modules NOT TESTABLE --")
  print(M[ANSWER %like% "NOT TESTABLE", .(method, partition, module, size, n_gene, n_TF, n_miRNA, n_tested)], nrows = 30)
}
