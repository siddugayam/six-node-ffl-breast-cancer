#!/usr/bin/env Rscript
# =============================================================================
# go_05_convergence.R
# Cross-method convergence of the ORA result.  A term "converges" if it is
# significant (BH < 0.05, Count >= 3) in at least one module from each of >= 2
# of the four methods.  Reported against universe_A (the pre-specified background)
# and, as the curation-matched control, universe_B.
# =============================================================================
suppressPackageStartupMessages({library(data.table)}); options(width = 280)
REV <- "/path/to/revision"; V7 <- file.path(REV, "results/v7")
LOGF <- file.path(REV, "logs", "go_05_convergence.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

ALL <- fread(file.path(V7, "goora_01_ORA_full_term_table.csv"))

for (U in c("A_DEtested_20250", "B_network364")) {
  S <- ALL[universe == U & sig == TRUE]
  lg("\n\n############ UNIVERSE ", U, " ############")
  conv <- S[, .(n_methods = uniqueN(method), methods = paste(sort(unique(method)), collapse = "+"),
                n_modules = .N, best_padj = min(p.adjust), best_Count = max(Count),
                median_Count = as.numeric(median(Count)),
                best_where = paste0(method[which.min(p.adjust)], ":", partition[which.min(p.adjust)],
                                    ":", module[which.min(p.adjust)]),
                best_GeneRatio = GeneRatio[which.min(p.adjust)],
                best_genes = geneID[which.min(p.adjust)]),
            by = .(ontology, ID, Description, theme)]
  setorder(conv, -n_methods, best_padj)
  fwrite(conv, file.path(V7, paste0("goora_05_convergent_terms_", U, ".csv")))
  lg("terms significant somewhere: ", nrow(conv),
     " | in >=2 methods: ", sum(conv$n_methods >= 2),
     " | in all 4: ", sum(conv$n_methods == 4))
  lg("\n-- terms significant in ALL FOUR methods, by theme --")
  print(conv[n_methods == 4, .N, by = theme])
  lg("\n-- top 30 terms significant in all four methods (by best BH p) --")
  print(conv[n_methods == 4][order(best_padj)][1:30,
        .(ontology, Description, theme, n_modules, best_padj = signif(best_padj, 2),
          best_GeneRatio, best_where)])
  lg("\n-- ECM/stroma terms in all four methods --")
  print(conv[n_methods == 4 & theme %in% c("ECM_stroma", "both")][order(best_padj)][1:20,
        .(ontology, Description, n_modules, best_padj = signif(best_padj, 2), best_GeneRatio, best_where)])
  lg("\n-- proliferation terms in all four methods --")
  pr <- conv[n_methods == 4 & theme %in% c("proliferation", "both")][order(best_padj)]
  if (nrow(pr)) print(pr[1:min(20, nrow(pr)), .(ontology, Description, n_modules,
        best_padj = signif(best_padj, 2), best_GeneRatio, best_where)]) else lg("  (none)")
}

## per-method top term, and the per-method strongest core-ECM and strict cell-cycle term
CORE_ECM <- c("GO:0030198","GO:0043062","GO:0045229","GO:0032963","GO:0032964","GO:0030199",
              "GO:0007160","hsa04512","hsa04974","hsa04510")
STRICT_CC <- c("GO:0007049","GO:0000278","GO:1903047","GO:0051301","GO:0006260","GO:0007059",
               "GO:0044843","GO:0000082","hsa04110","hsa03030","hsa04115")
S <- ALL[universe == "A_DEtested_20250"]
out <- rbind(
  S[sig == TRUE][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "OVERALL TOP TERM"],
  S[ID %in% CORE_ECM][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "BEST CORE-ECM TERM"],
  S[ID %in% STRICT_CC][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "BEST STRICT CELL-CYCLE TERM"],
  fill = TRUE)
setorder(out, method, axis, ontology)
fwrite(out[, .(method, axis, ontology, partition, module, module_size, n_tested, ID, Description,
               GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID, sig)],
       file.path(V7, "goora_06_per_method_headline.csv"))
lg("\n\n############ PER-METHOD HEADLINE (universe A) ############")
print(out[, .(method, axis, ontology, module = paste0(sub("^(MCODE_|cm2_|BioNet_|jAM_)", "", partition), ":", module),
              N = module_size, Description, GeneRatio, padj = signif(p.adjust, 2), Count, sig)], nrows = 60)
