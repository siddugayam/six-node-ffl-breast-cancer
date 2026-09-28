#!/usr/bin/env Rscript
# =============================================================================
# v7_go_03_summarise.R
# Classify every SIGNIFICANT term (BH < 0.05, Count >= 3) into three themes and
# answer: do the recovered modules annotate to ECM/collagen/stroma biology, to
# cell cycle / proliferation, or to something else?
#
# The keyword classifier is declared here in full and applied unchanged to
# every module.  A term matching neither ECM nor proliferation is "other" and
# is reported verbatim, not discarded.
# =============================================================================
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; V7 <- file.path(REV, "results/v7")
LOGF <- file.path(REV, "logs", "v7_go_03_summarise.log")
lg <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOGF, append = TRUE) }
cat("", file = LOGF)

ECM_RE <- paste0("extracellular matrix|extracellular structure|collagen|ECM-receptor|",
  "basement membrane|elastic fib|connective tissue|fibroblast|",
  "cell-substrate adhes|cell-matrix adhes|focal adhesion|integrin|",
  "wound healing|proteoglycan|matrix organi|external encapsulating|",
  "protein digestion and absorption|relaxin signaling|AGE-RAGE|",
  "epithelial to mesenchymal|epithelial-mesenchymal|myofibroblast|",
  "cartilage|ossification|bone|skeletal system dev|osteoblast|",
  "extracellular structure organi")
PROLIF_RE <- paste0("cell cycle|mitotic|mitosis|cell division|DNA replication|",
  "chromosome segregation|spindle|G1/S|G2/M|S phase|nuclear division|",
  "cell population proliferation|DNA repair|DNA damage|centrosome|kinetochore|",
  "sister chromatid|cyclin|E2F|p53 signaling|telomer")

theme_of <- function(x) {
  e <- grepl(ECM_RE, x, ignore.case = TRUE); p <- grepl(PROLIF_RE, x, ignore.case = TRUE)
  ifelse(e & !p, "ECM_stroma", ifelse(p & !e, "proliferation",
    ifelse(e & p, "both", "other")))
}

ALL <- fread(file.path(V7, "goora_01_ORA_full_term_table.csv"))
ALL[, theme := theme_of(Description)]
fwrite(ALL, file.path(V7, "goora_01_ORA_full_term_table.csv"))

inv <- fread(file.path(V7, "goora_00_module_inventory_size5.csv"))

for (U in c("A_DEtested_20250", "B_network364")) {
  A <- ALL[universe == U]
  lg("\n================ UNIVERSE ", U, " ================")

  ## per module: counts by theme, and the single best term of each theme
  per <- A[sig == TRUE, .(n_sig = .N,
      n_ECM = sum(theme %in% c("ECM_stroma", "both")),
      n_prolif = sum(theme %in% c("proliferation", "both")),
      n_other = sum(theme == "other")),
      by = .(method, partition, tier, module, module_size, n_geneTF, n_tested, ontology)]

  best <- A[order(p.adjust, pvalue)][, .SD[1], by = .(method, partition, module, ontology)]
  setnames(best, c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count","geneID","theme","sig"),
           paste0("top_", c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count","geneID","theme","sig")))

  bestE <- A[theme %in% c("ECM_stroma","both")][order(p.adjust, pvalue)][, .SD[1], by = .(method, partition, module, ontology)]
  setnames(bestE, c("ID","Description","GeneRatio","p.adjust","Count","geneID","sig"),
           paste0("bestECM_", c("ID","Description","GeneRatio","p.adjust","Count","geneID","sig")))
  bestP <- A[theme %in% c("proliferation","both")][order(p.adjust, pvalue)][, .SD[1], by = .(method, partition, module, ontology)]
  setnames(bestP, c("ID","Description","GeneRatio","p.adjust","Count","geneID","sig"),
           paste0("bestPROLIF_", c("ID","Description","GeneRatio","p.adjust","Count","geneID","sig")))

  key <- c("method","partition","module","ontology")
  M <- merge(best[, c(key, grep("^top_", names(best), value = TRUE)), with = FALSE],
             bestE[, c(key, grep("^bestECM_", names(bestE), value = TRUE)), with = FALSE], by = key, all.x = TRUE)
  M <- merge(M, bestP[, c(key, grep("^bestPROLIF_", names(bestP), value = TRUE)), with = FALSE], by = key, all.x = TRUE)
  base <- unique(A[, .(method, partition, tier, module, module_size, n_geneTF, n_tested, ontology)])
  M <- merge(base, M, by = key, all.x = TRUE)
  M <- merge(M, per[, c(key, "n_sig","n_ECM","n_prolif","n_other"), with = FALSE], by = key, all.x = TRUE)
  for (j in c("n_sig","n_ECM","n_prolif","n_other")) set(M, which(is.na(M[[j]])), j, 0L)
  M[, verdict := ifelse(n_sig == 0, "NOTHING at BH<0.05 & >=3 genes",
                 ifelse(n_ECM > 0 & n_prolif == 0, "ECM/stroma only",
                 ifelse(n_prolif > 0 & n_ECM == 0, "proliferation only",
                 ifelse(n_ECM > 0 & n_prolif > 0, "both ECM and proliferation", "other only"))))]
  setorder(M, method, partition, -module_size, ontology)
  fwrite(M, file.path(V7, paste0("goora_02_per_module_verdict_", U, ".csv")))

  lg("\n-- modules with NO significant term at the 3-gene floor --")
  z <- M[n_sig == 0]
  if (nrow(z)) print(z[, .(method, partition, module, module_size, n_tested, ontology, top_Description, top_p.adjust)], nrows = 300) else lg("  (none)")

  lg("\n-- verdict tally --")
  print(M[, .N, by = .(method, ontology, verdict)][order(method, ontology, -N)], nrows = 200)
}

## ---- headline: the strongest ECM term and strongest proliferation term per method
A <- ALL[universe == "A_DEtested_20250" & sig == TRUE]
hl <- rbind(
  A[theme %in% c("ECM_stroma","both")][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "ECM/stroma"],
  A[theme %in% c("proliferation","both")][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "proliferation"],
  A[theme == "other"][order(p.adjust)][, .SD[1], by = .(method, ontology)][, axis := "other"], fill = TRUE)
setorder(hl, method, axis, ontology)
fwrite(hl[, .(method, axis, ontology, partition, module, module_size, n_tested,
              ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)],
       file.path(V7, "goora_03_headline_best_term_per_method_per_axis.csv"))
lg("\n\n================ HEADLINE (universe A) ================")
print(hl[, .(method, axis, ontology, module, Description, GeneRatio, p.adjust, Count)], nrows = 200)
