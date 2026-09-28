#!/usr/bin/env Rscript
## 34_collagenfree_variants.R
## CIRCULARITY GUARD. Several published fibroblast signatures literally contain the
## outcome variable: MCP-counter's Fibroblasts score is the mean of 8 markers of which
## COL1A1, COL3A1, COL6A1 and COL6A2 are four; the two highest-loading signature genes
## of EPIC's CAF reference are COL3A1 and COL1A1; ConsensusTME's Fibroblasts set
## contains COL3A1 and COL14A1; the ESTIMATE stromal signature contains 8 COL genes.
## Adjusting COL1A1 for such a score is partly adjusting COL1A1 for itself.
## Here each of those methods is recomputed with every COL* gene removed, and a
## collagen-content audit is written out.
suppressPackageStartupMessages({library(matrixStats); library(data.table); library(GSVA)})
BASE <- "/path/to/revision"
LOG  <- file.path(BASE, "logs/v3/34_collagenfree.log"); cat("", file = LOG)
logf <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
D <- file.path(BASE, "cache/v3/deconv")
LOG2 <- readRDS(file.path(D, "bulk_log2.rds")); TPM <- readRDS(file.path(D, "bulk_tpm.rds"))
SS <- colnames(LOG2)
FIB <- readRDS(file.path(D, "fibroblast_estimates.rds"))
AUD <- list(); NEW <- list()

## ---------------------------------------------------------------- MCP-counter ----
mg <- read.delim("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt",
                 check.names = FALSE, stringsAsFactors = FALSE)
fb <- mg[mg[["Cell population"]] == "Fibroblasts", "HUGO symbols"]
fbc <- grep("^COL", fb, value = TRUE); fbn <- setdiff(fb, fbc)
logf("MCP-counter Fibroblasts markers (", length(fb), "): ", paste(fb, collapse = ", "))
logf("  collagens (", length(fbc), "): ", paste(fbc, collapse = ", "),
     "  -> collagen-free score built from: ", paste(fbn, collapse = ", "))
AUD$MCPcounter_Fibroblasts <- data.table(method = "MCPcounter", feature = "Fibroblasts",
  n_signature_genes = length(fb), n_collagen = length(fbc), collagens = paste(fbc, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% fb, contains_COL3A1 = "COL3A1" %in% fb)
gi <- intersect(fb, rownames(LOG2))
NEW$MCPcounter_Fibroblasts_recomputed <- colMeans(LOG2[gi, , drop = FALSE])
NEW$MCPcounter_Fibroblasts_nocol      <- colMeans(LOG2[intersect(fbn, rownames(LOG2)), , drop = FALSE])
logf("  reproduction check: Spearman(recomputed, MCPcounter package) = ",
     round(cor(NEW$MCPcounter_Fibroblasts_recomputed, FIB[, "MCPcounter_Fibroblasts"], method = "spearman"), 5))

## ------------------------------------------------------------------------ EPIC ----
suppressPackageStartupMessages(library(EPIC))
sg <- EPIC::TRef$sigGenes; sgc <- grep("^COL", sg, value = TRUE)
logf("\nEPIC TRef signature genes: ", length(sg), " ; collagens: ", paste(sgc, collapse = ", "))
rp <- EPIC::TRef$refProfiles
logf("  CAF reference top 8 signature genes by abundance: ",
     paste(names(sort(rp[sg, "CAFs"], decreasing = TRUE))[1:8], collapse = ", "))
AUD$EPIC_CAFs <- data.table(method = "EPIC", feature = "CAFs", n_signature_genes = length(sg),
  n_collagen = length(sgc), collagens = paste(sgc, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% sg, contains_COL3A1 = "COL3A1" %in% sg)
refNC <- EPIC::TRef; refNC$sigGenes <- setdiff(sg, sgc)
e2 <- EPIC::EPIC(bulk = TPM, reference = refNC, withOtherCells = TRUE)$cellFractions
NEW$EPIC_CAFs_nocol <- e2[SS, "CAFs"]
logf("  Spearman(EPIC CAFs, EPIC CAFs collagen-free) = ",
     round(cor(FIB[, "EPIC_CAFs"], NEW$EPIC_CAFs_nocol, method = "spearman"), 4))

## ---------------------------------------------------------------- ConsensusTME ----
suppressPackageStartupMessages(library(ConsensusTME))
cs <- ConsensusTME::consensusGeneSets[["BRCA"]]
fs <- cs$Fibroblasts; fsc <- grep("^COL", fs, value = TRUE)
logf("\nConsensusTME BRCA Fibroblasts set (", length(fs), "): ", paste(fs, collapse = ", "))
logf("  collagens: ", paste(fsc, collapse = ", "))
AUD$ConsensusTME_Fibroblasts <- data.table(method = "ConsensusTME", feature = "Fibroblasts",
  n_signature_genes = length(fs), n_collagen = length(fsc), collagens = paste(fsc, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% fs, contains_COL3A1 = "COL3A1" %in% fs)
SL <- list(Fibroblasts = intersect(fs, rownames(LOG2)),
           Fibroblasts_nocol = intersect(setdiff(fs, fsc), rownames(LOG2)))
sg2 <- t(GSVA::gsva(ssgseaParam(as.matrix(LOG2), SL, normalize = TRUE)))
NEW$ConsensusTME_Fibroblasts_nocol <- sg2[SS, "Fibroblasts_nocol"]
logf("  Spearman(ConsensusTME Fibroblasts, collagen-free) = ",
     round(cor(FIB[, "ConsensusTME_Fibroblasts"], NEW$ConsensusTME_Fibroblasts_nocol, method = "spearman"), 4))

## -------------------------------------------------------------------- ESTIMATE ----
gmt <- readLines(system.file("extdata", "SI_geneset.gmt", package = "estimate"))
strom <- strsplit(gmt[grep("StromalSignature", gmt)], "\t")[[1]][-c(1, 2)]
sc <- grep("^COL", strom, value = TRUE)
logf("\nESTIMATE stromal signature: ", length(strom), " genes ; collagens (", length(sc), "): ",
     paste(sc, collapse = ", "))
AUD$ESTIMATE_StromalScore <- data.table(method = "ESTIMATE", feature = "StromalScore",
  n_signature_genes = length(strom), n_collagen = length(sc), collagens = paste(sc, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% strom, contains_COL3A1 = "COL3A1" %in% strom)

## ----------------------------------------------------------------------- xCell ----
sgn <- xCell::xCell.data$signatures
fbs <- grep("^Fibroblasts", names(sgn), value = TRUE)
allg <- unique(unlist(lapply(fbs, function(f) GSEABase::geneIds(sgn[[f]]))))
xc <- grep("^COL", allg, value = TRUE)
logf("\nxCell: 9 Fibroblasts signature sets, union of ", length(allg), " genes ; collagens: ",
     ifelse(length(xc), paste(xc, collapse = ", "), "NONE"))
AUD$xCell_Fibroblasts <- data.table(method = "xCell", feature = "Fibroblasts",
  n_signature_genes = length(allg), n_collagen = length(xc), collagens = paste(xc, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% allg, contains_COL3A1 = "COL3A1" %in% allg)

## --------------------------------------------------------------- Wu-derived sigs ---
SIG  <- read.delim(file.path(D, "wu2021_signature_matrix.txt"), row.names = 1, check.names = FALSE)
sc2  <- grep("^COL", rownames(SIG), value = TRUE)
AUD$CBSX_Wu_CAFs <- data.table(method = "CBSX_Wu", feature = "CAFs", n_signature_genes = nrow(SIG),
  n_collagen = length(sc2), collagens = paste(sc2, collapse = ";"),
  contains_COL1A1 = "COL1A1" %in% rownames(SIG), contains_COL3A1 = "COL3A1" %in% rownames(SIG))
AUD$InstaPrism_Wu_CAFs <- data.table(method = "InstaPrism_Wu", feature = "CAFs",
  n_signature_genes = 16732L, n_collagen = 43L, collagens = "all COL* in the 16,732-gene reference",
  contains_COL1A1 = TRUE, contains_COL3A1 = TRUE)

A <- rbindlist(AUD, fill = TRUE)
fwrite(A, file.path(BASE, "results/v3/deconv_collagen_content_audit.csv"))
logf("\n=== COLLAGEN CONTENT OF EACH FIBROBLAST SIGNATURE ===")
for (i in seq_len(nrow(A))) logf(sprintf("%-16s %-12s genes=%5d collagens=%2d  COL1A1=%-5s COL3A1=%-5s  %s",
  A$method[i], A$feature[i], A$n_signature_genes[i], A$n_collagen[i], A$contains_COL1A1[i],
  A$contains_COL3A1[i], A$collagens[i]))

N <- do.call(cbind, NEW); rownames(N) <- SS
saveRDS(N, file.path(D, "fibroblast_estimates_collagenfree.rds"))
FULL <- cbind(FIB, N[SS, , drop = FALSE])
saveRDS(FULL, file.path(D, "fibroblast_estimates_all.rds"))
fwrite(data.table(sample = SS, as.data.table(FULL)),
       file.path(BASE, "results/v3/deconv_fibroblast_estimates.csv"))
logf("\nfibroblast estimates now: ", ncol(FULL), " (", ncol(N), " newly added collagen-free variants)")
gcol <- readRDS(file.path(BASE, "data/brca_gene_expr.rds"))
for (j in colnames(FULL)) logf(sprintf("  %-36s rho with COL1A1 = %+.3f ; with COL3A1 = %+.3f",
  j, cor(FULL[, j], gcol["COL1A1", SS], method = "spearman"),
  cor(FULL[, j], gcol["COL3A1", SS], method = "spearman")))
logf("DONE 34")
