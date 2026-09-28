#!/usr/bin/env Rscript
## 22b) Per-NODE drug-sensitivity association, so that the integrated matrix (Part 1E) can
## carry a real per-node drug column rather than set membership.
## For every protein-coding network node: Spearman(expression in breast lines, drug response)
## across all PRISM compounds and all GDSC2 drugs; BH within gene across compounds.
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"; OUT <- file.path(REV,"results/v3"); CA <- file.path(REV,"cache/v7")
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

E <- fread(file.path(REV,"data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv")); setnames(E,1,"ModelID")
gs <- sub(" \\(\\d+\\)$","", setdiff(names(E),"ModelID")); stopifnot(!any(duplicated(gs)))
mod <- fread(file.path(REV,"data/depmap/Model_24Q4.csv"))
bre <- mod[OncotreeLineage=="Breast", ModelID]
Eb <- E[ModelID %in% bre]; X <- as.matrix(Eb[,-1]); colnames(X) <- gs; rownames(X) <- Eb$ModelID
nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
prot <- intersect(nodes[type %in% c("Gene","TF"), name], gs)
msg("protein-coding network nodes with breast expression:", length(prot), "of", nodes[type %in% c("Gene","TF"), .N])

assoc <- function(W, idcol, meta, label) {
  ids <- intersect(setdiff(names(W), idcol), rownames(X))
  msg(label, ": breast lines shared with expression:", length(ids), "| compounds:", nrow(W))
  D <- as.matrix(W[, ..ids])
  Xi <- X[ids, prot, drop=FALSE]
  ## rank once
  Rx <- apply(Xi, 2, rank)
  out <- rbindlist(lapply(seq_len(nrow(D)), function(i){
    y <- D[i,]; ok <- is.finite(y)
    if (sum(ok) < 15) return(NULL)
    ry <- rank(y[ok]); ryc <- ry - mean(ry)
    rxc <- sweep(Rx[ok, , drop=FALSE], 2, colMeans(Rx[ok, , drop=FALSE]))
    rho <- as.numeric(crossprod(rxc, ryc) / sqrt(colSums(rxc^2) * sum(ryc^2)))
    n <- sum(ok); tt <- rho*sqrt((n-2)/(1-rho^2)); p <- 2*pt(-abs(tt), n-2)
    data.table(compound=W[[idcol]][i], gene=prot, n_lines=n, rho=rho, p=p)
  }))
  out[, resource := label]
  out
}
P <- fread(file.path(CA,"prism_secondary_dose_response.csv"),
           select=c("broad_id","depmap_id","auc","name","moa"))
P <- P[depmap_id %in% bre & is.finite(auc)]
P <- P[, .(auc=median(auc), name=name[1], moa=moa[1]), by=.(broad_id, depmap_id)]
P <- P[broad_id %in% P[, .N, by=broad_id][N>=15, broad_id]]
WP <- dcast(P, broad_id ~ depmap_id, value.var="auc")
metaP <- unique(P[, .(broad_id, name, moa)])
A1 <- assoc(WP, "broad_id", metaP, "PRISM_secondary_AUC")
A1 <- merge(A1, metaP, by.x="compound", by.y="broad_id", all.x=TRUE)

suppressPackageStartupMessages(library(readxl))
G <- as.data.table(read_excel(file.path(CA,"GDSC2_fitted_dose_response.xlsx")))
map <- mod[SangerModelID!="" & !is.na(SangerModelID), .(SANGER_MODEL_ID=SangerModelID, ModelID)]
G <- merge(G, map, by="SANGER_MODEL_ID")
G <- G[ModelID %in% bre & is.finite(LN_IC50)]
G <- G[, .(LN_IC50=median(LN_IC50), name=DRUG_NAME[1], moa=PUTATIVE_TARGET[1]), by=.(DRUG_ID, ModelID)]
G <- G[DRUG_ID %in% G[, .N, by=DRUG_ID][N>=15, DRUG_ID]]
WG <- dcast(G, DRUG_ID ~ ModelID, value.var="LN_IC50")
metaG <- unique(G[, .(DRUG_ID, name, moa)])
A2 <- assoc(WG, "DRUG_ID", metaG, "GDSC2_LN_IC50")
A2 <- merge(A2, metaG, by.x="compound", by.y="DRUG_ID", all.x=TRUE)

A <- rbindlist(list(A1, A2), use.names=TRUE)
A[, fdr := p.adjust(p, "BH"), by=.(resource, gene)]     # BH within gene across compounds
fwrite(A[fdr < 0.25], file.path(OUT,"screens_node_drug_associations_fdr25.csv"))
msg("total gene x compound tests:", nrow(A), "| with FDR<0.05:", sum(A$fdr<0.05))
SUM <- A[, .(n_tests=.N, n_fdr05=sum(fdr<0.05), best_fdr=min(fdr),
             best_rho=rho[which.min(fdr)], best_compound=name[which.min(fdr)],
             best_moa=moa[which.min(fdr)]), by=.(gene, resource)]
SUM <- dcast(SUM, gene ~ resource, value.var=c("n_fdr05","best_fdr","best_rho","best_compound","best_moa"))
fwrite(SUM, file.path(OUT,"screens_node_drug_association_summary.csv"))
msg("wrote per-node drug-association summary:", nrow(SUM), "genes")
msg("\n=== nodes with the most FDR<0.05 drug associations (PRISM) ===")
print(SUM[order(-n_fdr05_PRISM_secondary_AUC)][1:20,
      .(gene, n_fdr05_PRISM_secondary_AUC, best_fdr_PRISM_secondary_AUC, best_rho_PRISM_secondary_AUC,
        best_compound_PRISM_secondary_AUC)], digits=3)
msg("\n=== focus genes ===")
print(SUM[gene %in% c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","MYC","STAT3","HIF1A","MKL1","VEGFA","ESR1")], digits=3)
msg("DONE 22b")
