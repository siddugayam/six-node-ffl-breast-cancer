# Complete interaction and supplementary tables.
suppressMessages({library(data.table)})
REV <- "/path/to/revision"; TB <- file.path(REV,"results/v5/tables")
dir.create(TB, showWarnings=FALSE, recursive=TRUE); lg <- function(...) cat(...,"\n")

nodes <- fread(file.path(REV,"data/canonical_nodes.tsv"))
edges <- fread(file.path(REV,"data/canonical_edges.tsv"))
ev    <- fread(file.path(REV,"data/edge_evidence_tier.tsv"))
setnames(ev, grep("^tier$|tier", names(ev), value=TRUE)[1], "tier")
ec    <- unique(fread(file.path(REV,"results/edge_correlation.csv")))
trr   <- fread(file.path(REV,"data/layer_TF_target.tsv"))
tmr   <- fread(file.path(REV,"data/layer_TF_miRNA.tsv"))

# ---- S1: complete interaction table ----------------------------------------------------
S1 <- copy(edges)
S1 <- merge(S1, ev[, .(source, target, tier, validated=if("validated" %in% names(ev)) validated else NA,
                       databases=if("databases" %in% names(ev)) databases else NA)],
            by=c("source","target"), all.x=TRUE)
cn <- intersect(c("source","target","rho","p","fdr","concordant","n_samples"), names(ec))
S1 <- merge(S1, ec[, ..cn], by=c("source","target"), all.x=TRUE)
S1 <- merge(S1, trr[, .(source, target, trrust_mode=mode)], by=c("source","target"), all.x=TRUE)
S1 <- merge(S1, tmr[, .(source, target, transmir_mode=mode)], by=c("source","target"), all.x=TRUE)
S1[, sign_source := fifelse(!is.na(trrust_mode) & trrust_mode %in% c("Activation","Repression"), "TRRUST",
                     fifelse(!is.na(transmir_mode) & transmir_mode %in% c("Activation","Repression"), "TransmiR",
                      fifelse(edge_type=="miRNA_target", "mechanistic (miRNA repression)", "assumed / unannotated")))]
setcolorder(S1, c("source","target","edge_type","source_type","target_type","sign","sign_source",
                  "trrust_mode","transmir_mode","tier"))
fwrite(S1, file.path(TB,"TableS1_all_interactions.csv")); lg("S1 all interactions:", nrow(S1))

# ---- S2: every 3-node FFL core ----------------------------------------------------------
S2 <- fread(file.path(REV,"results/v2/ffl_cores_coherence_corrected.csv"))
for (arm in list(c("regulator","intermediate","RM"), c("intermediate","target","MT"), c("regulator","target","RT"))) {
  e <- ev[, .(source, target, tier)]; setnames(e, c(arm[1], arm[2], paste0("tier_",arm[3])))
  S2 <- merge(S2, e, by=c(arm[1],arm[2]), all.x=TRUE)
}
fwrite(S2, file.path(TB,"TableS2_ffl_cores.csv")); lg("S2 FFL cores:", nrow(S2))

# ---- T2: network composition -------------------------------------------------------------
T2 <- edges[, .(edges=.N), by=.(edge_type, source_type, target_type)][order(-edges)]
T2 <- rbind(T2, data.table(edge_type="TOTAL", source_type="", target_type="", edges=nrow(edges)))
fwrite(T2, file.path(TB,"Table2_network_composition.csv"))
tierT <- ev[, .(edges=.N), by=tier][order(-edges)]
tierT[, pct := round(100*edges/sum(edges),1)]
fwrite(tierT, file.path(TB,"Table2b_evidence_tiers.csv")); lg("T2 written"); print(T2); print(tierT)

# ---- T2c: sign provenance ----------------------------------------------------------------
sg <- S1[, .(edges=.N), by=sign_source][order(-edges)]
sg[, pct := round(100*edges/sum(edges),1)]
fwrite(sg, file.path(TB,"Table2c_sign_provenance.csv")); print(sg)

# ---- T3: FFL census ------------------------------------------------------------------------
T3 <- data.table(n_nodes=3:7,
  ffl_modules=c("6,037","~1.2e5","~1.5e6","~1.9e7","~2.7e8"),
  method=c("exhaustive", rep("RAND-ESU", 4)),
  max_edge_classes=c(3,5,6,6,6),
  note=c("1,434 composite / 206 miRNA-FFL / 9 TF-FFL","","replicate SD 2.1e4",
         "replicate SD 5.5e4","class diversity saturated"))
fwrite(T3, file.path(TB,"Table3_ffl_census.csv"))

# ---- T5: prioritisation --------------------------------------------------------------------
file.copy(file.path(REV,"results/v5/node_prioritisation_full.csv"),
          file.path(TB,"TableS5_node_prioritisation_full.csv"), overwrite=TRUE)
file.copy(file.path(REV,"results/v5/node_compendium_table.csv"),
          file.path(TB,"Table5_prioritised_30.csv"), overwrite=TRUE)
file.copy(file.path(REV,"results/v5/node_top_partners.csv"),
          file.path(TB,"TableS3_prioritised_node_partners.csv"), overwrite=TRUE)
lg("\nwrote", length(list.files(TB)), "tables to results/v5/tables/")
print(data.table(file=list.files(TB), rows=sapply(list.files(TB, full.names=TRUE),
      function(f) tryCatch(nrow(fread(f)), error=function(e) NA))))
