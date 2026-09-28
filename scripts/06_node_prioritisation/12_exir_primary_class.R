#!/usr/bin/env Rscript
## ============================================================================
## 12_exir_primary_class.R
## ExIR does NOT return a partition: the Driver table and the Biomarker table
## each contain EVERY significant DE feature (same 4311 features, two different
## rankings). Only the nonDE-mediator table is disjoint. This script therefore
## derives an explicit, documented mutually-exclusive assignment and adds it to
## results/exir_classification.csv as the column `primary_class`.
##
##   nonDE_mediator            : in the nonDE-mediator table (disjoint by construction)
##   Driver+Biomarker          : DE feature, ExIR P.adj < 0.05 in BOTH tables
##   Driver (Accel/Decel)      : DE feature, ExIR P.adj < 0.05 in the Driver table only
##   Biomarker (Up/Down)       : DE feature, ExIR P.adj < 0.05 in the Biomarker table only
##   DE_unranked               : DE feature, not ExIR-significant in either table
## ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- "/path/to/revision"
LOG  <- file.path(ROOT, "logs", "exir_agent.log")
say <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | 12 | ", paste0(..., collapse = ""))
  cat(msg, "\n", file = LOG, append = TRUE); cat(msg, "\n"); flush.console()
}
say("=== 12 START ===")

cls <- fread(file.path(ROOT, "results", "exir_classification.csv"))

sigD <- cls[class == "Driver"        & exir_significant == TRUE, feature]
sigB <- cls[class == "Biomarker"     & exir_significant == TRUE, feature]
sigM <- cls[class == "nonDE_mediator"& exir_significant == TRUE, feature]
allM <- cls[class == "nonDE_mediator", feature]
allD <- cls[class == "Driver", feature]

say("significant Drivers=", length(sigD), "  significant Biomarkers=", length(sigB),
    "  significant nonDE-mediators=", length(sigM))
say("sigDriver AND sigBiomarker = ", length(intersect(sigD, sigB)),
    " ; sigDriver only = ", length(setdiff(sigD, sigB)),
    " ; sigBiomarker only = ", length(setdiff(sigB, sigD)))
say("Driver table and Biomarker table membership identical? ",
    identical(sort(allD), sort(cls[class == "Biomarker", feature])))
say("Driver/Biomarker set intersect nonDE-mediator set = ", length(intersect(allD, allM)),
    " (must be 0 by construction)")

dsub <- cls[class == "Driver", .(feature, driver_subtype = subtype)]
bsub <- cls[class == "Biomarker", .(feature, biomarker_dir = subtype)]
key  <- unique(cls[, .(feature, node_type, in_network, logFC, DE_adj_P)])
key  <- merge(key, dsub, by = "feature", all.x = TRUE)
key  <- merge(key, bsub, by = "feature", all.x = TRUE)
key[, primary_class := fifelse(feature %in% allM, "nonDE_mediator",
                        fifelse(feature %in% sigD & feature %in% sigB, "Driver+Biomarker",
                        fifelse(feature %in% sigD, "Driver",
                        fifelse(feature %in% sigB, "Biomarker", "DE_unranked"))))]
key[primary_class %in% c("Driver","Driver+Biomarker"),
    primary_class := paste0(primary_class, " (", driver_subtype, ")")]
pc <- key[, .N, by = primary_class][order(-N)]
print(pc)
say("primary_class counts: ", paste(sprintf("%s=%d", pc$primary_class, pc$N), collapse = "; "))

cls <- merge(cls, key[, .(feature, primary_class)], by = "feature", all.x = TRUE)
setorder(cls, class, rank)
fwrite(cls, file.path(ROOT, "results", "exir_classification.csv"))
say("re-wrote results/exir_classification.csv with primary_class, rows=", nrow(cls))
fwrite(key, file.path(ROOT, "results", "exir_primary_class.csv"))
say("wrote results/exir_primary_class.csv rows=", nrow(key))

## the 4 network nodes among the significant drivers
nn <- key[primary_class %like% "Driver" & in_network == TRUE]
say("canonical-network nodes among the ExIR-significant drivers: ", nrow(nn), " -> ",
    paste(sprintf("%s(%s,logFC %.2f)", nn$feature, nn$node_type, nn$logFC), collapse = "; "))
say("=== 12 DONE ===")
