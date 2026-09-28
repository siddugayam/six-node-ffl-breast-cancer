#!/usr/bin/env Rscript
## B) DGIdb druggability of network hubs
suppressPackageStartupMessages({library(data.table)})
REV <- "/path/to/revision"
OUT <- file.path(REV,"results/multiomics")
RAW <- "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data"
msg <- function(...) cat(format(Sys.time(),"%H:%M:%S"),"|",...,"\n")

ix <- fread(file.path(RAW,"interactions.tsv"), quote="", na.strings=c("NULL","NA",""))
dr <- fread(file.path(RAW,"drugs.tsv"), quote="", na.strings=c("NULL","NA",""))
msg("interactions rows:", nrow(ix), " drugs rows:", nrow(dr))
msg("DGIdb source dbs:", length(unique(ix$interaction_source_db_name)))
msg("distinct genes in DGIdb:", uniqueN(ix$gene_name), " distinct drugs:", uniqueN(ix$drug_name))

topo  <- fread(file.path(REV,"results/network_topology_hubs.csv"))
mods  <- readRDS(file.path(REV,"data/ffl_module_sets.rds"))
ess   <- fread(file.path(OUT,"depmap_essentiality_nodes.csv"))
prot_hubs <- topo[is_hub==TRUE & type %in% c("TF","Gene"), name]
msg("protein-coding hubs:", length(prot_hubs))
mirna_hubs <- topo[is_hub==TRUE & type=="miRNA", name]
msg("miRNA hubs (not druggable via DGIdb):", length(mirna_hubs))

named <- c("NFKB1","RELA","SP1","ETS1","COL1A1","COL3A1","VEGFA","CCND2","MYC",
           "E2F1","TP53","EZH2","STAT3","HIF1A","TGFBR2")
query <- union(prot_hubs, named)
msg("genes queried against DGIdb:", length(query))

sub <- ix[gene_name %in% query]
msg("interaction rows matching queried genes:", nrow(sub))

## drug-level approval / anti-neoplastic status from drugs.tsv (union over rows)
drstat <- dr[, .(drug_approved = any(approved %in% TRUE, na.rm=TRUE),
                 drug_antineoplastic = any(anti_neoplastic %in% TRUE, na.rm=TRUE),
                 drug_immunotherapy = any(immunotherapy %in% TRUE, na.rm=TRUE)),
             by = .(drug_name)]
sub <- merge(sub, drstat, by="drug_name", all.x=TRUE)
sub[, approved_flag := (approved %in% TRUE) | (drug_approved %in% TRUE)]
sub[, antineo_flag  := (anti_neoplastic %in% TRUE) | (drug_antineoplastic %in% TRUE)]
sub[, clinical_stage_source := interaction_source_db_name %in%
      c("TdgClinicalTrial","MyCancerGenomeClinicalTrial","ClearityFoundationClinicalTrial",
        "FDA","NCI","JAX-CKB","CIViC","OncoKB","MyCancerGenome","CancerCommons","CGI","TTD")]

## per gene-drug pair, collapse sources / interaction types
pair <- sub[, .(interaction_types = paste(sort(unique(na.omit(interaction_type))), collapse="|"),
                sources = paste(sort(unique(interaction_source_db_name)), collapse="|"),
                n_source_records = .N,
                max_interaction_score = suppressWarnings(max(as.numeric(interaction_score), na.rm=TRUE)),
                approved = any(approved_flag),
                anti_neoplastic = any(antineo_flag),
                immunotherapy = any(drug_immunotherapy %in% TRUE),
                clinical_stage_evidence = any(clinical_stage_source)),
            by = .(gene_name, drug_name)]
pair[!is.finite(max_interaction_score), max_interaction_score := NA_real_]
pair[interaction_types=="", interaction_types := "unspecified"]

## annotate gene-level network + essentiality attributes
pair <- merge(pair, topo[, .(gene_name=name, node_type=type, degree_tot, betweenness,
                             hub_by_degree, hub_by_betweenness, hub_named_in_MS, is_hub)],
              by="gene_name", all.x=TRUE)
pair <- merge(pair, ess[, .(gene_name=gene, mean_chronos_breast, frac_essential_breast,
                            essential_class, selective_breast)],
              by="gene_name", all.x=TRUE)
pair[, in_MS_module := gene_name %in% mods$MS_MODULE]
pair[, named_in_task := gene_name %in% named]
setorder(pair, -approved, -anti_neoplastic, gene_name, drug_name)
fwrite(pair, file.path(OUT,"hub_druggability.csv"))
msg("wrote hub_druggability.csv rows:", nrow(pair))

## gene-level summary
gsum <- pair[, .(n_drugs = uniqueN(drug_name),
                 n_approved_drugs = uniqueN(drug_name[approved]),
                 n_antineoplastic_drugs = uniqueN(drug_name[anti_neoplastic]),
                 n_clinical_stage_drugs = uniqueN(drug_name[clinical_stage_evidence]),
                 interaction_types = paste(sort(unique(unlist(strsplit(interaction_types,"\\|")))), collapse="|"),
                 example_approved = paste(head(sort(unique(drug_name[approved])),6), collapse="; ")),
             by = .(gene_name, node_type, degree_tot, is_hub, in_MS_module, named_in_task,
                    mean_chronos_breast, frac_essential_breast, essential_class)]
## add hubs with zero drugs
nodrug <- setdiff(query, pair$gene_name)
if(length(nodrug)){
  add <- topo[name %in% nodrug, .(gene_name=name, node_type=type, degree_tot, is_hub)]
  add <- merge(add, ess[, .(gene_name=gene, mean_chronos_breast, frac_essential_breast, essential_class)],
               by="gene_name", all.x=TRUE)
  add[, `:=`(in_MS_module = gene_name %in% mods$MS_MODULE, named_in_task = gene_name %in% named,
             n_drugs=0L, n_approved_drugs=0L, n_antineoplastic_drugs=0L, n_clinical_stage_drugs=0L,
             interaction_types="", example_approved="")]
  gsum <- rbind(gsum, add, fill=TRUE)
}
setorder(gsum, -n_approved_drugs, -n_drugs)
fwrite(gsum, file.path(OUT,"hub_druggability_gene_summary.csv"))
msg("wrote hub_druggability_gene_summary.csv rows:", nrow(gsum))
msg("hubs with >=1 approved drug:", gsum[is_hub==TRUE & n_approved_drugs>0, .N], "of", gsum[is_hub==TRUE, .N])
msg("hubs with 0 drug records:", gsum[is_hub==TRUE & n_drugs==0, .N])
print(gsum[n_drugs>0][order(-n_approved_drugs)][1:25,
   .(gene_name,node_type,degree_tot,n_drugs,n_approved_drugs,n_antineoplastic_drugs,
     essential_class)])
cat("\n--- the 15 named genes ---\n")
print(gsum[named_in_task==TRUE][order(-n_approved_drugs),
   .(gene_name,node_type,n_drugs,n_approved_drugs,n_antineoplastic_drugs,n_clinical_stage_drugs,
     mean_chronos_breast, example_approved)])
cat("\n--- MS module genes: all drug records ---\n")
print(pair[in_MS_module==TRUE, .(gene_name, drug_name, interaction_types, approved,
    anti_neoplastic, clinical_stage_evidence, sources)][order(gene_name,-approved)])
msg("DONE")
