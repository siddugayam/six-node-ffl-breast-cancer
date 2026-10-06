# Reactive stroma and neoadjuvant

## Scripts in this folder

| Script | What it does |
|---|---|
| `40_farmer_extract_sources.py` | Extract published stromal / ECM / wound-response gene lists from PRIMARY sources that were downloaded to cache/v5/farmer/. |
| `40b_chang_wound_lists.py` | Extract the tumour-level wound-response gene lists of Chang et al. 2004 (PLoS Biol, PMID 14737219) from Dataset S2 sheets 10 (Sorlie cohort) and 11 (van't Veer cohort): the SAM significant genes separating 'wound-activated' from … |
| `41_build_signature_collection.R` | Assemble every stromal / ECM / wound-response / EMT gene set used in this analysis, update deprecated symbols to current HGNC symbols with org.Hs.eg.db, and record how many members are measurable in TCGA-BRCA. |
| `42_farmer_scores_tcga.R` | Per-sample Farmer stromal score (and every comparison signature) in TCGA-BRCA primary tumours by ssGSEA (Barbie), GSVA (Gaussian kcdf) and the mean-of-z score that Farmer et al. actually used ("the average of the 50 genes"), then … |
| `43_farmer_module_overlap.R` | THE KEY TEST How much of the higher-order FFL module IS the Farmer reactive-stroma programme? |
| `44_farmer_mediation.R` | Repeat the causal-mediation analysis of the TF -> collagen and miR-29 -\| collagen edges using the *independently published, clinically validated* Farmer stromal signature as the mediator, alongside our own CAF estimates. |
| `45_farmer_survival.R` | Do the FFL module score and the Farmer stromal score predict survival in TCGA-BRCA and METABRIC, and does the FFL module carry any prognostic information that the published stromal programme does not? |
| `46_helpers.R` | Shared GEO series-matrix parsing helpers (extracted from 46_neoadjuvant_response.R) |
| `46_neoadjuvant_response.R` | Does the Farmer stromal score -- and does the FFL module score -- predict pathological complete response to neoadjuvant chemotherapy? |
| `47_neoadj_survival_and_agreement.R` | (a) GSE25066 distant relapse-free survival (the only neoadjuvant cohort with outcome time) (b) duplicate-patient detection across the GPL96 neoadjuvant series (they come from overlapping MDACC/USO/LAB FNA collections, so the … |
| `48_neoadj_stratified.R` | Farmer's original result was obtained in ER-NEGATIVE tumours receiving FEC. |
| `49_neoadj_meanz_sensitivity.R` | Sensitivity analysis: repeat the pCR association using the *mean-of-z* score that Farmer et al. actually used ("the average of the 50 genes"), instead of ssGSEA, and check that the small five-gene exemplar module behaves the same … |
| `50_stroma_axis_and_figures.R` | (a) how much of each FFL module belongs to the UNION of published stromal/ECM programmes (b) a single "reactive-stroma axis" (PC1 of eight independent stromal scores) and the fraction of each module score's variance that this … |
