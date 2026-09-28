# Script run order (reconstructed)

This note was reconstructed from script file-name ordering and each script's own header comment.
It was NOT found as a pipeline document in the project: no Makefile, run-order note or pipeline
script exists. Treat it as a guide, not an authoritative record. Scripts are grouped by the
directory they live in; within a directory the numeric prefix gives the order they were run.

442 scripts across 6 directories.


## scripts/

| script | purpose (from its header) |
|---|---|
| `00_mirna_canon.R` | Shared helper: canonicalise a miRBase name (precursor or mature, any case) |
| `01_build_canonical_network.py` | 01_build_canonical_network.py |
| `01_tcga_brca_prep_de.R` | TCGA-BRCA expression matrices + differential expression (tumour vs normal) |
| `01b_fix_entrez_to_symbol.R` | Repair: the pan-cancer gene matrix is keyed by Entrez Gene ID, not HGNC symbol. |
| `01b_fix_gene_de.R` | The pan-cancer gene matrix is keyed by HGNC SYMBOL (a minority of rows are bare Entrez |
| `02_ffl_census.py` | 02_ffl_census.py -- formal definition, exact census and Alon coherence typing of |
| `02_layer_TF_target.R` | A) TF-TF and TF-gene layer with literature-derived signs, from TRRUST v2 (human). |
| `02b_ffl_class_coverage.py` | 02b_ffl_class_coverage.py -- tests the manuscript's rationale for stopping at six |
| `02c_candidate_space.py` | """02c_candidate_space.py -- computes, exactly, the number of candidate vertex |
| `02d_handworked_checks.py` | """02d_handworked_checks.py -- the hand-worked sanity checks of the formal FFL |
| `02e_exemplar_audit.py` | """02e_exemplar_audit.py -- audits the authors' deposited exemplar circuits |
| `02f_crossvalidate.py` | """02f_crossvalidate.py -- independent cross-validation of scripts/ffl_enum.c. |
| `02g_nineclass_search.py` | """02g_nineclass_search.py -- exact targeted search for modules realising ALL NINE |
| `02h_exemplar_permissive.py` | 02h_exemplar_permissive.py -- re-tests the authors' deposited exemplars on the |
| `02x_ffl_core_crosscheck.py` | """Independent cross-check of the 3-node FFL census and coherence typing. |
| `02y_ffl_per_network.py` | """Enumerate 3-node FFL cores WITHIN each deposited network separately, after applying |
| `03_ffl_census.py` | n-node feed-forward loop census by connected-induced-subgraph enumeration. |
| `03_ffl_enumerate.py` | Formal enumeration of n-node feed-forward loops. |
| `03_layer_TF_miRNA.R` | B) TF-miRNA layer with literature-derived signs, from TransmiR v2 (human). |
| `04_layer_miRNA_miRNA.R` | C) miRNA-miRNA layer defined as POLYCISTRONIC CO-TRANSCRIPTION (genomic clustering), |
| `05_layer_gene_gene.R` | D) Gene-gene layer: two clearly separated evidence tiers. |
| `05a_string_map_and_fetch.R` | Exact symbol -> STRING v12 protein-ID mapping using the official STRING |
| `06_edge_evidence_tier.R` | E) Evidence tier for every miRNA_target edge in the canonical network. |
| `07_layer_summary.R` | B <- "/path/to/revision" |
| `08_exir_prep_run.R` | ============================================================================ |
| `08_expression_validation.R` | 08_expression_validation.R |
| `08_motif_node_sets.py` | 08_motif_node_sets.py |
| `08_motif_significance.py` | 08 -- Motif significance testing of the canonical breast-cancer regulatory network |
| `08a_control_pool_fetch.R` | 08a: build the node pools for the DISEASE-UNRELATED matched control networks and |
| `08b_mixing_check.py` | """08b -- convergence check on the Maslov-Sneppen mixing time. |
| `08c_reciprocity_audit.py` | """08c -- why the Composite-FFL result depends entirely on the null model. |
| `09_anova_ffl_named_axes.R` | 09_anova_ffl_named_axes.R |
| `09_enrichment_ora.R` | ============================================================================= |
| `09_ffl_hubs.R` | ============================================================================ |
| `09_hub_bootstrap.py` | 09 -- Bootstrap stability of the hub list. |
| `09b_claimed_terms_hypergeometric.R` | ============================================================================= |
| `09c_exemplar_reproduction.R` | REV<-"/path/to/revision" |
| `10_directional_evidence_and_anova_control.R` | 10_directional_evidence_and_anova_control.R |
| `10_exir_classify_overlap.R` | ============================================================================ |
| `10_motif_class_comparison.R` | ============================================================================= |
| `10_motif_specificity.py` | 10 -- SPECIFICITY control: is the FFL structure specific to breast cancer, or does it merely |
| `10_survival_cox_hubs.R` | 10_survival_cox_hubs.R |
| `10a_export_control_data.R` | 10a: export, as flat TSVs, the three edge sources used to rebuild BOTH the breast-cancer |
| `11_exir_mediator_enrichment.R` | ============================================================================ |
| `11_ffl_module_membership.R` | 11_ffl_module_membership.R |
| `11_gsea_signed_influence.R` | ============================================================================= |
| `11_specificity_enrichment.py` | 11 -- The scale-free half of the specificity control. |
| `12_exir_primary_class.R` | ============================================================================ |
| `12_module_score_survival.R` | 12_module_score_survival.R |
| `12_summarise_findings.py` | """12_summarise_findings.py - pull the response-letter numbers out of the result CSVs.""" |
| `13_circuit_level_3node_vs_6node.R` | 13_circuit_level_3node_vs_6node.R |
| `13_hub_rank_validation.R` | ============================================================================ |
| `14_external_met500.R` | 14_external_met500.R |
| `14_overlap_sweep.R` | ============================================================================ |
| `15_external_metabric.R` | 15_external_metabric.R |
| `16_external_geo_de.R` | 16_external_geo_de.R |
| `17_external_cohort_validation.R` | 17_external_cohort_validation.R |
| `18_metabric_treatment.R` | 18_metabric_treatment.R |
| `20_cptac_prep.R` | 20_cptac_prep.R -- harmonise CPTAC-BRCA proteome / phosphoproteome / RNA / miRNA |
| `20_figures_core.R` | Core revision figures: global network, the FFL census |
| `21_cptac_protein_validation.R` | 21_cptac_protein_validation.R |
| `21_figure_exemplars.R` | Redraw the three published higher-order circuits exactly as deposited, but DIRECTED and |
| `22_cptac_axes_mrnaprot_phospho.R` | 22_cptac_axes_mrnaprot_phospho.R |
| `23_rppa_second_platform.R` | 23_rppa_second_platform.R -- TCGA RPPA as an independent protein platform |
| `24_nfkb_activity_score.R` | 24_nfkb_activity_score.R -- NF-kB TRANSCRIPTIONAL ACTIVITY (target-gene signature) vs collagen, |
| `30_motif_significance.py` | Motif over-representation against randomised networks. |
| `ct_01_hpa.py` | """Human Protein Atlas single-cell RNA: cell-type source of network genes. |
| `ct_02_scrna_extract.sh` | Single streaming pass over the 178M-entry GSE176078 sparse matrix. |
| `ct_02b_pseudobulk_all.sh` | Second streaming pass: per-(gene x celltype_major) UMI sums and n-cells-expressing, |
| `ct_03_stromal_scores.R` | Stromal / fibroblast / epithelial content scores for TCGA-BRCA primary tumours. |
| `ct_04_celltype_source.py` | """Per-cell-type expression of the network's protein-coding nodes in 26 primary breast |
| `ct_05_caf_signature.py` | """Derive a data-driven CAF signature from GSE176078 that is deliberately |
| `ct_06_partial_corr.R` | Does the collagen signal survive adjustment for stromal / CAF content? |
| `ct_07_mirna_hostgenes.py` | """Cellular source of the miRNA loci, approached through their pri-miRNA HOST GENES. |
| `ct_08_network_stromal.R` | Network-wide: does the sign concordance of the edges survive adjustment for CAF content? |
| `ct_09_external_corr.R` | Independent-cohort replication of (i) the collagen/TF correlations and |
| `ct_10_metabric_survival.R` | C) Survival replication in METABRIC: hub set + collagen/TF module scores, |
| `ct_11_power_downsample.R` | Is the TCGA "no hub is prognostic" null a power artefact? |
| `ct_12_replication_summary.R` | Consolidated independent-cohort replication table. |
| `ct_13_stratified.R` | Stratified check: partial correlation on a covariate as collinear as the CAF score |
| `mo_01_mutation.R` | Part B: somatic non-silent mutation frequency in TCGA-BRCA for network genes/TFs |
| `mo_02_copynumber.R` | Part C: GISTIC2 thresholded copy number for network nodes + CN-expression correlation |
| `mo_03_probeselect.py` | RV="/path/to/revision"; CA=RV+"/cache/multiomics" |
| `mo_04_probemap.R` | RV<-"/path/to/revision"; CA<-file.path(RV,"cache/multiomics") |
| `mo_05_methylation.R` | Part A: TCGA-BRCA Illumina 450k promoter methylation of network nodes |
| `mo_06_mut_lengthmatched.R` | Length-matched control for the "network genes are more mutated" claim |
| `mo_07_integrate.R` | Part D: integrated multi-omic summary for the manuscript's featured hubs |

## scripts/v2/

| script | purpose (from its header) |
|---|---|
| `00_counting_convention.py` | """Resolve the composite-FFL counting ambiguity definitively. |
| `01_coherence_corrected.py` | """Coherent / incoherent typing (Mangan & Alon, PNAS 2003) on the CORRECTLY counted core set. |
| `01_de_limma.R` | Independent re-derivation of the TCGA-BRCA differential expression tables (task A) |
| `02_edge_corr_caf.R` | ========================================================================== |
| `03_null_constructions.R` | Sensitivity of the N8 null to how the decile-matched random pairs are drawn. |
| `03b_null_within_network.R` | 4th null construction: random pairs drawn from the NETWORK'S OWN node set, |
| `04_caf_adjusted_vs_matched_null.R` | N9 under the SAME (within-network-node, decile-matched) null used by the paper. |
| `05_mediation.R` | Task D(ii): formal causal-mediation analysis -- how much of TF -> COL1A1 is |
| `06_cptac.R` | Task E: CPTAC-BRCA replication of the miR-29 -> collagen axis at PROTEIN level, |
| `07_consolidate.R` | Consolidated claim-by-claim verification table -> results/v2/verify_expression.csv |
| `20_pam50_subtype.R` | ============================================================================ |
| `21_pancancer_axis.py` | TASK B: pan-cancer generalisation of the miR-29 -> collagen axis and of the |
| `22_cptac_pancancer.R` | ============================================================================ |
| `23_consolidate_generalisation.py` | Consolidation for the subtype / pan-cancer / CPTAC-pan-cancer task. |
| `24_within_subtype_pooled.R` | ============================================================================ |
| `30_pam50_subtype_v3.R` | =========================================================================== |
| `31_pancancer_axis_v3.py` | PART B. Pan-cancer generalisation of the miR-29 -> collagen and TF -> collagen axes. |
| `32_cptac_pancancer_v3.R` | =========================================================================== |
| `33_subtype_contrasts_v3.R` | Pairwise subtype contrasts + within-subtype pooled estimates (TCGA-BRCA, METABRIC) |
| `34_consolidate_v3.py` | """Consolidation / reporting for Parts A-C (fresh run 2026-09-09).""" |
| `35_pancancer_ecological_product_v3.py` | Ecological decomposition: per TCGA cohort compute rho(TF, CAF), rho(CAF, COL1A1), |
| `36_meta_v3.py` | """Random-effects (DerSimonian-Laird) meta-analysis across TCGA cohorts and CPTAC cohorts, |
| `37_final_checks_v3.py` | OUT="/path/to/revision/results/v2" |
| `v2_analyse_null.py` | RES = "/path/to/revision/results/v2" |
| `v2_assert_classical.py` | """Assert that at n=3 the formal definition D1-D4 reduces EXACTLY to the classical FFL |
| `v2_build_graph.py` | """v2 independent graph builder for motif verification. |
| `v2_build_network.py` | Independent re-derivation of the canonical miRNA-TF-gene FFL network from the |
| `v2_chip_atlas.py` | ChIP-seq evidence for the TF->target edges of the breast-cancer miRNA-TF FFL network. |
| `v2_chip_download.sh` | Download ChIP-Atlas hg38 "Target Genes" tables for a list of TF symbols at a given |
| `v2_chip_null.py` | Is the ChIP-seq support of the 770 published TF->target edges better than chance? |
| `v2_cl_00_prep.R` | Cell-line "stroma-free test": data preparation |
| `v2_cl_01_expr_corr.R` | ============================================================================ |
| `v2_cl_02_mirna.R` | ============================================================================ |
| `v2_cl_03_controls_nulls.R` | ============================================================================ |
| `v2_cl_04_crispr.R` | ============================================================================ |
| `v2_cl_05_drug.R` | ============================================================================ |
| `v2_cl_06_consolidate.R` | ============================================================================ |
| `v2_cl_07_checks.R` | REV<-"/path/to/revision"; DM<-file.path(REV,"data/depmap24q4") |
| `v2_claims.py` | """(c) case duplicates, (e) TRRUST, (f) manuscript claims, (g) diff vs existing canonical.""" |
| `v2_classdiv.py` | """Edge-class diversity per module for n=3..7 from the class-mask histograms.""" |
| `v2_core3.py` | """Exhaustive n=3 FFL census under the formal definition (D1-D4), plus the two |
| `v2_diag3.py` | """Diagnostics: n=3 core counts under every plausible convention, per graph variant.""" |
| `v2_dl_chunked.sh` | N-way parallel HTTP range download (targetscan.org throttles single connections hard). |
| `v2_encode_query.py` | """Third resource: ENCODE portal REST API - how many human TF ChIP-seq experiments exist for |
| `v2_export_graph.py` | sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))) |
| `v2_extra.py` | SIF="/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files" |
| `v2_ffl_grid.py` | """Brute-force grid over plausible 3-node FFL definitions, searching for which |
| `v2_ffl_variants.py` | """FFL-core counting under several explicit definitional variants, so the |
| `v2_graph_final.py` | """Final graph definitions used by all v2 motif work. |
| `v2_make_census_csv.py` | sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))) |
| `v2_pernet.py` | sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))) |
| `v2_remap_promoters.py` | Reverse query: which TFs are actually bound at the COL1A1 and COL3A1 promoters? |
| `v2_sweep3.py` | """Sweep additional augmentation variants (incl. layer_TF_miRNA) looking for the |
| `v2_targetscan.py` | TargetScan 8.0 sequence-level evidence for the miRNA->target edges of the network. |
| `v2_write_verdicts.py` | OUT="/path/to/revision/results/v2/verify_network.csv" |
| `v2d_chip_analyse.py` | """ChIP-seq support for the canonical network's 770 TF_target edges, from |
| `v2d_chip_fetch.py` | """Download ChIP-Atlas 'Target Genes' tables (hg38) for every TF that is a source |
| `v2d_chip_merge.py` | """Merge the three ChIP resources into one per-edge evidence table for Table 1, |
| `v2d_collagen_context.py` | """Which factors occupy the COL1A1 / COL3A1 promoters in FIBROBLAST, STROMAL or |
| `v2d_collagen_specificity.py` | """Which TFs are SPECIFICALLY bound at the COL1A1 / COL3A1 promoters? |
| `v2d_edge_rho.R` | Independent re-computation (this run) of edge-level Spearman rho in TCGA-BRCA |
| `v2d_encode.py` | """ENCODE portal: number of released human TF ChIP-seq experiments per network TF, |
| `v2d_integrate.R` | Does TargetScan site quality predict which network miRNA-target edges are |
| `v2d_remap_analyse.py` | """ReMap 2022 (hg38, all 8,103 non-merged datasets) promoter occupancy. |
| `v2d_remap_bg.py` | """Genome-wide background: ReMap 2022 promoter occupancy at 600 randomly chosen |
| `v2d_remap_bgnull.py` | """Genome-wide promoter null for ReMap edge support: for each network TF_target |
| `v2d_remap_fetch.py` | """Query ReMap 2022 (hg38, all non-merged datasets) for every ChIP-seq peak whose |
| `v2d_summary.py` | """Consolidated verification table for the direct-molecular-evidence task.""" |
| `v2d_targetscan.py` | """TargetScan 8.0 site-level evidence for the canonical network's miRNA_target edges. |
| `v2d_targetscan2.py` | """TargetScan 8.0 sequence-level evidence for the canonical network's miRNA_target |
| `v2d_ts_null.py` | """Are the network's 4,819 miRNA-target edges enriched for TargetScan sites over |
| `v2d_tss.py` | """Build hg38 TSS table for every gene/TF node of the canonical network from |

## scripts/v3/

| script | purpose (from its header) |
|---|---|
| `00_mir130a_direction.R` | Which direction does the data support for miR-130a? |
| `00_verify_references.py` | 00_verify_references.py |
| `01_architecture.py` | """(A) ARCHITECTURE of the canonical miRNA-TF-gene network. |
| `01_integrator_validation.py` | 01_integrator_validation.py |
| `01_mir130a_oncomir_test.R` | REV <- "/path/to/revision" |
| `01_mir130a_target_set.R` | A) DEFINE THE REAL miR-130a TARGET SET |
| `01b_architecture_nulls.py` | """(A cont.) Are the architectural statistics anything other than a consequence of the |
| `02_characterise_3node.py` | 02_characterise_3node.py  -- Part A + B. |
| `02_mir130a_decisive.R` | The decisive test of the oncomiR model for miR-130a. |
| `02_mir130a_tcga_correlations.R` | B) TEST EVERY miR-130a TARGET IN TCGA-BRCA (n = 1066 paired primary tumours) |
| `02_powerlaw_csn.R` | (A) Is the degree distribution scale-free? Clauset, Shalizi & Newman (2009) |
| `02b_summary_table.py` | 02b_summary_table.py |
| `02c_characterise_3node_N16384.py` | 02c_characterise_3node_N16384.py |
| `02d_adaptive_solver_crosscheck.py` | 02d_adaptive_solver_crosscheck.py |
| `02e_paired_tests_and_textbook_validation.py` | 02e_paired_tests_and_textbook_validation.py |
| `02f_pulse_longwindow.py` | 02f_pulse_longwindow.py |
| `02g_fcd_recomputed.py` | 02g_fcd_recomputed.py |
| `02h_stochastic_noise_buffering.py` | 02h_stochastic_noise_buffering.py |
| `02i_paper_table1.py` | 02i_paper_table1.py |
| `02j_noise_at_matched_mean.py` | 02j_noise_at_matched_mean.py |
| `03_controllability.py` | """(B) STRUCTURAL CONTROLLABILITY of the canonical network. |
| `03_higher_order_sweep.py` | 03_higher_order_sweep.py  --  Part C: DOES ORDER ADD FUNCTION? |
| `03_mir130a_adjusted.R` | The positive miR-130a~TSG correlations could be compositional (both higher in normal-like |
| `03_mir130a_positive_control_and_subtype.R` | B2) Is the NULL RESULT for miR-130a a failure of the assay, or a property of miR-130a? |
| `03c_wellposedness_check.py` | 03c_wellposedness_check.py |
| `03d_higher_order_compI1.py` | 03_higher_order_sweep.py  --  Part C: DOES ORDER ADD FUNCTION? |
| `03e_higher_order_convergence_audit.py` | 03e_higher_order_convergence_audit.py |
| `03f_n6_limit_cycle_certification.py` | 03f_n6_limit_cycle_certification.py |
| `03g_higher_order_empirical_support.py` | 03g_higher_order_empirical_support.py |
| `04_information_flow.py` | """(C) INFORMATION FLOW through the regulatory network. |
| `04_mir130a_clinical.R` | REV <- "/path/to/revision" |
| `04_mir130a_cptac_protein.R` | C) miR-130a target repression at PROTEIN level, CPTAC-BRCA |
| `04a_extract_mir29_data.R` | 04a_extract_mir29_data.R -- Part D data extraction. |
| `04b_mir29_circuit_fit.py` | 04b_mir29_circuit_fit.py  --  Part D: parameterise one concrete circuit from the data. |
| `04c_mir29_modelfree_prediction.py` | 04c_mir29_modelfree_prediction.py |
| `05_figures.py` | """05_figures.py -- publication figures for the dynamical analysis. |
| `05_mir130a_celllines.R` | REV <- "/path/to/revision"; R3 <- file.path(REV,"results/v3") |
| `05_mir130a_esr1_and_composition.R` | B4) Is the only strong anti-correlation (ESR1) a real repression signal, or a subtype marker? |
| `05_perturbation.py` | """(D) PERTURBATION AND CONTROL: in-silico single- and double-node knockouts. |
| `06_empirical_grounding.py` | 06_empirical_grounding.py |
| `06_knockout_landscape.py` | """(D) Complete single- and double-knockout landscape of the canonical network. |
| `06_mir130a_enrichment.R` | D) FUNCTIONAL COHERENCE of the anti-correlated validated miR-130a target set |
| `06b_composite_showcase.py` | 06b_composite_showcase.py |
| `06c_composite_bistability_export.py` | 06c_composite_bistability_export.py |
| `06d_composite_ringing_export.py` | 06d_composite_ringing_export.py |
| `07_consolidate.py` | 07_consolidate.py -- pull every headline number of the dynamical analysis into one |
| `07_dynamics.py` | """(E) DYNAMIC RESPONSE OF THE REAL NETWORK. |
| `07_mir130a_depmap.R` | E-i) Are miR-130a targets ESSENTIAL in breast cancer cell lines?  (DepMap CRISPR / Chronos) |
| `07b_dynamics_nulls.py` | """(E cont.) Nulls and self-consistency checks for the dynamic-response analysis. |
| `07c_consolidate_additions.py` | 07c_consolidate_additions.py |
| `08_mir130a_survival.R` | E-ii) PROGNOSIS: miR-130a itself, and a miR-130a target-activity score, in TCGA-BRCA |
| `08_mir_vs_txn_paired.py` | 08_mir_vs_txn_paired.py |
| `08_verify_citations.py` | """(F) Verify every reference used in the systems-pharmacology analysis against the live |
| `09_analytic_bounds.py` | 09_analytic_bounds.py |
| `09_combination_stats.py` | """(D cont.) Does multi-target perturbation actually beat single-target inhibition here? |
| `09_mir130a_methylation.R` | E-iii) Is the silencing REVERSIBLE in principle?  Methylation <-> expression at the miR-130a locus. |
| `09b_composite_oscillation.py` | 09b_composite_oscillation.py |
| `09c_relaxation_rate_theorem.py` | 09c_relaxation_rate_theorem.py |
| `10_figures.py` | """Figures for the systems-pharmacology / network-control analysis.""" |
| `10_mir130a_dac_geo.R` | E-iii b) Does DNA-demethylation (5-aza-2'-deoxycytidine, DAC) re-express miR-130a in breast lines? |
| `10b_dynamics_extra_figures.py` | 10b_dynamics_extra_figures.py |
| `11_master_summary.py` | """Consolidate every headline number produced by the v3 systems-pharmacology analysis into |
| `11_mir130a_celltype.py` | """F) CELL-TYPE SOURCE of miR-130a and of its targets. |
| `12_controllability_degree_control.py` | """(B cont.) Is the driver-node depletion among FFL hubs anything more than a degree effect? |
| `12_mir130a_gse19783.R` | EXTERNAL COHORT: GSE19783 (Enerly et al.) - matched miRNA (GPL8227) + mRNA (GPL6480) breast tumours. |
| `13_matching_crosscheck.py` | """(B cont.) Independent cross-check of the maximum-matching / driver-node computation. |
| `13_mir130a_consistency.R` | Which individual miR-130a targets hold up ACROSS all three independent tests? |
| `14_mir130a_literature.py` | """G) LITERATURE on miR-130a in breast cancer (PubMed E-utilities).""" |
| `14_regulator_subnetwork.py` | """(A cont.) Is the "hierarchical regulatory structure" real, or an artefact of how the |
| `15_igraph_crosscheck.R` | (A cont.) Independent cross-check, in R/igraph, of the topology numbers that the Python |
| `15_mir130a_extras.R` | Odds and ends: (a) arm-specific tumour-vs-normal test, (b) ESR1 relation in NORMAL breast, |
| `16_flow_vs_degree_sinkcontrol.py` | """(C cont.) The flow-vs-degree comparison, controlled for the structural sinks. |
| `16_mir130a_summary.R` | Consolidated headline table. Every number is read back out of the result files |
| `19_mir29_target_set.R` | 19) Build the miR-29 family target set (mirrors 01_mir130a_target_set.R) and its |
| `20_deconv_install.R` | 20_deconv_install.R -- install every deconvolution package we intend to use. |
| `20_screens_demeter2.R` | PART 1A) DepMap RNAi (DEMETER2 v6 combined: Achilles + DRIVE + Marcotte) as an |
| `20_seqreg_regions.py` | """20_seqreg_regions.py -- define hg38 loci for sequence-level regulatory analysis |
| `21_screens_projectscore.R` | PART 1B) Project Score (Sanger CRISPR-Cas9, Score2 release: Sanger v2 + Broad 21Q2), |
| `21_scrna_reference_build.py` | """21_scrna_reference_build.py |
| `21_seqreg_alphagenome.py` | """21_seqreg_alphagenome.py -- AlphaGenome (Google DeepMind, 2025) sequence-to-function |
| `21b_seqreg_alphagenome_access.py` | """21b -- record, reproducibly, exactly what was attempted to obtain AlphaGenome access, |
| `22_deconv_prep.R` | 22_deconv_prep.R -- build the matrices every deconvolution method needs. |
| `22_screens_drug_sensitivity.R` | PART 1C) Does a miR-29 target score, a miR-130a target score or a collagen score |
| `22_seqreg_enformer.py` | """22_seqreg_enformer.py -- Enformer (Avsec et al. 2021) sequence-to-function predictions |
| `22a_celline_scores.R` | 22a) Cell-line module scores used for the drug-sensitivity analysis (Part 1C). |
| `22b_node_drug_associations.R` | 22b) Per-NODE drug-sensitivity association, so that the integrated matrix (Part 1E) can |
| `22b_seqreg_enformer_ism.py` | """22b_seqreg_enformer_ism.py -- sequence attribution at the COL1A1 and COL3A1 promoters |
| `22c_seqreg_ism_motifmap.py` | """22c -- map the Enformer in-silico-mutagenesis peaks onto the JASPAR/HOCOMOCO motif sites |
| `23_cibersort_impl.R` | 23_cibersort_impl.R -- faithful open re-implementation of the CIBERSORT |
| `23_screens_orcs.py` | """PART 1D) Published phenotypic CRISPR screens (BioGRID ORCS 2.0.18, human). |
| `23_seqreg_borzoi.py` | """23_seqreg_borzoi.py -- Borzoi (Linder et al., Nat Genet 2025) sequence-to-function |
| `23b_screens_literature.py` | """PART 1D, second half) Literature search for published CRISPR screens in BREAST models |
| `24_screens_integrated_matrix.R` | PART 1E) One integrated essentiality matrix over the 587 network nodes: |
| `24_seqreg_chipatlas.py` | """24_seqreg_chipatlas.py -- ChIP-Atlas (Zou et al., NAR 2022) evidence that a TF is bound |
| `24_wu_signature_matrix.R` | 24_wu_signature_matrix.R |
| `25_run_deconvolution.R` | 25_run_deconvolution.R -- run every deconvolution method that installed, on the |
| `25_seqreg_remap.py` | """25_seqreg_remap.py -- ReMap 2022 (Hammal et al., NAR 2022) experimentally mapped TF peaks |
| `25b_seqreg_remap_celltypes.py` | """25b -- recompute the ReMap 2022 region tables with a CURATED cell-type (biotype) whitelist. |
| `25c_seqreg_remap_background.py` | """25c -- ReMap 2022 random-promoter background: for 300 random RefSeq-Select promoters, |
| `25d_seqreg_remap_background_counts.py` | """25d -- as 25c but storing the NUMBER of ReMap peaks per TF per background promoter, so |
| `26_instaprism_bayesprism.R` | 26_instaprism_bayesprism.R |
| `26_seqreg_screen_ccre.py` | """26_seqreg_screen_ccre.py -- ENCODE SCREEN candidate cis-regulatory elements (cCREs) |
| `27_assemble_scores.R` | 27_assemble_scores.R -- assemble every per-sample cell-type estimate into one |
| `27_seqreg_gtex.py` | """27_seqreg_gtex.py -- GTEx v8/v10 cis-eQTL and sQTL evidence for COL1A1, COL3A1 and the |
| `27b_seqreg_gtex_extra.py` | """27b -- GTEx addendum: (a) independent cis-eQTLs for the miRNA genes themselves, |
| `28_mediation_by_method.R` | 28_mediation_by_method.R -- THE DECISIVE RE-TEST. |
| `28_seqreg_motifscan.py` | """28_seqreg_motifscan.py -- FIMO-equivalent PWM scan of the COL1A1 / COL3A1 (and miRNA-locus) |
| `28b_seqreg_motif_calibration.py` | """28b_seqreg_motif_calibration.py -- higher-resolution GC-matched calibration of the |
| `28c_seqreg_motif_bh.py` | """28c -- Benjamini-Hochberg correction of the focus-motif calibration (25 motifs x 10 promoters |
| `29_correlations.R` | 29_correlations.R |
| `29_seqreg_cpg_architecture.py` | """29_seqreg_cpg_architecture.py -- promoter sequence architecture at the COL1A1/COL3A1 and |
| `30_sc_atlas2_pseudobulk.py` | """PART 2F) Second (and much larger) breast atlas: the Human Breast Cancer Single Cell Atlas |
| `30_seqreg_mir130a_probes.py` | """30_seqreg_mir130a_probes.py -- place the six Illumina 450k probes used in the project's |
| `30_simulation_benchmark.R` | 30_simulation_benchmark.R |
| `30a_wu_pseudobulk.sh` | Stream the Wu et al. 2021 (GSE176078) 29,733 x 100,064 count matrix once and accumulate |
| `31_gene_lengths.py` | """31_gene_lengths.py -- union-exon length per gene symbol from GENCODE v36 basic, |
| `31_sc_caf_subtypes.py` | """PART 2G) CAF subtypes in the Human Breast Cancer Single Cell Atlas stromal compartment. |
| `31_seqreg_enhancers.py` | """31_seqreg_enhancers.py -- candidate enhancers of COL1A1 and COL3A1: ENCODE cCREs within |
| `31b_wu_caf_subtypes.py` | """PART 2G, cross-check in the atlas the project already uses (Wu et al. 2021, GSE176078), |
| `31c_caf_summary_table.py` | """PART 2G summary: assign each CAF subtype to myCAF / iCAF / apCAF by marker score, and |
| `32_sc_celltype_enrichment.py` | """PART 2F (confirmation) and 2H (do the two miRNA target sets segregate by cell type?). |
| `32_seqreg_rela_enhancers.py` | """32_seqreg_rela_enhancers.py -- the only fibroblast-context TF-binding evidence at the |
| `32_tpm_length_sensitivity.R` | 32_tpm_length_sensitivity.R |
| `32b_sc_setcontrols.py` | """PART 2H controls.  (1) Is the CAF-weighting of the miR-29 target set just the collagen |
| `32c_wu_compartment_check.py` | """PART 2H, independent replication in the Wu et al. 2021 atlas (the project's original atlas), |
| `33_cibersort_lm22_variants.R` | 33_cibersort_lm22_variants.R -- CIBERSORT/LM22 goodness-of-fit under both |
| `33_sc_donor_pseudobulk.py` | """PART 2J) Donor x cell-type pseudobulk from the global atlas, so that the |
| `33_seqreg_synthesis.py` | """33_seqreg_synthesis.py -- assemble every sequence-level result into one tidy table |
| `34_collagenfree_variants.R` | 34_collagenfree_variants.R |
| `34_sc_tf_collagen_within_caf.R` | PART 2J) Does the TF -> COL1A1 association exist WITHIN the CAF compartment alone? |
| `34_seqreg_network_tf_edges.py` | """34_seqreg_network_tf_edges.py -- put every TF->collagen edge that the manuscript's own |
| `34b_caf_composition_control.py` | """PART 2J, control: is the surviving within-CAF TF -> COL1A1 association just CAF-subtype |
| `35_consensus_mediator.R` | 35_consensus_mediator.R -- build two method-agnostic consensus fibroblast axes and |
| `35_sc_cellchat_liana.py` | """PART 2I) Cell-cell communication on the Human Breast Cancer Single Cell Atlas, |
| `35_seqreg_nfkb_context.py` | """35_seqreg_nfkb_context.py -- every NF-kB-family ReMap 2022 peak within +/-100 kb of COL1A1 |
| `36_sc_cellchat.R` | PART 2I, named-method run: CellChat v2 on the Human Breast Cancer Single Cell Atlas |
| `36_seqreg_key_results_addendum.py` | """36 -- append the remaining headline numbers to results/v3/seqreg_key_results.csv.""" |
| `36a_cellchat_input.py` | """Export a subsampled, log-normalised expression matrix + cell labels for CellChat.""" |
| `36b_sc_cellchat_small.R` | PART 2I, named-method run: CellChat v2 on the Human Breast Cancer Single Cell Atlas |
| `40_seqreg_ext_alphagenome_access.py` | """40_seqreg_ext_alphagenome_access.py |
| `40_verify_screens_rnai.py` | """INDEPENDENT verification of PART 1A (DEMETER2 RNAi). |
| `41_deconv2_verify_existing.R` | 41_deconv2_verify_existing.R |
| `41_seqreg_ext_occlusion_null.py` | """41_seqreg_ext_occlusion_null.py |
| `41_verify_setlevel_null.py` | """Independent re-derivation of the DEMETER2 and Project SCORE set-level decile-matched null tests.""" |
| `42_deconv2_external_refs.R` | 42_deconv2_external_refs.R |
| `42_seqreg_ext_chipatlas_fibroblast.py` | """42_seqreg_ext_chipatlas_fibroblast.py |
| `42_verify_projectscore.py` | """Independent verification of PART 1B (Project SCORE / Sanger).""" |
| `42b_verify_projectscore_replicates.py` | """Project SCORE: handle replicate screens of the same model explicitly.""" |
| `43_deconv2_newmethods.R` | 43_deconv2_newmethods.R -- ADD independent deconvolution methods that were not in the |
| `43_seqreg_ext_gtex_singletissue.py` | """43_seqreg_ext_gtex_singletissue.py |
| `43_verify_drug.py` | """Independent verification of PART 1C (PRISM + GDSC2 drug-sensitivity associations).""" |
| `43b_verify_prism.py` | REV="/path/to/revision"; OUT=f"{REV}/results/v3"; CA=f"{REV}/cache/v7" |
| `43c_verify_orcs.py` | """Independent verification of PART 1D (BioGRID ORCS 2.0.18).""" |
| `44_deconv2_pathology_stroma.py` | """44_deconv2_pathology_stroma.py |
| `44_seqreg_ext_synthesis.py` | """44_seqreg_ext_synthesis.py -- assemble the extension results into one tidy table. |
| `44_verify_sc_atlas.py` | """Independent verification of PART 2F/2G: streams raw counts once and accumulates |
| `44b_verify_integrated_matrix.py` | """Independent verification of PART 1E (integrated essentiality matrix).""" |
| `45_circularity_check.R` | Quantify the circularity in PART 2H: do miR-29 and miR-130a-3p have OPPOSITE bulk |
| `45_deconv2_assemble.R` | 45_deconv2_assemble.R |
| `45_seqreg_ext_splicedonor_control.py` | """45_seqreg_ext_splicedonor_control.py |
| `46_compartment_decircularised.R` | PART 2H, decisive control.  Rebuild the "anti-correlated target" sets for miR-29 and |
| `46_deconv2_mediation.R` | 46_deconv2_mediation.R -- THE DECISIVE RE-TEST, extended. |
| `47_compartment_test_decirc.py` | """Compartment test (Part 2H) run on RAW vs CAF-ADJUSTED anti-correlated target sets, |
| `47_deconv2_bayesprism.R` | 47_deconv2_bayesprism.R |
| `48_deconv2_celltype_correlations.R` | 48_deconv2_celltype_correlations.R -- part (D). |
| `48_verify_within_caf.R` | Independent re-derivation of PART 2J: within-CAF TF -> COL1A1, with and without |
| `49_deconv2_benchmark_newmethods.R` | 49_deconv2_benchmark_newmethods.R |
| `49_verification_summary.py` | """Consolidated record of the independent verification pass over PART 1 (screens) and |
| `50_deconv2_matched_collagen_control.R` | 50_deconv2_matched_collagen_control.R |
| `51_deconv2_contradictions.R` | 51_deconv2_contradictions.R -- part (E). |
| `52_deconv2_methylation_epidish.R` | 52_deconv2_methylation_epidish.R |
| `53_deconv2_consolidate.R` | 53_deconv2_consolidate.R -- one compact, manuscript-ready summary of the multi-method |
| `90_independent_verification.py` | """INDEPENDENT re-derivation of the v3 systems-pharmacology headline numbers. |
| `91_verify_part2.py` | """Part 2 of the independent verification: (B) exhaustive criticality, (D) perturbation, |
| `92_verify_part3.py` | """Part 3: double-knockout matrix, combination null, dynamics, power-law cross-check.""" |
| `93_powerlaw_crosscheck.R` | Cross-implementation check of the Clauset-Shalizi-Newman fit: is the R/poweRlaw result |
| `94_verify_part4.py` | """Part 4: degree-adjusted controllability control, Boolean attractor, C-section spot checks.""" |
| `95_verify_finalise.py` | """Part 5: re-run the (A) architecture checks so they are written to file, verify the |
| `96_powerlaw_figure.R` | RES <- "/path/to/revision/results/v3" |
| `99_dynamics_driver.sh` | 99_dynamics_driver.sh |
| `99b_fcd_driver.sh` | waits for the enlarged 3-node ensemble, then recomputes fold-change detection |
| `99c_dynamics_driver2.sh` | 99c_dynamics_driver2.sh -- runs the analyses that were written but never executed, |

## scripts/v5/

| script | purpose (from its header) |
|---|---|
| `01_assemble_node_compendium.py` | """Assemble the quantitative compendium for the 30 prioritised nodes. |
| `01_prioritisation.R` | ============================================================================== |
| `02_dgidb_query.py` | """Query DGIdb v5 GraphQL for the 20 protein-coding prioritised nodes. |
| `10_theme.R` | Shared publication theme for all revision figures. |
| `11_fig_census_motif.R` | Figure: FFL census across module sizes, edge-class saturation, and motif significance |
| `12_fig_validation.R` | Figure: coherence typology, expression validation by evidence tier, and the mediation result. |
| `13_fig_prioritisation.R` | Figure: the 30 prioritised nodes, their evidence-domain profiles, and the divergence between |
| `14_fig_mir130a.R` | Figure: hsa-miR-130a - expression, target landscape, composition-adjusted relationships, |
| `15_fig_mir130a_benchmark.R` | Figure: the miRNA benchmark that makes the miR-130a negative interpretable. |
| `20_compendium.R` | REV <- "/path/to/revision"; R5 <- file.path(REV,"results/v5") |
| `21_tables.R` | Complete interaction and supplementary tables. |
| `31_reference_list.R` | Assemble the reference list in Springer (Nature-style author-year) format from the verified |
| `40_farmer_extract_sources.py` | 40_farmer_extract_sources.py |
| `40_node_compendium_assemble.py` | 40_node_compendium_assemble.py |
| `40b_chang_wound_lists.py` | """Extract the tumour-level wound-response gene lists of Chang et al. 2004 (PLoS Biol, |
| `41_build_signature_collection.R` | 41_build_signature_collection.R |
| `41_pubmed_literature.py` | """Literature volume + candidate references for the 30 prioritised nodes. |
| `42_farmer_scores_tcga.R` | 42_farmer_scores_tcga.R |
| `42_pubmed_targeted.py` | """Targeted PubMed lookups for specific claims; prints PMID/year/journal/title for selection.""" |
| `43_farmer_module_overlap.R` | 43_farmer_module_overlap.R  -- THE KEY TEST |
| `43_verify_pmids.py` | """Verify every PMID cited in node_literature_text.py against live NCBI esummary records. |
| `44_farmer_mediation.R` | 44_farmer_mediation.R |
| `44_write_compendium.py` | """Write results/v5/NODE_COMPENDIUM.md from node_compendium_table.csv, |
| `45_farmer_survival.R` | 45_farmer_survival.R |
| `45_pubmed_all_nodes.py` | """PubMed record counts for all 587 network nodes, to test whether network centrality |
| `45b_pubmed_all_nodes_par.py` | """PubMed record counts for all 587 network nodes (threaded, globally rate-limited to <3 req/s), |
| `46_helpers.R` | shared GEO series-matrix parsing helpers (extracted from 46_neoadjuvant_response.R) |
| `46_neoadjuvant_response.R` | 46_neoadjuvant_response.R |
| `46_synthesis_stats.py` | """Programme assignment and cross-cutting statistics for the compendium synthesis.""" |
| `47_literature_vs_topology.py` | """Does network centrality track how heavily a node has been studied? |
| `47_neoadj_survival_and_agreement.R` | 47_neoadj_survival_and_agreement.R |
| `48_neoadj_stratified.R` | 48_neoadj_stratified.R |
| `49_neoadj_meanz_sensitivity.R` | 49_neoadj_meanz_sensitivity.R |
| `50_stroma_axis_and_figures.R` | 50_stroma_axis_and_figures.R |
| `51_build_docx_condensed.sh` | set -euo pipefail |
| `60_census_nostring.py` | Sensitivity re-run of the RAND-ESU n-node FFL census (scripts/03_ffl_census.py) |
| `60_census_sensitivity.py` | 60_census_sensitivity.py |
| `61_tables_T1_T2_T4.py` | Publication tables T1 (data sources / filtering cascade), T2 (network composition) |

## scripts/v6/

| script | purpose (from its header) |
|---|---|
| `01_fetch_pdmr.py` | """Fetch NCI PDMR (cBioPortal study pancan_pdmr_2025) RNA-seq V2 RSEM expression |
| `01_parse_hpa_xml.py` | """Parse HPA v25.1 per-gene XML entries for the 26 genes of interest. |
| `02_parse_hpa_sc_and_tissue.py` | """HPA v25.1 per-gene XML: single-cell-type RNA (nCPM), consensus tissue RNA (nTPM), |
| `02_spatial_analysis.py` | Spatial transcriptomics test of the compartment argument. |
| `03_cox_crosscheck.R` | Cross-check HPA pathology-atlas prognostic calls against our own TCGA-BRCA and |
| `03_pdx_analysis.py` | PDX natural experiment. |
| `04_build_hpa_tables.py` | """Assemble the deliverable HPA tables for the FIG revision.""" |
| `04_pdx_foldchange_and_caf.py` | (a) Matched PDX-vs-originator log2 fold change, corrected for the compositional |
| `05_breast_singlecell.py` | """HPA breast-tissue-specific single-cell clusters (rna_single_cell_type_tissue.tsv, v23): |
| `05_figures_springer.py` | """05_figures.py -- publication figures for the dynamical analysis. |
| `05_spatial_full.py` | Spatial transcriptomics test of the compartment argument (full version). |
| `06_crosscheck_table.py` | """Part D: HPA Pathology-Atlas breast prognostic calls vs our own TCGA-BRCA and |
| `06_spatial_mediation_summary.py` | (1) Spot-level mediation: TF -> stromal content -> COL1A1/COL3A1/FN1, run inside |
| `07_figures.py` | """Figures for the PDX and spatial natural experiments.""" |
| `07_genomewide_direction_concordance.R` | Is the direction disagreement between HPA's breast prognostic calls and our own |
| `08_file_inventory.py` | """Part A: inventory of the HPA files actually obtained, with release, URL, bytes, rows.""" |
| `08_pdx_passage_trend.py` | """Dose-response: does the human collagen signal fall progressively with PDX passage, |
| `09_compartment_table.py` | """Part C: the compartment test at protein level. One row per gene, combining |
| `09_provenance.py` | BASE="/path/to/revision"; RES=f"{BASE}/results/v6" |
| `10_master_table.py` | """One master table: the 20 protein-coding prioritised nodes (10 TFs + 10 genes), |
| `10_theme.R` | Shared publication theme for all revision figures. |
| `11_fig_census_motif.R` | Figure: FFL census across module sizes, edge-class saturation, and motif significance |
| `12_fig_validation.R` | Figure: coherence typology, expression validation by evidence tier, and the mediation result. |
| `13_fig_prioritisation.R` | Figure: the 30 prioritised nodes, their evidence-domain profiles, and the divergence between |
| `14_fig_mir130a.R` | Figure: hsa-miR-130a - expression, target landscape, composition-adjusted relationships, |
| `15_fig_mir130a_benchmark.R` | Figure: the miRNA benchmark that makes the miR-130a negative interpretable. |
| `21_gwas_fetch_traits.py` | """Resolve breast-cancer EFO trait hierarchy via OLS4 and export the term list.""" |
| `21b_node_coordinates.py` | """GRCh38 coordinates for all 30 prioritised nodes + comparators. |
| `21c_gwas_breast_loci.py` | """Do the 30 prioritised nodes sit at genome-wide significant breast cancer risk loci? |
| `60_fig8_independent_validation.R` | Figure 4: spatial compartment occupancy, patient-to-xenograft change, bulk versus within-stroma |
| `61_fig5_reactive_stroma.R` | Figure 5: FFL module scores against the Farmer stroma-related signature, and the pCR meta-analysis. |
| `62_fig6_gain_loss.R` | Figure 6: behaviours gained and lost by higher-order modules relative to their embedded three-node |

## scripts/v7/

| script | purpose (from its header) |
|---|---|
| `00_build_graph.R` | v7 / 00 : build the undirected simple graph used by every community-detection |
| `01_cm2_mcl.R` | v7 / 01 : Markov clustering (van Dongen 2000), the algorithm clusterMaker2 |
| `02_clustermaker2_cyrest.py` | v7 / 02 : drive the REAL clusterMaker2 (v2.3.4) through Cytoscape 3.10.0 + CyREST. |
| `03_ffl_census_nolegacymirna.py` | n-node feed-forward loop census by connected-induced-subgraph enumeration. |
| `03_igraph_battery.R` | v7 / 03 : the same battery run independently in igraph 2.x, so that every |
| `04_analysis.R` | v7 / 04 : harmonise every partition (clusterMaker2 2.3.4 via CyREST + igraph), |
| `05_consensus_and_stability.R` | v7 / 05 : degeneracy filtering, stochastic-seed stability, consensus |
| `06_null_and_prioritisation.R` | v7 / 06 : (a) characterise the consensus modules; |
| `07_stroma_module_test.R` | v7 / 07 : the direct test of the paper's claim. |
| `08_degree_matched_null.R` | v7 / 08 : the ECM-concentration test of script 07 repeated with a DEGREE-MATCHED |
| `09_subnetwork_sensitivity.R` | v7 / 09 : sensitivity to network composition. 4,819 of the 6,859 directed edges |
| `10_compile_headline.R` | v7 / 10 : compile the headline numbers quoted in the write-up. |
| `11_figure.R` | v7 / 11 : one summary figure for the convergent-validity test. |
| `jam_01_search.R` | jam_01_search.R -- jActiveModules (Ideker et al. 2002) active-subnetwork |
| `jam_02_characterise.R` | jam_02_characterise.R -- what the jActiveModules search recovered: |
| `jam_03_landscape.R` | jam_03_landscape.R -- the shape of the Ideker objective on this network. |
| `jam_04_permutation_control.R` | jam_04_permutation_control.R |
| `jam_05_circularity_check.R` | jam_05_circularity_check.R -- the agreement between the active module and the |
| `v7_go_00_assemble_modules.R` | ============================================================================= |
| `v7_go_01_ora.R` | ============================================================================= |
| `v7_go_02_summarise.R` | ============================================================================= |
| `v7_go_03_named_terms.R` | ============================================================================= |
| `v7_go_04_convergence.R` | ============================================================================= |
| `v7_go_05_answer.R` | v7_go_05_answer.R -- consolidate the ORA into the answer table. |
| `v7_mcode_characterise.R` | --------------------------------------------------------------------------- |
| `v7_mcode_control_sweep.R` | v7 -- the negative controls, swept across the full node-score-cutoff range, |
| `v7_mcode_controls.R` | --------------------------------------------------------------------------- |
| `v7_mcode_run.R` | --------------------------------------------------------------------------- |
| `v7_mcode_validate.R` | Validation of the R port against the reference Java implementation's own |
