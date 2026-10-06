# Module detection (Supplementary Note S7; Table S9)

Community detection: `01_build_graph.R` → `02_cm2_mcl.R`, `03_clustermaker2_cyrest.py`, `04_igraph_battery.R` → `05`–`12`.
MCODE: `mcode.R` → `mcode_01` … `mcode_05`. BioNet: `bionet_01` … `bionet_07`. jActiveModules: `jam_common.R` with
`jam_engine.cpp` → `jam_01` … `jam_05`. GO over-representation: `go_01` … `go_06`. Re-runs on the analysed network:
`analyses/analysed_network_reruns/v7`, `v7b`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_build_graph.R` | Build the undirected simple graph used by every community-detection algorithm in the clusterMaker2 battery, and export it for Cytoscape. |
| `02_cm2_mcl.R` | Markov clustering (van Dongen 2000), the algorithm clusterMaker2 implements as its "MCL Cluster" command. |
| `03_clustermaker2_cyrest.py` | Drive the REAL clusterMaker2 (v2.3.4) through Cytoscape 3.10.0 + CyREST. |
| `04_igraph_battery.R` | The same battery run independently in igraph 2.x, so that every clusterMaker2 result has a second implementation behind it, and so that walktrap (absent from clusterMaker2 2.3.4) is covered. |
| `05_analysis.R` | Harmonise every partition (clusterMaker2 2.3.4 via CyREST + igraph), then the four questions the convergent-validity test asks: (a) how many clusters, what modularity, what size distribution (b) which cluster holds COL1A1 / … |
| `06_consensus_and_stability.R` | Degeneracy filtering, stochastic-seed stability, consensus co-clustering, and the two headline tests: does the stroma module exist, and does it match the prioritised 30? |
| `07_null_and_prioritisation.R` | (a) characterise the consensus modules; (b) degree-preserving rewiring null for the COL1A1/COL3A1/miR-29a triad; (c) agreement between every informative partition and the 30 prioritised nodes (hypergeometric, and a concentration … |
| `08_stroma_module_test.R` | Do the ECM and reactive-stroma genes concentrate in one module, and does that module also hold miR-29a and the TF arm (ETS1, NFKB1, RELA, SP1)? |
| `09_degree_matched_null.R` | The ECM-concentration test of script 07 repeated with a DEGREE-MATCHED null. |
| `10_subnetwork_sensitivity.R` | Sensitivity to network composition. 4,819 of the 6,859 directed edges are miRNA->target, so the community structure of the whole graph is dominated by shared-miRNA-target bipartite structure. |
| `11_compile_headline.R` | Compile the headline numbers of the module-detection analysis into one table. |
| `12_figure.R` | One summary figure for the convergent-validity test. |
| `bionet_01_bum_and_scan.R` | BioNet (Beisser et al. 2010, Bioconductor 1.68.0) Step 1: build undirected graph, attach DE p-values, fit beta-uniform mixture, scan FDR -> module size. |
| `bionet_02_main.R` | BioNet step 2: primary + sensitivity maximum-scoring modules |
| `bionet_03_sensitivity.R` | BioNet step 3: sensitivity of the maximum-scoring module to every analyst choice (BUM fitting set, miRNA handling, missing p-values, size rule) |
| `bionet_04_stability_and_nulls.R` | BioNet step 4: (a) edge-resampling stability of module membership (b) is the overlap with the prioritised 30 anything more than shared dependence on the same DE p-values? |
| `bionet_05_figures.R` | BioNet step 5: supplementary figure (print, light surface) |
| `bionet_06_compartment.R` | BioNet step 6: is the maximum-scoring module a fibroblast-compartment module? |
| `bionet_07_consolidate.R` | BioNet step 7: one consolidated headline table |
| `go_01_assemble_modules.R` | Assemble EVERY module of size >= 5 recovered by the four v7 module-detection runs into one tidy inventory, ready for BiNGO-equivalent ORA. |
| `go_02_ora.R` | BiNGO-equivalent GO-BP + KEGG over-representation for every module of size >= 5 recovered by the four v7 methods. |
| `go_03_summarise.R` | Classify every SIGNIFICANT term (BH < 0.05, Count >= 3) into three themes and answer: do the recovered modules annotate to ECM/collagen/stroma biology, to cell cycle / proliferation, or to something else? |
| `go_04_named_terms.R` | p and BH-adjusted p of the named GO terms (extracellular matrix and collagen terms, among others) in every module, whether or not they pass. |
| `go_05_convergence.R` | Cross-method convergence of the ORA result. |
| `go_06_answer_per_module.R` | Consolidate the over-representation results into one table per module. |
| `jam_01_search.R` | jActiveModules (Ideker et al. 2002) active-subnetwork search on the canonical miRNA-TF network, primary run + sensitivity grid. |
| `jam_02_characterise.R` | What the jActiveModules search recovered: reproducibility across seeds, overlap with the paper's prioritised nodes and module, and whether the active modules are stromal or proliferative. |
| `jam_03_landscape.R` | The shape of the Ideker objective on this network. |
| `jam_04_permutation_control.R` | (a) the component structure of the final annealed state (why ranks 2-5 of a single run are not modules), and (b) the control the corrected score does not provide: re-run the whole search on node scores permuted over the network. |
| `jam_05_circularity_check.R` | Circularity check: correlation of the jActiveModules node scores and selection frequencies with each domain of the prioritisation score, whose differential-expression domain comes from the same TCGA-BRCA contrast. |
| `jam_common.R` | Shared machinery for the jActiveModules re-implementation Ideker, Ozier, Schwikowski & Siegel (2002) Bioinformatics 18:S233-S240. |
| `jam_engine.cpp` | Simulated-annealing active-subnetwork search Faithful re-implementation of the jActiveModules search of Ideker, Ozier, Schwikowski & Siegel (2002) Bioinformatics 18:S233-S240. |
| `mcode.R` | a faithful R/igraph port of the Cytoscape MCODE app (Bader & Hogue, BMC Bioinformatics 2003;4:2) |
| `mcode_01_validate.R` | Validation of the R port against the reference Java implementation's own unit tests (MCODEAlgorithmTest.java) and against ProNet::mcode. |
| `mcode_02_run.R` | MCODE (Bader and Hogue 2003) on the network, as an independent test of the higher-order reactive-stroma module. |
| `mcode_03_characterise.R` | Characterisation of the MCODE result: agreement with the paper's prioritisation, placement of the TF arm, edge-evidence provenance of the top module, and a degree-preserving null for the top module's MCODE score. |
| `mcode_04_controls.R` | Controls for the MCODE result: does MCODE still find the module without the 30 miRNA–miRNA edges among the nine miRNAs of the exemplar module? |
| `mcode_05_control_sweep.R` | The negative controls of mcode_04_controls.R, swept across the full node-score-cutoff range. |
