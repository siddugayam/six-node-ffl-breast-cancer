# Feed-forward loops: three-node cores and the n-node census (Methods 2.2, Fig. S2, Supplementary Note S2)

- Definitions and three-node cores: `01a_counting_convention.py`, `01b_coherence_corrected.py` (Mangan–Alon coherence
  sub-types, Table S4); `ffl_def.py`, `ffl_graph.py`, `ffl_cores.py` → `02_ffl_census.py` and the follow-up analyses
  `02b` … `02y` (edge-class coverage, candidate space, the deposited exemplar circuits, cores per source network).
- n-node census: `03_ffl_census.py`, `03_ffl_enumerate.py` and the C enumerators `ffl_enum.c` and `ffl_enum2.c` (first
  round); `census_01` … `census_12`, the first census under the formal definition (graph building, the C program
  `census_09_esu_classmask.c` with its run scripts `run_*.sh`, class diversity and the count table); `reconcile_01` …
  `reconcile_09`, the reconciliation of these counts with the conventions of the original submission and the effect of
  the STRING arcs; and `03_ffl_census_nolegacymirna.py` (the census graph without the 30 exemplar edges).
- `11_ffl_module_membership.R`: node sets of the 3- to 6-node modules. `60_census_*`: STRING sensitivity.
The census of the paper was re-run in `analyses/census_and_motif_nulls/`.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01a_counting_convention.py` | Composite-FFL counting convention: one motif instance or two. |
| `01b_coherence_corrected.py` | Coherent / incoherent typing (Mangan & Alon, PNAS 2003) on the unique cores (01a convention). |
| `02_ffl_census.py` | Formal definition, exact census and Alon coherence typing of n-node feed-forward loops in the canonical breast-cancer miRNA-TF-gene network. |
| `02b_ffl_class_coverage.py` | Smallest module size that admits one edge of every regulatory type (TF–TF, TF–miRNA, miRNA–miRNA, miRNA–gene, TF–gene, gene–gene), and the edge classes of larger modules. |
| `02c_candidate_space.py` | Computes, exactly, the number of candidate vertex sets the exhaustive enumeration must examine for each module size n (the figures in section 8 of results/ffl_definition_check.md). |
| `02e_exemplar_audit.py` | Tests the exemplar circuits of the original submission (SIF_files/4-TF.sif, 5-TF.sif, 6-TF.sif) against the formal definition: edge support, whole-structure size, exhaustive test of every induced n-subset, and the subsets … |
| `02g_nineclass_search.py` | Exact targeted search for modules realising ALL NINE interaction classes present in G. |
| `02h_exemplar_permissive.py` | Re-tests the exemplar circuits of the original submission on the PERMISSIVE graph, in which the 941 undirected STRING associations are admitted as orientable edges (so COL1A1--COL3A1, which has no directed evidence anywhere, does … |
| `02i_exemplar_whole_circuit_check.py` | Test the whole six-node exemplar circuit (its 10 molecules) against conditions D1-D4 (Methods 2.2) in two graphs, from the deposited files of the public repository only: (1) the analysed network (Table S2 rows with … |
| `02x_ffl_core_crosscheck.py` | Independent cross-check of the 3-node FFL census and coherence typing. |
| `02y_ffl_per_network.py` | Enumerate 3-node FFL cores WITHIN each of the six networks of the original submission separately, after applying only the identifier harmonisation. |
| `03_ffl_census.py` | n-node feed-forward loop census by connected-induced-subgraph enumeration. |
| `03_ffl_census_nolegacymirna.py` | n-node feed-forward loop census by connected-induced-subgraph enumeration. |
| `03_ffl_enumerate.py` | Formal enumeration of n-node feed-forward loops. |
| `11_ffl_module_membership.R` | Enumerate the node sets of the three- to six-node FFL modules (composite core; plus a gene–gene edge at four nodes, a miRNA–miRNA edge at five and a TF–TF edge at six) and write data/ffl_module_sets.rds. |
| `60_census_nostring.py` | Sensitivity re-run of the RAND-ESU n-node FFL census (scripts/03_ffl_census/03_ffl_census.py) with STRING co-functional associations EXCLUDED from the gene-gene layer. |
| `60_census_sensitivity.py` | Sensitivity of the n-node feed-forward-loop census to three construction choices of the census graph. |
| `census_01_core_definition_variants.py` | FFL-core counting under several explicit definitional variants, to match the counts of the original submission to a definition. |
| `census_02_definition_grid.py` | Brute-force grid over plausible 3-node FFL definitions, compared with the 6,037 three-node FFLs of the original submission. |
| `census_03_build_graph.py` | Independent builder of the census graph, used to cross-check the motif counts. |
| `census_04_n3_graph_variants.py` | Exhaustive n=3 FFL census under the formal definition (D1-D4), plus the two alternative 'core' counting conventions, for several graph-augmentation variants. |
| `census_05_n3_per_network.py` | Three-node FFL count (conditions D1-D4) for each of the six networks of the original submission, selected by the in_networks column of data/canonical_edges.tsv; prints the counts. |
| `census_06_graph_definitions.py` | Definitions of the census graphs used by the census and motif scripts. |
| `census_07_export_graphs.py` | Writes a census graph in the plain-text format that the C census programs read (counts of nodes, arcs and edge classes; node types; one line per arc with its class) and a .meta.json with node names and classes. |
| `census_09_esu_classmask.c` | Independent n-node FFL census under the formal definition. |
| `census_11_class_diversity.py` | Edge-class diversity per module for n=3..7 from the class-mask histograms. |
| `census_12_census_table.py` | Collects the counts of the C census runs from their logs (logs/v2) into results/v2/verify_census.csv: count per graph and module size, counting method, seeds, maximum number of edge classes and the ratio to the counts of the … |
| `ffl_cores.py` | Exact enumeration and Mangan & Alon (2003) coherence typing of every 3-node FFL core (R->M, R->T, M->T) in the canonical directed signed graph. |
| `ffl_def.py` | INDEPENDENT reference implementation of the formal n-node FFL definition. |
| `ffl_enum2.c` | Ffl_enum.c -- exact enumeration of n-node feed-forward-loop (FFL) modules |
| `ffl_enum.c` | Exact enumeration of n-node feed-forward-loop (FFL) modules |
| `ffl_graph.py` | Construction of the canonical DIRECTED, SIGNED, TYPED regulatory graph G used for the n-node feed-forward-loop (FFL) census. |
| `reconcile_01_graph_and_n3.py` | Reconciliation pass: rebuild the census graph exactly as 03_ffl_census.py does, report its size, and re-enumerate n=3 FFLs exhaustively under D1-D4 on both the augmented graph and the deposited-only graph. |
| `reconcile_02_n3_breakdown.py` | Break the n=3 FFL census down by node-type composition, on both graphs. |
| `reconcile_04_export_graph.py` | Export the exact graph 03_ffl_census.py builds, for the C enumerator. |
| `reconcile_05_census_greedy_vs_maxflow.c` | Exhaustive / sampled ESU census of n-node FFLs under D1-D4, replicating |
| `reconcile_06_published_convention.py` | Reproduce the counting convention behind the published 1,434 / 206 / 9 / 1,649, directly from data/canonical_edges.tsv. |
| `reconcile_08_string_impact.py` | How much of the n=3 census depends on STRING-orientation arcs? |
| `reconcile_09_no_string_census.py` | Three-node census with the undirected STRING tier excluded. |
| `run_dep67.sh` | Census of the deposited-network census graph (results/v2/graph_dep_fine.txt) with the C program census_09_esu_classmask: exhaustive at six nodes, RAND-ESU sampling with three seeds at seven nodes. |
| `run_n6s.sh` | RAND-ESU census at six nodes of the published-network census graph (results/v2/graph_pub_fine.txt) with census_09_esu_classmask, five seeds. |
| `run_n7.sh` | RAND-ESU census at seven nodes of the published-network census graph with census_09_esu_classmask, five seeds. |
| `run_n7b.sh` | RAND-ESU census at seven nodes of the published-network census graph with census_09_esu_classmask, three further seeds with a second set of sampling probabilities. |
