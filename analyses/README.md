# Analyses

Each folder is one analysis added after the main rounds, with its code, the inputs it reads and the outputs it wrote,
in the layout in which it was run. Scripts refer to other analyses as `analyses/<name>/…` under the analysis root.

| Folder | Paper | What it holds |
|---|---|---|
| [`census_and_motif_nulls/`](census_and_motif_nulls) | Table 3, Fig. 1 | The census of 3- to 7-node FFL modules and the three-node motif null models, on four versions of the network |
| [`six_node_pattern/`](six_node_pattern) | Note S2 | Tests of the six-node composite pattern: over-representation, dynamics, signed circuits, robustness and survival |
| [`six_node_pattern_networks/`](six_node_pattern_networks) | Note S2 | The three- to six-node composite FFL networks rebuilt on the census graph, with topology and MCODE |
| [`six_node_followups/`](six_node_followups) | Notes S1, S2 | Follow-up questions on the miRNA–TF pair filter and the six-node pattern |
| [`analysed_network_reruns/`](analysed_network_reruns) | Notes S4, S7 | Re-runs without the 30 exemplar miRNA–miRNA edges: architecture, hubs, Tables S2, S3 and S6, module detection |
| [`dynamics_checks/`](dynamics_checks) | Results 3.7 | Checks of the dynamical models and the four-node TF–TF control |
| [`metabric_power_check/`](metabric_power_check) | Results 3.8 | METABRIC downsampled to the TCGA event count |
| [`perturbation_tests/`](perturbation_tests) | Note S1 | Tests with public perturbation data, purified-cell miRNA atlases and a second sequence model |
| [`perturbation_checks/`](perturbation_checks) | Notes S1, S2 | Three later checks of the perturbation tests |
| [`original_submission_code/`](original_submission_code) | — | The network files and R scripts of the original submission (MIT licence) |
