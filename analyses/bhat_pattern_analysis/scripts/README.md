# Scripts

## Scripts in this folder

| Script | What it does |
|---|---|
| `a1_analyse.py` | A1 of SETTINGS.md: prevalences and comparisons for the Bhat pattern modules. |
| `a1_certify.py` | A1 of SETTINGS.md: certification of sustained oscillation for the A1 modules, exactly as analyses/six_node_pattern S2/s2_certify_limit_cycles.py does it (after scripts/11_dynamics/03f_n6_limit_cycle_certification.py): every set … |
| `a1_config_coverage.py` | A1 (added 2026-09-29; SETTINGS.md, Changes and decisions, 11:24): how many instances each chosen six-node sign configuration covers. a1_configurations.csv gives the greedy cover's newly-covered count: instances compatible with … |
| `a1_configurations.py` | A1 of SETTINGS.md: the sign configurations of the Bhat patterns to model, and the module list. |
| `a1_nested_certify.py` | A1 nested check, certification part (SETTINGS.md, Changes and decisions, 2026-09-28 20:05): the sets that the nested run (extended model, TF2 -> TF1 absent) flags is_sustained_osc are certified exactly as a1_certify.py does it … |
| `a1_nested_compare.py` | A1 nested check (SETTINGS.md): the extended model (dyn_models_bhat.py) with TF2 -> TF1 absent reruns the paper's stored six-node modules (COMP_C2_toggle and COMP_I1_negfeedback, n6 = GG+MM+TT, 16,384 sets, seed 20260908). |
| `a1_nested_table.py` | A1 nested check, table (SETTINGS.md, Changes and decisions, 2026-09-28 20:05). |
| `a1_run.py` | A1 of SETTINGS.md: runs modules of the extended model (dyn_models_bhat.py) with the paper's run_module() (scripts/11_dynamics/03_higher_order_sweep.py: integration, scoring and thresholds unchanged) and the chunked runner of … |
| `a1b_analyse.py` | A1b of SETTINGS.md: the claim guard against the full layer factorial. |
| `a1b_certify.py` | A1b of SETTINGS.md: certification of the sustained oscillation flagged in the A1b modules. |
| `a1b_modules.py` | A1b of SETTINGS.md: the module lists. configurations the nine six-node configurations that passed A1's claim guard for memory AND pulse and count for the verdict (TF2 -> TF1 as activation or repression; miRNA-FFL without the … |
| `a1b_nested_check.py` | A1b nested check (SETTINGS.md A1b; it gates the A1b runs). |
| `a1b_run.py` | A1b of SETTINGS.md: a1_run.py with the four layers of S2's factorial that A1 did not need. |
| `bhat_common.py` | Shared loader for the analyses of the twelve Bhat FFL networks (see SETTINGS.md). |
| `bhat_null.c` | The S1 null program of analyses/six_node_pattern (s1_null6.c), with counters added for the Bhat three-, four- and five-node patterns in the three classes (see "ADDED" below). |
| `cluster_nulls.py` | A2 of SETTINGS.md, item 2: the twelve Bhat networks under the two cluster variants, on the observed graph and on every replicate graph of NULL-B, NULL-C and NULL-L (section O's program and settings; each replicate's arc list is … |
| `concordance.R` | Section X of SETTINGS.md: sign concordance of the twelve Bhat networks' edges in TCGA-BRCA and CPTAC, against matched random-pair nulls, with the methods of … |
| `dyn_models_bhat.py` | A copy of scripts/11_dynamics/dyn_models.py with two additions for the Bhat patterns (A1 of SETTINGS.md); nothing else is changed: (1) TF2 -> TF1 (Topology.s_T21: +1 activation, -1 repression, 0 absent, the default). |
| `enrichment.R` | Section E of SETTINGS.md: enrichment of the twelve Bhat networks against the tumour vs normal signature. |
| `evidence_nulls.py` | A3 of SETTINGS.md, item 2: the twelve Bhat networks under NULL-B, NULL-C and NULL-L on two restricted graphs, with section O's program (bhat_null.c) and settings (1,000 replicates, seed 20250908, 100 swap attempts per rewired … |
| `knocktf_tf_arcs.py` | A4 of SETTINGS.md: KnockTF 2.0 support for the TF-TF arcs of the six-node Bhat patterns, with the rule of P3b (analyses/perturbation_tests), unchanged. |
| `nulls.py` | Section O of SETTINGS.md: over-representation of the twelve Bhat networks under NULL-B, NULL-C and NULL-L. |
| `ora_overlap.py` | ORA term overlap (SETTINGS.md, Changes and decisions, 2026-09-28 20:14); descriptive, part of no decision line. |
| `prioritisation.R` | Section P of SETTINGS.md: the paper's prioritisation (scripts/06_node_prioritisation/01_prioritisation.R) with Bhat FFL participation. |
| `run_logged.py` | Runs one analysis step and appends a line to RUN_LOG.tsv (next to SETTINGS.md): analysis, command, start, end, seconds, exit status, and the md5 of SETTINGS.md at the start. |
| `survival.R` | Section S of SETTINGS.md: each circuit of the four-, five- and six-node Bhat networks against its three-node core. |
| `topology.py` | Section T of SETTINGS.md: topology of the twelve Bhat FFL networks. |
| `variant_networks.py` | Variant networks for the reruns of A2 and A3 (SETTINGS.md), written in the network deposit's file format so that enrichment.R, concordance.R, survival.R and prioritisation.R can read them: <network>_nodes.tsv node, type, n_ffl … |
