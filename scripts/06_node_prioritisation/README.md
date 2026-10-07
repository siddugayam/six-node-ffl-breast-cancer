# Hubs and node prioritisation (Results 3.3; Notes S3 and S5)

Hubs and ExIR: `09_ffl_hubs.R`, `09_hub_bootstrap.py` (first-round bootstrap; Table S3 comes from
`analyses/analysed_network_reruns/hubboot/31_hub_bootstrap_nolegacy.py`), `08_exir_prep_run.R`,
`10`–`12_exir_*`. Prioritisation: `01_prioritisation.R` (Note S3, Fig. S4, Table S1). Node compendium and literature (Note S5): `40_node_compendium_assemble.py`, `41`–`45_pubmed_*`,
`43_verify_pmids.py`, `46_synthesis_stats.py`, `47_literature_vs_topology.py`. `20_compendium.R` writes the partner table
(Table S5) and an earlier, smaller `node_compendium_table.csv`; run it before `40` and `46`, which write the deposited
version of that table.

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_prioritisation.R` | Pre-specified composite prioritisation of network nodes. |
| `02_dgidb_query.py` | Query DGIdb v5 GraphQL for the 20 protein-coding prioritised nodes. |
| `08_exir_prep_run.R` | ExIR (Experimental-data-based Integrative Ranking) on TCGA-BRCA Package: influential 2.3.2 (Salavaty, Ramialison & Currie) classification. |
| `09_ffl_hubs.R` | Define "high-centrality hubs in the FFLs" on the canonical network |
| `09_hub_bootstrap.py` | Bootstrap stability of the hub list. |
| `10_exir_classify_overlap.R` | B) build results/exir_classification.csv C) ExIR-driver vs FFL-hub overlap (hypergeometric, BOTH directions) D) mediator \|logFC\| distribution vs drivers vs biomarkers |
| `11_exir_mediator_enrichment.R` | Enrichment of the ExIR mediators, with every gene tested for differential expression as the background universe. |
| `12_exir_primary_class.R` | ExIR does NOT return a partition: the Driver table and the Biomarker table each contain EVERY significant DE feature (same 4311 features, two different rankings). |
| `20_compendium.R` | Assembles per-node evidence for the node compendium (results/v5/node_compendium_table.csv). |
| `40_node_compendium_assemble.py` | Assemble every quantitative field for the 30 prioritised nodes (10 TF, 10 gene, 10 miRNA) into results/v5/node_compendium_table.csv, plus full and top-N partner tables. |
| `41_pubmed_literature.py` | Literature volume + candidate references for the 30 prioritised nodes. |
| `42_pubmed_targeted.py` | Targeted PubMed lookups for statements about the prioritised nodes; prints PMID, year, journal and title. |
| `43_verify_pmids.py` | Verify every PMID cited in node_literature_text.py against live NCBI esummary records. |
| `45_pubmed_all_nodes.py` | PubMed record counts for all 587 network nodes, to test whether network centrality tracks how heavily a node has been studied (the core argument of PRIORITISATION.md). |
| `45b_pubmed_all_nodes_par.py` | PubMed record counts for all 587 network nodes (threaded, globally rate-limited to <3 req/s), to test whether network centrality tracks how heavily a node has been studied. |
| `46_synthesis_stats.py` | Programme assignment and cross-cutting statistics for the compendium synthesis. |
| `47_literature_vs_topology.py` | Does network centrality track how heavily a node has been studied? |
| `node_literature_text.py` | Curated literature notes and one-line summaries for the 30 prioritised nodes. |
