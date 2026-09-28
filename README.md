# Higher-order feed-forward loops in the breast cancer regulatory network

Code, network and result tables for:

> Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad. *Deciphering the
> Transcriptional Regulation of Breast Cancer Genes through a Novel Six-Node Feed-Forward Loop
> Analysis.* Functional & Integrative Genomics (under revision).

Corresponding author: Fayaz Shaik Mahammad, Manipal Institute of Technology, Manipal Academy of
Higher Education, Manipal, Karnataka 576104, India (fayaz.shaik@manipal.edu).

## Contents

| Folder | What it holds |
|---|---|
| `ANALYSIS_GUIDE.ipynb` | Start here: the analyses in the order of the paper, with the scripts of each (in run order) and the deposited results that hold the reported numbers. It reads only files in this repository and needs only pandas. |
| `data/network/` | The harmonised, directed, signed, evidence-tiered regulatory network (587 nodes, 6,859 edges): `canonical_nodes.tsv`, `canonical_edges.tsv`, the four layer files, `edge_evidence_tier.tsv` (per-edge multiMiR evidence tier) and the identifier maps. |
| `data/pair_filter/` | The inputs and pair-level output of the miRNA–TF pair filter (Methods 2.1), with a README that describes each file and how its *P* value was computed. |
| `scripts/` | Every analysis script, in the layout in which it was run: the top level is the first analysis round and `v2/` … `v7/` are the later rounds. Where two rounds write a file of the same name, the later round's output is the one reported. `six_node_27b/` and `followups_27c/` hold the six-node analyses and their follow-ups in the layout of their analysis folders, with the run records (`S1/jobs.txt`, `F/F4/f4_jobs.txt`, `S7/s7c_jobs.txt`), the graphs and instance lists that later steps read, the per-parameter-set dynamics outputs (`S2/perset/`, `S6/perset/`) and the randomisation summaries; `perturbation_27d/` holds the tests with public perturbation data, and `perturbation_checks_28a/` the three later checks of those tests (series-level and leave-one-series-out analyses; Notes S1 and S2) with the five tables of the tests that they read, so that they run unchanged. `upstream_26series/` holds the earlier runs whose outputs the six-node scripts read (census graphs and null runs, the census and dynamics re-runs, the six-node pattern networks, the legacy-edge audit and the control-role audit), with their inputs, outputs and run records; its `MANIFEST_26series.tsv`, `MANIFEST_missing_scripts.tsv` and `MANIFEST_results.tsv` give for each file the script that reads or writes it and the reported values it backs. and `original_submission_code/` holds the network files and six R scripts of the original submission that `followups_27c/Q5` compares with the revised network (MIT licence, as in its `LICENSE`). |
| `results/` | The analysis outputs the paper reports: feed-forward-loop census and module membership, motif-significance randomisations, differential expression, enrichment (ORA and GSEA, with gene ratios, background ratios, raw and adjusted *p*), edge-level correlation and concordance, survival and meta-analysis tables, mediation, compartment (single-cell, spatial, xenograft, HPA), sequence-level and dynamical-model results, and the module-detection battery. `results/multiomics/` holds the mutation, copy-number, methylation, CPTAC, RPPA, DepMap, stromal-adjustment, external-cohort and metabolomics outputs, and `results/v4/` the metabolic-pathway analyses and the gene-set enrichment of the feed-forward-loop classes. Files above 50 MB are gzip-compressed (`.gz`). `results/from_analysis_machine/` holds files delivered by the analysis machine, among them `six_node_pattern_survival_all_circuits.csv` (every six-node circuit compared with its three-node core in Cox models: one row per circuit, scoring method and endpoint, with cohort-level counts only) and `six_node_pattern_os_circuits_q_below_0.05.csv` (the 16 circuits with *q* < 0.05 for overall survival; its column `source_line_in_s8_circuits` is the line of the circuit in the former file, counting the header as line 1). The files `results/from_analysis_machine/perturbation_tests_*` hold the tests with public perturbation experiments, purified-cell miRNA atlases and a second sequence-to-function model (Supplementary Note S1): every dataset found with the reason for any exclusion (`perturbation_tests_datasets.tsv`; `perturbation_tests_p3a_triage_note.txt` explains how the 40 miRNA contrasts were chosen), the settings fixed in advance with their logged changes (`perturbation_tests_settings.txt`) and the summary of each test. Paths cited inside these files (for example `P5/p5_summary.txt`) are those of the analysis folder; the files are deposited here with the prefix `perturbation_tests_`. The cell-atlas files are given as corrected (settings change 16) and, with the suffix `_as_delivered`, as first computed. In `perturbation_tests_settings.txt` the parties to the logged decisions are named "the authors"; the log is otherwise as written. `scripts/perturbation_27d/` holds the scripts that produced them. Besides public downloads, they read the network tables in `data/network/`, `results/v3/seqreg_ext_occlusion_allsites.csv`, TRRUST v2 and, for `P2/p2_analyse.py`, the composite six-node instances of the six-node analysis. The per-gene differential-expression tables of these tests are not included; the scripts regenerate them from GEO. |
| `logs/` | Run logs of the census, class-diversity and dynamics runs, in the layout of the analysis folder (`logs/`, `logs/v2/`, `logs/v3/`). |
| `supplementary_tables/` | Supplementary Tables S1–S9 as CSV (identical to the journal's Excel file). |
| `docs/files_cited_in_supplementary_notes.tsv` | Every file named in the Supplementary Notes, with its path in this repository. |
| `environment/` | R and Python session information. |
| `MANIFEST.tsv` | Every file with its size and MD5 checksum. |

## Data that are not redistributed

The analyses use public resources whose terms do not allow redistribution of patient-level or
licensed data, so only derived results are deposited here. To re-run the scripts, download:

- TCGA-BRCA gene and miRNA expression, clinical and survival data: UCSC Xena (https://xenabrowser.net/), TCGA Pan-Cancer clinical resource.
- CPTAC-BRCA proteome: Proteomic Data Commons. METABRIC: cBioPortal.
- GEO series: GSE176078, GSE19783, GSE96058, GSE42568, GSE45827, GSE10780, GSE25066, GSE20194, GSE41998, GSE22093, GSE23988, GSE32646, GSE42822, GSE22216, GSE37405, GSE59829, GSE78870, GSE28969, GSE73002.
- DepMap (CRISPR, CCLE): https://depmap.org. Human Protein Atlas (v25.1; v23.0 archive): https://www.proteinatlas.org.
- NCI Patient-Derived Models Repository: https://pdmr.cancer.gov. 10x Genomics public Visium breast sections.
- GWAS Catalog (v1.0.2, GRCh38): https://www.ebi.ac.uk/gwas.
- Interaction and gene-disease resources: DisGeNET v7.0, GeneCards, miRTarBase v9.0, miRWalk v3, TarBase, miRecords (via multiMiR), TRRUST v2, TransmiR v2.0, hTFtarget, TcoF-DB v2, dbCoRC, STRING v12, HMDD v4.0, miR2Disease, PhenomiR 2.0, miRBase v22.
- OncoKB and COSMIC gene lists (licensed; not included).

The miRNA–TF pair filter (cumulative hypergeometric statistic; Methods 2.1 of the paper) was run in the
first analysis round outside the scripts collected here. Its inputs and pair-level output are included in
`data/pair_filter/`. The filter's *P* value is the lower tail of the hypergeometric distribution, so the filter did
not select pairs that share more targets than expected, and the filter did not define the network, whose own files
are in `data/network/` (see `data/pair_filter/README.md`).

## Files to read with care

- The `sign` column of `data/network/canonical_edges.tsv` is the original deposit's placeholder (+1 for every edge that
  is not miRNA→target). Literature-derived signs are in the layer files and in Table S2 (`sign`, `trrust_mode`,
  `transmir_mode`).
- `results/ffl_census_summary.csv`, `results/ffl_census.csv` and `results/ffl_higher_order.csv.gz` are the first
  census, on the graph that still holds the 30 exemplar miRNA–miRNA edges. The census the paper reports is
  `results/v2/legacy_rerun_2026-09-26b/census_all_graphs.csv`.
- The rank columns of `results/v5/node_prioritisation_full.csv` place the four unrankable miRNAs first; the ranks
  the paper reports are those of Table S1, recomputed from its `priority` column in `ANALYSIS_GUIDE.ipynb`.
- `results/v5/tables/Table1_filtering_cascade.csv`, `Table2_network_composition.csv`, `Table2c_sign_provenance.csv`
  and `Table3_ffl_census.csv` are copies of the main-text tables, assembled from the files named in the notebook.
- The published Table S2 (`supplementary_tables/TableS2_all_interactions.csv`; `results/v5/tables/TableS1_all_interactions.csv`
  is the same table) differs from the output of `scripts/v5/21_tables.R` by three later corrections. Edges without an
  annotated sign have an empty `sign`, where the script copied the placeholder of `canonical_edges.tsv`. The 521 repression
  edges (136 TRRUST, 385 TransmiR) have `sign` -1, as their `trrust_mode` and `transmir_mode` state and as every analysis
  used. The `in_analysed_network` column, written by
  `scripts/upstream_26series/INBOX_2026-09-27/legacy_audit/tables/flag_tableS1.py`, is FALSE for the 30 exemplar
  miRNA–miRNA edges.
- `results/v7/mcode_19_headline_summary.csv` was compiled by hand from the files `mcode_00` to `mcode_18` of the same
  folder, which hold its values.
- No saved script wrote `results/v6/atac_encode_compartment_specificity_tests.csv`;
  `scripts/v6/atac_compartment_specificity_tests_rebuild.py` rebuilds it byte for byte from
  `results/v6/atac_encode_region_accessible.csv` and writes the rebuilt copy and a summary next to itself.

## Running the scripts

The analysis was run on Ubuntu 22.04.5 LTS with R 4.5.1 and Python 3.13.13 (package versions in
`environment/`). Scripts use absolute paths from the original project root, which appears in this copy
as the placeholder `/path/to/revision`, with inputs under `data/` and outputs under `results/`; a few
scripts read raw inputs from sibling folders of the project, shown as `/path/to/home/Desktop/DD/R_GPR/`.
Replace the placeholders with your own paths, or set up the same layout. Random seeds are fixed in the scripts; randomisation and bootstrap procedures use
1,000 iterations unless stated otherwise. The NCBI E-utilities scripts need your own e-mail address
in place of `your.email@example.org`.

Compiled programs are not included; build each from its `.c` source with gcc. `scripts/upstream_26series/INBOX_2026-09-26b/run_all.sh`
and `INBOX_2026-09-26/rerun/B_census/run_census.sh` contain their compile commands; the flags used for `scripts/ffl_enum.c`,
`ffl_enum2.c` and the `scripts/v2/` census programs were not recorded (`v2/v2_null.c` is the source that `run_all.sh` compiles
with `gcc -O2`). In `scripts/six_node_27b/`, the C programs were compiled with `gcc -O2` (`S7/v2_null.c`) and `gcc -O3 -march=native`
(`S7/ffl_census_composition.c`); for `S1/s1_null6.c` and `F/F4/f4_node_union.c`, `gcc -O2 -o <name> <name>.c -lm`
reproduces the programs used. The TCGA-BRCA R objects that `E/E6` and `S8` read (`data/brca_*.rds`) are built by
`scripts/01_tcga_brca_prep_de.R` and `scripts/10_survival_cox_hubs.R` from the downloads listed above;
`S7/s7b_restricted_gene_gene.py` reads the STRING v12 downloads and `followups_27c/Q7` the TargetScan 8.0 miRNA family
file. `E/E5/e5a_hypergeometric_filter.py` and `followups_27c/Q123`, `Q5` and `Q9` take the authors' filter archive as an
argument; its files are in `data/pair_filter/`. The two `S5` scripts read `results/v5/tables/TableS1_all_interactions.csv`
(the published Table S2). They were run on an earlier copy that differs only in the sign columns of 497 edges without an
annotated sign, which these scripts do not read.

These scripts address files by their paths under the original analysis root (`/path/to/revision`). `INBOX_2026-09-27b/`,
`INBOX_2026-09-27c/` and `INBOX_2026-09-27d/` are `scripts/six_node_27b/`, `scripts/followups_27c/` and
`scripts/perturbation_27d/`, and `INBOX_2026-09-28a/` is `scripts/perturbation_checks_28a/INBOX_2026-09-28a/`. `INBOX_2026-09-26/`,
`INBOX_2026-09-26b/`, `INBOX_2026-09-26c/`, `INBOX_2026-09-27/`, `HANDOVER_2026-09-27/figure_fixes/` and `cache/` are under
`scripts/upstream_26series/` with the same subpaths (`HANDOVER_2026-09-27/inboxes/INBOX_2026-09-27/` is its `INBOX_2026-09-27/`);
its MANIFEST files list the script that reads or writes each file, and its `scripts/05a_string_map_and_fetch.R` is an identical copy of `scripts/05a_string_map_and_fetch.R`. The authors'
earlier code folder (`/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR`) is `scripts/original_submission_code/`.

The other paths under the analysis root map as follows. `results/`, `scripts/` and `logs/` have the same relative paths here;
`results/from_analysis_machine/` holds further files delivered from the analysis folders. The network files that the scripts read from `data/` (`canonical_nodes.tsv`,
`canonical_edges.tsv`, the four layer files, `edge_evidence_tier.tsv`, `mirna_id_map.tsv` and `name_map.json`) are in
`data/network/`. The following are not deposited:
- `data/brca_*.rds`: patient-level TCGA data, built as described above.
- `data/ffl_module_sets.rds`: built by `scripts/11_ffl_module_membership.R`.
- `data/db/`: the TRRUST v2 and TransmiR v2.0 downloads.
- `cache/`: public downloads. These are the STRING v12 protein aliases and information files for human (9606); the
  TargetScan 8.0 miRNA family file (`miR_Family_Info.txt`); and, for the sequence-model tests, 600-kb hg38 windows
  around *COL1A1* and *COL3A1*, the Enformer and Borzoi target tables, UCSC RefSeq Select and the Borzoi weights
  (Hugging Face `johahi/borzoi-replicate-0` to `-3`). The exception is `cache/string_symbol_map_exact.tsv`, which is in
  `scripts/upstream_26series/`.

## Licence

Code: MIT licence (`LICENSE`). Derived data and result tables: CC BY 4.0.
