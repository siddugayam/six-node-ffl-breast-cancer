<div align="center">

# Six-node feed-forward loops in the breast cancer regulatory network

**Code, regulatory network and result tables for**

*Deciphering the Transcriptional Regulation of Breast Cancer Genes through a Novel Six-Node Feed-Forward Loop Analysis*

Gayam Prasanna Kumar Reddy · Jesil Mathew A · Fayaz Shaik Mahammad

*Functional & Integrative Genomics* (under revision)

[![Code: MIT](https://img.shields.io/badge/code-MIT-2e6db4.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-4f9a6a.svg)](#licence)
[![R 4.5.1](https://img.shields.io/badge/R-4.5.1-276dc3.svg)](environment)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-3776ab.svg)](environment)

</div>

---

## At a glance

| | |
|---|---|
| **Regulatory network** | 587 nodes (157 TFs, 223 miRNAs, 207 genes); 6,829 directed edges in the analysed network |
| **Edge classes** | miRNA → gene 3,266 · miRNA → TF 1,553 · TF → miRNA 1,239 · TF → gene 681 · TF → TF 89 · gene–gene 1 |
| **Three-node FFL cores** | 1,649 typed cores: 1,434 composite, 206 miRNA-FFL and 9 TF-FFL (Table S4) |
| **Higher-order census** | 116,505 four-node, 1,506,263 five-node and 19,401,840 six-node modules (Table 3) |
| **Prioritised nodes** | 30: the ten highest-ranked TFs, genes and miRNAs (Table 4; all nodes in Table S1) |

<p align="center">
  <img src="docs/images/ffl_types.svg" width="900" alt="The three classes of three-node feed-forward loop (TF-FFL, miRNA-FFL, composite) and an n-node module">
</p>

<p align="center">
  <img src="docs/images/census_counts.svg" width="720" alt="Number of n-node feed-forward loop modules in the census graph, n = 3 to 7">
</p>

## Start here

1. **[`ANALYSIS_GUIDE.ipynb`](ANALYSIS_GUIDE.ipynb)** follows the paper section by section. For each analysis it lists the
   scripts in run order and loads the file behind each reported number. It reads only files in this repository and
   needs only pandas.
2. **[`data/network/`](data/network)** holds the regulatory network, and
   **[`supplementary_tables/`](supplementary_tables)** holds Tables S1–S9 as CSV (identical to the journal's Excel file).
3. Every folder has a README that explains its files.

## Repository map

| Folder | Contents |
|---|---|
| [`data/`](data) | The regulatory network ([`network/`](data/network)), the miRNA–TF pair filter ([`pair_filter/`](data/pair_filter)) and the curation rules ([`curation/`](data/curation)) |
| [`supplementary_tables/`](supplementary_tables) | Supplementary Tables S1–S9 |
| [`scripts/`](scripts) | The code of the main analysis, by round: the first round at the top level, then `v2/` … `v7/` |
| [`analyses/`](analyses) | Self-contained later analyses, each with its code, inputs and outputs: the census and motif null models, the six-node pattern, the perturbation-data tests and the re-runs on the analysed network |
| [`results/`](results) | The outputs of the scripts, in the same layout as `scripts/` |
| [`logs/`](logs) | Run logs of the census, class-diversity and dynamics runs |
| [`environment/`](environment) | R and Python session information |
| [`docs/`](docs) | Every file named in the Supplementary Notes with its path here, and the reconstructed run order |
| [`MANIFEST.tsv`](MANIFEST.tsv) | Every file with its size and MD5 checksum |

## Workflow

```mermaid
flowchart TD
    A["Regulatory network<br/>TRRUST · TransmiR · multiMiR · STRING<br/>587 nodes · 6,829 edges"]
    A --> B["Typed three-node FFL cores<br/>Table S4"]
    A --> C["Census of 3- to 7-node modules<br/>Table 3"]
    B --> D["Motif significance<br/>three null models"]
    C --> E["Six-node composite pattern<br/>Supplementary Note S2"]
    C --> I["Dynamics of<br/>higher-order modules"]
    B --> F["Node prioritisation<br/>Table 4, Table S1"]
    F --> G["Expression, protein and survival<br/>TCGA · GEO · CPTAC · METABRIC"]
    G --> H["Stromal mediation and<br/>compartment tests"]
    G --> K["Tests with public perturbation data<br/>Supplementary Note S1"]
    A --> J["Module detection<br/>Supplementary Note S7"]
```

## Reproducing the analyses

<details>
<summary><b>Software, paths and random seeds</b></summary>

<br>

The analysis was run on Ubuntu 22.04.5 LTS with R 4.5.1 and Python 3.13.13; package versions are in
[`environment/`](environment). Random seeds are fixed in the scripts, and randomisation and bootstrap procedures use
1,000 iterations unless stated otherwise.

Scripts use absolute paths under the original analysis root, shown as the placeholder `/path/to/revision`. Set it to
the root of this repository: `scripts/`, `results/`, `analyses/` and `logs/` then resolve as they are here. Two
differences remain. The network files that the scripts read from `data/` are in `data/network/`, and a few scripts read
raw inputs from sibling folders of the project, shown as `/path/to/home/Desktop/DD/R_GPR/`. Where two rounds write a file
of the same name, the later round's output is the one reported. The NCBI E-utilities scripts need your own e-mail
address in place of `your.email@example.org`.

</details>

<details>
<summary><b>Compiled programs</b></summary>

<br>

The C programs are given as source; build each with gcc. The recorded commands are:

| Program | Command |
|---|---|
| `analyses/census_and_motif_nulls/census/ffl_census_composition.c` and `nulls/v2_null.c` | in `analyses/census_and_motif_nulls/run_all.sh` (`gcc -O3 -march=native` and `gcc -O2`) |
| `analyses/census_and_motif_nulls/first_rerun/ffl_census_composition.c` | in `first_rerun/run_census.sh` |
| `analyses/six_node_pattern/S7/v2_null.c` | `gcc -O2` |
| `analyses/six_node_pattern/S7/ffl_census_composition.c` | `gcc -O3 -march=native` |
| `analyses/six_node_pattern/S1/s1_null6.c`, `F/F4/f4_node_union.c` | `gcc -O2 -o <name> <name>.c -lm` |

The flags used for `scripts/ffl_enum.c`, `scripts/ffl_enum2.c` and the census programs of `scripts/v2/` were not
recorded; `scripts/v2/v2_null.c` is the program that `run_all.sh` compiles with `gcc -O2`.

</details>

<details>
<summary><b>Inputs of particular scripts</b></summary>

<br>

- The TCGA-BRCA R objects that `analyses/six_node_pattern/E/E6` and `S8` read (`data/brca_*.rds`) are built by
  `scripts/01_tcga_brca_prep_de.R` and `scripts/10_survival_cox_hubs.R` from the downloads listed below.
- `analyses/six_node_pattern/S7/s7b_restricted_gene_gene.py` reads the STRING v12 downloads, and
  `analyses/six_node_followups/Q7` the TargetScan 8.0 miRNA family file.
- `analyses/six_node_pattern/E/E5/e5a_hypergeometric_filter.py` and `analyses/six_node_followups/Q123`, `Q5` and `Q9`
  take the authors' pair-filter archive as an argument; its files are in `data/pair_filter/`.
- The two scripts of `analyses/six_node_pattern/S5` read `results/v5/tables/TableS1_all_interactions.csv` (the published
  Table S2). They were run on an earlier copy that differs only in the sign columns of 497 edges without an annotated
  sign, which these scripts do not read.
- The perturbation-data tests read, besides public downloads, the network tables, `results/v3/seqreg_ext_occlusion_allsites.csv`,
  TRRUST v2 and, for `P2/p2_analyse.py`, the six-node composite instances of `analyses/six_node_pattern/S1`.

</details>

<details>
<summary><b>Data that are not redistributed</b></summary>

<br>

The analyses use public resources whose terms do not allow redistribution of patient-level or licensed data, so only
derived results are deposited. To re-run the scripts, download:

- TCGA-BRCA gene and miRNA expression, clinical and survival data: UCSC Xena (https://xenabrowser.net/), TCGA Pan-Cancer clinical resource.
- CPTAC-BRCA proteome: Proteomic Data Commons. METABRIC: cBioPortal.
- GEO series: GSE176078, GSE19783, GSE96058, GSE42568, GSE45827, GSE10780, GSE25066, GSE20194, GSE41998, GSE22093, GSE23988, GSE32646, GSE42822, GSE22216, GSE37405, GSE59829, GSE78870, GSE28969, GSE73002.
- DepMap (CRISPR, CCLE): https://depmap.org. Human Protein Atlas (v25.1; v23.0 archive): https://www.proteinatlas.org.
- NCI Patient-Derived Models Repository: https://pdmr.cancer.gov. 10x Genomics public Visium breast sections.
- GWAS Catalog (v1.0.2, GRCh38): https://www.ebi.ac.uk/gwas.
- Interaction and gene-disease resources: DisGeNET v7.0, GeneCards, miRTarBase v9.0, miRWalk v3, TarBase, miRecords (via multiMiR), TRRUST v2, TransmiR v2.0, hTFtarget, TcoF-DB v2, dbCoRC, STRING v12, HMDD v4.0, miR2Disease, PhenomiR 2.0, miRBase v22.
- OncoKB and COSMIC gene lists (licensed; not included).

Not deposited from the analysis root:

- `data/brca_*.rds`: patient-level TCGA data, built as described above.
- `data/ffl_module_sets.rds`: built by `scripts/11_ffl_module_membership.R`.
- `data/db/`: the TRRUST v2 and TransmiR v2.0 downloads.
- `cache/`: public downloads. These are the STRING v12 protein aliases and information files for human (9606); the
  TargetScan 8.0 miRNA family file (`miR_Family_Info.txt`); and, for the sequence-model tests, 600-kb hg38 windows around
  *COL1A1* and *COL3A1*, the Enformer and Borzoi target tables, UCSC RefSeq Select and the Borzoi weights (Hugging Face
  `johahi/borzoi-replicate-0` to `-3`). The exception is `string_symbol_map_exact.tsv`, which is in `data/network/`.

</details>

<details>
<summary><b>The miRNA–TF pair filter</b></summary>

<br>

The miRNA–TF pair filter (cumulative hypergeometric statistic; Methods 2.1 of the paper) was run in the first analysis
round outside the scripts collected here. Its inputs and pair-level output are in [`data/pair_filter/`](data/pair_filter).
The filter's *P* value is the lower tail of the hypergeometric distribution, so the filter did not select pairs that
share more targets than expected, and the filter did not define the network, whose own files are in `data/network/`
(see `data/pair_filter/README.md`).

</details>

<details>
<summary><b>Files to read with care</b></summary>

<br>

- The `sign` column of `data/network/canonical_edges.tsv` is the original deposit's placeholder (+1 for every edge that
  is not miRNA→target). Literature-derived signs are in the layer files and in Table S2 (`sign`, `trrust_mode`,
  `transmir_mode`).
- `results/ffl_census_summary.csv`, `results/ffl_census.csv` and `results/ffl_higher_order.csv.gz` are the first census,
  on the graph that still holds the 30 exemplar miRNA–miRNA edges. The census the paper reports is
  `results/v2/census_rerun/census_all_graphs.csv`.
- The rank columns of `results/v5/node_prioritisation_full.csv` place the four unrankable miRNAs first; the ranks the
  paper reports are those of Table S1, recomputed from its `priority` column in `ANALYSIS_GUIDE.ipynb`.
- `results/v5/tables/` holds the tables under the analysis numbering; its README maps them to the published tables.
- The published Table S2 (`supplementary_tables/TableS2_all_interactions.csv`) differs from the output of
  `scripts/v5/21_tables.R` by three later corrections. Edges without an annotated sign have an empty `sign`, where the
  script copied the placeholder of `canonical_edges.tsv`. The 521 repression edges (136 TRRUST, 385 TransmiR) have
  `sign` −1, as their `trrust_mode` and `transmir_mode` state and as every analysis used. The `in_analysed_network`
  column, written by `analyses/analysed_network_reruns/tables/flag_tableS1.py`, is FALSE for the 30 exemplar
  miRNA–miRNA edges.
- `results/v7/mcode_19_headline_summary.csv` was compiled by hand from the files `mcode_00` to `mcode_18` of the same
  folder, which hold its values.
- No saved script wrote `results/v6/atac_encode_compartment_specificity_tests.csv`;
  `scripts/v6/atac_compartment_specificity_tests_rebuild.py` rebuilds it byte for byte from
  `results/v6/atac_encode_region_accessible.csv` and writes the rebuilt copy and a summary next to itself.

</details>

## Citation

If you use this code or these data, please cite:

> Gayam Prasanna Kumar Reddy, Jesil Mathew A, Fayaz Shaik Mahammad. Deciphering the Transcriptional Regulation of Breast
> Cancer Genes through a Novel Six-Node Feed-Forward Loop Analysis. *Functional & Integrative Genomics* (under revision).

## Licence

Code: MIT licence ([`LICENSE`](LICENSE)). Derived data and result tables: CC BY 4.0.

## Contact

Corresponding author: Fayaz Shaik Mahammad, Manipal Institute of Technology, Manipal Academy of Higher Education,
Manipal, Karnataka 576104, India (fayaz.shaik@manipal.edu).
