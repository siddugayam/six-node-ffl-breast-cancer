# Scripts

The code of the main analysis, in the layout in which it was run. Scripts address their inputs and outputs by these
paths, so the folders keep the names of the analysis rounds. `ANALYSIS_GUIDE.ipynb` (repository root) lists the scripts
of each analysis in the order of the paper. The later, self-contained analyses are in `analyses/`.

| Folder | What it holds |
|---|---|
| top level | First round, numbered in run order: network assembly and its layers (`00`–`07`), the three-node census (`02`, `03`, `ffl_*`), motif significance and node sets (`08`, `30`), enrichment, hubs and survival (`09`–`13`, `31`), external cohorts (`14`–`18`) and CPTAC and RPPA (`20`–`24`). Letter prefixes: `A_` DepMap essentiality, `B_` druggability, `C1`–`C4` ChIP-seq and LINCS, `D_` miR-29 literature, `E1`–`E4` metabolic hypotheses, `ct_` cell-type and stromal analyses, `mo_` mutation, copy number and methylation. |
| `v2/` | Coherence and the census re-analysed; differential expression; CAF-adjusted correlations and mediation; subtype and pan-cancer analyses; cell lines (`celllines_*`, `v2_cl_*`); external cohorts (`N*`). |
| `v3/` | Dynamical models of the feed-forward loops; architecture, controllability and information flow; sequence-level regulation (`seqreg_*`); deconvolution and mediation; single-cell and screen analyses. |
| `v4/` | Gene-set enrichment of the FFL classes; metabolic analyses; survival cohorts (`M*`) and web-portal cross-checks (`portal_*`). |
| `v5/` | Node prioritisation, node compendium and literature counts; supplementary tables; reactive stroma and neoadjuvant response (Farmer signatures). |
| `v6/` | Spatial transcriptomics, patient-derived xenografts, Human Protein Atlas, ENCODE accessibility (`atac_*`), miRNA survival meta-analysis and liquid biopsy (`c*`), GWAS loci and figures. |
| `v7/` | Module detection: community detection, MCODE, BioNet, jActiveModules and GO over-representation. |

The C programs are given as source (`.c`); compile them with gcc before running the scripts that call them (see the
repository README).
