# The twelve networks of Bhat et al. (2024), analysed with the settings of the paper (Supplementary Note S2)

The miRNA-FFL, TF-FFL and composite patterns of Bhat et al. (2024) at three to six nodes, rebuilt on the census graph as
twelve networks (`data/ffl_networks/`). Each analysis of the paper is repeated on them with the paper's own settings.

- `SETTINGS.md`: every setting, written before its run, and every change made later, with its time.
- `RUN_LOG.tsv`: every run, with its start and end time, exit status and the md5 of `SETTINGS.md` at its start. The logged md5 values are those of the file as written; the copy here differs from it only in folder and script names, which were changed to those of this repository. The md5 at the start of the A1b certification and analysis steps differs from that at the start of the A1b module runs for one reason: on 29 Sep at 16:34, before any A1b result existed, the name of the nested-check script in the A1b list of scripts was corrected from `a1b_nested_compare.py` to `a1b_nested_check.py`, with a dated note on that line.
- `DOWNLOADS.tsv`: every file fetched from outside the analysis machine.

| Folder | Section | What it holds |
|---|---|---|
| `topology/` | T | Composition, degree and hubs of the twelve networks; overlap with the typed cores; D1–D4 modules |
| `nulls/` | O | Over-representation under NULL-B, NULL-C and NULL-L, and the one-layer variants (1,000 randomisations each) |
| `enrichment/` | E | GSEA against the tumour-versus-normal signature, the network-background test and over-representation analysis |
| `survival/` | S | Each four- to six-node circuit against its three-node core (TCGA-BRCA, overall and progression-free interval) |
| `concordance/` | X | Sign concordance in TCGA-BRCA and CPTAC against matched nulls |
| `prioritisation/` | P | The prioritisation with FFL participation counted from these networks; full rankings of the 583 rankable nodes |
| `knocktf/` | A4 | KnockTF 2.0 support for the TF–TF arcs of the six-node patterns |
| `dynamics/` | A1, A1b | The dynamical model of the main text with a TF2→TF1 arc added (A1), the comparison of each passing sign configuration with the seven smaller modules of the layer factorial (A1b), the nested checks against the stored modules, the sign configurations, per-parameter-set results and prevalences |
| `cluster/` | A2 | Cluster sensitivity: without miR-17~92 instances, and with each co-transcribed cluster counted once |
| `evidence/` | A3 | Evidence sensitivity: physical STRING links only, and validated miRNA targets only |
| `scripts/` | | The code of all sections |

`scripts/run_logged.py` runs one step and appends it to `RUN_LOG.tsv`. `scripts/bhat_null.c` is compiled by `nulls.py`,
`cluster_nulls.py` and `evidence_nulls.py` with `gcc -O2`. Every script names its inputs in its header. The scripts
read the analysis root (`/path/to/revision`, the root of this repository), the twelve networks, the original submission's
SIF files, the stored runs of `analyses/six_node_pattern/` and the KnockTF tables of `analyses/perturbation_tests/`.
