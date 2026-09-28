# Feed-forward-loop networks, one file set per class and size

Twelve networks: three FFL classes (miRNA-FFL, TF-FFL, composite FFL) at three to six nodes. They are the revised
counterparts of the six SIF files of the original submission (`analyses/original_submission_code/SIF_files/`).
`make_ffl_networks.py` writes every file here; `docs/MANIFEST.tsv` lists their md5.

## Files

For each network `<name>` (for example `6node_composite_FFL`):

- `<name>.sif`: one line per edge, `source <tab> interacts with <tab> target`, the format of the original files.
  Undirected links (gene-gene, miRNA pair) appear once, with the two names sorted.
- `<name>_edges.tsv`: `source`, `target`, `edge_type`, `sign`, `sign_source`, `evidence_tier`, `layer`, `n_ffl`.
- `<name>_nodes.tsv`: `node`, `type` (TF, miRNA or Gene), `n_ffl`.

`n_ffl` is the number of FFL instances of that network that contain the edge or node (defined below).

## Graph

All twelve networks are built on the census graph (the graph of Table 3), 9,226 arcs:
- the analysed network: the interaction table (Table S2) without the 30 miRNA-miRNA edges of the original circuits,
  587 nodes and 6,829 edges;
- 604 TF-to-target arcs from TRRUST;
- 833 gene-gene arcs from STRING;
- 960 arcs joining 480 pairs of miRNAs that lie within 10 kb of each other.

The `layer` column of every edge file names the source of each edge. The 3-node networks use only edges of the
analysed network.

## Definition: Bhat et al. (2024), at every size

Each size follows a fixed pattern of compulsory edges. M1 and M2 are miRNAs, T1 and T2 are TFs, G1 and G2 are genes,
and all are distinct nodes:

| size | compulsory edges |
|---|---|
| 3 | M1-T1 link, M1 -> G1, T1 -> G1 |
| 4 | M1-T1 link, M1 -> G1, T1 -> G2, G1 - G2 |
| 5 | M1-T1 link, M1 -> G1, M2 -> G1, M1 - M2, T1 -> G2, G1 - G2 |
| 6 | M1-T1 link, M1 -> G1, M2 -> G1, M1 - M2, T1-T2 link, T2 -> G2, G1 - G2 |

The class sets the direction of the links:
- **miRNA-FFL**: M1 -> T1, and at six nodes also T1 -> T2.
- **TF-FFL**: T1 -> M1, and at six nodes also T2 -> T1.
- **Composite FFL**: both directions.

G1 - G2 (gene-gene) and M1 - M2 (miRNA pair) are undirected. Extra edges among the nodes are allowed. An instance is
one assignment of nodes to the roles.

The classes overlap: every composite instance is also a miRNA-FFL and a TF-FFL instance. In this network almost every
TF -> miRNA edge is matched by a miRNA -> TF edge (1,223 of 1,239), so most TF-FFL and many miRNA-FFL instances are
composite. The six-node composite FFL is the pattern of Note S2.

## Counts

| network | FFL instances | distinct node sets | nodes (TF / miRNA / gene) | edges | edges by type |
|---|---|---|---|---|---|
| 3node_miRNA_FFL | 1,355 | 1,355 | 364 (83 / 153 / 128) | 1,760 | TF_target 351, miRNA_target 1,409 |
| 3node_TF_FFL | 1,182 | 1,182 | 348 (76 / 148 / 124) | 1,576 | TF_miRNA 556, TF_target 334, miRNA_target 686 |
| 3node_composite_FFL | 1,173 | 1,173 | 346 (76 / 146 / 124) | 2,112 | TF_miRNA 551, TF_target 331, miRNA_target 1,230 |
| 4node_miRNA_FFL | 2,580 | 2,554 | 365 (93 / 166 / 106) | 2,036 | TF_target 295, gene_gene 141, miRNA_target 1,600 |
| 4node_TF_FFL | 2,090 | 2,065 | 352 (88 / 159 / 105) | 1,769 | TF_miRNA 652, TF_target 282, gene_gene 140, miRNA_target 695 |
| 4node_composite_FFL | 2,054 | 2,029 | 350 (88 / 157 / 105) | 2,389 | TF_miRNA 642, TF_target 282, gene_gene 140, miRNA_target 1,325 |
| 5node_miRNA_FFL | 2,272 | 1,456 | 246 (66 / 87 / 93) | 1,335 | TF_target 223, gene_gene 116, miRNA_miRNA 112, miRNA_target 884 |
| 5node_TF_FFL | 1,719 | 1,168 | 239 (61 / 85 / 93) | 1,191 | TF_miRNA 339, TF_target 210, gene_gene 116, miRNA_miRNA 102, miRNA_target 424 |
| 5node_composite_FFL | 1,719 | 1,168 | 239 (61 / 85 / 93) | 1,530 | TF_miRNA 339, TF_target 210, gene_gene 116, miRNA_miRNA 102, miRNA_target 763 |
| 6node_miRNA_FFL | 14,003 | 8,606 | 330 (122 / 96 / 112) | 2,329 | TF_target 720, gene_gene 142, miRNA_miRNA 142, miRNA_target 1,325 |
| 6node_TF_FFL | 15,701 | 9,950 | 328 (117 / 99 / 112) | 2,203 | TF_miRNA 543, TF_target 694, gene_gene 141, miRNA_miRNA 169, miRNA_target 656 |
| 6node_composite_FFL | 3,898 | 2,498 | 210 (45 / 79 / 86) | 1,488 | TF_miRNA 290, TF_target 266, gene_gene 113, miRNA_miRNA 96, miRNA_target 723 |

Edges by source (the `layer` column):

| network | analysed network | TRRUST layer | STRING layer | 10 kb miRNA-pair layer |
|---|---|---|---|---|
| 3node_miRNA_FFL | 1,760 | 0 | 0 | 0 |
| 3node_TF_FFL | 1,576 | 0 | 0 | 0 |
| 3node_composite_FFL | 2,112 | 0 | 0 | 0 |
| 4node_miRNA_FFL | 1,896 | 0 | 140 | 0 |
| 4node_TF_FFL | 1,630 | 0 | 139 | 0 |
| 4node_composite_FFL | 2,250 | 0 | 139 | 0 |
| 5node_miRNA_FFL | 1,108 | 0 | 115 | 112 |
| 5node_TF_FFL | 974 | 0 | 115 | 102 |
| 5node_composite_FFL | 1,313 | 0 | 115 | 102 |
| 6node_miRNA_FFL | 1,707 | 339 | 141 | 142 |
| 6node_TF_FFL | 1,555 | 339 | 140 | 169 |
| 6node_composite_FFL | 1,193 | 87 | 112 | 96 |

## How these counts relate to the paper

- **Note S2:** 3,898 six-node composite instances. The file reproduces this, and the
  script stops otherwise.
- **Table S4 and Table 3, three-node row:** the paper counts 1,649 typed three-node cores with exclusive classes:
  1,434 composite, 206 miRNA-FFL and 9 TF-FFL. There, a composite core counts only as composite.
  - 285 of those cores have a TF as target; the Bhat pattern requires a gene target. The other 1,364 have a gene target:
    1,173 composite, 182 miRNA-FFL and 9 TF-FFL.
  - The 3-node files are exactly these cores, with overlapping classes:
    - composite: 1,173 = 1,173;
    - TF-FFL: 1,182 = 1,173 + 9;
    - miRNA-FFL: 1,355 = 1,173 + 182.
  - The script checks this correspondence core by core, and stops otherwise.
- **Table 3, other counts:** 5,833 at n = 3 and the modules at four to seven nodes are D1-D4 modules counted
  irrespective of node type. They are a different object from these typed networks.

## The original submission's files

Counts are taken from the files as deposited. Many miRNAs occur twice with different capitalisation (`hsa-miR-` and
`hsa-mir-`), so nodes are given both as written and with capitalisation merged.

| file | lines | distinct edges | nodes as written | nodes, capitalisation merged | revised counterpart |
|---|---|---|---|---|---|
| 3-miR.sif | 6,236 | 5,493 | 782 | 575 | 3node_miRNA_FFL |
| 3-TF.sif | 5,741 | 5,360 | 783 | 575 | 3node_TF_FFL |
| 3-Comp.sif | 6,526 | 6,526 | 770 | 567 | 3node_composite_FFL |
| 4-TF.sif | 18 | 18 | 11 | 11 | 4node_* (exemplar circuit) |
| 5-TF.sif | 41 | 41 | 14 | 14 | 5node_* (exemplar circuit) |
| 6-TF.sif | 30 | 30 | 10 | 10 | 6node_* (exemplar circuit) |

The original files are of a different kind from the revised ones:
- The original 3-node files are the full networks of the original analysis (5,360 to 6,526 distinct edges each).
- The original 4-, 5- and 6-node files are single exemplar circuits. The 5- and 6-node files hold 20 and 10
  miRNA-miRNA edges: the 30 edges the revision excludes.
- Each revised file contains only the edges that the FFL instances of its class use.

## Signs

`sign` and `sign_source` follow the published interaction table:
- miRNA -> target: -1 (mechanistic repression).
- TF edges: +1 or -1 from the TRRUST or TransmiR mode (Activation or Repression).
- All others: blank, with sign_source `unannotated (excluded from sign-dependent analyses)`.

Census-layer arcs: TRRUST arcs follow the same rule; STRING and miRNA-pair links are unsigned. `evidence_tier` is
the table's tier for miRNA -> target edges. When the table read is the published one, the script also checks that
these values equal its `sign` and `sign_source` columns for every edge, and stops otherwise.

## Rebuild

`python3 make_ffl_networks.py <analysis root> <folder of the original SIF files> <output folder>`. The script reads:
- `data/canonical_nodes.tsv` and `data/canonical_edges.tsv`;
- the three layer files `data/layer_TF_target.tsv`, `data/layer_gene_gene.tsv` and `data/layer_miRNA_miRNA.tsv`;
- the interaction table `results/v5/tables/TableS1_all_interactions.csv`;
- the typed-core table `results/v5/tables/TableS2_ffl_cores.csv`, for the 3-node check.
