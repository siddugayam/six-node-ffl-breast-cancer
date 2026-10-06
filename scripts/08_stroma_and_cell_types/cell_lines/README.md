# Cell lines

Two independent implementations of the stroma-free test, which asks whether the collagen and miR-29 associations hold
in breast cancer cell lines, where no stroma is present: the Python series `celllines_01` … `celllines_12` and the R
series `stroma_free_01` … `stroma_free_07` (DepMap Public 24Q4). Each wrote its own deposited outputs in `results/v2/`.
The R files that start with `_celllines_` are written and run by the Python scripts that hold their code.

## Scripts in this folder

| Script | What it does |
|---|---|
| `_celllines_full_pull.R` | Written and run by celllines_11_purity_axis.py, which holds this code; running that script regenerates this file. |
| `_celllines_full_write.R` | Written and run by celllines_11_purity_axis.py, which holds this code; running that script regenerates this file. |
| `_celllines_tcga_pull.R` | Written and run by celllines_02_corr.py, which holds this code; running that script regenerates this file. |
| `_celllines_tcga_targets.R` | Written and run by celllines_09_network_sign_tcga.py, which holds this code; running that script regenerates this file. |
| `celllines_01_prep.py` | Build the cell-line data objects for the "stroma-free test". |
| `celllines_02_corr.py` | The stroma-free test, parts B and C. |
| `celllines_03_controls.py` | Controls for celllines_02_corr.py: assay positive controls for the CCLE miRNA panel, dynamic range, network-level sign tests against a permutation null, epithelial–mesenchymal confounding and absolute expression by lineage. |
| `celllines_04_crispr.py` | CRISPR dependency of the network TFs against collagen expression in breast cancer lines (part D of the stroma-free test). |
| `celllines_05_drug.py` | Collagen and miR-29-target module scores against drug sensitivity in breast cancer lines (part E of the stroma-free test; exploratory). |
| `celllines_06_sc_epithelium.py` | Control for part C: are breast cancer CELL LINES collagen-negative because they are cultured, or because tumour epithelial cells are collagen-negative in vivo too? |
| `celllines_07_drug_partial.py` | Follow-up to part E. Every GDSC2 hit for the collagen module score is also a hit for a generic mesenchymal score (the two scores correlate rho = 0.75 across the 69 breast cancer lines). |
| `celllines_08_power_and_tissue_gradient.py` | Two things needed to read the cell-line null correctly. |
| `celllines_09_network_sign_tcga.py` | The same network-level sign test that celllines_03_controls.py ran in 50 breast cancer cell lines, run in 1,066 paired TCGA-BRCA tumours, so the two are directly comparable: over ALL canonical targets of a miRNA, what fraction … |
| `celllines_10_mir29_context.py` | Where does miR-29 sit relative to the stromal / mesenchymal axis in each system? - TCGA-BRCA bulk: miR-29a/b/c vs the 4-marker CAF score and vs COL1A1/COL3A1 - CCLE breast lines: miR-29a/b/c vs the mesenchymal score and vs … |
| `celllines_11_purity_axis.py` | Decisive test of whether the bulk miR-29 -> collagen correlation is anything more than a tumour-composition (purity) effect. |
| `celllines_12_summary.py` | Consolidate the stroma-free test into one table, plus the essentiality re-check on ALL breast-lineage lines that have a CRISPR screen (not only the ones that also have expression). |
| `stroma_free_01_prep.R` | Cell-line "stroma-free test": data preparation DepMap Public 24Q4 (figshare article 27993248), md5-verified. |
| `stroma_free_02_expr_corr.R` | STROMA-FREE TEST -- part B & C B) TF->collagen and miR-29->collagen correlations in breast cancer lines C) absolute collagen expression: cell lines vs bulk tumours vs fibroblasts Data: DepMap Public 24Q4 (figshare 27993248; md5 … |
| `stroma_free_03_mirna.R` | STROMA-FREE TEST -- part B (miRNA arm) + sensitivity analyses CCLE miRNA = Nanostring nCounter panel, CCLE_miRNA_20181103.gct (734 miRNAs x 954 lines), from data.broadinstitute.org/ccle/ |
| `stroma_free_04_controls_nulls.R` | STROMA-FREE TEST -- power, positive controls, matched nulls, pan-cancer |
| `stroma_free_05_crispr.R` | STROMA-FREE TEST -- part D: dependency x collagen-expression relationships CRISPRGeneEffect.csv = DepMap Public 24Q4 (md5 6edf7ade09b9b34199210b559d4745d3) |
| `stroma_free_06_drug.R` | STROMA-FREE TEST -- part E (EXPLORATORY): does a collagen / miR-29-target module score predict drug sensitivity in breast cancer cell lines? |
| `stroma_free_07_consolidate.R` | STROMA-FREE TEST -- consolidation, EMT adjustment, summary table |
