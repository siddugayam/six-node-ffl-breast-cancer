# Protein atlas

## Scripts in this folder

| Script | What it does |
|---|---|
| `01_parse_hpa_xml.py` | Parse HPA v25.1 per-gene XML entries for the 26 genes of interest. |
| `02_parse_hpa_sc_and_tissue.py` | HPA v25.1 per-gene XML: single-cell-type RNA (nCPM), consensus tissue RNA (nTPM), and cell-type specificity calls, for the 26 genes. |
| `04_build_hpa_tables.py` | Assemble the Human Protein Atlas tables (tissue validation). |
| `09_compartment_table.py` | The compartment test at protein level. One row per gene, combining HPA IHC annotator summaries, normal-breast/soft-tissue cell-type staining, breast-cancer tumour-cell staining and its recorded subcellular location, and the HPA … |
