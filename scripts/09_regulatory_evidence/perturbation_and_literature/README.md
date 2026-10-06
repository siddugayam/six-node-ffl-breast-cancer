# Perturbation and literature

## Scripts in this folder

| Script | What it does |
|---|---|
| `C1_lincs_fetch.py` | C-i) SigCom LINCS: does knocking down / over-expressing NFKB1, RELA, SP1, ETS1 change COL1A1 / COL3A1 in the L1000 CD-coefficient signatures? |
| `C2_chipseq_libs.py` | C-ii/iii) Is COL1A1 / COL3A1 a target of NFKB1, RELA, SP1 or ETS1 according to ChIP-seq consensus (ChEA/ENCODE), curated TRRUST, ARCHS4 co-expression, GEO TF perturbation signatures, and LINCS L1000 CRISPR-KO consensus up/down … |
| `C3_perturbation_summary.R` | C) Assemble perturbation evidence for TF -> COL1A1 / COL3A1 |
| `C4_lincs_mcf7.R` | LINCS transcription-factor perturbation signatures in breast cancer cell lines: change of COL1A1 and COL3A1 per TF and perturbation type; writes results/multiomics/lincs_TF_collagen_breastlines.csv. |
| `D2_mir29_evidence_table.py` | D) Curated PubMed literature table: miR-29 perturbation -> collagen genes. |
| `D_mir29_literature.py` | D) Literature evidence (PubMed) for miR-29 family perturbation -> collagen genes. |
