# Drug and crispr screens

## Scripts in this folder

| Script | What it does |
|---|---|
| `20_screens_demeter2.R` | DepMap RNAi (DEMETER2 v6 combined: Achilles + DRIVE + Marcotte) as an orthogonal loss-of-function modality to the CRISPR/Chronos analysis already done. |
| `21_screens_projectscore.R` | Project Score (Sanger CRISPR-Cas9, Score2 release: Sanger v2 + Broad 21Q2), a second independent CRISPR screen with its own cell-line panel and analysis pipeline. |
| `22_screens_drug_sensitivity.R` | Does a miR-29 target score or a collagen score predict sensitivity to any compound class in BREAST cell lines? |
| `22a_celline_scores.R` | 22a) Cell-line module scores used for the drug-sensitivity analysis (Part 1C). |
| `22b_node_drug_associations.R` | 22b) Per-NODE drug-sensitivity association, so that the integrated matrix (Part 1E) can carry a real per-node drug column rather than set membership. |
| `23_screens_orcs.py` | Published phenotypic CRISPR screens (BioGRID ORCS 2.0.18, human). |
| `23b_screens_literature.py` | Literature search for published CRISPR screens in breast models with migration, invasion, metastasis or ECM read-outs (BioGRID ORCS 2.0.18 has none). |
| `24_screens_integrated_matrix.R` | One integrated essentiality matrix over the 587 network nodes: CRISPR (DepMap 24Q4 Chronos, breast lines) RNAi (DEMETER2 v6 combined, breast lines) SCORE (Sanger Project Score 2, breast models) drug-sensitivity association (best … |
