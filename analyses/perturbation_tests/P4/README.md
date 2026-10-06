# P4

## Scripts in this folder

| Script | What it does |
|---|---|
| `p4_borzoi_occlusion.py` | P4 of analyses/perturbation_tests: replicate the Enformer all-site promoter occlusion with Borzoi. |
| `p4_download_weights.py` | P4: fetch the Borzoi replicate weights 1-3 (borzoi-pytorch port of the official Calico human weights) into the 27d cache, and record each file in DOWNLOADS.tsv. |
| `p4_family_matrices.py` | P4 checks (a) and (b): the matrices of each motif family at each collagen promoter, and the COL1A1 calls ranked 9–18 by Enformer effect with their position relative to the exon-1 splice donor. |
| `p4_summarise.py` | P4 summary: applies the revised decision rule (SETTINGS.md, change 1) to the deposited Enformer values (sanity check) and, when the Borzoi scores exist, to Borzoi. |
