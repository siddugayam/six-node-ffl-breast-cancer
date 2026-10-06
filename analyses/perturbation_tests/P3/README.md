# P3

## Scripts in this folder

| Script | What it does |
|---|---|
| `p3a_analyse.py` | P3a: validate the analysed network's miRNA->target edges by evidence tier with miRNA gain/loss datasets (SETTINGS.md, P3). |
| `p3a_arms.py` | P3a arms: the contrasts read from the P3a hits that reach the 40 largest (SETTINGS P3; see P3a/p3a_triage_note.txt). |
| `p3a_prep_counts.py` | P3a: GSM-labelled count matrices from the submitters' files for the P3a series without NCBI-generated counts. |
| `p3a_search.py` | P3a search: the P1b search pattern applied to every miRNA node of the analysed network (223), Homo sapiens, expression profiling; then the first-pass screen (samples naming a network miRNA near a perturbation word). |
| `p3a_series_sensitivity.py` | P3a sensitivity analysis (descriptive): 25 of the 40 contrasts come from one series (GSE115646), so contrasts are first averaged within each series. |
| `p3b_analyse.py` | P3b: validate the analysed network's TF->target edges with KnockTF 2.0 human knockdown datasets (SETTINGS.md, P3). |
| `p3b_datasets.py` | P3b: KnockTF 2.0 human datasets for the network TFs (sources of the analysed network's TF->target edges), with the number of control and treated samples read from each dataset's KnockTF page. |
| `p3b_fetch_missing.py` | P3b: KnockTF's bulk file (knocktf_v2_main_human.txt) holds only the DataSet_01/02 series and 4 DataSet_03 series. |
