# S1

## Scripts in this folder

| Script | What it does |
|---|---|
| `make_s1_inputs.py` | S1 inputs. Reads the project read-only; writes only into this folder. |
| `s1_cotranscription_overlap.py` | S1 context: do the co-transcribed miRNA pairs (10 kb layer) share targets and TF inputs more than other miRNA pairs? |
| `s1_null6.c` | S1 of analyses/six_node_pattern: are six-node patterns over-represented? |
| `s1_summarise.py` | S1 summary (analyses/six_node_pattern). Reads runs/null_{B,C,L,LT,LG,LM}.tsv (s1_null6) and obs/ (observed). |
