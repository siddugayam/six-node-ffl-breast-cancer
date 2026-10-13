# Census

## Scripts in this folder

| Script | What it does |
|---|---|
| `aggregate_census.py` | Aggregate ffl_census_composition partition outputs into one summary (JSON + readable text). |
| `ffl_census_composition.c` | Counting-only n-node FFL census (C), used for the census re-runs. |
| `summarise_census.py` | One table across the four census graphs: exact counts n = 3-6, RAND-ESU n = 7 (8 seeds), maximum edge classes, class sets at the maximum, share with 2 TF + 2 miRNA + 2 Gene, and the count of modules containing a miRNA-miRNA arc. |
