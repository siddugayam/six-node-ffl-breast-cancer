# Over representation

## Scripts in this folder

| Script | What it does |
|---|---|
| `09_enrichment_ora.R` | Over-representation analysis (ORA) for every FFL motif class, with an EXPLICIT background universe, full statistics, and a low-support flag. |
| `09b_claimed_terms_hypergeometric.R` | Hypergeometric tests, without filtering, of each named functional term for every motif set, with both background universes. |
| `09c_exemplar_enrichment.R` | GO (cellular component, biological process) and KEGG enrichment of the five protein-coding members of the six-node exemplar circuit against the tested universe; counts the terms supported by fewer than three genes. |
