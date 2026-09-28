# miRNA–TF pair filter (Methods 2.1)

`hypergeometric_pair_filter.py` recomputes the filter from `data/pair_filter/value.txt`: the lower-tail hypergeometric
*P* in a 256-gene universe, Benjamini–Hochberg at 5% over all 2,676 rows. It reproduces the retained pairs of the authors'
workbook (1,566 pairs, 1,547 of which share no target) and writes `results/pair_filter/hypergeometric_pair_filter.tsv`.
Standard library only: `python3 scripts/02_pair_filter/hypergeometric_pair_filter.py`.
`data/pair_filter/README.md` describes the inputs; `analyses/six_node_followups/Q123`, `Q5` and `Q9` compare the
filter with the network.
