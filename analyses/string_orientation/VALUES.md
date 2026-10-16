# VALUES: orientation of the STRING associations in the census graph

Lines count the header as line 1. Tables are written by `string_orientation_graphs.py` (graphs) and
`summarise_string_orientation.py` (counts from `counts/`, run by `run_string_orientation.sh`). The census graph
built here reproduces `census_noexemplar_fine.txt` arc for arc; the published counts at n = 3–6 and the published
n = 7 samples seed by seed are reproduced (the script stops otherwise).

## (a) Random orientation of the 941 associations (`string_bootstrap_summary.csv`)

Columns replicates, published_count, median, p2_5, p97_5, mean, sd, cv, published_percentile (where the published
count lies in the bootstrap distribution), recovery_median / p2_5 / p97_5 (share of the published modules found
again), jaccard_median / p2_5 / p97_5.

| line | n | replicates | published | median | 95 % interval | mean | SD | CV | published percentile | recovery median (95 %) | Jaccard median (95 %) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | 500 | 5,833 | 5,804 | 5,592–6,001 | 5,803.4 | 104.3 | 0.0180 | 60.1 | 0.8730 (0.8547–0.8936) | 0.7786 (0.7544–0.8029) |
| 3 | 4 | 500 | 116,505 | 115,632 | 111,604–120,003 | 115,689.4 | 2,120.9 | 0.0183 | 65.2 | 0.8398 (0.8240–0.8598) | 0.7284 (0.7060–0.7522) |
| 4 | 5 | 500 | 1,506,263 | 1,473,116 | 1,402,622–1,555,410 | 1,474,100.2 | 40,449.6 | 0.0274 | 78.6 | 0.7279 (0.7033–0.7582) | 0.5812 (0.5527–0.6107) |
| 5 | 6 | 500 | 19,401,840 | 18,536,388 | 17,382,861–19,943,487 | 18,566,338.7 | 682,996.9 | 0.0368 | 88.8 | 0.6154 (0.5832–0.6547) | 0.4587 (0.4270–0.4933) |

STRING arcs entering the graph per replicate (`string_bootstrap_graphs.csv`, column string_arcs_added): 799–839 (the listed orientation gives 833). Per-replicate counts: `string_bootstrap_counts.csv`.

## (b) Gene–gene layer variants (`gene_gene_layer_variants.csv`, column modules; sd for n = 7)

| line | graph | n | method | modules | SD |
|---|---|---|---|---|---|
| 2 | census graph (STRING in the listed orientation) | 3 | exhaustive | 5,833 |  |
| 3 | census graph (STRING in the listed orientation) | 4 | exhaustive | 116,505 |  |
| 4 | census graph (STRING in the listed orientation) | 5 | exhaustive | 1,506,263 |  |
| 5 | census graph (STRING in the listed orientation) | 6 | exhaustive | 19,401,840 |  |
| 6 | census graph (STRING in the listed orientation) | 7 | RAND-ESU, 8 seeds | 225,761,602 | 939186.5 |
| 7 | TRRUST-only gene-gene links (census graph minus STRING arcs) | 3 | exhaustive | 3,432 |  |
| 8 | TRRUST-only gene-gene links (census graph minus STRING arcs) | 4 | exhaustive | 81,104 |  |
| 9 | TRRUST-only gene-gene links (census graph minus STRING arcs) | 5 | exhaustive | 791,085 |  |
| 10 | TRRUST-only gene-gene links (census graph minus STRING arcs) | 6 | exhaustive | 8,751,860 |  |
| 11 | TRRUST-only gene-gene links (census graph minus STRING arcs) | 7 | RAND-ESU, 8 seeds | 87,313,549 | 477782.8 |
| 12 | no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs) | 3 | exhaustive | 1,774 |  |
| 13 | no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs) | 4 | exhaustive | 57,494 |  |
| 14 | no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs) | 5 | exhaustive | 232,488 |  |
| 15 | no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs) | 6 | exhaustive | 2,467,228 |  |
| 16 | no gene-gene layer (minus STRING and the 604 TRRUST-layer arcs) | 7 | RAND-ESU, 8 seeds | 15,473,009 | 146704.1 |

## (c) Modules that traverse no STRING arc (`orientation_robust_core.csv`)

no_string_arc: no STRING arc in the module's induced subgraph. no_string_pair: no two of its vertices form
a STRING association, so the module is unchanged under every orientation. n = 7: mean over the 8 RAND-ESU seeds.

| line | n | method | modules | no STRING arc (fraction) | no STRING pair (fraction) |
|---|---|---|---|---|---|
| 2 | 3 | exhaustive | 5,833 | 2,767 (0.4744) | 2,342 (0.4015) |
| 3 | 4 | exhaustive | 116,505 | 66,634 (0.5719) | 61,889 (0.5312) |
| 4 | 5 | exhaustive | 1,506,263 | 524,397 (0.3481) | 422,344 (0.2804) |
| 5 | 6 | exhaustive | 19,401,840 | 4,711,285 (0.2428) | 3,568,213 (0.1839) |
| 6 | 7 | RAND-ESU, 8 seeds | 225,761,602 | 36,586,055 (0.1621) | 26,093,870 (0.1156) |

