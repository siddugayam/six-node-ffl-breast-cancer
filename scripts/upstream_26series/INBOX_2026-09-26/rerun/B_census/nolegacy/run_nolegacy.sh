#!/bin/bash
# Same enumerator on the no-legacy-miRNA graph (scripts/v7 build_graph, 8,003 arcs).
set -euo pipefail
cd "$(dirname "$0")"
BIN=../ffl_census_composition
F=graph_nolegacy_fine.txt; R=graph_nolegacy_R4.txt; O=graph_nolegacy_orig.txt
for k in 3 4 5; do $BIN $F $R $O $k 1 0 42 1 1 1 1 1 1 > nl_n$k.txt; python3 ../aggregate_census.py census_nolegacy_n$k nl_n$k.txt > /dev/null; done
mkdir -p n6_parts
seq 0 63 | xargs -P 3 -I{} sh -c "$BIN $F $R $O 6 64 {} 42 1 1 1 1 1 1 > n6_parts/part_{}.txt"
python3 ../aggregate_census.py census_nolegacy_n6_exhaustive n6_parts/part_*.txt > /dev/null
