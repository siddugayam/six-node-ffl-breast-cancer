#!/bin/bash
# Census of the deposited-network census graph (results/v2/graph_dep_fine.txt) with the C program census_09_esu_classmask:
# exhaustive at six nodes, RAND-ESU sampling with three seeds at seven nodes.
cd /path/to/revision
./scripts/03_ffl_census/census_09_esu_classmask results/v2/graph_dep_fine.txt 6 1 1 1 1 1 42
for s in 11 22 33; do ./scripts/03_ffl_census/census_09_esu_classmask results/v2/graph_dep_fine.txt 7 1 1 0.3 0.3 0.3 0.3 $s; done
