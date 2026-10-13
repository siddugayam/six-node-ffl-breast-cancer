#!/bin/bash
# RAND-ESU census at six nodes of the published-network census graph (results/expression_validation/graph_pub_fine.txt) with
# census_09_esu_classmask, five seeds.
cd /path/to/six-node-ffl-breast-cancer
for s in 11 22 33 44 55; do
  ./scripts/03_ffl_census/census_09_esu_classmask results/expression_validation/graph_pub_fine.txt 6 1 1 0.3 0.3 0.3 $s
done
