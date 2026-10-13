#!/bin/bash
# RAND-ESU census at seven nodes of the published-network census graph with census_09_esu_classmask, five seeds.
cd /path/to/six-node-ffl-breast-cancer
for s in 11 22 33 44 55; do
  ./scripts/03_ffl_census/census_09_esu_classmask results/expression_validation/graph_pub_fine.txt 7 1 1 0.2 0.2 0.2 0.2 $s
done
