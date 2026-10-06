#!/bin/bash
# RAND-ESU census at seven nodes of the published-network census graph with census_09_esu_classmask, three further seeds
# with a second set of sampling probabilities.
cd /path/to/revision
for s in 71 72 73; do ./scripts/03_ffl_census/census_09_esu_classmask results/v2/graph_pub_fine.txt 7 1 1 1 0.1 0.15 0.15 $s; done
