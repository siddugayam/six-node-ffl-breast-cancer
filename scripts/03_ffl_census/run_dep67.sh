#!/bin/bash
cd /path/to/revision
./scripts/03_ffl_census/v2_census2 results/v2/graph_dep_fine.txt 6 1 1 1 1 1 42
for s in 11 22 33; do ./scripts/03_ffl_census/v2_census2 results/v2/graph_dep_fine.txt 7 1 1 0.3 0.3 0.3 0.3 $s; done
