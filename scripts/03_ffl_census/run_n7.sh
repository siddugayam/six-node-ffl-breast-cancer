#!/bin/bash
cd /path/to/revision
for s in 11 22 33 44 55; do
  ./scripts/03_ffl_census/v2_census2 results/v2/graph_pub_fine.txt 7 1 1 0.2 0.2 0.2 0.2 $s
done
