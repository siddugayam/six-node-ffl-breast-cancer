#!/bin/bash
cd /path/to/revision
for s in 11 22 33 44 55; do
  ./scripts/v2/v2_census2 results/v2/graph_pub_fine.txt 6 1 1 0.3 0.3 0.3 $s
done
