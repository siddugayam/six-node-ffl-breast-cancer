#!/bin/bash
cd /path/to/revision
for s in 71 72 73; do ./scripts/v2/v2_census2 results/v2/graph_pub_fine.txt 7 1 1 1 0.1 0.15 0.15 $s; done
