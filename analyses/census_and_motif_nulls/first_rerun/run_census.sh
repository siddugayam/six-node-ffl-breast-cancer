#!/bin/bash
# Part B re-runs for the first re-runs.  Reads the v2 census graphs read-only; writes only here.
set -euo pipefail
cd "$(dirname "$0")"
G=/path/to/revision/results/v2
BIN=./ffl_census_composition
gcc -O3 -march=native -o $BIN ffl_census_composition.c

# validation against logs/v2/census2_n{3,4,5}.log (exhaustive, single partition)
for k in 3 4 5; do
  $BIN $G/graph_pub_fine.txt $G/R4_graph_aug.txt $G/orig_pub.txt $k 1 0 42 1 1 1 1 1 1 > val_n$k.txt
done

# n = 6, exhaustive, 64 disjoint source partitions on 10 workers
mkdir -p n6_parts
seq 0 63 | xargs -P 10 -I{} sh -c "$BIN $G/graph_pub_fine.txt $G/R4_graph_aug.txt $G/orig_pub.txt 6 64 {} 42 1 1 1 1 1 1 > n6_parts/part_{}.txt"
python3 aggregate_census.py census_n6_exhaustive n6_parts/part_*.txt > /dev/null

# n = 7, RAND-ESU with the eight seeds and q-profiles of scripts/03_ffl_census/run_n7.sh and run_n7b.sh
mkdir -p n7_seeds
for s in 11 22 33 44 55; do echo "$s 1 1 0.2 0.2 0.2 0.2"; done >  n7_jobs.txt
for s in 71 72 73;       do echo "$s 1 1 1 0.1 0.15 0.15"; done >> n7_jobs.txt
cat n7_jobs.txt | xargs -P 8 -L 1 sh -c "$BIN $G/graph_pub_fine.txt $G/R4_graph_aug.txt $G/orig_pub.txt 7 1 0 \$0 \$1 \$2 \$3 \$4 \$5 \$6 > n7_seeds/seed_\$0.txt"
for s in 11 22 33 44 55 71 72 73; do python3 aggregate_census.py n7_seeds/census_n7_seed_$s n7_seeds/seed_$s.txt > /dev/null; done
