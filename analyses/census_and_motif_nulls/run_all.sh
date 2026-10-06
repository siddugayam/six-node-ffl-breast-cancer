#!/bin/bash
# analyses/census_and_motif_nulls re-runs.  Project read-only; all outputs here.
set -euo pipefail
cd "$(dirname "$0")"
gcc -O3 -march=native -o census/ffl_census_composition census/ffl_census_composition.c
gcc -O2 -o nulls/v2_null nulls/v2_null.c
: > jobs.txt
# ---- motif nulls: 1,000 randomisations, seed 20250908, 100 swaps per edge (the stored runs' settings)
for g in dep pub dep_nolegacy pub_nolegacy pub_nostring pub_nolegacy_nostring; do
  for m in A B C; do
    echo "nulls/v2_null graphs/null_$g.txt $m 1000 20250908 100 > nulls/null_${m}_${g}.tsv 2> nulls/null_${m}_${g}.err" >> jobs.txt
  done
done
# ---- census: exhaustive n = 3-5 (single partition), n = 6 (64 partitions), n = 7 RAND-ESU (8 seeds)
for v in table3 nolegacy table3_nostring nolegacy_nostring; do
  G="graphs/census_${v}_fine.txt graphs/census_${v}_R4.txt graphs/census_${v}_orig.txt"
  mkdir -p census/$v/n6_parts census/$v/n7_seeds
  for k in 3 4 5; do echo "census/ffl_census_composition $G $k 1 0 42 1 1 1 1 1 1 > census/$v/n$k.txt" >> jobs.txt; done
  for p in $(seq 0 63); do echo "census/ffl_census_composition $G 6 64 $p 42 1 1 1 1 1 1 > census/$v/n6_parts/part_$p.txt" >> jobs.txt; done
  for s in 11 22 33 44 55; do echo "census/ffl_census_composition $G 7 1 0 $s 1 1 0.2 0.2 0.2 0.2 > census/$v/n7_seeds/seed_$s.txt" >> jobs.txt; done
  for s in 71 72 73;       do echo "census/ffl_census_composition $G 7 1 0 $s 1 1 1 0.1 0.15 0.15 > census/$v/n7_seeds/seed_$s.txt" >> jobs.txt; done
done
xargs -P 11 -I{} sh -c "{}" < jobs.txt
# ---- aggregate
for v in table3 nolegacy table3_nostring nolegacy_nostring; do
  for k in 3 4 5; do python3 census/aggregate_census.py census/$v/census_n$k census/$v/n$k.txt > /dev/null; done
  python3 census/aggregate_census.py census/$v/census_n6 census/$v/n6_parts/part_*.txt > /dev/null
  for s in 11 22 33 44 55 71 72 73; do python3 census/aggregate_census.py census/$v/n7_seeds/census_n7_seed_$s census/$v/n7_seeds/seed_$s.txt > /dev/null; done
done
echo ALL_DONE
