#!/bin/bash
# the topological definition census of the census graph under four treatments of reciprocal TF<->miRNA pairs (see SETTINGS.md):
# base (contracted, the published census), unc (uncontracted), scc (graph-level SCC condensation), mod (module-level
# condensation). Exhaustive at n = 3-6 (n = 6 in 64 partitions of the root vertices), RAND-ESU at n = 7 with the 8 seeds
# and keep-probabilities of analyses/census_and_motif_nulls/run_all.sh. The condensation DAG is enumerated exhaustively
# at n = 3-7.
# usage: FFL_ROOT=<repository checkout> bash run_reciprocal_variants.sh [parallel jobs, default 12]
set -euo pipefail
cd "$(dirname "$0")"
: "${FFL_ROOT:?set FFL_ROOT to a checkout of the repository}"
P=${1:-12}
G="$FFL_ROOT/analyses/census_and_motif_nulls/graphs"
mkdir -p work/bin work/keys counts
gcc -O3 -march=native -o work/bin/crv census_reciprocal_variants.c
python3 make_scc_file.py > counts/make_scc_file.log
ARGS="$G/census_noexemplar_fine.txt $G/census_noexemplar_uncontracted.txt work/scc_uncontracted.txt"
J=work/jobs.txt; : > $J
for m in base unc mod; do
  mkdir -p counts/$m work/keys/$m
  for k in 3 4 5; do echo "work/bin/crv $m $ARGS $k 1 0 42 1 1 1 1 1 1 work/keys/$m/n$k > counts/$m/n$k.txt" >> $J; done
  for p in $(seq 0 63); do echo "work/bin/crv $m $ARGS 6 64 $p 42 1 1 1 1 1 1 work/keys/$m/n6_part_$p > counts/$m/n6_part_$p.txt" >> $J; done
  for s in 11 22 33 44 55; do echo "work/bin/crv $m $ARGS 7 1 0 $s 1 1 0.2 0.2 0.2 0.2 - > counts/$m/n7_seed_$s.txt" >> $J; done
  for s in 71 72 73;       do echo "work/bin/crv $m $ARGS 7 1 0 $s 1 1 1 0.1 0.15 0.15 - > counts/$m/n7_seed_$s.txt" >> $J; done
done
mkdir -p counts/scc work/keys/scc
for k in 3 4 5 6; do echo "work/bin/crv scc $ARGS $k 1 0 42 1 1 1 1 1 1 work/keys/scc/n$k > counts/scc/n$k.txt" >> $J; done
echo "work/bin/crv scc $ARGS 7 1 0 42 1 1 1 1 1 1 - > counts/scc/n7.txt" >> $J
# longest jobs first
grep " mod " $J > $J.sorted || true; grep -v " mod " $J >> $J.sorted
xargs -P "$P" -I{} sh -c "{}" < $J.sorted
echo ALL_DONE
