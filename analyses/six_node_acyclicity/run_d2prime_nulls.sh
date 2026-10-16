#!/bin/bash
# D2' six-node counts on the census graph (observed) and under NULL-C and NULL-L (see SETTINGS.md): s1_null6.c with the
# added D2' count, seed 20250908, 100 swaps or trades per item, R replicates per model split over P processes (each
# process randomises every replicate and counts those with r % P == its index, as s1_null6.c does).
# usage: FFL_ROOT=<repository checkout> bash run_d2prime_nulls.sh [replicates, default 1000] [processes per model, default 4]
set -euo pipefail
cd "$(dirname "$0")"
: "${FFL_ROOT:?set FFL_ROOT to a checkout of the repository}"
R=${1:-1000}; P=${2:-4}
S1="$FFL_ROOT/analyses/six_node_pattern/S1"
mkdir -p work/bin runs
gcc -O2 -o work/bin/d2p six_node_d2prime_nulls.c
work/bin/d2p "$S1/graph_noexemplar_null.txt" "$S1/graph_noexemplar_labels.txt" OBS 0 20250908 100 1 1 0 > runs/obs.tsv 2> runs/obs.err
J=work/jobs.txt; : > $J
for m in C L; do for p in $(seq 0 $((P-1))); do
  echo "work/bin/d2p $S1/graph_noexemplar_null.txt $S1/graph_noexemplar_labels.txt $m $R 20250908 100 $R $P $p > runs/null_${m}_part$p.tsv 2> runs/null_${m}_part$p.err" >> $J
done; done
xargs -P $((2*P)) -I{} sh -c "{}" < $J
echo ALL_DONE
