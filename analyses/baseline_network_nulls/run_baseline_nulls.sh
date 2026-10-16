#!/bin/bash
# Three-node FFL motif test on the baseline networks (see SETTINGS.md): the degree-preserving, reciprocity-preserving
# and bipartite-layer null models of motif_null_models.c (program arguments A, B and C), 1,000 randomised networks
# each, seed 20250908, 100 swaps or trades per arc; two networks x two typings.
# usage: FFL_ROOT=<repository checkout> bash run_baseline_nulls.sh [parallel jobs, default 4]
set -euo pipefail
cd "$(dirname "$0")"
: "${FFL_ROOT:?set FFL_ROOT to a checkout of the repository}"
P=${1:-4}
mkdir -p work runs
gcc -O2 -o work/mnm motif_null_models_large.c
Rscript build_baseline_network.R > runs/build_baseline_network.log 2>&1
J=work/jobs.txt; : > $J
for n in validated validated_plus_predicted; do for t in published retyped; do
  for m in A:degree_preserving B:reciprocity_preserving C:bipartite_layer; do
    a=${m%%:*}; o=runs/null_${m#*:}_${n}_${t}
    echo "work/mnm network_${n}_${t}.txt $a 1000 20250908 100 > $o.tsv 2> $o.err" >> $J
  done
done; done
xargs -P "$P" -I{} sh -c "{}" < $J
echo ALL_DONE
