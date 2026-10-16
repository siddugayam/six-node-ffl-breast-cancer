#!/bin/bash
# Three-node FFL motif test under NULL-A/B/C and the tier- and sign-restricted null models (see SETTINGS.md):
# 8 models x 2 graphs (analysed network, census graph) x 2 typings (published, 22 co-regulators retyped), 1,000
# replicates each, seed 20250908, 100 swaps or trades per arc (the published settings).
# usage: FFL_ROOT=<repository checkout> bash run_tier_sign_nulls.sh [parallel jobs, default 4]
set -euo pipefail
cd "$(dirname "$0")"
: "${FFL_ROOT:?set FFL_ROOT to a checkout of the repository}"
P=${1:-4}
mkdir -p work/bin runs
gcc -O2 -o work/bin/mnts motif_null_tier_sign.c
python3 make_tier_sign_inputs.py > runs/make_tier_sign_inputs.log
J=work/jobs.txt; : > $J
for g in analysed census; do for t in published retyped; do for m in A B C BT BS BTS CT CTS; do
  echo "work/bin/mnts work/${g}_${t}.txt $m 1000 20250908 100 work/${g}_${t}.attr > runs/null_${m}_${g}_${t}.tsv 2> runs/null_${m}_${g}_${t}.err" >> $J
done; done; done
xargs -P "$P" -I{} sh -c "{}" < $J
echo ALL_DONE
