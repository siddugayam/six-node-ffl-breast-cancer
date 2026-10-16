#!/bin/bash
# Census counts for the STRING orientation analysis (see SETTINGS.md).
# Phase 1: the census graph as published (listed orientation; module keys kept; STRING arc and pair counts) and the
#          two gene-gene layer variants, exhaustive at n = 3-6 and RAND-ESU at n = 7 (published seeds).
# Phase 2: the 500 re-oriented graphs, exhaustive at n = 3-6, each module looked up among the published modules.
# Jobs whose output file is already complete are skipped, so the script can be restarted.
# usage: FFL_ROOT=<repository checkout> bash run_string_orientation.sh [parallel jobs, default 5] [replicates, default 500]
set -euo pipefail
cd "$(dirname "$0")"
: "${FFL_ROOT:?set FFL_ROOT to a checkout of the repository}"
P=${1:-5}; R=${2:-500}
mkdir -p work/bin work/keys counts/listed counts/nostring counts/nogenegene counts/boot
gcc -O3 -march=native -o work/bin/csv census_string_variants.c
[ -f work/graphs/boot_$(printf %03d $((R-1))).txt ] || python3 string_orientation_graphs.py "$R" > counts/string_orientation_graphs.log
PAIRS=work/string_pairs.txt
J=work/jobs_phase1.txt; : > $J
for g in listed nostring nogenegene; do
  A=$PAIRS; S=work/string_arcs_listed.txt; K=work/keys/$g
  [ $g = listed ] || { A=-; S=-; K=-; }
  for k in 3 4 5; do echo "work/bin/csv work/graphs/$g.txt $k 1 0 42 1 1 1 1 1 1 $A $S - $( [ $K = - ] && echo - || echo ${K}_n$k ) > counts/$g/n$k.txt" >> $J; done
  for p in $(seq 0 63); do echo "work/bin/csv work/graphs/$g.txt 6 64 $p 42 1 1 1 1 1 1 $A $S - $( [ $K = - ] && echo - || echo ${K}_n6_part_$p ) > counts/$g/n6_part_$p.txt" >> $J; done
  for s in 11 22 33 44 55; do echo "work/bin/csv work/graphs/$g.txt 7 1 0 $s 1 1 0.2 0.2 0.2 0.2 $A $S - - > counts/$g/n7_seed_$s.txt" >> $J; done
  for s in 71 72 73;       do echo "work/bin/csv work/graphs/$g.txt 7 1 0 $s 1 1 1 0.1 0.15 0.15 $A $S - - > counts/$g/n7_seed_$s.txt" >> $J; done
done
todo(){ while read -r line; do out=${line##*> }; grep -q "^$1" "$out" 2>/dev/null || echo "$line"; done; }
todo no_string_pair < $J > $J.todo
xargs -P "$P" -I{} sh -c "{}" < $J.todo
# sorted key files of the published modules, one per n
[ -f work/keys/listed_n6_sorted.bin ] || python3 - <<'EOF'
import glob, numpy as np
for k in (3, 4, 5, 6):
    fs = sorted(glob.glob(f"work/keys/listed_n{k}.keys")) if k < 6 else sorted(glob.glob("work/keys/listed_n6_part_*.keys"))
    a = np.concatenate([np.fromfile(f, dtype=np.uint64) for f in fs]); u = np.unique(a)
    assert len(u) == len(a), "duplicate module keys"
    u.tofile(f"work/keys/listed_n{k}_sorted.bin"); print(f"n={k}: {len(u)} module keys")
EOF
J=work/jobs_phase2.txt; : > $J
for r in $(seq 0 $((R-1))); do
  b=$(printf %03d $r)
  echo "work/bin/csv work/graphs/boot_$b.txt 6 1 0 42 1 1 1 1 1 1 - - work/keys/listed_n6_sorted.bin - > counts/boot/boot_${b}_n6.txt" >> $J
  for k in 3 4 5; do echo "work/bin/csv work/graphs/boot_$b.txt $k 1 0 42 1 1 1 1 1 1 - - work/keys/listed_n${k}_sorted.bin - > counts/boot/boot_${b}_n$k.txt" >> $J; done
done
todo in_baseline < $J > $J.todo
xargs -P "$P" -I{} sh -c "{}" < $J.todo
echo ALL_DONE
