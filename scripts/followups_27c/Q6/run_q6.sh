#!/bin/bash
# Q6: S1's exact command (INBOX_2026-09-27b/S1/jobs.txt: R 1000, seed 20250908, 100 swaps per edge, census_reps 1000,
# nproc 1, proc 0) on the three variant graphs; NULL-C and NULL-L, plus OBS.  One process per job.
cd "$(dirname "$0")"
B=/path/to/revision/INBOX_2026-09-27b
BIN=$B/S1/s1_null6
declare -A G=( [retyped22]="graph_retyped22_reclassed.txt $B/S1/graph_nolegacy_labels.txt"
               [validated]="$B/S7/null_validated.txt $B/S7/labels_validated.txt"
               [physical900]="$B/S7/graph_physical_ge_0.900.txt $B/S7/labels_physical_ge_0.900.txt" )
for k in retyped22 validated physical900; do
  set -- ${G[$k]}
  $BIN $1 $2 OBS 0 20250908 100 1000 1 0 > runs/obs_$k.tsv 2> runs/obs_$k.err
  for M in C L; do
    ( $BIN $1 $2 $M 1000 20250908 100 1000 1 0 > runs/null_${M}_$k.tsv 2> runs/null_${M}_$k.err; echo "DONE $k $M $(date +%H:%M:%S)" >> runs/q6_driver.log ) &
  done
done
wait
echo "Q6_ALL_DONE $(date)" >> runs/q6_driver.log
