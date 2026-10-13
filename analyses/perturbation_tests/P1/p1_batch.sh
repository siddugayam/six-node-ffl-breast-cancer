#!/bin/bash
# run p1_run_dataset.py over a list file: part<TAB>arms<TAB>gse<TAB>tf<TAB>cell<TAB>tag<TAB>extra-args
cd "$(dirname "$0")/.."
while IFS=$'\t' read -r part arms gse tf cell tag extra; do
  [ -z "$gse" ] && continue
  echo "=== $gse $tf $cell $tag $(date +%H:%M:%S)"
  python3 P1/p1_run_dataset.py "$part" "$arms" "$gse" "$tf" "$cell" --tag "$tag" $extra 2>&1 | tail -2 | cut -c1-400
done < "$1"
echo BATCH_DONE
