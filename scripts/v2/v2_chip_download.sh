#!/bin/bash
# Download ChIP-Atlas hg38 "Target Genes" tables for a list of TF symbols at a given
# TSS-distance threshold (1, 5 or 10 kb). Records HTTP status for every request.
SP="$1"      # output dir
LIST="$2"    # file of TF symbols, one per line
THR="$3"     # 1 | 5 | 10
mkdir -p "$SP/$THR"
fetch() {
  tf="$1"; sp="$2"; thr="$3"
  out="$sp/$thr/${tf}.${thr}.tsv"
  if [ -s "$out" ]; then echo -e "${tf}\t${thr}\tCACHED\t$(stat -c%s "$out")"; return; fi
  code=$(curl -sS -w '%{http_code}' --max-time 300 -o "$out.part" \
     "https://chip-atlas.dbcls.jp/data/hg38/target/${tf}.${thr}.tsv")
  sz=$(stat -c%s "$out.part" 2>/dev/null || echo 0)
  if [ "$code" = "200" ] && [ "$sz" -gt 100 ]; then mv "$out.part" "$out"; else rm -f "$out.part"; fi
  echo -e "${tf}\t${thr}\t${code}\t${sz}"
}
export -f fetch
cat "$LIST" | xargs -P 8 -I{} bash -c 'fetch "$@"' _ {} "$SP" "$THR"
