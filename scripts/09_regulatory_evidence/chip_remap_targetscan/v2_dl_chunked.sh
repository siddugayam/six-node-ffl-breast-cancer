#!/bin/bash
# N-way parallel HTTP range download (targetscan.org throttles single connections hard).
URL="$1"; OUT="$2"; N="${3:-12}"
LEN=$(curl -sSI "$URL" | awk 'tolower($1)=="content-length:"{print $2+0}')
echo "total bytes: $LEN in $N chunks"
CH=$(( (LEN + N - 1) / N ))
for i in $(seq 0 $((N-1))); do
  S=$((i*CH)); E=$((S+CH-1)); [ $E -ge $LEN ] && E=$((LEN-1))
  curl -sS --retry 5 -r ${S}-${E} -o "${OUT}.part${i}" "$URL" &
done
wait
cat $(for i in $(seq 0 $((N-1))); do echo "${OUT}.part${i}"; done) > "$OUT"
rm -f ${OUT}.part*
echo "expected $LEN got $(stat -c%s "$OUT")"
[ "$(stat -c%s "$OUT")" = "$LEN" ] && echo SIZE_OK || echo SIZE_MISMATCH
