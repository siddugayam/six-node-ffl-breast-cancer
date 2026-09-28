#!/usr/bin/env python3
"""Download ChIP-Atlas 'Target Genes' tables (hg38) for every TF that is a source
of a TF_target edge in the canonical network.  Tables give, per target gene,
the average MACS2 -10*log10(p) score across all ChIP-seq experiments for that TF,
plus the per-experiment score (0 = no peak in the window).
Distances: 1 = TSS +/-1 kb, 5 = +/-5 kb, 10 = +/-10 kb."""
import csv, os, sys, time, urllib.request, urllib.error, json
from concurrent.futures import ThreadPoolExecutor

BASE = "/path/to/revision"
OUT  = os.path.join(BASE, "cache/direct/chipatlas/target")
os.makedirs(OUT, exist_ok=True)

rows = list(csv.DictReader(open(os.path.join(BASE,"data/canonical_edges.tsv")), delimiter="\t"))
tfs = sorted({r["source"] for r in rows if r["edge_type"] == "TF_target"})
sys.stderr.write("TFs with TF_target edges: %d\n" % len(tfs))

def fetch(args):
    tf, dist = args
    path = os.path.join(OUT, "%s.%s.tsv" % (tf, dist))
    if os.path.exists(path) and os.path.getsize(path) > 100:
        return (tf, dist, "cached", os.path.getsize(path))
    url = "https://chip-atlas.dbcls.jp/data/hg38/target/%s.%s.tsv" % (tf, dist)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                data = r.read()
            with open(path, "wb") as fh:
                fh.write(data)
            return (tf, dist, "ok", len(data))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return (tf, dist, "404", 0)
            time.sleep(3)
        except Exception:
            time.sleep(5)
    return (tf, dist, "fail", 0)

jobs = [(tf, d) for tf in tfs for d in ("1", "10")]
res = []
with ThreadPoolExecutor(max_workers=6) as ex:
    for i, r in enumerate(ex.map(fetch, jobs)):
        res.append(r)
        if (i+1) % 40 == 0:
            sys.stderr.write("  %d/%d\n" % (i+1, len(jobs))); sys.stderr.flush()

from collections import Counter
sys.stderr.write("status: %s\n" % Counter(r[2] for r in res))
json.dump([{"tf":a,"dist":b,"status":c,"bytes":d} for a,b,c,d in res],
          open(os.path.join(BASE,"cache/direct/chipatlas/fetch_status.json"),"w"), indent=1)
missing = sorted({r[0] for r in res if r[2] != "ok" and r[2] != "cached"})
sys.stderr.write("TFs with any failure/404: %d -> %s\n" % (len(missing), missing[:40]))
