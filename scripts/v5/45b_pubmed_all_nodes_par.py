#!/usr/bin/env python3
"""PubMed record counts for all 587 network nodes (threaded, globally rate-limited to <3 req/s),
to test whether network centrality tracks how heavily a node has been studied."""
import json, time, threading, queue, urllib.parse, urllib.request, csv, os

EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "&tool=ffl_revision&email=your.email@example.org"
ROOT = "/path/to/revision"
OUT = os.path.join(ROOT, "results/v5/network_literature_volume.csv")
LOCK = threading.Lock(); LAST = [0.0]; WLOCK = threading.Lock()


def throttle():
    with LOCK:
        dt = time.time() - LAST[0]
        if dt < 0.36:
            time.sleep(0.36 - dt)
        LAST[0] = time.time()


def get(u, tries=5):
    for i in range(tries):
        throttle()
        try:
            with urllib.request.urlopen(u, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(1.5 * (i + 1))


def count(term):
    return int(get(EUT + "esearch.fcgi?db=pubmed&retmode=json&retmax=0&term=%s%s"
                   % (urllib.parse.quote(term), TOOL))["esearchresult"]["count"])


nodes = []
with open(os.path.join(ROOT, "data/canonical_nodes.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        nodes.append((r["name"], r["type"]))

done = set()
if os.path.exists(OUT):
    with open(OUT) as f:
        for r in csv.DictReader(f):
            if r.get("name"):
                done.add(r["name"])
todo = [(n, t) for n, t in nodes if n not in done]
print(f"{len(done)} already done, {len(todo)} to fetch", flush=True)

BC = '"Breast Neoplasms"[Mesh]'


def term(name, typ):
    if typ == "miRNA":
        s = name.replace("hsa-", "")
        return f'("{s}"[tiab] OR "micro{s}"[tiab])'
    return f'"{name}"[tiab]'


fh = open(OUT, "a" if done else "w", newline="")
w = csv.DictWriter(fh, fieldnames=["name", "type", "pubmed_total", "pubmed_breast"])
if not done:
    w.writeheader(); fh.flush()

q = queue.Queue()
for x in todo:
    q.put(x)
cnt = [0]


def worker():
    while True:
        try:
            n, t = q.get_nowait()
        except queue.Empty:
            return
        tm = term(n, t)
        try:
            tot = count(tm); bc = count(tm + " AND " + BC)
        except Exception as e:
            print("FAIL", n, repr(e)[:80], flush=True)
            continue
        with WLOCK:
            w.writerow(dict(name=n, type=t, pubmed_total=tot, pubmed_breast=bc)); fh.flush()
            cnt[0] += 1
            if cnt[0] % 50 == 0:
                print(f"{cnt[0]}/{len(todo)} last={n}", flush=True)


ths = [threading.Thread(target=worker, daemon=True) for _ in range(6)]
for t_ in ths:
    t_.start()
for t_ in ths:
    t_.join()
fh.close()
print("done", flush=True)
