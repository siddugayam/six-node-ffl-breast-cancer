#!/usr/bin/env python3
# Collects the counts of the C census runs from their logs (logs/v2) into results/v2/verify_census.csv: count
# per graph and module size, counting method, seeds, maximum number of edge classes and the ratio to the
# counts of the original submission.
import os, csv, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from census_11_class_diversity import parse, TAX, diversity
from collections import Counter

BASE = "/path/to/revision"
L = os.path.join(BASE, "logs", "v2")
PUB = {3: 6037, 4: 1.2e5, 5: 1.5e6, 6: 1.9e7, 7: 2.7e8}

JOBS = [
 ("pub_augmented", 3, ["census2_n3.log"], "exhaustive"),
 ("pub_augmented", 4, ["census2_n4.log"], "exhaustive"),
 ("pub_augmented", 5, ["census2_n5.log"], "exhaustive"),
 ("pub_augmented", 6, ["census2_n6_exhaustive.log"], "exhaustive"),
 ("pub_augmented", 5, ["census_n5_sampled.log"], "RAND-ESU q=(1,1,.3,.3)"),
 ("pub_augmented", 6, ["census_n6_sampled.log"], "RAND-ESU q=(1,1,.3,.3,.3)"),
 ("pub_augmented", 7, ["census_n7_sampled.log", "census_n7_sampled_altq.log"], "RAND-ESU (2 q-profiles)"),
 ("deposited_only", 3, ["census2_dep_n3.log"], "exhaustive"),
 ("deposited_only", 4, ["census2_dep_n4.log"], "exhaustive"),
 ("deposited_only", 5, ["census2_dep_n5.log"], "exhaustive"),
 ("deposited_only", 6, ["census_dep_n67.log"], "exhaustive"),
 ("deposited_only", 7, ["census_dep_n67.log"], "RAND-ESU q=(1,1,.3,.3,.3,.3)"),
]
rows = []
for g, n, files, meth in JOBS:
    runs = []
    for f in files:
        p = os.path.join(L, f)
        if not os.path.exists(p): continue
        runs += [r for r in parse(p) if r["K"] == n and (("exhaustive" in meth) == (set(r["q"].split(",")) <= {"1", ""}))]
    if not runs: continue
    est = [r["est"] for r in runs]
    m = statistics.mean(est); sd = statistics.stdev(est) if len(est) > 1 else 0.0
    mask = Counter()
    for r in runs:
        for k, v in r["mask"].items(): mask[k] += v
    div = {}
    for tname, tax in TAX.items():
        b, c, h = diversity(mask, tax); div[tname] = b
    rows.append(dict(graph=g, n=n, method=meth, n_independent_seeds=len(runs),
                     mean_count=round(m, 1), sd_count=round(sd, 1),
                     min_count=round(min(est), 1), max_count=round(max(est), 1),
                     published_claim=(PUB[n] if g == "pub_augmented" else ""),
                     ratio_to_published=(round(m / PUB[n], 4) if g == "pub_augmented" else ""),
                     max_classes_fine8=div["fine8 (source_type -> target_type)"],
                     max_classes_7cls=div["7-class (all Gene-> edges pooled as gene_gene)"],
                     max_classes_6cls=div["6-class coarse edge_type (miRNA_target pools miRNA->Gene/TF)"],
                     total_visited=sum(r["visited"] for r in runs)))
out = os.path.join(BASE, "results", "v2", "verify_census.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow(r)
print("wrote", out)
for r in rows:
    print("%-15s n=%d %-26s seeds=%d mean=%16.1f sd=%12.1f  pub=%-8s ratio=%-7s  maxcls fine8=%d 7cls=%d 6cls=%d"
          % (r["graph"], r["n"], r["method"], r["n_independent_seeds"], r["mean_count"], r["sd_count"],
             r["published_claim"], r["ratio_to_published"], r["max_classes_fine8"],
             r["max_classes_7cls"], r["max_classes_6cls"]))
