#!/usr/bin/env python3
"""One table across the four census graphs: exact counts n = 3-6, RAND-ESU n = 7 (8 seeds),
maximum edge classes, class sets at the maximum, share with 2 TF + 2 miRNA + 2 Gene, and the
count of modules containing a miRNA-miRNA arc.  Writes census/census_all_graphs.csv and .json."""
import os, json, csv, glob, statistics as st, collections
D = os.path.dirname(os.path.abspath(__file__))
VAR = {"table3": "legacy kept, STRING kept (Table 3 graph, 8,031 arcs)",
       "nolegacy": "legacy dropped, STRING kept (8,003 arcs)",
       "table3_nostring": "legacy kept, STRING dropped (7,198 arcs)",
       "nolegacy_nostring": "legacy dropped, STRING dropped (7,170 arcs)"}
rows = []; full = {}
for v, desc in VAR.items():
    full[v] = {}
    for k in (3, 4, 5, 6):
        d = json.load(open(f"{D}/{v}/census_n{k}.json"))
        mm = 0
        files = [f"{D}/{v}/n{k}.txt"] if k < 6 else glob.glob(f"{D}/{v}/n6_parts/part_*.txt")
        for fn in files:
            for line in open(fn):
                w = line.split()
                if w and w[0] == "mask" and (int(w[1]) >> 6) & 1: mm += int(w[2])
        rows.append(dict(graph=v, graph_desc=desc, n=k, method="exhaustive", modules=d["found"], modules_sd="",
                         max_edge_classes=d["max_classes"], n_modules_at_max=sum(d["class_sets_at_max"].values()),
                         class_sets_at_max="; ".join(f"{a}: {b}" for a, b in d["class_sets_at_max"].items()),
                         pct_2TF_2miR_2Gene=round(d["flags_pct"]["a_2TF_2miR_2Gene"], 4),
                         modules_with_miRNA_miRNA_arc=mm, visited=d["visited"]))
        full[v][k] = d
    seeds = []
    for s in (11, 22, 33, 44, 55, 71, 72, 73):
        d = json.load(open(f"{D}/{v}/n7_seeds/census_n7_seed_{s}.json")); seeds.append(d)
    est = [d["estimated_modules"] for d in seeds]; mx = max(d["max_classes"] for d in seeds)
    at7 = [sum(v2 for k2, v2 in d["n_classes_hist"].items() if int(k2) == mx) for d in seeds]
    pct = [d["flags_pct"]["a_2TF_2miR_2Gene"] for d in seeds]
    rows.append(dict(graph=v, graph_desc=desc, n=7, method="RAND-ESU, 8 seeds (q 1,1,.2,.2,.2,.2 x5; 1,1,1,.1,.15,.15 x3)",
                     modules=round(st.mean(est), 1), modules_sd=round(st.stdev(est), 1), max_edge_classes=mx,
                     n_modules_at_max=f"{min(at7)}-{max(at7)} sampled per seed", class_sets_at_max="",
                     pct_2TF_2miR_2Gene=round(st.mean(pct), 3), modules_with_miRNA_miRNA_arc="", visited=""))
    full[v][7] = dict(estimates=est, found=[d["found"] for d in seeds], visited=[d["visited"] for d in seeds])
with open(f"{D}/census_all_graphs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
json.dump(full, open(f"{D}/census_all_graphs.json", "w"), indent=1, default=str)
for r in rows:
    print(f"{r['graph']:18s} n={r['n']}  modules {r['modules']:>16,}  sd {r['modules_sd']!s:>10}  max classes {r['max_edge_classes']}  "
          f"at max {r['n_modules_at_max']!s:>10}  2/2/2 {r['pct_2TF_2miR_2Gene']}%  miR-miR {r['modules_with_miRNA_miRNA_arc']}")
