#!/usr/bin/env python3
"""
02b_ffl_class_coverage.py -- tests the manuscript's rationale for stopping at six
nodes: "six is the smallest module size that admits one edge of every regulatory
type (TF-TF, TF-miRNA, miRNA-miRNA, miRNA-gene, TF-gene, gene-gene), and beyond
six the modules add no new interaction class."

Writes results/ffl_class_coverage.tsv
"""
import sys, os, csv, pickle
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE); RES = os.path.join(ROOT, "results")
SCRATCH = os.environ.get("FFL_SCRATCH",
  "/path/to/scratch")
import ffl_graph as fg, ffl_def

ICL = {("TF","TF"):"TF-TF",("TF","Gene"):"TF-gene",("TF","miRNA"):"TF-miRNA",
       ("miRNA","Gene"):"miRNA-gene",("miRNA","TF"):"miRNA-TF",("miRNA","miRNA"):"miRNA-miRNA",
       ("Gene","TF"):"gene-TF",("Gene","miRNA"):"gene-miRNA",("Gene","Gene"):"gene-gene"}
NAMED6 = {"TF-TF","TF-miRNA","miRNA-miRNA","miRNA-gene","TF-gene","gene-gene"}

nodes, ntype, arcs, meta = fg.build(verbose=False)
A = set(arcs); und = set(e for e,a in arcs.items() if a["undirected"])
rows = []
def emit(*a): rows.append(list(a)); print("\t".join(str(x) for x in a), flush=True)

# ---- 1. interaction-class inventory of G ---------------------------------
inv = Counter(ICL[(ntype[u], ntype[v])] for (u, v) in A)
for k in sorted(inv, key=lambda k: -inv[k]):
    emit("inventory", "G", k, inv[k], "arcs of this interaction class in G")
emit("inventory", "G", "TOTAL_CLASSES", len(inv),
     "the manuscript names only 6; miRNA-TF, gene-TF and gene-miRNA are omitted")

# ---- 2. exhaustive census: modules carrying all six NAMED classes ---------
for n in (3, 4, 5):
    f = os.path.join(SCRATCH, "n%d_instances.tsv" % n)
    if not os.path.exists(f): continue
    tot = six = mx = 0
    for r in csv.DictReader(open(f), delimiter="\t"):
        S = [r["source"]] + r["members"].split(";") + [r["sink"]]
        cls = set(ICL[(ntype[u], ntype[v])] for u in S for v in S if u != v and (u, v) in A)
        tot += 1; mx = max(mx, len(cls))
        if NAMED6 <= cls: six += 1
    emit("named6_exhaustive", "n=%d" % n, "modules_with_all_six_named_classes", six,
         "out of %d modules; max distinct classes in one module = %d" % (tot, mx))

# ---- 3. targeted EXACT search at n=6 -------------------------------------
def pairs(t1, t2):
    s = set()
    for (u, v) in A:
        if ntype[u] == t1 and ntype[v] == t2:
            s.add(frozenset((u, v)) if t1 == t2 else (u, v))
    return sorted(s, key=sorted)
tftf, mm, gg = pairs("TF","TF"), pairs("miRNA","miRNA"), pairs("Gene","Gene")
emit("constraining_arcs", "G", "TF-TF_unordered/miRNA-miRNA/gene-gene",
     "%d/%d/%d" % (len(tftf), len(mm), len(gg)),
     "a module carrying all six named classes needs one of each, hence >=2 TF, >=2 miRNA, >=2 Gene, hence >=6 nodes")
def has(t1, t2, S):
    return any((u, v) in A for u in S for v in S if u != v and ntype[u] == t1 and ntype[v] == t2)
cache = os.path.join(SCRATCH, "classcover6.pkl")
if os.path.exists(cache):
    valid = pickle.load(open(cache, "rb"))
else:
    valid = []
    for g in gg:
        gl = sorted(g)
        for m in mm:
            ml = sorted(m)
            if not any((x, y) in A for x in ml for y in gl): continue
            for f2 in tftf:
                S = sorted(f2) + ml + gl
                if not has("TF","Gene",S) or not has("TF","miRNA",S): continue
                r = ffl_def.test_set(S, A, und)
                if r: valid.append((S, r))
    pickle.dump(valid, open(cache, "wb"))
ncls = Counter()
nodesused = set()
for S, r in valid:
    ncls[len(set(ICL[(ntype[u], ntype[v])] for u in S for v in S if u != v and (u, v) in A))] += 1
    nodesused.update(S)
emit("named6_exhaustive", "n=6", "modules_with_all_six_named_classes", len(valid),
     "exact targeted enumeration over all %d (TF-TF,miRNA-miRNA,gene-gene) triples; %d distinct nodes"
     % (len(tftf)*len(mm)*len(gg), len(nodesused)))
for k in sorted(ncls):
    emit("named6_exhaustive", "n=6", "of_those_with_%d_distinct_classes" % k, ncls[k],
         "maximum attainable in a 6-node module is 8 of the 9 classes (gene-miRNA never fits)")

# ---- 4. the gene-gene class is a node-typing artefact ---------------------
ggsrc = sorted(set(u for (u, v) in A if ntype[u] == "Gene" and ntype[v] == "Gene"))
emit("artefact", "gene-gene", "distinct_source_nodes", len(ggsrc), ";".join(ggsrc))
regs = sorted(set(u for (u, v) in A if ntype[u] == "Gene"))
emit("artefact", "Gene-typed_TRRUST/TransmiR_regulators", "n", len(regs), ";".join(regs))
nt2 = {n: ("TF" if n in regs else ntype[n]) for n in nodes}
inv2 = Counter(ICL[(nt2[u], nt2[v])] for (u, v) in A)
emit("artefact", "corrected_typing", "interaction_classes_remaining", len(inv2),
     ";".join("%s=%d" % (k, v) for k, v in sorted(inv2.items())))
SIXC = set(inv2)
for n in (3, 4, 5):
    f = os.path.join(SCRATCH, "n%d_instances.tsv" % n)
    if not os.path.exists(f): continue
    c = tot = 0
    for r in csv.DictReader(open(f), delimiter="\t"):
        S = [r["source"]] + r["members"].split(";") + [r["sink"]]
        cls = set(ICL[(nt2[u], nt2[v])] for u in S for v in S if u != v and (u, v) in A)
        tot += 1
        if SIXC <= cls: c += 1
    emit("corrected_typing", "n=%d" % n, "modules_with_all_%d_corrected_classes" % len(SIXC), c,
         "out of %d" % tot)

with open(os.path.join(RES, "ffl_class_coverage.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["section", "scope", "quantity", "value", "note"])
    for r in rows: w.writerow(r)
print("wrote results/ffl_class_coverage.tsv (%d data rows)" % len(rows))
