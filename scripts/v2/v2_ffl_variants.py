#!/usr/bin/env python3
"""FFL-core counting under several explicit definitional variants, so the
published numbers can be matched to a definition rather than assumed."""
import csv, os, re, collections, itertools, json
import sys

SIF_DIR  = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files"
ATTR_DIR = "/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/node_attributes"
NETS = ["3-miR","3-TF","3-Comp","4-TF","5-TF","6-TF"]

raw_edges={}
for net in NETS:
    rows=[]
    for line in open(os.path.join(SIF_DIR,net+".sif")):
        line=line.rstrip("\n")
        if not line.strip(): continue
        f=line.split("\t"); rows.append((f[0].strip(),f[2].strip()))
    raw_edges[net]=rows
attr_type=collections.defaultdict(set)
for net in NETS:
    for r in csv.DictReader(open(os.path.join(ATTR_DIR,net+".csv"),newline="")):
        attr_type[r["name"].strip()].add(r["Type"].strip())
RAW=set(attr_type)|{x for net in NETS for e in raw_edges[net] for x in e}
MIR=re.compile(r"^hsa-(mir|miR|let)-")
tfs={n for n,t in attr_type.items() if "TF" in t}
rtype={l:("TF" if l in tfs else ("miRNA" if MIR.match(l) else "Gene")) for l in RAW}
def cn(l): return "hsa-miR-"+l[8:] if l.startswith("hsa-mir-") else l
cns={cn(l) for l in RAW}
P1=re.compile(r"^(hsa-(?:miR|let)-.+?)-([12])$")
def col(l):
    m=P1.match(l)
    return m.group(1) if (m and m.group(1) in cns) else l
canon={l:col(cn(l)) for l in RAW}
ctype={}
for l in RAW:
    c=canon[l]; ctype.setdefault(c,set()).add(rtype[l])
ctype={c:("TF" if "TF" in s else ("miRNA" if "miRNA" in s else "Gene")) for c,s in ctype.items()}

def ffl(pairs, types, target_types, exclusive):
    E={(s,t) for s,t in pairs if s!=t}
    out=collections.defaultdict(set)
    for s,t in E: out[s].add(t)
    nodes={x for p in E for x in p}
    tf_n=[n for n in nodes if types.get(n)=="TF"]
    mi_n=[n for n in nodes if types.get(n)=="miRNA"]
    tf=mi=cp=0
    for a in tf_n:
        for m in out.get(a,()):
            if types.get(m)!="miRNA": continue
            mutual=(m,a) in E
            sh={g for g in (out.get(a,set())&out.get(m,set())) if g not in (a,m) and types.get(g) in target_types}
            k=len(sh)
            if mutual:
                cp+=k
                if not exclusive: tf+=k
            else: tf+=k
    for m in mi_n:
        for a in out.get(m,()):
            if types.get(a)!="TF": continue
            mutual=(a,m) in E
            sh={g for g in (out.get(a,set())&out.get(m,set())) if g not in (a,m) and types.get(g) in target_types}
            k=len(sh)
            if mutual:
                if not exclusive: mi+=k
            else: mi+=k
    return tf,mi,cp

def ffl_agnostic(pairs):
    """Type-agnostic exhaustive FFL: ordered triples a->b, b->c, a->c (a!=b!=c)."""
    E={(s,t) for s,t in pairs if s!=t}
    out=collections.defaultdict(set)
    for s,t in E: out[s].add(t)
    n=0
    for a in out:
        oa=out[a]
        for b in oa:
            n+=len({c for c in out.get(b,()) if c in oa and c!=a})
    return n

print("VARIANT TABLE  (harmonised labels)")
for tt,ttl in [(("Gene","TF"),"g in {Gene,TF}"), (("Gene","TF","miRNA"),"g any type"), (("Gene",),"g Gene only")]:
    for excl in (False,True):
        print(f"\n--- target set: {ttl} | composites {'EXCLUDED from TF/miR classes' if excl else 'ALSO counted in TF & miR classes'}")
        print(f"{'network':9s} {'TF-FFL':>8s} {'miRNA-FFL':>10s} {'Composite':>10s} {'TF+miR':>8s} {'TF+miR+Comp':>12s}")
        for net in NETS+["UNION"]:
            if net=="UNION":
                pairs={(canon[s],canon[t]) for nn in NETS for s,t in raw_edges[nn]}
            else:
                pairs={(canon[s],canon[t]) for s,t in raw_edges[net]}
            a,b,c=ffl(pairs,ctype,set(tt),excl)
            print(f"{net:9s} {a:8d} {b:10d} {c:10d} {a+b:8d} {a+b+c:12d}")

print("\n\nTYPE-AGNOSTIC EXHAUSTIVE 3-node FFL (any a->b->c with a->c), harmonised:")
for net in NETS+["UNION"]:
    if net=="UNION":
        pairs={(canon[s],canon[t]) for nn in NETS for s,t in raw_edges[nn]}
    else:
        pairs={(canon[s],canon[t]) for s,t in raw_edges[net]}
    print(f"   {net:9s} {ffl_agnostic(pairs):8d}")
