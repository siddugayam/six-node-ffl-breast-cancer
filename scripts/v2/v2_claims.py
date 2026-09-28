#!/usr/bin/env python3
"""(c) case duplicates, (e) TRRUST, (f) manuscript claims, (g) diff vs existing canonical."""
import csv,os,re,collections,json,sys
SIF="/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/SIF_files"
ATT="/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/miRNA_Github_GPR/node_attributes"
REV="/path/to/revision"
OUT=f"{REV}/results/v2"
NETS=["3-miR","3-TF","3-Comp","4-TF","5-TF","6-TF"]
L=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); L.append(s)

raw={n:[(l.rstrip('\n').split("\t")[0].strip(),l.rstrip('\n').split("\t")[2].strip())
        for l in open(f"{SIF}/{n}.sif") if l.strip()] for n in NETS}
at=collections.defaultdict(set); deg={}
for n in NETS:
    d={}
    for r in csv.DictReader(open(f"{ATT}/{n}.csv",newline="")):
        at[r["name"].strip()].add(r["Type"].strip()); d[r["name"].strip()]=int(r["Degree"])
    deg[n]=d
RAW=set(at)|{x for n in NETS for e in raw[n] for x in e}
M=re.compile(r"^hsa-(mir|miR|let)-"); TFS={k for k,v in at.items() if "TF" in v}
rt={l:("TF" if l in TFS else "miRNA" if M.match(l) else "Gene") for l in RAW}
def cn(l): return "hsa-miR-"+l[8:] if l.startswith("hsa-mir-") else l
cns={cn(l) for l in RAW}; PR=re.compile(r"^(hsa-(?:miR|let)-.+?)-[12]$")
def co(l):
    m=PR.match(l); return m.group(1) if m and m.group(1) in cns else l
canon={l:co(cn(l)) for l in RAW}
ct={}
for l in RAW: ct.setdefault(canon[l],set()).add(rt[l])
ct={c:("TF" if "TF" in s else "miRNA" if "miRNA" in s else "Gene") for c,s in ct.items()}
E=collections.defaultdict(set)
for n in NETS:
    for s,t in raw[n]: E[(canon[s],canon[t])].add(n)
edges=set(E)
CLASS={("miRNA","Gene"):"miRNA_target",("miRNA","TF"):"miRNA_target",("TF","miRNA"):"TF_miRNA",
       ("TF","Gene"):"TF_target",("TF","TF"):"TF_target",("miRNA","miRNA"):"miRNA_miRNA",
       ("Gene","Gene"):"gene_gene",("Gene","TF"):"gene_gene",("Gene","miRNA"):"gene_miRNA"}
ecls={e:CLASS[(ct[e[0]],ct[e[1]])] for e in edges}

# ---------------------------------------------------------------- (c) detail
P("="*70); P("(c) miRBase case duplication detail")
mir_raw={l for l in RAW if rt[l]=="miRNA"}
groups=collections.defaultdict(set)
for l in mir_raw: groups[l.lower()].add(l)
both={k:v for k,v in groups.items() if len(v)>1}
P(f"  raw miRNA labels                                  = {len(mir_raw)}")
P(f"  distinct case-insensitive miRNA names             = {len(groups)}")
P(f"  names present in BOTH miRBase cases (raw labels)  = {len(both)}")
# what if the -1/-2 precursor forms are first collapsed?
gc=collections.defaultdict(set)
for l in mir_raw: gc[canon[l].lower()].add(l)
both_c={k:v for k,v in gc.items() if len({x.lower() for x in v})>1}
P(f"  ...after -1/-2 precursor collapse                 = {len(both_c)}")
P(f"  the extra name vs 214: {sorted(set(both)-{k.lower() for k in both_c})}")
P(f"  '129' group members: {sorted(groups.get('hsa-mir-129',set()))} ; '129-2': {sorted(groups.get('hsa-mir-129-2',set()))}")

outt=collections.defaultdict(set)
for n in NETS:
    for s,t in raw[n]: outt[s].add(t)
P("\n  Jaccard of the two copies' RAW target sets:")
jac={}
for b in ["130a","21","124","29a","34a"]:
    a1,a2=f"hsa-mir-{b}",f"hsa-miR-{b}"
    A,B=outt.get(a1,set()),outt.get(a2,set())
    j=len(A&B)/len(A|B) if A|B else float('nan')
    jac[b]=j
    P(f"   miR-{b:5s} lower n={len(A):3d} upper n={len(B):3d} shared={len(A&B):2d} {sorted(A&B)} union={len(A|B):3d} J={j:.4f}")

# ------------------------------------------------------------- (e) TRRUST
P("\n"+"="*70); P("(e) TRRUST coverage of published TF_target edges")
tr=collections.defaultdict(set)
nline=0
with open(f"{REV}/data/db/trrust_human.tsv") as fh:
    for line in fh:
        line=line.rstrip("\n")
        if not line.strip(): continue
        f=line.split("\t")
        if f[0]=="TF" or len(f)<3: continue
        nline+=1
        tr[(f[0].strip(),f[1].strip())].add(f[2].strip())
P(f"  TRRUST rows parsed = {nline}; distinct TF-target pairs = {len(tr)}")
P(f"  mode vocabulary = {sorted({m for v in tr.values() for m in v})}")
multi=sum(1 for v in tr.values() if len(v)>1)
P(f"  TRRUST pairs with >1 distinct mode = {multi}")

tf_edges=sorted([e for e in edges if ecls[e]=="TF_target"])
tf_edges_nsl=[e for e in tf_edges if e[0]!=e[1]]
P(f"  TF_target edges in my network = {len(tf_edges)} (excluding the 1 self-loop: {len(tf_edges_nsl)})")
for label,ES in [("incl. self-loop",tf_edges),("excl. self-loop",tf_edges_nsl)]:
    cnt=collections.Counter()
    for e in ES:
        modes=tr.get(e)
        if modes is None: cnt["not_in_TRRUST"]+=1
        elif len(modes)>1: cnt["Ambiguous"]+=1
        else: cnt[next(iter(modes))]+=1
    cov=len(ES)-cnt["not_in_TRRUST"]
    P(f"  [{label}] n={len(ES)}  covered={cov} ({100*cov/len(ES):.1f}%)  "
      f"Activation={cnt['Activation']} Repression={cnt['Repression']} "
      f"({100*cnt['Repression']/len(ES):.1f}% of all) Unknown={cnt['Unknown']} "
      f"Ambiguous={cnt['Ambiguous']} not_in_TRRUST={cnt['not_in_TRRUST']}")

# ------------------------------------------------------ (f) manuscript claims
P("\n"+"="*70); P("(f) manuscript-specific claims")
def targets_of(node, netlist=NETS, harmon=True):
    s=set()
    for n in netlist:
        for a,b in raw[n]:
            if (canon[a] if harmon else a)==node: s.add(canon[b] if harmon else b)
    return s
t130=targets_of("hsa-miR-130a")
P(f"  miR-130a distinct targets (harmonised, all 6 networks) = {len(t130)}  [claim: 109]")
P(f"    of which Gene/TF = {sum(1 for x in t130 if ct[x]!='miRNA')}, miRNA = {sum(1 for x in t130 if ct[x]=='miRNA')}")
for g in ["VEGFA","MMP2","COL1A1"]:
    P(f"    {g} a target of miR-130a? {g in t130}   ({g} in network at all: {g in ct})")
P(f"    raw-label copies: hsa-miR-130a n={len(outt.get('hsa-miR-130a',set()))}, hsa-mir-130a n={len(outt.get('hsa-mir-130a',set()))}")
P(f"    3-Comp only, harmonised: {len(targets_of('hsa-miR-130a',['3-Comp']))}")

P(f"\n  HK2 present as any node?  {'HK2' in RAW}   in any edge? "
  f"{any('HK2' in e for e in edges)}")
P(f"  let-7b nodes: {sorted(l for l in RAW if 'let-7b' in l)}")

s35=["miR-124","miR-133a","miR-133b","miR-143","miR-218","miR-3619","miR-412","miR-4776",
     "miR-4800","miR-516a","miR-516b","miR-7162","miR-92a-1"]
P(f"\n  Section 3.5 miRNAs claimed to target COL1A1 (n={len(s35)}):")
present=[];absent=[];targets_col1a1=[]
allnodes=set(ct)
for m in s35:
    cands=[l for l in RAW if l.lower() in (f"hsa-{m}".lower(), f"hsa-{m}".lower().replace("mir-","mir-"))]
    # robust: case-insensitive match on hsa-<m>
    cands=[l for l in RAW if l.lower()==f"hsa-{m}".lower()]
    cnodes={canon[l] for l in cands}
    # also allow the -1 precursor collapse (e.g. miR-92a-1 -> hsa-miR-92a)
    if not cands:
        alt=[l for l in RAW if l.lower()==f"hsa-{m}".lower().rsplit('-',1)[0]] if re.search(r"-[12]$",m) else []
        cands=alt; cnodes={canon[l] for l in cands}
    hit=bool(cnodes)
    tcol = any(("COL1A1" in targets_of(c)) for c in cnodes) if hit else False
    (present if hit else absent).append(m)
    if tcol: targets_col1a1.append(m)
    P(f"    {m:10s} in network: {str(hit):5s} node(s)={sorted(cnodes)}  ->COL1A1: {tcol}")
P(f"  present={len(present)} absent={len(absent)}  absent list={absent}")
P(f"  of those present, actually have an edge to COL1A1: {len(targets_col1a1)} {targets_col1a1}")

P("\n  Degrees in 3-Comp.sif (raw label vs de-duplicated harmonised node):")
P(f"  {'node':16s} {'rawlines':>9s} {'raw_uniq':>9s} {'harm_uniq':>10s} {'attrCSV':>8s}")
def degs(node):
    rawlines=sum(1 for a,b in raw["3-Comp"] if a==node or b==node)
    rawu=len({(a,b) for a,b in raw["3-Comp"] if a==node or b==node})
    return rawlines,rawu
for nd in ["VEGFA","CCND2","TP53","MYC","RELA","SP1","hsa-miR-130a","hsa-miR-21","hsa-miR-124","hsa-miR-34a"]:
    rl,ru=degs(nd)
    hn=canon.get(nd,nd)
    hu=len({(canon[a],canon[b]) for a,b in raw["3-Comp"] if canon[a]==hn or canon[b]==hn})
    csvdeg=deg["3-Comp"].get(nd,"NA")
    P(f"  {nd:16s} {rl:9d} {ru:9d} {hu:10d} {str(csvdeg):>8s}")
P("  (for miRNAs the 'harm_uniq' column merges the hsa-mir-/hsa-miR- copies)")
for nd in ["hsa-mir-130a","hsa-mir-21","hsa-mir-124","hsa-mir-34a"]:
    rl,ru=degs(nd)
    P(f"  lower-case copy {nd:16s} rawlines={rl} raw_uniq={ru} attrCSV={deg['3-Comp'].get(nd,'NA')}")
# top degree genes in 3-Comp for context
hd=collections.Counter()
for a,b in {(canon[a],canon[b]) for a,b in raw["3-Comp"]}:
    hd[a]+=1; hd[b]+=1
P(f"  top-6 harmonised-degree nodes in 3-Comp: {hd.most_common(6)}")
rd=collections.Counter()
for a,b in set(raw["3-Comp"]):
    rd[a]+=1; rd[b]+=1
P(f"  top-6 RAW-label-degree nodes in 3-Comp:  {rd.most_common(6)}")

# ------------------------------------------------------------------ (g) diff
P("\n"+"="*70); P("(g) DIFF vs existing revision/data/canonical_*.tsv")
ex_e=set(); ex_cls={}
for r in csv.DictReader(open(f"{REV}/data/canonical_edges.tsv"),delimiter="\t"):
    ex_e.add((r["source"],r["target"])); ex_cls[(r["source"],r["target"])]=r["edge_type"]
ex_n={}
for r in csv.DictReader(open(f"{REV}/data/canonical_nodes.tsv"),delimiter="\t"):
    ex_n[r["name"]]=r["type"]
P(f"  existing: {len(ex_e)} edges, {len(ex_n)} nodes | mine: {len(edges)} edges, {len(ct)} nodes")
P(f"  edges in mine only  ({len(edges-ex_e)}): {sorted(edges-ex_e)}")
P(f"  edges in theirs only({len(ex_e-edges)}): {sorted(ex_e-edges)}")
diffcls=[(e,ecls[e],ex_cls[e]) for e in edges&ex_e if ecls[e]!=ex_cls[e]]
P(f"  edges present in both but with a DIFFERENT edge_type: {len(diffcls)} {diffcls[:10]}")
P(f"  nodes in mine only  ({len(set(ct)-set(ex_n))}): {sorted(set(ct)-set(ex_n))}")
P(f"  nodes in theirs only({len(set(ex_n)-set(ct))}): {sorted(set(ex_n)-set(ct))}")
difft=[(n,ct[n],ex_n[n]) for n in set(ct)&set(ex_n) if ct[n]!=ex_n[n]]
P(f"  nodes with a DIFFERENT type: {len(difft)} {difft[:10]}")
P(f"  existing per-class counts: {dict(collections.Counter(ex_cls.values()))}")
P(f"  my       per-class counts: {dict(collections.Counter(ecls.values()))}")

open(f"{REV}/logs/v2/v2_claims.log","w").write("\n".join(L)+"\n")
