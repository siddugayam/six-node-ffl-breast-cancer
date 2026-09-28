#!/usr/bin/env python3
import csv,os,re,collections
SIF="/path/to/revision/analyses/original_submission_code/SIF_files"
ATT="/path/to/revision/analyses/original_submission_code/node_attributes"
REV="/path/to/revision"
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
CL={("miRNA","Gene"):"miRNA_target",("miRNA","TF"):"miRNA_target",("TF","miRNA"):"TF_miRNA",
    ("TF","Gene"):"TF_target",("TF","TF"):"TF_target",("miRNA","miRNA"):"miRNA_miRNA",
    ("Gene","Gene"):"gene_gene",("Gene","TF"):"gene_gene"}
# finer classes matching manuscript wording
def fine(a,b,types): return f"{types[a]}->{types[b]}"

P("== per-network fine-grained edge counts (RAW labels, unique pairs) vs manuscript ==")
for n in NETS:
    u=set(raw[n])
    c=collections.Counter(fine(a,b,rt) for a,b in u)
    P(f"  {n:7s} uniq={len(u):5d} {dict(sorted(c.items()))}")
P("\n== per-network fine-grained edge counts (HARMONISED, unique pairs) ==")
for n in NETS:
    u={(canon[a],canon[b]) for a,b in raw[n]}
    c=collections.Counter(fine(a,b,ct) for a,b in u)
    P(f"  {n:7s} uniq={len(u):5d} {dict(sorted(c.items()))}")

P("\n== self-loops in the raw SIFs ==")
for n in NETS:
    sl=[e for e in set(raw[n]) if e[0]==e[1]]
    if sl: P(f"  {n}: {sl}")

P("\n== degree of selected nodes in EVERY network (attribute CSV Degree column) ==")
for nd in ["VEGFA","CCND2","ADAMTS5","TGFBR2","TP53","MYC","RELA","SP1",
           "hsa-miR-130a","hsa-mir-130a","hsa-miR-21","hsa-mir-21","hsa-miR-124","hsa-mir-124",
           "hsa-miR-34a","hsa-mir-34a"]:
    P(f"  {nd:15s} " + "  ".join(f"{n}={deg[n].get(nd,'-')}" for n in NETS))

P("\n== harmonised (merged-copy) degree in each network ==")
for nd in ["VEGFA","CCND2","ADAMTS5","TGFBR2","TP53","MYC","RELA","SP1",
           "hsa-miR-130a","hsa-miR-21","hsa-miR-124","hsa-miR-34a"]:
    row=[]
    for n in NETS:
        u={(canon[a],canon[b]) for a,b in raw[n]}
        row.append(f"{n}={sum(1 for a,b in u if a==nd or b==nd)}")
    P(f"  {nd:15s} " + "  ".join(row))

# ---- TRRUST with the alternative 'resolve Unknown' rule
P("\n== TRRUST, alternative resolution rules ==")
tr=collections.defaultdict(set)
for line in open(f"{REV}/data/db/trrust_human.tsv"):
    f=line.rstrip("\n").split("\t")
    if len(f)<3 or f[0]=="TF": continue
    tr[(f[0].strip(),f[1].strip())].add(f[2].strip())
E=set()
for n in NETS:
    for a,b in raw[n]: E.add((canon[a],canon[b]))
tft=sorted(e for e in E if CL[(ct[e[0]],ct[e[1]])]=="TF_target")
tft_ns=[e for e in tft if e[0]!=e[1]]
def rule_strict(m):
    if len(m)>1: return "Ambiguous"
    return next(iter(m))
def rule_resolve(m):
    if "Activation" in m and "Repression" in m: return "Ambiguous"
    if "Activation" in m: return "Activation"
    if "Repression" in m: return "Repression"
    return "Unknown"
for name,f in [("strict (any multi-mode = Ambiguous)",rule_strict),
               ("resolve (Unknown absorbed; Act+Rep = Ambiguous)",rule_resolve)]:
    for lab,ES in [("771 incl self-loop",tft),("770 excl self-loop",tft_ns)]:
        c=collections.Counter(f(tr[e]) if e in tr else "not_in_TRRUST" for e in ES)
        P(f"  {name:48s} [{lab}] {dict(c)}  Rep%={100*c['Repression']/len(ES):.1f}")
P(f"  TF_target edges NOT in TRRUST: {[e for e in tft if e not in tr]}")
# and TF_miRNA / other classes coverage for context
open(f"{REV}/logs/v2/v2_extra.log","w").write("\n".join(L)+"\n")
