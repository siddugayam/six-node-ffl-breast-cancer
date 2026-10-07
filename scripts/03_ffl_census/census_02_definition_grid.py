#!/usr/bin/env python3
"""Brute-force grid over plausible 3-node FFL definitions, compared with the 6,037
three-node FFLs of the original submission."""
import csv,os,re,collections,itertools
SIF="/path/to/revision/analyses/original_submission_code/SIF_files"
ATT="/path/to/revision/analyses/original_submission_code/node_attributes"
NETS=["3-miR","3-TF","3-Comp","4-TF","5-TF","6-TF"]
raw={n:[(l.split("\t")[0].strip(),l.split("\t")[2].strip()) for l in open(f"{SIF}/{n}.sif") if l.strip()] for n in NETS}
at=collections.defaultdict(set)
for n in NETS:
    for r in csv.DictReader(open(f"{ATT}/{n}.csv",newline="")): at[r["name"].strip()].add(r["Type"].strip())
RAW=set(at)|{x for n in NETS for e in raw[n] for x in e}
M=re.compile(r"^hsa-(mir|miR|let)-"); TFS={k for k,v in at.items() if "TF" in v}
rt={l:("TF" if l in TFS else "miRNA" if M.match(l) else "Gene") for l in RAW}
def cn(l): return "hsa-miR-"+l[8:] if l.startswith("hsa-mir-") else l
cns={cn(l) for l in RAW}; P=re.compile(r"^(hsa-(?:miR|let)-.+?)-[12]$")
def co(l):
    m=P.match(l); return m.group(1) if m and m.group(1) in cns else l
canon={l:co(cn(l)) for l in RAW}
ct={}
for l in RAW: ct.setdefault(canon[l],set()).add(rt[l])
ct={c:("TF" if "TF" in s else "miRNA" if "miRNA" in s else "Gene") for c,s in ct.items()}

def graphs():
    H={n:{(canon[s],canon[t]) for s,t in raw[n]} for n in NETS}
    R={n:set(raw[n]) for n in NETS}
    H["UNION"]=set().union(*H.values()); R["UNION"]=set().union(*R.values())
    return {"harmonised":(H,ct),"raw":(R,rt)}

def count(E,types,srcpair,tgt,mode):
    """srcpair: allowed (type(a),type(b)) for a->b regulator chain; tgt: allowed type(c);
       mode: 'ordered' counts every (a,b,c); 'mutual_once' counts mutual pairs once."""
    E={(s,t) for s,t in E if s!=t}
    out=collections.defaultdict(set)
    for s,t in E: out[s].add(t)
    n=0
    for a in out:
        for b in out[a]:
            if (types.get(a),types.get(b)) not in srcpair: continue
            if mode=="mutual_once" and (b,a) in E and a>b: continue
            n+=sum(1 for c in (out[a]&out.get(b,set())) if c not in (a,b) and types.get(c) in tgt)
    return n

ALL={"TF","miRNA","Gene"}
srcopts={
 "TFmiR_and_miRTF":{("TF","miRNA"),("miRNA","TF")},
 "any":{(x,y) for x in ALL for y in ALL},
 "TFmiR_only":{("TF","miRNA")},
 "no_gene_reg":{(x,y) for x in ("TF","miRNA") for y in ("TF","miRNA")},
}
tgtopts={"GeneTF":{"Gene","TF"},"any":ALL,"Gene":{"Gene"}}
G=graphs()
hits=[]
print(f"{'labels':11s} {'graph':8s} {'srcpair':18s} {'tgt':7s} {'mode':12s} {'count':>8s}")
for lab,(D,types) in G.items():
    for gname in NETS+["UNION"]:
        for sn,sv in srcopts.items():
            for tn,tv in tgtopts.items():
                for mode in ("ordered","mutual_once"):
                    c=count(D[gname],types,sv,tv,mode)
                    if gname=="UNION":
                        print(f"{lab:11s} {gname:8s} {sn:18s} {tn:7s} {mode:12s} {c:8d}")
                    if c in (6037,6036,6038): hits.append((lab,gname,sn,tn,mode,c))
print("\nMATCHES to ~6037:",hits if hits else "NONE")
# sum-over-networks variants
print("\nSum over the 6 deposited networks (not the union):")
for lab,(D,types) in G.items():
    for sn,sv in srcopts.items():
        for tn,tv in tgtopts.items():
            for mode in ("ordered","mutual_once"):
                s=sum(count(D[n],types,sv,tv,mode) for n in NETS)
                flag=" <== 6037" if s in (6036,6037,6038) else ""
                print(f"  {lab:11s} {sn:18s} {tn:7s} {mode:12s} {s:8d}{flag}")
