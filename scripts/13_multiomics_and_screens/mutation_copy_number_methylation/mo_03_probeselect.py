# Selects the 450k methylation probes for the promoter analysis of the network genes and of a background gene
# set, and writes the miRNA loci (hg19) to cache/multiomics/.
import gzip, random, sys
RV="/path/to/revision"; CA=RV+"/cache/multiomics"
# network nodes
genes=set(); mirs=[]
with open(RV+"/data/canonical_nodes.tsv") as f:
    next(f)
    for l in f:
        p=l.rstrip("\n").split("\t")
        if p[1] in ("Gene","TF"): genes.add(p[0])
        else: mirs.append(p[0])
print("network protein-coding:",len(genes)," miRNA:",len(mirs))

# miRNA precursor coords hg19 from wgRna
import re
pre=[]
with gzip.open(CA+"/wgRna_hg19.txt.gz","rt") as f:
    for l in f:
        p=l.rstrip("\n").split("\t")
        if p[9]=="miRNA": pre.append((p[4],p[1],int(p[2]),int(p[3]),p[6]))
def mat2pre(m):
    stem=re.sub(r"-(3p|5p)$","",m.replace("hsa-","").lower())
    rx=re.compile(r"^hsa-"+re.escape(stem)+r"(-[0-9]+)?$")
    return [x for x in pre if rx.match(x[0])]
mirloci={}
for m in mirs:
    hits=mat2pre(m)
    if hits: mirloci[m]=hits
print("miRNAs with hg19 precursor:",len(mirloci),"/",len(mirs))
def mk(n):
    y=n.replace("hsa-","",1)
    y="MIR"+y[4:] if y.startswith("mir-") else ("MIRLET"+y[4:] if y.startswith("let-") else y)
    return y.replace("-","").upper()
mirsyms=set(mk(h[0]) for v in mirloci.values() for h in v)
print("miRNA precursor symbols:",len(mirsyms))

# interval index for miRNA windows (+/- 10kb around precursor)
WIN=10000
byc={}
for m,hits in mirloci.items():
    for nm,ch,st,en,strand in hits:
        byc.setdefault(ch,[]).append((st-WIN,en+WIN,m,nm))

# background genes: 2000 random symbols from probeMap not in network
allsym=set()
rows=[]
with open(CA+"/probeMap_illuminaMethyl450_hg19") as f:
    hdr=next(f)
    for l in f:
        p=l.rstrip("\n").split("\t")
        rows.append(p)
        for g in p[1].split(","):
            if g and g!=".": allsym.add(g)
print("probeMap probes:",len(rows)," symbols:",len(allsym))
random.seed(42)
bg=set(random.sample(sorted(allsym-genes-mirsyms),2000))

sel=set(); reason={}
for p in rows:
    pid,gf,ch,st,en,strand=p[0],p[1],p[2],p[3],p[4],p[5]
    gs=set(g for g in gf.split(",") if g and g!=".")
    hit=gs & genes
    if hit: sel.add(pid); reason.setdefault(pid,set()).add("net_gene")
    if gs & mirsyms: sel.add(pid); reason.setdefault(pid,set()).add("mir_symbol")
    if gs & bg: sel.add(pid); reason.setdefault(pid,set()).add("background")
    try: pos=int(st)
    except: continue
    for a,b,m,nm in byc.get(ch,[]):
        if a<=pos<=b: sel.add(pid); reason.setdefault(pid,set()).add("mir_window"); break
print("probes selected:",len(sel))
with open(CA+"/probes_selected.txt","w") as o:
    for p in sorted(sel): o.write(p+"\n")
with open(CA+"/background_genes.txt","w") as o:
    for g in sorted(bg): o.write(g+"\n")
# also dump miRNA loci table for R
with open(CA+"/mirna_loci_hg19.tsv","w") as o:
    o.write("mature\tprecursor\tchrom\tstart\tend\tstrand\tsym\n")
    for m,hits in mirloci.items():
        for nm,ch,st,en,strand in hits:
            o.write(f"{m}\t{nm}\t{ch}\t{st}\t{en}\t{strand}\t{mk(nm)}\n")
print("done")
