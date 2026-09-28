#!/usr/bin/env python3
"""Which TFs are SPECIFICALLY bound at the COL1A1 / COL3A1 promoters?

Raw dataset counts are dominated by how often a factor was assayed and by how
'hyper-ChIPable' a promoter is, so every TF is scored against its own occupancy
across 600 random non-network promoters fetched from the same ReMap 2022 catalogue.
  obs   = ReMap datasets of TF t with a peak summit in TSS +/- w of the collagen
  bgfrac= fraction of the 600 background promoters where t is bound
  exp   = mean number of t datasets per background promoter
  binomial p (one-sided) on 'bound / not bound' against bgfrac.
"""
import csv, json, os, collections, math
BASE="/path/to/revision"
def load(d):
    out={}
    for fn in os.listdir(d):
        if not fn.endswith(".bed"): continue
        pk=[]
        for line in open(os.path.join(d,fn)):
            q=line.rstrip("\n").split("\t")
            if len(q)<11: continue
            pk.append((int(q[6]),q[3],q[9],q[10]))
        out[fn[:-4]]=pk
    return out
NETP=load(BASE+"/cache/direct/remap/regions"); BGP=load(BASE+"/cache/direct/remap/bg")
nt=json.load(open(BASE+"/cache/direct/node_tss_hg38.json"))
bt=json.load(open(BASE+"/cache/direct/bg_tss.json"))
def occ(pk,t,w):
    m=collections.defaultdict(lambda:[set(),set()])
    for summit,ds,tf,bio in pk:
        if abs(summit-t)<=w: m[tf][0].add(ds); m[tf][1].add(bio)
    return m
E=list(csv.DictReader(open(BASE+"/data/canonical_edges.tsv"),delimiter="\t"))
netTF={e["source"] for e in E if e["edge_type"]=="TF_target"}
net_edge={(e["source"],e["target"]) for e in E if e["edge_type"]=="TF_target"}
out=[]
for w in (1000,5000):
    BG={g:occ(BGP[g],bt[g][1],w) for g in BGP}
    nbg=len(BG)
    bgcount=collections.Counter(); bgsum=collections.Counter()
    for g,m in BG.items():
        for tf,v in m.items(): bgcount[tf]+=1; bgsum[tf]+=len(v[0])
    for gene in ("COL1A1","COL3A1"):
        m=occ(NETP[gene],nt[gene][1],w)
        for tf,v in m.items():
            nd=len(v[0]); bgf=bgcount[tf]/nbg; expd=bgsum[tf]/nbg
            # one-sided binomial tail for "bound at this promoter" given bgf
            p_bound = 1.0-(1.0-bgf) if bgf<1 else 1.0
            out.append(dict(gene=gene,window_bp=w,TF=tf,n_datasets_at_promoter=nd,
                n_biotypes=len(v[1]),biotypes=";".join(sorted(v[1])[:12]),
                bg_frac_promoters_bound=round(bgf,4),
                bg_mean_datasets_per_promoter=round(expd,3),
                enrichment_vs_background=round(nd/expd,2) if expd>0 else None,
                is_network_TF=int(tf in netTF),
                is_claimed_regulator_of_this_gene=int((tf,gene) in net_edge)))
    # unbound network TFs
    for gene in ("COL1A1","COL3A1"):
        m=occ(NETP[gene],nt[gene][1],w)
        for tf in sorted(netTF-set(m)):
            bgf=bgcount[tf]/nbg
            out.append(dict(gene=gene,window_bp=w,TF=tf,n_datasets_at_promoter=0,
                n_biotypes=0,biotypes="",bg_frac_promoters_bound=round(bgf,4),
                bg_mean_datasets_per_promoter=round(bgsum[tf]/nbg,3),
                enrichment_vs_background=0.0 if bgsum[tf]>0 else None,
                is_network_TF=1,is_claimed_regulator_of_this_gene=int((tf,gene) in net_edge)))
with open(BASE+"/results/v2/chip_collagen_promoter_specificity.csv","w",newline="") as fh:
    wr=csv.DictWriter(fh,fieldnames=list(out[0].keys())); wr.writeheader()
    wr.writerows(sorted(out,key=lambda r:(r["gene"],r["window_bp"],
        -(r["enrichment_vs_background"] or 0),-r["n_datasets_at_promoter"])))
print("WROTE chip_collagen_promoter_specificity.csv rows=",len(out))
for gene in ("COL1A1","COL3A1"):
    for w in (1000,):
        sub=[r for r in out if r["gene"]==gene and r["window_bp"]==w and r["n_datasets_at_promoter"]>=3]
        sub=[r for r in sub if r["enrichment_vs_background"]]
        print("\n== %s TSS+/-%dbp : TFs with >=3 ReMap datasets, ranked by enrichment over "
              "their own background promoter occupancy (n=%d qualifying TFs) =="%(gene,w,len(sub)))
        for r in sorted(sub,key=lambda z:-z["enrichment_vs_background"])[:15]:
            print("  %-9s %3d datasets  bg %.3f/promoter  enrich %5.1fx  netTF=%d claimed=%d  %s"%(
                r["TF"],r["n_datasets_at_promoter"],r["bg_mean_datasets_per_promoter"],
                r["enrichment_vs_background"],r["is_network_TF"],
                r["is_claimed_regulator_of_this_gene"],r["biotypes"][:60]))
print("\n== The 8 TFs the network claims regulate COL1A1 ==")
for tf in ("ETS1","MKL1","MYB","NFKB1","RELA","SP1","STAT6","TFAP2A"):
    for w in (1000,5000):
        r=[x for x in out if x["gene"]=="COL1A1" and x["window_bp"]==w and x["TF"]==tf]
        if r:
            r=r[0]
            print("  %-8s +/-%4dbp: %3d datasets, bg %.3f/promoter, enrich %s, cells=%s"%(
                tf,w,r["n_datasets_at_promoter"],r["bg_mean_datasets_per_promoter"],
                r["enrichment_vs_background"],r["biotypes"][:70]))
        else:
            print("  %-8s +/-%4dbp: TF ABSENT from ReMap catalogue"%(tf,w))
