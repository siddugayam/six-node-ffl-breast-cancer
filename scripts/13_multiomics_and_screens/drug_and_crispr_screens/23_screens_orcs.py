#!/usr/bin/env python3
"""PART 1D) Published phenotypic CRISPR screens (BioGRID ORCS 2.0.18, human).
This script asks (i) what such screens actually exist in breast models, (ii) whether any
network hub scores as a hit in them, and (iii) what the closest migration/invasion
screens are anywhere in the database.
"""
import os, glob, gzip, sys
import pandas as pd, numpy as np
from scipy import stats

REV="/path/to/revision"
ORCS=f"{REV}/cache/v7/orcs"; OUT=f"{REV}/results/v3"
def msg(*a): print(*a, flush=True)

idx = pd.read_csv(f"{ORCS}/BIOGRID-ORCS-SCREEN_INDEX-2.0.18.index.tab.txt", sep="\t", dtype=str)
idx = idx.rename(columns={'#SCREEN_ID':'SCREEN_ID'})
msg(f"BioGRID ORCS 2.0.18 human screens in index: {len(idx)}")
msg("phenotype counts:\n" + idx.PHENOTYPE.value_counts().to_string())

is_breast = idx.CELL_TYPE.str.lower().str.contains('breast|mammary', na=False)
msg(f"\nbreast/mammary screens: {is_breast.sum()}")
msg(idx[is_breast].PHENOTYPE.value_counts().to_string())

# the only migration/invasion-phenotype screens in the whole human database
mig = idx[idx.PHENOTYPE.str.contains('migration', case=False, na=False)]
msg(f"\nscreens with a MIGRATION phenotype anywhere in human ORCS: {len(mig)}")
for _,r in mig.iterrows():
    msg(f"   screen {r.SCREEN_ID} | {r.AUTHOR} | {r.CELL_LINE} ({r.CELL_TYPE}) | PMID {r.SOURCE_ID} | setup: {str(r.EXPERIMENTAL_SETUP)[:80]} | condition: {str(r.CONDITION_NAME)[:90]}")

# invasion / ECM / matrigel mentioned anywhere in the free-text fields
free = (idx.CONDITION_NAME.fillna('')+' '+idx.NOTES.fillna('')+' '+idx.SCREEN_RATIONALE.fillna('')+' '+idx.SCREEN_NAME.fillna(''))
inv = idx[free.str.contains('invasi|matrigel|transwell|extracellular matrix|collagen|metasta', case=False, na=False)]
msg(f"\nscreens whose free text mentions invasion/matrigel/transwell/ECM/collagen/metastasis: {len(inv)}")
for _,r in inv.iterrows():
    msg(f"   screen {r.SCREEN_ID} | {r.AUTHOR} | {r.CELL_LINE} ({r.CELL_TYPE}) | {r.PHENOTYPE} | PMID {r.SOURCE_ID}")

nodes = pd.read_csv(f"{REV}/data/canonical_nodes.tsv", sep="\t")
ctrl  = pd.read_csv(f"{OUT}/systems_controllability_nodes.csv")
hubs  = set(ctrl.loc[ctrl.is_FFL_hub==True,'node'])
prot  = set(nodes.loc[nodes.type.isin(['Gene','TF']),'name'])
focus = ['NFKB1','RELA','SP1','ETS1','COL1A1','COL3A1','MYC','TP53','STAT3','HIF1A','MKL1','VEGFA','MMP2','TGFB1','ESR1','XIAP']
msg(f"\nnetwork protein-coding nodes: {len(prot)} | FFL hubs: {len(hubs)} "
    f"| protein-coding FFL hubs: {len(hubs & prot)}")

def load_screen(sid):
    p=f"{ORCS}/BIOGRID-ORCS-SCREEN_{sid}-2.0.18.screen.tab.txt"
    if not os.path.exists(p): return None
    d=pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    d=d.rename(columns={'#SCREEN_ID':'SCREEN_ID'})
    return d

# ---- which screens to interrogate ----
sel = idx[is_breast].copy()
sel['group']=np.where(sel.PHENOTYPE.str.contains('proliferation|viability',case=False,na=False),'proliferation',
              np.where(sel.PHENOTYPE.str.contains('chemicals|radiation',case=False,na=False),'drug_response',
              np.where(sel.PHENOTYPE.str.contains('tumorigenicity',case=False,na=False),'tumorigenicity','other')))
extra = idx[idx.SCREEN_ID.isin(mig.SCREEN_ID)]
todo = pd.concat([sel, extra.assign(group='migration_nonbreast')]).drop_duplicates('SCREEN_ID')
msg(f"\nscreens to interrogate: {len(todo)}  ({todo.group.value_counts().to_dict()})")

rows=[]; hitrows=[]
for _,r in todo.iterrows():
    d = load_screen(r.SCREEN_ID)
    if d is None:
        msg(f"   screen {r.SCREEN_ID}: file missing"); continue
    d['sym']=d.OFFICIAL_SYMBOL.astype(str)
    n_tested=d.sym.nunique(); n_hit=d.loc[d.HIT=='YES','sym'].nunique()
    tested=set(d.sym); hitset=set(d.loc[d.HIT=='YES','sym'])
    ph = prot & tested; hh = (hubs & prot) & tested
    ph_hit = ph & hitset; hh_hit = hh & hitset
    # Fisher: are network nodes / FFL hubs enriched among hits?
    def fisher(sub, subhit):
        a=len(subhit); b=len(sub)-a; c=n_hit-a; dd=n_tested-n_hit-b
        if min(a+b, c+dd)==0: return np.nan, np.nan
        orr,p = stats.fisher_exact([[a,b],[c,dd]], alternative='greater')
        return orr,p
    or_n,p_n = fisher(ph, ph_hit); or_h,p_h = fisher(hh, hh_hit)
    rows.append(dict(screen_id=r.SCREEN_ID, pmid=r.SOURCE_ID, author=r.AUTHOR, cell_line=r.CELL_LINE,
                     cell_type=r.CELL_TYPE, phenotype=r.PHENOTYPE, group=r.group,
                     screen_type=r.SCREEN_TYPE, condition=r.CONDITION_NAME,
                     n_genes_tested=n_tested, n_hits=n_hit, hit_rate=n_hit/max(n_tested,1),
                     network_nodes_tested=len(ph), network_nodes_hit=len(ph_hit),
                     ffl_hubs_tested=len(hh), ffl_hubs_hit=len(hh_hit),
                     OR_network=or_n, p_network=p_n, OR_ffl_hub=or_h, p_ffl_hub=p_h,
                     focus_hits=";".join(sorted(set(focus)&hitset)),
                     ffl_hub_hits=";".join(sorted(hh_hit))))
    for g in sorted(set(focus)|(hubs&prot)):
        if g in tested:
            sc = d.loc[d.sym==g,'SCORE.1'].iloc[0]
            hitrows.append(dict(screen_id=r.SCREEN_ID, group=r.group, cell_line=r.CELL_LINE,
                                phenotype=r.PHENOTYPE, gene=g, is_focus=g in focus,
                                is_ffl_hub=g in hubs, score1=sc, hit=(g in hitset)))
R=pd.DataFrame(rows); H=pd.DataFrame(hitrows)
R.to_csv(f"{OUT}/screens_orcs_screen_summary.csv", index=False)
H.to_csv(f"{OUT}/screens_orcs_gene_level.csv", index=False)
msg(f"\nwrote {len(R)} screen summaries and {len(H)} gene-level rows")

msg("\n=== per-group summary ===")
for g,sub in R.groupby('group'):
    msg(f"{g}: {len(sub)} screens | median hit rate {sub.hit_rate.median():.4f} "
        f"| FFL hubs hit in >=1 screen: {sorted(set(sum([s.split(';') for s in sub.ffl_hub_hits if s],[])))}")

msg("\n=== how often is each focus gene a hit, by group ===")
if len(H):
    t = H[H.is_focus].groupby(['gene','group']).agg(n_screens=('hit','size'), n_hit=('hit','sum')).reset_index()
    t['frac']=t.n_hit/t.n_screens
    print(t.pivot(index='gene',columns='group',values='frac').round(3).to_string())
    print("\n(counts)")
    print(t.pivot(index='gene',columns='group',values='n_hit').to_string())
    t.to_csv(f"{OUT}/screens_orcs_focus_gene_hitrates.csv", index=False)

msg("\n=== FFL hubs ranked by how often they are a hit in breast screens ===")
if len(H):
    hb = H[(H.is_ffl_hub) & (H.group!='migration_nonbreast')]
    tt = hb.groupby('gene').agg(n_screens=('hit','size'), n_hit=('hit','sum')).reset_index()
    tt['frac_hit']=tt.n_hit/tt.n_screens
    tt=tt.sort_values('frac_hit',ascending=False)
    print(tt.head(25).to_string(index=False))
    tt.to_csv(f"{OUT}/screens_orcs_ffl_hub_hitrates.csv", index=False)
    # background hit rate over all genes in the same screens
    bg = R[R.group!='migration_nonbreast'].hit_rate.mean()
    msg(f"\nmean genome-wide hit rate in those breast screens: {bg:.4f}")
    msg(f"mean FFL-hub hit rate: {tt.n_hit.sum()/tt.n_screens.sum():.4f}")
    a=int(tt.n_hit.sum()); b=int(tt.n_screens.sum())-a
    tot_hit=int(R[R.group!='migration_nonbreast'].n_hits.sum()); tot=int(R[R.group!='migration_nonbreast'].n_genes_tested.sum())
    orr,p=stats.fisher_exact([[a,b],[tot_hit-a, tot-tot_hit-b]], alternative='greater')
    msg(f"pooled Fisher (FFL hubs vs all genes, all breast screens): OR={orr:.3f} p={p:.3g}")
msg("DONE 23")
