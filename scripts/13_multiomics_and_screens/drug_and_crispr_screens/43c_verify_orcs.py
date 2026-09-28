#!/usr/bin/env python3
"""Independent verification of PART 1D (BioGRID ORCS 2.0.18)."""
import pandas as pd, numpy as np
REV="/path/to/revision"; ORCS=f"{REV}/cache/v7/orcs"; OUT=f"{REV}/results/v3"
idx=pd.read_csv(f"{ORCS}/BIOGRID-ORCS-SCREEN_INDEX-2.0.18.index.tab.txt",sep="\t",dtype=str).rename(columns={'#SCREEN_ID':'SCREEN_ID'})
print("total human screens in index:",len(idx),"| organisms:",idx.ORGANISM_OFFICIAL.value_counts().to_dict())
isb=idx.CELL_TYPE.str.lower().str.contains("breast|mammary",na=False)
print("breast/mammary screens:",int(isb.sum()))
print(idx[isb].PHENOTYPE.value_counts().to_string())
mig=idx[idx.PHENOTYPE.str.contains("migration",case=False,na=False)]
print("\nmigration-phenotype screens anywhere in human ORCS:",len(mig))
print(mig[["SCREEN_ID","AUTHOR","CELL_LINE","CELL_TYPE","PHENOTYPE","SOURCE_ID"]].to_string(index=False))
print("\nof these, breast:",int(mig.CELL_TYPE.str.lower().str.contains('breast|mammary',na=False).sum()))
free=(idx.CONDITION_NAME.fillna('')+' '+idx.NOTES.fillna('')+' '+idx.SCREEN_RATIONALE.fillna('')+' '+idx.SCREEN_NAME.fillna(''))
inv=idx[free.str.contains("invasi|matrigel|transwell|extracellular matrix|metasta",case=False,na=False)]
print("free-text invasion/matrigel/transwell/ECM/metastasis screens:",len(inv),"| breast among them:",
      int(inv.CELL_TYPE.str.lower().str.contains('breast|mammary',na=False).sum()))
dep=pd.read_csv(f"{OUT}/screens_orcs_screen_summary.csv")
rows=[]
for sid in [188,220,222]:
    d=pd.read_csv(f"{ORCS}/BIOGRID-ORCS-SCREEN_{sid}-2.0.18.screen.tab.txt",sep="\t",dtype=str,low_memory=False).rename(columns={'#SCREEN_ID':'SCREEN_ID'})
    hit=d["HIT"].str.upper()=="YES"; r=dep[dep.screen_id==sid].iloc[0]
    print(f"\nscreen {sid}: unique genes {d.OFFICIAL_SYMBOL.nunique()} hits {int(hit.sum())} | deposited genes {r.n_genes_tested} hits {r.n_hits}")
    for g in ["MYC","NFKB1","COL1A1","COL3A1","SP1","RELA","ETS1","STAT3","HIF1A"]:
        s=d[d.OFFICIAL_SYMBOL==g]
        if len(s): rows.append(dict(screen=sid,gene=g,HIT=s.HIT.iloc[0],score1=s["SCORE.1"].iloc[0]))
V=pd.DataFrame(rows); print("\nspot-check hit calls:"); print(V.pivot(index="gene",columns="screen",values="HIT").to_string())
V.to_csv(f"{OUT}/screens_verify_orcs_spotcheck.csv",index=False)
