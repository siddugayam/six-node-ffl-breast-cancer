#!/usr/bin/env python3
"""22c -- map the Enformer in-silico-mutagenesis peaks onto the JASPAR/HOCOMOCO motif sites
called by 28_, so that the model's most load-bearing promoter positions are named."""
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
ism=pd.read_csv(f"{RES}/seqreg_enformer_ism.csv")
sites=pd.read_csv(f"{RES}/seqreg_motif_sites.csv")
sites["len"]=sites.seq.astype(str).str.len()
rows=[]
for reg in ism.region.unique():
    s=ism[ism.region==reg]
    pos=s.groupby("rel_to_TSS").agg(mean_d_fib_cage=("d_fib_cage","mean"),
                                    min_d_fib_cage=("d_fib_cage","min"),
                                    max_abs=("d_fib_cage",lambda x: float(np.abs(x).max())),
                                    mean_d_fib_dnase=("d_fib_dnase","mean")).reset_index()
    ss=sites[sites.region==reg]
    for _,p in pos.iterrows():
        r=int(p.rel_to_TSS)
        ov=ss[(ss.rel_to_TSS<=r)&(ss.rel_to_TSS+ss["len"]>=r)]
        rows.append(dict(region=reg, rel_to_TSS=r,
                         mean_d_fib_cage=round(p.mean_d_fib_cage,4),
                         min_d_fib_cage=round(p.min_d_fib_cage,4),
                         max_abs_d_fib_cage=round(p.max_abs,4),
                         mean_d_fib_dnase=round(p.mean_d_fib_dnase,4),
                         n_motifs_covering=len(ov),
                         best_motif=(ov.nsmallest(1,"pvalue").motif.iloc[0] if len(ov) else ""),
                         best_motif_p=(float(ov.pvalue.min()) if len(ov) else np.nan),
                         motifs=";".join(sorted(set(ov.motif.str.split("|").str[-1].str.split("_HUMAN").str[0])))[:200]))
P=pd.DataFrame(rows)
P.to_csv(f"{RES}/seqreg_enformer_ism_by_position.csv", index=False)
pd.set_option("display.width",250)
for reg in P.region.unique():
    s=P[P.region==reg]
    print(f"\n===== {reg}: 20 most load-bearing promoter positions (Enformer fibroblast CAGE)")
    print(s.nlargest(20,"max_abs_d_fib_cage")[["rel_to_TSS","mean_d_fib_cage","min_d_fib_cage",
          "mean_d_fib_dnase","best_motif","best_motif_p","motifs"]].to_string(index=False))
