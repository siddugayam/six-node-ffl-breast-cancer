#!/usr/bin/env python3
"""28c -- Benjamini-Hochberg correction of the focus-motif calibration (25 motifs x 8 promoters
= 200 tests), plus a note of which called sites are the same genomic site detected by several
matrices (they are not independent tests)."""
import pandas as pd, numpy as np
from scipy.stats import false_discovery_control
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
d=pd.read_csv(f"{RES}/seqreg_motif_focus_calibrated_5000.csv")
d["p_used"]=d[["emp_p_gcmatched","poisson_p_ge"]].max(axis=1)   # conservative
d["BH_q"]=false_discovery_control(d.p_used.values, method="bh")
d=d.sort_values("p_used")
d.to_csv(f"{RES}/seqreg_motif_focus_calibrated_BH.csv", index=False)
pd.set_option("display.width",240)
print("n tests:",len(d))
print(d.head(15)[["region","db","tf","n_sites","bg_gcmatched_mean","emp_p_gcmatched","poisson_p_ge","p_used","BH_q"]].to_string(index=False))
print("\nany BH q<0.10:", int((d.BH_q<0.10).sum()), "  any BH q<0.25:", int((d.BH_q<0.25).sum()))
# collapse to distinct genomic sites for the NF-kB family at MIR29B2CHG
s=pd.read_csv(f"{RES}/seqreg_motif_sites.csv")
s["tf"]=s.motif.str.split("|").str[-1].str.split("_HUMAN").str[0]
nk=s[(s.region=="MIR29B2CHG")&(s.tf.isin(["NFKB1","NFKB2","RELA","RELB","REL"]))]
print("\nNF-kB-family sites called at the MIR29B2CHG promoter (distinct positions collapse):")
print(nk[["db","motif","rel_to_TSS","strand","score","pvalue","seq"]].sort_values("rel_to_TSS").to_string(index=False))
print("distinct start positions:", sorted(nk.rel_to_TSS.unique()))
