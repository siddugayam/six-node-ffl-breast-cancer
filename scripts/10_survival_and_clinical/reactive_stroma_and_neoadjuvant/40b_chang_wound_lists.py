#!/usr/bin/env python3
"""Extract the tumour-level wound-response gene lists of Chang et al. 2004 (PLoS Biol,
PMID 14737219) from Dataset S2 sheets 10 (Sorlie cohort) and 11 (van't Veer cohort):
the SAM significant genes separating 'wound-activated' from 'quiescent' breast tumours."""
import pandas as pd, re, warnings, os
warnings.filterwarnings("ignore")
BASE = "/path/to/revision"
F = os.path.join(BASE, "cache/v5/farmer/sup/chang2004/pbio.0020007.sd002.xls")
OUT = os.path.join(BASE, "cache/v5/farmer/lists")
xl = pd.ExcelFile(F)

def parse(sheet, symbol_from):
    d = xl.parse(sheet, header=None)
    pos, neg, mode = [], [], None
    for i in range(d.shape[0]):
        c0 = str(d.iloc[i, 0])
        if "Positive Significant Genes" in c0: mode = "pos"; continue
        if "Negative Significant Genes" in c0: mode = "neg"; continue
        gid = d.iloc[i, 3]
        if mode is None or not isinstance(gid, str): continue
        if symbol_from == "second_token":
            toks = gid.split()
            sym = toks[1] if len(toks) > 1 else None
        else:
            m = re.findall(r"\(([A-Za-z0-9_\-\.]+)\)\s*$", gid.strip())
            sym = m[-1] if m else None
        if not sym or sym in ("ID", "Gene"): continue
        (pos if mode == "pos" else neg).append(sym)
    return list(dict.fromkeys(pos)), list(dict.fromkeys(neg))

sp, sn = parse("10. SAM_Sorlie", "second_token")
vp, vn = parse("11. SAM_vant_Veer", "paren")
for nm, g in [("CHANG_WOUND_ACTIVATED_UP_SORLIE", sp), ("CHANG_WOUND_ACTIVATED_DN_SORLIE", sn),
              ("CHANG_WOUND_ACTIVATED_UP_VANTVEER", vp), ("CHANG_WOUND_ACTIVATED_DN_VANTVEER", vn)]:
    open(os.path.join(OUT, nm + ".txt"), "w").write("\n".join(g) + "\n")
    print(f"{nm:36s} n={len(g):4d}  e.g. {g[:8]}")
