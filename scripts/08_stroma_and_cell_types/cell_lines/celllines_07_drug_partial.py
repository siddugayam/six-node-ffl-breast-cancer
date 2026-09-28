#!/usr/bin/env python3
"""
celllines_07_drug_partial.py

Follow-up to part E.  Every GDSC2 hit for the collagen module score is also a hit for a
generic mesenchymal score (the two scores correlate rho = 0.75 across the 69 breast
cancer lines).  This script asks whether the collagen module adds anything once the
mesenchymal phenotype is held fixed, and vice versa.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
GDSC = "/path/to/home/Desktop/DD/R_GPR/data/depmap/GDSC2_dose_response.xlsx"

scores = pd.read_csv(os.path.join(RES, "celllines_E_module_scores.csv"), index_col=0)
hits = pd.read_csv(os.path.join(RES, "celllines_E_drug_hits_FDR10.csv"))
drugs = sorted(set(hits.loc[hits.panel == "GDSC2", "compound"]))
print("GDSC2 drugs to re-test:", drugs)

gd = pd.read_excel(GDSC, usecols=["SANGER_MODEL_ID", "DRUG_NAME", "LN_IC50", "AUC"])
mdl = pd.read_csv(os.path.join(DATA, "depmap", "Model_24Q4.csv"), low_memory=False)
s2m = mdl.dropna(subset=["SangerModelID"]).set_index("SangerModelID")["ModelID"].to_dict()
gd["ModelID"] = gd["SANGER_MODEL_ID"].map(s2m)
gd = gd[gd.ModelID.isin(scores.index) & gd.DRUG_NAME.isin(drugs)]


def partial_spear(x, y, z):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[m], y[m], z[m]
    rx, ry, rz = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    Z = np.column_stack([np.ones(len(rz)), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p), int(m.sum())


rows = []
for metric in ["AUC", "LN_IC50"]:
    mat = gd.pivot_table(index="ModelID", columns="DRUG_NAME", values=metric, aggfunc="median")
    common = [i for i in mat.index if i in scores.index]
    mat = mat.loc[common]
    coll = scores.loc[common, "collagen_module"].values
    mes = scores.loc[common, "mesenchymal_score"].values
    for d in mat.columns:
        y = mat[d].values
        m = np.isfinite(y)
        if m.sum() < 20:
            continue
        r_c, p_c = stats.spearmanr(coll[m], y[m])
        r_m, p_m = stats.spearmanr(mes[m], y[m])
        rc_pm, pc_pm, n = partial_spear(coll, y, mes)
        rm_pc, pm_pc, _ = partial_spear(mes, y, coll)
        rows.append(dict(metric=metric, drug=d, n=n,
                         rho_collagen=r_c, p_collagen=p_c,
                         rho_mesenchymal=r_m, p_mesenchymal=p_m,
                         rho_collagen_given_MES=rc_pm, p_collagen_given_MES=pc_pm,
                         rho_MES_given_collagen=rm_pc, p_MES_given_collagen=pm_pc))
out = pd.DataFrame(rows).sort_values(["metric", "p_collagen"])
out.to_csv(os.path.join(RES, "celllines_E_drug_partial_vs_mesenchymal.csv"), index=False)
print(out.round(4).to_string(index=False))
print("\ncollagen retains p<0.05 after adjusting for mesenchymal score:",
      int((out.p_collagen_given_MES < 0.05).sum()), "of", len(out))
print("mesenchymal retains p<0.05 after adjusting for collagen score:",
      int((out.p_MES_given_collagen < 0.05).sum()), "of", len(out))
print("DONE")
