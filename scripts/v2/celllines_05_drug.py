#!/usr/bin/env python3
"""
celllines_05_drug.py -- part E of the stroma-free test.  EXPLORATORY.

Does a collagen module score, or a miR-29-target module score, predict drug sensitivity
across breast cancer cell lines?  Two independent panels:
  PRISM secondary screen (AUC; depmap_id joins directly)
  GDSC2 fitted dose response (LN_IC50 and AUC; SANGER_MODEL_ID -> ModelID via Model.csv)
Per-compound Spearman, BH-FDR across compounds within each panel and each score.
Only FDR-corrected results are reported.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
PRISM = "/path/to/home/Desktop/DD/R_GPR/data/depmap/prism_secondary_dose_response.csv"
GDSC = "/path/to/home/Desktop/DD/R_GPR/data/depmap/GDSC2_dose_response.xlsx"
MIN_N = 20

expr = pd.read_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", index_col=0)
models = pd.read_csv(os.path.join(RES, "celllines_breast_models.csv")).set_index("ModelID")
cancer_ids = [i for i in expr.index if models.loc[i, "OncotreePrimaryDisease"] != "Non-Cancerous"]
E = expr.loc[cancer_ids]
edges = pd.read_csv(os.path.join(DATA, "canonical_edges.tsv"), sep="\t")

z = lambda df: (df - df.mean()) / df.std(ddof=1)
COLLAGEN = [g for g in ["COL1A1", "COL3A1", "COL1A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
                        "COL6A3", "COL4A1", "COL11A1", "FN1", "SPARC", "LOX"] if g in E.columns]
mir29_targets = sorted(set(edges.loc[edges.source.isin(["hsa-miR-29a", "hsa-miR-29b", "hsa-miR-29c"])
                                     & (edges.edge_type == "miRNA_target"), "target"]) & set(E.columns))
MES = [g for g in ["VIM", "ZEB1", "ZEB2", "SNAI2", "TWIST1", "CDH2"] if g in E.columns]
scores = pd.DataFrame({
    "collagen_module": z(E[COLLAGEN]).mean(axis=1),
    "miR29_target_module": z(E[mir29_targets]).mean(axis=1),
    "mesenchymal_score": z(E[MES]).mean(axis=1),
    "COL1A1": E["COL1A1"], "COL3A1": E["COL3A1"]})
scores.to_csv(os.path.join(RES, "celllines_E_module_scores.csv"))
print(f"collagen module {len(COLLAGEN)} genes; miR-29 target module {len(mir29_targets)} genes; "
      f"mesenchymal {len(MES)} genes; {len(scores)} breast cancer lines")
print("score intercorrelation (Spearman):")
print(scores.corr(method="spearman").round(3).to_string())

out = []


def test_panel(panel, mat, score_names, metric):
    """mat: DataFrame index=ModelID, columns=compound, values=metric."""
    common = [i for i in mat.index if i in scores.index]
    mat = mat.loc[common]
    print(f"{panel}/{metric}: {len(common)} breast lines, {mat.shape[1]} compounds")
    for sname in score_names:
        s = scores.loc[common, sname].values
        recs = []
        for c in mat.columns:
            y = mat[c].values
            m = np.isfinite(y) & np.isfinite(s)
            if m.sum() < MIN_N:
                continue
            r, p = stats.spearmanr(s[m], y[m])
            if not np.isfinite(p):
                continue
            recs.append((c, int(m.sum()), float(r), float(p)))
        if not recs:
            continue
        d = pd.DataFrame(recs, columns=["compound", "n", "rho", "p"])
        d["FDR_BH"] = stats.false_discovery_control(d["p"].values, method="bh")
        d.insert(0, "metric", metric)
        d.insert(0, "score", sname)
        d.insert(0, "panel", panel)
        d["n_compounds_tested"] = len(d)
        out.append(d)
        nsig = int((d.FDR_BH < 0.10).sum())
        print(f"   {sname}: {len(d)} compounds tested, {int((d.FDR_BH<0.05).sum())} at FDR<0.05, "
              f"{nsig} at FDR<0.10")


# ------------------------------------------------------------------ PRISM
pr = pd.read_csv(PRISM, usecols=["depmap_id", "name", "auc", "ic50", "moa", "screen_id"],
                 low_memory=False)
pr = pr[pr.depmap_id.isin(scores.index)]
print("PRISM rows for breast lines:", len(pr), "| distinct compounds:", pr.name.nunique())
moa_map = pr.dropna(subset=["moa"]).drop_duplicates("name").set_index("name")["moa"].to_dict()
auc = pr.pivot_table(index="depmap_id", columns="name", values="auc", aggfunc="median")
test_panel("PRISM_secondary", auc, ["collagen_module", "miR29_target_module",
                                    "mesenchymal_score", "COL1A1"], "AUC")

# ------------------------------------------------------------------ GDSC2
gd = pd.read_excel(GDSC, usecols=["SANGER_MODEL_ID", "DRUG_NAME", "DRUG_ID", "PUTATIVE_TARGET",
                                  "PATHWAY_NAME", "LN_IC50", "AUC", "TCGA_DESC"])
mdl = pd.read_csv(os.path.join(DATA, "depmap", "Model_24Q4.csv"), low_memory=False)
s2m = mdl.dropna(subset=["SangerModelID"]).set_index("SangerModelID")["ModelID"].to_dict()
gd["ModelID"] = gd["SANGER_MODEL_ID"].map(s2m)
gd = gd[gd.ModelID.isin(scores.index)]
print("GDSC2 rows for breast lines:", len(gd), "| distinct drugs:", gd.DRUG_NAME.nunique(),
      "| distinct lines:", gd.ModelID.nunique())
gd_target = gd.drop_duplicates("DRUG_NAME").set_index("DRUG_NAME")[
    ["PUTATIVE_TARGET", "PATHWAY_NAME"]]
for metric in ["LN_IC50", "AUC"]:
    mat = gd.pivot_table(index="ModelID", columns="DRUG_NAME", values=metric, aggfunc="median")
    test_panel("GDSC2", mat, ["collagen_module", "miR29_target_module",
                              "mesenchymal_score", "COL1A1"], metric)

res = pd.concat(out, ignore_index=True)
res["moa_or_target"] = res["compound"].map(moa_map)
res.loc[res.panel == "GDSC2", "moa_or_target"] = res.loc[res.panel == "GDSC2", "compound"].map(
    gd_target["PUTATIVE_TARGET"])
res["gdsc_pathway"] = res["compound"].map(gd_target["PATHWAY_NAME"])
res = res.sort_values(["panel", "score", "metric", "p"])
res.to_csv(os.path.join(RES, "celllines_E_drug_associations.csv"), index=False)

sig = res[res.FDR_BH < 0.10]
print("\n=== FDR<0.10 hits ===")
print(sig.to_string(index=False) if len(sig) else "NONE at FDR<0.10 in any panel/score")
sig.to_csv(os.path.join(RES, "celllines_E_drug_hits_FDR10.csv"), index=False)

# summary row per panel/score
summ = (res.groupby(["panel", "score", "metric"])
           .agg(n_compounds=("compound", "size"),
                n_FDR05=("FDR_BH", lambda v: int((v < 0.05).sum())),
                n_FDR10=("FDR_BH", lambda v: int((v < 0.10).sum())),
                min_FDR=("FDR_BH", "min"))
           .reset_index())
summ.to_csv(os.path.join(RES, "celllines_E_drug_summary.csv"), index=False)
print("\n", summ.to_string(index=False))
print("DONE")
