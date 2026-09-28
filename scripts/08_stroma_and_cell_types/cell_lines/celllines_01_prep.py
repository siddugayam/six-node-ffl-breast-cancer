#!/usr/bin/env python3
"""
celllines_01_prep.py  --  build the cell-line data objects for the "stroma-free test".

Inputs (all read-only):
  data/depmap/Model_24Q4.csv                                   DepMap 24Q4 Public model annotation
  data/depmap/OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv   DepMap 24Q4 log2(TPM+1)
  data/ccle/CCLE_miRNA_20181103.gct                            CCLE nanoString miRNA, 734 x 954
  data/ccle/Cell_lines_annotations_20181226.txt                CCLE_Name -> DepMap_ID map
  ../../TCGA_PAN_CAN/.../CRISPRGeneEffect.csv                  local DepMap CRISPR (release recorded, not assumed)

Outputs (results/v2/):
  celllines_00_provenance.csv
  celllines_breast_models.csv
  celllines_expr_breast.tsv.gz      (breast models x protein-coding genes, log2(TPM+1))
  celllines_mirna_breast.tsv        (breast models x miRNA, CCLE nanoString log2)
Nothing outside results/v2, scripts/08_stroma_and_cell_types/cell_lines, logs/v2, data/depmap, data/ccle is written.
"""
import os, sys, gzip, hashlib, json
import numpy as np
import pandas as pd

ROOT = "/path/to/revision"
RES = os.path.join(ROOT, "results", "v2")
DATA = os.path.join(ROOT, "data")
CRISPR = "/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data/CRISPRGeneEffect.csv"
os.makedirs(RES, exist_ok=True)

prov = []


def note(k, v):
    print(f"[prov] {k} = {v}", flush=True)
    prov.append({"item": k, "value": str(v)})


# ---------------------------------------------------------------- models
model = pd.read_csv(os.path.join(DATA, "depmap", "Model_24Q4.csv"), low_memory=False)
note("depmap_release", "DepMap 24Q4 Public (figshare article 27993248, published 2024-12-10)")
note("Model.csv_n_models", len(model))
breast = model[model["OncotreeLineage"] == "Breast"].copy()
note("Model.csv_n_breast_lineage", len(breast))
# carcinoma-only subset (drop the non-cancerous / immortalized normal-breast models)
breast_ca = breast[breast["OncotreePrimaryDisease"] != "Non-Cancerous"].copy()
note("n_breast_cancer_models(excl Non-Cancerous)", len(breast_ca))

# ---------------------------------------------------------------- expression
exp_path = os.path.join(DATA, "depmap", "OmicsExpressionProteinCodingGenesTPMLogp1_24Q4.csv")
sz = os.path.getsize(exp_path)
note("expression_file_bytes", sz)
assert sz == 506628654, f"expression file incomplete/unexpected size {sz} (expected 506628654)"

expr = pd.read_csv(exp_path, index_col=0, low_memory=False)
note("expression_n_models_total", expr.shape[0])
note("expression_n_genes_total", expr.shape[1])
# column names are "SYMBOL (ENTREZ)"
sym = pd.Index([c.split(" (")[0] for c in expr.columns])
note("expression_n_duplicated_symbols", int(sym.duplicated().sum()))
expr.columns = sym

br_ids = [m for m in breast["ModelID"] if m in expr.index]
note("n_breast_models_with_expression", len(br_ids))
br_ca_ids = [m for m in breast_ca["ModelID"] if m in expr.index]
note("n_breast_CANCER_models_with_expression", len(br_ca_ids))

expr_br = expr.loc[br_ids]
expr_br.to_csv(os.path.join(RES, "celllines_expr_breast.tsv.gz"), sep="\t", compression="gzip")

# per-lineage COL1A1/COL3A1 summary needs the whole matrix; compute now, then free
lin = model.set_index("ModelID")["OncotreeLineage"]
for g in ["COL1A1", "COL3A1"]:
    assert g in expr.columns, g
sub = expr[["COL1A1", "COL3A1"]].copy()
sub["lineage"] = lin.reindex(sub.index)
lin_summary = (sub.groupby("lineage")
                 .agg(n=("COL1A1", "size"),
                      COL1A1_median=("COL1A1", "median"),
                      COL1A1_max=("COL1A1", "max"),
                      COL3A1_median=("COL3A1", "median"),
                      COL3A1_max=("COL3A1", "max"))
                 .sort_values("COL1A1_median", ascending=False))
lin_summary.to_csv(os.path.join(RES, "celllines_C_lineage_collagen.csv"))
note("lineages_summarised", len(lin_summary))

# within-sample percentile rank of the two collagens, all DepMap models
rk = expr.rank(axis=1, pct=True)
pct = rk[["COL1A1", "COL3A1"]].copy()
pct["lineage"] = lin.reindex(pct.index)
pct.to_csv(os.path.join(RES, "celllines_C_percentile_all_models.csv"))
del rk

# ---------------------------------------------------------------- miRNA
gct = os.path.join(DATA, "ccle", "CCLE_miRNA_20181103.gct")
with open(gct) as fh:
    fh.readline()
    dims = fh.readline().split()
    hdr = fh.readline().rstrip("\n").split("\t")
note("ccle_mirna_file", "CCLE_miRNA_20181103.gct (data.broadinstitute.org/ccle, nanoString nCounter)")
note("ccle_mirna_dims", f"{dims[0]} miRNA x {dims[1]} cell lines")
mir = pd.read_csv(gct, sep="\t", skiprows=2)
mir = mir.set_index("Description").drop(columns=["Name"])
note("ccle_mirna_n_probes", mir.shape[0])
note("ccle_mirna_n_lines", mir.shape[1])
note("ccle_mirna_value_range", f"{np.nanmin(mir.values):.3f} .. {np.nanmax(mir.values):.3f}")

ann = pd.read_csv(os.path.join(DATA, "ccle", "Cell_lines_annotations_20181226.txt"),
                  sep="\t", low_memory=False)
# CCLEName -> ModelID from the 24Q4 Model file (authoritative), fall back to legacy annotation
ccle2id = model.dropna(subset=["CCLEName"]).set_index("CCLEName")["ModelID"].to_dict()
note("Model.csv_CCLEName_mapped", len(ccle2id))

mir_cols = list(mir.columns)
mir_breast_cols = [c for c in mir_cols if c.endswith("_BREAST")]
note("ccle_mirna_n_breast_columns", len(mir_breast_cols))
mapped = {c: ccle2id[c] for c in mir_breast_cols if c in ccle2id}
note("ccle_mirna_breast_mapped_to_ModelID", len(mapped))
unmapped = [c for c in mir_breast_cols if c not in ccle2id]
note("ccle_mirna_breast_unmapped", ";".join(unmapped) if unmapped else "none")

mir_br = mir[list(mapped.keys())].T
mir_br.index = [mapped[c] for c in mir_br.index]
mir_br.index.name = "ModelID"
# keep only breast-lineage models that are also in Model.csv Breast set
keep = [i for i in mir_br.index if i in set(breast["ModelID"])]
mir_br = mir_br.loc[keep]
note("n_breast_models_with_miRNA", mir_br.shape[0])
both = [i for i in mir_br.index if i in set(expr_br.index)]
note("n_breast_models_with_miRNA_AND_expression", len(both))
mir_br.to_csv(os.path.join(RES, "celllines_mirna_breast.tsv"), sep="\t")

# ---------------------------------------------------------------- CRISPR provenance
with open(CRISPR) as fh:
    h = fh.readline()
ncol = h.count(",") + 1
nrow = sum(1 for _ in open(CRISPR))
note("crispr_file_path", CRISPR)
note("crispr_file_bytes", os.path.getsize(CRISPR))
note("crispr_n_models", nrow - 1)
note("crispr_n_genes", ncol - 1)
note("crispr_release_note",
     "local copy; byte size does not match the 22Q4/23Q2/23Q4/24Q2/24Q4 public "
     "CRISPRGeneEffect.csv files, so the exact release is not identifiable from the file itself")

# ---------------------------------------------------------------- breast model table
cols = ["ModelID", "StrippedCellLineName", "CCLEName", "OncotreeLineage",
        "OncotreePrimaryDisease", "OncotreeSubtype", "PrimaryOrMetastasis", "GrowthPattern"]
bm = breast[cols].copy()
bm["has_expression"] = bm["ModelID"].isin(expr_br.index)
bm["has_miRNA"] = bm["ModelID"].isin(mir_br.index)
crispr_ids = set(pd.read_csv(CRISPR, usecols=[0]).iloc[:, 0])
bm["has_CRISPR"] = bm["ModelID"].isin(crispr_ids)
bm["COL1A1_log2TPM1"] = expr_br["COL1A1"].reindex(bm["ModelID"]).values
bm["COL3A1_log2TPM1"] = expr_br["COL3A1"].reindex(bm["ModelID"]).values
bm.to_csv(os.path.join(RES, "celllines_breast_models.csv"), index=False)
note("n_breast_models_with_CRISPR", int(bm["has_CRISPR"].sum()))
note("n_breast_models_expr_and_CRISPR", int((bm["has_CRISPR"] & bm["has_expression"]).sum()))

pd.DataFrame(prov).to_csv(os.path.join(RES, "celllines_00_provenance.csv"), index=False)
print("DONE")
