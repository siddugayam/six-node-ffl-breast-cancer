#!/usr/bin/env python3
"""Assemble the deliverable HPA tables (tissue validation)."""
import csv, os, collections, math

R  = "/path/to/revision/results/v6"
D  = "/path/to/revision/data/hpa"

TFS   = ["E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"]
GENES = ["CCND2","COL1A1","STAT5A","MYBL2","FN1","PDGFRB","MET","CXCL12","MMP14","PLAU"]
EXTRA = ["COL3A1","POSTN","NFKB1","RELA","SP1","ETS1"]
ROLE  = {g:"TF (prioritised)" for g in TFS}
ROLE.update({g:"gene (prioritised)" for g in GENES})
ROLE.update({g:"comparator (compartment test)" for g in EXTRA})
ALL = TFS + GENES + EXTRA

def rd(fn, base=R):
    return list(csv.DictReader(open(os.path.join(base, fn), newline="", encoding="utf-8")))

def wr(fn, rows, fields=None):
    if not rows:
        print("EMPTY", fn); return
    fields = fields or list(rows[0].keys())
    with open(os.path.join(R, fn), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(fn, len(rows))

# ---------------------------------------------------------------- inputs
cnt   = rd("hpa_cancer_ihc_counts_v25.csv")
norm  = rd("hpa_normal_tissue_ihc_v25.csv")
surv  = rd("hpa_pathology_survival_v25.csv")
subl  = rd("hpa_subcellular_v25.csv")
abrel = rd("hpa_antibody_reliability_v25.csv")
summ  = {r["gene"]: r for r in rd("hpa_gene_summaries_v25.csv")}
scspec= {r["gene"]: r for r in rd("hpa_singlecell_specificity_v25.csv")}
sc    = rd("hpa_singlecell_rna_v25.csv")

# v23 pathology.tsv = HPA's own per-gene selected-antibody breast counts
p23 = {}
with open(os.path.join(D, "pathology.tsv"), newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["Cancer"] == "breast cancer" and r["Gene name"] in ALL:
            p23[r["Gene name"]] = r

# ---------------------------------------------------- B) breast-cancer IHC
bc = [r for r in cnt if r["cancer"] == "Breast cancer"]
by_gene = collections.defaultdict(list)
for r in bc: by_gene[r["gene"]].append(r)

rows = []
for g in ALL:
    abs_ = by_gene.get(g, [])
    sel = None
    ref = p23.get(g)
    if ref:
        for a in abs_:
            if (int(a["high"]), int(a["medium"]), int(a["low"]), int(a["not_detected"])) == \
               (int(ref["High"]), int(ref["Medium"]), int(ref["Low"]), int(ref["Not detected"])):
                sel = a; break
    if sel is None and abs_:
        sel = max(abs_, key=lambda a: int(a["n_patients"]))
    tot = {k: sum(int(a[k]) for a in abs_) for k in ("high","medium","low","not_detected")} if abs_ else None
    rows.append(dict(
        gene=g, role=ROLE[g],
        n_antibodies_with_breast_IHC=len(abs_),
        selected_antibody=(sel["antibody"] if sel else ""),
        selected_matches_HPA_pathology_table=("n/a (no tissue IHC)" if not abs_ else "yes" if (ref and sel and
            (int(sel["high"]),int(sel["medium"]),int(sel["low"]),int(sel["not_detected"])) ==
            (int(ref["High"]),int(ref["Medium"]),int(ref["Low"]),int(ref["Not detected"]))) else ("no" if ref else "n/a")),
        high=(sel["high"] if sel else ""), medium=(sel["medium"] if sel else ""),
        low=(sel["low"] if sel else ""), not_detected=(sel["not_detected"] if sel else ""),
        n_patients=(sel["n_patients"] if sel else 0),
        pct_high_or_medium=(round(100*(int(sel["high"])+int(sel["medium"]))/int(sel["n_patients"]),1) if sel and int(sel["n_patients"]) else ""),
        all_antibody_pooled_high=(tot["high"] if tot else ""),
        all_antibody_pooled_medium=(tot["medium"] if tot else ""),
        all_antibody_pooled_low=(tot["low"] if tot else ""),
        all_antibody_pooled_not_detected=(tot["not_detected"] if tot else ""),
        ihc_reliability=summ[g]["ihc_reliability"] or "no IHC tissue data",
        note=("no antibody-based tissue IHC in HPA v25.1" if not abs_ else ""),
        hpa_cancer_summary=summ[g]["cancer_summary"]))
wr("hpa_breast_cancer_ihc_nodes.csv", rows)

# ------------------------------------------- B/C) normal breast by cell type
keys = [("Breast","Glandular cells"), ("Breast","Myoepithelial cells"), ("Breast","Adipocytes"),
        ("Soft tissue 1","Fibroblasts"), ("Soft tissue 2","Fibroblasts"), ("Adipose tissue","Adipocytes")]
lvl = collections.defaultdict(dict)
for r in norm: lvl[r["gene"]][(r["tissue"], r["cell_type"])] = r["level"]
rows = []
for g in ALL:
    d = lvl.get(g, {})
    rows.append(dict(gene=g, role=ROLE[g],
        breast_glandular=d.get(keys[0], "no data"),
        breast_myoepithelial=d.get(keys[1], "no data"),
        breast_adipocytes=d.get(keys[2], "no data"),
        soft_tissue1_fibroblasts=d.get(keys[3], "no data"),
        soft_tissue2_fibroblasts=d.get(keys[4], "no data"),
        adipose_adipocytes=d.get(keys[5], "no data"),
        ihc_reliability=summ[g]["ihc_reliability"] or "no IHC tissue data",
        hpa_normal_tissue_summary=summ[g]["normal_tissue_summary"]))
wr("hpa_normal_breast_by_celltype.csv", rows)

# --------------------------------------- B/D) HPA pathology-atlas prognostics
sv = collections.defaultdict(dict)
for r in surv:
    if r["cancer"].startswith("Breast Invasive Carcinoma"):
        sv[r["gene"]]["validation" if "validation" in r["cancer"] else "TCGA"] = r
rows = []
for g in ALL:
    t = sv[g].get("TCGA"); v = sv[g].get("validation")
    def fmt(x, k): return x[k] if x else ""
    rows.append(dict(gene=g, role=ROLE[g],
        hpa_TCGA_call=fmt(t,"prognostic"), hpa_TCGA_direction=fmt(t,"prognostic_type"),
        hpa_TCGA_p=fmt(t,"p_value"),
        hpa_validation_call=fmt(v,"prognostic"), hpa_validation_direction=fmt(v,"prognostic_type"),
        hpa_validation_p=fmt(v,"p_value")))
wr("hpa_prognostic_breast.csv", rows)

# -------------------------------------------------- C) compartment, scRNA
scv = collections.defaultdict(dict)
for r in sc: scv[r["gene"]][r["cell_type"]] = float(r["value"])
epi = ["Breast glandular cells","Breast hormone-responsive cells","Breast secretory cells",
       "Breast myoepithelial cells","Breast lactating cells"]
rows = []
for g in ALL:
    d = scv.get(g, {})
    fib = d.get("Fibroblasts")
    epis = {k: d[k] for k in epi if k in d}
    mx = max(epis.values()) if epis else None
    rows.append(dict(gene=g, role=ROLE[g],
        fibroblasts_nCPM=fib,
        **{k.replace(" ","_")+"_nCPM": epis.get(k) for k in epi},
        max_breast_epithelial_nCPM=mx,
        fibroblast_over_max_breast_epithelial=(round(fib/mx,1) if fib and mx and mx>0 else ""),
        sc_specificity_category=scspec[g]["sc_specificity_category"],
        sc_enhanced_cell_types=scspec[g]["enhanced_cell_types"],
        sc_expression_cluster=scspec[g]["sc_expression_cluster"]))
wr("hpa_compartment_singlecell.csv", rows)

# -------------------------------------------------- E) reliability / subcellular
sub = {r["gene"]: r for r in subl}
ab = collections.defaultdict(list)
for r in abrel: ab[r["gene"]].append(r)
rows = []
for g in ALL:
    a = ab.get(g, [])
    ihc = summ[g]["ihc_reliability"]
    rows.append(dict(gene=g, role=ROLE[g],
        ihc_reliability=ihc or "no IHC tissue data",
        ihc_reliability_note=summ[g]["ihc_reliability_note"],
        if_reliability=sub[g]["if_reliability"] or "no ICC/IF data",
        subcellular_main=sub[g]["main_location"], subcellular_additional=sub[g]["additional_location"],
        n_antibodies=len(a), antibodies=";".join(x["antibody"] for x in a),
        antibody_ihc_validation=";".join(f'{x["antibody"]}:{x["ihc_validation"] or "-"}' for x in a),
        antibody_western_blot=";".join(f'{x["antibody"]}:{x["western_blot"] or "-"}' for x in a),
        FLAG_gene_level_uncertain=("YES - do not over-read staining" if (ihc or "").lower()=="uncertain" else "no"),
        FLAG_some_antibody_uncertain=("yes (at least one antibody scored uncertain on IHC or WB; "
                                      "gene-level IHC reliability is still " + (ihc or "n/a") + ")"
                                      if any((x["ihc_validation"] or "").lower()=="uncertain" for x in a)
                                      or any((x["western_blot"] or "").lower()=="uncertain" for x in a) else "no"),
        FLAG_no_ihc=("YES - no antibody-based tissue IHC in HPA v25.1" if not ihc else "no")))
wr("hpa_antibody_reliability_flags.csv", rows)
