#!/usr/bin/env python3
"""Parse HPA v25.1 per-gene XML entries for the 26 genes of interest.
Emits: normal tissue IHC (all tissues/cell types), cancer IHC per-antibody staining
counts + per-patient rows, pathology-atlas survival calls, subcellular location,
antibody reliability."""
import os, csv, glob
import xml.etree.ElementTree as ET

XD  = "/path/to/revision/data/hpa/xml"
OUT = "/path/to/revision/results/v6"
os.makedirs(OUT, exist_ok=True)

normal_rows, cancer_count_rows, cancer_patient_rows = [], [], []
surv_rows, subloc_rows, ab_rows, summary_rows = [], [], [], []

for path in sorted(glob.glob(os.path.join(XD, "*.xml"))):
    tree = ET.parse(path); root = tree.getroot()
    entry = root.find("entry")
    gene = entry.findtext("name")
    ver  = entry.get("version")
    ident = entry.find("identifier"); ens = ident.get("id") if ident is not None else ""

    # ---- top-level normal tissue IHC (consensus across antibodies) ----
    ihc_rel, ihc_rel_desc, ihc_summary = "", "", ""
    for te in entry.findall("tissueExpression"):
        if te.get("technology") == "IHC" and te.get("assayType") == "tissue":
            v = te.find("verification")
            if v is not None and v.get("type") == "reliability":
                ihc_rel = (v.text or "").strip(); ihc_rel_desc = v.get("description") or ""
            s = te.find("summary")
            if s is not None: ihc_summary = (s.text or "").strip()
            for d in te.findall("data"):
                t = d.find("tissue")
                tissue = (t.text or "").strip(); organ = t.get("organ") or ""
                tl = d.find("level")
                tlevel = (tl.text or "").strip() if tl is not None else ""
                cells = d.findall("tissueCell")
                if not cells:
                    normal_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, organ=organ,
                                            tissue=tissue, cell_type="(tissue level)", level=tlevel,
                                            ihc_reliability=ihc_rel))
                for tc in cells:
                    ct = tc.findtext("cellType") or ""
                    lv = tc.find("level"); lvv = (lv.text or "").strip() if lv is not None else ""
                    normal_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, organ=organ,
                                            tissue=tissue, cell_type=ct, level=lvv,
                                            ihc_reliability=ihc_rel))

    # ---- pathology-atlas survival (RNA / TCGA) ----
    ce = entry.find("cancerExpression")
    if ce is not None:
        for d in ce.findall("data"):
            t = d.find("tissue"); cancer = (t.text or "").strip()
            for sa in d.findall("survivalAnalysis"):
                surv_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, cancer=cancer,
                                      is_prognostic=sa.get("isPrognostic"),
                                      prognostic=sa.get("prognostic"),
                                      prognostic_type=sa.get("prognosticType"),
                                      p_value=sa.get("pValue"),
                                      data_source=sa.get("dataSource")))

    # ---- subcellular location (ICC/IF) ----
    cellx = entry.find("cellExpression")
    if_rel = ""
    if cellx is not None:
        v = cellx.find("verification")
        if v is not None: if_rel = (v.text or "").strip()
        locs_main, locs_add = [], []
        for d in cellx.findall("data"):
            for L in d.findall("location"):
                (locs_main if L.get("status") == "main" else locs_add).append((L.text or "").strip())
        subloc_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver,
                                if_reliability=if_rel,
                                main_location=";".join(locs_main),
                                additional_location=";".join(locs_add),
                                summary=(cellx.findtext("summary") or "").strip()))
    else:
        subloc_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, if_reliability="",
                                main_location="", additional_location="", summary=""))

    # ---- per-antibody: cancer IHC + validation ----
    cancer_summary = ""
    for ab in entry.findall("antibody"):
        abid = ab.get("id")
        wb, pa, ihc_val = "", "", ""
        for v in ab.findall("verification"):
            pass
        # verifications sit inside westernBlot / proteinArray / tissueExpression children
        for child in ab:
            if child.tag == "westernBlot":
                v = child.find("verification")
                if v is not None: wb = (v.text or "").strip()
            if child.tag == "proteinArray":
                v = child.find("verification")
                if v is not None: pa = (v.text or "").strip()
            if child.tag == "tissueExpression" and child.get("assayType") == "tissue":
                v = child.find("verification")
                if v is not None: ihc_val = (v.text or "").strip()
        ab_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, antibody=abid,
                            release_version=ab.get("releaseVersion"), release_date=ab.get("releaseDate"),
                            ihc_validation=ihc_val, western_blot=wb, protein_array=pa,
                            gene_ihc_reliability=ihc_rel, gene_if_reliability=if_rel))
        for te in ab.findall("tissueExpression"):
            if te.get("assayType") != "cancer":
                continue
            s = te.find("summary")
            if s is not None and not cancer_summary:
                cancer_summary = (s.text or "").strip()
            for d in te.findall("data"):
                t = d.find("tissue")
                cancer = (t.text or "").strip(); organ = t.get("organ") or ""
                for tc in d.findall("tissueCell"):
                    ct = tc.findtext("cellType") or ""
                    counts = {}
                    for L in tc.findall("level"):
                        if L.get("type") == "staining":
                            counts[(L.text or "").strip().lower()] = int(L.get("count") or 0)
                    cancer_count_rows.append(dict(
                        gene=gene, ensembl=ens, hpa_version=ver, antibody=abid, organ=organ,
                        cancer=cancer, cell_type=ct,
                        high=counts.get("high", 0), medium=counts.get("medium", 0),
                        low=counts.get("low", 0), not_detected=counts.get("not detected", 0),
                        n_patients=sum(counts.values()),
                        ihc_reliability=ihc_rel))
                for p in d.findall("patient"):
                    st, inten = "", ""
                    for L in p.findall("level"):
                        if L.get("type") == "staining":   st = (L.text or "").strip()
                        if L.get("type") == "intensity":  inten = (L.text or "").strip()
                    snomeds = []
                    for smp in p.findall("sample"):
                        sp = smp.find("snomedParameters")
                        if sp is not None:
                            for sn in sp.findall("snomed"):
                                snomeds.append(sn.get("tissueDescription") or "")
                    cancer_patient_rows.append(dict(
                        gene=gene, ensembl=ens, hpa_version=ver, antibody=abid, cancer=cancer,
                        patient_id=p.findtext("patientId"), sex=p.findtext("sex"), age=p.findtext("age"),
                        staining=st, intensity=inten, quantity=p.findtext("quantity"),
                        location=p.findtext("location"),
                        snomed=";".join(sorted(set(x for x in snomeds if x)))))

    summary_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver,
                             ihc_reliability=ihc_rel, ihc_reliability_note=ihc_rel_desc,
                             if_reliability=if_rel,
                             normal_tissue_summary=ihc_summary, cancer_summary=cancer_summary))

def w(fn, rows, fields=None):
    if not rows:
        print("EMPTY", fn); return
    fields = fields or list(rows[0].keys())
    with open(os.path.join(OUT, fn), "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields); wr.writeheader(); wr.writerows(rows)
    print(fn, len(rows))

w("hpa_normal_tissue_ihc_v25.csv", normal_rows)
w("hpa_cancer_ihc_counts_v25.csv", cancer_count_rows)
w("hpa_cancer_ihc_patients_v25.csv", cancer_patient_rows)
w("hpa_pathology_survival_v25.csv", surv_rows)
w("hpa_subcellular_v25.csv", subloc_rows)
w("hpa_antibody_reliability_v25.csv", ab_rows)
w("hpa_gene_summaries_v25.csv", summary_rows)
