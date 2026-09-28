#!/usr/bin/env python3
"""HPA v25.1 per-gene XML: single-cell-type RNA (nCPM), consensus tissue RNA (nTPM),
and cell-type specificity calls, for the 26 genes."""
import os, csv, glob
import xml.etree.ElementTree as ET
XD="/path/to/revision/data/hpa/xml"
OUT="/path/to/revision/results/v6"
sc_rows, spec_rows, tis_rows = [], [], []
for path in sorted(glob.glob(os.path.join(XD,"*.xml"))):
    e=ET.parse(path).getroot().find("entry")
    gene=e.findtext("name"); ver=e.get("version")
    ens=e.find("identifier").get("id")
    cte=e.find("cellTypeExpression")
    if cte is not None:
        spec=cte.find("cellTypeSpecificity")
        cat=spec.get("category") if spec is not None else ""
        cells=[c.text for c in spec.findall("cellType")] if spec is not None else []
        clus=cte.find("cellTypeExpressionCluster")
        spec_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, sc_specificity_category=cat,
                              enhanced_cell_types=";".join(cells),
                              sc_distribution=cte.findtext("cellTypeDistribution") or "",
                              sc_expression_cluster=(clus.text if clus is not None else "")))
        for s in cte.findall("singleCellTypeExpression"):
            sc_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver, cell_type=s.get("name"),
                                unit=s.get("unitRNA"), value=s.get("expRNA")))
    for rx in e.findall("rnaExpression"):
        if rx.get("assayType")=="consensusTissue":
            for d in rx.findall("data"):
                t=d.find("tissue")
                lv=d.find("level")
                tis_rows.append(dict(gene=gene, ensembl=ens, hpa_version=ver,
                                     organ=t.get("organ") or "", tissue=(t.text or "").strip(),
                                     unit=lv.get("unitRNA") if lv is not None else "",
                                     nTPM=lv.get("expRNA") if lv is not None else ""))
def w(fn,rows):
    if not rows: print("EMPTY",fn); return
    with open(os.path.join(OUT,fn),"w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    print(fn,len(rows))
w("hpa_singlecell_rna_v25.csv",sc_rows)
w("hpa_singlecell_specificity_v25.csv",spec_rows)
w("hpa_consensus_tissue_rna_v25.csv",tis_rows)
