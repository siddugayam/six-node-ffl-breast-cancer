#!/usr/bin/env python3
"""43_seqreg_ext_gtex_singletissue.py

EXTENSION of 27_seqreg_gtex.py.

27_seqreg_gtex.py tested GTEx v8 independent cis-eQTLs, which are the conditionally
independent lead signals only; a gene can have many significant eQTLs and no independent
signal. This script therefore uses the
full significant single-tissue cis-eQTL set (association/singleTissueEqtl) across all 54
GTEx v8 tissues, which answers the question "is there an eQTL for COL1A1 in
breast or fibroblast?".

The design is self-controlling: the same query is run in all 54 tissues, so tissues WITH a
signal act as the positive control for the query itself, and a zero in breast/fibroblast can
be reported as a real negative rather than a failed request.

Also run: the proposed TFs (NFKB1, RELA, SP1, ETS1), so the report can say whether GTEx
even offers genetic variation in these factors in the relevant tissues, and the miRNA genes.
"""
import time, sys
import numpy as np, pandas as pd, requests

ROOT = "/path/to/revision"
RES = f"{ROOT}/results/v3"
API = "https://gtexportal.org/api/v2"
sess = requests.Session()

GENES = {
    "COL1A1": "ENSG00000108821.13", "COL3A1": "ENSG00000168542.14",
    "NFKB1": "ENSG00000109320.11", "RELA": "ENSG00000173039.18",
    "SP1": "ENSG00000185591.9", "ETS1": "ENSG00000134954.14",
    "MIR29A": "ENSG00000284032.1", "MIR29B1": "ENSG00000283797.1",
    "MIR29C": "ENSG00000284214.1", "LINC-PINT": "ENSG00000231721.6",
}
TSS = {"COL1A1": ("chr17", 50201632, "-"), "COL3A1": ("chr2", 188974320, "+")}
FOCUS_TISSUES = ["Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts"]

def get(path, **params):
    for attempt in range(4):
        try:
            r = sess.get(f"{API}/{path}", params=params, timeout=120)
            if r.status_code == 200:
                return r.json()
            time.sleep(3)
        except Exception:
            time.sleep(3)
    return None

tis = get("dataset/tissueSiteDetail", page=0, itemsPerPage=100)
TISSUES = [x["tissueSiteDetailId"] for x in tis["data"]]
print(f"{len(TISSUES)} GTEx v8 tissues", flush=True)

count_rows, var_rows = [], []
for sym, gid in GENES.items():
    for t in TISSUES:
        j = get("association/singleTissueEqtl", gencodeId=gid, tissueSiteDetailId=t,
                datasetId="gtex_v8", page=0, itemsPerPage=250)
        if j is None:
            count_rows.append(dict(gene=sym, tissue=t, n_significant_eqtl=np.nan,
                                   status="REQUEST_FAILED"))
            print(f"  FAILED {sym} {t}", flush=True)
            continue
        n = j["paging_info"]["totalNumberOfItems"]
        count_rows.append(dict(gene=sym, tissue=t, n_significant_eqtl=int(n), status="OK"))
        for v in j["data"]:
            d = None
            if sym in TSS:
                chrom, tss, strand = TSS[sym]
                p = v.get("pos")
                if p is not None:
                    d = (p - tss) if strand == "+" else (tss - p)
            var_rows.append(dict(gene=sym, tissue=t, variantId=v.get("variantId"),
                                 snpId=v.get("snpId"), chromosome=v.get("chromosome"),
                                 pos=v.get("pos"), pValue=v.get("pValue"),
                                 nes=v.get("nes"), maf=v.get("maf"),
                                 dist_to_TSS=d))
    done = sum(1 for r in count_rows if r["gene"] == sym and r["status"] == "OK")
    tot = int(np.nansum([r["n_significant_eqtl"] for r in count_rows if r["gene"] == sym]))
    print(f"{sym}: {done}/{len(TISSUES)} tissues queried OK, {tot} significant eQTLs total",
          flush=True)

counts = pd.DataFrame(count_rows)
variants = pd.DataFrame(var_rows)
counts.to_csv(f"{RES}/seqreg_ext_gtex_singletissue_counts.csv", index=False)
variants.to_csv(f"{RES}/seqreg_ext_gtex_singletissue_variants.csv", index=False)

# headline table: is the breast / fibroblast zero a real zero?
summ = []
for sym in GENES:
    c = counts[(counts.gene == sym) & (counts.status == "OK")]
    pos = c[c.n_significant_eqtl > 0]
    for t in FOCUS_TISSUES:
        row = c[c.tissue == t]
        n = int(row.n_significant_eqtl.iloc[0]) if len(row) else np.nan
        summ.append(dict(gene=sym, tissue=t, n_significant_eqtl=n,
                         n_tissues_queried=len(c),
                         n_tissues_with_any_eqtl=len(pos),
                         max_n_in_any_tissue=int(c.n_significant_eqtl.max()),
                         tissue_with_most=c.loc[c.n_significant_eqtl.idxmax(), "tissue"]
                         if len(c) else "",
                         interpretation=("REAL NEGATIVE (query works elsewhere)"
                                         if n == 0 and len(pos) > 0 else
                                         ("eQTL present" if n and n > 0 else
                                          "no eQTL in ANY tissue - uninformative"))))
summ = pd.DataFrame(summ)
summ.to_csv(f"{RES}/seqreg_ext_gtex_focus_summary.csv", index=False)
print("\n=== breast / cultured-fibroblast eQTL status ===")
print(summ.to_string(index=False))

# how close does any collagen eQTL get to the promoter?
if len(variants):
    for sym in TSS:
        v = variants[(variants.gene == sym) & variants.dist_to_TSS.notna()]
        if len(v):
            a = v.dist_to_TSS.abs()
            print(f"\n{sym}: {len(v)} significant eQTL records across tissues; "
                  f"nearest to TSS {int(a.min())} bp; within 2 kb {int((a<2000).sum())}; "
                  f"within 10 kb {int((a<10000).sum())}")
print("\nDONE")
