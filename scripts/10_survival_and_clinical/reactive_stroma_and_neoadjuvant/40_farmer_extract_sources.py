#!/usr/bin/env python3
"""
40_farmer_extract_sources.py
Extract published stromal / ECM / wound-response gene lists from PRIMARY sources
that were downloaded to cache/v5/farmer/.

PRIMARY SOURCES (all verified against PubMed before use):
 * Farmer P, Bonnefoi H, Anderle P, Cameron D, Wirapati P, Becette V, et al.
   "A stroma-related gene signature predicts resistance to neoadjuvant chemotherapy
   in breast cancer." Nat Med. 2009;15(1):68-74. PMID 19122658. doi:10.1038/nm.1908
   -> Supplementary Table 4 (MOESM20): the 50-gene stromal metagene (DCN-anchored).
   -> Supplementary Table 5 (MOESM21): published gene signatures used in that study
      (DTF desmoid-type fibromatosis, SFT, EMT, WNT, TGF-beta, Finak E/I reference lists).
 * Triulzi T, et al. PLoS One. 2013;8(2):e56761. PMID 23441215 -> Table S1 ECM gene list.
 * Winslow S, et al. Breast Cancer Res. 2015;17:23. PMID 25848820 -> Table 2 gene signature 1.
Outputs plain one-symbol-per-line .txt files into cache/v5/farmer/lists/.
"""
import os, re, html, csv

BASE = "/path/to/revision"
CACHE = os.path.join(BASE, "cache/v5/farmer")
OUT = os.path.join(CACHE, "lists")
os.makedirs(OUT, exist_ok=True)

def rd(p):
    return [l.split("\t") for l in open(p, newline="").read().replace("\r", "\n").split("\n")]

def write(name, genes, prov):
    genes = [g for g in dict.fromkeys([g.strip() for g in genes]) if g and g.lower() != "nan"]
    with open(os.path.join(OUT, name + ".txt"), "w") as fh:
        fh.write("\n".join(genes) + "\n")
    print(f"{name:34s} n={len(genes):4d}   {prov}")
    return genes

inv = []

# ---------------------------------------------------------------- Farmer 2009 Sup Table 4
rows = rd(os.path.join(CACHE, "farmer2009_MOESM20.txt"))
assert rows[5][0].startswith("Representative gene used in the linear model"), rows[5][0]
assert rows[57][0].startswith("Number of genes from the original signature")
assert rows[57][2].strip() == "50", "column 1 must retain all 50 original genes"
farmer50 = [rows[i][2].strip() for i in range(5, 55)]
assert len(farmer50) == 50 and farmer50[0] == "DCN"
g = write("FARMER2009_STROMAL_50", farmer50,
          "Farmer 2009 Nat Med Suppl Table 4, signature column 1 (DCN-anchored)")
inv.append(("FARMER2009_STROMAL_50", len(g), "Farmer et al. Nat Med 2009", "19122658",
            "10.1038/nm.1908", "Supplementary Table 4, column 1 (top 50 genes of the "
            "DCN-anchored multilinear model = the stroma-related signature)"))

# the 49 leave-one-out variants (columns 2..50) -- used to show robustness of overlap
for j in range(3, 52):
    rep = rows[5][j].strip()
    var = [rows[i][j].strip() for i in range(5, 55)]
    with open(os.path.join(CACHE, "farmer2009_variant_%02d_%s.txt" % (j - 1, rep)), "w") as fh:
        fh.write("\n".join(var) + "\n")
print("wrote 49 Farmer leave-one-out variant metagenes (Suppl Table 4 columns 2-50)")

# ---------------------------------------------------------------- Farmer 2009 Sup Table 5
rows5 = rd(os.path.join(CACHE, "farmer2009_MOESM21.txt"))
hdr = [h.strip() for h in rows5[0]]
cols = {j: [] for j in range(len(hdr))}
for r in rows5[1:]:
    for j in range(min(len(hdr), len(r))):
        v = r[j].strip()
        if v:
            cols[j].append(v)
want = {
    0: ("FINAK_REFERENCE_EPITHELIUM_E", "Finak 2006 reference epithelial list, via Farmer 2009 Suppl Table 5 col 1"),
    1: ("FINAK_REFERENCE_STROMA_I", "Finak 2006 reference stromal list, via Farmer 2009 Suppl Table 5 col 2"),
    2: ("WEST_DTF_FIBROMATOSIS", "desmoid-type fibromatosis (DTF) top genes of West 2005 / Beck 2008, via Farmer 2009 Suppl Table 5 col 3"),
    3: ("MAMMOSPHERE_MSP", "mammosphere signature, via Farmer 2009 Suppl Table 5 col 4"),
    4: ("FARMER_EMT_SIGNATURE", "EMT gene set used by Farmer 2009, Suppl Table 5 col 5"),
    5: ("FARMER_WNT_SIGNATURE", "WNT gene set used by Farmer 2009, Suppl Table 5 col 6"),
    6: ("FARMER_TGFBETA_TARGETS", "TGF-beta target gene set used by Farmer 2009, Suppl Table 5 col 7"),
    7: ("WEST_SFT_SIGNATURE", "solitary fibrous tumour (SFT) signature of West 2005, via Farmer 2009 Suppl Table 5 col 8"),
}
for j, (nm, prov) in want.items():
    gg = write(nm, cols[j], prov)
    inv.append((nm, len(gg), "Farmer et al. Nat Med 2009 Suppl Table 5", "19122658",
                "10.1038/nm.1908", prov))

# ---------------------------------------------------------------- Triulzi 2013 ECM list
try:
    import pandas as pd, warnings
    warnings.filterwarnings("ignore")
    d = pd.read_excel(os.path.join(CACHE, "sup/triulzi2013/pone.0056761.s007.xls"),
                      sheet_name="ECMlist", header=None)
    tri = [str(x).strip() for x in d.iloc[3:, 0].dropna().tolist()]
    gg = write("TRIULZI_ECM_GENELIST", tri,
               "Triulzi 2013 PLoS One Table S1 (a priori ECM gene list used to define ECM3)")
    inv.append(("TRIULZI_ECM_GENELIST", len(gg), "Triulzi et al. PLoS One 2013", "23441215",
                "10.1371/journal.pone.0056761", "Table S1 ECM gene list"))
except Exception as e:
    print("Triulzi extraction FAILED:", e)

# ---------------------------------------------------------------- Winslow 2015 signature 1
t = open(os.path.join(CACHE, "winslow2015.xml")).read()
tw = [b for b in re.findall(r"<table-wrap\b.*?</table-wrap>", t, re.S)
      if "Stromal gene signatures of highly correlating genes" in b][0]
grid = []
for r in re.findall(r"<tr\b.*?</tr>", tw, re.S):
    grid.append([html.unescape(re.sub("<[^>]+>", "", c)).strip()
                 for c in re.findall(r"<t[hd]\b.*?</t[hd]>", r, re.S)])
# Column 1 of the first block: the collagen/CAF gene signature. Rows are ragged because
# empty cells are dropped, so we take column 0 only for as long as it stays in signature 1
# (rows 1..14; row 15 onwards the first cell belongs to a later signature).
wins1 = [grid[i][0] for i in range(1, 15)]
assert wins1[0] == "COL1A2" and wins1[-1] == "VCAN", wins1
gg = write("WINSLOW_STROMAL_SIG1", wins1,
           "Winslow 2015 Breast Cancer Res Table 2, gene signature 1 (collagen/CAF set; "
           "multivariable HR 1.79, p=0.0039 for new tumour event)")
inv.append(("WINSLOW_STROMAL_SIG1", len(gg), "Winslow et al. Breast Cancer Res 2015",
            "25848820", "10.1186/s13058-015-0530-2", "Table 2, gene signature 1"))

with open(os.path.join(BASE, "results/v5/farmer_signature_sources.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["signature", "n_genes_as_published", "source_paper", "pmid", "doi", "provenance"])
    w.writerows(inv)
print("\nWROTE results/v5/farmer_signature_sources.csv")
