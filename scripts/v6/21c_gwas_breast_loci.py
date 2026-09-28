#!/usr/bin/env python3
"""Do the 30 prioritised nodes sit at genome-wide significant breast cancer risk loci?
GWAS Catalog v1.0.2 ontology-annotated full association file (FTP release 'latest')."""
import sys, math
import numpy as np, pandas as pd

DATA = "/path/to/revision/data/gwas"
RES  = "/path/to/revision/results/v6"
GWS  = 5e-8

assoc = pd.read_csv(f"{DATA}/gwas-catalog-download-associations-alt-full.tsv",
                    sep="\t", dtype=str, low_memory=False)
print("catalog associations:", len(assoc))

# --- build-check on a canonical breast-cancer SNP -------------------------
chk = assoc[assoc["SNPS"] == "rs2981582"]
if len(chk):
    print("build check rs2981582 (FGFR2): CHR_ID/CHR_POS seen =",
          sorted(set(zip(chk["CHR_ID"].fillna("NA"), chk["CHR_POS"].fillna("NA"))))[:4])

terms = pd.read_csv(f"{RES}/../../data/gwas/efo_breast_carcinoma_terms.csv")
desc  = set(terms["short_form"].dropna())
strict = desc | {"EFO_0000305"}                       # EFO_0000305 is obsolete -> MONDO_0004989
inclusive = strict | {"MONDO_0007254", "MONDO_0021100"}   # 'breast cancer', 'breast neoplasm'

def uri_terms(u):
    if not isinstance(u, str): return []
    return [x.strip().rsplit("/", 1)[-1] for x in u.split(",") if x.strip()]

assoc["_terms"] = assoc["MAPPED_TRAIT_URI"].map(uri_terms)
assoc["hit_strict"]    = assoc["_terms"].map(lambda t: bool(set(t) & strict))
assoc["hit_inclusive"] = assoc["_terms"].map(lambda t: bool(set(t) & inclusive))

def pnum(x):
    try: return float(x)
    except Exception: return np.nan
assoc["P"] = assoc["P-VALUE"].map(pnum)

br = assoc[assoc["hit_inclusive"]].copy()
print(f"breast-cancer associations: strict={int(assoc['hit_strict'].sum())} "
      f"inclusive={len(br)}")
brg = br[br["P"] <= GWS].copy()
print("genome-wide significant (P<=5e-8):", len(brg),
      "| unique SNPs:", brg["SNPS"].nunique(),
      "| unique studies:", brg["STUDY ACCESSION"].nunique(),
      "| unique PMIDs:", brg["PUBMEDID"].nunique())

# --- positions -----------------------------------------------------------
def expand(df):
    rows = []
    for _, r in df.iterrows():
        cids = str(r["CHR_ID"]).split(";") if pd.notna(r["CHR_ID"]) else []
        cps  = str(r["CHR_POS"]).split(";") if pd.notna(r["CHR_POS"]) else []
        snps = str(r["SNPS"]).split("; ") if pd.notna(r["SNPS"]) else [""]
        n = max(len(cids), len(cps))
        for i in range(n):
            try: pos = int(float(cps[i]))
            except Exception: continue
            ch = cids[i].strip()
            if not ch: continue
            rows.append(dict(chrom="chr"+ch, pos=pos,
                             snp=(snps[i] if i < len(snps) else snps[0]).strip(),
                             p=r["P"], trait=r["DISEASE/TRAIT"],
                             mapped_trait=r["MAPPED_TRAIT"],
                             reported_gene=r["REPORTED GENE(S)"],
                             mapped_gene=r["MAPPED_GENE"],
                             risk_allele=r["STRONGEST SNP-RISK ALLELE"],
                             or_beta=r["OR or BETA"], ci=r["95% CI (TEXT)"],
                             pmid=r["PUBMEDID"], first_author=r["FIRST AUTHOR"],
                             date=r["DATE"], study=r["STUDY ACCESSION"],
                             initial_n=r["INITIAL SAMPLE SIZE"],
                             strict=r["hit_strict"]))
    return pd.DataFrame(rows)

pos = expand(brg)
pos = pos[pos["pos"] > 0]
print("positioned GWS breast associations:", len(pos), "| unique loci(SNP):", pos["snp"].nunique())
pos.to_csv(f"{RES}/gwas_breast_gws_associations.csv", index=False)

# --- node proximity ------------------------------------------------------
nodes = pd.read_csv(f"{RES}/gwas_node_coordinates_grch38.csv")
WINDOWS = [0, 50_000, 250_000, 500_000, 1_000_000]
recs = []
for _, nd in nodes.iterrows():
    sub = pos[pos["chrom"] == nd["chrom"]]
    if not len(sub): continue
    d = np.where(sub["pos"] < nd["start"], nd["start"] - sub["pos"],
        np.where(sub["pos"] > nd["end"], sub["pos"] - nd["end"], 0))
    sub = sub.assign(dist_bp=d, dist_tss=np.abs(sub["pos"] - nd["tss"]))
    keep = sub[sub["dist_bp"] <= max(WINDOWS)]
    for _, s in keep.iterrows():
        recs.append(dict(node=nd["node"], node_type=nd["node_type"], locus=nd["locus"],
                         node_chrom=nd["chrom"], node_start=nd["start"], node_end=nd["end"],
                         **{k: s[k] for k in ["snp","chrom","pos","p","trait","mapped_trait",
                                              "reported_gene","mapped_gene","risk_allele",
                                              "or_beta","ci","pmid","first_author","date",
                                              "study","initial_n","strict"]},
                         dist_bp=int(s["dist_bp"]), dist_tss=int(s["dist_tss"])))
prox = pd.DataFrame(recs)
prox = prox.sort_values(["node","dist_bp","p"])
prox.to_csv(f"{RES}/gwas_node_proximity_all.csv", index=False)

# --- reported / mapped gene text match -----------------------------------
def genehit(field, sym):
    if not isinstance(field, str): return False
    toks = [t.strip().upper() for t in field.replace(" - ", ",").replace(";", ",").replace(" x ", ",").split(",")]
    return sym.upper() in toks
gene_nodes = nodes[nodes["node_type"].isin(["TF","gene","comparator"])]["node"].unique()
grecs = []
for sym in gene_nodes:
    m = pos[pos["mapped_gene"].map(lambda x: genehit(x, sym)) |
            pos["reported_gene"].map(lambda x: genehit(x, sym))]
    for _, s in m.iterrows():
        grecs.append(dict(node=sym, match="MAPPED_or_REPORTED_GENE",
                          **{k: s[k] for k in ["snp","chrom","pos","p","trait","mapped_gene",
                                               "reported_gene","pmid","first_author","study"]}))
gm = pd.DataFrame(grecs)
gm.to_csv(f"{RES}/gwas_node_genename_matches.csv", index=False)

# --- summary table -------------------------------------------------------
srows = []
for _, nd in nodes.drop_duplicates(subset=["node"])[["node","node_type"]].iterrows():
    n = nd["node"]
    sub = prox[prox["node"] == n] if len(prox) else pd.DataFrame()
    r = dict(node=n, node_type=nd["node_type"])
    for w in WINDOWS:
        k = "in_gene" if w == 0 else f"within_{w//1000}kb"
        s2 = sub[sub["dist_bp"] <= w] if len(sub) else sub
        r[f"n_gws_snps_{k}"] = 0 if not len(s2) else s2["snp"].nunique()
    if len(sub):
        best = sub.loc[sub["p"].idxmin()]
        nearest = sub.loc[sub["dist_bp"].idxmin()]
        r.update(min_p=best["p"], min_p_snp=best["snp"], min_p_dist_bp=int(best["dist_bp"]),
                 min_p_trait=best["trait"], min_p_mapped_gene=best["mapped_gene"],
                 min_p_reported_gene=best["reported_gene"], min_p_or_beta=best["or_beta"],
                 min_p_pmid=best["pmid"], min_p_study=best["study"],
                 nearest_snp=nearest["snp"], nearest_dist_bp=int(nearest["dist_bp"]),
                 nearest_p=nearest["p"])
    else:
        r.update(min_p=np.nan, min_p_snp="", min_p_dist_bp=np.nan, min_p_trait="",
                 min_p_mapped_gene="", min_p_reported_gene="", min_p_or_beta="",
                 min_p_pmid="", min_p_study="", nearest_snp="", nearest_dist_bp=np.nan,
                 nearest_p=np.nan)
    gsub = gm[gm["node"] == n] if len(gm) else pd.DataFrame()
    r["n_gws_snps_named_to_node"] = 0 if not len(gsub) else gsub["snp"].nunique()
    r["named_snps"] = "" if not len(gsub) else ";".join(sorted(set(gsub["snp"]))[:8])
    srows.append(r)
summ = pd.DataFrame(srows)
summ.to_csv(f"{RES}/gwas_node_summary.csv", index=False)

print()
print(summ[["node","node_type","n_gws_snps_in_gene","n_gws_snps_within_50kb",
            "n_gws_snps_within_500kb","n_gws_snps_named_to_node","min_p","nearest_dist_bp"]]
      .to_string(index=False))
