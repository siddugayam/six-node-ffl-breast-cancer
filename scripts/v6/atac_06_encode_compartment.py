#!/usr/bin/env python3
"""(a, supporting) ENCODE chromatin accessibility in breast carcinoma cells vs
fibroblasts at the COL1A1/COL3A1 promoters, and the extent to which the
TCGA-BRCA tumour ATAC peak set captures fibroblast-specific elements.
"""
import json, os, subprocess, gzip, tempfile, numpy as np, pandas as pd

ROOT = "/path/to/revision"
C = f"{ROOT}/cache/v6/atac"; E = f"{C}/encode"; R = f"{ROOT}/results/v6"
man = pd.DataFrame(json.load(open(f"{E}/manifest.json")))
EPI = ["MCF-7", "MCF 10A", "breast epithelium"]
FIB = ["fibroblast of mammary gland", "fibroblast of dermis", "fibroblast of lung", "IMR-90"]
man["compartment"] = np.where(man.biosample.isin(EPI), "epithelial",
                       np.where(man.biosample.isin(FIB), "fibroblast", "other"))
man["label"] = man.biosample + " (" + man.assay.str.replace("-seq", "") + ", " + man.file + ")"
print(man[["biosample", "assay", "experiment", "file", "output_type", "n_peaks", "compartment"]].to_string(index=False))

tss = pd.read_csv(f"{C}/refseq_select_tss.tsv", sep="\t").drop_duplicates("symbol").set_index("symbol")
GENES = ["COL1A1", "COL3A1", "COL1A2", "COL5A1", "FN1", "DCN", "LUM", "POSTN", "FAP", "PDGFRB",
         "EPCAM", "KRT8", "ESR1", "GATA3", "GAPDH", "ACTB", "PTPRC"]
rows = []
for g in GENES:
    if g not in tss.index: continue
    t = tss.loc[g]
    rows.append((t.chrom, max(0, t.tss - 1000), t.tss + 1000, f"{g}_prom"))
    rows.append((t.chrom, max(0, t.tss - 100), t.tss + 100, f"{g}_TSS"))
    lo, hi = min(t.txStart, t.txEnd), max(t.txStart, t.txEnd)
    rows.append((t.chrom, lo, hi, f"{g}_body"))
q = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"]).sort_values(["chrom", "start"])
qbed = f"{C}/query_regions.bed"; q.to_csv(qbed, sep="\t", header=False, index=False)

def peaks(f):
    p = f"{E}/{f}.bed.gz"
    d = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    d = d[d.chrom.str.match(r"^chr[0-9XY]+$")].drop_duplicates()      # ENCODE beds repeat rows
    return d.sort_values(["chrom", "start"])

# ---- 1. promoter/body overlap per biosample --------------------------------
res = {}
tmpd = tempfile.mkdtemp()
for _, r in man.iterrows():
    pk = peaks(r.file); pf = f"{tmpd}/{r.file}.bed"; pk.to_csv(pf, sep="\t", header=False, index=False)
    out = subprocess.run(["bedtools", "intersect", "-a", qbed, "-b", pf, "-c"],
                         capture_output=True, text=True).stdout
    cnt = {ln.split("\t")[3]: int(ln.split("\t")[4]) for ln in out.strip().split("\n")}
    res[r.label] = cnt
O = pd.DataFrame(res)
O.index.name = "region"
O.to_csv(f"{R}/atac_encode_region_overlap.csv")
Ob = (O > 0).astype(int)
Ob.to_csv(f"{R}/atac_encode_region_accessible.csv")
for suf, lab in [("_prom", "TSS +/-1 kb"), ("_TSS", "TSS +/-100 bp")]:
    sub = Ob.loc[[i for i in Ob.index if i.endswith(suf)]]
    sub.index = [i.rsplit("_", 1)[0] for i in sub.index]
    sub.columns = [c.split(" (")[0] + "|" + c.split("(")[1].split(",")[0] for c in sub.columns]
    print(f"\n== accessible (1/0) at {lab} ==")
    print(sub.to_string())
    epi = [c for c in sub.columns if c.split("|")[0] in EPI]
    fib = [c for c in sub.columns if c.split("|")[0] in FIB]
    print("  n_epithelial_accessible / n_fibroblast_accessible")
    for g in sub.index:
        print(f"   {g:10s} {int(sub.loc[g, epi].sum())}/{len(epi)}   {int(sub.loc[g, fib].sum())}/{len(fib)}")

# ---- 2. compartment-specific ENCODE elements vs the TCGA-BRCA peak set -----
def multi(files, tag):
    fs = [f"{tmpd}/{f}.bed" for f in files]
    cat = f"{tmpd}/{tag}_cat.bed"
    with open(cat, "w") as fh:
        subprocess.run(["cat"] + fs, stdout=fh)
    srt = f"{tmpd}/{tag}_srt.bed"
    with open(srt, "w") as fh:
        subprocess.run(["bedtools", "sort", "-i", cat], stdout=fh)
    mrg = f"{tmpd}/{tag}_mrg.bed"
    with open(mrg, "w") as fh:
        subprocess.run(["bedtools", "merge", "-i", srt], stdout=fh)
    return mrg

epi_files = man[man.compartment == "epithelial"].file.tolist()
fib_files = man[man.compartment == "fibroblast"].file.tolist()
epi_u = multi(epi_files, "epi"); fib_u = multi(fib_files, "fib")

def subtract(a, b, tag):
    out = f"{tmpd}/{tag}.bed"
    with open(out, "w") as fh:
        subprocess.run(["bedtools", "intersect", "-a", a, "-b", b, "-v"], stdout=fh)
    return out
fib_only = subtract(fib_u, epi_u, "fibonly")
epi_only = subtract(epi_u, fib_u, "epionly")
shared_f = f"{tmpd}/shared.bed"
with open(shared_f, "w") as fh:
    subprocess.run(["bedtools", "intersect", "-a", fib_u, "-b", epi_u, "-u"], stdout=fh)

brca = pd.read_csv(f"{C}/brca_peaks.bed", sep="\t", header=None, names=["chrom", "start", "end", "name"])
bb = f"{tmpd}/brca_srt.bed"; brca.sort_values(["chrom", "start"]).to_csv(bb, sep="\t", header=False, index=False)

summ = []
for tag, f in [("fibroblast_specific", fib_only), ("epithelial_specific", epi_only), ("shared", shared_f)]:
    n = sum(1 for _ in open(f))
    hit = subprocess.run(["bedtools", "intersect", "-a", f, "-b", bb, "-u"], capture_output=True, text=True).stdout
    nh = len([l for l in hit.strip().split("\n") if l])
    summ.append(dict(set=tag, n_regions=n, n_in_TCGA_BRCA_peakset=nh, frac=nh / n if n else np.nan))
Sm = pd.DataFrame(summ)
Sm.to_csv(f"{R}/atac_encode_compartment_capture.csv", index=False)
print("\n== ENCODE compartment-specific elements captured by the TCGA-BRCA ATAC peak set ==")
print(Sm.to_string(index=False))
print("DONE")
