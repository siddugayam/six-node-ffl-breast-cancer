#!/usr/bin/env python3
"""
Reverse query: which TFs are actually bound at the COL1A1 and COL3A1 promoters?

Source: ReMap 2022 (https://remap.univ-amu.fr), all-peaks catalogue, hg38, queried through the
UCSC bigBed mirror https://hgdownload.soe.ucsc.edu/gbdb/hg38/reMap/reMap2022.bb with
bigBedToBed range queries (no full download; the file is 4.76 GB).
BED fields: chrom,start,end,name(GSE.TF.biotype),score,strand,thickStart(=peak summit),
thickEnd,itemRgb,TF,biotype.

TSS (hg38, NCBI RefSeq curated):
  COL1A1 NM_000088.4  chr17:50,184,100-50,201,631  (-)  -> TSS 50,201,631
  COL3A1 NM_000090.4  chr2:188,974,372-189,012,746 (+)  -> TSS 188,974,372
"""
import os, sys, csv
from collections import defaultdict

REV = "/path/to/revision"
SCR = "/path/to/scratch"
OUT = os.path.join(REV, "results", "v2")
TSS = {"COL1A1": ("chr17", 50201631), "COL3A1": ("chr2", 188974372)}
BREAST = ("MCF-7", "MCF7", "MDA-MB", "SUM_", "SUM1", "T-47D", "T47D", "ZR-75", "BT-", "HMEC",
          "breast", "MCF_10A", "MCF-10A", "SK-BR", "CAMA", "HCC1", "MCF")

focus = ["NFKB1", "RELA", "SP1", "ETS1"]
rows, summary = [], []
for gene, (chrom, tss) in TSS.items():
    bed = os.path.join(SCR, "remap", f"{gene}_10kb.bed")
    peaks = []
    with open(bed) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            peaks.append(dict(start=int(p[1]), end=int(p[2]), name=p[3],
                              summit=int(p[6]), tf=p[9], biotype=p[10]))
    for win in (1000, 5000):
        lo, hi = tss - win, tss + win
        # summit-in-window (strict) and peak-overlap (lenient)
        d_sum = defaultdict(lambda: dict(ds=set(), bt=set()))
        d_ovl = defaultdict(lambda: dict(ds=set(), bt=set()))
        for pk in peaks:
            if lo <= pk["summit"] <= hi:
                d_sum[pk["tf"]]["ds"].add(pk["name"]); d_sum[pk["tf"]]["bt"].add(pk["biotype"])
            if pk["start"] < hi and pk["end"] > lo:
                d_ovl[pk["tf"]]["ds"].add(pk["name"]); d_ovl[pk["tf"]]["bt"].add(pk["biotype"])
        for tf, v in sorted(d_ovl.items(), key=lambda kv: -len(kv[1]["ds"])):
            s = d_sum.get(tf, dict(ds=set(), bt=set()))
            bts = sorted(v["bt"])
            br = [b for b in bts if any(k.lower() in b.lower() for k in BREAST)]
            rows.append(dict(gene=gene, resource="ReMap2022_hg38", window_bp=win, TF=tf,
                             n_datasets_summit_in_window=len(s["ds"]),
                             n_biotypes_summit_in_window=len(s["bt"]),
                             n_datasets_peak_overlap=len(v["ds"]),
                             n_biotypes_peak_overlap=len(v["bt"]),
                             reproducible=str(len(s["ds"]) >= 2 and len(s["bt"]) >= 2).upper(),
                             breast_biotypes=";".join(br[:8]),
                             biotypes=";".join(bts[:15])))
        summary.append(dict(gene=gene, window_bp=win,
                            n_peaks_in_window_overlap=sum(1 for pk in peaks if pk["start"] < hi and pk["end"] > lo),
                            n_peaks_summit_in_window=sum(1 for pk in peaks if lo <= pk["summit"] <= hi),
                            n_TFs_overlap=len(d_ovl), n_TFs_summit=len(d_sum),
                            n_TFs_summit_reproducible=sum(1 for tf, v in d_sum.items()
                                                          if len(v["ds"]) >= 2 and len(v["bt"]) >= 2)))
        if win == 1000:
            print(f"\n### {gene} promoter TSS+/-{win} bp (ReMap 2022, hg38)")
            print(f"    peaks with summit in window: {summary[-1]['n_peaks_summit_in_window']}; "
                  f"distinct TFs: {summary[-1]['n_TFs_summit']}; "
                  f"reproducible (>=2 datasets AND >=2 biotypes): {summary[-1]['n_TFs_summit_reproducible']}")
            top = sorted(d_sum.items(), key=lambda kv: -len(kv[1]["ds"]))[:30]
            for tf, v in top:
                print(f"      {tf:12s} datasets={len(v['ds']):3d} biotypes={len(v['bt']):3d}  "
                      f"{';'.join(sorted(v['bt'])[:6])}")
            print("    FOCUS TFs:")
            for tf in focus:
                v = d_sum.get(tf); o = d_ovl.get(tf)
                print(f"      {tf:8s} summit-in-window datasets={len(v['ds']) if v else 0}"
                      f"  peak-overlap datasets={len(o['ds']) if o else 0}"
                      + (f"  biotypes={';'.join(sorted(v['bt'])[:6])}" if v else ""))

with open(os.path.join(OUT, "chip_reverse_collagen_promoters.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT, "chip_reverse_promoter_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
print("\nwrote chip_reverse_collagen_promoters.csv rows:", len(rows))
for s in summary: print(s)
