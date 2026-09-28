#!/usr/bin/env python3
"""34_seqreg_network_tf_edges.py -- put every TF->collagen edge that the manuscript's own
network asserts to a sequence-level test.

data/canonical_edges.tsv contains 8 TF_target edges into COL1A1 (ETS1, MKL1, MYB, NFKB1,
RELA, SP1, STAT6, TFAP2A) and NONE into COL3A1. For each TF this script reports:
  - ReMap 2022 peaks at the promoter (TSS-1000..+500) and across the +/-100 kb locus, with
    the number in fibroblast/mesenchymal biotypes, calibrated against 289 random promoters;
  - whether ChIP-Atlas hg38 contains any ChIP-seq for that factor in a fibroblast/mesenchymal
    cell type at all (so that absence of evidence can be distinguished from evidence of absence);
  - the number of JASPAR/HOCOMOCO motif matches at p<1e-4 in the promoter;
  - the Enformer and Borzoi predicted ChIP signal percentile at the promoter, where the model
    has a track for that factor.
"""
import numpy as np, pandas as pd
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"; CACHE=f"{ROOT}/cache/seqreg"
edges=pd.read_csv(f"{ROOT}/data/canonical_edges.tsv", sep="\t")
tf_edges=edges[(edges.target.isin(["COL1A1","COL3A1"]))&(edges.edge_type=="TF_target")]
print("TF->collagen edges in the network:"); print(tf_edges[["source","target","sign"]].to_string(index=False))
TFS=sorted(tf_edges.source.unique())
ALIAS={"MKL1":["MKL1","MRTFA"]}

rm=pd.read_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv")
cur=pd.read_csv(f"{RES}/seqreg_remap_tf_by_region_curated.csv")
fb=pd.read_csv(f"{RES}/seqreg_chipatlas_fibroblast_experiments.csv")
sites=pd.read_csv(f"{RES}/seqreg_motif_sites.csv")
sites["tfname"]=sites.motif.str.split("|").str[-1].str.split("_HUMAN").str[0]
en=pd.read_csv(f"{RES}/seqreg_enformer_promoter_tracks.csv.gz")
bz=pd.read_csv(f"{RES}/seqreg_borzoi_promoter_tracks.csv.gz")

rows=[]
for g in ["COL1A1","COL3A1"]:
    for tf in TFS:
        names=ALIAS.get(tf,[tf])
        p=rm[(rm.region==g)&(rm.tf.isin(names))]
        loc=cur[(cur.region==g)&(cur.window=="locus_100kb")&(cur.tf.isin(names))]
        f=fb[fb.antigen.isin(names)]
        m=sites[(sites.region==g)&(sites.tfname.isin(names))]
        e=en[(en.region==g)&(en.assay=="CHIP")&(en.tf.isin(names))]
        b=bz[(bz.region==g)&(bz.assay=="CHIP")&(bz.tf.isin(names))]
        rows.append(dict(
            gene=g, tf=tf, in_network_edge=bool(((tf_edges.source==tf)&(tf_edges.target==g)).any()),
            remap_promoter_peaks=int(p.n_peaks.iloc[0]) if len(p) else 0,
            remap_promoter_fibroblast_peaks=int(p.n_peaks_fibroblast.iloc[0]) if len(p) else 0,
            remap_bg_mean_peaks=float(p.bg_mean_peaks.iloc[0]) if len(p) else np.nan,
            remap_percentile=float(p.percentile_vs_bg.iloc[0]) if len(p) else np.nan,
            remap_emp_p=float(p.emp_p_ge.iloc[0]) if len(p) else np.nan,
            remap_locus100kb_peaks=int(loc.n_peaks.iloc[0]) if len(loc) else 0,
            remap_locus100kb_fibroblast_peaks=int(loc.n_peaks_fibroblast.iloc[0]) if len(loc) else 0,
            chipatlas_fibroblast_experiments=int(f.n_experiments.sum()),
            chipatlas_fibroblast_celltypes=";".join(sorted(f.celltype.unique())),
            motif_sites_promoter=len(m),
            motif_best_p=float(m.pvalue.min()) if len(m) else np.nan,
            enformer_tracks=len(e),
            enformer_median_pct=round(float(e.bg_percentile.median()),1) if len(e) else np.nan,
            enformer_max_pct=round(float(e.bg_percentile.max()),1) if len(e) else np.nan,
            borzoi_tracks=len(b),
            borzoi_median_pct=round(float(b.bg_percentile.median()),1) if len(b) else np.nan,
            borzoi_max_pct=round(float(b.bg_percentile.max()),1) if len(b) else np.nan))
D=pd.DataFrame(rows)
D.to_csv(f"{RES}/seqreg_network_tf_collagen_edges.csv", index=False)
pd.set_option("display.width",280)
print("\n", D.to_string(index=False))
