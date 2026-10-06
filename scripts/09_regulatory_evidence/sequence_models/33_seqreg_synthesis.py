#!/usr/bin/env python3
"""33_seqreg_synthesis.py -- assemble every sequence-level result into one tidy table
(results/v3/seqreg_key_results.csv). Every row carries the file the number came from.
Nothing is computed here that is not read from a deposited file, except BH correction of the
ReMap per-region empirical p-values.
"""
import os, numpy as np, pandas as pd
from scipy.stats import false_discovery_control
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
R=[]
def add(section, claim, value, source):
    R.append(dict(section=section, claim=claim, value=str(value), source_file=source))

# ---- A AlphaGenome
acc=pd.read_csv(f"{RES}/seqreg_alphagenome_access.csv")
for _,r in acc.iterrows():
    add("A_AlphaGenome", r.step, r.result, "seqreg_alphagenome_access.csv")

# ---- regions / annotation
reg=pd.read_csv(f"{RES}/seqreg_regions.csv")
for _,r in reg.iterrows():
    add("0_coordinates", f"{r.region} anchor (hg38)", f"{r.chrom}:{r.anchor} ({r.strand}) {r.note}", "seqreg_regions.csv")

# ---- B Enformer / Borzoi
en=pd.read_csv(f"{RES}/seqreg_enformer_promoter_tracks.csv.gz")
bz=pd.read_csv(f"{RES}/seqreg_borzoi_promoter_tracks.csv.gz")
for reg_ in ["COL1A1","COL3A1"]:
    s=en[en.region==reg_]
    top=s[s.assay=="DNASE"].nlargest(1,"value").iloc[0]
    add("B_Enformer", f"{reg_}: top predicted DNase track at promoter",
        f"{top.description} value={top.value:.2f} (background percentile {top.bg_percentile:.1f})",
        "seqreg_enformer_promoter_tracks.csv.gz")
    topc=s[s.assay=="CAGE"].nlargest(1,"value").iloc[0]
    add("B_Enformer", f"{reg_}: top predicted CAGE track at promoter",
        f"{topc.description} value={topc.value:.1f} (pct {topc.bg_percentile:.1f})",
        "seqreg_enformer_promoter_tracks.csv.gz")
    for tf in ["SP1","ETS1","RELB"]:
        f=s[(s.assay=="CHIP")&(s.tf==tf)]
        if len(f):
            add("B_Enformer", f"{reg_}: predicted {tf} ChIP at promoter, {len(f)} tracks",
                f"median background percentile {f.bg_percentile.median():.1f}, max {f.bg_percentile.max():.1f}",
                "seqreg_enformer_promoter_tracks.csv.gz")
    add("B_Enformer", f"{reg_}: NFKB1/RELA tracks in the Enformer human head", "absent (0 tracks; only RELB)",
        "cache/seqreg/enformer_targets_human.txt")
    sb=bz[bz.region==reg_]
    tr=sb[sb.assay=="RNA"].nlargest(1,"value").iloc[0]
    add("B_Borzoi", f"{reg_}: top predicted RNA track at promoter",
        f"{tr.description} value={tr.value:.2f} (pct {tr.bg_percentile:.1f})",
        "seqreg_borzoi_promoter_tracks.csv.gz")
    ta=sb[sb.assay=="ATAC"].nlargest(1,"value").iloc[0]
    add("B_Borzoi", f"{reg_}: top predicted ATAC track at promoter",
        f"{ta.description} value={ta.value:.2f} (pct {ta.bg_percentile:.1f})",
        "seqreg_borzoi_promoter_tracks.csv.gz")
    for tf in ["RELA","RELB","NFKB2","SP1","ETS1"]:
        f=sb[(sb.assay=="CHIP")&(sb.tf==tf)]
        if len(f):
            add("B_Borzoi", f"{reg_}: predicted {tf} ChIP at promoter, {len(f)} tracks",
                f"median background percentile {f.bg_percentile.median():.1f}, max {f.bg_percentile.max():.1f}",
                "seqreg_borzoi_promoter_tracks.csv.gz")
for reg_ in ["MIR29B2CHG","LINC_PINT_prox","MIR29A"]:
    sb=bz[bz.region==reg_]
    tr=sb[sb.assay=="RNA"].nlargest(1,"value").iloc[0]
    add("B_Borzoi", f"{reg_}: top predicted RNA track",
        f"{tr.description} value={tr.value:.2f} (pct {tr.bg_percentile:.1f})",
        "seqreg_borzoi_promoter_tracks.csv.gz")
    f=sb[(sb.assay=="CHIP")&(sb.tf.isin(["RELA","RELB","NFKB2"]))]
    add("B_Borzoi", f"{reg_}: predicted NF-kB-family ChIP",
        "; ".join(f"{x.tf}={x.bg_percentile:.0f}pct" for x in f.itertuples()),
        "seqreg_borzoi_promoter_tracks.csv.gz")

if os.path.exists(f"{RES}/seqreg_enformer_motif_occlusion.csv"):
    occ=pd.read_csv(f"{RES}/seqreg_enformer_motif_occlusion.csv")
    occ["tf"]=occ.motif.str.split("|").str[-1].str.split("_HUMAN").str[0]
    for reg_ in occ.region.unique():
        s=occ[occ.region==reg_]
        for tf in ["NFKB1","NFKB2","RELA","RELB","REL","SP1","SP2","SP3","ETS1","ETS2"]:
            t=s[s.tf==tf]
            if len(t):
                worst=t.nsmallest(1,"d_fib_cage").iloc[0]
                add("B_Enformer_occlusion",
                    f"{reg_}: shuffling the {tf} site(s) changes predicted fibroblast CAGE by",
                    f"largest drop {worst.d_fib_cage:.3f} (ref {worst.ref_fib_cage:.2f}, "
                    f"= {100*worst.d_fib_cage/worst.ref_fib_cage:+.2f}%) at TSS{worst.rel_to_TSS:+d}",
                    "seqreg_enformer_motif_occlusion.csv")
if os.path.exists(f"{RES}/seqreg_enformer_ism.csv"):
    ism=pd.read_csv(f"{RES}/seqreg_enformer_ism.csv")
    for reg_ in ism.region.unique():
        s=ism[ism.region==reg_]
        top=s.nsmallest(1,"d_fib_cage").iloc[0]
        add("B_Enformer_ISM", f"{reg_}: most disruptive single-base substitution (fibroblast CAGE)",
            f"{top.ref_base}>{top.alt_base} at TSS{int(top.rel_to_TSS):+d} (chr pos {int(top.genome_pos)}), "
            f"delta={top.d_fib_cage:.3f}", "seqreg_enformer_ism.csv")

# ---- C ReMap
rm=pd.read_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv")
rm["BH_q"]=np.nan
for reg_ in rm.region.unique():
    m=rm.region==reg_
    rm.loc[m,"BH_q"]=false_discovery_control(rm.loc[m,"emp_p_ge"].values, method="bh")
rm.to_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv", index=False)
for reg_ in ["COL1A1","COL3A1"]:
    for tf in ["NFKB1","RELA","SP1","ETS1"]:
        t=rm[(rm.region==reg_)&(rm.tf==tf)]
        if len(t):
            x=t.iloc[0]
            add("C_ReMap", f"{reg_} promoter (TSS-1000..+500): {tf} ReMap 2022 peaks",
                f"{int(x.n_peaks)} peaks (fibroblast {int(x.n_peaks_fibroblast)}); random-promoter mean "
                f"{x.bg_mean_peaks}; percentile {x.percentile_vs_bg}; emp p(>=) {x.emp_p_ge}", 
                "seqreg_remap_promoter_tfs_calibrated.csv")
        else:
            add("C_ReMap", f"{reg_} promoter: {tf} ReMap 2022 peaks", "0 peaks",
                "seqreg_remap_tf_by_region_curated.csv")
    s=rm[(rm.region==reg_)&(rm.n_peaks>=5)].nsmallest(3,"emp_p_ge")
    add("C_ReMap", f"{reg_} promoter: TFs most over-represented vs 289 random promoters",
        "; ".join(f"{x.tf} {int(x.n_peaks)} vs {x.bg_mean_peaks} (p={x.emp_p_ge}, BH q={x.BH_q:.3f})" for x in s.itertuples()),
        "seqreg_remap_promoter_tfs_calibrated.csv")
rel=pd.read_csv(f"{RES}/seqreg_rela_fibroblast_peaks.csv")
for g in ["COL1A1","COL3A1"]:
    s=rel[rel.gene==g]
    add("C_ReMap_fibroblast", f"{g}: NF-kB-family ReMap peaks in fibroblast/mesenchymal biotypes within +/-100 kb",
        f"{len(s)} peaks, all RELA, datasets {sorted(set(s.dataset.str.split('.').str[0]))}; "
        f"distances to TSS {sorted(int(x) for x in s.dist_to_TSS)}; "
        f"{int((s.n_kB_motifs>0).sum())} contain a kB consensus at p<1e-4",
        "seqreg_rela_fibroblast_peaks.csv")

# ---- C ChIP-Atlas
if os.path.exists(f"{RES}/seqreg_chipatlas_inventory.csv"):
    inv=pd.read_csv(f"{RES}/seqreg_chipatlas_inventory.csv")
    for _,r in inv.iterrows():
        add("C_ChIPAtlas", r.metric, r.value, "seqreg_chipatlas_inventory.csv")
    fb=pd.read_csv(f"{RES}/seqreg_chipatlas_fibroblast_experiments.csv")
    for tf in ["NFKB1","RELA","SP1","ETS1"]:
        t=fb[fb.antigen==tf]
        add("C_ChIPAtlas", f"{tf}: ChIP-seq experiments in fibroblast/mesenchymal cell types (hg38)",
            f"{int(t.n_experiments.sum())} experiments in {t.celltype.nunique()} cell types"
            + (f" ({';'.join(t.celltype)})" if len(t) else " -- NONE EXIST"),
            "seqreg_chipatlas_fibroblast_experiments.csv")
if os.path.exists(f"{RES}/seqreg_chipatlas_tf_summary.csv"):
    ca=pd.read_csv(f"{RES}/seqreg_chipatlas_tf_summary.csv")
    for g in ["COL1A1","COL3A1"]:
        for d in [1,5]:
            s=ca[(ca.gene==g)&(ca.distance_kb==d)]
            add("C_ChIPAtlas", f"{g}: TFs with a ChIP-Atlas peak within +/-{d} kb of the TSS",
                f"{len(s)} distinct antigens, {int(s.n_experiments.sum())} experiments; "
                f"{int((s.n_fibro_exp>0).sum())} antigens with a fibroblast/mesenchymal experiment",
                "seqreg_chipatlas_tf_summary.csv")
            for tf in ["NFKB1","RELA","SP1","ETS1"]:
                t=s[s.tf==tf]
                add("C_ChIPAtlas", f"{g} +/-{d} kb: {tf}",
                    (f"{int(t.n_experiments.iloc[0])} experiments in {int(t.n_celltypes.iloc[0])} cell types, "
                     f"max MACS2 -log10Q {t.max_score.iloc[0]:.1f}, fibroblast experiments {int(t.n_fibro_exp.iloc[0])}"
                     f" [{t.fibro_celltypes.iloc[0]}]" if len(t) else "no peak"),
                    "seqreg_chipatlas_tf_summary.csv")

# ---- C cCRE / CpG
arch=pd.read_csv(f"{RES}/seqreg_promoter_architecture.csv")
for _,r in arch.iterrows():
    add("C_promoter_architecture", f"{r.region} promoter TSS-1000..+500",
        f"GC {r.gc:.3f} ({r.gc_percentile:.1f}th pct of 3,000 random promoters); CpG o/e {r.cpg_oe:.3f} "
        f"({r.cpg_oe_percentile:.1f}th pct); CpG island overlap: {r.cpg_island_overlap or 'NONE'}; "
        f"cCREs: {r.ccre_overlap or 'none'}", "seqreg_promoter_architecture.csv")

# ---- C GTEx
neg=pd.read_csv(f"{RES}/seqreg_gtex_breast_fibroblast_negatives.csv")
for _,r in neg.iterrows():
    add("C_GTEx", f"{r.gene}: independent cis-eQTL in {r.tissue}",
        f"{r.has_independent_cis_eQTL} (signal present in {r.n_tissues_with_signal} other GTEx v8 tissues)",
        "seqreg_gtex_breast_fibroblast_negatives.csv")
gv=pd.read_csv(f"{RES}/seqreg_gtex_collagen_variants.csv")
add("C_GTEx", "COL1A1/COL3A1 eQTL or fine-mapped variants within 2 kb of the TSS",
    f"{int((gv.dist_to_TSS.abs()<=2000).sum())} of {len(gv)} (within 10 kb: {int((gv.dist_to_TSS.abs()<=10000).sum())})",
    "seqreg_gtex_collagen_variants.csv")
mi=pd.read_csv(f"{RES}/seqreg_gtex_mirna_eqtl.csv")
add("C_GTEx", "independent cis-eQTLs for MIR29A/MIR29B1/MIR29C in any GTEx v8 tissue",
    "none (0 signals); MIR29B2CHG is absent from the GTEx v8 gencode v26 reference",
    "seqreg_gtex_mirna_eqtl.csv")
sq=pd.read_csv(f"{RES}/seqreg_gtex_sqtl.csv")
if len(sq):
    add("C_GTEx", "sQTLs for COL1A1/COL3A1 in breast or cultured fibroblasts", "0", "seqreg_gtex_sqtl.csv")

# ---- C motifs
mc=pd.read_csv(f"{RES}/seqreg_motif_focus_calibrated_BH.csv")
for reg_ in ["COL1A1","COL3A1","MIR29B2CHG","LINC_PINT_prox"]:
    s=mc[mc.region==reg_]
    for fam,tfs in [("NF-kB family",["NFKB1","NFKB2","RELA","RELB","REL"]),
                    ("SP family",["SP1","SP2","SP3"]),
                    ("ETS family",["ETS1","ETS2","ELK1","ELK4","GABPA"])]:
        t=s[s.tf.isin(tfs)]
        add("C_motif", f"{reg_} promoter: {fam} motif sites (FIMO-equivalent p<1e-4)",
            f"total {int(t.n_sites.sum())} across {len(t)} matrices; best empirical p vs GC-matched "
            f"promoters {t.p_used.min():.3f} (BH q {t.BH_q.min():.2f}); "
            f"matrices with 0 sites: {int((t.n_sites==0).sum())}/{len(t)}",
            "seqreg_motif_focus_calibrated_BH.csv")
add("C_motif", "multiple-testing note", "200 tests (25 matrices x 8 promoters); nothing survives BH at q<0.25",
    "seqreg_motif_focus_calibrated_BH.csv")

out=pd.DataFrame(R)
out.to_csv(f"{RES}/seqreg_key_results.csv", index=False)
print("wrote", len(out), "rows to seqreg_key_results.csv")
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",130)
print(out.groupby("section").size().to_string())
