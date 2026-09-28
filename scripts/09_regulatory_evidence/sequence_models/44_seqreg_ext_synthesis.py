#!/usr/bin/env python3
"""44_seqreg_ext_synthesis.py -- assemble the extension results into one tidy table.

Produces:
  seqreg_ext_occlusion_calibrated.csv  focus-TF occlusion effects expressed as a percentile
                                       of the promoter's OWN all-site occlusion distribution
  seqreg_ext_key_results.csv           tidy claim/value/source summary of the extension
"""
import numpy as np, pandas as pd

ROOT = "/path/to/revision"
RES = f"{ROOT}/results/v3"

FAM = {"NF-kB": ["NFKB1", "NFKB2", "RELA", "RELB", "REL"],
       "SP": ["SP1", "SP2", "SP3"],
       "ETS": ["ETS1", "ETS2", "ELK1", "ELK4", "GABPA"]}
key = []
def K(section, claim, value, src):
    key.append(dict(section=section, claim=claim, value=value, source_file=src))

# ---------------------------------------------------------------- occlusion calibration
occ = pd.read_csv(f"{RES}/seqreg_ext_occlusion_allsites.csv")
cal = []
for reg_name, g in occ.groupby("region"):
    ref = float(g.ref_fib_cage.iloc[0])
    n = len(g)
    worst = g.nsmallest(1, "d_fib_cage").iloc[0]
    K("B_occlusion_calibration", f"{reg_name}: sites occluded (all JASPAR/HOCOMOCO calls p<1e-4)",
      f"{n} sites, {int(g.n_shuffles.iloc[0])} dinucleotide shuffles each; reference "
      f"fibroblast CAGE {ref:.3f}", "seqreg_ext_occlusion_allsites.csv")
    K("B_occlusion_calibration",
      f"{reg_name}: MOST load-bearing site in the promoter (de-facto positive control)",
      f"{worst.motif} at TSS{int(worst.rel_to_TSS):+d}, "
      f"d={worst.d_fib_cage:.3f} ({worst.pct_change_fib_cage:+.2f}%)",
      "seqreg_ext_occlusion_allsites.csv")

    for fam, tfs in FAM.items():
        s = g[g.tf_name.isin(tfs)]
        if not len(s):
            cal.append(dict(region=reg_name, family=fam, n_sites=0, best_d=np.nan,
                            best_pct_change=np.nan, percentile_among_all_sites=np.nan,
                            rank=np.nan, n_sites_total=n, best_motif="", best_rel_to_TSS=np.nan))
            K("B_occlusion_calibration", f"{reg_name}: {fam} family sites at p<1e-4", "NONE",
              "seqreg_ext_occlusion_allsites.csv")
            continue
        b = s.nsmallest(1, "d_fib_cage").iloc[0]
        rank = int((g.d_fib_cage < b.d_fib_cage).sum()) + 1
        pct = 100.0 * (g.d_fib_cage < b.d_fib_cage).mean()
        cal.append(dict(region=reg_name, family=fam, n_sites=len(s), best_d=float(b.d_fib_cage),
                        best_pct_change=float(b.pct_change_fib_cage),
                        percentile_among_all_sites=pct, rank=rank, n_sites_total=n,
                        best_motif=b.motif, best_rel_to_TSS=int(b.rel_to_TSS)))
        K("B_occlusion_calibration",
          f"{reg_name}: best {fam} site, calibrated against all {n} occluded sites",
          f"{b.motif} at TSS{int(b.rel_to_TSS):+d}: d={b.d_fib_cage:.3f} "
          f"({b.pct_change_fib_cage:+.2f}%), rank {rank}/{n} "
          f"({pct:.1f}th percentile; {rank-1} sites are MORE load-bearing)",
          "seqreg_ext_occlusion_allsites.csv")
        for t in sorted(s.tf_name.unique()):
            bt = s[s.tf_name == t].nsmallest(1, "d_fib_cage").iloc[0]
            K("B_occlusion_calibration", f"{reg_name}: {t} best single site",
              f"TSS{int(bt.rel_to_TSS):+d} d={bt.d_fib_cage:.3f} ({bt.pct_change_fib_cage:+.2f}%)",
              "seqreg_ext_occlusion_allsites.csv")
    top = g.nsmallest(10, "d_fib_cage")
    K("B_occlusion_calibration", f"{reg_name}: top-10 most load-bearing motifs",
      "; ".join(f"{r.motif}@{int(r.rel_to_TSS):+d}({r.pct_change_fib_cage:+.1f}%)"
                for _, r in top.iterrows()), "seqreg_ext_occlusion_allsites.csv")
cal = pd.DataFrame(cal)
cal.to_csv(f"{RES}/seqreg_ext_occlusion_calibrated.csv", index=False)

# ---------------------------------------------------------------- distal RELA peaks
try:
    pk = pd.read_csv(f"{RES}/seqreg_ext_rela_peak_occlusion.csv")
    for reg_name, g in pk.groupby("region"):
        sig = g[(g.z_vs_ctrl < -3)]
        K("B_distal_RELA_peaks", f"{reg_name}: fibroblast NF-kB ReMap peaks occluded",
          f"{len(g)} peaks; {len(sig)} with z < -3 vs width-matched control windows "
          f"(i.e. the model treats them as functional for fibroblast CAGE at the TSS)",
          "seqreg_ext_rela_peak_occlusion.csv")
        for _, r in g.sort_values("z_vs_ctrl").iterrows():
            K("B_distal_RELA_peaks",
              f"{reg_name} peak {int(r.peak_index)} ({r.biotype}, {int(r.dist_to_TSS):+d} bp, "
              f"{int(r.n_kB_motifs) if not pd.isna(r.n_kB_motifs) else 0} kB motifs)",
              f"d={r.d_fib_cage:+.3f} ({r.pct_change_fib_cage:+.2f}%), "
              f"control {r.ctrl_mean_d:+.3f}+/-{r.ctrl_sd_d:.3f} (n={int(r.n_ctrl)}), "
              f"z={r.z_vs_ctrl:+.2f}, emp p={r.emp_p_more_negative:.3f}",
              "seqreg_ext_rela_peak_occlusion.csv")
except FileNotFoundError:
    K("B_distal_RELA_peaks", "distal peak occlusion", "NOT AVAILABLE", "")

# ---------------------------------------------------------------- splice-donor control
try:
    sd = pd.read_csv(f"{RES}/seqreg_ext_splicedonor_control.csv")
    for reg_name, g in sd.groupby("region"):
        don = g[g.variant.str.startswith("donor_kill")]
        ets = g[g.variant.str.startswith("ets_core")]
        ctl = g[g.variant == "exonic_control_4bp"]
        K("B_splice_donor_confound",
          f"{reg_name}: killing the 2-bp exon1 5' splice donor (GT->AA / GT->CC)",
          "; ".join(f"{r.description}: {r.pct_change_fib_cage:+.2f}%" for _, r in don.iterrows()),
          "seqreg_ext_splicedonor_control.csv")
        if len(ets):
            K("B_splice_donor_confound",
              f"{reg_name}: killing the 4-bp exonic GGAA 'ETS core' instead",
              "; ".join(f"{r.description}: {r.pct_change_fib_cage:+.2f}%" for _, r in ets.iterrows()),
              "seqreg_ext_splicedonor_control.csv")
        K("B_splice_donor_confound", f"{reg_name}: neutral 4-bp exonic control",
          f"{ctl.pct_change_fib_cage.iloc[0]:+.2f}%", "seqreg_ext_splicedonor_control.csv")
    K("B_splice_donor_confound", "CORRECTION to the first pass",
      "The deposited claim 'COL3A1: shuffling the ETS1 site changes predicted fibroblast CAGE "
      "by -47.60% at TSS+186' is confounded: COL3A1 exon 1 ends at TSS+195, the called ETS1 "
      "site (TSS+186..+198, CAACAGGAAGGTG) straddles the exon1/intron1 junction, and its GGAA "
      "is the last exonic bases while the rest is the canonical donor AG|GTGAGT. Destroying "
      "only the 2-bp donor reproduces -27% to -29%, while destroying only the 4-bp GGAA gives "
      "-5% to -8% (and even that mutates two bases of the donor's exonic consensus). The same "
      "geometry explains the COL1A1 TSS+214 sites (exon 1 ends at TSS+221; donor kill -9% to "
      "-10%). These are splice-donor effects, not evidence of ETS1 binding.",
      "seqreg_ext_splicedonor_control.csv")
except FileNotFoundError:
    pass

# ---------------------------------------------------------------- ChIP-Atlas fibroblast
foc = pd.read_csv(f"{RES}/seqreg_ext_chipatlas_focus_tfs.csv")
summ = pd.read_csv(f"{RES}/seqreg_ext_chipatlas_fibroblast_tf_summary.csv")
elist = pd.read_csv(f"{RES}/seqreg_ext_chipatlas_fibroblast_experiment_list.csv")
K("C_ChIPAtlas_direct", "fibroblast/mesenchymal TF ChIP-seq experiments downloaded and "
  "intersected at peak level (hg38, MACS2 q<1e-5)",
  f"{elist.srx.nunique()} experiments, {elist.celltype.nunique()} cell types, "
  f"{elist.antigen.nunique()} antigens", "seqreg_ext_chipatlas_fibroblast_experiment_list.csv")
for rg in ["COL1A1", "COL3A1"]:
    s = summ[(summ.region == rg) & (summ.window == "promoter_1kb_500")]
    K("C_ChIPAtlas_direct", f"{rg} promoter (TSS-1000..+500): TFs with a FIBROBLAST peak",
      f"{len(s)} TFs: " + "; ".join(f"{r.antigen}({r.n_experiments_with_peak}/"
                                    f"{r.n_fibroblast_experiments})" for _, r in s.iterrows()),
      "seqreg_ext_chipatlas_fibroblast_tf_summary.csv")
    for ag in ["NFKB1", "RELA", "SP1", "ETS1", "MKL1", "CTCF"]:
        r = foc[(foc.antigen == ag) & (foc.region == rg)]
        p = r[r.window == "promoter_1kb_500"]
        w5 = r[r.window == "pm5kb"]; w100 = r[r.window == "pm100kb"]
        ne = int(p.n_fibroblast_experiments.iloc[0]) if len(p) else 0
        K("C_ChIPAtlas_direct", f"{rg}: {ag} in fibroblast/mesenchymal experiments",
          (f"{ne} such experiments exist genome-wide; peaks: promoter "
           f"{int(p.n_experiments_with_peak.iloc[0])}, +/-5 kb {int(w5.n_experiments_with_peak.iloc[0])}, "
           f"+/-100 kb {int(w100.n_experiments_with_peak.iloc[0])}")
          if ne else "NO fibroblast/mesenchymal ChIP-seq experiment exists for this factor "
                     "-- untestable, absence of evidence not evidence of absence",
          "seqreg_ext_chipatlas_focus_tfs.csv")

# ---------------------------------------------------------------- GTEx
g = pd.read_csv(f"{RES}/seqreg_ext_gtex_focus_summary.csv")
for _, r in g.iterrows():
    K("C_GTEx_singletissue", f"{r.gene}: significant cis-eQTLs in {r.tissue}",
      f"{r.n_significant_eqtl} (this gene has eQTLs in {r.n_tissues_with_any_eqtl}/"
      f"{r.n_tissues_queried} tissues; max {r.max_n_in_any_tissue} in {r.tissue_with_most}) "
      f"-- {r.interpretation}", "seqreg_ext_gtex_focus_summary.csv")
v = pd.read_csv(f"{RES}/seqreg_ext_gtex_singletissue_variants.csv")
for sym in ["COL1A1", "COL3A1"]:
    s = v[(v.gene == sym) & v.dist_to_TSS.notna()]
    if len(s):
        a = s.dist_to_TSS.abs()
        K("C_GTEx_singletissue", f"{sym}: distance from significant eQTLs to the TSS",
          f"{len(s)} eQTL records across all tissues; nearest {int(a.min())} bp; "
          f"within 2 kb {int((a < 2000).sum())}; within 10 kb {int((a < 10000).sum())}",
          "seqreg_ext_gtex_singletissue_variants.csv")

# ---------------------------------------------------------------- AlphaGenome access
acc = pd.read_csv(f"{RES}/seqreg_ext_alphagenome_access.csv")
for _, r in acc.iterrows():
    K("A_AlphaGenome_access", r.step, r.result, "seqreg_ext_alphagenome_access.csv")

kdf = pd.DataFrame(key)
kdf.to_csv(f"{RES}/seqreg_ext_key_results.csv", index=False)
print(f"wrote seqreg_ext_key_results.csv with {len(kdf)} rows")
print(cal.to_string(index=False))
