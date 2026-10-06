#!/usr/bin/env python3
"""36 -- append the remaining headline numbers to results/v3/seqreg_key_results.csv."""
import pandas as pd, numpy as np
ROOT="/path/to/revision"; RES=f"{ROOT}/results/v3"
K=pd.read_csv(f"{RES}/seqreg_key_results.csv")
R=[]
def add(section, claim, value, source): R.append(dict(section=section, claim=claim, value=str(value), source_file=source))

# ISM -> motif attribution
P=pd.read_csv(f"{RES}/seqreg_enformer_ism_by_position.csv")
for reg in P.region.unique():
    s=P[P.region==reg].nlargest(6,"max_abs_d_fib_cage")
    add("B_Enformer_ISM", f"{reg}: the six most load-bearing promoter positions and the motifs covering them",
        "; ".join(f"TSS{int(x.rel_to_TSS):+d} (max|d|={x.max_abs_d_fib_cage:.1f}, best motif {x.best_motif} "
                  f"p={x.best_motif_p:.1e})" for x in s.itertuples()),
        "seqreg_enformer_ism_by_position.csv")

# NF-kB binding context
NC=pd.read_csv(f"{RES}/seqreg_nfkb_context_summary.csv")
for _,r in NC.iterrows():
    add("C_ReMap_NFkB_context", f"{r.gene}: NF-kB-family ReMap peaks within +/-100 kb",
        f"{int(r.total_peaks)} peaks over {int(r.n_biotypes)} biotypes; {int(r.promoter_peaks)} in the promoter "
        f"(+/-1 kb), of which {int(r.promoter_fibroblast)} fibroblast and {int(r.promoter_mesenchymal)} "
        f"mesenchymal-lineage; {int(r.fibroblast_peaks)} fibroblast peaks in total, "
        f"{int(r.mesenchymal_peaks)} mesenchymal-lineage; {100*r.stimulated_fraction:.1f}% of all peaks come "
        f"from a cytokine/genotoxin-stimulated condition", "seqreg_nfkb_context_summary.csv")

# cCRE summary
C=pd.read_csv(f"{RES}/seqreg_ccre_summary.csv")
for _,r in C.iterrows():
    add("C_ENCODE_cCRE", f"{r.region}: ENCODE cCREs within +/-100 kb",
        f"{int(r.n_ccre_100kb)} total (PLS {int(r.PLS)}, pELS {int(r.pELS)}, dELS {int(r.dELS)}); "
        f"promoter-like (PLS) element within 2 kb: {r.PLS_accessions if isinstance(r.PLS_accessions,str) and r.PLS_accessions else 'NONE'}",
        "seqreg_ccre_summary.csv")

# BH-corrected ReMap
RM=pd.read_csv(f"{RES}/seqreg_remap_promoter_tfs_calibrated.csv")
for reg in ["COL1A1","COL3A1","LINC_PINT_prox","MIR29B2CHG","MIR29A","MIR29B1"]:
    s=RM[RM.region==reg]
    sig=s[s.BH_q<0.25].sort_values("BH_q")
    add("C_ReMap_BH", f"{reg} promoter: TFs over-represented vs 289 random promoters after BH (q<0.25), "
                      f"{len(s)} TFs tested",
        ("; ".join(f"{x.tf} {int(x.n_peaks)} peaks (fibroblast {int(x.n_peaks_fibroblast)}) vs bg mean "
                   f"{x.bg_mean_peaks} q={x.BH_q:.3f}" for x in sig.itertuples()) if len(sig)
         else f"NONE (minimum BH q = {s.BH_q.min():.3f})"),
        "seqreg_remap_promoter_tfs_calibrated.csv")

# network TF edges
N=pd.read_csv(f"{RES}/seqreg_network_tf_collagen_edges.csv")
for _,r in N[N.in_network_edge].iterrows():
    add("D_network_edges", f"network edge {r.tf} -> {r.gene}: sequence-level support",
        f"ReMap promoter peaks {int(r.remap_promoter_peaks)} (fibroblast {int(r.remap_promoter_fibroblast_peaks)}), "
        f"locus +/-100 kb {int(r.remap_locus100kb_peaks)} (fibroblast {int(r.remap_locus100kb_fibroblast_peaks)}); "
        f"ChIP-Atlas fibroblast/mesenchymal experiments for this factor anywhere in the genome: "
        f"{int(r.chipatlas_fibroblast_experiments)}; promoter motif sites p<1e-4: {int(r.motif_sites_promoter)}; "
        f"Enformer tracks {int(r.enformer_tracks)} (median pct {r.enformer_median_pct}); "
        f"Borzoi tracks {int(r.borzoi_tracks)} (median pct {r.borzoi_median_pct})",
        "seqreg_network_tf_collagen_edges.csv")
add("D_network_edges", "TF_target edges into COL3A1 in data/canonical_edges.tsv", "0 (all 8 TF->collagen edges are into COL1A1)",
    "data/canonical_edges.tsv")

# caveats
add("E_caveats", "ChIP-Atlas target-gene files that could not be fetched", "203 of 3,702 TF x distance requests (115 distinct antigens), including MKL1/MRTFA (HTTP 404)",
    "seqreg_chipatlas_fetch_errors.csv")
add("E_caveats", "ReMap peaks served by UCSC carry score 0", "peak scores are not available in the UCSC-hosted reMap2022.bb; only presence, TF and biotype were used",
    "cache/seqreg/remap2022_hg38_loci.bed")
add("E_caveats", "ISM allele convention", "ref_base/alt_base are plus-strand alleles at genome_pos; rel_to_TSS is in transcription orientation; all 1,000 positions verified against hg38",
    "seqreg_enformer_ism_conventions.csv")

K2=pd.concat([K, pd.DataFrame(R)], ignore_index=True)
K2.to_csv(f"{RES}/seqreg_key_results.csv", index=False)
print("total rows now", len(K2))
pd.set_option("display.width",300); pd.set_option("display.max_colwidth",200)
print(pd.DataFrame(R)[["section","claim","value"]].to_string(index=False))
