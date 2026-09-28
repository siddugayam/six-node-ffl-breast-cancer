#!/usr/bin/env python3
"""Assemble the quantitative compendium for the 30 prioritised nodes.
Every field is read from an existing verified results file; nothing is imputed."""
import numpy as np, pandas as pd, json, os

R = "results"
OUT = "results/v5"
os.makedirs(OUT, exist_ok=True)

NODES = [("TF", n) for n in ['E2F1','EZH2','GATA3','BRCA1','JUN','EGR2','ESR1','SREBF1','DNMT1','E2F3']] + \
        [("Gene", n) for n in ['CCND2','COL1A1','STAT5A','MYBL2','FN1','PDGFRB','MET','CXCL12','MMP14','PLAU']] + \
        [("miRNA", n) for n in ['hsa-miR-21','hsa-miR-195','hsa-miR-204','hsa-miR-383','hsa-miR-124',
                                'hsa-miR-155','hsa-miR-429','hsa-miR-141','hsa-miR-34a','hsa-miR-101']]
names = [n for _, n in NODES]

# ---------- a) network position ----------
top = pd.read_csv(f"{R}/network_node_topology.csv").set_index("name")
edges = pd.read_csv("data/canonical_edges.tsv", sep="\t")
cores = pd.read_csv(f"{R}/v2/ffl_cores_coherence_corrected.csv")
mods = pd.read_csv(f"{R}/ffl_module_membership.csv").set_index("node")
msets = pd.read_csv(f"{R}/motif_node_sets.tsv", sep="\t")
prior = pd.read_csv(f"{R}/v5/node_prioritisation_full.csv").set_index("name")

# per-node core role counts + class + coherence breakdown
def core_stats(n):
    r = cores[cores.regulator == n]; m = cores[cores.intermediate == n]; t = cores[cores.target == n]
    sub = pd.concat([r, m, t])
    cls = sub["class"].value_counts().to_dict()
    coh = sub["coherence"].value_counts().to_dict()
    return dict(ffl_cores=len(sub), cores_as_regulator=len(r), cores_as_intermediate=len(m),
                cores_as_target=len(t),
                cores_composite=cls.get("Composite-FFL", 0), cores_miRNA_FFL=cls.get("miRNA-FFL", 0),
                cores_TF_FFL=cls.get("TF-FFL", 0),
                cores_coherent=sum(v for k, v in coh.items() if k.startswith("C")),
                cores_incoherent=sum(v for k, v in coh.items() if k.startswith("I")),
                cores_unresolved=coh.get("unresolved", 0),
                coherence_profile="; ".join(f"{k}={v}" for k, v in sorted(coh.items()) if k != "unresolved"))

# ---------- c) expression ----------
de = pd.read_csv(f"{R}/BRCA_DEX_ALL_nodes.csv").set_index("Gene")
geo = pd.read_csv(f"{R}/external_GEO_DE_vs_TCGA.csv")
cptac = pd.read_csv(f"{R}/multiomics/cptac_mRNA_protein_concordance.csv").set_index("gene")

def geo_stats(n, logfc_tcga):
    g = geo[geo.feature == n]
    if len(g) == 0 or pd.isna(logfc_tcga):
        return dict(geo_n_cohorts=0, geo_n_same_dir=np.nan, geo_n_same_dir_sig=np.nan,
                    geo_detail="")
    same = np.sign(g.logFC_ext) == np.sign(logfc_tcga)
    sig = g.q_ext < 0.05
    det = "; ".join(f"{row.gse}: logFC {row.logFC_ext:+.2f}, q={row.q_ext:.2g}" for row in g.itertuples())
    return dict(geo_n_cohorts=len(g), geo_n_same_dir=int(same.sum()),
                geo_n_same_dir_sig=int((same & sig).sum()), geo_detail=det)

# ---------- d) regulation of the node ----------
meth = pd.read_csv(f"{R}/multiomics/brca_methylation_nodes.csv").set_index("node")
cn = pd.read_csv(f"{R}/multiomics/brca_copynumber_nodes.csv").set_index("node")
mut = pd.read_csv(f"{R}/multiomics/brca_mutation_frequency.csv").set_index("gene")

# ---------- e) clinical ----------
cox = pd.read_csv(f"{R}/survival_cox_hubs.csv")
mb = pd.read_csv(f"{R}/external_METABRIC_cox.csv")

def cox_get(df, feat, endpoint, hrcol="HR_per_SD", pcol="p_value", qcol="q_value"):
    s = df[(df.iloc[:, 0] == feat) & (df.endpoint == endpoint)]
    if len(s) == 0:
        return (np.nan,) * 5
    s = s.iloc[0]
    return (s[hrcol], s["CI_low"], s["CI_high"], s[pcol], s.get(qcol, np.nan))

# ---------- f) function ----------
dep = pd.read_csv(f"{R}/multiomics/depmap_essentiality_nodes.csv").set_index("gene")
dgi = pd.read_csv(f"{OUT}/node_druggability_dgidb_summary.csv").set_index("gene")

rows = []
for typ, n in NODES:
    d = dict(node=n, node_class=typ)
    # topology
    t = top.loc[n]
    d.update(degree=int(t.degree), indegree=int(t.indegree), outdegree=int(t.outdegree),
             betweenness_raw=float(t.betweenness), clustering=float(t.clustering),
             closeness=float(t.closeness))
    d.update(core_stats(n))
    ms = sorted(msets[msets.node == n].motif_set.unique())
    d["motif_sets"] = ";".join(ms)
    mm = mods.loc[n]
    d.update(in_3node=bool(mm.in_3node), in_4node_set=bool(mm.in_4node_set),
             in_5node_set=bool(mm.in_5node_set), in_6node_set=bool(mm.in_6node_set),
             higher_order_only=bool(mm.higher_order_only),
             in_MS_exemplar_module=bool(mm.in_MS_exemplar_module))
    # ranks
    p = prior.loc[n]
    d.update(priority_score=p.priority, rank_within_type=int(p.rank_within_type),
             rank_overall=int(p.rank_overall), n_domains=int(p.n_domains),
             n_mirna_target_edges=p.n_edges, frac_strong_evidence=p.frac_strong)
    # expression
    lf = de.loc[n, "logFC"] if n in de.index else np.nan
    fdr = de.loc[n, "adj.P.Val"] if n in de.index else np.nan
    d.update(tcga_logFC=lf, tcga_FDR=fdr)
    d.update(geo_stats(n, lf))
    if n in cptac.index:
        d.update(cptac_mrna_prot_rho=cptac.loc[n, "rho"], cptac_mrna_prot_p=cptac.loc[n, "p"],
                 cptac_mrna_prot_fdr=cptac.loc[n, "fdr"], cptac_n=cptac.loc[n, "n"])
    else:
        d.update(cptac_mrna_prot_rho=np.nan, cptac_mrna_prot_p=np.nan,
                 cptac_mrna_prot_fdr=np.nan, cptac_n=np.nan)
    # regulation
    if n in meth.index:
        m = meth.loc[n]
        d.update(meth_n_probes=m.n_probes, meth_beta_normal=m.beta_normal, meth_beta_tumour=m.beta_tumour,
                 meth_delta_beta=m.delta_beta, meth_fdr=m.wilcox_fdr,
                 meth_expr_rho=m.meth_expr_rho, meth_expr_fdr=m.meth_expr_fdr)
    if n in cn.index:
        c = cn.loc[n]
        d.update(cn_locus=c.locus, cn_chrom=c.chrom, cn_frac_amp=c.frac_amp, cn_frac_del=c.frac_del,
                 cn_frac_highamp=c.frac_highamp, cn_frac_deepdel=c.frac_deepdel,
                 cn_expr_rho=c.cn_expr_rho, cn_expr_fdr=c.cn_expr_fdr, cn_all_loci=c.all_loci)
    if n in mut.index:
        d.update(mut_n=mut.loc[n, "n_mutated"], mut_n_samples=mut.loc[n, "n_samples"],
                 mut_freq=mut.loc[n, "mut_freq"])
    else:
        d.update(mut_n=np.nan, mut_n_samples=np.nan, mut_freq=np.nan)
    # clinical
    for ep in ["OS", "PFI", "DSS"]:
        hr, lo, hi, pv, qv = cox_get(cox, n, ep)
        d.update({f"tcga_{ep}_HR": hr, f"tcga_{ep}_CIlow": lo, f"tcga_{ep}_CIhigh": hi,
                  f"tcga_{ep}_p": pv, f"tcga_{ep}_q": qv})
    for ep in ["OS", "RFS"]:
        hr, lo, hi, pv, qv = cox_get(mb, n, ep)
        d.update({f"metabric_{ep}_HR": hr, f"metabric_{ep}_CIlow": lo, f"metabric_{ep}_CIhigh": hi,
                  f"metabric_{ep}_p": pv, f"metabric_{ep}_q": qv})
    # function
    if n in dep.index:
        dd = dep.loc[n]
        d.update(depmap_n_breast=dd.n_breast_screened, depmap_mean_chronos=dd.mean_chronos_breast,
                 depmap_frac_essential=dd.frac_essential_breast,
                 depmap_n_essential=dd.n_essential_breast, depmap_class=dd.essential_class,
                 depmap_selective_breast=dd.selective_breast,
                 depmap_mean_chronos_other=dd.mean_chronos_other,
                 depmap_wilcox_p_breast_vs_other=dd.wilcox_p_breast_vs_other)
    if n in dgi.index:
        g = dgi.loc[n]
        d.update(dgidb_n_drugs=g.n_drugs, dgidb_n_approved=g.n_approved,
                 dgidb_n_antineoplastic=g.n_antineoplastic,
                 dgidb_n_approved_antineoplastic=g.n_approved_antineoplastic,
                 dgidb_interaction_types=g.interaction_types,
                 dgidb_example_approved_antineoplastic=g.top_approved_antineoplastic)
    rows.append(d)

out = pd.DataFrame(rows)
out.to_csv(f"{OUT}/node_compendium_table.csv", index=False)
print(out.shape)
print(out[["node","node_class","degree","indegree","outdegree","ffl_cores","cores_as_regulator",
           "cores_as_intermediate","cores_as_target","motif_sets","tcga_logFC"]].to_string(index=False))
