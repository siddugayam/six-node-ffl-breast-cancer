#!/usr/bin/env python3
"""
40_node_compendium_assemble.py
Assemble every quantitative field for the 30 prioritised nodes (10 TF, 10 gene, 10 miRNA)
into results/v5/node_compendium_table.csv, plus full and top-N partner tables.
All values are read from existing verified result files or computed here from
data/canonical_edges.tsv; nothing is hard-coded.
"""
import os, sys, json, math
import pandas as pd, numpy as np, networkx as nx

ROOT = "/path/to/revision"
D   = os.path.join(ROOT, "data")
R   = os.path.join(ROOT, "results")
V2  = os.path.join(R, "v2")
MO  = os.path.join(R, "multiomics")
V5  = os.path.join(R, "v5")
os.makedirs(V5, exist_ok=True)

log = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.append(s)

# ---------------------------------------------------------------- 0. node set
top10 = pd.read_csv(os.path.join(V5, "top10_per_class.csv"))
# Recompute ranks among RANKABLE nodes only. The deposited rank columns in
# node_prioritisation_full.csv sort the four nodes with <3 evidence domains (priority = NA)
# to the top, which offsets every miRNA rank by 4; PRIORITISATION.md quotes the corrected
# ranks, so we recompute them here rather than propagate the offset.
FULL = pd.read_csv(os.path.join(V5, "node_prioritisation_full.csv"))
RK = FULL.dropna(subset=["priority"]).copy()
RK["rank_overall_fixed"] = RK["priority"].rank(ascending=False, method="min").astype(int)
RK["rank_within_type_fixed"] = RK.groupby("type")["priority"].rank(ascending=False,
                                                                  method="min").astype(int)
say("rankable nodes:", len(RK), RK.groupby("type").size().to_dict())
RK = RK.set_index("name")
NODES = top10[["name", "type", "priority"]].copy()
NODES["rank_overall"] = NODES["name"].map(RK["rank_overall_fixed"])
NODES["rank_within_type"] = NODES["name"].map(RK["rank_within_type_fixed"])
say("nodes:", len(NODES), NODES["type"].value_counts().to_dict())
nodeset = set(NODES["name"])

# ---------------------------------------------------------------- 1. network position
E = pd.read_csv(os.path.join(D, "canonical_edges.tsv"), sep="\t")
say("edges:", len(E), "unique nodes:", len(set(E.source) | set(E.target)))
G = nx.DiGraph()
G.add_nodes_from(pd.read_csv(os.path.join(D, "canonical_nodes.tsv"), sep="\t")["name"])
for s, t in zip(E.source, E.target):
    G.add_edge(s, t)
say("graph:", G.number_of_nodes(), "nodes,", G.number_of_edges(), "directed edges")
btw_raw  = nx.betweenness_centrality(G, normalized=False)
btw_norm = nx.betweenness_centrality(G, normalized=True)
outdeg = dict(G.out_degree()); indeg = dict(G.in_degree())
U = nx.Graph(); U.add_nodes_from(G.nodes()); U.add_edges_from(G.edges())
btw_undir = nx.betweenness_centrality(U, normalized=False)
say("undirected simple graph:", U.number_of_nodes(), "nodes,", U.number_of_edges(), "edges")
# betweenness/degree as used by the prioritisation composite (reconstructed network,
# results/network_topology_hubs.csv) - reported alongside for full transparency
TOPO = pd.read_csv(os.path.join(R, "network_topology_hubs.csv")).set_index("name")

# ---------------------------------------------- 1b. higher-order module membership
HO = pd.read_csv(os.path.join(R, "ffl_higher_order.csv"),
                 usecols=["n_nodes", "members", "source", "sink", "ffl_class", "listing"])
say("higher-order module records:", len(HO),
    HO.groupby(["n_nodes", "listing"]).size().to_dict())
ho_count = {}   # (node, n) -> count ; (node, n, role) -> count ; (node, n, "cls", class)
for n_nodes, members, src, snk, cls in zip(HO.n_nodes, HO.members, HO.source, HO.sink, HO.ffl_class):
    for m in members.split(";"):
        if m in nodeset:
            ho_count[(m, n_nodes)] = ho_count.get((m, n_nodes), 0) + 1
            role = "source" if m == src else ("sink" if m == snk else "internal")
            ho_count[(m, n_nodes, role)] = ho_count.get((m, n_nodes, role), 0) + 1
            ho_count[(m, n_nodes, "cls", cls)] = ho_count.get((m, n_nodes, "cls", cls), 0) + 1
del HO

# --------------------------------------------------------- 2. FFL core roles
FFL = pd.read_csv(os.path.join(V2, "ffl_cores_coherence_corrected.csv"))
say("FFL 3-node cores:", len(FFL), FFL["class"].value_counts().to_dict())
COH_COHERENT   = {"C1", "C2", "C3", "C4"}
COH_INCOHERENT = {"I1", "I2", "I3", "I4"}

def ffl_stats(n):
    reg = FFL[FFL.regulator == n]; imd = FFL[FFL.intermediate == n]; tgt = FFL[FFL.target == n]
    allc = pd.concat([reg, imd, tgt])
    cls = allc["class"].value_counts().to_dict()
    coh = allc["coherence"].value_counts().to_dict()
    nc = sum(v for k, v in coh.items() if k in COH_COHERENT)
    ni = sum(v for k, v in coh.items() if k in COH_INCOHERENT)
    nu = coh.get("unresolved", 0)
    resolved = nc + ni
    top = sorted(((k, v) for k, v in coh.items() if k != "unresolved"),
                 key=lambda x: -x[1])
    return dict(
        ffl_total=len(allc), ffl_as_regulator=len(reg), ffl_as_intermediate=len(imd),
        ffl_as_target=len(tgt),
        ffl_composite=cls.get("Composite-FFL", 0), ffl_mirnaFFL=cls.get("miRNA-FFL", 0),
        ffl_tfFFL=cls.get("TF-FFL", 0),
        coh_n_coherent=nc, coh_n_incoherent=ni, coh_n_unresolved=nu,
        coh_frac_coherent_of_resolved=(nc / resolved) if resolved else np.nan,
        coh_dominant_type=(top[0][0] if top else ""),
        coh_type_counts="; ".join(f"{k}:{v}" for k, v in top) if top else "")

# ------------------------------------------------- 3. motif sets and modules
MS = pd.read_csv(os.path.join(R, "motif_node_sets.tsv"), sep="\t")
motif_of = MS.groupby("node")["motif_set"].apply(lambda s: sorted(set(s))).to_dict()
MM = pd.read_csv(os.path.join(R, "ffl_module_membership.csv")).set_index("node")

# ------------------------------------------------------------- 4. partner table
EC = pd.read_csv(os.path.join(R, "edge_correlation.csv"))
ET = pd.read_csv(os.path.join(D, "edge_evidence_tier.tsv"), sep="\t")
ET["key"] = ET.source + "|" + ET.target
tier_of = dict(zip(ET.key, ET.tier))
dbs_of  = dict(zip(ET.key, ET.databases))
sup_of  = dict(zip(ET.key, ET.support_types))
EC["key"] = EC.source + "|" + EC.target
SA = pd.read_csv(os.path.join(MO, "stromal_adjusted_edge_rho.csv"))
SA["key"] = SA.source + "|" + SA.target
sa_of = dict(zip(SA.key, SA.rho_adj))
say("stromal(CAF)-adjusted rho available for", len(SA), "edges by type:",
    SA.edge_type.value_counts().to_dict())
EC["rho_caf_adj"] = EC["key"].map(sa_of)
EC["tier2"] = EC["key"].map(tier_of)
EC["databases"] = EC["key"].map(dbs_of)
EC["support_types"] = EC["key"].map(sup_of)

rows = []
for _, r in EC.iterrows():
    for node, partner, role in ((r.source, r.target, "regulates"),
                                (r.target, r.source, "regulated_by")):
        if node in nodeset:
            rows.append(dict(
                node=node, partner=partner, role=role, edge_type=r.edge_type,
                partner_type=(r.target_type if role == "regulates" else r.source_type),
                annotation=r.annotation, predicted_sign=r.predicted_sign,
                sign_source=r.sign_source, evidence_tier=r.tier2,
                databases=r.databases, support_types=r.support_types,
                rho=r.rho, p=r.p, fdr=r.fdr, n_samples=r.n_samples,
                concordant=r.concordant, rho_pcadj=r.rho_pcadj, p_pcadj=r.p_pcadj,
                rho_caf_adj=r.rho_caf_adj))
P = pd.DataFrame(rows)
P.to_csv(os.path.join(V5, "node_compendium_partners_full.csv"), index=False)
say("partner rows (all edges of the 30 nodes):", len(P))

tier_rank = {"strong": 0, "weak": 1, "predicted_only": 2}
P["_tr"] = P["evidence_tier"].map(tier_rank).fillna(1.5)
P["_ar"] = P["rho"].abs()
TOPN = 8   # per direction, so up to 16 partners per node
tops = (P.sort_values(["node", "role", "_tr", "_ar"], ascending=[True, True, True, False])
          .groupby(["node", "role"]).head(TOPN)
          .drop(columns=["_tr", "_ar"]))
tops.to_csv(os.path.join(V5, "node_compendium_top_partners.csv"), index=False)
say("top-partner rows written:", len(tops))

def partner_stats(n):
    d = P[P.node == n]
    ds = d.drop_duplicates(subset=["partner", "edge_type", "role"])
    mt = ds[ds.edge_type == "miRNA_target"]
    nstrong = int((mt.evidence_tier == "strong").sum())
    rhos = ds["rho"].dropna()
    sig  = ds[(ds.fdr < 0.05)]
    cc = ds["concordant"].dropna()
    ad = ds.dropna(subset=["rho_caf_adj", "rho"])
    return dict(n_partner_edges=len(ds), n_unique_partners=ds.partner.nunique(),
                n_mirna_target_edges=len(mt), n_strong_edges=nstrong,
                frac_strong=(nstrong / len(mt)) if len(mt) else np.nan,
                median_abs_rho=float(rhos.abs().median()) if len(rhos) else np.nan,
                n_edges_fdr05=len(sig),
                n_edges_sign_tested=len(cc),
                frac_edges_sign_concordant=(float(cc.mean()) if len(cc) else np.nan),
                n_edges_caf_adjusted=len(ad),
                median_abs_rho_caf_adj=(float(ad["rho_caf_adj"].abs().median()) if len(ad) else np.nan))

# --------------------------------------------------------------- 5. expression
DEg = pd.read_csv(os.path.join(V2, "v2_DE_genes.csv")).set_index("feature")
DEm = pd.read_csv(os.path.join(V2, "v2_DE_mirnas.csv")).set_index("feature")
GEO = pd.read_csv(os.path.join(R, "external_GEO_DE.csv"))
COHORTS = ["GSE42568", "GSE45827", "GSE10780"]
geo_n = GEO.drop_duplicates("gse").set_index("gse")[["n_normal", "n_tumour"]].to_dict("index")
say("GEO cohorts:", geo_n)
CPTAC = pd.read_csv(os.path.join(MO, "cptac_mRNA_protein_concordance.csv")).set_index("gene")
RPPA  = pd.read_csv(os.path.join(MO, "rppa_protein_vs_mRNA.csv"))
RPPA_nophos = RPPA[~RPPA.is_phospho].drop_duplicates("gene").set_index("gene")

def expr_stats(n, typ):
    out = {}
    src = DEm if typ == "miRNA" else DEg
    if n in src.index:
        out["logFC_TCGA"] = float(src.loc[n, "logFC"])
        out["fdr_TCGA"]   = float(src.loc[n, "adj.P.Val"])
        out["P_TCGA"]     = float(src.loc[n, "P.Value"])
    g = GEO[GEO.feature == n]
    same = 0; tot = 0; samesig = 0
    for c in COHORTS:
        gg = g[g.gse == c]
        if len(gg):
            lfc = float(gg["logFC"].iloc[0]); pv = float(gg["adj.P.Val"].iloc[0])
            out[f"logFC_{c}"] = lfc; out[f"fdr_{c}"] = pv
            tot += 1
            if not np.isnan(out.get("logFC_TCGA", np.nan)) and np.sign(lfc) == np.sign(out["logFC_TCGA"]):
                same += 1
                if pv < 0.05: samesig += 1
    out["geo_n_cohorts"] = tot; out["geo_n_same_dir"] = same if tot else np.nan
    out["geo_n_same_dir_sig"] = samesig if tot else np.nan
    out["geo_frac_same_dir"] = (same / tot) if tot else np.nan
    if n in CPTAC.index:
        out["cptac_mrna_prot_rho"] = float(CPTAC.loc[n, "rho"])
        out["cptac_mrna_prot_p"]   = float(CPTAC.loc[n, "p"])
        out["cptac_mrna_prot_fdr"] = float(CPTAC.loc[n, "fdr"])
        out["cptac_n"]             = int(CPTAC.loc[n, "n"])
    if n in RPPA_nophos.index:
        out["rppa_rho"] = float(RPPA_nophos.loc[n, "rho"])
        out["rppa_n"]   = int(RPPA_nophos.loc[n, "n"])
        out["rppa_antibody"] = RPPA_nophos.loc[n, "antibody"]
    return out

# ------------------------------------------------- 6. regulation of the node
METH = pd.read_csv(os.path.join(MO, "brca_methylation_nodes.csv")).set_index("node")
CN   = pd.read_csv(os.path.join(MO, "brca_copynumber_nodes.csv")).set_index("node")
MUT  = pd.read_csv(os.path.join(MO, "brca_mutation_frequency.csv")).set_index("gene")

def reg_stats(n):
    o = {}
    if n in METH.index:
        m = METH.loc[n]
        o.update(meth_n_probes=int(m.n_probes), meth_beta_tumour=float(m.beta_tumour),
                 meth_beta_normal=float(m.beta_normal), meth_delta_beta=float(m.delta_beta),
                 meth_fdr=float(m.wilcox_fdr),
                 meth_expr_rho=(float(m.meth_expr_rho) if pd.notna(m.meth_expr_rho) else np.nan),
                 meth_expr_fdr=(float(m.meth_expr_fdr) if pd.notna(m.meth_expr_fdr) else np.nan))
    if n in CN.index:
        c = CN.loc[n]
        o.update(cn_frac_amp=float(c.frac_amp), cn_frac_del=float(c.frac_del),
                 cn_frac_highamp=float(c.frac_highamp), cn_frac_deepdel=float(c.frac_deepdel),
                 cn_locus=str(c.locus), cn_chrom=str(c.chrom), cn_n_samples=int(c.n_cn_samples),
                 cn_expr_rho=(float(c.cn_expr_rho) if pd.notna(c.cn_expr_rho) else np.nan),
                 cn_expr_fdr=(float(c.cn_expr_fdr) if pd.notna(c.cn_expr_fdr) else np.nan))
    if n in MUT.index:
        u = MUT.loc[n]
        o.update(mut_n=int(u.n_mutated), mut_n_samples=int(u.n_samples), mut_freq=float(u.mut_freq))
    return o

# -------------------------------------------------------------- 7. clinical
SV = pd.read_csv(os.path.join(R, "survival_cox_hubs.csv"))
MB = pd.read_csv(os.path.join(R, "external_METABRIC_cox.csv"))
MBr = pd.read_csv(os.path.join(MO, "metabric_survival_replication.csv"))

def surv_stats(n):
    o = {}
    for ep in ["OS", "PFI", "DSS"]:
        s = SV[(SV.feature == n) & (SV.endpoint == ep)]
        if len(s):
            s = s.iloc[0]
            o[f"tcga_{ep}_HR"] = float(s.HR_per_SD); o[f"tcga_{ep}_CIlow"] = float(s.CI_low)
            o[f"tcga_{ep}_CIhigh"] = float(s.CI_high); o[f"tcga_{ep}_p"] = float(s.p_value)
            o[f"tcga_{ep}_q"] = float(s.q_value); o[f"tcga_{ep}_n"] = int(s.n)
            o[f"tcga_{ep}_nevent"] = int(s.n_event)
    for ep in ["OS", "RFS"]:
        s = MB[(MB.feature == n) & (MB.endpoint == ep)]
        if len(s):
            s = s.iloc[0]
            o[f"metabric_{ep}_HR"] = float(s.HR_per_SD); o[f"metabric_{ep}_CIlow"] = float(s.CI_low)
            o[f"metabric_{ep}_CIhigh"] = float(s.CI_high); o[f"metabric_{ep}_p"] = float(s.p_value)
            o[f"metabric_{ep}_q"] = float(s.q_value); o[f"metabric_{ep}_n"] = int(s.n)
        s2 = MBr[(MBr.feature == n) & (MBr.endpoint == ep) & (MBr.model == "adjusted_for_CAF")]
        if len(s2):
            o[f"metabric_{ep}_HR_CAFadj"] = float(s2.iloc[0].HR_per_SD)
            o[f"metabric_{ep}_p_CAFadj"]  = float(s2.iloc[0].p)
    return o

# -------------------------------------------------------------- 8. function
DEP = pd.read_csv(os.path.join(MO, "depmap_essentiality_nodes.csv")).set_index("gene")
DG  = pd.read_csv(os.path.join(V5, "node_druggability_dgidb_summary.csv")).set_index("gene")
HD  = pd.read_csv(os.path.join(MO, "hub_druggability_gene_summary.csv")).set_index("gene_name")

def func_stats(n):
    o = {}
    if n in DEP.index:
        d = DEP.loc[n]
        o.update(chronos_mean_breast=float(d.mean_chronos_breast),
                 chronos_median_breast=float(d.median_chronos_breast),
                 chronos_sd_breast=float(d.sd_chronos_breast),
                 n_breast_lines=int(d.n_breast_screened),
                 n_essential_breast=int(d.n_essential_breast),
                 frac_essential_breast=float(d.frac_essential_breast),
                 essential_class=str(d.essential_class),
                 chronos_mean_other=float(d.mean_chronos_other),
                 selective_breast=bool(d.selective_breast))
    if n in DG.index:
        g = DG.loc[n]
        o.update(dgidb_in=bool(g.in_dgidb), dgidb_n_drugs=int(g.n_drugs),
                 dgidb_n_approved=int(g.n_approved),
                 dgidb_n_antineoplastic=int(g.n_antineoplastic),
                 dgidb_n_approved_antineoplastic=int(g.n_approved_antineoplastic),
                 dgidb_interaction_types=(str(g.interaction_types) if pd.notna(g.interaction_types) else ""),
                 dgidb_top_approved_antineoplastic=(str(g.top_approved_antineoplastic) if pd.notna(g.top_approved_antineoplastic) else ""))
    if n in HD.index:
        o["dgidb_hubtable_n_clinical_stage"] = int(HD.loc[n, "n_clinical_stage_drugs"])
    return o

# ------------------------------------------------ 8b. ExIR
EX  = pd.read_csv(os.path.join(R, "exir_classification.csv"))
EXP = pd.read_csv(os.path.join(R, "exir_primary_class.csv")).set_index("feature")
say("ExIR classes present for prioritised nodes:",
    EX[EX.feature.isin(nodeset)]["class"].value_counts().to_dict())

def exir_stats(n):
    o = {}
    if n in EXP.index:
        o["exir_primary_class"] = str(EXP.loc[n, "primary_class"])
    d = EX[EX.feature == n]
    for cl, tag in (("Driver", "driver"), ("Biomarker", "biomarker"),
                    ("Mediator", "mediator")):
        dd = d[d["class"] == cl]
        if len(dd):
            o[f"exir_{tag}_rank"] = int(dd.iloc[0]["rank"])
            o[f"exir_{tag}_score"] = float(dd.iloc[0]["score"])
            o[f"exir_{tag}_subtype"] = str(dd.iloc[0]["subtype"])
            o[f"exir_{tag}_padj"] = float(dd.iloc[0]["P_adj"])
            o[f"exir_{tag}_significant"] = bool(dd.iloc[0]["exir_significant"])
    return o

# -------------------------------------------------------------- 9. assemble
recs = []
for _, r in NODES.iterrows():
    n, typ = r["name"], r["type"]
    rec = dict(name=n, type=typ, composite_priority=r["priority"],
               rank_within_type=int(r["rank_within_type"]), rank_overall=int(r["rank_overall"]))
    rec.update(degree_total=indeg.get(n, 0) + outdeg.get(n, 0),
               in_degree=indeg.get(n, 0), out_degree=outdeg.get(n, 0),
               betweenness_directed=btw_raw.get(n, np.nan),
               betweenness_directed_norm=btw_norm.get(n, np.nan),
               betweenness_undirected=btw_undir.get(n, np.nan))
    if n in TOPO.index:
        rec.update(recon_degree_total=int(TOPO.loc[n, "degree_tot"]),
                   recon_degree_out=int(TOPO.loc[n, "degree_out"]),
                   recon_degree_in=int(TOPO.loc[n, "degree_in"]),
                   recon_betweenness=float(TOPO.loc[n, "betweenness"]),
                   hub_by_degree=bool(TOPO.loc[n, "hub_by_degree"]),
                   hub_by_betweenness=bool(TOPO.loc[n, "hub_by_betweenness"]),
                   hub_named_in_MS=bool(TOPO.loc[n, "hub_named_in_MS"]))
    for k in (4, 5, 6):
        rec[f"ho_n{k}_modules"] = ho_count.get((n, k), 0)
        for role in ("source", "internal", "sink"):
            rec[f"ho_n{k}_{role}"] = ho_count.get((n, k, role), 0)
    rec.update(ffl_stats(n))
    rec["motif_subnetworks_containing_node"] = "; ".join(motif_of.get(n, []))
    n3cls = []
    for lab, k in (("Composite-FFL", "ffl_composite"), ("miRNA-FFL", "ffl_mirnaFFL"),
                   ("TF-FFL", "ffl_tfFFL")):
        if rec.get(k, 0): n3cls.append(f"{lab}:{rec[k]}")
    rec["motif_classes_n3"] = "; ".join(n3cls)
    n4cls = [(c, ho_count.get((n, 4, "cls", c), 0))
             for c in ("Composite-FFL", "miRNA-FFL", "TF-FFL", "gene-FFL")]
    rec["motif_classes_n4"] = "; ".join(f"{c}:{v}" for c, v in n4cls if v)
    if n in MM.index:
        m = MM.loc[n]
        rec.update(in_3node=bool(m.in_3node), in_4node_set=bool(m.in_4node_set),
                   in_5node_set=bool(m.in_5node_set), in_6node_set=bool(m.in_6node_set),
                   higher_order_only=bool(m.higher_order_only),
                   in_MS_exemplar_module=bool(m.in_MS_exemplar_module))
    rec.update(partner_stats(n))
    rec.update(expr_stats(n, typ))
    rec.update(reg_stats(n))
    rec.update(surv_stats(n))
    rec.update(func_stats(n))
    rec.update(exir_stats(n))
    recs.append(rec)

T = pd.DataFrame(recs)
order = ["Gene", "TF", "miRNA"]
T["_o"] = T["type"].map({t: i for i, t in enumerate(order)})
T = T.sort_values(["_o", "rank_within_type"]).drop(columns="_o")
T.to_csv(os.path.join(V5, "node_compendium_table.csv"), index=False)
say("wrote node_compendium_table.csv:", T.shape)
say("columns:", len(T.columns))

with open(os.path.join(ROOT, "logs", "v5", "40_node_compendium_assemble.log"), "w") as f:
    f.write("\n".join(log) + "\n")
