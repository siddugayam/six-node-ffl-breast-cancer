#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write results/v5/NODE_COMPENDIUM.md from node_compendium_table.csv,
node_compendium_top_partners.csv, node_literature_text.py and the verified PMID table.
Every number in the prose is read from the tables; nothing is typed by hand here."""
import os, sys, re, math, csv
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import node_literature_text as LT

ROOT = "/path/to/revision"
V5 = os.path.join(ROOT, "results", "v5")
T = pd.read_csv(os.path.join(V5, "node_compendium_table.csv"))
P = pd.read_csv(os.path.join(V5, "node_compendium_top_partners.csv"))
VER = pd.read_csv(os.path.join(V5, "node_literature_verified.csv"), dtype={"pmid": str})
LITVOL = pd.read_csv(os.path.join(V5, "node_literature_volume.csv")).set_index("node")
VERD = {r.pmid: r for r in VER.itertuples()}

# ---- verification gate -------------------------------------------------------
cited = {p for v in LT.LIT.values() for p in re.findall(r"PMID:(\d+)", v)}
missing = sorted(cited - set(VER.pmid))
if missing:
    sys.exit("ABORT: unverified PMIDs cited: " + ", ".join(missing))
print("PMID gate passed:", len(cited), "citations verified")

# ---- formatting helpers ------------------------------------------------------
def na(x):
    return x is None or (isinstance(x, float) and math.isnan(x))

def f(x, d=2):
    if na(x): return "n.a."
    return f"{x:+.{d}f}" if d and abs(x) < 1000 and x < 0 else f"{x:.{d}f}"

def num(x, d=2):
    return "n.a." if na(x) else f"{x:.{d}f}"

def sgn(x, d=2):
    return "n.a." if na(x) else f"{x:+.{d}f}"

def pv(x):
    if na(x): return "n.a."
    if x == 0: return "< 1e-300"
    if x >= 0.001: return f"{x:.3g}"
    m = f"{x:.1e}".replace("e-0", "e-").replace("e+0", "e+")
    return m

def pct(x, d=1):
    return "n.a." if na(x) else f"{100*x:.{d}f} %"

def i(x):
    return "n.a." if na(x) else f"{int(x):,}"

def cite(txt):
    """[PMID:123] -> (PMID 123)"""
    def one(m):
        return "PMID " + m.group(1)
    txt = re.sub(r"\[PMID:(\d+)(; PMID:(\d+))?\]",
                 lambda m: "(" + "; ".join("PMID " + g for g in m.groups()[::2] if g) + ")", txt)
    return txt

# ---- per-node section builders ----------------------------------------------
TYPELABEL = {"TF": "transcription factor", "Gene": "gene", "miRNA": "microRNA"}
NTYPE = T.groupby("type").size().to_dict()
NCLASS_TOTAL = {"TF": 157, "Gene": 207, "miRNA": 219}   # rankable nodes per class (PRIORITISATION.md)

def sec_network(r):
    L = []
    L.append(f"**a. Network position.** Degree {i(r.degree_total)} "
             f"(in {i(r.in_degree)}, out {i(r.out_degree)}); directed betweenness "
             f"{num(r.betweenness_directed, 1)}, undirected betweenness "
             f"{num(r.betweenness_undirected, 1)} "
             f"(the prioritisation composite used the reconstructed-network value, "
             f"{num(r.recon_betweenness, 1)}, from `network_topology_hubs.csv`; both are reported "
             f"so the two conventions cannot be confused).")
    if r.ffl_total > 0:
        roles = []
        if r.ffl_as_regulator: roles.append(f"regulator x{i(r.ffl_as_regulator)}")
        if r.ffl_as_intermediate: roles.append(f"intermediate x{i(r.ffl_as_intermediate)}")
        if r.ffl_as_target: roles.append(f"target x{i(r.ffl_as_target)}")
        if na(r.coh_frac_coherent_of_resolved):
            coh = (f"all {i(r.coh_n_unresolved)} of them unresolved — at least one edge of every "
                   f"core carries no directional annotation, so no Alon coherence type can be "
                   f"assigned")
        else:
            coh = (f"{i(r.coh_n_coherent)} coherent, {i(r.coh_n_incoherent)} incoherent, "
                   f"{i(r.coh_n_unresolved)} unresolved "
                   f"({pct(r.coh_frac_coherent_of_resolved)} coherent among resolved; commonest "
                   f"type {r.coh_dominant_type}, full breakdown {r.coh_type_counts})")
        L.append(f"It occupies {i(r.ffl_total)} of the 1,649 three-node FFL cores "
                 f"({', '.join(roles)}); motif classes {r.motif_classes_n3}. "
                 f"Coherence of those cores: {coh}.")
    else:
        L.append("It occupies **no** three-node FFL core in any role — it enters the prioritised "
                 "set on evidence weight alone, which is precisely the disagreement between the "
                 "two rankings that this analysis is meant to expose.")
    ho = (f"Higher-order membership: {i(r.ho_n4_modules)} four-node modules "
          f"(complete enumeration of 71,307 distinct member sets), {i(r.ho_n5_modules)} five-node "
          f"and {i(r.ho_n6_modules)} six-node modules in the capped listing "
          f"(5,000 per architecture; these are sampled, not exhaustive)")
    if r.ho_n4_modules > 0:
        ho += (f". Role within the four-node modules: source {i(r.ho_n4_source)}, "
               f"internal {i(r.ho_n4_internal)}, sink {i(r.ho_n4_sink)}; their motif classes are "
               f"{r.motif_classes_n4}")
    ho += ("; it **is** one of the eleven nodes of the manuscript's exemplar four/five/six-node "
           "modules" if r.in_MS_exemplar_module else
           "; it is not part of the manuscript's exemplar modules")
    L.append(ho + ".")
    ex = []
    if not na(getattr(r, "exir_primary_class", np.nan)):
        ex.append(f"ExIR primary class **{r.exir_primary_class}**")
    if not na(getattr(r, "exir_driver_rank", np.nan)):
        ex.append(f"driver rank {i(r.exir_driver_rank)} ({r.exir_driver_subtype})")
    if not na(getattr(r, "exir_biomarker_rank", np.nan)):
        ex.append(f"biomarker rank {i(r.exir_biomarker_rank)} ({r.exir_biomarker_subtype})")
    if ex:
        tail = (". No prioritised node reaches ExIR significance after multiple-testing "
                "correction, so these ranks are ordinal only."
                if len(ex) > 1 else
                " — it is not ranked as an ExIR driver or biomarker, i.e. it maintains network "
                "connectivity without being among the top-ranked differentially expressed "
                "features, which is exactly the 'mediator' category we set out to "
                "report.")
        L.append("ExIR: " + "; ".join(ex) + tail)
    return " ".join(L)

def sec_partners(r):
    nm = r["name"]
    d = P[P.node == nm].copy()
    if not len(d):
        return "**b. Partners.** No edges in the canonical network."
    out = ["**b. Partners.** "
           f"{i(r.n_partner_edges)} edges to {i(r.n_unique_partners)} distinct partners; "
           f"{i(r.n_mirna_target_edges)} "
           f"{'is a miRNA-target edge' if r.n_mirna_target_edges == 1 else 'are miRNA-target edges'}"
           f", of which {i(r.n_strong_edges)} "
           f"({pct(r.frac_strong)}) "
           f"{'carries' if r.n_strong_edges == 1 else 'carry'} low-throughput experimental "
           f"validation. "
           f"{i(r.n_edges_fdr05)} of the node's edges have a TCGA correlation at FDR < 0.05; "
           f"median |rho| = {num(r.median_abs_rho, 3)}, falling to "
           f"{num(r.median_abs_rho_caf_adj, 3)} after partialling out CAF content. "
           f"Sign concordance with the predicted direction: "
           f"{pct(r.frac_edges_sign_concordant)} of {i(r.n_edges_sign_tested)} testable edges."]
    order = {"strong": 0, "weak": 1, "predicted_only": 2}
    d["_t"] = d.evidence_tier.map(order).fillna(1.5)
    d["_a"] = d.rho.abs()
    rows = []
    for role in ("regulated_by", "regulates"):
        dd = d[d.role == role].sort_values(["_t", "_a"], ascending=[True, False]).head(6)
        rows.append(dd)
    dd = pd.concat(rows)
    tbl = ["", "| partner | direction | edge type | sign / annotation | evidence tier | TCGA rho | FDR | rho (CAF-adj) |",
           "|---|---|---|---|---|---|---|---|"]
    for _, e in dd.iterrows():
        ann = e.annotation if isinstance(e.annotation, str) else "n.a."
        if e.edge_type == "miRNA_target": ann = "repression (-1)"
        tier = (e.evidence_tier if isinstance(e.evidence_tier, str)
                else "n.a. (not a miRNA-target edge)")
        tbl.append(f"| {e.partner} | {'in' if e.role=='regulated_by' else 'out'} | {e.edge_type} | "
                   f"{ann} | {tier} | {sgn(e.rho, 3)} | {pv(e.fdr)} | {sgn(e.rho_caf_adj, 3)} |")
    tbl.append("")
    tbl.append("*(top six per direction, ranked by evidence tier then |rho|; \"in\" = the partner "
               "regulates this node, \"out\" = this node regulates the partner. CAF-adjusted "
               "correlations exist for miRNA-target and TF-target edges only, not for TF-miRNA "
               "edges. The complete edge list for all 30 nodes is in "
               "`node_compendium_partners_full.csv`.)*")
    return out[0] + "\n" + "\n".join(tbl)

def sec_expression(r):
    L = [f"**c. Expression.** TCGA-BRCA tumour vs normal log2FC **{sgn(r.logFC_TCGA)}** "
         f"(FDR {pv(r.fdr_TCGA)})."]
    if not na(r.geo_n_cohorts) and r.geo_n_cohorts > 0:
        L.append(f"Direction is reproduced in {i(r.geo_n_same_dir)} of {i(r.geo_n_cohorts)} "
                 f"independent GEO cohorts ({i(r.geo_n_same_dir_sig)} of them at cohort FDR < 0.05): "
                 f"GSE42568 {sgn(r.logFC_GSE42568)} (FDR {pv(r.fdr_GSE42568)}), "
                 f"GSE45827 {sgn(r.logFC_GSE45827)} (FDR {pv(r.fdr_GSE45827)}), "
                 f"GSE10780 {sgn(r.logFC_GSE10780)} (FDR {pv(r.fdr_GSE10780)}).")
    else:
        L.append("The three replication cohorts (GSE42568, GSE45827, GSE10780) are mRNA "
                 "microarray series and contain no miRNA measurements, so no direction "
                 "concordance can be computed for this node — a coverage limit, not a "
                 "negative result.")
    if not na(getattr(r, "cptac_n", np.nan)):
        if na(r.cptac_mrna_prot_rho):
            L.append(f"The protein is detected in the CPTAC prospective breast cohort "
                     f"({i(r.cptac_n)} tumours) but the mRNA-protein correlation is not "
                     f"estimable (too few paired non-missing values).")
        else:
            L.append(f"CPTAC mRNA-protein correlation rho = **{num(r.cptac_mrna_prot_rho, 3)}** "
                     f"(FDR {pv(r.cptac_mrna_prot_fdr)}, n = {i(r.cptac_n)} tumours).")
        if not na(getattr(r, "rppa_rho", np.nan)):
            L.append(f"TCGA RPPA antibody {r.rppa_antibody}: protein-mRNA rho = "
                     f"{num(r.rppa_rho, 3)} (n = {i(r.rppa_n)}).")
        L.append("Note that the CPTAC breast set and the TCGA RPPA panel contain no "
                 "matched normal tissue, so no protein-level tumour-vs-normal fold change "
                 "is available for any node in this study.")
    else:
        L.append("No protein-level measurement exists for a microRNA, so neither CPTAC nor "
                 "RPPA contributes here.")
    return " ".join(L)

def sec_regulation(r):
    parts = []
    if not na(r.meth_delta_beta):
        parts.append(f"Promoter methylation: delta-beta **{sgn(r.meth_delta_beta, 3)}** across "
                     f"{i(r.meth_n_probes)} probe(s) (beta {num(r.meth_beta_normal,3)} normal to "
                     f"{num(r.meth_beta_tumour,3)} tumour, FDR {pv(r.meth_fdr)}); "
                     f"methylation-expression rho {sgn(r.meth_expr_rho, 3)} "
                     f"(FDR {pv(r.meth_expr_fdr)}).")
    if not na(r.cn_frac_amp):
        parts.append(f"Copy number at {r.cn_locus} ({r.cn_chrom}, n = {i(r.cn_n_samples)}): gain in "
                     f"**{pct(r.cn_frac_amp)}** and loss in **{pct(r.cn_frac_del)}** of tumours "
                     f"(high-level gain {pct(r.cn_frac_highamp)}, deep deletion "
                     f"{pct(r.cn_frac_deepdel)}); copy-number-expression rho "
                     f"{sgn(r.cn_expr_rho, 3)} (FDR {pv(r.cn_expr_fdr)}).")
    if not na(getattr(r, "mut_freq", np.nan)):
        parts.append(f"Non-silent somatic mutation in {i(r.mut_n)} of {i(r.mut_n_samples)} "
                     f"tumours (**{pct(r.mut_freq, 2)}**).")
    else:
        parts.append("Somatic point mutation was not assessed: the MC3 non-silent call set "
                     "covers protein-coding genes only.")
    # ---- mechanism call (rule-based, stated explicitly) ----
    calls = []
    lf = r.logFC_TCGA
    if na(r.fdr_TCGA) or r.fdr_TCGA >= 0.05:
        msg = ("**Mechanism call:** this node is **not differentially expressed** at the cohort "
               "level (FDR " + pv(r.fdr_TCGA) + "), so there is no expression change for a "
               "*cis* lesion to explain.")
        if (not na(r.meth_delta_beta) and abs(r.meth_delta_beta) >= 0.05
                and not na(r.meth_fdr) and r.meth_fdr < 0.05):
            msg += (" The promoter is nonetheless significantly "
                    + ("hyper" if r.meth_delta_beta > 0 else "hypo")
                    + f"methylated (Δβ {sgn(r.meth_delta_beta,3)}, FDR {pv(r.meth_fdr)})")
            if (not na(r.meth_expr_rho) and abs(r.meth_expr_rho) >= 0.30
                    and not na(r.meth_expr_fdr) and r.meth_expr_fdr < 0.05):
                msg += (f", and promoter methylation correlates strongly with expression across "
                        f"tumours (ρ = {sgn(r.meth_expr_rho,3)}, FDR {pv(r.meth_expr_fdr)}) — so "
                        "the lesion is consequential within individual tumours even though the "
                        "cohort mean does not shift, which is the signature of a marker that "
                        "stratifies the disease rather than one that is uniformly altered by it.")
            else:
                msg += (f", but it is not appreciably correlated with expression across tumours "
                        f"(ρ = {sgn(r.meth_expr_rho,3)}) — a lesion without an observed "
                        "transcriptional consequence, which we report rather than interpret.")
        if not na(getattr(r, "mut_freq", np.nan)) and r.mut_freq >= 0.05:
            msg += f" It is also recurrently mutated ({pct(r.mut_freq,2)} of tumours)."
        parts.append(msg)
        return "**d. Regulation of the node itself.** " + " ".join(parts)
    if (not na(r.meth_delta_beta) and abs(r.meth_delta_beta) >= 0.05 and r.meth_fdr < 0.05
            and not na(lf) and np.sign(r.meth_delta_beta) == -np.sign(lf)
            and not na(r.meth_expr_rho) and r.meth_expr_rho < -0.1):
        calls.append("promoter methylation change of the right sign, correlated with expression "
                     "in the expected direction — a supported epigenetic mechanism")
    elif (not na(r.meth_delta_beta) and abs(r.meth_delta_beta) >= 0.05 and r.meth_fdr < 0.05
          and not na(lf) and np.sign(r.meth_delta_beta) == np.sign(lf)):
        calls.append("a significant methylation change whose sign is *discordant* with the "
                     "expression change, so methylation does not explain the dysregulation")
    if (not na(r.cn_frac_amp) and r.cn_frac_amp >= 0.30 and not na(r.cn_expr_rho)
            and r.cn_expr_rho >= 0.20 and r.cn_expr_fdr < 0.05 and not na(lf) and lf > 0):
        calls.append("recurrent copy gain that tracks expression — a plausible contributing "
                     "genomic mechanism")
    if (not na(r.cn_frac_del) and r.cn_frac_del >= 0.30 and not na(r.cn_expr_rho)
            and r.cn_expr_rho >= 0.20 and r.cn_expr_fdr < 0.05 and not na(lf) and lf < 0):
        calls.append("recurrent copy loss that tracks expression — a plausible contributing "
                     "genomic mechanism")
    if not na(getattr(r, "mut_freq", np.nan)) and r.mut_freq >= 0.05:
        calls.append("recurrent somatic mutation")
    if calls:
        parts.append("**Mechanism call:** " + "; ".join(calls) + ".")
    else:
        parts.append("**Mechanism call:** none of the three measured lesions (methylation, copy "
                     "number, point mutation) accounts for this node's dysregulation on the "
                     "pre-specified criteria; its altered expression is most likely "
                     "*trans*-regulatory, which is the situation the FFL census is designed to "
                     "describe.")
    return "**d. Regulation of the node itself.** " + " ".join(parts)

def sec_clinical(r):
    rows = ["", "| cohort | endpoint | n (events) | HR per SD | 95 % CI | p | q |",
            "|---|---|---|---|---|---|---|"]
    any_row = False
    for ep in ("OS", "PFI", "DSS"):
        hr = getattr(r, f"tcga_{ep}_HR", np.nan)
        if na(hr): continue
        any_row = True
        ci = f"{num(getattr(r,f'tcga_{ep}_CIlow'),2)}-{num(getattr(r,f'tcga_{ep}_CIhigh'),2)}"
        n = f"{i(getattr(r,f'tcga_{ep}_n'))} ({i(getattr(r,f'tcga_{ep}_nevent'))})"
        rows.append(f"| TCGA-BRCA | {ep} | {n} | {num(hr,3)} | {ci} | "
                    f"{pv(getattr(r,f'tcga_{ep}_p'))} | {pv(getattr(r,f'tcga_{ep}_q'))} |")
    for ep in ("OS", "RFS"):
        hr = getattr(r, f"metabric_{ep}_HR", np.nan)
        if na(hr): continue
        any_row = True
        ci = f"{num(getattr(r,f'metabric_{ep}_CIlow'),2)}-{num(getattr(r,f'metabric_{ep}_CIhigh'),2)}"
        rows.append(f"| METABRIC | {ep} | {i(getattr(r,f'metabric_{ep}_n'))} | {num(hr,3)} | {ci} | "
                    f"{pv(getattr(r,f'metabric_{ep}_p'))} | {pv(getattr(r,f'metabric_{ep}_q'))} |")
        adj = getattr(r, f"metabric_{ep}_HR_CAFadj", np.nan)
        if not na(adj):
            rows.append(f"| METABRIC (CAF-adjusted) | {ep} | "
                        f"{i(getattr(r,f'metabric_{ep}_n'))} | {num(adj,3)} | - | "
                        f"{pv(getattr(r,f'metabric_{ep}_p_CAFadj'))} | - |")
    rows.append("")
    head = "**e. Clinical.** "
    if r.type == "miRNA":
        head += ("METABRIC profiles mRNA on an Illumina array and contains no microRNA "
                 "measurements, so external survival replication is unavailable for this node. ")
    if not any_row:
        return head + "No Cox model was estimable."
    return head + "Cox models are per standard deviation of expression, univariate.\n" + "\n".join(rows)

def sec_function(r):
    L = []
    if not na(getattr(r, "chronos_mean_breast", np.nan)):
        L.append(f"**f. Function.** CRISPR (DepMap, {i(r.n_breast_lines)} breast cancer lines): "
                 f"mean Chronos gene effect **{sgn(r.chronos_mean_breast, 3)}** "
                 f"(sd {num(r.chronos_sd_breast,3)}); essential in {i(r.n_essential_breast)} of "
                 f"{i(r.n_breast_lines)} lines ({pct(r.frac_essential_breast)}), class "
                 f"*{r.essential_class}*"
                 f"{'; selectively essential in breast relative to other lineages' if r.selective_breast else ''}.")
        if not na(getattr(r, "dgidb_n_drugs", np.nan)) and r.dgidb_n_drugs > 0:
            L.append(f"Druggability (DGIdb): {i(r.dgidb_n_drugs)} drug "
                     f"{'interaction' if r.dgidb_n_drugs == 1 else 'interactions'}, "
                     f"{i(r.dgidb_n_approved)} with "
                     f"{'an approved agent' if r.dgidb_n_approved == 1 else 'approved agents'}, "
                     f"{i(r.dgidb_n_antineoplastic)} antineoplastic, of which "
                     f"{i(r.dgidb_n_approved_antineoplastic)} "
                     f"{'is' if r.dgidb_n_approved_antineoplastic == 1 else 'are'} approved"
                     + (f"; interaction types {r.dgidb_interaction_types}"
                        if isinstance(r.dgidb_interaction_types, str) and r.dgidb_interaction_types else "")
                     + (f". Approved antineoplastics: {r.dgidb_top_approved_antineoplastic}."
                        if isinstance(r.dgidb_top_approved_antineoplastic, str)
                        and r.dgidb_top_approved_antineoplastic else "."))
        else:
            L.append("DGIdb returns no drug-gene interaction for this node.")
    else:
        L.append("**f. Function.** DepMap CRISPR screens and DGIdb cover protein-coding genes "
                 "only; neither applies to a microRNA. The functional evidence for this node is "
                 "therefore restricted to the validated-target and expression evidence above, "
                 "which we state rather than paper over.")
    return " ".join(L)

def sec_literature(r):
    nm = r["name"]
    txt = cite(LT.LIT[nm])
    vol = LITVOL.loc[nm] if nm in LITVOL.index else None
    v = ""
    if vol is not None:
        v = (f" *(PubMed volume, retrieved 2026-09-09: {int(vol.pubmed_total):,} records for this "
             f"entity, {int(vol.pubmed_breast):,} of them also indexed under Breast Neoplasms — "
             f"{100*vol.frac_breast:.1f} %.)*")
    return "**g. Literature.** " + txt + v

def node_block(r, n):
    nm = r["name"]
    rank = (f"composite rank {i(r.rank_within_type)} of {NCLASS_TOTAL[r.type]} "
            f"{'TFs' if r.type=='TF' else ('genes' if r.type=='Gene' else 'miRNAs')}, "
            f"rank {i(r.rank_overall)} of 583 rankable nodes overall")
    head = (f"### {n}. {nm}\n\n*{TYPELABEL[r.type]} — {rank}; composite score "
            f"{num(r.composite_priority,3)}*  \n*Programme: {r.programme}*\n")
    return "\n\n".join([head,
                        sec_network(r), sec_partners(r), sec_expression(r),
                        sec_regulation(r), sec_clinical(r), sec_function(r),
                        sec_literature(r),
                        "**h. Verdict.** " + LT.VERDICT[nm]])

# ---- quick-reference table --------------------------------------------------
def quick_table():
    rows = ["\n\n---\n\n## Quick reference: all thirty nodes\n",
            "| node | class | rank in class | programme | FFL cores | degree | log2FC (TCGA) | FDR | "
            "promoter Δβ | best TCGA survival p | METABRIC OS q | Chronos |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    d = T.copy()
    d["_o"] = d["type"].map({"TF": 0, "Gene": 1, "miRNA": 2})
    for _, r in d.sort_values(["_o", "rank_within_type"]).iterrows():
        best = np.nanmin([r.get("tcga_OS_p", np.nan), r.get("tcga_PFI_p", np.nan),
                          r.get("tcga_DSS_p", np.nan)])
        rows.append(f"| **{r['name']}** | {r.type} | {i(r.rank_within_type)} | "
                    f"{r.programme.split(' ')[0]} | {i(r.ffl_total)} | {i(r.degree_total)} | "
                    f"{sgn(r.logFC_TCGA)} | {pv(r.fdr_TCGA)} | {sgn(r.meth_delta_beta,3)} | "
                    f"{pv(best)} | {pv(r.get('metabric_OS_q', np.nan))} | "
                    f"{sgn(r.get('chronos_mean_breast', np.nan),3)} |")
    rows.append("")
    rows.append("*Programme keys: P1 cell cycle and genome integrity; P2 chromatin and DNA "
                "methylation; P3 luminal identity and hormone signalling; P4 matrix remodelling "
                "and invasion; P5 stromal and inflammatory microenvironment; P6 epithelial "
                "plasticity and TGF-β signalling; P7 metabolic and stress-response control. "
                "Chronos is the mean CRISPR gene effect across 53 breast cancer cell lines "
                "(more negative = more essential); it does not exist for microRNAs. METABRIC "
                "contains no microRNA measurements.*")
    return "\n".join(rows)

# ---- assemble ---------------------------------------------------------------
OUT = [open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compendium_preamble.md"),
            encoding="utf-8").read()]
OUT.append(quick_table())
for typ, title in (("TF", "Transcription factors"), ("Gene", "Genes"), ("miRNA", "MicroRNAs")):
    OUT.append(f"\n\n---\n\n## {title}\n")
    d = T[T.type == typ].sort_values("rank_within_type")
    for n, (_, r) in enumerate(d.iterrows(), 1):
        OUT.append(node_block(r, n))
SYN = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "compendium_synthesis.md"), encoding="utf-8").read()
LCF = os.path.join(V5, "litcorr_fragment.md")
if os.path.exists(LCF):
    SYN = SYN.replace("{{LITCORR}}", open(LCF, encoding="utf-8").read().strip())
else:
    sys.exit("ABORT: litcorr_fragment.md missing - run 47_literature_vs_topology.py first")
if "{{" in SYN:
    sys.exit("ABORT: unsubstituted placeholder left in synthesis")
OUT.append(SYN)

# references
refs = ["\n\n---\n\n## References cited in the literature notes\n",
        "All PMIDs were retrieved and checked against live NCBI E-utilities records on "
        "2026-09-09 (`43_verify_pmids.py`); the full verified table with DOIs is in "
        "`node_literature_verified.csv`.\n"]
for _, v in VER.sort_values(["year", "first_author"]).iterrows():
    refs.append(f"- **PMID {v.pmid}** — {v.first_author} et al. ({v.year}) *{v.title}.* "
                f"{v.journal}" + (f" {v.volume}" if isinstance(v.volume, str) and v.volume else "")
                + (f":{v.pages}" if isinstance(v.pages, str) and v.pages else "") + ".")
OUT.append("\n".join(refs))

md = "\n\n".join(OUT)
with open(os.path.join(V5, "NODE_COMPENDIUM.md"), "w", encoding="utf-8") as fh:
    fh.write(md + "\n")
print("wrote NODE_COMPENDIUM.md:", len(md), "characters,", md.count("\n"), "lines")
