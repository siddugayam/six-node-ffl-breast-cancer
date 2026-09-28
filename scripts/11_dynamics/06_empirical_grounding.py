#!/usr/bin/env python3
"""
06_empirical_grounding.py
Tie every modelled topology to its count in the BRCA network, and record the two
structural facts that motivate the modelling choices:
  (i)  no core in the network is Alon's C1 (all-activating): the miRNA->target edge
       is repressive by construction, so C1 is unrealisable for a miRNA-mediated FFL;
  (ii) the composite class splits into a delayed-negative-feedback subclass
       (TF -> miRNA -| TF) and a mutual-repression / toggle subclass (TF -| miRNA -| TF).
Out: dynamics_empirical_grounding.csv, dynamics_composite_feedback_split.csv
"""
import pandas as pd, numpy as np

REV = "/path/to/revision"
RES = f"{REV}/results/v3"
D = pd.read_csv(f"{REV}/results/v2/ffl_cores_coherence_corrected.csv")
E = pd.read_csv(f"{REV}/data/canonical_edges.tsv", sep="\t")

assert len(D) == 1649, len(D)
res = D[D.coherence != "unresolved"]
coh = res.coherence.isin(["C1", "C2", "C3", "C4"])
print(f"cores {len(D)}; resolved {len(res)}; coherent {coh.sum()} ({100*coh.mean():.1f}%); "
      f"incoherent {(~coh).sum()} ({100*(~coh).mean():.1f}%)")
print("C1 cores anywhere in the corrected set:", int((D.coherence == 'C1').sum()))

rows = []
for cls, g in D.groupby("class"):
    r = g[g.coherence != "unresolved"]
    rows.append(dict(ffl_class=cls, n_cores=len(g), n_resolved=len(r),
                     n_coherent=int(r.coherence.isin(["C1", "C2", "C3", "C4"]).sum()),
                     pct_coherent=100 * r.coherence.isin(["C1", "C2", "C3", "C4"]).mean()
                     if len(r) else np.nan,
                     **{f"n_{k}": int((r.coherence == k).sum())
                        for k in ["C1", "C2", "C3", "C4", "I1", "I2", "I3", "I4"]}))
rows.append(dict(ffl_class="ALL", n_cores=len(D), n_resolved=len(res),
                 n_coherent=int(coh.sum()), pct_coherent=100 * coh.mean(),
                 **{f"n_{k}": int((res.coherence == k).sum())
                    for k in ["C1", "C2", "C3", "C4", "I1", "I2", "I3", "I4"]}))
G = pd.DataFrame(rows)
G.to_csv(f"{RES}/dynamics_empirical_grounding.csv", index=False)
print(G.to_string(index=False))

# ---- composite feedback split ------------------------------------------------
C = D[D["class"] == "Composite-FFL"]
Cr = C[C.coherence != "unresolved"]
split = pd.DataFrame([
    dict(subclass="TF -> miRNA -| TF  (delayed NEGATIVE FEEDBACK inside the FFL)",
         sign_TF_to_miRNA=+1, embedded_loop="negative feedback (activation-repression)",
         n_all=int((C.sign_RM == 1).sum()), n_resolved=int((Cr.sign_RM == 1).sum()),
         pct_of_resolved_composite=100 * (Cr.sign_RM == 1).mean(),
         coherence_types="; ".join(f"{k}={int((Cr[Cr.sign_RM==1].coherence==k).sum())}"
                                   for k in ["C4", "I1"]),
         Alon_counterpart="none in the 3-node transcriptional catalogue"),
    dict(subclass="TF -| miRNA -| TF  (MUTUAL REPRESSION / toggle inside the FFL)",
         sign_TF_to_miRNA=-1, embedded_loop="double-negative (positive) feedback",
         n_all=int((C.sign_RM == -1).sum()), n_resolved=int((Cr.sign_RM == -1).sum()),
         pct_of_resolved_composite=100 * (Cr.sign_RM == -1).mean(),
         coherence_types="; ".join(f"{k}={int((Cr[Cr.sign_RM==-1].coherence==k).sum())}"
                                   for k in ["C2", "I3"]),
         Alon_counterpart="none in the 3-node transcriptional catalogue"),
])
split.to_csv(f"{RES}/dynamics_composite_feedback_split.csv", index=False)
print()
print(split.to_string(index=False))

# ---- sanity: every miRNA->target / miRNA->TF edge is repressive --------------
m = E[E.edge_type.isin(["miRNA_target"])]
print(f"\nmiRNA-outgoing edges: {len(m)}, all sign -1: {bool((m.sign == -1).all())}")
print(f"TF_target edges: {(E.edge_type=='TF_target').sum()}, all sign +1 in canonical file: "
      f"{bool((E[E.edge_type=='TF_target'].sign==1).all())}")
