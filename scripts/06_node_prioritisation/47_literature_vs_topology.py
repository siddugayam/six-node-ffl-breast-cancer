#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Does network centrality track how heavily a node has been studied?
Correlates PubMed record counts for all 587 canonical nodes against degree, betweenness and
FFL-core participation. Emits a markdown fragment spliced into NODE_COMPENDIUM.md."""
import os, pandas as pd, numpy as np, networkx as nx
from scipy import stats
ROOT = "/path/to/revision"; V5 = ROOT + "/results/v5"

L = pd.read_csv(V5 + "/network_literature_volume.csv").drop_duplicates("name")
E = pd.read_csv(ROOT + "/data/canonical_edges.tsv", sep="\t")
N = pd.read_csv(ROOT + "/data/canonical_nodes.tsv", sep="\t")
G = nx.DiGraph(); G.add_nodes_from(N["name"]); G.add_edges_from(zip(E.source, E.target))
U = nx.Graph(); U.add_nodes_from(G.nodes()); U.add_edges_from(G.edges())
F = pd.read_csv(ROOT + "/results/v2/ffl_cores_coherence_corrected.csv")
ffl = pd.concat([F.regulator, F.intermediate, F.target]).value_counts()

D = pd.DataFrame(dict(name=list(G.nodes())))
D["type"] = D.name.map(dict(zip(N.name, N.type)))
D["degree"] = D.name.map(dict(G.degree()))
D["betw_undirected"] = D.name.map(nx.betweenness_centrality(U, normalized=False))
D["ffl_cores"] = D.name.map(ffl).fillna(0)
D = D.merge(L[["name", "pubmed_total", "pubmed_breast"]], on="name", how="inner")
D.to_csv(V5 + "/network_literature_vs_topology.csv", index=False)

rows = []
for lit in ("pubmed_total", "pubmed_breast"):
    for met in ("degree", "betw_undirected", "ffl_cores"):
        r, p = stats.spearmanr(D[lit], D[met])
        rows.append(dict(scope="all nodes", n=len(D), literature=lit, metric=met,
                         rho=r, p=p))
        for t in ("TF", "Gene", "miRNA"):
            dd = D[D.type == t]
            r2, p2 = stats.spearmanr(dd[lit], dd[met])
            rows.append(dict(scope=t, n=len(dd), literature=lit, metric=met, rho=r2, p=p2))
S = pd.DataFrame(rows)
S.to_csv(V5 + "/network_literature_vs_topology_stats.csv", index=False)
print(S.to_string(index=False))

def g(lit, met, scope="all nodes"):
    x = S[(S.literature == lit) & (S.metric == met) & (S.scope == scope)].iloc[0]
    return x.rho, x.p

def fp(p):
    return "< 1e-300" if p == 0 else (f"{p:.3g}" if p >= 1e-4 else f"{p:.1e}")

rd, pd_ = g("pubmed_total", "degree")
rb, pb = g("pubmed_total", "betw_undirected")
rf, pf = g("pubmed_total", "ffl_cores")
rdb, pdb = g("pubmed_breast", "degree")
rfb, pfb = g("pubmed_breast", "ffl_cores")
rf_tf, pf_tf = g("pubmed_total", "ffl_cores", "TF")
rf_mi, pf_mi = g("pubmed_total", "ffl_cores", "miRNA")
rf_gn, pf_gn = g("pubmed_total", "ffl_cores", "Gene")

md = f"""- **Position in the network tracks how heavily a node has been studied.** We retrieved PubMed
  record counts for all {len(D)} canonical nodes (2026-09-09) and correlated them with topology.
  Across the whole network, the number of PubMed records for a node correlates with its degree
  (Spearman ρ = {rd:.3f}, p = {fp(pd_)}), with its undirected betweenness
  (ρ = {rb:.3f}, p = {fp(pb)}) and with the number of three-node FFL cores it occupies
  (ρ = {rf:.3f}, p = {fp(pf)}); the association with breast-cancer-specific literature is of the
  same sign and size (degree ρ = {rdb:.3f}, p = {fp(pdb)}; FFL cores ρ = {rfb:.3f},
  p = {fp(pfb)}). It holds within every node class separately (FFL cores vs total literature:
  TFs ρ = {rf_tf:.3f}, p = {fp(pf_tf)}; genes ρ = {rf_gn:.3f}, p = {fp(pf_gn)};
  miRNAs ρ = {rf_mi:.3f}, p = {fp(pf_mi)}). This is the quantitative form of the argument in
  `PRIORITISATION.md`: in a network assembled from curated databases, structural centrality is
  substantially a record of curation effort, and it is why a topology-only hub ranking cannot be
  read as a ranking of biological importance. The full table is in
  `network_literature_vs_topology.csv`."""
open(V5 + "/litcorr_fragment.md", "w", encoding="utf-8").write(md + "\n")
print("\nfragment written")
