#!/usr/bin/env python3
"""scripts/06_node_prioritisation/47_literature_vs_topology.py statistics (PubMed volume vs degree / undirected betweenness /
FFL cores, Spearman, all nodes) recomputed with and without the 30 legacy miRNA-miRNA edges.
Same code for the graph and the correlations; outputs only here."""
import pandas as pd, networkx as nx
from scipy import stats
ROOT = "/path/to/revision"; V5 = ROOT + "/results/v5"
L = pd.read_csv(V5 + "/network_literature_volume.csv").drop_duplicates("name")
E0 = pd.read_csv(ROOT + "/data/canonical_edges.tsv", sep="\t"); N = pd.read_csv(ROOT + "/data/canonical_nodes.tsv", sep="\t")
F = pd.read_csv(ROOT + "/results/v2/ffl_cores_coherence_corrected.csv")
ffl = pd.concat([F.regulator, F.intermediate, F.target]).value_counts()
stored = pd.read_csv(V5 + "/network_literature_vs_topology_stats.csv")
rows = []
for mode, E in (("full", E0), ("nolegacy", E0[E0.edge_type != "miRNA_miRNA"])):
    G = nx.DiGraph(); G.add_nodes_from(N["name"]); G.add_edges_from(zip(E.source, E.target))
    U = nx.Graph(); U.add_nodes_from(G.nodes()); U.add_edges_from(G.edges())
    D = pd.DataFrame(dict(name=list(G.nodes())))
    D["degree"] = D.name.map(dict(G.degree())); D["betw_undirected"] = D.name.map(nx.betweenness_centrality(U, normalized=False))
    D["ffl_cores"] = D.name.map(ffl).fillna(0)
    D = D.merge(L[["name", "pubmed_total", "pubmed_breast"]], on="name", how="inner")
    for lit in ("pubmed_total", "pubmed_breast"):
        for met in ("degree", "betw_undirected", "ffl_cores"):
            r, p = stats.spearmanr(D[lit], D[met]); rows.append(dict(network=mode, n=len(D), literature=lit, metric=met, rho=round(r, 4), p=p))
R = pd.DataFrame(rows); R.to_csv("literature_vs_topology_both.csv", index=False); print(R.to_string(index=False))
print("\nstored (all nodes):"); print(stored[stored.scope == "all nodes"].to_string(index=False))
