#!/usr/bin/env python3
"""ORA term overlap (SETTINGS.md, Changes and decisions, 2026-09-28 20:14); descriptive, part of no decision line.
For each of the twelve networks of one enrichment run (section E, or a rerun of A2 or A3), per universe and ontology:
  n_sig          significant terms here (p.adjust < 0.05, E3_ora_full.csv)
  paper          the paper's significant terms for the matching motif set (results/enrichment_all_motifs_FULL.csv,
                 column significant_BH_0.05): 3node_miRNA_FFL -> 3-miR, 3node_TF_FFL -> 3-TF, 3node_composite_FFL ->
                 3-Comp, four- to six-node networks -> 4-node, 5-node, 6-node (the paper has one set per size above 3)
  shared_paper   terms significant both here and in the paper's set (by term ID); Jaccard = shared / union
  sectionE       (reruns only) the same network's significant terms in section E, and the terms shared with them
KEGG is downloaded when enrichKEGG runs (release in E3_kegg_release.txt), so KEGG overlaps include release drift.
usage: python3 ora_overlap.py <analysis root> <run folder> [<section E folder>]
Output: <run folder>/E3_ora_overlap.csv"""
import os, sys
import pandas as pd

NETS = [f'{n}node_{c}' for n in (3, 4, 5, 6) for c in ('miRNA_FFL', 'TF_FFL', 'composite_FFL')]
PAPER = {'3node_miRNA_FFL': '3-miR', '3node_TF_FFL': '3-TF', '3node_composite_FFL': '3-Comp'}


def sig(F, key):
    return {k: set(g.ID) for k, g in F[F.significant_BH_0_05].groupby(key)}


def load(p):
    F = pd.read_csv(p, usecols=['set', 'universe', 'ontology', 'ID', 'significant_BH_0.05'])
    F = F.rename(columns={'significant_BH_0.05': 'significant_BH_0_05'})
    F['significant_BH_0_05'] = F['significant_BH_0_05'].astype(str).str.upper().eq('TRUE')
    return F


def main(root, run, efolder=None):
    F = load(os.path.join(run, 'E3_ora_full.csv'))
    P = pd.read_csv(os.path.join(root, 'results/enrichment_all_motifs_FULL.csv'), usecols=['motif_set', 'universe', 'ontology', 'ID', 'significant_BH_0.05'])
    P = P.rename(columns={'significant_BH_0.05': 'significant_BH_0_05'})
    P['significant_BH_0_05'] = P['significant_BH_0_05'].astype(str).str.upper().eq('TRUE')
    here, paper = sig(F, ['set', 'universe', 'ontology']), sig(P, ['motif_set', 'universe', 'ontology'])
    tested = F.groupby(['set', 'universe', 'ontology']).size()
    E = sig(load(os.path.join(efolder, 'E3_ora_full.csv')), ['set', 'universe', 'ontology']) if efolder else None
    rows = []
    for net in NETS:
        ps = PAPER.get(net, f'{net[0]}-node')
        for uni in ('TCGA_tested_plus_network', 'network_protein_coding'):
            for ont in ('GO_BP', 'GO_MF', 'GO_CC', 'KEGG', 'Reactome'):
                if (net, uni, ont) not in tested.index: continue
                h = here.get((net, uni, ont), set()); p = paper.get((ps, uni, ont), set())
                d = dict(set=net, universe=uni, ontology=ont, n_tested=int(tested[(net, uni, ont)]), n_sig=len(h), paper_motif_set=ps,
                         n_sig_paper=len(p), n_shared_paper=len(h & p), jaccard_paper=(len(h & p) / len(h | p)) if (h | p) else float('nan'))
                if E is not None:
                    e = E.get((net, uni, ont), set())
                    d.update(n_sig_sectionE=len(e), n_shared_sectionE=len(h & e), jaccard_sectionE=(len(h & e) / len(h | e)) if (h | e) else float('nan'))
                rows.append(d)
    R = pd.DataFrame(rows); R.to_csv(os.path.join(run, 'E3_ora_overlap.csv'), index=False)
    print(f'{os.path.basename(os.path.normpath(run))}: {len(R)} rows')


if __name__ == '__main__':
    main(*sys.argv[1:4])
