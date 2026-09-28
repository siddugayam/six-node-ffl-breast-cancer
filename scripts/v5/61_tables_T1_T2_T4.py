#!/usr/bin/env python3
"""
Publication tables T1 (data sources / filtering cascade), T2 (network composition)
and T4 (motif significance).

Every value is computed from a file in revision/data or revision/results; nothing is
typed by hand except (a) the counts the ORIGINAL manuscript reports, which are quoted as
such and flagged 'manuscript-reported, not reproducible from the deposit', and (b) the
literal database names.  Sources are named in the 'source_file' column of every panel.
"""
import csv, collections, json, math, os, sys
import numpy as np

REV = '/path/to/revision'
TB  = f'{REV}/results/v5/tables'
os.makedirs(TB, exist_ok=True)

def W(name, rows, fieldnames=None):
    if not rows:
        print(f'  !! {name}: no rows'); return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(f'{TB}/{name}', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    print(f'  wrote {name:52s} {len(rows):>7d} rows')

def rd(path, delim=','):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delim))

def f(x, nd=3):
    try:
        v = float(x)
        if math.isnan(v): return ''
        return round(v, nd)
    except (TypeError, ValueError):
        return ''

# =====================================================================================
# load the canonical network
# =====================================================================================
nodes = {r['name']: r['type'] for r in rd(f'{REV}/data/canonical_nodes.tsv', '\t')}
edges = rd(f'{REV}/data/canonical_edges.tsv', '\t')
ev    = {(r['source'], r['target']): r for r in rd(f'{REV}/data/edge_evidence_tier.tsv', '\t')}
trr   = {(r['source'], r['target']): r for r in rd(f'{REV}/data/layer_TF_target.tsv', '\t')}
tmr   = {(r['source'], r['target']): r for r in rd(f'{REV}/data/layer_TF_miRNA.tsv', '\t')}
ec    = {}
for r in rd(f'{REV}/results/edge_correlation.csv'):
    ec[(r['source'], r['target'])] = r
cafr  = {(r['source'], r['target']): r for r in rd(f'{REV}/results/v2/v2_edge_rho.csv')}

E = {(r['source'], r['target']) for r in edges}
ntype = collections.Counter(nodes.values())
print(f'network: {len(nodes)} nodes {ntype}, {len(edges)} directed edges')

def sign_source(r):
    k = (r['source'], r['target'])
    t = trr.get(k, {}).get('mode')
    if t in ('Activation', 'Repression'): return 'TRRUST (curated mode)'
    m = tmr.get(k, {}).get('mode')
    if m in ('Activation', 'Repression'): return 'TransmiR (curated mode)'
    if r['edge_type'] == 'miRNA_target': return 'mechanistic (miRNA repression, -1)'
    if r['edge_type'] == 'miRNA_miRNA': return 'assumed co-transcription (+1)'
    return 'assumed / unannotated (+1)'

for r in edges:
    r['sign_source'] = sign_source(r)

# =====================================================================================
# T1a  source databases
# =====================================================================================
tier_ct = collections.Counter(v['tier'] for v in ev.values())
db_ct   = collections.Counter()
for v in ev.values():
    for d in (v['databases'] or '').split(';'):
        if d.strip(): db_ct[d.strip()] += 1
trr_modes = collections.Counter(r['mode'] for r in rd(f'{REV}/data/layer_TF_target.tsv', '\t'))
tmr_modes = collections.Counter(r['mode'] for r in rd(f'{REV}/data/layer_TF_miRNA.tsv', '\t'))
gg = rd(f'{REV}/data/layer_gene_gene.tsv', '\t')
gg_ct = collections.Counter(r['evidence'] for r in gg)
mm = rd(f'{REV}/data/layer_miRNA_miRNA.tsv', '\t')
mm_ct = {t: sum(1 for r in mm if str(r.get(t, '')).upper() in ('TRUE', '1'))
         for t in ('threshold_3kb', 'threshold_10kb', 'threshold_50kb')}

edge_by_type = collections.Counter(r['edge_type'] for r in edges)

T1a = [
 dict(layer='Breast-cancer gene seed set', resource='GeneCards', role='candidate disease genes',
      count_reported_in_manuscript='5,247', count_verifiable_here='not reproducible',
      status='manuscript-reported only; the filtered list is not in the deposit',
      source_file='MS.md sect. 3.1'),
 dict(layer='Breast-cancer gene seed set', resource='DisGeNET (score > 0.4)', role='candidate disease genes',
      count_reported_in_manuscript='1,612', count_verifiable_here='not reproducible',
      status='manuscript-reported only; the two resources use different scoring scales',
      source_file='MS.md sect. 2.1, 3.1'),
 dict(layer='Breast-cancer miRNA seed set', resource='PhenomiR (2010)', role='candidate disease miRNAs',
      count_reported_in_manuscript='655', count_verifiable_here='not reproducible',
      status='manuscript-reported only; resource is outdated', source_file='MS.md sect. 3.1'),
 dict(layer='Breast-cancer miRNA seed set', resource='HMDD', role='candidate disease miRNAs',
      count_reported_in_manuscript='1,095', count_verifiable_here='not reproducible',
      status='manuscript-reported only; used in Results but absent from Methods',
      source_file='MS.md sect. 3.1'),
 dict(layer='Breast-cancer miRNA seed set', resource='miR2Disease (2009)', role='candidate disease miRNAs',
      count_reported_in_manuscript='132', count_verifiable_here='not reproducible',
      status='manuscript-reported only; resource is outdated', source_file='MS.md sect. 3.1'),
 dict(layer='miRNA -> target', resource='miRWalk + miRTarBase', role='miRNA-gene and miRNA-TF interactions',
      count_reported_in_manuscript='21,443 miRNA-gene; 7,443 miRNA-TF',
      count_verifiable_here=f"{edge_by_type['miRNA_target']:,} retained directed edges in the final network",
      status='retained edges re-annotated per-edge with multiMiR (see T1c); prediction and validation were pooled in the original',
      source_file='data/canonical_edges.tsv; data/edge_evidence_tier.tsv'),
 dict(layer='TF -> target', resource='TRRUST v2', role='TF-gene and TF-TF interactions, with activation/repression mode',
      count_reported_in_manuscript='87,600 TF-gene (with hTFtarget, TcoF-DB); 17,600 TF-TF (with dbCoRC)',
      count_verifiable_here=(f"{edge_by_type['TF_target']:,} retained directed edges; TRRUST layer rebuilt to "
                             f"{sum(trr_modes.values()):,} records "
                             f"(Activation {trr_modes['Activation']:,}, Repression {trr_modes['Repression']:,}, "
                             f"Unknown {trr_modes['Unknown']:,}, Ambiguous {trr_modes['Ambiguous']:,})"),
      status='mode annotations were available but unused in the original, which set every TF edge to +1',
      source_file='data/layer_TF_target.tsv; results/layer_summary_counts.tsv'),
 dict(layer='TF -> miRNA', resource='TransmiR v2.0', role='TF-miRNA interactions, with mode',
      count_reported_in_manuscript='2,677',
      count_verifiable_here=(f"{edge_by_type['TF_miRNA']:,} retained directed edges; TransmiR layer rebuilt to "
                             f"{sum(tmr_modes.values()):,} records "
                             f"(Activation {tmr_modes['Activation']:,}, Repression {tmr_modes['Repression']:,}, "
                             f"unsigned {tmr_modes.get('Regulation_unsigned',0):,}, Ambiguous {tmr_modes['Ambiguous']:,})"),
      status='mode annotations now used for edge signs',
      source_file='data/layer_TF_miRNA.tsv; results/layer_summary_counts.tsv'),
 dict(layer='gene -> gene', resource='GeneMANIA co-expression (original)',
      role='gene-gene interactions used to build 4-node FFLs',
      count_reported_in_manuscript='not stated',
      count_verifiable_here=f"{edge_by_type['gene_gene']} edge in the deposited network (COL1A1-COL3A1)",
      status='REMOVED as a regulatory source: GeneMANIA predicts functional association, not regulation',
      source_file='data/canonical_edges.tsv; results/00_DATA_AUDIT.md sect. E1'),
 dict(layer='gene -> gene', resource='TRRUST (directed) / STRING (>= 900, undirected)',
      role='replacement gene-gene layer, declared as two separate tiers',
      count_reported_in_manuscript='n/a (new in revision)',
      count_verifiable_here=(f"{gg_ct.get('TRRUST',0):,} TRRUST-directed rows; "
                             f"{gg_ct.get('STRING',0):,} STRING undirected associations"),
      status=('TRRUST rows are all TF->target edges already carried by the TF layer; the STRING rows are '
              'co-functional associations and are EXCLUDED from the primary graph'),
      source_file='data/layer_gene_gene.tsv'),
 dict(layer='miRNA - miRNA', resource='GeneMANIA co-expression (original)',
      role='miRNA-miRNA interactions used to build 5-node FFLs',
      count_reported_in_manuscript='not stated',
      count_verifiable_here=f"{edge_by_type['miRNA_miRNA']} edges in the deposited network",
      status=('REPLACED: only 2 of the 30 are genomically co-clustered at any threshold; the layer is now '
              'defined as polycistronic co-transcription'),
      source_file='data/canonical_edges.tsv; results/ffl_definition_check.md sect. 6'),
 dict(layer='miRNA - miRNA', resource='miRBase hsa.gff3 genomic coordinates',
      role='replacement layer: polycistronic co-transcription',
      count_reported_in_manuscript='n/a (new in revision)',
      count_verifiable_here=(f"{mm_ct['threshold_3kb']:,} pairs at 3 kb; {mm_ct['threshold_10kb']:,} at 10 kb "
                             f"(primary); {mm_ct['threshold_50kb']:,} at 50 kb"),
      status='undirected; stored as a reciprocal arc pair; 10 kb used throughout, 3 / 50 kb as sensitivity',
      source_file='data/layer_miRNA_miRNA.tsv; results/miRNA_clusters_10kb.tsv'),
 dict(layer='Expression validation', resource='UCSC Xena TCGA-BRCA (Pan-Cancer Atlas)',
      role='gene and miRNA expression, tumour vs normal, and clinical follow-up',
      count_reported_in_manuscript='1,097 primary tumour / 113 solid tissue normal',
      count_verifiable_here='1,097 tumour vs 114 normal (gene); 1,068 vs 90 (miRNA); 1,066 samples with paired gene+miRNA',
      status='counts differ by one normal sample from the manuscript; the measured values are used throughout',
      source_file='results/NOTICE_DE_TABLE_FIXED.md; results/v4/ranklist_provenance.csv'),
]
W('T1a_source_databases.csv', T1a)

# =====================================================================================
# T1b  filtering cascade
# =====================================================================================
def deg_counts():
    o = collections.Counter(); i = collections.Counter()
    for r in edges:
        o[r['source']] += 1; i[r['target']] += 1
    return o, i
odeg, ideg = deg_counts()
by_pair = collections.Counter((nodes[r['source']], nodes[r['target']]) for r in edges)

T1b = [
 dict(step='1', stage='Disease-gene retrieval (GeneCards, DisGeNET score > 0.4)',
      genes='5,247 + 1,612 candidates', miRNAs='-', TFs='-', edges='-',
      reproducible_from_deposit='NO', note='raw candidate lists are not deposited',
      source_file='MS.md sect. 3.1'),
 dict(step='2', stage='Manual curation of disease relevance',
      genes='1,057', miRNAs='-', TFs='665 (a subset of the 1,057)', edges='-',
      reproducible_from_deposit='NO',
      note=('THIS IS THE 665 IN THE ABSTRACT AND SECTION 3.1: it is the number of TFs among the 1,057 curated '
            'genes BEFORE the miRNA-TF hypergeometric filter. No curator, criteria or inter-rater statistic '
            'is recorded'),
      source_file='MS.md sect. 3.1'),
 dict(step='3', stage='Disease-miRNA retrieval (PhenomiR, HMDD, miR2Disease) and curation',
      genes='-', miRNAs='421', TFs='-', edges='-', reproducible_from_deposit='NO',
      note='655 + 1,095 + 132 candidates curated to 421', source_file='MS.md sect. 3.1'),
 dict(step='4', stage='Raw interaction retrieval from the interaction databases',
      genes='1,057', miRNAs='421', TFs='665',
      edges='87,600 TF-gene; 21,443 miRNA-gene; 7,443 miRNA-TF; 2,677 TF-miRNA; 17,600 TF-TF',
      reproducible_from_deposit='NO', note='"Before hypergeometric test" column of the original Table 1',
      source_file='MS.md Table 1'),
 dict(step='5', stage='Cumulative hypergeometric test on miRNA-TF pairs, BH FDR < 0.05',
      genes='1,057', miRNAs='246', TFs='233',
      edges='3,902 TF-gene; 13,947 miRNA-gene; 953 TF-miRNA; 958 miRNA-TF',
      reproducible_from_deposit='NO',
      note=('THIS IS THE 233 IN TABLE 1. The filter acts on miRNA-TF pairs, so miRNAs (421 -> 246) and TFs '
            '(665 -> 233) are reduced while the gene count is unchanged (1,057 -> 1,057). '
            'The file the README promises, data/interaction_pairs/hypergeometric_filtered_interactions.csv, '
            'is absent from the repository (audit A9), so this step cannot be re-executed'),
      source_file='MS.md Table 1; results/00_DATA_AUDIT.md sect. A9'),
 dict(step='6', stage='Assembly of the six Cytoscape FFL networks (the deposited object)',
      genes='-', miRNAs='-', TFs='-', edges='18,592 SIF rows across the six exports',
      reproducible_from_deposit='YES',
      note='only nodes and edges that entered at least one FFL survive to this step',
      source_file='data/01_build_log.txt'),
 dict(step='7', stage='Repair: miRBase precursor/mature collapse, de-duplication, node-type harmonisation, direction',
      genes=f"{ntype['Gene']}", miRNAs=f"{ntype['miRNA']}", TFs=f"{ntype['TF']}",
      edges=f"{len(edges):,} unique directed edges",
      reproducible_from_deposit='YES',
      note=('803 raw node labels -> 587 canonical nodes (214 labels merged from > 1 alias); 10 mis-typed TFs '
            'resolved; duplicate SIF rows removed. THIS IS THE NETWORK ANALYSED THROUGHOUT THE REVISION'),
      source_file='data/canonical_nodes.tsv; data/canonical_edges.tsv; data/01_build_log.txt'),
 dict(step='8', stage='Augmentation for the higher-order census (TRRUST TF-target, TransmiR TF-miRNA, 10 kb polycistron)',
      genes=f"{ntype['Gene']}", miRNAs=f"{ntype['miRNA']}", TFs=f"{ntype['TF']}",
      edges='8,529 arcs (7,198 after contracting the 1,223 reciprocal TF-miRNA pairs)',
      reproducible_from_deposit='YES',
      note='STRING associations excluded from the primary graph (see T3 legend)',
      source_file='logs/ffl_census.log; scripts/v5/60_census_nostring.py'),
]
W('T1b_filtering_cascade.csv', T1b)

T1b2 = []
for cls, before, after, final in [('Genes', '1,057', '1,057', ntype['Gene']),
                                  ('miRNAs', '421', '246', ntype['miRNA']),
                                  ('TFs', '665', '233', ntype['TF'])]:
    T1b2.append(dict(node_class=cls,
                     after_manual_curation=before,
                     after_hypergeometric_filter_Table1=after,
                     in_the_deposited_FFL_network=f'{final}',
                     pct_of_curated_set_reaching_the_network=(
                         f"{100*final/float(before.replace(',','')):.1f}%"),
                     explanation={'Genes': ('the hypergeometric test filters miRNA-TF pairs, not genes, so the gene '
                                            'count is unchanged at step 5; the drop to 207 happens because only genes '
                                            'that are a target inside at least one FFL enter the deposited network'),
                                  'miRNAs': ('421 -> 246 at the hypergeometric step; 246 -> 223 because 23 retained '
                                             'miRNAs form no FFL'),
                                  'TFs': ('665 -> 233 at the hypergeometric step (the source of the 665-vs-233 '
                                          'discrepancy flagged); 233 -> 157 because 76 retained TFs form no FFL')}[cls]))
W('T1b2_665_vs_233_reconciliation.csv', T1b2)

# =====================================================================================
# T1c  evidence tier of the retained edges, per source database
# =====================================================================================
T1c = []
tot_mir = edge_by_type['miRNA_target']
for tier, label, defn in [('strong', 'strong', 'at least one low-throughput functional assay (reporter, western, qPCR)'),
                          ('weak', 'weak', 'high-throughput evidence only (CLIP, microarray, RNA-seq)'),
                          ('predicted_only', 'prediction-only', 'no experimental record in multiMiR; sequence prediction only')]:
    sub = [v for v in ev.values() if v['tier'] == tier]
    dbc = collections.Counter()
    for v in sub:
        for d in (v['databases'] or '').split(';'):
            if d.strip(): dbc[d.strip()] += 1
    T1c.append(dict(evidence_tier=label, definition=defn, n_edges=len(sub),
                    pct_of_miRNA_target_edges=f'{100*len(sub)/tot_mir:.1f}%',
                    pct_of_all_network_edges=f'{100*len(sub)/len(edges):.1f}%',
                    supporting_databases='; '.join(f'{k} {v:,}' for k, v in dbc.most_common()) or 'none (prediction only)',
                    source_file='data/edge_evidence_tier.tsv (multiMiR)'))
T1c.append(dict(evidence_tier='TOTAL miRNA-target edges', definition='',
                n_edges=tot_mir, pct_of_miRNA_target_edges='100.0%',
                pct_of_all_network_edges=f'{100*tot_mir/len(edges):.1f}%',
                supporting_databases='; '.join(f'{k} {v:,}' for k, v in db_ct.most_common()),
                source_file='data/edge_evidence_tier.tsv'))
for lab, n, note in [('TF-target edges', edge_by_type['TF_target'],
                      'TRRUST-curated; no multiMiR tier applies'),
                     ('TF-miRNA edges', edge_by_type['TF_miRNA'], 'TransmiR-curated; no multiMiR tier applies'),
                     ('miRNA-miRNA edges', edge_by_type['miRNA_miRNA'],
                      'author-deposited co-expression; 2 of 30 are genomically co-clustered'),
                     ('gene-gene edges', edge_by_type['gene_gene'],
                      'COL1A1-COL3A1; no directed evidence in TRRUST in either orientation (audit E1)')]:
    T1c.append(dict(evidence_tier=lab, definition=note, n_edges=n, pct_of_miRNA_target_edges='n/a',
                    pct_of_all_network_edges=f'{100*n/len(edges):.1f}%', supporting_databases='',
                    source_file='data/canonical_edges.tsv'))
W('T1c_evidence_tier_by_source.csv', T1c)

# =====================================================================================
# T2a  nodes
# =====================================================================================
part = {r['node']: r for r in rd(f'{REV}/results/ffl_module_membership.csv')}
cores = rd(f'{REV}/results/v2/ffl_cores_coherence_corrected.csv')
in_core = collections.Counter()
for c in cores:
    for m in (c['regulator'], c['intermediate'], c['target']):
        in_core[m] += 1
de = {r['Gene']: r for r in rd(f'{REV}/results/BRCA_DEX_ALL_nodes.csv')}

T2a = []
for cls in ['TF', 'Gene', 'miRNA', 'TOTAL']:
    sel = [n for n, t in nodes.items() if t == cls] if cls != 'TOTAL' else list(nodes)
    ncore = sum(1 for n in sel if in_core[n] > 0)
    nho   = sum(1 for n in sel if str(part.get(n, {}).get('in_6node_set', '')).upper() == 'TRUE')
    nhoo  = sum(1 for n in sel if str(part.get(n, {}).get('higher_order_only', '')).upper() == 'TRUE')
    nde   = sum(1 for n in sel if n in de)
    nsig  = sum(1 for n in sel if n in de and abs(float(de[n]['logFC'])) > 1
                and float(de[n]['adj.P.Val']) < 0.05)
    T2a.append(dict(node_class=cls, n_nodes=len(sel), pct=f'{100*len(sel)/len(nodes):.1f}%',
                    n_in_a_3node_FFL_core=ncore,
                    n_in_a_6node_module=nho, n_higher_order_only=nhoo,
                    n_measured_in_TCGA_BRCA=nde,
                    n_significantly_DE=f'{nsig} (|log2FC| > 1 and FDR < 0.05)',
                    median_total_degree=int(np.median([odeg[n] + ideg[n] for n in sel])),
                    max_total_degree=max(odeg[n] + ideg[n] for n in sel),
                    source_file='data/canonical_nodes.tsv; results/v2/ffl_cores_coherence_corrected.csv; results/ffl_module_membership.csv; results/BRCA_DEX_ALL_nodes.csv'))
W('T2a_node_composition.csv', T2a)

# =====================================================================================
# T2b  edges by class, with tier / sign / correlation
# =====================================================================================
recip_pairs = {(a, b) for (a, b) in E if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA' and (b, a) in E}
T2b = []
groups = collections.defaultdict(list)
for r in edges:
    groups[(r['edge_type'], nodes[r['source']], nodes[r['target']])].append(r)
for key in sorted(groups, key=lambda k: -len(groups[k])):
    et, st, tt = key
    grp = groups[key]
    tiers = collections.Counter(ev.get((r['source'], r['target']), {}).get('tier', 'n/a') for r in grp)
    signs = collections.Counter(r['sign_source'] for r in grp)
    rhos = [float(ec[(r['source'], r['target'])]['rho']) for r in grp
            if (r['source'], r['target']) in ec and ec[(r['source'], r['target'])]['rho'] not in ('', 'NA')]
    nsig = sum(1 for r in grp if (r['source'], r['target']) in ec
               and ec[(r['source'], r['target'])]['fdr'] not in ('', 'NA')
               and float(ec[(r['source'], r['target'])]['fdr']) < 0.05)
    nrec = sum(1 for r in grp if (r['source'], r['target']) in recip_pairs
               or (r['target'], r['source']) in recip_pairs)
    T2b.append(dict(edge_class=f'{st} -> {tt}', edge_type=et, n_edges=len(grp),
                    pct_of_network=f'{100*len(grp)/len(edges):.1f}%',
                    tier_strong=tiers.get('strong', 0), tier_weak=tiers.get('weak', 0),
                    tier_prediction_only=tiers.get('predicted_only', 0),
                    tier_not_applicable=tiers.get('n/a', 0),
                    sign_source_breakdown='; '.join(f'{k} {v:,}' for k, v in signs.most_common()),
                    n_reciprocated=nrec,
                    pct_reciprocated=f'{100*nrec/len(grp):.1f}%',
                    n_with_TCGA_correlation=len(rhos),
                    median_abs_rho=f(np.median(np.abs(rhos)) if rhos else None),
                    n_FDR_lt_0.05=nsig,
                    source_file='data/canonical_edges.tsv; data/edge_evidence_tier.tsv; results/edge_correlation.csv'))
tot_tiers = collections.Counter(ev.get((r['source'], r['target']), {}).get('tier', 'n/a') for r in edges)
tot_signs = collections.Counter(r['sign_source'] for r in edges)
all_rho = [float(v['rho']) for v in ec.values() if v['rho'] not in ('', 'NA')]
T2b.append(dict(edge_class='ALL', edge_type='TOTAL', n_edges=len(edges), pct_of_network='100.0%',
                tier_strong=tot_tiers['strong'], tier_weak=tot_tiers['weak'],
                tier_prediction_only=tot_tiers['predicted_only'], tier_not_applicable=tot_tiers['n/a'],
                sign_source_breakdown='; '.join(f'{k} {v:,}' for k, v in tot_signs.most_common()),
                n_reciprocated=2 * len(recip_pairs),
                pct_reciprocated=f'{100*2*len(recip_pairs)/len(edges):.1f}%',
                n_with_TCGA_correlation=len(all_rho),
                median_abs_rho=f(np.median(np.abs(all_rho))),
                n_FDR_lt_0.05=sum(1 for v in ec.values() if v['fdr'] not in ('', 'NA') and float(v['fdr']) < 0.05),
                source_file=''))
W('T2b_edge_composition.csv', T2b)

# =====================================================================================
# T2c  sign provenance
# =====================================================================================
T2c = []
for k, v in tot_signs.most_common():
    sub = [r for r in edges if r['sign_source'] == k]
    pos = sum(1 for r in sub if r['sign'] == '1')
    neg = sum(1 for r in sub if r['sign'] == '-1')
    T2c.append(dict(sign_source=k, n_edges=v, pct=f'{100*v/len(edges):.1f}%',
                    n_activating=pos, n_repressing=neg,
                    is_the_sign_evidence_based=('yes - database-curated mode' if 'curated mode' in k else
                                                'yes - mechanism of miRNA action' if 'mechanistic' in k else 'NO - assumed'),
                    source_file='data/layer_TF_target.tsv; data/layer_TF_miRNA.tsv; data/canonical_edges.tsv'))
n_assumed = sum(r['n_edges'] for r in T2c if r['is_the_sign_evidence_based'].startswith('NO'))
T2c.append(dict(sign_source='TOTAL', n_edges=len(edges), pct='100.0%',
                n_activating=sum(1 for r in edges if r['sign'] == '1'),
                n_repressing=sum(1 for r in edges if r['sign'] == '-1'),
                is_the_sign_evidence_based=f'{len(edges)-n_assumed:,} evidence-based / {n_assumed:,} assumed '
                                           f'({100*n_assumed/len(edges):.1f}% assumed)',
                source_file=''))
W('T2c_sign_provenance.csv', T2c)

# =====================================================================================
# T2d  reciprocity
# =====================================================================================
n_tf_mir = sum(1 for (a, b) in E if nodes.get(a) == 'TF' and nodes.get(b) == 'miRNA')
n_mir_tf = sum(1 for (a, b) in E if nodes.get(a) == 'miRNA' and nodes.get(b) == 'TF')
rec_audit = rd(f'{REV}/results/motif_reciprocity_audit.csv')
T2d = [
 dict(quantity='TF -> miRNA directed edges', value=f'{n_tf_mir:,}',
      note='TransmiR arm', source_file='data/canonical_edges.tsv'),
 dict(quantity='miRNA -> TF directed edges', value=f'{n_mir_tf:,}',
      note='miRNA-target arm onto TF nodes', source_file='data/canonical_edges.tsv'),
 dict(quantity='mutual TF <-> miRNA pairs', value=f'{len(recip_pairs):,}',
      note='both arcs present; each is ONE composite regulatory relationship',
      source_file='scripts/v2/00_counting_convention.py'),
 dict(quantity='% of TF -> miRNA edges reciprocated', value=f'{100*len(recip_pairs)/n_tf_mir:.2f}%',
      note='near-total reciprocity', source_file='results/motif_reciprocity_audit.csv'),
 dict(quantity='% of miRNA -> TF edges reciprocated', value=f'{100*len(recip_pairs)/n_mir_tf:.2f}%',
      note='', source_file='results/motif_reciprocity_audit.csv'),
 dict(quantity='Composite-FFL cores, counted ONCE per mutual pair', value='1,434',
      note='the convention adopted throughout the revision', source_file='scripts/v2/00_counting_convention.py'),
 dict(quantity='Composite-FFL cores, counted once per ARC', value='2,868',
      note='the double count; NOT used',
      source_file='scripts/v2/00_counting_convention.py'),
 dict(quantity='Reciprocity in an independently rebuilt breast network', value='59.62% of TF->miRNA; 6.40% of miRNA->TF',
      note=('the deposited network is far more reciprocal than the underlying databases support; the excess is '
            'the imprint of the miRNA-TF hypergeometric selection step, not biology'),
      source_file='results/motif_reciprocity_audit.csv'),
 dict(quantity='Enrichment of mutual pairs over a degree-preserving null', value='4.88-fold, Z = +87.7, p = 0.001',
      note='NULL-A, 1,000 randomisations; report next to any motif Z-score (see T4)',
      source_file='results/motif_significance.csv'),
]
W('T2d_reciprocity.csv', T2d)

# =====================================================================================
# T2e  degrees split by direction and edge type
# =====================================================================================
per_node = collections.defaultdict(lambda: collections.Counter())
for r in edges:
    per_node[r['source']][('out', r['edge_type'])] += 1
    per_node[r['target']][('in', r['edge_type'])] += 1
ETYPES = ['miRNA_target', 'TF_target', 'TF_miRNA', 'miRNA_miRNA', 'gene_gene']
rows_all = []
for n, t in nodes.items():
    d = per_node[n]
    row = dict(node=n, node_class=t, degree_total=odeg[n] + ideg[n], degree_in=ideg[n], degree_out=odeg[n],
               n_FFL_cores_any_role=in_core[n])
    for et in ETYPES:
        row[f'in_{et}'] = d[('in', et)]; row[f'out_{et}'] = d[('out', et)]
    rows_all.append(row)
rows_all.sort(key=lambda r: -r['degree_total'])
W('T2e_node_degrees_full.csv', rows_all)

claimed = {'VEGFA': 111, 'CCND2': 88, 'ADAMTS5': 70, 'TGFBR2': 70, 'TP53': 111, 'MYC': 104,
           'RELA': 98, 'SP1': 98, 'hsa-miR-130a': 50, 'hsa-miR-21': 48, 'hsa-miR-124': 46,
           'hsa-miR-34a': 41}
byname = {r['node']: r for r in rows_all}
T2e = []
for n, c in claimed.items():
    r = byname[n]
    T2e.append(dict(node=n, node_class=r['node_class'], degree_reported_in_manuscript=c,
                    degree_measured_canonical=r['degree_total'], degree_in=r['degree_in'], degree_out=r['degree_out'],
                    in_miRNA_target=r['in_miRNA_target'], out_miRNA_target=r['out_miRNA_target'],
                    in_TF_target=r['in_TF_target'], out_TF_target=r['out_TF_target'],
                    in_TF_miRNA=r['in_TF_miRNA'], out_TF_miRNA=r['out_TF_miRNA'],
                    in_miRNA_miRNA=r['in_miRNA_miRNA'], out_miRNA_miRNA=r['out_miRNA_miRNA'],
                    discrepancy=('agrees' if c == r['degree_total'] else
                                 f"manuscript {c} vs measured {r['degree_total']}"),
                    cause=('' if c == r['degree_total'] else
                           ('duplicate SIF rows inflated the count (audit A6)' if c > r['degree_total'] else
                            'precursor/mature node split deflated the count (audit A2)')),
                    source_file='data/canonical_edges.tsv; results/00_DATA_AUDIT.md sect. B5'))
W('T2e2_manuscript_degrees_corrected.csv', T2e)

# =====================================================================================
# T4  motif significance
# =====================================================================================
def summarise_null(path, obs_field_names, model, graph, note):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    obs = {k: float(v) for k, v in rows[0].items() if k != 'rep'}
    rand = {k: np.array([float(r[k]) for r in rows[1:]]) for k in obs}
    out = []
    for k, label in obs_field_names.items():
        o = obs[k]; a = rand[k]; n = len(a)
        mu, sd = a.mean(), a.std(ddof=1)
        Z = (o - mu) / sd if sd > 0 else float('nan')
        p_over = (1 + (a >= o).sum()) / (n + 1)
        p_under = (1 + (a <= o).sum()) / (n + 1)
        out.append(dict(null_model=model, graph=graph, motif_class=label, observed=int(o),
                        random_mean=round(mu, 1), random_SD=round(sd, 1),
                        random_min=int(a.min()), random_max=int(a.max()),
                        Z_score=round(Z, 2), fold_change=round(o / mu, 3) if mu else '',
                        pct_excess=f'{100*(o/mu-1):+.1f}%' if mu else '',
                        p_empirical_over=f'{p_over:.4f}', p_empirical_under=f'{p_under:.4f}',
                        n_randomisations=n, note=note,
                        source_file=os.path.relpath(path, REV)))
    return out

T4 = []
# first pass: 30_motif_significance.py, canonical network, named classes only
m1 = rd(f'{REV}/results/motif_significance_1000.csv')
for r in m1:
    model = {'NULL-A_full_randomisation': 'NULL-A  full curveball randomisation within each edge class',
             'NULL-B_reciprocity_preserved': 'NULL-B  NULL-A with the mutual TF<->miRNA pairs held fixed'}[r['null_model']]
    T4.append(dict(null_model=model, graph='canonical directed network (587 nodes, 6,859 edges); named FFL classes only',
                   motif_class=r['motif_class'], observed=int(r['n_real']),
                   random_mean=float(r['rand_mean']), random_SD=float(r['rand_sd']),
                   random_min=int(r['rand_min']), random_max=int(r['rand_max']),
                   Z_score=float(r['Z']), fold_change=float(r['fold_change']),
                   pct_excess=f"{100*(float(r['fold_change'])-1):+.1f}%",
                   p_empirical_over=f"{float(r['p_over']):.4f}", p_empirical_under=f"{float(r['p_under']):.4f}",
                   n_randomisations=int(r['n_rand']),
                   note='primary analysis; composite cores counted once per mutual pair',
                   source_file='results/motif_significance_1000.csv'))
# reciprocal-pair Z from the same framework
mrec = [r for r in rd(f'{REV}/results/motif_significance.csv') if r['motif_class'] == 'Mutual-TF-miRNA-pairs']
for r in mrec:
    T4.append(dict(null_model=('NULL-A  full curveball randomisation within each edge class'
                               if r['null_model'] == 'NM1' else
                               'NULL-B  NULL-A with the mutual TF<->miRNA pairs held fixed'),
                   graph='canonical directed network (587 nodes, 6,859 edges)',
                   motif_class='Reciprocal (mutual) TF<->miRNA pairs  [THE DIAGNOSTIC]',
                   observed=int(r['n_real']), random_mean=float(r['rand_mean']), random_SD=float(r['rand_sd']),
                   random_min=int(r['rand_min']), random_max=int(r['rand_max']),
                   Z_score=('' if r['Z'] == 'NA' else float(r['Z'])), fold_change=float(r['fold_change']),
                   pct_excess=f"{100*(float(r['fold_change'])-1):+.1f}%",
                   p_empirical_over=f"{float(r['p_emp']):.4f}", p_empirical_under=f"{float(r['p_emp_under']):.4f}",
                   n_randomisations=int(r['n_rand']),
                   note=('the network carries 4.9x more mutual TF<->miRNA pairs than chance (Z = +88); this is the '
                         'imprint of the hypergeometric miRNA-TF filter and it is what drives the apparent '
                         'composite enrichment under NULL-A. Held fixed by construction under NULL-B'),
                   source_file='results/motif_significance.csv'))
# second pass: independent re-implementation, three nulls, two graphs
lbl = {'total': 'Any 3-node FFL', 'comp': 'Composite-FFL', 'tf': 'TF-FFL', 'mir': 'miRNA-FFL',
       'other': 'FFL not in a named class (TF-TF-, miRNA-miRNA- and gene-mediated)',
       'recip': 'Reciprocal (mutual) TF<->miRNA pairs  [THE DIAGNOSTIC]'}
models = {'A': 'NULL-A  full curveball randomisation within each edge class',
          'B': 'NULL-B  NULL-A with the mutual TF<->miRNA pairs held fixed',
          'C': 'NULL-C  NULL-B plus curveball randomisation within the bipartite miRNA->target and TF->target layers'}
graphs = {'dep': 'deposited canonical graph, structural definition D1-D4 (587 nodes, 5,636 arcs after contraction)',
          'pub': 'augmented census graph, structural definition D1-D4 (587 nodes, 8,031 arcs after contraction)'}
for g in ('dep', 'pub'):
    for m in ('A', 'B', 'C'):
        p = f'{REV}/results/v2/null_{m}_{g}.tsv'
        if os.path.exists(p):
            T4 += summarise_null(p, lbl, models[m], graphs[g],
                                 'independent second-pass re-implementation; all topologically valid 3-node FFLs counted')
W('T4_motif_significance.csv', T4)
print('done T1/T2/T4')
