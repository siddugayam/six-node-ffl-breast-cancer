#!/usr/bin/env python3
"""(D cont.) Does multi-target perturbation actually beat single-target inhibition here?

Uses the complete 171,991-pair double-knockout landscape from 06_knockout_landscape.py.

  * how far the best pair exceeds the best single node
  * whether the top-damage pairs are super- or sub-additive
  * the empirical null over all random pairs, and the z of the most synergistic pair
  * targeted multi-node knockouts of biologically defined sets (the miR-29 family, the
    NF-kB module, the top FFL hubs, the miR-29/collagen module)
"""
import os, sys, csv, json, collections, itertools
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, ffl_cores
from perturblib import Kernel

DGI = '/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data'
nt = load_nodes(); E = load_edges()
P = [(e['source'], e['target']) for e in E]
sg = {(e['source'], e['target']): e['sign_published'] for e in E}
NODES = sorted(nt)
cores = ffl_cores(set(P), nt)
K = Kernel(NODES, P, sg, nt, cores)

arr = np.load(f'{RES}/systems_double_knockout_matrix.npy')
ii = arr[:, 0].astype(int); jj = arr[:, 1].astype(int); comp = arr[:, 5]
single = {n: K.damage([n]) for n in NODES}
s_comp = np.array([single[n]['composite'] for n in NODES])
syn = comp - (s_comp[ii] + s_comp[jj])

drugs = collections.defaultdict(set); appr = collections.defaultdict(set)
with open(f'{DGI}/interactions.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['gene_name'] in nt:
            drugs[r['gene_name']].add(r['drug_name'])
            if r['approved'] == 'TRUE':
                appr[r['gene_name']].add(r['drug_name'])

out = {}
best_single_i = int(np.argmax(s_comp))
out['best_single'] = dict(node=NODES[best_single_i], damage=float(s_comp[best_single_i]))
top100 = np.argsort(-comp)[:100]
out['top100_pairs'] = dict(
    n_super_additive=int((syn[top100] > 0).sum()),
    n_sub_additive=int((syn[top100] < 0).sum()),
    mean_synergy=float(syn[top100].mean()),
    mean_ratio_to_best_single=float((comp[top100] / s_comp[best_single_i]).mean()),
    max_ratio_to_best_single=float((comp[top100] / s_comp[best_single_i]).max()))
k_syn = int(np.argmax(syn))
out['most_synergistic_pair'] = dict(
    A=NODES[ii[k_syn]], B=NODES[jj[k_syn]], damage=float(comp[k_syn]),
    sum_of_singles=float(s_comp[ii[k_syn]] + s_comp[jj[k_syn]]),
    synergy=float(syn[k_syn]),
    z_vs_all_pairs=float((syn[k_syn] - syn.mean()) / syn.std(ddof=1)),
    percent_of_best_single_damage=float(100 * syn[k_syn] / s_comp[best_single_i]))
out['synergy_distribution'] = dict(
    mean=float(syn.mean()), sd=float(syn.std(ddof=1)),
    frac_super_additive=float((syn > 0).mean()),
    q999=float(np.quantile(syn, 0.999)), max=float(syn.max()))

# how much does a second target buy you on average, conditional on the first being the best?
mask_b = (ii == best_single_i) | (jj == best_single_i)
part = comp[mask_b]
bi = ii[mask_b]; bj = jj[mask_b]
kbest = int(np.argmax(part))
partner = NODES[bj[kbest] if bi[kbest] == best_single_i else bi[kbest]]
out['adding_a_second_target_to_the_best_single'] = dict(
    first_target=NODES[best_single_i], n_partners=int(mask_b.sum()),
    mean_damage=float(part.mean()), max_damage=float(part.max()),
    best_partner=partner,
    mean_gain_over_best_single=float(part.mean() - s_comp[best_single_i]),
    max_gain_over_best_single=float(part.max() - s_comp[best_single_i]))

# ---------------------------------------------------------------- targeted sets
SETS = {
    'miR-29 family': ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c'],
    'miR-29/collagen module': ['hsa-miR-29a', 'hsa-miR-29b', 'hsa-miR-29c', 'COL1A1', 'COL3A1'],
    'NF-kB module (RELA+NFKB1)': ['RELA', 'NFKB1'],
    'top3 FFL hubs (SP1,RELA,NFKB1)': ['SP1', 'RELA', 'NFKB1'],
    'top5 FFL hubs': ['SP1', 'RELA', 'NFKB1', 'VEGFA', 'CCND2'],
    'best druggable pair (RELA+SP1)': ['RELA', 'SP1'],
    'MKL1 + miR-143 (most synergistic)': ['MKL1', 'hsa-miR-143'],
    'miR-17-92 seed cluster': ['hsa-miR-17', 'hsa-miR-18a', 'hsa-miR-19a', 'hsa-miR-19b',
                               'hsa-miR-20a'],
}
rows = []
for lab, s in SETS.items():
    s = [x for x in s if x in nt]
    d = K.damage(s)
    add = sum(single[x]['composite'] for x in s)
    rows.append(dict(set_name=lab, members=';'.join(s), n=len(s),
                     dFFL=round(d['dFFL'], 6), dReach=round(d['dReach'], 6),
                     dCOL_signed=('' if np.isnan(d['dCOL']) else round(d['dCOL'], 6)),
                     composite_damage=round(d['composite'], 6),
                     sum_of_single_damages=round(add, 6),
                     synergy=round(d['composite'] - add, 6),
                     ratio_to_best_single=round(d['composite'] / s_comp[best_single_i], 4),
                     all_druggable=all(x in drugs for x in s),
                     n_druggable=sum(1 for x in s if x in drugs)))
    print(f"{lab:<36} n={len(s)} comp={d['composite']:.5f} sum_singles={add:.5f} "
          f"syn={d['composite']-add:+.5f} ratio={d['composite']/s_comp[best_single_i]:.2f}x")
with open(f'{RES}/systems_targeted_set_knockouts.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---------------------------------------------------------------- best druggable combos with
# both partners carrying an APPROVED drug
both_appr = np.array([(NODES[a] in appr) and (NODES[b] in appr) for a, b in zip(ii, jj)])
ka = np.where(both_appr)[0]
ka = ka[np.argsort(-comp[ka])]
arows = []
for k in ka[:200]:
    a, b = NODES[ii[k]], NODES[jj[k]]
    arows.append(dict(node_A=a, node_B=b, type_A=nt[a], type_B=nt[b],
                      composite_damage=round(float(comp[k]), 6),
                      synergy=round(float(syn[k]), 6),
                      ratio_to_best_single=round(float(comp[k] / s_comp[best_single_i]), 4),
                      approved_drugs_A=';'.join(sorted(appr[a])[:4]),
                      approved_drugs_B=';'.join(sorted(appr[b])[:4])))
with open(f'{RES}/systems_druggable_combinations_approved.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(arows[0].keys())); w.writeheader(); w.writerows(arows)
out['n_pairs_both_with_approved_drug'] = int(both_appr.sum())
out['best_approved_drug_pair'] = arows[0]
print(f"\npairs where BOTH nodes have an approved drug: {both_appr.sum():,}")
print(f"  best: {arows[0]['node_A']} + {arows[0]['node_B']} "
      f"comp={arows[0]['composite_damage']} ({arows[0]['ratio_to_best_single']}x best single)")
for a in arows[:8]:
    print(f"   {a['node_A']:<10}+{a['node_B']:<10} {a['composite_damage']:.5f} "
          f"syn={a['synergy']:+.5f}  [{a['approved_drugs_A'][:28]} | {a['approved_drugs_B'][:28]}]")

json.dump(out, open(f'{RES}/systems_combination_stats.json', 'w'), indent=2)
print('\n' + json.dumps(out, indent=2))
