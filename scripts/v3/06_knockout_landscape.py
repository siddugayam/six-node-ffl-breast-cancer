#!/usr/bin/env python3
"""(D) Complete single- and double-knockout landscape of the canonical network.

All 587 single knockouts and all 587*586/2 = 171,791 unordered double knockouts are
evaluated exactly (no sampling).  For each pair the EXCESS OVER ADDITIVITY

    synergy(u,v) = damage(u,v) - [damage(u) + damage(v)]

is reported.  Positive synergy = the pair disrupts more than the sum of its parts, which is
the precise form of Iyengar's claim that multi-target perturbation can beat single-target
inhibition (Azeloglu & Iyengar 2015 CSH Perspect Biol 7:a005934; Berger & Iyengar 2009
Bioinformatics 25:2466).  The null distribution is the empirical distribution over all
random pairs, so no resampling assumption is needed.

Druggability: DGIdb interactions.tsv / drugs.tsv.
"""
import os, sys, csv, json, gzip, collections, itertools, time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, load_edges, RES, ffl_cores
from perturblib import Kernel

DGI = '/path/to/home/Desktop/DD/R_GPR/TCGA_PAN_CAN/CPTAC_PanCancer_Analysis/00_Raw_Data'
NPROC = int(os.environ.get('NPROC', '12'))

nt = load_nodes(); E = load_edges()
P = [(e['source'], e['target']) for e in E]
sg = {(e['source'], e['target']): e['sign_published'] for e in E}
NODES = sorted(nt)
cores = ffl_cores(set(P), nt)
K = Kernel(NODES, P, sg, nt, cores)
_, n_iter = K.influence(K.alive_all, return_iters=True)
print(f'WT: FFL={K.FFL_WT} reach_pairs={K.R_WT_all} influence_to_collagens={K.INF_WT_all:.4f} '
      f'(Neumann iterations to 1e-12: {n_iter})')

# ---------------------------------------------------------------- druggability
drugs = collections.defaultdict(set); appr = collections.defaultdict(set)
anti = collections.defaultdict(set)
with open(f'{DGI}/interactions.tsv') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        g = r['gene_name']
        if g in nt:
            drugs[g].add(r['drug_name'])
            if r['approved'] == 'TRUE':
                appr[g].add(r['drug_name'])
            if r['anti_neoplastic'] == 'TRUE':
                anti[g].add(r['drug_name'])
DRUGGABLE = sorted(drugs)
APPROVED = sorted(appr)
print(f'druggable network nodes: {len(DRUGGABLE)}  with an approved drug: {len(APPROVED)}  '
      f'with an antineoplastic: {len(anti)}')

part = collections.Counter()
for R_, M_, T_, c in cores:
    part[R_] += 1; part[M_] += 1; part[T_] += 1
deg = collections.Counter()
for u, v in P:
    deg[u] += 1; deg[v] += 1

# ---------------------------------------------------------------- singles
t0 = time.time()
single = {n: K.damage([n]) for n in NODES}
print(f'singles done {time.time()-t0:.1f}s')
srows = []
for n in NODES:
    d = single[n]
    srows.append(dict(node=n, type=nt[n], degree=deg[n], ffl_participation=part[n],
                      dFFL=round(d['dFFL'], 6), dReach=round(d['dReach'], 6),
                      dCOL_signed=('' if np.isnan(d['dCOL']) else round(d['dCOL'], 6)),
                      abs_dCOL=('' if np.isnan(d['dCOL']) else round(abs(d['dCOL']), 6)),
                      composite_damage=round(d['composite'], 6),
                      druggable=n in drugs, n_drugs=len(drugs.get(n, ())),
                      n_approved_drugs=len(appr.get(n, ())),
                      n_antineoplastic_drugs=len(anti.get(n, ())),
                      example_drugs=';'.join(sorted(appr.get(n, ()))[:5])))
srows.sort(key=lambda r: -r['composite_damage'])
with open(f'{RES}/systems_single_knockout.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(srows[0].keys())); w.writeheader(); w.writerows(srows)
print('\nTop 12 single knockouts (all nodes):')
for r in srows[:12]:
    print(f"  {r['node']:<14} {r['type']:<6} dFFL={r['dFFL']:.3f} dReach={r['dReach']:.4f} "
          f"dCOL={r['dCOL_signed']} comp={r['composite_damage']:.4f} drugs={r['n_drugs']}")
best_single = srows[0]
best_drug_single = next(r for r in srows if r['druggable'])
print(f"\nBEST SINGLE overall : {best_single['node']} comp={best_single['composite_damage']:.5f}")
print(f"BEST SINGLE druggable: {best_drug_single['node']} comp={best_drug_single['composite_damage']:.5f}")

# ---------------------------------------------------------------- all pairs
ALLP = list(itertools.combinations(range(len(NODES)), 2))
print(f'\nevaluating {len(ALLP):,} double knockouts on {NPROC} processes ...')

def work(chunk):
    out = []
    for i, j in chunk:
        d = K.damage([NODES[i], NODES[j]])
        out.append((i, j, d['dFFL'], d['dReach'],
                    -9.0 if np.isnan(d['dCOL']) else d['dCOL'], d['composite']))
    return out

CH = 4000
chunks = [ALLP[k:k + CH] for k in range(0, len(ALLP), CH)]
t0 = time.time(); res = []
with Pool(NPROC) as pool:
    for k, r in enumerate(pool.imap_unordered(work, chunks)):
        res.extend(r)
        if (k + 1) % 10 == 0:
            print(f'  {len(res):,}/{len(ALLP):,}  {time.time()-t0:.0f}s', flush=True)
print(f'doubles done {time.time()-t0:.1f}s')
arr = np.array(res, dtype=float)
order = np.lexsort((arr[:, 1], arr[:, 0]))
arr = arr[order]
np.save(f'{RES}/systems_double_knockout_matrix.npy', arr)

ii = arr[:, 0].astype(int); jj = arr[:, 1].astype(int)
comp = arr[:, 5]
s_comp = np.array([single[n]['composite'] for n in NODES])
add = s_comp[ii] + s_comp[jj]
syn = comp - add

# ---------------------------------------------------------------- summaries
def name(i):
    return NODES[i]

hdr = ['node_A', 'node_B', 'type_A', 'type_B', 'dFFL', 'dReach', 'dCOL_signed',
       'composite_damage', 'sum_of_singles', 'synergy_excess_over_additive',
       'druggable_A', 'druggable_B', 'both_druggable', 'drugs_A', 'drugs_B']
def row_for(k):
    a, b = name(ii[k]), name(jj[k])
    return [a, b, nt[a], nt[b], round(arr[k, 2], 6), round(arr[k, 3], 6),
            ('' if arr[k, 4] == -9.0 else round(arr[k, 4], 6)),
            round(comp[k], 6), round(add[k], 6), round(syn[k], 6),
            a in drugs, b in drugs, (a in drugs) and (b in drugs),
            ';'.join(sorted(appr.get(a, ()))[:3]), ';'.join(sorted(appr.get(b, ()))[:3])]

top_damage = np.argsort(-comp)[:2000]
top_syn = np.argsort(-syn)[:2000]
with open(f'{RES}/systems_double_knockout_top.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['rank_by'] + hdr)
    for r, k in enumerate(top_damage):
        w.writerow(['damage'] + row_for(k))
    for r, k in enumerate(top_syn):
        w.writerow(['synergy'] + row_for(k))

both_drug = np.array([(NODES[a] in drugs) and (NODES[b] in drugs) for a, b in zip(ii, jj)])
dk = np.where(both_drug)[0]
dorder = dk[np.argsort(-comp[dk])]
with open(f'{RES}/systems_double_knockout_druggable.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(hdr)
    for k in dorder[:5000]:
        w.writerow(row_for(k))
print(f'\ndruggable pairs: {both_drug.sum():,}')

with gzip.open(f'{RES}/systems_double_knockout_all.csv.gz', 'wt', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['node_A', 'node_B', 'dFFL', 'dReach', 'dCOL_signed', 'composite_damage',
                'sum_of_singles', 'synergy'])
    for k in range(len(comp)):
        w.writerow([name(ii[k]), name(jj[k]), round(arr[k, 2], 6), round(arr[k, 3], 6),
                    ('' if arr[k, 4] == -9.0 else round(arr[k, 4], 6)),
                    round(comp[k], 6), round(add[k], 6), round(syn[k], 6)])

# ---------------------------------------------------------------- the key comparison
best_pair_k = int(np.argmax(comp))
best_drug_pair_k = int(dk[np.argmax(comp[dk])])
bs = best_single['composite_damage']
bds = best_drug_single['composite_damage']
q = np.quantile(comp, [0.5, 0.9, 0.99, 0.999])
n_beat = int((comp > bs).sum())
n_beat_drug = int((comp[dk] > bds).sum())

summary = dict(
    n_nodes=len(NODES), n_pairs=int(len(comp)),
    best_single=dict(node=best_single['node'], composite=bs, dFFL=best_single['dFFL'],
                     dReach=best_single['dReach'], dCOL=best_single['dCOL_signed']),
    best_single_druggable=dict(node=best_drug_single['node'], composite=bds),
    best_pair=dict(A=name(ii[best_pair_k]), B=name(jj[best_pair_k]),
                   composite=float(comp[best_pair_k]),
                   sum_of_singles=float(add[best_pair_k]),
                   synergy=float(syn[best_pair_k]),
                   ratio_to_best_single=float(comp[best_pair_k] / bs)),
    best_druggable_pair=dict(A=name(ii[best_drug_pair_k]), B=name(jj[best_drug_pair_k]),
                             composite=float(comp[best_drug_pair_k]),
                             sum_of_singles=float(add[best_drug_pair_k]),
                             synergy=float(syn[best_drug_pair_k]),
                             ratio_to_best_druggable_single=float(comp[best_drug_pair_k] / bds)),
    random_pair_null=dict(mean=float(comp.mean()), sd=float(comp.std(ddof=1)),
                          median=float(q[0]), q90=float(q[1]), q99=float(q[2]),
                          q999=float(q[3]), max=float(comp.max())),
    n_pairs_beating_best_single=n_beat,
    frac_pairs_beating_best_single=n_beat / len(comp),
    n_druggable_pairs=int(both_drug.sum()),
    n_druggable_pairs_beating_best_druggable_single=n_beat_drug,
    synergy=dict(mean=float(syn.mean()), sd=float(syn.std(ddof=1)),
                 n_positive=int((syn > 0).sum()), frac_positive=float((syn > 0).mean()),
                 max=float(syn.max()), min=float(syn.min()),
                 q99=float(np.quantile(syn, 0.99))),
    z_of_best_pair_vs_random_pairs=float((comp[best_pair_k] - comp.mean()) / comp.std(ddof=1)),
    empirical_p_best_pair=float((np.sum(comp >= comp[best_pair_k])) / len(comp)),
)
json.dump(summary, open(f'{RES}/systems_combination_summary.json', 'w'), indent=2)
print(json.dumps(summary, indent=2))

print('\nTop 15 pairs by total damage:')
for k in top_damage[:15]:
    print(f'  {name(ii[k]):<14}+{name(jj[k]):<14} comp={comp[k]:.5f} '
          f'sum_singles={add[k]:.5f} synergy={syn[k]:+.5f}')
print('\nTop 15 pairs by SYNERGY (excess over additivity):')
for k in top_syn[:15]:
    print(f'  {name(ii[k]):<14}+{name(jj[k]):<14} comp={comp[k]:.5f} '
          f'sum_singles={add[k]:.5f} synergy={syn[k]:+.5f}')
print('\nTop 15 DRUGGABLE pairs by total damage:')
for k in dorder[:15]:
    print(f'  {name(ii[k]):<14}+{name(jj[k]):<14} comp={comp[k]:.5f} synergy={syn[k]:+.5f}')
