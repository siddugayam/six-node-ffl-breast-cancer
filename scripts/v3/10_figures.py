#!/usr/bin/env python3
"""Figures for the systems-pharmacology / network-control analysis."""
import os, sys, csv, json, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlib import load_nodes, RES, FIG

plt.rcParams.update({'font.size': 8, 'axes.linewidth': .8, 'figure.dpi': 160,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'font.family': 'DejaVu Sans'})
C = dict(tf='#2166ac', mir='#b2182b', gene='#4d9221', grey='#9e9e9e', hi='#e08214')
tcol = {'TF': C['tf'], 'miRNA': C['mir'], 'Gene': C['gene']}
nt = load_nodes()


def save(fig, name):
    fig.savefig(f'{FIG}/{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


# ============================================================ Fig 1: architecture
arch = json.load(open(f'{RES}/systems_architecture_summary.json'))
nulls = json.load(open(f'{RES}/systems_architecture_nulls.json'))
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

# 1A bow-tie
bt = arch['bowtie']; bty = arch['bowtie_types']
a = ax[0]; a.set_xlim(0, 10); a.set_ylim(0, 6); a.axis('off')
a.add_patch(Polygon([[0.3, 3], [2.6, 5.0], [2.6, 1.0]], closed=True,
                    fc='#cfe3f5', ec='#4a708b'))
a.add_patch(Circle((4.3, 3), 1.35, fc='#f6c6c6', ec='#a33'))
a.add_patch(Polygon([[6.0, 5.0], [6.0, 1.0], [8.3, 3]], closed=True,
                    fc='#cfe8cf', ec='#4a8b4a'))
for x0, x1 in [(2.7, 2.9), (5.7, 5.9)]:
    a.add_patch(FancyArrowPatch((x0, 3), (x1 + .25, 3), arrowstyle='-|>',
                                mutation_scale=9, lw=1.1, color='k'))
a.text(1.55, 3.0, f"IN\n{bt['IN']}", ha='center', va='center', fontsize=8.5, fontweight='bold')
a.text(4.3, 3.25, f"CORE\n{bt['CORE']}", ha='center', va='center', fontsize=9, fontweight='bold')
a.text(4.3, 2.35, f"{bty['CORE']['TF']} TF\n{bty['CORE']['miRNA']} miRNA", ha='center',
       va='center', fontsize=6.4, color='#7a1f1f')
a.text(6.95, 3.25, f"OUT\n{bt['OUT']}", ha='center', va='center', fontsize=8.5,
       fontweight='bold')
a.text(6.95, 2.4, f"{bty['OUT']['Gene']} genes", ha='center', va='center', fontsize=6.4,
       color='#2d5a2d')
a.text(5, 0.35, f"tubes {bt['TUBES']} · tendrils "
                f"{bt['TENDRILS_from_IN']+bt['TENDRILS_to_OUT']} · "
                f"{arch['edges_within_core_frac']*100:.0f}% of edges inside the core",
       ha='center', fontsize=6.3, color='#555')
a.set_title('A  Bow-tie decomposition (587 nodes)', loc='left', fontsize=9)

# 1B observed vs null for the architecture statistics
a = ax[1]
mets = [('LSCC', 'core size'), ('flow_hierarchy', 'flow hierarchy'),
        ('GRC', 'global reaching\ncentrality'), ('F0', 'trophic\nincoherence'),
        ('transitivity', 'transitivity'), ('ffl_unique', '3-node FFLs')]
zs_a = [nulls['NULL_A'][m]['z'] for m, _ in mets]
zs_b = [nulls['NULL_B'][m]['z'] for m, _ in mets]
xp = np.arange(len(mets))
a.bar(xp - .2, zs_a, .38, label='NULL-A (degree preserved)', color='#bbbbbb', ec='k', lw=.4)
a.bar(xp + .2, zs_b, .38, label='NULL-B (+ reciprocity preserved)', color=C['tf'], ec='k', lw=.4)
a.axhline(0, color='k', lw=.7)
for yv, ls in [(1.96, ':'), (-1.96, ':')]:
    a.axhline(yv, color='#c33', lw=.7, ls=ls)
a.set_xticks(xp); a.set_xticklabels([l for _, l in mets], rotation=35, ha='right', fontsize=6.6)
a.set_ylabel('z vs randomised networks')
a.legend(fontsize=5.9, frameon=False, loc='lower left')
a.set_title('B  Architecture vs randomised networks', loc='left', fontsize=9)

# 1C trophic levels by node class
a = ax[2]
tl = collections.defaultdict(list)
for r in csv.DictReader(open(f'{RES}/systems_trophic_levels.csv')):
    tl[r['type']].append(float(r['trophic_level']))
for i, k in enumerate(['miRNA', 'TF', 'Gene']):
    v = np.array(tl[k])
    a.scatter(np.full(len(v), i) + np.random.uniform(-.22, .22, len(v)), v, s=3.5,
              color=tcol[k], alpha=.45, lw=0)
    a.plot([i - .32, i + .32], [v.mean()] * 2, color='k', lw=1.6)
    a.text(i, v.mean() + .12, f'{v.mean():.2f}', ha='center', fontsize=6.5)
a.set_xticks([0, 1, 2]); a.set_xticklabels(['miRNA', 'TF', 'Gene'])
a.set_ylabel('trophic level')
a.set_title(f"C  Trophic levels (incoherence F0={arch['trophic_incoherence_F0']:.2f})",
            loc='left', fontsize=9)
fig.tight_layout()
save(fig, 'Fig_v3_A_architecture')

# ============================================================ Fig 2: controllability
ctrl = list(csv.DictReader(open(f'{RES}/systems_controllability_nodes.csv')))
summ = json.load(open(f'{RES}/systems_controllability_summary.json'))
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

a = ax[0]
ffl = np.array([float(r['ffl_participation']) for r in ctrl])
dfq = np.array([float(r['driver_frequency']) for r in ctrl])
ish = np.array([r['is_FFL_hub'] == 'True' for r in ctrl])
a.scatter(ffl[~ish] + .8, dfq[~ish], s=7, color=C['grey'], alpha=.55, lw=0, label='other nodes')
a.scatter(ffl[ish] + .8, dfq[ish], s=13, color=C['hi'], ec='k', lw=.25, label='FFL hub')
a.set_xscale('log'); a.set_xlabel('FFL participation (+1)')
a.set_ylabel(f"driver frequency\n({summ['n_matchings_sampled']} maximum matchings)")
a.legend(fontsize=6, frameon=False)
r, p = summ['spearman_driverfreq_vs_fflparticipation']
a.set_title(f'A  Driver frequency falls with FFL participation\n'
            f'(Spearman {r:+.2f}, p={p:.1e}; degree-confounded)',
            loc='left', fontsize=8.5)

a = ax[1]
enr = list(csv.DictReader(open(f'{RES}/systems_controllability_enrichment.csv')))
sel = [(g, 'driver_node') for g in ('FFL_hub (hub_ffl==TRUE)', 'top20_FFL_participation',
                                    'degree_hub', 'miR-29/collagen module')]
labs = ['FFL hubs\n(n=60)', 'top-20 FFL\nparticipation', 'degree hubs', 'miR-29/collagen\nmodule (n=5)']
fg = []; fr = []
for (g, pr) in sel:
    e = next(x for x in enr if x['group'] == g and x['property'] == pr)
    fg.append(float(e['frac_group'])); fr.append(float(e['frac_rest']))
xp = np.arange(len(sel))
a.bar(xp - .19, fg, .36, color=C['hi'], ec='k', lw=.4, label='group')
a.bar(xp + .19, fr, .36, color=C['grey'], ec='k', lw=.4, label='rest of network')
for i, (g, pr) in enumerate(sel):
    e = next(x for x in enr if x['group'] == g and x['property'] == pr)
    pv = float(e['p'])
    a.text(i, max(fg[i], fr[i]) + .015, ('p=%.3g' % pv) if pv < .05 else 'ns',
           ha='center', fontsize=6)
a.set_xticks(xp); a.set_xticklabels(labs, fontsize=6.2)
a.set_ylabel('fraction that are driver nodes')
a.legend(fontsize=6, frameon=False)
dgc = json.load(open(f'{RES}/systems_controllability_degree_control.json'))
_ph = dgc['logistic_driver']['is_FFL_hub']['p']
_pm = dgc['degree_matched']['p_one_sided']
a.set_title(f'B  Raw depletion is a DEGREE effect\nadjusted for degree: p={_ph:.2f}; '
            f'degree-matched: p={_pm:.2f}', loc='left', fontsize=8.5)

a = ax[2]
dc = collections.Counter(r['deletion_class'] for r in ctrl)
cr = collections.Counter(r['control_role'] for r in ctrl)
ks = ['indispensable', 'neutral', 'dispensable']
a.barh([2, 1, 0], [dc[k] for k in ks], color=['#b2182b', '#cccccc', '#4d9221'], ec='k', lw=.4)
for i, k in enumerate(ks[::-1]):
    a.text(dc[k] + 6, i, f'{dc[k]}', va='center', fontsize=7)
a.set_yticks([2, 1, 0]); a.set_yticklabels(ks, fontsize=7)
a.set_xlabel('nodes')
_indisp = [r for r in ctrl if r['deletion_class'] == 'indispensable']
_ind_hub = sum(1 for r in _indisp if r['is_FFL_hub'] == 'True')
_ind_mod = sum(1 for r in _indisp if r['is_mir29_module'] == 'True')
a.set_title(f"C  N_D = {summ['N_D']} drivers ({summ['n_D_over_N']*100:.1f}% of nodes)\n"
            f"{_ind_hub}/{len(_indisp)} indispensable nodes are FFL hubs, "
            f"{_ind_mod}/{len(_indisp)} are miR-29 module members",
            loc='left', fontsize=8.5)
fig.tight_layout()
save(fig, 'Fig_v3_B_controllability')

# ============================================================ Fig 3: information flow
inf = list(csv.DictReader(open(f'{RES}/systems_information_flow.csv')))
deg = np.array([float(r['degree_total']) for r in inf])
cfb = np.array([float(r['current_flow_betweenness']) for r in inf])
rwt = np.array([float(r['rw_throughput_betweenness']) for r in inf])
nm = [r['node'] for r in inf]
ty = [r['type'] for r in inf]
cmpd = list(csv.DictReader(open(f'{RES}/systems_flow_vs_degree_comparison.csv')))
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

for pfx, a, v, lab, key in [('A  Undirected current flow', ax[0], cfb,
                             'current-flow betweenness\n(undirected)',
                             'current_flow_betweenness'),
                            ('B  Directed signal throughput', ax[1], rwt,
                             'directed random-walk\nthroughput',
                             'rw_throughput_betweenness')]:
    for t in ('Gene', 'TF', 'miRNA'):
        m = np.array([x == t for x in ty])
        a.scatter(deg[m], v[m], s=8, color=tcol[t], alpha=.6, lw=0, label=t)
    rho = float(next(c for c in cmpd if c['measure'] == key)['spearman_vs_degree'])
    ov = next(c for c in cmpd if c['measure'] == key)['top20_overlap']
    a.set_xlabel('degree'); a.set_ylabel(lab)
    a.set_title(f'{pfx}\nvs degree: Spearman rho = {rho:+.3f}; top-20 overlap {ov}/20',
                loc='left', fontsize=8)
    a.legend(fontsize=6, frameon=False)
    for n in ('VEGFA', 'CCND2', 'SP1', 'TP53', 'hsa-miR-21'):
        i = nm.index(n)
        a.annotate(n, (deg[i], v[i]), fontsize=5.6, xytext=(3, 2),
                   textcoords='offset points')

a = ax[2]
t20 = list(csv.DictReader(open(f'{RES}/systems_top20_by_measure.csv')))
meas = ['degree_total', 'current_flow_betweenness', 'information_centrality',
        'rw_throughput_betweenness', 'ffl_participation']
short = ['degree', 'current-flow\nbetweenness', 'information\ncentrality',
         'directed RW\nthroughput', 'FFL\nparticipation']
sets = {m: {r[m] for r in t20} for m in meas}
J = np.zeros((len(meas), len(meas)))
for i, m1 in enumerate(meas):
    for j, m2 in enumerate(meas):
        J[i, j] = len(sets[m1] & sets[m2]) / len(sets[m1] | sets[m2])
im = a.imshow(J, cmap='Blues', vmin=0, vmax=1)
a.set_xticks(range(len(meas))); a.set_xticklabels(short, rotation=45, ha='right', fontsize=5.8)
a.set_yticks(range(len(meas))); a.set_yticklabels(short, fontsize=5.8)
for i in range(len(meas)):
    for j in range(len(meas)):
        a.text(j, i, f'{J[i,j]:.2f}', ha='center', va='center', fontsize=5.6,
               color='white' if J[i, j] > .5 else 'k')
a.set_title('C  Top-20 Jaccard overlap', loc='left', fontsize=8.5)
a.spines[:].set_visible(True)
fig.tight_layout()
save(fig, 'Fig_v3_C_information_flow')

# ============================================================ Fig 4: perturbation
sk = list(csv.DictReader(open(f'{RES}/systems_single_knockout.csv')))
cs = json.load(open(f'{RES}/systems_combination_stats.json'))
arr = np.load(f'{RES}/systems_double_knockout_matrix.npy')
comp = arr[:, 5]
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))

a = ax[0]
d = np.array([float(r['composite_damage']) for r in sk])
col = [tcol[r['type']] for r in sk]
a.bar(range(30), d[:30], color=col[:30], ec='k', lw=.3)
a.set_xticks(range(30))
a.set_xticklabels([r['node'] for r in sk[:30]], rotation=90, fontsize=5.3)
a.set_ylabel('composite damage')
for i, r in enumerate(sk[:30]):
    if r['druggable'] == 'True':
        a.text(i, d[i] + .0012, '*', ha='center', fontsize=7, color='#333')
a.set_title('A  Top single knockouts (* = druggable, DGIdb)', loc='left', fontsize=8.5)

a = ax[1]
a.hist(comp, bins=np.linspace(0, comp.max(), 120), color=C['grey'], log=True)
bs = cs['best_single']['damage']
a.axvline(bs, color=C['tf'], lw=1.2)
a.text(bs, 3e4, f"  best single\n  ({cs['best_single']['node']})", fontsize=6, color=C['tf'])
cmb = json.load(open(f'{RES}/systems_combination_summary.json'))
bp = comp.max()
a.axvline(bp, color='#b2182b', lw=1.2)
a.text(bp, 3e3, f"  best pair\n ({cmb['best_pair']['A']}+{cmb['best_pair']['B']})",
       fontsize=6, color='#b2182b', ha='right')
a.set_xlabel('composite damage'); a.set_ylabel('number of node pairs')
a.set_title(f"B  All {len(comp):,} double knockouts", loc='left', fontsize=8.5)

a = ax[2]
syn = comp - (arr[:, 5] * 0)  # placeholder replaced below
NODES = [r['node'] for r in sorted(sk, key=lambda x: x['node'])]
sd = {r['node']: float(r['composite_damage']) for r in sk}
nodes_sorted = sorted(sd)
sc = np.array([sd[n] for n in nodes_sorted])
ii = arr[:, 0].astype(int); jj = arr[:, 1].astype(int)
syn = comp - (sc[ii] + sc[jj])
a.hist(syn, bins=200, color=C['grey'], log=True)
a.axvline(0, color='k', lw=.8)
mx = cs['most_synergistic_pair']
a.axvline(mx['synergy'], color='#b2182b', lw=1.2)
a.text(mx['synergy'], 2e3, f"  {mx['A']}+{mx['B']}\n  z={mx['z_vs_all_pairs']:.1f}",
       fontsize=6, color='#b2182b', ha='right')
a.set_xlabel('synergy = damage(A,B) - [damage(A)+damage(B)]')
a.set_ylabel('number of node pairs')
a.set_title(f"C  {cs['synergy_distribution']['frac_super_additive']*100:.0f}% of pairs "
            f"super-additive;\ntop-100 damage pairs: "
            f"{cs['top100_pairs']['n_sub_additive']}/100 SUB-additive",
            loc='left', fontsize=8.5)
fig.tight_layout()
save(fig, 'Fig_v3_D_perturbation')

# ============================================================ Fig 5: dynamics
dyn = json.load(open(f'{RES}/systems_dynamics_summary.json'))
dn = json.load(open(f'{RES}/systems_dynamics_nulls.json'))
R = np.load(f'{RES}/systems_response_matrix_published.npy')
de = {}
for p in ('BRCA_DEX_genes.csv', 'BRCA_DEX_mirnas.csv'):
    for r in csv.DictReader(open(f'/path/to/revision/results/{p}')):
        de[r['feature']] = float(r['logFC'])
NN = sorted(nt); ix = {n: i for i, n in enumerate(NN)}
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

best = dyn['published']['best_single']
a = ax[0]
s = ix[best['node']]; sgn_ = int(best['direction'])
xs = []; ys = []; cs2 = []
for n in NN:
    if n in de and n != best['node']:
        xs.append(sgn_ * R[ix[n], s]); ys.append(de[n]); cs2.append(tcol[nt[n]])
a.scatter(xs, ys, s=7, c=cs2, alpha=.6, lw=0)
a.axhline(0, color='k', lw=.5); a.axvline(0, color='k', lw=.5)
a.set_xscale('symlog', linthresh=1e-4)
a.set_xlabel(f"model steady-state response to {best['node']} ({sgn_:+d})")
a.set_ylabel('observed tumour vs normal logFC')
a.set_title(f"A  Best single perturbation\nSpearman {best['spearman']:+.3f}, "
            f"FWER p={dyn['published']['fwer_p_best_single']:.3f}", loc='left', fontsize=8.5)

a = ax[1]
for tag, c, lab in [('published', C['tf'], 'published signs (all TF->target = +1)'),
                    ('curated', '#b2182b', 'TRRUST/TransmiR curated signs')]:
    lp = list(csv.DictReader(open(f'{RES}/systems_dynamics_lasso_path_{tag}.csv')))
    a.plot([int(r['order']) for r in lp], [float(r['cumulative_R2']) for r in lp],
           'o-', ms=2.6, lw=1.1, color=c, label=lab)
a.set_xlabel('number of simultaneous upstream perturbations')
a.set_ylabel('cumulative $R^2$ of the observed phenotype')
a.legend(fontsize=5.9, frameon=False)
a.set_title('B  How many drivers to reproduce the phenotype', loc='left', fontsize=8.5)

a = ax[2]
b = dn['boolean']
a.bar([0], [b['observed_agreement']], .55, color=C['hi'], ec='k', lw=.4, label='observed network')
a.bar([1], [b['null_permuted_seed_mean']], .55, yerr=[b['null_permuted_seed_sd']],
      color=C['grey'], ec='k', lw=.4, capsize=3, label='permuted seed')
a.bar([2], [b['null_rewired_mean']], .55, yerr=[b['null_rewired_sd']],
      color='#cccccc', ec='k', lw=.4, capsize=3, label='rewired network')
a.axhline(.5, color='#c33', lw=.7, ls=':')
a.set_xticks([0, 1, 2]); a.set_xticklabels(['observed', 'permuted\nseed', 'rewired\nnetwork'],
                                           fontsize=6.6)
a.set_ylabel('attractor / phenotype sign agreement')
a.set_ylim(0, 1)
a.text(0, b['observed_agreement'] + .03, f"{b['observed_agreement']:.3f}", ha='center', fontsize=6.5)
a.text(1, .70, f"p={b['p_permuted_seed']:.1e}", ha='center', fontsize=6)
a.text(2, .78, f"p={b['p_rewired']:.1e}", ha='center', fontsize=6)
a.set_title(f"C  Boolean threshold model\n(limit cycle, period {b['period']}, "
            f"reached in {b['steps']} steps)", loc='left', fontsize=8.5)
fig.tight_layout()
save(fig, 'Fig_v3_E_dynamics')
# ============================================================ Fig 6: hierarchy test
reg = list(csv.DictReader(open(f'{RES}/systems_regulator_subnetwork.csv')))
regn = list(csv.DictReader(open(f'{RES}/systems_regulator_subnetwork_nulls.csv')))
full = next(r for r in reg if r['network'].startswith('full'))
ronly = next(r for r in reg if r['network'].startswith('regulators'))
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

a = ax[0]
labs = ['full network\n(587 nodes)', 'regulators only\n(TF + miRNA)']
core = [float(r['frac_nodes_in_largest_SCC']) * 100 for r in (full, ronly)]
a.bar([0, 1], core, .55, color=[C['grey'], C['tf']], ec='k', lw=.4)
for i, v in enumerate(core):
    a.text(i, v + 1.5, f'{v:.1f}%', ha='center', fontsize=7)
a.set_xticks([0, 1]); a.set_xticklabels(labs, fontsize=6.6)
a.set_ylabel('% of nodes in the largest SCC')
a.set_ylim(0, 105)
a.set_title('A  Where the recurrence is', loc='left', fontsize=8.5)

a = ax[1]
ks = ['flow_hierarchy', 'global_reaching_centrality', 'trophic_incoherence_F0', 'transitivity']
sl = ['flow\nhierarchy', 'global reaching\ncentrality', 'trophic\nincoherence', 'transitivity']
zs = []
for k in ks:
    zs.append(float(next(r for r in regn if r['metric'] == k)['z']))
a.bar(range(len(ks)), zs, .55, color=C['tf'], ec='k', lw=.4)
a.axhline(0, color='k', lw=.7)
for yv in (1.96, -1.96):
    a.axhline(yv, color='#c33', lw=.7, ls=':')
a.set_xticks(range(len(ks))); a.set_xticklabels(sl, rotation=35, ha='right', fontsize=6.4)
a.set_ylabel('z vs degree-preserving null')
a.set_title('B  Regulator subnetwork vs its own null', loc='left', fontsize=8.5)

a = ax[2]
dep = [int(full['condensation_depth_layers']), int(ronly['condensation_depth_layers'])]
a.bar([0, 1], dep, .55, color=[C['grey'], C['tf']], ec='k', lw=.4)
for i, v in enumerate(dep):
    a.text(i, v + .05, str(v), ha='center', fontsize=7)
a.set_xticks([0, 1]); a.set_xticklabels(labs, fontsize=6.6)
a.set_ylabel('layers in the condensation DAG')
a.set_title('C  How many genuine hierarchical layers', loc='left', fontsize=8.5)
fig.tight_layout()
save(fig, 'Fig_v3_F_hierarchy_test')

print('\nall figures written to', FIG)
