#!/usr/bin/env python3
"""Fig. S5 (Fig_v3_C_information_flow) only: the Fig 3 block of scripts/05_network_architecture/10_figures.py, copied verbatim,
with RES/FIG pointed at out_<mode>_hash0/ (the other v3 figures need analyses not re-run here).
usage: python3 fig_s5_information_flow.py <full|nolegacy>"""
import os, sys, csv, json, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch

sys.dont_write_bytecode = True
sys.path.insert(0, '/path/to/revision/scripts/05_network_architecture')
from netlib import load_nodes
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, f'out_{sys.argv[1]}_hash0'); FIG = os.path.join(RES, 'fig'); os.makedirs(FIG, exist_ok=True)

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
