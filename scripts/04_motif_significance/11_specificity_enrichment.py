#!/usr/bin/env python3
"""
11 -- The scale-free half of the specificity control.

Raw FFL density (script 10) cannot by itself separate "breast-cancer-specific architecture"
from "these nodes simply have more known interactions", because FFL counts grow super-linearly
with edge count.  The scale-free question is:

    within its OWN degree- and type-preserving randomised ensemble, is the breast-cancer
    network more FFL-enriched than a disease-unrelated network of the same size is within ITS
    own ensemble?

For every network built in script 10 this script therefore computes
    fold = N_FFL(real) / mean(N_FFL(randomised))       and     Z = (N-mu)/sigma
using the SAME null model NM2 as script 08 (Maslov-Sneppen switching confined to each
edge_type x node-type layer, with mutual TF<->miRNA pairs swapped as units so reciprocity is
preserved).  If the BRCA fold-enrichment sits inside the distribution of control
fold-enrichments, the reported structures are NOT specific to breast cancer.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import csv, collections, importlib.util, json, random, sys, time
from multiprocessing import Pool
import numpy as np

REV = '/path/to/revision'
LOG = f'{REV}/logs/specificity_enrichment.log'
spec = importlib.util.spec_from_file_location('m08', f'{REV}/scripts/04_motif_significance/08_motif_significance.py')
m08 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m08)

N_RAND_NET  = int(os.environ.get('N_RAND_NET', 100))
SWAP_FACTOR = int(os.environ.get('SWAP_FACTOR', 100))
NPROC       = int(os.environ.get('NPROC', 12))
PER_CLASS   = int(os.environ.get('PER_CLASS', 30))
SEED = 20260908
KEYS = ['miRNA-FFL', 'TF-FFL', 'Composite-FFL', 'TF-TF-FFL', 'Any-3node-FFL',
        '4N-multi-output-FFL', '4N-multi-input-FFL', 'Mutual-TF-miRNA-pairs']

def say(*a):
    m = time.strftime('%H:%M:%S') + ' | ' + ''.join(str(x) for x in a)
    print(m, flush=True)
    with open(LOG, 'a') as fh: fh.write(m + '\n')

def rd(p, sep='\t'):
    with open(p) as fh: return list(csv.DictReader(fh, delimiter=sep))

D = f'{REV}/data/control'
TRR = collections.defaultdict(set)
for r in rd(f'{D}/trrust_pairs.tsv'): TRR[r['TF']].add(r['target'])
TMR = collections.defaultdict(set)
for r in rd(f'{D}/transmir_pairs.tsv'): TMR[r['TF']].add(r['mirna'])
VAL = collections.defaultdict(set)
for r in rd(f'{D}/validated_mirna_target_pairs.tsv'): VAL[r['mirna']].add(r['gene'])

def layers_for(ns):
    tfs, genes, mirs = set(ns['TF']), set(ns['Gene']), set(ns['miRNA'])
    prot = tfs | genes
    nodes = {}
    for v in tfs: nodes[v] = 'TF'
    for v in genes: nodes[v] = 'Gene'
    for v in mirs: nodes[v] = 'miRNA'
    L = collections.defaultdict(list)
    if ns.get('edges') is not None:
        for a, b in ns['edges']:
            t = ('TF_miRNA' if nodes[b] == 'miRNA' and nodes[a] == 'TF'
                 else 'miRNA_target' if nodes[a] == 'miRNA'
                 else 'TF_target')
            L[t].append((a, b))
    else:
        for t in tfs:
            for g in TRR.get(t, ()):
                if g in prot and g != t: L['TF_target'].append((t, g))
            for mi in TMR.get(t, ()):
                if mi in mirs: L['TF_miRNA'].append((t, mi))
        for mi in mirs:
            for g in VAL.get(mi, ()):
                if g in prot: L['miRNA_target'].append((mi, g))
    for k in L: L[k] = list(dict.fromkeys(L[k]))
    return nodes, dict(L)

_S = {}
def _init(payload): _S.update(payload)

def _task(arg):
    label, seed = arg
    nodes, L = _S['nets'][label]
    enc = m08.Encoded(nodes, L)
    groups, eset0 = m08.build_groups(enc, 'NM2')
    UV, _, _ = m08.rewire_groups(groups, eset0, enc.n, SWAP_FACTOR, random.Random(seed), 'NM2')
    c, _ = m08.count_motifs(enc, UV, want_4node=False)
    return label, [c[k] for k in KEYS]

def main():
    open(LOG, 'w').close()
    say('=== 11_specificity_enrichment :: START ===')
    say(f'N_RAND_NET={N_RAND_NET}  SWAP_FACTOR={SWAP_FACTOR}  PER_CLASS={PER_CLASS}  NPROC={NPROC}')
    NS = json.load(open(f'{REV}/cache/specificity_nodesets.json'))
    # the published network keeps its own edge list
    pub = [(r['source'], r['target']) for r in rd(f'{REV}/data/canonical_edges.tsv')
           if r['edge_type'] in ('miRNA_target', 'TF_miRNA', 'TF_target')]
    NS['BRCA_published']['edges'] = pub
    byclass = collections.defaultdict(list)
    for lab, v in NS.items(): byclass[v['cls']].append(lab)
    chosen = []
    for cls, labs in byclass.items():
        labs = sorted(labs)
        chosen += labs if cls.startswith('BRCA_') and len(labs) <= 2 else labs[:PER_CLASS]
    say(f'networks under test: {len(chosen)}  ' +
        str({c: sum(1 for l in chosen if NS[l]['cls'] == c) for c in byclass}))

    nets = {}; real = {}
    for lab in chosen:
        nodes, L = layers_for(NS[lab])
        if sum(len(v) for v in L.values()) < 10: continue
        nets[lab] = (nodes, L)
        enc = m08.Encoded(nodes, L)
        U0 = np.concatenate([u for u, _ in enc.layers]); V0 = np.concatenate([v for _, v in enc.layers])
        c, _ = m08.count_motifs(enc, (U0, V0), want_4node=False)
        real[lab] = ([c[k] for k in KEYS], enc.m)
    say(f'{len(nets)} networks retained (>=10 edges); total edges = '
        f'{sum(v[1] for v in real.values())}')

    tasks = [(lab, SEED + 1 + i) for lab in nets for i in range(N_RAND_NET)]
    say(f'{len(tasks)} randomisation tasks')
    acc = collections.defaultdict(list)
    t0 = time.time()
    with Pool(NPROC, initializer=_init, initargs=({'nets': nets},)) as pool:
        for i, (lab, v) in enumerate(pool.imap_unordered(_task, tasks, chunksize=2), 1):
            acc[lab].append(v)
            if i % 1000 == 0:
                el = time.time() - t0
                say(f'   {i}/{len(tasks)}  elapsed {el/60:.1f} min  '
                    f'eta {el/i*(len(tasks)-i)/60:.1f} min')
    say(f'done in {(time.time()-t0)/60:.1f} min')

    out = []
    for lab in nets:
        A = np.array(acc[lab], float)
        nr, ne = real[lab]
        for j, k in enumerate(KEYS):
            col = A[:, j]; mu = col.mean(); sd = col.std(ddof=1)
            out.append(dict(network=lab, network_class=NS[lab]['cls'], motif_class=k,
                            n_edges=ne, n_real=int(nr[j]), rand_mean=round(float(mu), 3),
                            rand_sd=round(float(sd), 3),
                            fold_enrichment=(round(nr[j] / mu, 4) if mu > 0 else 'Inf'),
                            Z=(round((nr[j] - mu) / sd, 3) if sd > 0 else 'NA'),
                            p_emp=(1 + int((col >= nr[j]).sum())) / (1 + len(col)),
                            n_rand=len(col)))
    with open(f'{REV}/results/motif_specificity_enrichment.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader()
        for r in out: w.writerow(r)
    say(f'wrote results/motif_specificity_enrichment.csv ({len(out)} rows)')

    say('FOLD-ENRICHMENT of Any-3node-FFL over each network OWN randomised ensemble:')
    for cls in ('BRCA_published', 'BRCA_rebuilt', 'BRCA_edgematched', 'CTRL_matched',
                'CTRL_topann', 'CTRL_random'):
        v = [r['fold_enrichment'] for r in out
             if r['network_class'] == cls and r['motif_class'] == 'Any-3node-FFL'
             and r['fold_enrichment'] != 'Inf']
        z = [r['Z'] for r in out if r['network_class'] == cls
             and r['motif_class'] == 'Any-3node-FFL' and r['Z'] != 'NA']
        if not v: continue
        v = np.array(v, float); z = np.array(z, float)
        say(f'   {cls:<18} n={len(v):<4} fold={v.mean():.3f} +-{v.std(ddof=1) if len(v)>1 else 0:.3f} '
            f'[{v.min():.3f},{v.max():.3f}]   Z={z.mean():.2f} +-{z.std(ddof=1) if len(z)>1 else 0:.2f}')
    say('=== 11 DONE ===')

if __name__ == '__main__':
    main()
