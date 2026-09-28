#!/usr/bin/env python3
"""08b -- convergence check on the Maslov-Sneppen mixing time.
Re-estimates the null means at several swap factors.  If 100x|E| is enough, the null means
must not drift as the swap factor is increased.  Writes results/motif_mixing_check.csv;
touches none of the deliverables."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import csv, importlib.util, random, sys, time
from multiprocessing import Pool
import numpy as np
REV = '/path/to/revision'
spec = importlib.util.spec_from_file_location('m08', f'{REV}/scripts/08_motif_significance.py')
m08 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m08)
KEYS = ['miRNA-FFL','TF-FFL','Composite-FFL','TF-TF-FFL','Any-3node-FFL',
        '4N-multi-output-FFL','4N-multi-input-FFL','Mutual-TF-miRNA-pairs']
NPROC = int(os.environ.get('NPROC', 3)); NR = int(os.environ.get('NR', 100))
_S = {}
def _init(nodes, layers, nm):
    enc = m08.Encoded(nodes, layers)
    _S['enc'] = enc; _S['g'] = m08.build_groups(enc, nm); _S['nm'] = nm
def _t(arg):
    sf, seed = arg
    enc = _S['enc']; groups, eset0 = _S['g']
    UV, _, _ = m08.rewire_groups(groups, eset0, enc.n, sf, random.Random(seed), _S['nm'])
    c, _ = m08.count_motifs(enc, UV, want_4node=False)
    return sf, [c[k] for k in KEYS]
if __name__ == '__main__':
    nodes, layers = m08.load_graph()
    out = []
    for nm in ('NM1', 'NM2'):
        tasks = [(sf, 777000 + i) for sf in (25, 100, 400) for i in range(NR)]
        acc = {}
        t0 = time.time()
        with Pool(NPROC, initializer=_init, initargs=(nodes, layers, nm)) as p:
            for sf, v in p.imap_unordered(_t, tasks, chunksize=2):
                acc.setdefault(sf, []).append(v)
        print(f'{nm}: {len(tasks)} randomisations in {(time.time()-t0)/60:.1f} min', flush=True)
        for sf in sorted(acc):
            A = np.array(acc[sf], float)
            for j, k in enumerate(KEYS):
                out.append(dict(null_model=nm, swap_factor=sf, motif_class=k, n_rand=len(A),
                                rand_mean=round(float(A[:, j].mean()), 3),
                                rand_sd=round(float(A[:, j].std(ddof=1)), 3)))
            print(f'  {nm} sf={sf:<4} ' + '  '.join(
                f'{k}={A[:,j].mean():.1f}' for j, k in enumerate(KEYS)), flush=True)
    with open(f'{REV}/results/motif_mixing_check.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader()
        for r in out: w.writerow(r)
    print('wrote results/motif_mixing_check.csv')
