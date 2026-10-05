#!/usr/bin/env python3
"""A3 of SETTINGS.md, item 2: the twelve Bhat networks under NULL-B, NULL-C and NULL-L on two restricted graphs, with
section O's program (bhat_null.c) and settings (1,000 replicates, seed 20250908, 100 swap attempts per rewired item):
  physical   gene-gene links restricted to STRING v12.0 physical score >= 0.900 (analyses/six_node_pattern S7(b) inputs)
  validated  miRNA -> target edges of the strong and weak tiers only, exemplar edges removed (S7(c) inputs)
Check: the observed counts equal S7's reported counts where S7 reports them. The observed instance lists are kept, to
rebuild the restricted networks. Statistics and decision rule as section O.

usage: python3 evidence_nulls.py <S7 folder> <output folder>
"""
import os, sys, shutil, hashlib, subprocess
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bhat_common as bc

R, SEED, SPE = 1000, 20250908, 100
NULLS = [('NULL-B', 'B'), ('NULL-C', 'C'), ('NULL-L', 'L')]
CK = {'miRNA_FFL': 'mirFFL', 'TF_FFL': 'TFFFL', 'composite_FFL': 'comp'}
VAR = {'physical': ('graph_physical_ge_0.900.txt', 'labels_physical_ge_0.900.txt'),
       'validated': ('null_validated.txt', 'labels_validated.txt')}


def col(net, what): return f"bhat{net[0]}_{CK[net.split('_', 1)[1]]}_{what}"


def main(s7, out):
    log = []; say = lambda *a: (print(*a, flush=True), log.append(' '.join(map(str, a))))
    os.makedirs(out, exist_ok=True)
    exe = os.path.join(out, 'bhat_null'); subprocess.run(['gcc', '-O2', '-o', exe, os.path.join(HERE, 'bhat_null.c')], check=True)
    ref = {'physical': pd.read_csv(os.path.join(s7, 's7b_restricted_gene_gene_counts.csv')).set_index('restriction').loc['physical_ge_0.900'],
           'validated': pd.read_csv(os.path.join(s7, 's7c_validated_only_summary.csv')).iloc[0]}
    allrows, alldec = [], []
    for v, (gf, lf) in VAR.items():
        d = os.path.join(out, v); os.makedirs(os.path.join(d, 'inputs'), exist_ok=True); os.makedirs(os.path.join(d, 'obs'), exist_ok=True)
        os.makedirs(os.path.join(d, 'runs'), exist_ok=True)
        for f in (gf, lf): shutil.copy2(os.path.join(s7, f), os.path.join(d, 'inputs', f))
        G, L = os.path.join(d, 'inputs', gf), os.path.join(d, 'inputs', lf)
        say(f'{v}: inputs {gf} md5 {hashlib.md5(open(G, "rb").read()).hexdigest()}, {lf} md5 {hashlib.md5(open(L, "rb").read()).hexdigest()}')
        with open(os.path.join(d, 'obs', 'observed.tsv'), 'w') as f:
            subprocess.run([exe, G, L, 'OBS', '0', str(SEED), str(SPE), '0', '1', '0', os.path.join(d, 'obs', 'observed')], stdout=f, check=True)
        for fn in os.listdir(os.path.join(d, 'obs')):
            if 'census6c' in fn or 'model6' in fn: os.remove(os.path.join(d, 'obs', fn))
        O = pd.read_csv(os.path.join(d, 'obs', 'observed.tsv'), sep='\t').iloc[0]
        chk = {c: (int(O[c]), int(ref[v][c])) for c in ref[v].index if c in O.index and c.startswith(('bhat', 'model6'))}
        ok = all(a == b for a, b in chk.values())
        say(f'CHECK {v}: observed counts equal S7\'s for ' + ', '.join(f'{c} {a}' for c, (a, b) in chk.items()) + f': {ok}')
        assert ok
        say(f'{v} observed: ' + ', '.join(f'{net} {int(O[col(net, "inst")])}/{int(O[col(net, "sets")])}' for net in bc.NETWORKS))
        procs = []
        for nm, m in NULLS:
            fo = open(os.path.join(d, 'runs', f'null_{m}.tsv'), 'w'); fe = open(os.path.join(d, 'runs', f'null_{m}.err'), 'w')
            procs.append((m, subprocess.Popen([exe, G, L, m, str(R), str(SEED), str(SPE), '0', '1', '0'], stdout=fo, stderr=fe), fo, fe))
        for m, p, fo, fe in procs:
            assert p.wait() == 0, (v, m); fo.close(); fe.close()
        rows = []
        for nm, m in NULLS:
            X = pd.read_csv(os.path.join(d, 'runs', f'null_{m}.tsv'), sep='\t'); assert len(X) == R
            for net in bc.NETWORKS:
                for w, lab in (('inst', 'instances'), ('sets', 'node sets')):
                    c = col(net, w); x = X[c].to_numpy(float); o = float(O[c]); mu = x.mean(); sd = x.std(ddof=1)
                    rows.append(dict(variant=v, network=net, count=lab, null=nm, observed=o, null_mean=mu, null_sd=sd,
                                     obs_over_mean=o / mu if mu > 0 else np.nan, Z=(o - mu) / sd if sd > 0 else np.nan,
                                     p_upper=(np.sum(x >= o) + 1) / (R + 1), p_lower=(np.sum(x <= o) + 1) / (R + 1), replicates=len(x)))
        T = pd.DataFrame(rows); allrows.append(T)
        for net in bc.NETWORKS:
            for lab in ('instances', 'node sets'):
                g = lambda nm: T[(T.network == net) & (T['count'] == lab) & (T.null == nm)].iloc[0]
                b, c, l = g('NULL-B'), g('NULL-C'), g('NULL-L')
                if c.observed == 0: dd = 'not tested (no instance on the observed graph)'
                elif c.p_upper < 0.05 and l.p_upper < 0.05: dd = 'over-represented'
                elif c.p_upper < 0.05: dd = 'explained by the three-node cores'
                elif c.p_lower < 0.05 and l.p_lower < 0.05: dd = 'depleted'
                else: dd = 'none of the three (not over-represented, not depleted)'
                alldec.append(dict(variant=v, network=net, count=lab, observed=c.observed, NULL_B_ratio=b.obs_over_mean, NULL_B_p_upper=b.p_upper,
                                   NULL_C_ratio=c.obs_over_mean, NULL_C_p_upper=c.p_upper, NULL_C_p_lower=c.p_lower,
                                   NULL_L_ratio=l.obs_over_mean, NULL_L_p_upper=l.p_upper, NULL_L_p_lower=l.p_lower, decision=dd))
    pd.concat(allrows).to_csv(os.path.join(out, 'evidence_null_summary.csv'), index=False)
    D = pd.DataFrame(alldec); D.to_csv(os.path.join(out, 'evidence_decision.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    say('\n' + D.to_string(index=False))
    os.remove(exe)
    open(os.path.join(out, 'evidence_nulls.log'), 'w').write('\n'.join(log) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:3])
