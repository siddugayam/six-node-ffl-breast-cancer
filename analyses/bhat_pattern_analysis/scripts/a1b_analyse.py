#!/usr/bin/env python3
"""A1b of SETTINGS.md: the claim guard against the full layer factorial.
For each A1b configuration (a1b_configurations.csv: the nine of the verdict and option (b)) and each joint behaviour
(memory AND pulse; memory AND oscillation), its prevalence against each of the seven smaller modules of its family and
signs (a1b_smaller_map.csv: core, GG and GG+MM from A1, run or stored; MM from S2 or A1b; TT, GG+TT and MM+TT from A1b),
exact two-sided McNemar paired by parameter set. Extended guard: higher than all seven, each at p < 0.05.
Cells and verdict (SETTINGS.md A1b: rule written 2026-09-29 11:20, confirmed by the authors 11:25): five cells on
memory AND pulse (miRNA-FFL; TF-FFL and composite with TF2 -> TF1 as activation and as repression); a cell passes if one
of its verdict configurations passes. Instance coverage of every configuration: a1_config_coverage.csv (compatible and
exact counts). Behaviour flags and certified sustained oscillation are those of a1_analyse.py.
usage: python3 a1b_analyse.py <analysis root> <dynamics folder> <S2 folder>
Outputs (dynamics folder): a1b_prevalence.csv, a1b_mcnemar.csv, a1b_claim_guard.csv, a1b_cells.csv, a1b_analyse.log"""
import os, sys
import numpy as np, pandas as pd
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_analyse import flags, mcn, cp, BEH, STORED_FILE, LAYER2N

JOINT = ('memory AND pulse', 'memory AND oscillation')
ORDER = ('core', 'GG', 'MM', 'TT', 'GG+MM', 'GG+TT', 'MM+TT')
CLASS = {'Bhat_miR_': ('miRNA-FFL', '6node_miRNA_FFL'), 'Bhat_TF_': ('TF-FFL', '6node_TF_FFL'), 'Bhat_comp_': ('composite', '6node_composite_FFL')}
SIGN = {'0': 'TF2 -> TF1 not in the pattern', '1': 'TF2 -> TF1 activation', '-1': 'TF2 -> TF1 repression'}


def cls_of(fam):
    return next(v for k, v in CLASS.items() if fam.startswith(k))


def main(root, dyn, s2):
    CF = pd.read_csv(os.path.join(dyn, 'a1b_configurations.csv'), dtype=str).fillna('')
    SM = pd.read_csv(os.path.join(dyn, 'a1b_smaller_map.csv'), dtype=str)
    COV = pd.read_csv(os.path.join(dyn, 'a1_config_coverage.csv'), dtype=str).fillna('')
    S2C = pd.read_csv(os.path.join(s2, 's2_limit_cycle_certification.csv'))
    A1C = pd.read_csv(os.path.join(dyn, 'a1_certification.csv'))
    BC = pd.read_csv(os.path.join(dyn, 'a1b_certification.csv'))
    cache, stored = {}, {}

    def certset(T, fam, mod):
        return set(T.loc[(T.family == fam) & (T.module == mod) & T.certified_limit_cycle.astype(bool), 'param_set'].astype(int))

    def certified_module(T, fam, mod):
        return bool(((T.family == fam) & (T.module == mod)).any())

    def load(fam, source):
        """per-set behaviour flags of one module, by its source label in a1b_smaller_map.csv"""
        if (fam, source) in cache: return cache[(fam, source)]
        kind, what = source.split(': ', 1)
        if kind == 'stored':                                        # the paper's stored module (results/v3), S2 certification
            sfam, lay = what.split(' ')
            if STORED_FILE[sfam] not in stored: stored[STORED_FILE[sfam]] = pd.read_csv(os.path.join(root, STORED_FILE[sfam]))
            S = stored[STORED_FILE[sfam]]; X = S[(S.family == sfam) & (S.module == LAYER2N[lay])]; T, tf, tm = S2C, sfam, lay
        elif kind == 'S2 stored':                                   # S2's stored MM
            X = pd.read_csv(os.path.join(s2, 'perset', what)); T, tf, tm = S2C, what.split('__')[0], 'MM'
        elif kind == 'A1 run':
            X = pd.read_csv(os.path.join(dyn, 'runs', 'perset', f'{fam}__{what}.csv.gz')); T, tf, tm = A1C, fam, what
        elif kind == 'A1b run':
            X = pd.read_csv(os.path.join(dyn, 'a1b', 'runs', 'perset', f'{fam}__{what}.csv.gz')); T, tf, tm = BC, fam, what
        else:
            raise ValueError(source)
        assert len(X) == 16384, (fam, source, len(X))
        cs = certset(T, tf, tm)
        assert int((X.is_sustained_osc > 0.5).sum()) == 0 or certified_module(T, tf, tm), (fam, source, 'flagged but not certified')
        cache[(fam, source)] = flags(X, cs)
        return cache[(fam, source)]

    prev, comp, guard = [], [], []
    for c in CF.to_dict('records'):
        fam, mod = c['family'], c['module']
        six = load(fam, f'A1 run: {mod}')
        small = {r['smaller']: r['source'] for r in SM[(SM.family == fam) & (SM.module == mod)].to_dict('records')}
        assert set(small) == set(ORDER), (fam, mod, small)
        for lay in ORDER:
            if small[lay].startswith('A1b run'):
                f = load(fam, small[lay]); d = dict(family=fam, module=small[lay].split(': ')[1], size=lay)
                for b in BEH:
                    k = int(f[b].sum()); lo, hi = cp(k, len(f)); d.update({f'{b} n': k, f'{b} %': 100 * k / len(f), f'{b} CI low': lo, f'{b} CI high': hi})
                prev.append(d)
        for b in JOINT:
            ok_all, fails = True, []
            for lay in ORDER:
                o = load(fam, small[lay]); bb, cc, p = mcn(six[b].values, o[b].values)
                higher = six[b].mean() > o[b].mean()
                comp.append(dict(family=fam, module=mod, role=c['role'], behaviour=b, smaller=lay, source=small[lay],
                                 pct_this=100 * six[b].mean(), pct_other=100 * o[b].mean(), this_only=bb, other_only=cc, mcnemar_p=p))
                if not (higher and p < 0.05): ok_all = False; fails.append(lay)
            guard.append(dict(family=fam, module=mod, role=c['role'], s_T21=c['s_T21'], behaviour=b, pct=100 * six[b].mean(), n_sets=int(six[b].sum()),
                              exceeds_all_seven_p_lt_005=ok_all, not_exceeded=', '.join(fails)))
    P = pd.DataFrame(prev).drop_duplicates(['family', 'module']); P.to_csv(os.path.join(dyn, 'a1b_prevalence.csv'), index=False)
    pd.DataFrame(comp).to_csv(os.path.join(dyn, 'a1b_mcnemar.csv'), index=False)
    G = pd.DataFrame(guard); G.to_csv(os.path.join(dyn, 'a1b_claim_guard.csv'), index=False)

    def coverage(fam, mod):
        cname, net = cls_of(fam); c = CF[(CF.family == fam) & (CF.module == mod)].iloc[0]
        for r in COV[COV.network == net].to_dict('records'):
            if all(r[a] == '' or r[a] == c[a] for a in ('s_TM', 's_T12', 's_T2Y')):
                return int(r['compatible']), int(r['exact']), int(r['instances']), int(r['config'])
        raise KeyError((fam, mod))
    cells = []
    for pre, (cname, net) in CLASS.items():
        for sg in (('0',) if pre == 'Bhat_miR_' else ('1', '-1')):
            for b in JOINT:
                g = G[G.family.str.startswith(pre) & (G.s_T21 == sg) & (G.behaviour == b)
                      & ((G.role == 'verdict') if b == 'memory AND pulse' else True)]
                okc = g[g.exceeds_all_seven_p_lt_005]
                det = []
                for r in okc.to_dict('records'):
                    comp_n, exact_n, tot, k = coverage(r['family'], r['module'])
                    det.append(f"{r['family'].replace('Bhat_', '')} {r['module'].replace('GG+MM+TT_', '')}: {r['pct']:.2f} %, {r['n_sets']} sets; "
                               f"configuration #{k}: {comp_n:,} of {tot:,} instances compatible, {exact_n:,} with exactly its signs")
                cells.append(dict(cell=f'{cname}, {SIGN[sg]}', behaviour=b, configurations=len(g), passing=len(okc), passing_detail=' | '.join(det)))
    CE = pd.DataFrame(cells); CE.to_csv(os.path.join(dyn, 'a1b_cells.csv'), index=False)
    mp = CE[CE.behaviour == 'memory AND pulse']; rep = int((mp.passing > 0).sum())
    verdict = 'supports' if rep == len(mp) else ('qualifies' if rep else 'not reproduced on the Bhat patterns')
    txt = (G.round(3).to_string(index=False) + '\n\n' + CE.to_string(index=False) + f'\n\nVERDICT (memory AND pulse, {rep} of {len(mp)} cells pass): {verdict}\n')
    open(os.path.join(dyn, 'a1b_analyse.log'), 'w').write(txt); print(txt)


if __name__ == '__main__':
    main(*sys.argv[1:4])
