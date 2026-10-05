#!/usr/bin/env python3
"""A1 of SETTINGS.md: prevalences and comparisons for the Bhat pattern modules.

Per module (51: 41 run here, 10 identical to the paper's stored modules and read from the stored per-set files):
prevalence (% of the 16,384 sets, 95 % Clopper-Pearson CI) of memory, pulse, damped oscillation, certified sustained
oscillation, oscillation (damped OR certified sustained), memory AND pulse, memory AND oscillation (definitions of
analyses/six_node_pattern S2; certification from a1_certify.py for the run modules and from S2's certification file, which
includes the stored 03f result, for the stored modules).
Comparisons, exact two-sided McNemar (binomial on the discordant sets), paired by parameter set:
  (i)   every six-node configuration against the paper's six-node module (COMP_C2_toggle GG+MM+TT)
  (ii)  every six-node configuration against the same class's core, GG and GG+MM modules with the same s_TM and s_TY = +1
  (iii) TF2 -> TF1 attribution: each six-node TF-FFL or composite configuration with the arc as activation and as
        repression against the arc absent (same class and other signs)
Claim guard: a six-node configuration "exceeds every smaller pattern" for a joint behaviour only if its prevalence is
higher than that of each of its core, GG and GG+MM modules with McNemar p < 0.05.
usage: python3 a1_analyse.py <analysis root> <dynamics folder> <S2 certification csv>
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta

BEH = ['memory', 'pulse', 'damped oscillation', 'certified sustained oscillation', 'oscillation', 'memory AND pulse', 'memory AND oscillation']
STORED_FILE = {'COMP_C2_toggle': 'results/v3/dynamics_higher_order_persetset.csv.gz', 'I1_miRNA_FFL': 'results/v3/dynamics_higher_order_persetset.csv.gz',
               'COMP_I1_negfeedback': 'results/v3/dynamics_higher_order_persetset_compI1.csv.gz'}
LAYER2N = {'core': 'n3', 'GG': 'n4', 'GG+MM': 'n5', 'GG+MM+TT': 'n6'}


def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1); hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return 100 * lo, 100 * hi


def flags(X, cert):
    X = X.sort_values('param_set').reset_index(drop=True)
    cs = X.param_set.isin(cert) & (X.is_sustained_osc > 0.5)
    F = pd.DataFrame({'param_set': X.param_set, 'memory': X.is_memory > 0.5, 'pulse': X.is_pulse > 0.5,
                      'damped oscillation': X.is_damped_osc > 0.5, 'certified sustained oscillation': cs})
    F['oscillation'] = F['damped oscillation'] | F['certified sustained oscillation']
    F['memory AND pulse'] = F.memory & F.pulse
    F['memory AND oscillation'] = F.memory & F.oscillation
    return F


def mcn(a, b):
    b_only = int((a & ~b).sum()); c_only = int((~a & b).sum())
    p = binomtest(b_only, b_only + c_only, 0.5).pvalue if b_only + c_only else 1.0
    return b_only, c_only, p


def main(root, dyn, s2cert):
    M = pd.read_csv(os.path.join(dyn, 'a1_modules.csv'), dtype=str).fillna('')
    certA = pd.read_csv(os.path.join(dyn, 'a1_certification.csv')) if os.path.exists(os.path.join(dyn, 'a1_certification.csv')) else pd.DataFrame(columns=['family', 'module', 'param_set', 'certified_limit_cycle'])
    certS = pd.read_csv(s2cert)
    stored = {}
    F = {}
    for r in M.to_dict('records'):
        key = (r['family'], r['module'])
        if r['stored_as']:
            fam, lay = r['stored_as'].split(' ')
            if STORED_FILE[fam] not in stored: stored[STORED_FILE[fam]] = pd.read_csv(os.path.join(root, STORED_FILE[fam]))
            S = stored[STORED_FILE[fam]]; X = S[(S.family == fam) & (S.module == LAYER2N[lay])]
            cset = set(certS.loc[(certS.family == fam) & (certS.module == lay) & certS.certified_limit_cycle.astype(bool), 'param_set'].astype(int))
            flagged = int((X.is_sustained_osc > 0.5).sum())
            assert flagged == 0 or ((certS.family == fam) & (certS.module == lay)).any(), (fam, lay, 'flagged but not certified in S2')
        else:
            X = pd.read_csv(os.path.join(dyn, 'runs', 'perset', f"{r['family']}__{r['module'].replace('|', '')}.csv.gz"))
            cset = set(certA.loc[(certA.family == r['family']) & (certA.module == r['module']) & certA.certified_limit_cycle.astype(bool), 'param_set'].astype(int))
            flagged = int((X.is_sustained_osc > 0.5).sum())
            assert flagged == 0 or ((certA.family == r['family']) & (certA.module == r['module'])).any(), (key, 'flagged but not certified')
        assert len(X) == 16384
        F[key] = flags(X, cset)
    paper = ('COMP_C2_toggle', 'GG+MM+TT')
    S = stored[STORED_FILE['COMP_C2_toggle']]; X = S[(S.family == 'COMP_C2_toggle') & (S.module == 'n6')]
    F[paper] = flags(X, set(certS.loc[(certS.family == 'COMP_C2_toggle') & (certS.module == 'GG+MM+TT') & certS.certified_limit_cycle.astype(bool), 'param_set'].astype(int)))
    info = {(r['family'], r['module']): r for r in M.to_dict('records')}
    info[paper] = dict(family='COMP_C2_toggle', module='GG+MM+TT', size_key='GG+MM+TT', s_TM='-1', s_TY='1', s_T12='1', s_T21='0', s_T2Y='1', stored_as='COMP_C2_toggle GG+MM+TT (the paper\'s six-node module)')
    rows = []
    for key, f in F.items():
        r = info[key]; d = dict(family=key[0], module=key[1], size=r['size_key'], s_TM=r['s_TM'], s_TY=r['s_TY'], s_T12=r['s_T12'], s_T21=r['s_T21'],
                                s_T2Y=r['s_T2Y'], source=r['stored_as'] or 'run here', n_sets=len(f))
        for b in BEH:
            k = int(f[b].sum()); lo, hi = cp(k, len(f))
            d[f'{b} n'] = k; d[f'{b} %'] = 100 * k / len(f); d[f'{b} CI low'] = lo; d[f'{b} CI high'] = hi
        rows.append(d)
    PR = pd.DataFrame(rows); PR.to_csv(os.path.join(dyn, 'a1_prevalence.csv'), index=False)
    comp, guard = [], []
    six = [k for k in F if info[k]['size_key'] == 'GG+MM+TT' and k != paper]
    for k in six:
        r = info[k]
        smaller = [(kk, info[kk]['size_key']) for kk in F if kk[0] == k[0] and info[kk]['size_key'] in ('core', 'GG', 'GG+MM')]
        for b in BEH:
            bb, cc, p = mcn(F[k][b].values, F[paper][b].values)
            comp.append(dict(comparison='six-node configuration vs the paper\'s six-node module', family=k[0], module=k[1], other=f'{paper[0]} {paper[1]}',
                             behaviour=b, pct_this=100 * F[k][b].mean(), pct_other=100 * F[paper][b].mean(), this_only=bb, other_only=cc, mcnemar_p=p))
            ok_all = True
            for kk, lay in sorted(smaller, key=lambda x: ('core', 'GG', 'GG+MM').index(x[1])):
                bb, cc, p = mcn(F[k][b].values, F[kk][b].values)
                comp.append(dict(comparison=f'six-node configuration vs the same class\'s {lay}', family=k[0], module=k[1], other=f'{kk[0]} {kk[1]}',
                                 behaviour=b, pct_this=100 * F[k][b].mean(), pct_other=100 * F[kk][b].mean(), this_only=bb, other_only=cc, mcnemar_p=p))
                ok_all &= (F[k][b].mean() > F[kk][b].mean()) and (p < 0.05)
            if b in ('memory AND pulse', 'memory AND oscillation'):
                guard.append(dict(family=k[0], module=k[1], behaviour=b, smaller_patterns=', '.join(lay for kk, lay in smaller),
                                  exceeds_every_smaller_pattern_p_lt_005=bool(ok_all and len(smaller) == 3)))
        if r['s_T21'] in ('1', '-1'):
            base = [kk for kk in F if kk[0] == k[0] and info[kk]['size_key'] == 'GG+MM+TT' and info[kk]['s_T21'] == '0'
                    and all(info[kk][x] == r[x] for x in ('s_T12', 's_T2Y'))]
            for kk in base:
                for b in BEH:
                    bb, cc, p = mcn(F[k][b].values, F[kk][b].values)
                    comp.append(dict(comparison=f"TF2->TF1 {'activation' if r['s_T21'] == '1' else 'repression'} vs absent", family=k[0], module=k[1],
                                     other=f'{kk[0]} {kk[1]}', behaviour=b, pct_this=100 * F[k][b].mean(), pct_other=100 * F[kk][b].mean(),
                                     this_only=bb, other_only=cc, mcnemar_p=p))
    pd.DataFrame(comp).to_csv(os.path.join(dyn, 'a1_mcnemar.csv'), index=False)
    pd.DataFrame(guard).to_csv(os.path.join(dyn, 'a1_claim_guard.csv'), index=False)
    txt = (PR[['family', 'module', 'memory %', 'pulse %', 'oscillation %', 'memory AND pulse %', 'memory AND oscillation %']].round(3).to_string(index=False)
           + '\n\n' + pd.DataFrame(guard).to_string(index=False) + '\n')
    open(os.path.join(dyn, 'a1_analyse.log'), 'w').write(txt); print(txt)


if __name__ == '__main__':
    main(*sys.argv[1:4])
