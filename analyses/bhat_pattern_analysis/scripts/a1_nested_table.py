#!/usr/bin/env python3
"""A1 nested check, table (SETTINGS.md, Changes and decisions, 2026-09-28 20:05). For each of the paper's six-node
behaviours (analyses/six_node_pattern S2's list, the Results 3.7 classifiers and joint behaviours), per stored module
(COMP_C2_toggle and COMP_I1_negfeedback, GG+MM+TT):
  paper value   S2's s2_prevalence_gain_loss.csv, column pct (its line given; header = line 1)
  stored        the same prevalence recomputed from the stored per-set file and the stored certification (must equal)
  nested        the nested run (extended model, TF2 -> TF1 absent) and its certification (a1_nested_certify.py)
  difference    nested - paper, in percentage points; sets whose flag differs between stored and nested
Tolerance (SETTINGS.md A1, nested check): every per-set flag equal, i.e. 0 sets differing.
usage: python3 a1_nested_table.py <analysis root> <nested perset folder> <nested certification csv> <S2 folder> <output csv>"""
import os, sys
import pandas as pd
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a1_nested_certify import stored_certification

STORED = {'COMP_C2_toggle': 'results/v3/dynamics_higher_order_persetset.csv.gz',
          'COMP_I1_negfeedback': 'results/v3/dynamics_higher_order_persetset_compI1.csv.gz'}
BEH = [('memory', 'is_memory'), ('ultrasensitivity', 'is_ultrasensitive'), ('bistability', 'is_bistable'),
       ('pulse generation', 'is_pulse'), ('noise rejection', 'is_noise_rejecting'), ('damped oscillation', 'is_damped_osc')]


def behaviours(X, cert):                                  # S2's definitions (s2_analyse.py lines 41-52)
    X = X.sort_values('param_set').reset_index(drop=True)
    b = pd.DataFrame({name: X[col] > 0.5 for name, col in BEH})
    b['certified sustained oscillation'] = X.param_set.isin(cert) & (X.is_sustained_osc > 0.5)
    b['flagged sustained oscillation (sweep)'] = X.is_sustained_osc > 0.5
    b['memory AND pulse'] = b['memory'] & b['pulse generation']
    b['memory AND damped oscillation'] = b['memory'] & b['damped oscillation']
    b['memory AND certified sustained oscillation'] = b['memory'] & b['certified sustained oscillation']
    b['memory AND either oscillation'] = b['memory'] & (b['damped oscillation'] | b['certified sustained oscillation'])
    return b


def main(root, folder, ncert, s2, out):
    P = pd.read_csv(os.path.join(s2, 's2_prevalence_gain_loss.csv'))
    P['line'] = P.index + 2
    st = stored_certification(root, os.path.join(s2, 's2_limit_cycle_certification.csv'))
    NC = pd.read_csv(ncert)
    rows = []
    for fam, rel in STORED.items():
        S = pd.read_csv(os.path.join(root, rel)); S = S[(S.family == fam) & (S.module == 'n6')]
        X = pd.read_csv(os.path.join(folder, f'{fam}__GG+MM+TT@nested.csv.gz'))
        assert len(S) == len(X) == 16384
        bs = behaviours(S, set(st[fam].loc[st[fam].certified_limit_cycle, 'param_set']))
        bn = behaviours(X, set(NC.loc[(NC.family == fam) & NC.certified_limit_cycle.astype(bool), 'param_set'].astype(int)))
        for b in bs.columns:
            p = P[(P.family == fam) & (P.module == 'GG+MM+TT') & (P.behaviour == b)]
            assert len(p) == 1, (fam, b)
            p = p.iloc[0]
            diff = int((bs[b] != bn[b]).sum())
            rows.append(dict(family=fam, module='GG+MM+TT', behaviour=b, paper_n=int(p.n_sets), paper_pct=p.pct,
                             paper_source=f's2_prevalence_gain_loss.csv line {p.line}, column pct',
                             stored_recomputed_pct=100 * bs[b].mean(), nested_n=int(bn[b].sum()), nested_pct=100 * bn[b].mean(),
                             abs_difference_pct_points=abs(100 * bn[b].mean() - p.pct), sets_differing=diff,
                             tolerance='0 sets differing', result='PASS' if diff == 0 else 'FAIL'))
            assert abs(rows[-1]['stored_recomputed_pct'] - p.pct) < 1e-9, (fam, b, 'stored per-set file does not reproduce S2')
    R = pd.DataFrame(rows); R.to_csv(out, index=False)
    print(R[['family', 'behaviour', 'paper_pct', 'nested_pct', 'abs_difference_pct_points', 'sets_differing', 'result']].to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:6])
