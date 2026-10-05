#!/usr/bin/env python3
"""A1b of SETTINGS.md: the module lists.
  configurations  the nine six-node configurations that passed A1's claim guard for memory AND pulse and count for
                  the verdict (TF2 -> TF1 as activation or repression; miRNA-FFL without the arc), plus option (b)
                  (authors, 2026-09-29 11:25): Bhat_miR_TM0_TY+ T12+ / T21 absent / T2Y+, reported for memory AND
                  oscillation only
  run list        TT, GG+TT and MM+TT of each configuration, with its own TF2 arcs; MM of each family whose smaller
                  modules A1 ran (families whose smaller modules are stored reuse S2's stored MM)
  nested list     S2's stored COMP_I1_negfeedback TT, GG+TT, MM and MM+TT, rebuilt through a1b_run.py's layers as the
                  Bhat composite with TM+, TY+, T12+, T21 absent, T2Y+ (the same model; A1 maps it to COMP_I1_negfeedback)
  smaller map     for each configuration, the source of each of its seven smaller modules
usage: python3 a1b_modules.py <dynamics folder>
Outputs: a1b_configurations.csv, a1b_run_modules.csv, a1b_nested_modules.csv, a1b_smaller_map.csv, a1b_modules.log"""
import os, sys
import pandas as pd

LAYERS_NEW = ['TT', 'GG+TT', 'MM+TT']
S2_MM = {'COMP_C2_toggle': 'COMP_C2_toggle__MM.csv.gz', 'COMP_I1_negfeedback': 'COMP_I1_negfeedback__MM.csv.gz',
         'I1_miRNA_FFL': 'I1_miRNA_FFL__MM.csv.gz'}
FIELDS = ['family', 'module', 'registry', 'size_key', 'mir_to_TF', 'input_node', 's_TM', 's_TY', 's_T12', 's_T21', 's_T2Y', 's_T2M']
SEC_PER_CHUNK = {'TT': 450, 'MM': 450, 'MM+TT': 700, 'GG+TT': 1100}      # SETTINGS.md A1b, estimate


def main(dyn):
    M = pd.read_csv(os.path.join(dyn, 'a1_modules.csv'), dtype=str).fillna('')
    G = pd.read_csv(os.path.join(dyn, 'a1_claim_guard.csv'), dtype=str)
    key = {(r['family'], r['module']): r for r in M.to_dict('records')}
    ok = G[(G.behaviour == 'memory AND pulse') & (G.exceeds_every_smaller_pattern_p_lt_005 == 'True')]
    conf = [dict(key[(f, m)], role='verdict') for f, m in zip(ok.family, ok.module) if key[(f, m)]['s_T21'] != '0' or f.startswith('Bhat_miR_')]
    assert len(conf) == 9, len(conf)
    conf.append(dict(key[('Bhat_miR_TM0_TY+', 'GG+MM+TT_T12+_T210_T2Y+')], role='option (b): memory AND oscillation only'))
    C = pd.DataFrame(conf)[FIELDS + ['role']]; C.to_csv(os.path.join(dyn, 'a1b_configurations.csv'), index=False)
    run, smap, mm_done = [], [], set()
    for c in conf:
        suffix = c['module'].replace('GG+MM+TT_', '')
        for lay in LAYERS_NEW:
            r = {k: c[k] for k in FIELDS}; r.update(module=f'{lay}_{suffix}', size_key=lay, registry=''); run.append(r)
        fam = c['family']
        small = {k: v for k, v in key.items() if k[0] == fam and v['size_key'] in ('core', 'GG', 'GG+MM')}
        stored = {v['size_key']: v['stored_as'] for v in small.values() if v['stored_as']}
        if stored:                                   # the family's smaller modules are stored S2/paper modules
            s2fam = next(iter(stored.values())).split(' ')[0]
            mm_src = f'S2 stored: {S2_MM[s2fam]}'
        else:
            mm_src = 'A1b run: MM'
            if fam not in mm_done:
                core = next(v for v in small.values() if v['size_key'] == 'core')
                r = dict(core); r.update(module='MM', size_key='MM', registry=''); run.append(r); mm_done.add(fam)
        for v in sorted(small.values(), key=lambda v: ('core', 'GG', 'GG+MM').index(v['size_key'])):
            smap.append(dict(family=fam, module=c['module'], role=c['role'], smaller=v['size_key'],
                             source=('stored: ' + v['stored_as']) if v['stored_as'] else f"A1 run: {v['module']}"))
        smap.append(dict(family=fam, module=c['module'], role=c['role'], smaller='MM', source=mm_src))
        for lay in LAYERS_NEW:
            smap.append(dict(family=fam, module=c['module'], role=c['role'], smaller=lay, source=f'A1b run: {lay}_{suffix}'))
    R = pd.DataFrame(run)[FIELDS]; R.to_csv(os.path.join(dyn, 'a1b_run_modules.csv'), index=False)
    pd.DataFrame(smap).to_csv(os.path.join(dyn, 'a1b_smaller_map.csv'), index=False)
    ref = key[('Bhat_comp_TM+_TY+', 'GG+MM+TT_T12+_T210_T2Y+')]
    assert ref['stored_as'] == 'COMP_I1_negfeedback GG+MM+TT', ref['stored_as']
    nested = []
    for lay in ('TT', 'GG+TT', 'MM', 'MM+TT'):
        r = dict(ref); r.update(family='COMP_I1_negfeedback', module=f'{lay}@nested', size_key=lay, registry=''); nested.append(r)
    pd.DataFrame(nested)[FIELDS].to_csv(os.path.join(dyn, 'a1b_nested_modules.csv'), index=False)
    cpu = lambda rows: sum(SEC_PER_CHUNK[r['size_key']] * 8 for r in rows) / 3600
    L = [f'configurations: {len(conf)} ({sum(c["role"] == "verdict" for c in conf)} for the verdict)', f'run modules: {len(run)} ({sum(r["size_key"] == "MM" for r in run)} MM), {8 * len(run)} chunks, '
         f'about {cpu(run):.0f} CPU-h; with 12 workers about {cpu(run) / 12:.1f} h plus the tail',
         f'nested check modules: {len(nested)}, {8 * len(nested)} chunks, about {cpu(nested):.0f} CPU-h']
    open(os.path.join(dyn, 'a1b_modules.log'), 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))


if __name__ == '__main__':
    main(sys.argv[1])
