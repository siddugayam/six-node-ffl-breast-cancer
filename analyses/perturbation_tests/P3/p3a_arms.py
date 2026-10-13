#!/usr/bin/env python3
"""P3a arms: the contrasts read from the P3a hits that reach the 40 largest (SETTINGS P3; see P3a/p3a_triage_note.txt).
Each contrast = one perturbed network miRNA (mature arm suffix removed -> network node) in one cell type against the
pooled controls of the same cell type and condition; latest time point; unstimulated where present.
Writes P3a/p3a_arms.tsv (one row per sample) in the p1_run_dataset.py format."""
import os, sys, re
import pandas as pd
INB = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, INB)
import _geo
N = pd.read_csv(f'{os.path.dirname(os.path.dirname(INB))}/data/canonical_nodes.tsv', sep='\t'); NODES = set(N[N.type == 'miRNA'].name)
def node(mature):
    m = re.sub(r'-(3p|5p)$', '', mature.strip())
    return m if m in NODES else None
rows = []
def ch(s): return dict(c.split(': ', 1) for c in s['characteristics'] if ': ' in c)
def add(gse, label, cell, direction, mirna, pert, ctrl, S, platform=None):
    assert len(pert) >= 2 and len(ctrl) >= 2, (gse, label, len(pert), len(ctrl))
    for arm, gl in (('pert', pert), ('ctrl', ctrl)):
        for g in gl:
            s = S[g]; rows.append(dict(gse=gse, tf=label, gsm=g, arm=arm, title=s['title'], platform=platform or s['platform'], strategy=s['strategy'],
                                       cell=cell, cls='P3a', direction=direction, mirna=mirna))
# GSE115646: HepG2, 48 h, each miRNA mimic vs the 17 mock transfections (miR-101-3p is the P1b contrast and is not repeated)
S = {s['gsm']: s for s in _geo.samples('GSE115646', 'P3a')}
hep = {g: ch(s) for g, s in S.items() if ch(s).get('cell line') == 'HepG2' and s['strategy'] == 'RNA-Seq' and 'treatment' in ch(s)}
mock = [g for g, c in hep.items() if c['treatment'].startswith('mock')]
by = {}
for g, c in hep.items():
    m = re.match(r'(hsa-\S+) transfection', c['treatment'])
    if m: by.setdefault(m.group(1), []).append(g)
for mat, gl in sorted(by.items()):
    nd = node(mat)
    if nd and mat != 'hsa-miR-101-3p' and len(gl) >= 2: add('GSE115646', mat, f'HepG2 ({mat})', 'gain', nd, gl, mock, S)
# GSE246676: DU145 radioresistant (RR), mimics vs miR control + untransfected RR (pooled controls)
S = {s['gsm']: s for s in _geo.samples('GSE246676', 'P3a')}
rr = {g: ch(s) for g, s in S.items() if ch(s).get('phenotype') == 'RR'}
ctrl = [g for g, c in rr.items() if c['transfection'] in ('none', 'miR control')]
for mat in ('miR-200c-3p', 'miR-92a-3p'):
    add('GSE246676', f'hsa-{mat}', f'DU145 radioresistant ({mat})', 'gain', node('hsa-' + mat), [g for g, c in rr.items() if c['transfection'] == mat], ctrl, S)
# GSE277200: hiPSC-derived NSCs and neurons, MIR155 over-expression vs control
S = {s['gsm']: s for s in _geo.samples('GSE277200', 'P3a')}
for stage, lab in (('NSC', 'NSCs'), ('NEURO', 'neurons')):
    add('GSE277200', 'hsa-miR-155', f'hiPSC-derived {lab}', 'gain', 'hsa-miR-155',
        [g for g, s in S.items() if s['title'].startswith(f'MIR155_{stage}')], [g for g, s in S.items() if s['title'].startswith(f'CONTR_{stage}')], S)
# GSE56268: P3HR1, doxycycline-induced miR-28 vs control vector, 24 h (latest)
S = {s['gsm']: s for s in _geo.samples('GSE56268', 'P3a')}
t24 = {g: ch(s) for g, s in S.items() if ch(s).get('time') == '24h'}
add('GSE56268', 'hsa-miR-28', 'P3HR1 Burkitt lymphoma (24 h)', 'gain', 'hsa-miR-28',
    [g for g, c in t24.items() if c.get('genotype/variation') == 'miR28'], [g for g, c in t24.items() if c.get('genotype/variation') == 'control'], S)
# GSE98848: IOMM-Lee, miR-16 vs control (miR-519 is not mapped: the series does not name the family member)
S = {s['gsm']: s for s in _geo.samples('GSE98848', 'P3a')}
add('GSE98848', 'hsa-miR-16', 'IOMM-Lee meningioma', 'gain', 'hsa-miR-16',
    [g for g, s in S.items() if s['title'].startswith('miR16.')], [g for g, s in S.items() if s['title'].startswith('Ctrl.')], S)
# GSE179200: LX-2 hepatic stellate cells, miR-34b + miR-34c mimics vs negative-control mimic + untreated (no TGF-beta)
S = {s['gsm']: s for s in _geo.samples('GSE179200', 'P3a')}
add('GSE179200', 'hsa-miR-34b/c', 'LX-2 hepatic stellate cells (no TGF-beta)', 'gain', 'hsa-miR-34c',
    [g for g, s in S.items() if re.match(r'miR_\d', s['title'])], [g for g, s in S.items() if re.match(r'(NC|NT)_\d', s['title'])], S)
# GSE252885: THP-1 macrophages, LPS-stimulated only (flagged): miR-122-5p mimic vs control miRNA
S = {s['gsm']: s for s in _geo.samples('GSE252885', 'P3a')}
add('GSE252885', 'hsa-miR-122-5p', 'THP-1 macrophages (LPS-stimulated)', 'gain', 'hsa-miR-122',
    [g for g, s in S.items() if 'miR-122-5p' in ch(s).get('treatment', '')], [g for g, s in S.items() if 'control miRNA' in ch(s).get('treatment', '')], S)
# GSE274664 (B cells, miR-99a antagomir) is not used: its only data file (TPM) covers 4 vs 4 samples (8), below the cut-off
# GSE206584: UT-12 cSCC cells, lentiviral miR-23b over-expression vs scrambled control
S = {s['gsm']: s for s in _geo.samples('GSE206584', 'P3a')}
add('GSE206584', 'hsa-miR-23b', 'UT-12 cutaneous SCC', 'gain', 'hsa-miR-23b',
    [g for g, s in S.items() if '23b OE' in s['title']], [g for g, s in S.items() if re.search(r', Ctl\d', s['title'])], S)
# GSE46039: Flp-in T-REx 293, miR-92a inhibitor vs control inhibitor.  Human Exon 1.0 ST arrays (GPL5188): RMA to
# transcript clusters (core), so the transcript-cluster annotation GPL5175 maps them to genes
S = {s['gsm']: s for s in _geo.samples('GSE46039', 'P3a')}
add('GSE46039', 'hsa-miR-92a', 'HEK293 (Flp-in T-REx)', 'loss', 'hsa-miR-92a',
    [g for g, s in S.items() if 'miR-92a inhibitor' in ch(s).get('treatment', '')], [g for g, s in S.items() if 'control inhibitor' in ch(s).get('treatment', '')], S,
    platform='GPL5175')
# GSE121892: primary osteoblasts (5 patients), miR-320a mimic vs control mimic; miR-320a inhibitor vs control inhibitor
S = {s['gsm']: s for s in _geo.samples('GSE121892', 'P3a')}
tr = {g: ch(s).get('treatment', '') for g, s in S.items()}
add('GSE121892', 'hsa-miR-320a mimic', 'primary osteoblasts (mimic)', 'gain', 'hsa-miR-320a',
    [g for g, t in tr.items() if t.endswith('miR-320a mimic (M)')], [g for g, t in tr.items() if t.endswith('control mimic (CM)')], S)
add('GSE121892', 'hsa-miR-320a inhibitor', 'primary osteoblasts (inhibitor)', 'loss', 'hsa-miR-320a',
    [g for g, t in tr.items() if t.endswith('miR-320a inhibitor (I)')], [g for g, t in tr.items() if t.endswith('control inhibitor (CI)')], S)
# GSE133527: M2 macrophages, non-stimulated: miR-221 mimic vs control
S = {s['gsm']: s for s in _geo.samples('GSE133527', 'P3a')}
ns = {g: ch(s) for g, s in S.items() if ch(s).get('stimulation') == 'non-stimulated'}
add('GSE133527', 'hsa-miR-221', 'M2 macrophages (non-stimulated)', 'gain', 'hsa-miR-221',
    [g for g, c in ns.items() if 'miR-221' in c['transfection']], [g for g, c in ns.items() if c['transfection'] == 'control'], S)
# GSE196161: HUVEC, TNF-stimulated only (flagged): miR-125a-5p vs negative control
S = {s['gsm']: s for s in _geo.samples('GSE196161', 'P3a')}
add('GSE196161', 'hsa-miR-125a-5p', 'HUVEC (TNF-stimulated)', 'gain', 'hsa-miR-125a',
    [g for g, s in S.items() if s['title'].startswith('HUVEC_miR-125a')], [g for g, s in S.items() if s['title'].startswith('HUVEC_NC')], S)
A = pd.DataFrame(rows); A.to_csv(f'{INB}/P3a/p3a_arms.tsv', sep='\t', index=False)
print(A.groupby(['gse', 'tf', 'arm']).size().unstack().to_string())
