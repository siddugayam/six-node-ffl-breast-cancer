#!/usr/bin/env python3
"""DATASETS.tsv: every dataset found, per part, with its status and reason (SETTINGS.md, "Every hit is listed").
One row per GEO series (or per contrast where a series gave several), per KnockTF dataset (P3b) and per atlas (P5).
Sources: the search hits (<part>_hits.tsv), the triage files, the re-reading reasons (P1/reread_reasons.tsv,
P2/p2_reasons.tsv, P3a/p3a_triage.tsv), the COL1A1 screens and the per-contrast result tables."""
import os, glob, json, re
import numpy as np, pandas as pd
INB = os.path.dirname(os.path.abspath(__file__))
rows = []
def add(**k): rows.append(k)
def meta_index(part):
    out = {}
    for mp in glob.glob(f'{INB}/{part}/de/*.meta.json'):
        m = json.load(open(mp)); out.setdefault((m['gse'], m['tf'], m['cell'], m.get('tag', 'primary')), []).append((mp, m))
    return out
RR = pd.read_csv(f'{INB}/P1/reread_reasons.tsv', sep='\t')
REREAD = {(r.part, r.gse): r.reason for r in RR.itertuples()}
NOMENTION_P1A = 'no human sample names ETS1, NFKB1 (p50/p105), RELA (p65) or SP1; the TF appears only in the series text'
NOMENTION_P1B = 'no human sample names a miR-29 or miR-101 perturbation; the miRNA appears only in the series text'

# ------------------------------------------------------------------ P1a
H = pd.read_csv(f'{INB}/P1a/P1a_hits.tsv', sep='\t')
T = pd.read_csv(f'{INB}/P1a/p1a_triage.tsv', sep='\t')
T2 = pd.read_csv(f'{INB}/P1a/p1a_triage_reread.tsv', sep='\t')
SCR = pd.concat([pd.read_csv(f'{INB}/P1a/p1a_col1a1_screen.tsv', sep='\t'), pd.read_csv(f'{INB}/P1a/p1a_col1a1_screen_reread.tsv', sep='\t')])
SC = {(r.gse, r.tf, r.cell): r for r in SCR.itertuples()}
CON = pd.read_csv(f'{INB}/P1a/p1a_contrasts.tsv', sep='\t')
CC = {(r.gse, r.tf, r.cell, r.tag): r for r in CON.itertuples()}
title = dict(zip(H.gse, H.title))
done = set()
def p1a_candidate(t, src):
    k = (t.gse, t.tf, t.cell); s = SC.get(k)
    base = dict(part='P1a', source=src, accession=t.gse, perturbation=f'{t.tf} knockdown', cell=t.cell, cls=t['class'] if isinstance(t, pd.Series) else t.cls,
                title=title.get(t.gse, ''))
    if s is None:
        add(**base, status='excluded', reason='candidate without a COL1A1 screen row (arms not assigned)'); return
    base.update(n_pert=s.n_pert, n_ctrl=s.n_ctrl, data=s.data)
    if s.col1a1_pass == 'no':
        add(**base, status='excluded', reason=f'COL1A1 not expressed in the controls ({s.col1a1_ctrl}; rule: {s.col1a1_rule})'); return
    if s.col1a1_pass in ('NOT CHECKED', 'NOT FOUND', 'ERROR'):
        why = {'NOT CHECKED': 'RNA-seq without NCBI-generated counts; the submitter files were not processed in this round',
               'NOT FOUND': f'COL1A1 could not be checked: {s.note}', 'ERROR': f'arms not found in the NCBI count table ({str(s.note)[:80]})'}[s.col1a1_pass]
        add(**base, status='not processed', reason=why); return
    for tag in ('primary', 'secondary', 'descriptive'):
        c = CC.get((t.gse, t.tf, t.cell, tag))
        if c is None: continue
        st = ('included' if c.eligible == 'yes' else 'excluded') if tag == 'primary' else f'reported ({tag})'
        why = (f'knockdown verified (TF log2FC {c.TF_log2FC:+.2f}) and COL1A1 expressed ({c.COL1A1_ctrl})' if c.eligible == 'yes'
               else f'knockdown not verified (TF log2FC {c.TF_log2FC:+.2f} > -0.74)' if c.kd_verified == 'no' else f'COL1A1 in controls {c.COL1A1_ctrl}')
        flags = '; '.join(x for x in ('processed' if c.route in ('matrix', 'matrix2') else '', 'stimulated' if 'stimulat' in str(t.reason).lower() or 'induced' in t.cell else '',
                                      'CRISPR knockout' if 'CRISPR' in str(t.reason) else '', 'parental control' if 'parental' in str(t.reason) else '') if x)
        add(**dict(base, n_pert=c.n_pert, n_ctrl=c.n_ctrl, data=c.route), status=st, reason=why, flags=flags,
            result_file=f'P1a/de/{c.gse}__{re.sub(r"[^A-Za-z0-9.-]+", "_", c.tf)}__{re.sub(r"[^A-Za-z0-9]+", "_", c.cell)[:40]}__{tag}.tsv')
    if not any(CC.get((t.gse, t.tf, t.cell, tag)) is not None for tag in ('primary', 'secondary', 'descriptive')):
        add(**base, status='not processed', reason='passed the COL1A1 screen but no contrast was computed')
for _, t in pd.concat([T[T.decision != 'check'], T2]).iterrows():
    done.add(t.gse)
    if (('P1a', t.gse) in REREAD) and t.decision == 'exclude' and t.gse not in set(T2.gse):
        pass
    if t.decision == 'exclude':
        add(part='P1a', source='GEO search', accession=t.gse, perturbation=t.tf, cell=t.cell, cls=t['class'], title=title.get(t.gse, ''), status='excluded', reason=t.reason)
    else:
        p1a_candidate(t, 'GEO search' if t.gse not in set(T2.gse) else 'GEO search (second reading)')
rf = lambda c: f'P1a/de/{c.gse}__{re.sub(r"[^A-Za-z0-9.-]+", "_", c.tf)}__{re.sub(r"[^A-Za-z0-9]+", "_", c.cell)[:40]}__{c.tag}.tsv'
listed = {r.get('result_file') for r in rows if r.get('part') == 'P1a'}
for c in CON.itertuples():                  # secondary (stimulated) and descriptive contrasts, listed with their own cell label
    if rf(c) in listed: continue
    add(part='P1a', source='GEO search', accession=c.gse, perturbation=f'{c.tf} knockdown', cell=c.cell, cls=c.cls, n_pert=c.n_pert, n_ctrl=c.n_ctrl, data=c.route,
        status=f'reported ({c.tag})' if c.tag != 'primary' else ('included' if c.eligible == 'yes' else 'excluded'),
        reason=f'TF log2FC {c.TF_log2FC:+.2f}; COL1A1 in controls {c.COL1A1_ctrl}' + ('; not part of the pools' if c.tag != 'primary' else ''), result_file=rf(c),
        flags='stimulated' if c.tag == 'secondary' else '')
for g in H.gse:
    if g in done: continue
    if ('P1a', g) in REREAD: add(part='P1a', source='GEO search', accession=g, title=title.get(g, ''), status='excluded', reason=REREAD[('P1a', g)])
    else: add(part='P1a', source='GEO search', accession=g, title=title.get(g, ''), status='excluded', reason=NOMENTION_P1A)
add(part='P1a', source='ENCODE portal', accession='ENCSR052ZPX', perturbation='RELA CRISPRi', cell='K562', cls='other', status='not processed',
    reason='ENCODE gene quantifications were not downloaded in this round; K562 is an "other" cell type (descriptive pool only)')
add(part='P1a', source='ENCODE portal', accession='ENCSR650QAK', perturbation='SP1 CRISPRi', cell='K562', cls='other', status='not processed',
    reason='ENCODE gene quantifications were not downloaded in this round; K562 is an "other" cell type (descriptive pool only)')
add(part='P1a', source='LINCS L1000 (GSE92742)', accession='GSE92742', perturbation='TF shRNA signatures', status='not used',
    reason='optional in the analysis plan; COL1A1 is a measured landmark and COL3A1 is not; of the four TFs only ETS1 is a landmark, so NFKB1, RELA and SP1 knockdowns could not be verified without inferred genes; not run')

# ------------------------------------------------------------------ P1b
H = pd.read_csv(f'{INB}/P1b/P1b_hits.tsv', sep='\t'); title = dict(zip(H.gse, H.title))
T = pd.read_csv(f'{INB}/P1b/p1b_triage.tsv', sep='\t')
CON = pd.read_csv(f'{INB}/P1b/p1b_contrasts.tsv', sep='\t')
done = set()
for _, t in T.iterrows():
    done.add(t.gse)
    base = dict(part='P1b', source='GEO search', accession=t.gse, perturbation=f'{t.mirna} {t.direction}'.strip(), cell=t.cell, cls=t['class'], title=title.get(t.gse, ''))
    if t.decision == 'exclude': add(**base, status='excluded', reason=t.reason); continue
    c = CON[(CON.gse == t.gse) & (CON.cell == t.cell)]
    NP = {'GSE13674': 'array raw-data table (GPL5799, GSE13674_Raw_Data.txt) not parsed in this round; positive-control dataset (miR-101 -> EZH2), which has k = 6 without it',
          'GSE137048': "the planned arms (the DHSM clone samples) are absent from NCBI's counts and the submitter table has 1,880 genes; the NCBI-count contrast used in P2 (control-ASO clones vs wild type, 24 h) was not run for P1b; 'other' cell type, outside the meta-regression",
          'GSE310550': "RNA-seq without NCBI-generated counts; the submitter RSEM table was not processed in this round; 'other' cell type, outside the meta-regression",
          'GSE159688': "RNA-seq without NCBI-generated counts; the submitter files were not processed in this round; 'other' cell type, outside the meta-regression"}
    if not len(c): add(**base, status='not processed', reason=f'candidate ({t.reason}); ' + NP.get(t.gse, 'no contrast computed')); continue
    for r in c.itertuples():
        st = 'included' if r.tag == 'primary' else f'reported ({r.tag})'
        add(**dict(base, n_pert=r.n_pert, n_ctrl=r.n_ctrl, data=r.route, cls=r.cls), status=st, reason=t.reason,
            flags='processed' if r.route in ('matrix', 'matrix2') else '')
for g in H.gse:
    if g in done: continue
    add(part='P1b', source='GEO search', accession=g, title=title.get(g, ''), status='excluded', reason=REREAD.get(('P1b', g), NOMENTION_P1B))

# ------------------------------------------------------------------ P2
H = pd.read_csv(f'{INB}/P2/P2_hits.tsv', sep='\t'); title = dict(zip(H.gse, H.title))
RS = pd.read_csv(f'{INB}/P2/p2_reasons.tsv', sep='\t')
PD = pd.read_csv(f'{INB}/P2/p2_datasets.tsv', sep='\t').drop_duplicates(['gse', 'cell'])
used = {}
for mp in glob.glob(f'{INB}/P2/de/*.meta.json'):
    m = json.load(open(mp)); used.setdefault(m['gse'], []).append(m)
FP = pd.read_csv(f'{INB}/P2/p2_datasets_first_pass.tsv', sep='\t'); FPK = set(zip(FP.gse, FP.cell))
done = set()
for g, ms in used.items():
    for m in ms:
        src = ('GEO search' if g in set(H.gse) else 'found in the P1b search') + ('' if (g, m['cell']) in FPK else ' (second reading)')
        add(part='P2', source=src, accession=g, perturbation=f"{m['tf']} {m.get('direction', '')}",
            cell=m['cell'], cls=m['cls'], title=title.get(g, ''), n_pert=m.get('n_pert'), n_ctrl=m.get('n_ctrl'), data=m.get('route', ''), status='included',
            reason=f"members {m.get('mirna', '')}", flags=('processed' if m.get('route') in ('matrix', 'processed DE', 'linear_matrix') or 'processed' in m.get('note', '') else ''))
    done.add(g)
for r in RS.itertuples():
    if r.gse in done: continue
    add(part='P2', source='GEO search', accession=r.gse, title=title.get(r.gse, ''), status=r.status, reason=r.reason); done.add(r.gse)
for g in H.gse:
    if g in done: continue
    add(part='P2', source='GEO search', accession=g, title=title.get(g, ''), status='excluded',
        reason='no sample names a perturbation of a listed cluster or of two or more of its members')

# ------------------------------------------------------------------ P3a
if os.path.exists(f'{INB}/P3a/p3a_triage.tsv'):
    H = pd.read_csv(f'{INB}/P3a/P3a_hits.tsv', sep='\t'); title = dict(zip(H.gse, H.title))
    T = pd.read_csv(f'{INB}/P3a/p3a_triage.tsv', sep='\t')
    U = pd.read_csv(f'{INB}/P3/p3a_datasets_used.tsv', sep='\t')
    done = set()
    for r in U.itertuples():
        add(part='P3a', source={'P3': 'GEO search', 'P1b': 'P1b dataset', 'P2': 'P2 dataset'}.get(r.source, r.source), accession=r.gse, perturbation=f'{r.members} {r.direction}',
            cell=r.cell, n_total=r.n_samples, title=title.get(r.gse, ''), status='included' if r.used else 'qualified, not used',
            reason='among the 40 contrasts with the most samples' if r.used else 'qualified; outside the 40 contrasts with the most samples')
        done.add(r.gse)
    for r in T.itertuples():
        if r.gse in done and r.decision != 'exclude': continue
        add(part='P3a', source='GEO search', accession=r.gse, perturbation=getattr(r, 'mirna', ''), cell=getattr(r, 'cell', ''), title=title.get(r.gse, ''),
            status='excluded' if r.decision == 'exclude' else r.decision, reason=r.reason)
        done.add(r.gse)
    for g in H.gse:
        if g in done: continue
        n = H.set_index('gse').loc[g]
        add(part='P3a', source='GEO search', accession=g, title=title.get(g, ''), status='not read',
            reason=f'screen counts {n.auto_perturbed} perturbed / {n.auto_control} control samples: below the size of the 40th contrast, so it could not enter the 40 (see P3a/p3a_triage_note.txt)')

# ------------------------------------------------------------------ P3b
R = pd.read_csv(f'{INB}/P3/p3b_per_dataset.tsv', sep='\t')
DS = pd.read_csv(f'{INB}/P3/p3b_datasets.tsv', sep='\t')
st = R.drop_duplicates('sample_id').set_index('sample_id')
for d in DS.itertuples():
    base = dict(part='P3b', source='KnockTF 2.0', accession=d.sample_id, perturbation=f'{d.tf} {d.knock_method}', cell=d.biosample, title=str(d.profile_id),
                n_pert=d.n_treat, n_ctrl=d.n_control, data='KnockTF differential expression (processed)', flags='processed')
    if not (d.n_control >= 2 and d.n_treat >= 2):
        add(**base, status='excluded', reason='fewer than 2 control or 2 treated samples on the KnockTF dataset page'); continue
    if d.sample_id not in st.index: add(**base, status='excluded', reason='not evaluated'); continue
    s = st.loc[d.sample_id]
    add(**base, status='included' if s.status == 'included' else 'excluded', reason=s.status if s.status != 'included' else f'TF log2FC {s.tf_log2fc:+.2f}')

# ------------------------------------------------------------------ P5
for acc, what, why in (('FANTOM5 miRNA atlas (de Rie et al. 2017)', 'mature miRNA CPM', 'whole-cell samples; 24 subcellular-fraction samples excluded (6 in the version first delivered; SETTINGS change 16)'),
                       ('microRNAome (McCall et al. 2017; Bioconductor microRNAome 1.30.0)', 'mature miRNA counts', 'all samples of the package; the 32 iPSC_fibroblast samples (Class "Stem") are not counted as fibroblasts (SETTINGS change 16)')):
    add(part='P5', source='miRNA atlas', accession=acc, perturbation=what, status='included', reason=why)

D = pd.DataFrame(rows)
cols = ['part', 'source', 'accession', 'perturbation', 'cell', 'cls', 'status', 'reason', 'n_pert', 'n_ctrl', 'n_total', 'data', 'flags', 'result_file', 'title']
for c in cols:
    if c not in D: D[c] = ''
D = D[cols].rename(columns={'cls': 'class', 'title': 'series_title'})
D['series_title'] = D.series_title.fillna('').astype(str).str.replace('\t', ' ').str[:160]
D.to_csv(f'{INB}/DATASETS.tsv', sep='\t', index=False)
print(D.groupby(['part', 'status']).size().to_string())
