#!/usr/bin/env python3
"""
INBOX_2026-09-27c, Q1-Q3: the miRNA-TF pair filter files supplied by the authors
(hypergeometrictestoutput.zip, 2026-09-27; the eight files whose md5s are lines 1-8 of
INBOX_2026-09-27b/E/E5/e5a_hypergeometric_filter.txt).  Everything is read from the zip itself.
Read-only on the project; reads REV/data/db/transmir_hsa.tsv and trrust_human.tsv for provenance only.

Q1  retained rows of the right-hand Benjamini-Hochberg block: distinct miRNA entries, and how "Yes" was decided.
Q2  can the N = 256 universe (the constant in HYPGEOM.DIST) be identified from the eight files?
Q3  are the eight files the complete inputs and output, is there restricted / patient-level content,
    and is there any other first-round file of the filter in the project?
Outputs (this folder): q123_filter_files.txt (one statement per line; cited by line in VALUES.md),
q1_retained_mirna_entries.txt, q2_universe_256_genes.txt.
usage: q123_filter_files.py <author zip>
"""
import sys, os, io, re, zipfile, hashlib, collections, subprocess
AUTHORS = {'Gayam Prasanna Kumar Reddy', 'Jesil Mathew A', 'Fayaz Shaik Mahammad'}
creator_txt = lambda c: f"'{c}'" if c in AUTHORS else 'a non-author'   # other people's names are not written (authors' rule, 2026-09-27)
import openpyxl
ZIP = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
REV = '/path/to/revision'
z = zipfile.ZipFile(ZIP)
B = {n: z.read(n) for n in z.namelist()}
out = []
P = lambda *a: (out.append(' '.join(str(x) for x in a)), print(*a))
for n in sorted(B): P(f'md5 {hashlib.md5(B[n]).hexdigest()}  {n}  ({len(B[n])} bytes)')
rd = lambda n: [l.split('\t') for l in B[n].decode('utf-8', 'replace').replace('\r', '').split('\n') if l.strip()]
wbF = openpyxl.load_workbook(io.BytesIO(B['Hypergeometric_Test_2.xlsx']), data_only=False).worksheets[0]
wbV = openpyxl.load_workbook(io.BytesIO(B['Hypergeometric_Test_2.xlsx']), data_only=True).worksheets[0]
ROWS = range(3, 2679)
Lf = [[wbF.cell(r, c).value for c in range(2, 8)] for r in ROWS]          # B-G formulas
L = [[wbV.cell(r, c).value for c in range(2, 8)] for r in ROWS]           # miRNA, TF, Nmir, Ntf, common, P
Rf = [[wbF.cell(r, c).value for c in range(9, 16)] for r in ROWS]         # I-O formulas
R = [[wbV.cell(r, c).value for c in range(9, 16)] for r in ROWS]          # rank, P, miRNA, TF, crit, adj, sig
num = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
assert wbV.cell(2679, 9).value is None and wbV.cell(2, 15).value.startswith('Significant')

# ------------------------------------------------------------------ Q1
m_formula = {re.sub(r'\d+', '#', str(x[4])) for x in Rf} | {re.sub(r'\d+(?=\))|(?<=I)\d+|(?<=J)\d+', '#', str(x[5])) for x in Rf}
P(f"Q1 formulas: P-value (column G) {sorted({re.sub(r'[A-Z]+[0-9]+', 'cell', str(x[5])) for x in Lf})}; critical value (column M) "
  f"{sorted({re.sub(r'I[0-9]+', 'rank', str(x[4])) for x in Rf})}; adjusted P (column N) {sorted({re.sub(r'[IJ][0-9]+', lambda k: 'P' if k.group(0)[0] == 'J' else 'rank', str(x[5])) for x in Rf})}; "
  f"'Significant using an FDR of 0.05?' (column O) is typed text in {sum(1 for x in Rf if isinstance(x[6], str) and not x[6].startswith('='))} of {len(Rf)} rows (no formula)")
yes = [x for x in R if x[6] == 'Yes']
ry = [x[0] for x in yes]
P(f'Q1 flags: Yes {len(yes)}, No {sum(1 for x in R if x[6] == "No")}; Yes rows are ranks {min(ry)}-{max(ry)}, contiguous: {ry == list(range(1, len(yes) + 1))}')
m = 2676
pv = [x[1] if num(x[1]) else None for x in R]
rules = {'P < (rank/2676)*0.05': [p is not None and p < x[0] / m * 0.05 for p, x in zip(pv, R)],
         'P <= (rank/2676)*0.05': [p is not None and p <= x[0] / m * 0.05 for p, x in zip(pv, R)],
         'adjusted P (P*2676/rank) < 0.05': [p is not None and p * m / x[0] < 0.05 for p, x in zip(pv, R)]}
k = max(i for i, (p, x) in enumerate(zip(pv, R)) if p is not None and p <= x[0] / m * 0.05)
rules['BH step-up (all ranks up to the largest rank with P <= critical value)'] = [i <= k for i in range(len(R))]
flag = [x[6] == 'Yes' for x in R]
P('Q1 "Yes" reproduced by: ' + '; '.join(f'{nm}: {sum(a != b for a, b in zip(r, flag))} disagreements' for nm, r in rules.items()))
nn = [x for x in R if not num(x[1])]
P(f'Q1 P values ranked: {len(R)} rows (m = 2676 is typed into every critical-value and adjusted-P formula); numeric P {sum(1 for p in pv if p is not None)}; '
  f"'{nn[0][1]}' {len(nn)} (ranks {nn[0][0]}-{nn[-1][0]}, flags {dict(collections.Counter(x[6] for x in nn))}); rows repeating an earlier (miRNA, TF) pair "
  f"{len(R) - len({(x[2], str(x[3])) for x in R})}")
Lp = collections.defaultdict(list)
for x in L: Lp[(x[0], str(x[1]))].append(x[5])
eq = sum(1 for x in R if num(x[1]) and any(num(q) and abs(q - x[1]) <= 1e-12 * max(1.0, abs(q)) for q in Lp.get((x[2], str(x[3])), [])))
odd = [x for x in R if num(x[1]) and not any(num(q) and abs(q - x[1]) <= 1e-12 * max(1.0, abs(q)) for q in Lp.get((x[2], str(x[3])), []))]
P(f'Q1 right-block P equals the left-block lower-tail P of the same (miRNA, TF) pair: {eq} of {sum(1 for p in pv if p is not None)} numeric rows; '
  f'not matched by name: {[(x[2], x[3]) for x in odd]} (left block holds that TF as the date {[str(y[1]) for y in L if y[0] == odd[0][2] and not isinstance(y[1], str)]}, i.e. DEC1 converted by Excel)' if odd else '')
ym = sorted({str(x[2]) for x in yes})
P(f'Q1 retained rows: distinct miRNA entries {len(ym)} (hsa-mir-* {sum(s.startswith("hsa-mir-") for s in ym)}, hsa-let-* {sum(s.startswith("hsa-let-") for s in ym)}, '
  f'other {sum(not s.startswith(("hsa-mir-", "hsa-let-")) for s in ym)}); case-insensitive duplicates {len(ym) - len({s.lower() for s in ym})}; mature arm names (-3p/-5p) '
  f'{sum(s.endswith(("-3p", "-5p")) for s in ym)}; distinct (miRNA, TF) pairs {len({(x[2], str(x[3])) for x in yes})}; distinct TF entries {len({str(x[3]) for x in yes})}')
P(f'Q1 all tested rows: distinct miRNA entries {len({str(x[0]) for x in L})}; distinct TF entries {len({str(x[1]) for x in L})}')
open(f'{HERE}/q1_retained_mirna_entries.txt', 'w').write('\n'.join(ym) + '\n')

# ------------------------------------------------------------------ Q2
mg, tg = rd('mirna-gene.txt'), rd('tf-gene.txt')
Tm, Tt = collections.defaultdict(set), collections.defaultdict(set)
for r in mg: Tm[r[0]].add(r[1].strip())
for r in tg: Tt[r[0]].add(r[1].strip())
Gm, Gt = set().union(*Tm.values()), set().union(*Tt.values())
U = sorted(Gm & Gt)
P(f'Q2 target genes: mirna-gene.txt {len(Gm)} (of {len(Tm)} miRNAs); tf-gene.txt {len(Gt)} (of {len(Tt)} TFs); in both lists {len(U)}; in either {len(Gm | Gt)}')
open(f'{HERE}/q2_universe_256_genes.txt', 'w').write('\n'.join(U) + '\n')
P(f'Q2 genes targeted in both mirna-gene.txt and tf-gene.txt: {len(U)} = the N typed into HYPGEOM.DIST: {len(U) == 256}; written to q2_universe_256_genes.txt '
  f'(no file names the universe; the identification rests on this count)')
Nmir = {r[0]: int(r[1]) for r in rd('Nmir.txt')}; Ntf = {r[0]: int(r[1]) for r in rd('Ntf.txt')}
vl = rd('value.txt')
P(f"Q2 the counts used are not restricted to those genes: Nmir in value.txt / column D takes the value(s) {sorted({int(r[2]) for r in vl})} in all {len(vl)} rows; "
  f"Ntf.txt maximum {max(Ntf.values())}, TFs with > 256 targets {dict(sorted((k, v) for k, v in Ntf.items() if v > 256))}; all {len(nn)} '#NUM!' rows have Ntf > 256: "
  f"{all(y[3] > 256 for y in L if not num(y[5]))}; Ntf restricted to the 256 genes would be at most {max(len(Tt[t] & set(U)) for t in Tt)}")

# ------------------------------------------------------------------ Q3
Lrows = [(str(x[0]), x[1], x[2], x[3], x[4]) for x in L]
same = sum(1 for a, b in zip(vl, Lrows) if a[0] == b[0] and (a[1] == b[1] or not isinstance(b[1], str)) and int(a[2]) == b[2] and int(a[3]) == b[3] and int(a[4]) == b[4])
P(f'Q3 value.txt = the workbook input block (columns B-F, same order): {same} of {len(vl)} rows identical (TF cells compared as text except the Excel date)')
mt = rd('mirna-tf.txt'); ci = rd('common-i.txt')
cd = [(a[0], a[1], int(a[2]), int(b[4])) for a, b in zip(ci, vl) if (a[0], a[1], a[2].strip()) != (b[0], b[1], b[4].strip())]
P(f'Q3 mirna-tf.txt pairs = value.txt pairs in the same order: {[(r[0], r[1]) for r in mt] == [(r[0], r[1]) for r in vl]}; common-i.txt rows differing from value.txt '
  f'(the workbook input) in the common-target count: {len(cd)} of {len(ci)} (miRNA, TF, common-i.txt, value.txt): {cd}')
tr = [l.rstrip('\n').split('\t') for l in open(f'{REV}/data/db/transmir_hsa.tsv') if l.strip()]
tpairs = {(r[1], r[0]) for r in tr}
dmt = {(r[0], r[1]) for r in mt}
P(f'Q3 mirna-tf.txt: {len(mt)} lines, {len(dmt)} distinct pairs; in TransmiR human (REV/data/db/transmir_hsa.tsv, {len(tr)} rows, read as TF -> miRNA): '
  f'{len(dmt & tpairs)}; not found: {sorted(dmt - tpairs)}')
trr = {(l.split('\t')[0], l.split('\t')[1]) for l in open(f'{REV}/data/db/trrust_human.tsv') if l.strip()}
dtg = {(r[0], r[1].strip()) for r in tg}
P(f'Q3 tf-gene.txt: {len(dtg)} distinct pairs; in TRRUST v2 human (REV/data/db/trrust_human.tsv): {len(dtg & trr)}; differing (Excel-converted target symbols): {sorted(dtg - trr)}')
tested_m = {r[0] for r in mt}; tested_t = {r[1] for r in mt}
P(f'Q3 value.txt Nmir derivable from Nmir.txt / mirna-gene.txt: tested miRNA names {len(tested_m)} (precursor or family names such as hsa-let-7a-1), of which found in Nmir.txt '
  f'{len(tested_m & set(Nmir))} and in mirna-gene.txt {len(tested_m & set(Tm))} (these use mature names such as hsa-let-7a-5p); Nmir is {sorted({int(r[2]) for r in vl})} in every row')
dv = {(r[0], r[1]): (int(r[2]), int(r[3]), int(r[4])) for r in vl}
P(f'Q3 value.txt Ntf equal to Ntf.txt: {sum(1 for (a, t), v in dv.items() if Ntf.get(t) == v[1])} of {len(dv)} distinct pairs; TF absent from Ntf.txt: '
  f'{sum(1 for (a, t) in dv if t not in Ntf)} pairs ({len({t for (a, t) in dv if t not in Ntf})} TF entries), which still carry an Ntf value')
P(f'Q3 common-i.txt equal to |targets(miRNA) & targets(TF)| from mirna-gene.txt and tf-gene.txt by name: {sum(1 for (a, t), v in dv.items() if v[2] == len(Tm.get(a, set()) & Tt.get(t, set())))} '
  f'of {len(dv)} distinct pairs (the tested miRNA names are absent from mirna-gene.txt, so this is 0 for every pair); pairs with common > 0: {sum(1 for v in dv.values() if v[2] > 0)}')
P(f'Q3 Nmir.txt equal to the target count in mirna-gene.txt: {sum(1 for k, v in Nmir.items() if v == len(Tm.get(k, ())))} of {len(Nmir)}; differing: '
  f'{[(k, v, len(Tm.get(k, ()))) for k, v in Nmir.items() if v != len(Tm.get(k, ()))]}; Ntf.txt equal to the target count in tf-gene.txt: '
  f'{sum(1 for k, v in Ntf.items() if v == len(Tt.get(k, ())))} of {len(Ntf)}; misspelt miRNA names in mirna-gene.txt: {sorted(k for k in Tm if not k.startswith(("hsa-miR-", "hsa-let-")))}')
# patient-level / restricted content scan
txt = {n: B[n].decode('utf-8', 'replace') for n in B if n.endswith('.txt')}
xl_strings = ' '.join(str(c.value) for row in wbV.iter_rows() for c in row if c.value is not None)
allt = ' '.join(txt.values()) + ' ' + xl_strings
pats = {'TCGA barcode': r'TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}', 'GEO sample (GSM)': r'\bGSM\d{3,}', 'e-mail address': r'[\w.+-]+@[\w-]+\.[\w.]+',
        'date of birth / age field': r'(?i)\b(dob|date of birth|age_at|patient|sample_id|barcode)\b'}
P('Q3 patient-level / identifier scan of all eight files: ' + '; '.join(f'{k} {len(re.findall(v, allt))}' for k, v in pats.items()) +
  f"; the tokens matched: {dict(collections.Counter(re.findall(r'(?i)[A-Za-z0-9-]*(?:dob|date of birth|age_at|patient|sample_id|barcode)[A-Za-z0-9-]*', allt)))}")
kinds = {}
for n, t in txt.items():
    rows = [l.split('\t') for l in t.replace('\r', '').split('\n') if l.strip()]
    ncol = collections.Counter(len(r) for r in rows)
    kinds[n] = f"{len(rows)} rows, columns {dict(ncol)}"
P('Q3 text files (rows, columns): ' + '; '.join(f'{n} {v}' for n, v in sorted(kinds.items())))
core = B['Hypergeometric_Test_2.xlsx']
zx = zipfile.ZipFile(io.BytesIO(core))
cp = zx.read('docProps/core.xml').decode()
g = lambda tag: (re.search(f'<{tag}[^>]*>(.*?)</{tag}>', cp) or [None, ''])[1]
parts = zx.namelist()
P(f"Q3 workbook: sheets {openpyxl.load_workbook(io.BytesIO(core)).sheetnames} (states {[ws.sheet_state for ws in openpyxl.load_workbook(io.BytesIO(core)).worksheets]}); "
  f"document properties: creator {creator_txt(g('dc:creator'))}, lastModifiedBy '{g('cp:lastModifiedBy')}', created {g('dcterms:created')}, modified {g('dcterms:modified')}; "
  f"comment parts {[p for p in parts if 'comment' in p.lower()]}; persons part empty: {'<person ' not in zx.read('xl/persons/person.xml').decode()}; "
  f"external links {[p for p in parts if 'externalLink' in p]}; note cells outside the two blocks: "
  f"{[(c.coordinate, c.value) for row in wbV.iter_rows() for c in row if c.value is not None and (c.row > 2678 or c.column > 15)]}")
P(f"Q3 Excel conversions: workbook TF cells that are not text {[(str(x[0]), str(x[1])) for x in L if not isinstance(x[1], str)]} "
  f"(right block {[(x[2], x[3]) for x in R if not isinstance(x[3], str)]}); mirna-tf.txt / value.txt TF entries that look like dates "
  f"{sorted({r[1] for r in mt if re.fullmatch(r'\d{1,2}-[A-Z][a-z]{2}|[A-Z][a-z]{2}-\d{2}|\d{5}', r[1])})}")
# other first-round files of the filter anywhere in the project
pat = ['*hypergeom*', 'Hypergeometric*', 'Nmir*', 'Ntf*', 'common-i*', 'mirna-tf*', 'tf-gene*', 'mirna-gene*', 'value.txt', '*interaction_pairs*', '*filtered_interactions*', '*pair_filter*']
cmd = ['find', '/path/to/home/Desktop/DD/R_GPR/miRNA_FFL', '-path', f'{REV}/INBOX_2026-09-27c', '-prune', '-o', '(']
for i, s in enumerate(pat): cmd += (['-o'] if i else []) + ['-iname', s]
cmd += [')', '-print']
hits = sorted(subprocess.run(cmd, capture_output=True, text=True).stdout.split())
P(f"Q3 project search for other filter files (find -iname {', '.join(pat)} under /path/to/home/Desktop/DD/R_GPR/miRNA_FFL): "
  f"{[h.replace('/path/to/home/Desktop/DD/R_GPR/miRNA_FFL/', '') for h in hits]}")
open(f'{HERE}/q123_filter_files.txt', 'w').write('\n'.join(out) + '\n')
