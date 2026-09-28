#!/usr/bin/env python3
"""12_summarise_findings.py - pull the response-letter numbers out of the result CSVs."""
import csv, collections, sys, os

REV = '/path/to/revision/results'
CLS = ['3-miR', '3-TF', '3-Comp', '4-node', '5-node', '6-node',
       'exemplar_4node', 'exemplar_5node', 'exemplar_6node']

csv.field_size_limit(10 ** 9)


def load(fn):
    p = os.path.join(REV, fn)
    if not os.path.exists(p):
        print(f'!! MISSING {fn}')
        return []
    return list(csv.DictReader(open(p)))


print('=' * 100)
print('A. LOW-SUPPORT SUMMARY')
print('=' * 100)
for r in load('enrichment_low_support_summary.csv'):
    print(f"  {r['motif_set']:<16}{r['universe']:<26}{r['ontology']:<10}"
          f"sig={r['n_sig_terms']:>5}  <3genes={r['n_sig_lt3_genes']:>5} "
          f"({r['pct_sig_lt3_genes']:>5}%)  1gene={r['n_sig_1gene']:>5}  "
          f"medianCount={r['median_Count']:>6}  min_padj={float(r['min_padj']):.2e}")

print()
print('=' * 100)
print('B. MANUSCRIPT-CLAIMED TERMS, targeted hypergeometric test (no size filter)')
print('=' * 100)
rows = load('enrichment_claimed_terms_hypergeometric.csv')
by = collections.defaultdict(dict)
for r in rows:
    by[(r['Description'], r['universe'])][r['motif_set']] = r
show = ['fibrillar collagen trimer', 'banded collagen fibril', 'Proteoglycans in cancer',
        'Antifolate resistance', 'Relaxin signaling pathway',
        'AGE-RAGE signaling pathway in diabetic complications', 'KEAP1-NFE2L2 pathway',
        'supramolecular fiber organization', 'collagen fibril organization',
        'collagen-containing extracellular matrix', 'extracellular structure organization',
        'Diabetic cardiomyopathy', 'Fluid shear stress and atherosclerosis',
        'MicroRNAs in cancer', 'Cellular senescence', 'Cell cycle',
        'p53 signaling pathway', 'Breast cancer',
        'Transcriptional misregulation in cancer', 'PI3K-Akt signaling pathway',
        'Nuclear events mediated by NFE2L2']
for uni in ['TCGA_tested_plus_network', 'network_protein_coding']:
    print(f'\n--- universe = {uni} ---')
    hdr = f"{'term':<50}" + ''.join(f'{c:<17}' for c in CLS)
    print(hdr)
    for t in show:
        d = by.get((t, uni))
        if not d:
            print(f'{t[:48]:<50}  ** term not in test table **')
            continue
        line = f'{t[:48]:<50}'
        for c in CLS:
            r = d.get(c)
            if r is None:
                line += f"{'--':<17}"
            else:
                pa = r['p.adjust']
                s = f"{r['Count']}g " + (f'{float(pa):.1e}' if pa else 'NA')
                if r['significant_BH_0.05'] != 'TRUE':
                    s += ' ns'
                line += f'{s:<17}'
        print(line)
    any_r = next((v for k, v in by.items() if k[1] == uni), None)

print('\nterm sizes / testability under clusterProfiler default filters:')
seen = set()
for r in rows:
    if r['Description'] in show and r['Description'] not in seen:
        seen.add(r['Description'])
        print(f"  {r['Description'][:52]:<54}{r['ontology']:<10}"
              f"n_total={r['term_size_total']:>5}  n_in_universeA/B varies  "
              f"tested_under_default_filters={r['tested_under_default_filters']}")

print()
print('=' * 100)
print('C. GENE SYMBOLS BEHIND THE EXEMPLAR TERMS (universe A)')
print('=' * 100)
for r in rows:
    if (r['universe'] == 'TCGA_tested_plus_network'
            and r['motif_set'].startswith('exemplar') and int(r['Count']) > 0):
        print(f"  {r['motif_set']:<16}{r['Description'][:46]:<48}"
              f"Count={r['Count']:>2}  p={float(r['pvalue']):.2e}  "
              f"padj={float(r['p.adjust']):.2e}  genes={r['geneID']}")
