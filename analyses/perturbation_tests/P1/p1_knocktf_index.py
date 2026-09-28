#!/usr/bin/env python3
"""P1a/P3b: parse the KnockTF 2.0 download table (one row per dataset) into knocktf_datasets.tsv."""
import os, re, sys, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _dl
p = _dl.fetch('http://www.licpathway.net/KnockTFv2/download.php', 'P1/knocktf', 'knocktf2_download_page.html', 'P1a', 'KnockTF 2.0', 'download page', 'not stated on the download page')
t = open(p, encoding='utf-8', errors='replace').read()
out = []
for r in re.findall(r'<tr[^>]*>(.*?)</tr>', t, flags=re.S):
    cells = [re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', ' ', c)).strip() for c in re.findall(r'<td[^>]*>(.*?)</td>', r, flags=re.S)]
    links = re.findall(r'href="(search/gene_profiles/[^"]+)"', r)
    if len(cells) >= 12 and links: out.append(cells[:12] + [links[0]])
H = ['knock_method', 'tf', 'species', 'tf_class', 'tf_superclass', 'biosample', 'tissue', 'biosample_type', 'data_source', 'profile_id', 'platform', 'pubmed', 'profile_file']
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'knocktf_datasets.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(H); w.writerows(out)
import collections
print(len(out), 'datasets;', collections.Counter(r[2] for r in out))
for tf in ('ETS1', 'NFKB1', 'RELA', 'SP1'):
    for r in out:
        if r[1].upper() == tf: print(tf, '|', ' | '.join(r[i] for i in (0, 5, 6, 7, 9, 10, 12)))
