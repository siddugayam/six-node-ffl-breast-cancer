#!/usr/bin/env python3
"""TargetScan 8.0 Summary Counts (all predictions) -> human rows only (Species ID 9606), one row per gene and miRNA
family, written to the cache as ts80_human_family_sites.tsv.gz (gene, family, cons8, cons7m8, cons7a1, cons_total,
noncons_total).  Used by P1b, P2 and P3a."""
import os, sys, zipfile, io, gzip, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import _dl
src = os.path.join(_dl.CACHE, 'common/targetscan/Summary_Counts.all_predictions.txt.zip')
dst = os.path.join(_dl.CACHE, 'common/targetscan/ts80_human_family_sites.tsv.gz')
z = zipfile.ZipFile(src); name = z.namelist()[0]
n = 0
with z.open(name) as f, gzip.open(dst, 'wt') as o:
    r = csv.reader(io.TextIOWrapper(f, encoding='utf-8', errors='replace'), delimiter='\t'); h = next(r)
    ix = {k: h.index(k) for k in ('Gene Symbol', 'miRNA family', 'Species ID', 'Total num conserved sites', 'Number of conserved 8mer sites',
                                  'Number of conserved 7mer-m8 sites', 'Number of conserved 7mer-1a sites', 'Total num nonconserved sites')}
    o.write('gene\tfamily\tcons8\tcons7m8\tcons7a1\tcons_total\tnoncons_total\n')
    for row in r:
        if row[ix['Species ID']] != '9606': continue
        o.write('\t'.join(row[ix[k]] for k in ('Gene Symbol', 'miRNA family', 'Number of conserved 8mer sites', 'Number of conserved 7mer-m8 sites',
                                              'Number of conserved 7mer-1a sites', 'Total num conserved sites', 'Total num nonconserved sites')) + '\n'); n += 1
print('human gene x family rows', n, 'header', h[:6])
