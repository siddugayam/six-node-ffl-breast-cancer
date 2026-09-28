#!/usr/bin/env python3
"""Export the exact graph 03_ffl_census.py builds, for the C enumerator."""
import sys, importlib.util, collections
REV='/path/to/revision'
spec=importlib.util.spec_from_file_location('cen', f'{REV}/scripts/03_ffl_census.py')
cen=importlib.util.module_from_spec(spec); sys.argv=['x','3']
try: spec.loader.exec_module(cen)
except SystemExit: pass
CLS=['TF->Gene','TF->TF','TF->miRNA','miRNA->Gene','miRNA->TF','Gene->Gene','miRNA-miRNA']
for tag,extra in (('aug',True),('dep',False)):
    nodes,edges,nrecip=cen.build_graph(extra)
    def cls(a,b):
        t=edges[(a,b)]
        if t=='TF_target': return 'TF->TF' if nodes[b]=='TF' else 'TF->Gene'
        if t=='miRNA_target': return 'miRNA->TF' if nodes[b]=='TF' else 'miRNA->Gene'
        return {'TF_miRNA':'TF->miRNA','gene_gene':'Gene->Gene','miRNA_miRNA':'miRNA-miRNA'}[t]
    V=sorted(nodes); idx={v:i for i,v in enumerate(V)}
    with open(f'{REV}/results/v2/R4_graph_{tag}.txt','w') as f:
        f.write(f"{len(V)} {len(edges)}\n")
        for v in V: f.write(f"{v}\t{nodes[v]}\n")
        for (a,b) in sorted(edges): f.write(f"{idx[a]} {idx[b]} {CLS.index(cls(a,b))}\n")
    print(tag,"nodes",len(V),"edges",len(edges),"recip_contracted",nrecip,
          "classes present",sorted(set(cls(a,b) for (a,b) in edges)))
