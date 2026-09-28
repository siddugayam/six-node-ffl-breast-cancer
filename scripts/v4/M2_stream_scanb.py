#!/usr/bin/env python3
"""Stream the 592 MB SCAN-B (GSE96058) expression CSV and keep the rows for the
1,901 symbols the multi-cohort module analysis needs.
Writes cache/v4/multicohort/GSE96058_subset_v4.csv"""
import gzip, os, sys
B="/path/to/revision"
CA=os.path.join(B,"cache/v4/multicohort")
want=set(l.strip() for l in open(os.path.join(CA,"scanb_wanted_symbols.txt")) if l.strip())
print("wanted symbols:", len(want), flush=True)
src=os.path.join(B,"cache/newcohorts/GSE96058_expr.csv.gz")
out=os.path.join(CA,"GSE96058_subset_v4.csv")
found=set(); n=0
with gzip.open(src,"rt") as fi, open(out,"w") as fo:
    hdr=fi.readline(); fo.write(hdr)
    print("columns in header:", hdr.count(",")+1, flush=True)
    for line in fi:
        n+=1
        i=line.find(",")
        g=line[:i].strip().strip('"')
        if g in want:
            found.add(g); fo.write(line)
print("rows scanned:", n, " symbols kept:", len(found), flush=True)
miss=sorted(want-found)
print("symbols absent from SCAN-B:", len(miss), flush=True)
open(os.path.join(CA,"scanb_missing_symbols.txt"),"w").write("\n".join(miss))
