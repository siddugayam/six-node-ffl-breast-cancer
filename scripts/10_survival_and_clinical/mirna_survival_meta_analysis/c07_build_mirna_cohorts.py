#!/usr/bin/env python3
"""Parse GEO series matrices for breast-cancer miRNA cohorts, map probes to the
12 target miRNAs, extract clinical/outcome, and write tidy per-cohort CSVs."""
import gzip, os, re, sys, csv, json, math
import numpy as np

BASE="/path/to/revision"
MX=os.path.join(BASE,"cache/v6/cohorts/matrix")
GMAP=os.path.join(BASE,"cache/v6/cohorts/gpl/maps")
OUTD=os.path.join(BASE,"cache/v6/cohorts/built"); os.makedirs(OUTD,exist_ok=True)

TARGETS = {
 "miR-21":   r"^hsa-mir-21(-[0-9])?(-[35]p)?\*?$",
 "miR-195":  r"^hsa-mir-195(-[35]p)?\*?$",
 "miR-204":  r"^hsa-mir-204(-[35]p)?\*?$",
 "miR-383":  r"^hsa-mir-383(-[35]p)?\*?$",
 "miR-124":  r"^hsa-mir-124a?(-[123])?(-[35]p)?\*?$",
 "miR-155":  r"^hsa-mir-155(-[35]p)?\*?$",
 "miR-429":  r"^hsa-mir-429(-[35]p)?\*?$",
 "miR-141":  r"^hsa-mir-141(-[35]p)?\*?$",
 "miR-34a":  r"^hsa-mir-34a(-[35]p)?\*?$",
 "miR-101":  r"^hsa-mir-101(-[12])?(-[35]p)?\*?$",
 "miR-130a": r"^hsa-mir-130a(-[35]p)?\*?$",
 "miR-29a":  r"^hsa-mir-29a(-[35]p)?\*?$",
}
TRE={k:re.compile(v,re.I) for k,v in TARGETS.items()}

def norm(nm):
    nm=(nm or "").strip()
    if not nm: return ""
    nm=nm.replace("miR","mir").replace("MIR","mir").lower()
    if not nm.startswith("hsa-") and nm.startswith("mir"): nm="hsa-"+nm
    return nm

def read_matrix(path):
    """returns (probes, samples, X ndarray, chars dict, meta dict)"""
    chars=[]; sample_ids=None; meta={}
    rows=[]; probes=[]
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as fh:
        intab=False
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"): intab=True; continue
            if line.startswith("!series_matrix_table_end"): break
            if intab:
                p=[x.strip('"') for x in line.split("\t")]
                if sample_ids is None: sample_ids=p[1:]; continue
                probes.append(p[0])
                rows.append([np.nan if v in ("","NA","null","NaN") else float(v) if _isnum(v) else np.nan for v in p[1:]])
            else:
                if line.startswith("!Sample_characteristics_ch1"):
                    chars.append([x.strip('"') for x in line.split("\t")[1:]])
                elif line.startswith("!Sample_geo_accession"):
                    meta["gsm"]=[x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_title"):
                    meta["title"]=[x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_platform_id"):
                    meta["gpl"]=[x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_source_name_ch1"):
                    meta["source"]=[x.strip('"') for x in line.split("\t")[1:]]
    X=np.array(rows,dtype=float)
    return probes, meta.get("gsm",sample_ids), X, chars, meta

def _isnum(v):
    try: float(v); return True
    except Exception: return False

def chars_to_dict(chars, n):
    """list of rows -> dict key -> list(n)"""
    out={}
    for row in chars:
        for j,cell in enumerate(row[:n]):
            if ":" in cell:
                k,v=cell.split(":",1); k=k.strip().lower(); v=v.strip()
            else:
                k="_free"; v=cell.strip()
            out.setdefault(k,[None]*n)
            if out[k][j] is None: out[k][j]=v
    return out

def load_map(gpl):
    p=os.path.join(GMAP,gpl+"_map.tsv")
    if not os.path.exists(p): return None
    m={}
    with open(p) as f:
        r=csv.DictReader(f,delimiter="\t")
        cols=[c for c in r.fieldnames if c!="probe"]
        for row in r:
            names=[]
            for c in cols:
                v=row.get(c) or ""
                for piece in re.split(r"[,;|/]", v):
                    piece=piece.strip()
                    if piece: names.append(piece)
            m[row["probe"]]=names
    return m

def match_targets(probes, gplmap):
    """probe list -> {target: [probe indices]} with the probe names used"""
    hits={t:[] for t in TARGETS}
    for i,p in enumerate(probes):
        cands=[p]
        if gplmap and p in gplmap: cands=cands+gplmap[p]
        for c in cands:
            nc=norm(c)
            for t,rx in TRE.items():
                if rx.match(nc):
                    hits[t].append((i,c)); break
            else:
                continue
            break
    return hits

def build(cohort, matrix_file, gpl, extra_note=""):
    path=os.path.join(MX,matrix_file)
    probes,gsm,X,chars,meta=read_matrix(path)
    n=len(gsm)
    cd=chars_to_dict(chars,n)
    gplmap=load_map(gpl)
    hits=match_targets(probes,gplmap)
    # dominant probe per target = highest mean
    expr={}; probe_used={}; probe_all={}
    for t,lst in hits.items():
        if not lst: continue
        # prefer the guide arm: names without a '*' (passenger/star) suffix
        guide=[(i,nm) for i,nm in lst if not nm.strip().endswith("*")]
        pool = guide if guide else lst
        best=None; bestmu=-1e18
        for i,nm in pool:
            col=X[i,:]
            if not np.any(np.isfinite(col)): continue
            mu=np.nanmean(col)
            if np.isfinite(mu) and mu>bestmu: bestmu=mu; best=(i,nm)
        if best is None: continue
        expr[t]=X[best[0],:]
        probe_used[t]="%s|%s"%(probes[best[0]],best[1])
        probe_all[t]=";".join(sorted(set("%s(%s)"%(probes[i],nm) for i,nm in lst)))
    return dict(cohort=cohort, gsm=gsm, n=n, probes=probes, X=X, chars=cd,
                expr=expr, probe_used=probe_used, probe_all=probe_all,
                gpl=gpl, note=extra_note, source=meta.get("source",[""]*n),
                title=meta.get("title",[""]*n))

def write_cohort(b, clin):
    """clin: dict of column -> list"""
    p=os.path.join(OUTD, b["cohort"]+".csv")
    cols=["gsm"]+list(clin.keys())+sorted(b["expr"].keys())
    with open(p,"w",newline="") as f:
        w=csv.writer(f); w.writerow(cols)
        for j in range(b["n"]):
            row=[b["gsm"][j]]+[clin[k][j] for k in clin]+[
                 ("" if not np.isfinite(b["expr"][t][j]) else "%.6f"%b["expr"][t][j]) for t in sorted(b["expr"])]
            w.writerow(row)
    print("  wrote", p, "targets:", sorted(b["expr"].keys()))
    return p
