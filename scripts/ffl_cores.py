# Sub-type numbering of Mangan and Alon (2003), signs (R->M, M->T, R->T); relabelled 2026-09-26 from an earlier non-standard numbering.
"""
ffl_cores.py -- exact enumeration and Mangan & Alon (2003) coherence typing of
every 3-node FFL core (R->M, R->T, M->T) in the canonical directed signed graph.

Coherence: a core is COHERENT when sign(R->T) == sign(R->M)*sign(M->T).
Mangan & Alon numbering (PNAS 2003, Fig. 1), (R->M, M->T, R->T):
   C1 (+,+,+)   C2 (-,+,-)   C3 (+,-,-)   C4 (-,-,+)
   I1 (+,-,+)   I2 (-,-,-)   I3 (+,+,-)   I4 (-,+,+)
"""
import csv, os
DATA = "/path/to/revision/data"

TYPE_TABLE = {
    (1, 1, 1): "C1", (-1, 1, -1): "C2", (1, -1, -1): "C3", (-1, -1, 1): "C4",
    (1, -1, 1): "I1", (-1, -1, -1): "I2", (1, 1, -1): "I3", (-1, 1, 1): "I4",
}


def load_edge_annotation():
    """per-arc evidence string + curated sign (None where not curated)."""
    ann = {}
    with open(os.path.join(DATA, "layer_TF_target.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ann[(r["source"], r["target"])] = dict(
                layer="TF_target", ev="TRRUST:%s:%s_PMID" % (r["mode"], r["n_pmid"]),
                curated_sign=(int(r["sign"]) if int(r["sign"]) != 0 else None),
                mode=r["mode"])
    with open(os.path.join(DATA, "layer_TF_miRNA.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ann[(r["source"], r["target"])] = dict(
                layer="TF_miRNA", ev="TransmiR:%s:%s" % (r["mode"], r["evidence"]),
                curated_sign=(int(r["sign"]) if int(r["sign"]) != 0 else None),
                mode=r["mode"])
    tier = {}
    with open(os.path.join(DATA, "edge_evidence_tier.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tier[(r["source"], r["target"])] = r["tier"]
    with open(os.path.join(DATA, "canonical_edges.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["edge_type"] == "miRNA_target":
                k = (r["source"], r["target"])
                ann[k] = dict(layer="miRNA_target",
                              ev="multiMiR:%s" % tier.get(k, "NA"),
                              curated_sign=-1, mode="Repression")
            elif r["edge_type"] in ("TF_target", "TF_miRNA") and (r["source"], r["target"]) not in ann:
                ann[(r["source"], r["target"])] = dict(
                    layer=r["edge_type"], ev="canonical_only", curated_sign=None, mode="NA")
    with open(os.path.join(DATA, "layer_miRNA_miRNA.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["threshold_10kb"] != "TRUE":
                continue
            e = "cotranscription:%s:%sbp" % (r["chrom"], r["distance_bp"])
            for k in ((r["miRNA_1"], r["miRNA_2"]), (r["miRNA_2"], r["miRNA_1"])):
                ann[k] = dict(layer="miRNA_miRNA", ev=e, curated_sign=None, mode="Co-transcription")
    return ann


def enumerate_cores(nodes, ntype, arcs):
    """all (R,M,T) with R->M, M->T, R->T and <=1 undirected arc on the R->M->T arm"""
    from collections import defaultdict
    out = defaultdict(set)
    inn = defaultdict(set)
    for (u, v) in arcs:
        out[u].add(v); inn[v].add(u)
    cores = []
    for M in nodes:
        for R in inn[M]:
            oR = out[R]
            for T in out[M]:
                if T == R or T not in oR:
                    continue
                nu = (1 if arcs[(R, M)]["undirected"] else 0) + \
                     (1 if arcs[(M, T)]["undirected"] else 0)
                if nu > 1:
                    continue
                cores.append((R, M, T))
    return cores


def core_class(R, M, T, ntype, arcs):
    """manuscript FFL class of a 3-node core"""
    tR, tM = ntype[R], ntype[M]
    mutual = (M, R) in arcs and (R, M) in arcs
    if mutual and {tR, tM} == {"TF", "miRNA"}:
        return "Composite-FFL"
    if tR == "TF" and tM == "miRNA":
        return "TF-FFL"
    if tR == "miRNA" and tM == "TF":
        return "miRNA-FFL"
    return "%s-%s-%s" % (tR, tM, ntype[T])
