#!/usr/bin/env python3
"""PWM scan of all TCGA-BRCA ATAC peaks (Corces 2018 BRCA peak set) with JASPAR 2024
CORE vertebrate PFMs.

Log-odds PWM, 0-order background estimated from the peak sequences themselves,
per-position p-value threshold obtained by exact dynamic programming over the
quantised score distribution (the FIMO convention, p < 1e-4). Both strands.
Outputs per-peak hit counts for each motif.
"""
import os, sys, time, numpy as np

ROOT = "/path/to/revision"
CACHE = os.path.join(ROOT, "cache/v6/atac")
PVAL = 1e-4
SCALE = 100

def read_jaspar(path):
    mot = {}
    mid = name = None; rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if mid: mot[mid] = (name, np.array(rows, dtype=float))
            parts = line[1:].split()
            mid, name = parts[0], parts[1] if len(parts) > 1 else parts[0]
            rows = []
        elif line.strip():
            nums = line.split("[")[1].split("]")[0] if "[" in line else line.split(None, 1)[1]
            rows.append([float(x) for x in nums.split()])
    if mid: mot[mid] = (name, np.array(rows, dtype=float))
    return {k: (v[0], v[1].T) for k, v in mot.items()}   # -> (w x 4) counts, ACGT

def logodds(counts, bg):
    N = counts.sum(axis=1, keepdims=True)
    pseudo = bg[None, :] * 1.0                 # 1 total pseudocount split by background
    p = (counts + pseudo) / (N + 1.0)
    return np.log2(p / bg[None, :])

def score_threshold(lo, bg, alpha=PVAL, scale=SCALE):
    ints = np.rint(lo * scale).astype(np.int64)
    offset = 0
    dist = np.array([1.0])
    for j in range(ints.shape[0]):
        mn = ints[j].min(); offset += mn
        sh = ints[j] - mn
        v = np.zeros(sh.max() + 1)
        for b in range(4): v[sh[b]] += bg[b]
        dist = np.convolve(dist, v)
    cum = np.cumsum(dist[::-1])[::-1]           # P(score >= idx)
    ok = np.where(cum <= alpha)[0]
    thr_idx = ok[0] if len(ok) else len(dist) - 1
    return (thr_idx + offset) / scale

def main():
    t0 = time.time()
    seqfile = os.path.join(CACHE, "brca_peak_seq.tab")
    names, seqs = [], []
    for line in open(seqfile):
        a, b = line.rstrip("\n").split("\t")
        names.append(a.split("::")[0]); seqs.append(b.upper())
    L = len(seqs[0])
    assert all(len(s) == L for s in seqs[:1000]), "non fixed-width peaks"
    print("peaks", len(seqs), "width", L, flush=True)

    code = np.full(256, 4, dtype=np.int8)
    for i, c in enumerate("ACGT"): code[ord(c)] = i
    S = np.frombuffer("".join(seqs).encode(), dtype=np.uint8)
    S = code[S].reshape(len(seqs), L)
    del seqs

    cnt = np.bincount(S.ravel(), minlength=5)[:4].astype(float)
    bg = cnt / cnt.sum()
    print("background ACGT", np.round(bg, 4), flush=True)

    mot = read_jaspar(os.path.join(CACHE, "JASPAR2024_CORE_vert_nr_pfms.txt"))
    panel = [l.strip() for l in open(os.path.join(CACHE, "motif_panel.txt")) if l.strip()]
    print("motifs to scan", len(panel), flush=True)

    out_hits = np.zeros((len(names), len(panel)), dtype=np.int16)
    meta = []
    CH = 20000
    for mi, mid in enumerate(panel):
        if mid not in mot:
            print("MISSING", mid, flush=True); meta.append((mid, "NA", 0, np.nan)); continue
        nm, counts = mot[mid]
        lo = logodds(counts, bg)
        thr = score_threshold(lo, bg)
        w = lo.shape[0]
        # reverse complement PWM
        lo_rc = lo[::-1, ::-1]
        tab = np.vstack([lo, np.full((1, 4), -1e4)])          # index 4 = N
        tabr = np.vstack([lo_rc, np.full((1, 4), -1e4)])
        npos = L - w + 1
        for st in range(0, len(names), CH):
            en = min(st + CH, len(names))
            sub = S[st:en]
            sc_f = np.zeros((en - st, npos), dtype=np.float32)
            sc_r = np.zeros((en - st, npos), dtype=np.float32)
            for j in range(w):
                col = sub[:, j:j + npos]
                sc_f += lo[j][np.minimum(col, 3)] * (col < 4) + (col == 4) * (-1e4)
                sc_r += lo_rc[j][np.minimum(col, 3)] * (col < 4) + (col == 4) * (-1e4)
            out_hits[st:en, mi] = (sc_f >= thr).sum(1) + (sc_r >= thr).sum(1)
        meta.append((mid, nm, w, thr))
        print(f"  {mid} {nm} w={w} thr={thr:.2f} hits/peak={out_hits[:, mi].mean():.3f} "
              f"frac>=1={np.mean(out_hits[:, mi] > 0):.3f} [{time.time()-t0:.0f}s]", flush=True)

    np.save(os.path.join(CACHE, "motif_hits.npy"), out_hits)
    with open(os.path.join(CACHE, "motif_meta.tsv"), "w") as fh:
        fh.write("matrix_id\tname\twidth\tthreshold_bits\n")
        for m in meta: fh.write(f"{m[0]}\t{m[1]}\t{m[2]}\t{m[3]}\n")
    with open(os.path.join(CACHE, "motif_peaknames.txt"), "w") as fh:
        fh.write("\n".join(names))
    print("DONE", time.time() - t0, flush=True)

if __name__ == "__main__":
    main()
