#!/usr/bin/env python3
"""41_seqreg_ext_occlusion_null.py

EXTENSION of 22b_seqreg_enformer_ism.py.

The deposited occlusion table (seqreg_enformer_motif_occlusion.csv) occluded ONLY the
NF-kB / SP / ETS focus sites, so a statement like "shuffling the NFKB1 site changes
predicted fibroblast CAGE by -0.75%" has no scale attached: we do not know what a
load-bearing site at that promoter looks like, nor what an irrelevant one looks like.

This script calibrates those effects, in two experiments.

EXPERIMENT 1 -- ALL-SITE PROMOTER OCCLUSION (a null distribution for the focus sites).
  Every JASPAR2024/HOCOMOCOv11 site called at p<1e-4 in the COL1A1 and COL3A1 promoters
  (618 and 1,008 sites) is occluded in place by N_SHUF dinucleotide-preserving shuffles.
  Each focus-TF site can then be reported as a PERCENTILE of the promoter's own
  occlusion distribution, and the most load-bearing sites (the de-facto positive
  controls) are identified rather than assumed.

EXPERIMENT 2 -- DISTAL FIBROBLAST RELA PEAK OCCLUSION.
  The 20 ReMap 2022 NF-kB-family peaks observed in fibroblast/mesenchymal biotypes
  within +/-100 kb of the collagen TSSs (seqreg_rela_fibroblast_peaks.csv) all fall
  inside Enformer's 196,608 bp receptive field (max |dist| = 90,025). Each whole peak is
  dinucleotide-shuffled and the effect on predicted fibroblast CAGE AT THE COLLAGEN TSS
  is measured, against a null of width-matched random control windows drawn from the
  same locus. This asks whether the model treats the RELA peaks that are actually
  observed in fibroblasts as functional for collagen transcription.

Conventions (window construction, rel_to_index, dinucleotide shuffle, scoring) are copied
verbatim from 22b so the numbers are directly comparable; the reference prediction is
asserted against the deposited value before anything else runs.
"""
import os, sys, time, random
import numpy as np, pandas as pd, torch
from enformer_pytorch import from_pretrained

ROOT = "/path/to/revision"
RES, CACHE = f"{ROOT}/results/v3", f"{ROOT}/cache/seqreg"
SEQ_LEN, N_BINS, BIN = 196_608, 896, 128
CROP = (SEQ_LEN - N_BINS * BIN) // 2
CENTRE_BIN = (SEQ_LEN // 2 - CROP) // BIN
HALF = 300_000
N_SHUF = 5            # shuffles per promoter motif site
N_SHUF_PEAK = 5       # shuffles per distal peak
N_CTRL = 20           # width-matched control windows per distal peak
N_SHUF_CTRL = 2       # shuffles per control window

targets = pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets["assay"] = targets.description.str.split(":").str[0]
FIB_CAGE = targets.index[(targets.assay == "CAGE") & targets.description.str.contains("ibroblast")].to_numpy()
FIB_DNASE = targets.index[(targets.assay == "DNASE") & targets.description.str.contains("ibroblast")].to_numpy()
ALL_CAGE = targets.index[targets.assay == "CAGE"].to_numpy()
print(f"fibroblast CAGE tracks: {len(FIB_CAGE)}  fibroblast DNase tracks: {len(FIB_DNASE)}", flush=True)

MAP = {"A": 0, "C": 1, "G": 2, "T": 3}
def onehot(s):
    a = np.zeros((len(s), 4), dtype=np.float32)
    idx = np.frombuffer(s.encode(), dtype=np.uint8)
    for b, i in MAP.items():
        a[idx == ord(b), i] = 1.0
    return a
def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

model = from_pretrained("EleutherAI/enformer-official-rough").eval().cuda()
BATCH = 2
@torch.no_grad()
def predict_batch(seqs):
    outs = []
    for k in range(0, len(seqs), BATCH):
        x = torch.from_numpy(np.stack([onehot(s) for s in seqs[k:k + BATCH]])).cuda()
        with torch.autocast("cuda", dtype=torch.bfloat16):
            o = model(x)["human"]
        outs.append(o.float().cpu().numpy())
        del x, o
    torch.cuda.empty_cache()
    return np.concatenate(outs, axis=0)
def score(pred):
    prom = pred[:, CENTRE_BIN - 1:CENTRE_BIN + 2, :].mean(axis=1)
    return dict(fib_cage=prom[:, FIB_CAGE].mean(axis=1),
                fib_dnase=prom[:, FIB_DNASE].mean(axis=1),
                all_cage=prom[:, ALL_CAGE].mean(axis=1))

reg = pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
def genomic_window(regname):
    r = reg.loc[regname]; a = int(r.anchor)
    s = open(f"{CACHE}/seq/{r.chrom}_{a - HALF}_{a + HALF}.txt").read().strip()
    return s[HALF - SEQ_LEN // 2: HALF + SEQ_LEN // 2], r.strand, r.chrom, a
def rel_to_index(rel, strand):
    return SEQ_LEN // 2 + (rel if strand == "+" else -rel)

rng = random.Random(20260910)
def dinuc_shuffle(s):
    if len(s) < 3:
        return s
    for _ in range(50):
        edges = {}
        for a, b in zip(s, s[1:]):
            edges.setdefault(a, []).append(b)
        for k in edges:
            rng.shuffle(edges[k])
        out = [s[0]]; cur = s[0]; ok = True
        for _ in range(len(s) - 1):
            if not edges.get(cur):
                ok = False; break
            nxt = edges[cur].pop(); out.append(nxt); cur = nxt
        if ok and len(out) == len(s):
            return "".join(out)
    return "".join(rng.sample(s, len(s)))

# deposited reference values, asserted so the pipeline is provably the same one
DEPOSITED_REF = {"COL1A1": 259.80950927734375, "COL3A1": 114.47}

# --------------------------------------------------------------------------------------
sites = pd.read_csv(f"{RES}/seqreg_motif_sites.csv")
sites["tf_name"] = sites.motif.str.split("|").str[-1].str.split("_HUMAN").str[0]
peaks = pd.read_csv(f"{RES}/seqreg_rela_fibroblast_peaks.csv")

occ_rows, peak_rows, ctrl_rows = [], [], []

for regname in ["COL1A1", "COL3A1"]:
    ref, strand, chrom, anchor = genomic_window(regname)
    b = score(predict_batch([ref]))
    ref_fc = float(b["fib_cage"][0])
    print(f"\n=== {regname} reference fib_cage={ref_fc:.5f} "
          f"(deposited {DEPOSITED_REF[regname]}) fib_dnase={b['fib_dnase'][0]:.4f}", flush=True)
    assert abs(ref_fc - DEPOSITED_REF[regname]) < 0.02 * abs(DEPOSITED_REF[regname]), \
        f"reference prediction does not reproduce deposited value for {regname}"

    # ---------- EXPERIMENT 1: occlude every called site ----------
    ss = sites[sites.region == regname].reset_index(drop=True)
    t0 = time.time()
    for n, (_, r) in enumerate(ss.iterrows()):
        rel = int(r.rel_to_TSS); L = len(str(r.seq))
        i = rel_to_index(rel, strand)
        lo, hi = (i, i + L) if strand == "+" else (i - L + 1, i + 1)
        wt = ref[lo:hi]
        shuffles = []
        for _ in range(N_SHUF):
            sh = dinuc_shuffle(wt if strand == "+" else revcomp(wt))
            if strand == "-":
                sh = revcomp(sh)
            shuffles.append(ref[:lo] + sh + ref[hi:])
        S = score(predict_batch(shuffles))
        occ_rows.append(dict(
            region=regname, db=r.db, motif=r.motif, tf_name=r.tf_name, rel_to_TSS=rel,
            motif_strand=r.strand, motif_score=r.score, motif_p=r.pvalue, site_len=L,
            site_seq=r.seq, n_shuffles=N_SHUF, ref_fib_cage=ref_fc,
            mean_shuffled_fib_cage=float(S["fib_cage"].mean()),
            d_fib_cage=float(S["fib_cage"].mean() - ref_fc),
            pct_change_fib_cage=float(100 * (S["fib_cage"].mean() - ref_fc) / ref_fc),
            sd_fib_cage=float(S["fib_cage"].std()),
            d_fib_dnase=float(S["fib_dnase"].mean() - b["fib_dnase"][0]),
            d_all_cage=float(S["all_cage"].mean() - b["all_cage"][0])))
        if n % 50 == 0:
            print(f"  [{regname}] site {n}/{len(ss)}  {time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(occ_rows).to_csv(f"{RES}/seqreg_ext_occlusion_allsites.csv", index=False)
    print(f"  {regname}: occluded {len(ss)} sites in {time.time()-t0:.0f}s", flush=True)

    # ---------- EXPERIMENT 2: distal fibroblast RELA peaks ----------
    pk = peaks[peaks.gene == regname].reset_index(drop=True)
    forbidden = []  # peak intervals, in window-index space, to avoid when drawing controls
    for _, r in pk.iterrows():
        lo = int(r.start) - (anchor - SEQ_LEN // 2)
        forbidden.append((lo, lo + int(r.width)))
    prom_i = SEQ_LEN // 2
    for _, r in pk.iterrows():
        w = int(r.width)
        lo = int(r.start) - (anchor - SEQ_LEN // 2)
        hi = lo + w
        if lo < 0 or hi > SEQ_LEN:
            print(f"  peak {r.peak_index} outside window -- skipped", flush=True)
            continue
        wt = ref[lo:hi]
        shuffles = [ref[:lo] + dinuc_shuffle(wt) + ref[hi:] for _ in range(N_SHUF_PEAK)]
        S = score(predict_batch(shuffles))
        d = float(S["fib_cage"].mean() - ref_fc)
        # width-matched controls from the same locus
        cds = []
        tries = 0
        while len(cds) < N_CTRL and tries < 2000:
            tries += 1
            c = rng.randrange(prom_i - 98_000, prom_i + 98_000 - w)
            if abs(c - prom_i) < 1500:
                continue
            if any(not (c + w <= a or c >= b2) for a, b2 in forbidden):
                continue
            cwt = ref[c:c + w]
            if cwt.count("N") > 0:
                continue
            csh = [ref[:c] + dinuc_shuffle(cwt) + ref[c + w:] for _ in range(N_SHUF_CTRL)]
            CS = score(predict_batch(csh))
            cd = float(CS["fib_cage"].mean() - ref_fc)
            cds.append(cd)
            ctrl_rows.append(dict(region=regname, peak_index=int(r.peak_index), width=w,
                                  ctrl_offset_from_TSS=int(c - prom_i), d_fib_cage=cd))
        cds = np.array(cds)
        peak_rows.append(dict(
            region=regname, peak_index=int(r.peak_index), tf=r.tf, biotype=r.biotype,
            dataset=r.dataset, chrom=r.chrom, start=int(r.start), end=int(r.end), width=w,
            dist_to_TSS=int(r.dist_to_TSS), n_kB_motifs=r.n_kB_motifs, ccre=r.ccre,
            ref_fib_cage=ref_fc, d_fib_cage=d,
            pct_change_fib_cage=float(100 * d / ref_fc),
            sd_fib_cage=float(S["fib_cage"].std()),
            d_fib_dnase=float(S["fib_dnase"].mean() - b["fib_dnase"][0]),
            n_ctrl=len(cds), ctrl_mean_d=float(cds.mean()), ctrl_sd_d=float(cds.std()),
            ctrl_min_d=float(cds.min()), ctrl_max_d=float(cds.max()),
            z_vs_ctrl=float((d - cds.mean()) / cds.std()) if cds.std() > 0 else np.nan,
            emp_p_more_negative=float((np.sum(cds <= d) + 1) / (len(cds) + 1))))
        print(f"  peak {int(r.peak_index)} d={d:+.4f} ctrl {cds.mean():+.4f}+/-{cds.std():.4f} "
              f"z={(d-cds.mean())/cds.std() if cds.std()>0 else float('nan'):+.2f}", flush=True)
    pd.DataFrame(peak_rows).to_csv(f"{RES}/seqreg_ext_rela_peak_occlusion.csv", index=False)
    pd.DataFrame(ctrl_rows).to_csv(f"{RES}/seqreg_ext_rela_peak_controls.csv", index=False)

print("\nDONE", flush=True)
