#!/usr/bin/env python3
"""45_seqreg_ext_splicedonor_control.py

A CONTROL THE FIRST PASS DID NOT RUN, and it overturns one of its headline numbers.

The deposited occlusion table reports, as the largest focus-TF effect anywhere:
    "COL3A1: shuffling the ETS1 site changes predicted fibroblast CAGE by
     largest drop -54.482 (ref 114.47, = -47.60%) at TSS+186"
which reads as strong sequence-level evidence for ETS1 at COL3A1.

But COL3A1 exon 1 is chr2:188,974,372-188,974,568, so the exon1/intron1 junction sits at
TSS+195/+196 and the local sequence is

        ... C A A C A G G A A G | G T G A G T A G ...
                        ^^^^                          the "ETS" GGAA core (exonic)
                               ^^^^^^^^               the canonical 5' splice donor

i.e. the called ETS1 site STRADDLES the splice donor: the GGAA that makes it look like an
ETS site is the last exonic bases, and the rest of the match is the donor consensus
AG|GTGAGT. Dinucleotide-shuffling that 13 bp destroys the donor as well as the GGAA.
COL1A1 has the same geometry (exon 1 = TSS+0..+221, and the -11% sites sit at TSS+214).

This script separates the two explanations with single-purpose point mutations instead of
shuffles:
    donor_kill : GT -> AA at the first two intronic bases (donor destroyed, GGAA intact)
    ets_kill   : GGAA -> CCTT in the exon    (ETS core destroyed, donor intact)
    exon_ctrl  : 4 exonic bases mutated well away from both (neutral control)
If the effect is the splice donor, donor_kill reproduces most of the drop and ets_kill
does not.
"""
import numpy as np, pandas as pd, torch
from enformer_pytorch import from_pretrained

ROOT = "/path/to/revision"
RES, CACHE = f"{ROOT}/results/v3", f"{ROOT}/cache/seqreg"
SEQ_LEN, N_BINS, BIN = 196_608, 896, 128
CROP = (SEQ_LEN - N_BINS * BIN) // 2
CENTRE_BIN = (SEQ_LEN // 2 - CROP) // BIN
HALF = 300_000

targets = pd.read_csv(f"{CACHE}/enformer_targets_human.txt", sep="\t")
targets["assay"] = targets.description.str.split(":").str[0]
FIB_CAGE = targets.index[(targets.assay == "CAGE") & targets.description.str.contains("ibroblast")].to_numpy()
ALL_CAGE = targets.index[targets.assay == "CAGE"].to_numpy()

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
@torch.no_grad()
def predict(seqs):
    outs = []
    for k in range(0, len(seqs), 2):
        x = torch.from_numpy(np.stack([onehot(s) for s in seqs[k:k + 2]])).cuda()
        with torch.autocast("cuda", dtype=torch.bfloat16):
            o = model(x)["human"]
        outs.append(o.float().cpu().numpy()); del x, o
    torch.cuda.empty_cache()
    p = np.concatenate(outs, 0)[:, CENTRE_BIN - 1:CENTRE_BIN + 2, :].mean(1)
    return p[:, FIB_CAGE].mean(1), p[:, ALL_CAGE].mean(1)

reg = pd.read_csv(f"{RES}/seqreg_regions.csv").set_index("region")
def window(rn):
    r = reg.loc[rn]; a = int(r.anchor)
    s = open(f"{CACHE}/seq/{r.chrom}_{a - HALF}_{a + HALF}.txt").read().strip().upper()
    return s[HALF - SEQ_LEN // 2: HALF + SEQ_LEN // 2], r.strand, a
def idx(rel, strand):
    return SEQ_LEN // 2 + (rel if strand == "+" else -rel)

def put(ref, rel, strand, newbases):
    """write newbases (given in transcription orientation) starting at TSS+rel"""
    s = list(ref)
    for k, b in enumerate(newbases):
        i = idx(rel + k, strand)
        s[i] = b if strand == "+" else revcomp(b)
    return "".join(s)

def read(ref, rel, strand, n):
    return "".join(ref[idx(rel + k, strand)] if strand == "+"
                   else revcomp(ref[idx(rel + k, strand)]) for k in range(n))

# exon-1 ends (transcription orientation), from ncbiRefSeqCurated
LOCI = {
    # region : (last exonic base rel to TSS, ETS/GGAA core start rel to TSS)
    "COL3A1": dict(last_exonic=195, ggaa_start=191),
    "COL1A1": dict(last_exonic=221, ggaa_start=None),
}

rows = []
for rn, cfg in LOCI.items():
    ref, strand, anchor = window(rn)
    le = cfg["last_exonic"]
    ctx = read(ref, le - 9, strand, 20)
    print(f"\n=== {rn} exon1/intron1 junction (transcription orientation) ===")
    print(f"    TSS+{le-9}..+{le+10}: {ctx[:10]} | {ctx[10:]}")
    base_fc, base_ac = predict([ref])
    print(f"    reference fibroblast CAGE = {base_fc[0]:.4f}")

    variants = {}
    donor = read(ref, le + 1, strand, 2)
    variants["donor_kill_GT_to_AA"] = (put(ref, le + 1, strand, "AA"), f"TSS+{le+1}..+{le+2} {donor}->AA")
    variants["donor_kill_GT_to_CC"] = (put(ref, le + 1, strand, "CC"), f"TSS+{le+1}..+{le+2} {donor}->CC")
    if cfg["ggaa_start"] is not None:
        gs = cfg["ggaa_start"]
        core = read(ref, gs, strand, 4)
        variants["ets_core_kill_GGAA_to_CCTT"] = (put(ref, gs, strand, "CCTT"), f"TSS+{gs}..+{gs+3} {core}->CCTT")
        variants["ets_core_kill_GGAA_to_TTTT"] = (put(ref, gs, strand, "TTTT"), f"TSS+{gs}..+{gs+3} {core}->TTTT")
    cs = le - 25
    ctrl = read(ref, cs, strand, 4)
    variants["exonic_control_4bp"] = (put(ref, cs, strand, "".join("A" if b != "A" else "C" for b in ctrl)),
                                      f"TSS+{cs}..+{cs+3} {ctrl}->neutral")

    names = list(variants)
    fc, ac = predict([variants[n][0] for n in names])
    for n, f, a in zip(names, fc, ac):
        d = float(f - base_fc[0])
        rows.append(dict(region=rn, variant=n, description=variants[n][1],
                         ref_fib_cage=float(base_fc[0]), var_fib_cage=float(f),
                         d_fib_cage=d, pct_change_fib_cage=100 * d / float(base_fc[0]),
                         d_all_cage=float(a - base_ac[0])))
        print(f"    {n:32s} {variants[n][1]:34s} d={d:+9.3f} ({100*d/float(base_fc[0]):+7.2f}%)")

out = pd.DataFrame(rows)
out.to_csv(f"{RES}/seqreg_ext_splicedonor_control.csv", index=False)
print(f"\nwrote {RES}/seqreg_ext_splicedonor_control.csv")
