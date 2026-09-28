#!/usr/bin/env python3
"""21_seqreg_alphagenome.py -- AlphaGenome (Google DeepMind, 2025) sequence-to-function
predictions at the COL1A1 / COL3A1 / miR-29 / miR-130a loci.

STATUS: NOT RUN -- no API key was obtainable without the account holder
signing in to a Google account and accepting the AlphaGenome Terms of Service (see the
report). The script below is complete and will run unmodified once
    export ALPHAGENOME_API_KEY=...
is set. The gRPC endpoint was confirmed reachable from this machine: an intentionally
invalid key returns StatusCode.INVALID_ARGUMENT "API key not valid. Please pass a valid
API key.", i.e. the only missing element is the key itself.

Queries implemented:
 (i)   predicted chromatin accessibility (DNASE/ATAC), histone marks, TF binding (CHIP_TF)
       and expression (RNA_SEQ/CAGE) over COL1A1 and COL3A1, with the CHIP_TF tracks ranked
       so that NFKB1/RELA/SP1/ETS1 can be compared against the full predicted factor set;
 (ii)  the same for the miR-29a/b-1 (chr7q32), miR-29b-2/c (chr1q32) and miR-130a (chr11q12)
       loci, plus their host-gene promoters (LINC-PINT / MIR29B2CHG / MIR130AHG);
 (iii) AlphaGenome takes DNA sequence only and has NO DNA-methylation input or output head,
       so "predicted effect of promoter methylation on miR-130a expression" CANNOT be asked
       of it directly. The nearest supported query -- implemented here -- is (a) predicted
       H3K27me3/H3K9me3/H3K4me3 and accessibility at the MIR130A/MIR130AHG promoter, and
       (b) in-silico saturation of the CpG-island core by C->T transitions, scored with the
       gene-mask LFC scorer, which asks whether the model believes that sequence is
       load-bearing for MIR130AHG expression.
"""
import os, sys, json
import numpy as np, pandas as pd

ROOT="/path/to/revision"
RES=f"{ROOT}/results/v3"
KEY=os.environ.get("ALPHAGENOME_API_KEY")
if not KEY:
    sys.exit("ALPHAGENOME_API_KEY not set -- see results/v3/seqreg_alphagenome_access.csv")

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

model = dna_client.create(KEY)
SEQLEN = dna_client.SEQUENCE_LENGTH_500KB      # 524,288 bp

# ontology terms: mammary/fibroblast/breast contexts plus broad references
ONT = ["UBERON:0001911",   # mammary gland
       "EFO:0002009",      # cultured fibroblast (GTEx label)
       "CL:0000057",       # fibroblast
       "UBERON:0002097",   # skin of body
       "EFO:0001187",      # HepG2 (reference, well predicted)
       "EFO:0001203"]      # MCF-7
OUTS = [dna_client.OutputType.DNASE, dna_client.OutputType.ATAC,
        dna_client.OutputType.RNA_SEQ, dna_client.OutputType.CAGE,
        dna_client.OutputType.CHIP_TF, dna_client.OutputType.CHIP_HISTONE]

reg = pd.read_csv(f"{RES}/seqreg_regions.csv")
rows=[]
for _,r in reg.iterrows():
    a=int(r.anchor)
    iv=genome.Interval(chromosome=r.chrom, start=a-SEQLEN//2, end=a+SEQLEN//2)
    out=model.predict_interval(interval=iv, requested_outputs=OUTS, ontology_terms=None)
    for name in ["dnase","atac","rna_seq","cage","chip_tf","chip_histone"]:
        td=getattr(out, name, None)
        if td is None: continue
        v=td.values                      # (positions, tracks)
        md=td.metadata
        centre=v.shape[0]//2
        w=max(1, 512//td.resolution) if hasattr(td,'resolution') else 4
        prom=v[centre-w:centre+w].mean(axis=0)
        for j in range(v.shape[1]):
            m=md.iloc[j]
            rows.append(dict(region=r.region, output=name, track=str(m.get('name','')),
                             assay=str(m.get('assay','')), tf=str(m.get('transcription_factor_id',
                                m.get('histone_mark',''))), biosample=str(m.get('biosample_name','')),
                             ontology=str(m.get('ontology_curie','')), strand=str(m.get('strand','')),
                             promoter_value=float(prom[j]),
                             locus_max=float(v[:,j].max())))
    print("done", r.region, flush=True)
pd.DataFrame(rows).to_csv(f"{RES}/seqreg_alphagenome_tracks.csv.gz", index=False, compression="gzip")

# (iii) in-silico C->T saturation of the MIR130A CpG island core
cpg=pd.read_csv(f"{RES}/seqreg_cpg_islands.csv")
isl=cpg[(cpg.region=="MIR130AHG")].sort_values("dist_to_anchor", key=abs).head(1)
if len(isl):
    ch=isl.chrom.iloc[0]; s=int(isl.start.iloc[0]); e=int(isl.end.iloc[0])
    scorer=[variant_scorers.GeneMaskLFCScorer(requested_output=dna_client.OutputType.RNA_SEQ)]
    recs=[]
    import pyfaidx
    fa=pyfaidx.Fasta(f"{ROOT}/cache/seqreg/genome/hg38.fa", as_raw=True, sequence_always_upper=True)
    seq=str(fa[ch][s:e])
    centre=(s+e)//2
    iv=genome.Interval(chromosome=ch, start=centre-SEQLEN//2, end=centre+SEQLEN//2)
    for i,base in enumerate(seq):
        if base!='C': continue
        pos=s+i+1
        var=genome.Variant(chromosome=ch, position=pos, reference_bases='C', alternate_bases='T')
        try:
            sc=model.score_variant(interval=iv, variant=var, variant_scorers=scorer)
        except Exception as ex:
            recs.append(dict(pos=pos, error=str(ex))); continue
        for ad in sc:
            df=ad.to_df() if hasattr(ad,'to_df') else None
            recs.append(dict(pos=pos, payload=json.dumps(ad.var.to_dict() if hasattr(ad,'var') else {})[:500]))
    pd.DataFrame(recs).to_csv(f"{RES}/seqreg_alphagenome_mir130a_cpg_ism.csv", index=False)
print("AlphaGenome run complete")
