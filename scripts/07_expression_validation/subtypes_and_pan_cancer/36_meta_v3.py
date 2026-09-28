#!/usr/bin/env python3
"""Random-effects (DerSimonian-Laird) meta-analysis across TCGA cohorts and CPTAC cohorts,
   plus BRCA's rank within the pan-cancer distribution."""
import csv, os, math
import numpy as np
from scipy import stats
OUT="/path/to/revision/results/v2"
def rd(f):
    with open(os.path.join(OUT,f)) as fh: return list(csv.DictReader(fh))
def fl(v):
    try:
        x=float(v); return x if math.isfinite(x) else float('nan')
    except: return float('nan')
def dl(r, n):
    k=np.isfinite(r)&np.isfinite(n)&(n>4)
    r=np.clip(r[k],-0.999999,0.999999); n=n[k]
    z=np.arctanh(r); v=1.0/(n-3); w=1/v
    zf=np.sum(w*z)/np.sum(w)
    Q=float(np.sum(w*(z-zf)**2)); df=len(z)-1
    C=np.sum(w)-np.sum(w**2)/np.sum(w)
    tau2=max(0.0,(Q-df)/C) if C>0 else 0.0
    w2=1/(v+tau2); zr=np.sum(w2*z)/np.sum(w2); se=1/math.sqrt(np.sum(w2))
    p=2*stats.norm.sf(abs(zr/se))
    return dict(k=int(len(z)), rho_RE=float(np.tanh(zr)),
                lo=float(np.tanh(zr-1.96*se)), hi=float(np.tanh(zr+1.96*se)),
                p=float(p), Q=Q, df=df, p_Q=float(stats.chi2.sf(Q,df)) if df>0 else float('nan'),
                I2=float(max(0,(Q-df)/Q)*100) if Q>0 else 0.0, tau2=float(tau2))

print("="*95); print("RANDOM-EFFECTS META-ANALYSIS ACROSS 31 TCGA COHORTS (pancancer_axis.csv)"); print("="*95)
pa=rd("pancancer_axis.csv")
axes=["hsa-miR-29a -> COL1A1","hsa-miR-29b -> COL1A1","hsa-miR-29c -> COL1A1",
      "hsa-miR-29a -> COL3A1","hsa-miR-29b -> COL3A1","hsa-miR-29c -> COL3A1",
      "hsa-miR-29b-3p -> COL3A1","hsa-miR-29c-3p -> COL3A1","hsa-miR-101 -> EZH2",
      "ETS1 -> COL1A1","NFKB1 -> COL1A1","SP1 -> COL1A1","RELA -> COL1A1","MYC -> COL1A1",
      "ETS1 -> COL3A1","NFKB1 -> COL3A1","COL1A1 -> COL3A1","CAF_A -> COL1A1"]
res=[]
print(f"{'axis':28s}{'k':>3s}{'rho_RE':>9s}{'95% CI':>18s}{'p':>10s}{'I2%':>7s}{'p_Q':>10s}"
      f"{'BRCA_rho':>10s}{'BRCA_rank':>11s}")
for ax in axes:
    d=[r for r in pa if r["axis"]==ax]
    if not d: continue
    r=np.array([fl(x["rho"]) for x in d]); n=np.array([fl(x["n"]) for x in d])
    m=dl(r,n)
    brca=[x for x in d if x["cohort"]=="breast invasive carcinoma"]
    br=fl(brca[0]["rho"]) if brca else float('nan')
    srt=np.sort(r[np.isfinite(r)])
    rank=int(np.searchsorted(srt,br)+1) if math.isfinite(br) else -1
    print(f"{ax:28s}{m['k']:>3d}{m['rho_RE']:>9.3f}  [{m['lo']:+.3f},{m['hi']:+.3f}]{m['p']:>10.1e}"
          f"{m['I2']:>7.1f}{m['p_Q']:>10.1e}{br:>10.3f}{f'{rank}/{m[chr(107)]}':>11s}")
    res.append(dict(dataset="TCGA_pancancer", axis=ax, **m, rho_BRCA=br, BRCA_rank_ascending=rank))

print("\n"+"="*95); print("RANDOM-EFFECTS META-ANALYSIS ACROSS 10 CPTAC COHORTS (protein level)"); print("="*95)
cp=rd("cptac_pancancer_mir29.csv")
print(f"{'miRNA -> target':28s}{'layer':9s}{'k':>3s}{'rho_RE':>9s}{'95% CI':>18s}{'p':>10s}{'I2%':>7s}{'p_Q':>10s}")
for m_ in ["hsa-miR-29a","hsa-miR-29b","hsa-miR-29c","hsa-miR-29b-3p","hsa-miR-29c-3p"]:
    for tgt in ["COL1A1","COL3A1"]:
        d=[x for x in cp if x["miRNA"]==m_ and x["target"]==tgt]
        if not d: continue
        for layer,rc,nc in (("protein","rho_protein","n_protein"),("mRNA","rho_mrna","n_mrna")):
            r=np.array([fl(x[rc]) for x in d]); n=np.array([fl(x[nc]) for x in d])
            mm=dl(r,n)
            print(f"{m_+' -> '+tgt:28s}{layer:9s}{mm['k']:>3d}{mm['rho_RE']:>9.3f}  "
                  f"[{mm['lo']:+.3f},{mm['hi']:+.3f}]{mm['p']:>10.1e}{mm['I2']:>7.1f}{mm['p_Q']:>10.1e}")
            res.append(dict(dataset=f"CPTAC_{layer}", axis=f"{m_} -> {tgt}", **mm,
                            rho_BRCA=fl([x[rc] for x in d if x["cohort"]=="BRCA"][0]) if any(x["cohort"]=="BRCA" for x in d) else float('nan'),
                            BRCA_rank_ascending=-1))

print("\n"+"="*95); print("SUBTYPE meta (TCGA-BRCA PAM50, within-subtype random effects)"); print("="*95)
ss=rd("subtype_stratified.csv")
print(f"{'axis':28s}{'k':>3s}{'rho_RE':>9s}{'95% CI':>18s}{'p':>10s}{'I2%':>7s}{'p_Q':>10s}{'rho_ALL':>9s}")
for ax in sorted({r["axis"] for r in ss if r["cohort"]=="TCGA-BRCA"}):
    d=[r for r in ss if r["cohort"]=="TCGA-BRCA" and r["axis"]==ax and r["stratum"]!="ALL"]
    a=[r for r in ss if r["cohort"]=="TCGA-BRCA" and r["axis"]==ax and r["stratum"]=="ALL"]
    r=np.array([fl(x["rho"]) for x in d]); n=np.array([fl(x["n"]) for x in d])
    mm=dl(r,n)
    print(f"{ax:28s}{mm['k']:>3d}{mm['rho_RE']:>9.3f}  [{mm['lo']:+.3f},{mm['hi']:+.3f}]"
          f"{mm['p']:>10.1e}{mm['I2']:>7.1f}{mm['p_Q']:>10.1e}{fl(a[0]['rho']) if a else float('nan'):>9.3f}")
    res.append(dict(dataset="TCGA_BRCA_PAM50", axis=ax, **mm,
                    rho_BRCA=fl(a[0]["rho"]) if a else float('nan'), BRCA_rank_ascending=-1))

with open(os.path.join(OUT,"meta_analysis_v3.csv"),"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(res[0].keys())); w.writeheader()
    for r in res: w.writerow(r)
print("\nwrote meta_analysis_v3.csv", len(res),"rows")
