#!/usr/bin/env python3
"""Independent verification of PART 2F/2G: streams raw counts once and accumulates
ONLY the genes of interest plus per-group total UMI. Different code path from the
deposited full-matrix pseudobulk."""
import h5py, numpy as np, pandas as pd, anndata as ad, time
REV="/path/to/revision"; CA=f"{REV}/cache/v7"; OUT=f"{REV}/results/v3"
H5=f"{CA}/hbca_global.h5ad"
A=ad.read_h5ad(H5,backed="r"); obs=A.obs
genes=A.var["feature_name"].astype(str).values
GOI=["COL1A1","COL1A2","COL3A1","COL5A1","COL11A1","EPCAM","KRT18","PTPRC","ACTA2","TAGLN","PDGFRA","PDGFRB",
     "FAP","POSTN","MMP11","CXCL12","PI16","DPT","CA12","LAMP5","HLA-DRA","CD74","TGFB1","TGFBR2",
     "NFKB1","RELA","SP1","ETS1","RUNX1","RUNX2","SRF","MYC","YAP1","WWTR1","CEBPB","JUN","FOS","EGR1","SMAD3","STAT3"]
cols={g:np.where(genes==g)[0] for g in GOI}
for g,v in cols.items():
    if len(v)!=1: print("WARNING gene",g,"maps to",len(v),"columns")
sel=np.concatenate([v for v in cols.values() if len(v)>0])
selname=np.concatenate([[g]*len(v) for g,v in cols.items() if len(v)>0])
lut=np.full(A.n_vars,-1,dtype=np.int64); lut[sel]=np.arange(len(sel))
ct=obs["cell_type"].astype(str).values; act=obs["author_cell_type"].astype(str).values; st=obs["batch"].astype(str).values
fine=np.where(pd.Series(act).str.startswith(("CAFs","Myofib","PCs","VSMC","Endo","Prolif Stromal")).values,act,ct)
grp=np.array([f"{a}||{b}" for a,b in zip(fine,st)])
gl=sorted(set(grp)); gi={g:i for i,g in enumerate(gl)}; gidx=np.array([gi[g] for g in grp])
print("cells",A.n_obs,"groups",len(gl))
f=h5py.File(H5,"r"); Xg=f["raw/X"]; indptr=Xg["indptr"][:]
acc=np.zeros((len(gl),len(sel))); tot=np.zeros(len(gl)); ndet=np.zeros((len(gl),len(sel)))
ncell=np.bincount(gidx,minlength=len(gl))
CH=50000; t0=time.time()
for s0 in range(0,A.n_obs,CH):
    s1=min(s0+CH,A.n_obs); lo,hi=indptr[s0],indptr[s1]
    dat=Xg["data"][lo:hi]; ind=Xg["indices"][lo:hi]; ip=indptr[s0:s1+1]-lo
    rows=np.repeat(np.arange(s1-s0),np.diff(ip))
    gsub=gidx[s0:s1][rows]
    np.add.at(tot,gsub,dat)
    m=lut[ind]>=0
    np.add.at(acc,(gsub[m],lut[ind[m]]),dat[m])
    np.add.at(ndet,(gsub[m],lut[ind[m]]),1.0)
    if (s0//CH)%3==0: print(f"  {s1}/{A.n_obs} {time.time()-t0:.0f}s",flush=True)
f.close()
CPM=pd.DataFrame(acc/tot[:,None]*1e6,index=gl,columns=selname)
CPM=CPM.T.groupby(level=0).sum().T
DET=pd.DataFrame(ndet/np.maximum(ncell,1)[:,None],index=gl,columns=selname).T.groupby(level=0).max().T
meta=pd.DataFrame({"group":gl,"cell":[g.split("||")[0] for g in gl],"study":[g.split("||")[1] for g in gl],
                   "n_cells":ncell,"total_umi":tot})
out=pd.concat([meta.set_index("group"),CPM],axis=1)
out.to_csv(f"{OUT}/sc_verify_hbca_goi_cpm_by_celltype_study.csv")
print("\ntotal UMI %.3e"%tot.sum())
# ---- COL1A1/COL3A1 fibroblast vs malignant per study ----
rows=[]
for s in sorted(set(meta.study)):
    sub=out[out.study==s]
    fib=sub[sub.cell.astype(str).str.contains("fibro|CAFs|Myofib",case=False,regex=True)]
    mal=sub[sub.cell=="malignant cell"]
    if len(fib)==0 or len(mal)==0: continue
    fu=fib.total_umi.sum(); mu=mal.total_umi.sum()
    r=dict(study=s,n_fib_cells=int(fib.n_cells.sum()),n_mal_cells=int(mal.n_cells.sum()))
    for g in ["COL1A1","COL3A1","EPCAM"]:
        fc=(fib[g]*fib.total_umi).sum()/fu; mc=(mal[g]*mal.total_umi).sum()/mu
        r[f"{g}_fib_CPM"]=fc; r[f"{g}_mal_CPM"]=mc; r[f"{g}_fold"]=fc/mc if mc>0 else np.inf
    rows.append(r)
R=pd.DataFrame(rows)
print("\n=== COL1A1/COL3A1 fibroblast-vs-malignant, per constituent study (independent) ===")
print(R.to_string(index=False,float_format=lambda x:"%.1f"%x))
nw=R[R.study!="wu_natgen_2021"]
print("\nexcluding Wu (n=%d studies): COL1A1 fold median %.1f (range %.1f-%.1f); COL3A1 median %.1f (range %.1f-%.1f)"%(
  len(nw),nw.COL1A1_fold.median(),nw.COL1A1_fold.min(),nw.COL1A1_fold.max(),
  nw.COL3A1_fold.median(),nw.COL3A1_fold.min(),nw.COL3A1_fold.max()))
pal=R[R.study=="pal_2021"]
if len(pal): print("Pal 2021 alone: COL1A1 %.1f vs %.1f CPM = %.1fx ; COL3A1 %.1f vs %.1f = %.1fx ; EPCAM %.1f vs %.1f"%(
  pal.COL1A1_fib_CPM.iloc[0],pal.COL1A1_mal_CPM.iloc[0],pal.COL1A1_fold.iloc[0],
  pal.COL3A1_fib_CPM.iloc[0],pal.COL3A1_mal_CPM.iloc[0],pal.COL3A1_fold.iloc[0],
  pal.EPCAM_fib_CPM.iloc[0],pal.EPCAM_mal_CPM.iloc[0]))
R.to_csv(f"{OUT}/sc_verify_collagen_caf_fold_by_study.csv",index=False)
