#!/usr/bin/env python3
"""
Spatial transcriptomics test of the compartment argument (full version).

Sections
  10x Genomics public Visium human breast cancer:
      V1_Breast_Cancer_Block_A_Section_1 / _2 (Space Ranger 1.1.0)
      Visium_FFPE_Human_Breast_Cancer        (Space Ranger 1.3.0)
  Wu et al. 2021 Nat Genet Visium sections (Zenodo 4739739), which carry the
  authors' own per-spot pathology annotation.

Compartment definitions
  (i) marker-score tertiles - stroma markers contain NO collagen/ECM outcome gene
  (ii) the authors' pathology annotation (Wu sections only)

Outputs: results/v6/spatial_*.csv
"""
import os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats, sparse, io as spio
warnings.filterwarnings("ignore")

BASE="/path/to/revision"
CACHE=f"{BASE}/cache/v6"; RES=f"{BASE}/results/v6"; FIG=f"{BASE}/figures/v6"
os.makedirs(RES,exist_ok=True); os.makedirs(FIG,exist_ok=True)

EPI=["EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","ELF3","CLDN4","CLDN3","MUC1",
     "KRT5","KRT14","SLPI","AGR2","TACSTD2"]
STR=["DCN","LUM","FBLN1","MMP2","THY1","PDGFRA","VCAN","SFRP2","C1S","C1R","CFD",
     "SERPINF1","FBN1","MFAP4","ISLR","PDPN","CALD1"]
OUT=["COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12"]
TFS=["NFKB1","ETS1","E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"]
EXTRA=["EPCAM","PTPRC","ACTA2","DCN","MKI67","KRT8"]

AUTH_STROMA={"Stroma","Stroma + adipose tissue"}
AUTH_EPI={"Invasive cancer","Invasive cancer + lymphocytes","DCIS","Normal duct"}

def wu_sections():
    out=[]
    cm=glob.glob(f"{CACHE}/spatial_wu/**/*matrix*", recursive=True)
    roots=sorted(glob.glob(f"{CACHE}/spatial_wu/filtered_count_matrices/*"))
    for r in roots:
        if not os.path.isdir(r): continue
        samp=os.path.basename(r).replace("_filtered_count_matrix","")
        cand=glob.glob(f"{CACHE}/spatial_wu/spatial/{samp}*")+glob.glob(f"{CACHE}/spatial_wu/**/spatial/{samp}*",recursive=True)
        sp=[c for c in cand if os.path.isdir(c)]
        md=glob.glob(f"{CACHE}/spatial_wu/metadata/{samp}*metadata.csv")
        out.append((samp, r, sp[0] if sp else None, md[0] if md else None))
    return out

def read_wu(samp, mtxdir, spdir, mdfile):
    import anndata as ad
    X = spio.mmread(os.path.join(mtxdir,"matrix.mtx")).T.tocsr().astype(np.float32)
    bc = pd.read_csv(os.path.join(mtxdir,"barcodes.tsv"), header=None)[0].astype(str).values
    gn = pd.read_csv(os.path.join(mtxdir,"features.tsv"), header=None, sep="\t")
    gn = (gn[1] if gn.shape[1]>1 else gn[0]).astype(str).values
    a = ad.AnnData(X)
    a.obs_names = bc; a.var_names = gn
    a.var_names_make_unique()
    if spdir:
        tp=[x for x in os.listdir(spdir) if x.startswith("tissue_positions")]
        if tp:
            pos=pd.read_csv(os.path.join(spdir,tp[0]), header=None,
                            names=["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"])
            pos=pos.set_index("barcode")
            common=[b for b in a.obs_names if b in pos.index]
            a=a[common].copy()
            a.obsm["spatial"]=pos.loc[a.obs_names,["pxl_col","pxl_row"]].values.astype(float)
            sfp=os.path.join(spdir,"scalefactors_json.json")
            img=os.path.join(spdir,"tissue_hires_image.png")
            if os.path.exists(sfp) and os.path.exists(img):
                import matplotlib.pyplot as plt
                a.uns["spatial"]={samp:{"images":{"hires":plt.imread(img)},
                                        "scalefactors":json.load(open(sfp)),
                                        "metadata":{}}}
    if mdfile:
        md=pd.read_csv(mdfile, index_col=0)
        a.obs["Classification"]=md.reindex(a.obs_names)["Classification"].values
        a.obs["subtype"]=md.reindex(a.obs_names)["subtype"].values
    a.obs["section"]=samp
    return a

def prep(a):
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None, log1p=False)
    a=a[a.obs.total_counts>=500,:].copy()
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a

def dense(a,g):
    x=a[:,g].X
    return np.asarray(x.todense()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()

def spear(x,y):
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<40: return np.nan,np.nan,int(ok.sum())
    r,p=stats.spearmanr(x[ok],y[ok]); return r,p,int(ok.sum())

def meta_z(rs,ns):
    d=[(r,n) for r,n in zip(rs,ns) if np.isfinite(r) and n>3]
    if not d: return np.nan,np.nan,0
    z=np.arctanh([np.clip(r,-0.999999,0.999999) for r,_ in d]); w=np.array([n-3 for _,n in d],float)
    zm=(z*w).sum()/w.sum(); se=1/np.sqrt(w.sum())
    return float(np.tanh(zm)), float(2*stats.norm.sf(abs(zm/se))), len(d)

sections=[]
for s in ["V1_Breast_Cancer_Block_A_Section_1","V1_Breast_Cancer_Block_A_Section_2",
          "Visium_FFPE_Human_Breast_Cancer"]:
    p=f"{CACHE}/spatial/{s}"
    if os.path.isdir(p): sections.append(("10x_public",s,p,None))
for samp,mtx,sp,md in wu_sections(): sections.append(("Wu2021",samp,mtx,(sp,md)))
print("sections:",[(a,b) for a,b,_,_ in sections], flush=True)

qc=[]; mk=[]; cr=[]; frames=[]
for src,lib,path,aux in sections:
    if src=="10x_public":
        a=sc.read_visium(path, count_file="filtered_feature_bc_matrix.h5", library_id=lib)
        a.var_names_make_unique(); a.obs["section"]=lib; a.obs["Classification"]=np.nan
    else:
        a=read_wu(lib, path, aux[0], aux[1])
    n0=a.n_obs; a=prep(a)
    ge=[g for g in EPI if g in a.var_names]; gs=[g for g in STR if g in a.var_names]
    sc.tl.score_genes(a, ge, score_name="epi_score", random_state=0)
    sc.tl.score_genes(a, gs, score_name="stroma_score", random_state=0)
    a.obs["delta"]=a.obs.stroma_score-a.obs.epi_score
    q1,q2=np.quantile(a.obs.delta,[1/3,2/3])
    a.obs["comp_marker"]=pd.cut(a.obs.delta,[-np.inf,q1,q2,np.inf],
                                labels=["epithelial","intermediate","stromal"]).astype(str)
    cls=a.obs.get("Classification")
    a.obs["comp_author"]=[("stromal" if c in AUTH_STROMA else
                           "epithelial" if c in AUTH_EPI else "other")
                          for c in (cls if cls is not None else [np.nan]*a.n_obs)]
    genes=[g for g in OUT+TFS+EXTRA if g in a.var_names]
    E=pd.DataFrame({g:dense(a,g) for g in genes}, index=a.obs_names)
    for c in ["comp_marker","comp_author","stroma_score","epi_score","delta"]:
        E[c]=a.obs[c].values
    E["Classification"]=(cls.values if cls is not None else np.nan)
    E["section"]=lib; E["source"]=src
    if "spatial" in a.obsm: E["x"]=a.obsm["spatial"][:,0]; E["y"]=a.obsm["spatial"][:,1]
    frames.append(E)
    qc.append(dict(source=src,section=lib,spots_raw=n0,spots_used=a.n_obs,genes=a.n_vars,
        median_counts=float(np.median(a.obs.total_counts)),
        median_genes=float(np.median(a.obs.n_genes_by_counts)),
        n_marker_stromal=int((E.comp_marker=="stromal").sum()),
        n_marker_epithelial=int((E.comp_marker=="epithelial").sum()),
        n_author_stromal=int((E.comp_author=="stromal").sum()),
        n_author_epithelial=int((E.comp_author=="epithelial").sum())))

    for scheme,col in [("marker_tertile","comp_marker"),("author_annotation","comp_author")]:
        st=E[col]=="stromal"; ep=E[col]=="epithelial"
        if st.sum()<20 or ep.sum()<20: continue
        for g in genes:
            v=E[g].values
            u,pu=stats.mannwhitneyu(v[st.values],v[ep.values],alternative="two-sided")
            mk.append(dict(source=src,section=lib,scheme=scheme,gene=g,
                n_stromal=int(st.sum()),n_epithelial=int(ep.sum()),
                mean_stromal=float(v[st.values].mean()),mean_epithelial=float(v[ep.values].mean()),
                log_diff=float(v[st.values].mean()-v[ep.values].mean()),
                auc=float(u/(st.sum()*ep.sum())),p=pu,pct_detected=float((v>0).mean())))
        for o in [g for g in OUT if g in E.columns]:
            for t in [g for g in TFS if g in E.columns]:
                ra,pa,na=spear(E[o].values,E[t].values)
                rs,ps,ns=spear(E.loc[st,o].values,E.loc[st,t].values)
                re_,pe,ne=spear(E.loc[ep,o].values,E.loc[ep,t].values)
                cr.append(dict(source=src,section=lib,scheme=scheme,outcome=o,tf=t,
                    rho_all=ra,p_all=pa,n_all=na,rho_stromal=rs,p_stromal=ps,n_stromal=ns,
                    rho_epithelial=re_,p_epithelial=pe,n_epithelial=ne))
    print("done",lib,a.shape,flush=True)

pd.DataFrame(qc).to_csv(f"{RES}/spatial_section_qc.csv",index=False)
MK=pd.DataFrame(mk); MK.to_csv(f"{RES}/spatial_compartment_marker_enrichment.csv",index=False)
CR=pd.DataFrame(cr); CR.to_csv(f"{RES}/spatial_tf_collagen_correlation_bysection.csv",index=False)
S=pd.concat(frames); S.to_csv(f"{RES}/spatial_spot_level_values.csv.gz",compression="gzip")

rows=[]
for (scheme,o,t),g in CR.groupby(["scheme","outcome","tf"]):
    r=dict(scheme=scheme,outcome=o,tf=t)
    for lbl,rc,nc in [("all","rho_all","n_all"),("stromal","rho_stromal","n_stromal"),
                      ("epithelial","rho_epithelial","n_epithelial")]:
        rr,pp,k=meta_z(g[rc].values,g[nc].values)
        r[f"rho_{lbl}_meta"]=rr; r[f"p_{lbl}_meta"]=pp; r[f"n_sections_{lbl}"]=k
        r[f"n_spots_{lbl}"]=int(g[nc].sum())
    r["delta_stromal_minus_all"]=r["rho_stromal_meta"]-r["rho_all_meta"]
    r["delta_epithelial_minus_all"]=r["rho_epithelial_meta"]-r["rho_all_meta"]
    rows.append(r)
META=pd.DataFrame(rows)
tb=pd.read_csv(f"{RES}/spatial_tcga_bulk_reference_correlations.csv")
META=META.merge(tb[["outcome","tf","rho_tcga_bulk","p_tcga_bulk"]],on=["outcome","tf"],how="left")
META["sign_flip_vs_tcga_bulk"]=np.sign(META.rho_stromal_meta)!=np.sign(META.rho_tcga_bulk)
META.sort_values(["scheme","outcome","tf"]).to_csv(f"{RES}/spatial_tf_collagen_correlation_meta.csv",index=False)

rows=[]
for (scheme,g),d in MK.groupby(["scheme","gene"]):
    rr,pp,k=meta_z(2*d.auc.values-1,d.n_stromal.values+d.n_epithelial.values)
    rows.append(dict(scheme=scheme,gene=g,n_sections=len(d),mean_log_diff=d.log_diff.mean(),
        min_log_diff=d.log_diff.min(),max_log_diff=d.log_diff.max(),mean_auc=d.auc.mean(),
        n_sections_auc_gt_0p5=int((d.auc>0.5).sum()),mean_pct_detected=d.pct_detected.mean()))
pd.DataFrame(rows).sort_values(["scheme","mean_auc"],ascending=[True,False]).to_csv(
    f"{RES}/spatial_compartment_marker_enrichment_meta.csv",index=False)
print("WROTE")
