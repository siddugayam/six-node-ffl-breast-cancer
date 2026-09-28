#!/usr/bin/env python3
"""Figures for the PDX and spatial natural experiments."""
import os, json, warnings
import numpy as np, pandas as pd, scanpy as sc
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy import sparse, stats
warnings.filterwarnings("ignore")

BASE="/path/to/revision"
CACHE=f"{BASE}/cache/v6"; RES=f"{BASE}/results/v6"; FIG=f"{BASE}/figures/v6"
os.makedirs(FIG, exist_ok=True)

BLUE="#2a78d6"; ORANGE="#eb6834"; AQUA="#1baf7a"; RED="#d03b3b"
INK="#0b0b0b"; INK2="#52514e"; GRID="#e3e2df"; SURF="#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seqblue",
      ["#f4f8fd","#cde2fb","#9ec5f4","#5598e7","#2a78d6","#184f95","#0d366b"])
DIV = LinearSegmentedColormap.from_list("divbr", [ORANGE,"#f0efec",BLUE])

mpl.rcParams.update({
 "figure.facecolor":SURF,"axes.facecolor":SURF,"savefig.facecolor":SURF,
 "font.family":"DejaVu Sans","font.size":8,
 "axes.edgecolor":GRID,"axes.linewidth":0.8,"axes.labelcolor":INK,
 "xtick.color":INK2,"ytick.color":INK2,"text.color":INK,
 "axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"grid.alpha":0.9,
 "axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False,
 "savefig.dpi":400,"savefig.bbox":"tight"})

def dense(a,g):
    x=a[:,g].X
    return np.asarray(x.todense()).ravel() if sparse.issparse(x) else np.asarray(x).ravel()

EPI=["EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","ELF3","CLDN4","CLDN3","MUC1",
     "KRT5","KRT14","SLPI","AGR2","TACSTD2"]
STR=["DCN","LUM","FBLN1","MMP2","THY1","PDGFRA","VCAN","SFRP2","C1S","C1R","CFD",
     "SERPINF1","FBN1","MFAP4","ISLR","PDPN","CALD1"]

def load(lib):
    a=sc.read_visium(f"{CACHE}/spatial/{lib}", count_file="filtered_feature_bc_matrix.h5",
                     library_id=lib)
    a.var_names_make_unique()
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None, log1p=False)
    a=a[a.obs.total_counts>=500,:].copy(); sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.tl.score_genes(a,[g for g in EPI if g in a.var_names],score_name="epi",random_state=0)
    sc.tl.score_genes(a,[g for g in STR if g in a.var_names],score_name="str",random_state=0)
    a.obs["delta"]=a.obs["str"]-a.obs["epi"]
    q1,q2=np.quantile(a.obs.delta,[1/3,2/3])
    a.obs["comp"]=pd.cut(a.obs.delta,[-np.inf,q1,q2,np.inf],
                         labels=["epithelial","intermediate","stromal"]).astype(str)
    return a

def tissue(ax, a, lib):
    im=np.asarray(a.uns["spatial"][lib]["images"]["hires"])
    g=im.mean(axis=2) if im.ndim==3 else im
    g=(g-g.min())/(g.max()-g.min()+1e-9)
    ax.imshow(1-0.35*(1-g), cmap="gray", vmin=0, vmax=1, zorder=0)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for sp in ax.spines.values(): sp.set_visible(False)
    sf=a.uns["spatial"][lib]["scalefactors"]["tissue_hires_scalef"]
    x=a.obsm["spatial"][:,0]*sf; y=a.obsm["spatial"][:,1]*sf
    mx=0.03*(x.max()-x.min()); my=0.03*(y.max()-y.min())
    ax.set_xlim(x.min()-mx, x.max()+mx); ax.set_ylim(y.max()+my, y.min()-my)
    return x, y

def gene_panel(ax,a,lib,gene,size=3.0,title=None,qmax=0.98):
    x,y=tissue(ax,a,lib); v=dense(a,gene)
    vmax=np.quantile(v,qmax); vmax=vmax if vmax>0 else max(v.max(),1e-6)
    sca=ax.scatter(x,y,c=v,cmap=SEQ,s=size,linewidths=0,
                   norm=Normalize(0,vmax),zorder=2)
    ax.set_title(title or gene, fontsize=9, color=INK, pad=4)
    cb=plt.colorbar(sca,ax=ax,fraction=0.040,pad=0.02)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6,length=2)
    cb.set_label("log norm. expr.",fontsize=6,color=INK2)
    pct=100*np.mean(v>0)
    ax.text(0.02,0.02,f"detected in {pct:.0f}% of spots",transform=ax.transAxes,
            fontsize=6,color=INK2,va="bottom")
    return sca

# ---------------------------------------------------------------- FIGURE 1
libs=["Visium_FFPE_Human_Breast_Cancer","V1_Breast_Cancer_Block_A_Section_1"]
A={l:load(l) for l in libs}
lab={"Visium_FFPE_Human_Breast_Cancer":"10x Visium FFPE human breast cancer",
     "V1_Breast_Cancer_Block_A_Section_1":"10x Visium breast cancer block A, section 1"}

fig,axes=plt.subplots(2,5,figsize=(16.0,6.8))
for r,l in enumerate(libs):
    a=A[l]
    ax=axes[r,0]; x,y=tissue(ax,a,l)
    cols={"stromal":ORANGE,"epithelial":BLUE,"intermediate":"#c9c8c4"}
    for k in ["intermediate","epithelial","stromal"]:
        m=(a.obs.comp==k).values
        ax.scatter(x[m],y[m],c=cols[k],s=5.0,linewidths=0,label=k,zorder=2)
    ax.set_title("compartment (marker tertiles)",fontsize=9,pad=4)
    if r==0: ax.legend(loc="lower left",fontsize=6,markerscale=2.2,handletextpad=0.2)
    ax.text(-0.06,0.5,lab[l],transform=ax.transAxes,rotation=90,va="center",
            ha="center",fontsize=8,color=INK2)
    for c,g in enumerate(["COL1A1","COL3A1","NFKB1","ETS1"],start=1):
        gene_panel(axes[r,c],a,l,g,size=5.0)
fig.suptitle("COL1A1 and COL3A1 mark the stromal compartment; NFKB1 does not",
             fontsize=11,y=0.985)
fig.tight_layout(rect=[0.01,0,1,0.96])
fig.savefig(f"{FIG}/fig_spatial_compartment_maps.png"); fig.savefig(f"{FIG}/fig_spatial_compartment_maps.pdf")
plt.close(fig)

# the single requested image: COL1A1 vs NFKB1 side by side
fig,axes=plt.subplots(1,2,figsize=(8.6,4.4))
a=A["Visium_FFPE_Human_Breast_Cancer"]; l=libs[0]
gene_panel(axes[0],a,l,"COL1A1",size=13.0,title="COL1A1")
gene_panel(axes[1],a,l,"NFKB1",size=13.0,title="NFKB1")
fig.suptitle("10x Visium FFPE human breast cancer: COL1A1 tracks the stroma, NFKB1 does not",
             fontsize=10,y=0.97)
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig(f"{FIG}/fig_spatial_COL1A1_vs_NFKB1.png"); fig.savefig(f"{FIG}/fig_spatial_COL1A1_vs_NFKB1.pdf")
plt.close(fig)

# ---------------------------------------------------------------- FIGURE 2
mk=pd.read_csv(f"{RES}/spatial_compartment_marker_enrichment.csv")
mk=mk[mk.scheme=="marker_tertile"]
g=mk.groupby("gene").auc.agg(["mean","min","max","size"]).reset_index().sort_values("mean")
OUTG={"COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12"}
TFG={"NFKB1","ETS1","E2F1","EZH2","GATA3","BRCA1","JUN","EGR2","ESR1","SREBF1","DNMT1","E2F3"}
def cat(x): return "collagen / ECM node" if x in OUTG else ("transcription factor" if x in TFG else "cell-type marker")
g["cat"]=g.gene.map(cat)
cmap={"collagen / ECM node":ORANGE,"transcription factor":BLUE,"cell-type marker":"#8d8c88"}
fig,ax=plt.subplots(figsize=(6.2,7.0))
yy=np.arange(len(g))
ax.hlines(yy,g["min"],g["max"],color=[cmap[c] for c in g.cat],lw=1.6,alpha=0.45)
ax.scatter(g["mean"],yy,c=[cmap[c] for c in g.cat],s=42,zorder=3,edgecolor=SURF,linewidth=1.2)
ax.axvline(0.5,color=INK2,lw=1.0,ls=(0,(4,3)))
ax.set_yticks(yy); ax.set_yticklabels(g.gene,fontsize=8)
ax.set_xlabel("AUC, stromal vs epithelial spots  (0.5 = no compartment preference)")
ax.set_title("Which genes live in the stroma?\n9 Visium breast cancer sections; dot = mean, bar = range",
             fontsize=10,loc="left")
ax.grid(axis="y",visible=False)
for k,v in cmap.items(): ax.scatter([],[],c=v,s=42,label=k)
ax.legend(loc="lower right",fontsize=7)
ax.text(0.99,0.985,"→ stromal",transform=ax.transAxes,ha="right",va="top",fontsize=7,color=INK2)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_spatial_compartment_enrichment.png")
fig.savefig(f"{FIG}/fig_spatial_compartment_enrichment.pdf"); plt.close(fig)

# ---------------------------------------------------------------- FIGURE 3
meta=pd.read_csv(f"{RES}/spatial_tf_collagen_correlation_meta.csv")
d=meta[(meta.scheme=="marker_tertile")&(meta.outcome=="COL1A1")].copy()
d=d.sort_values("rho_tcga_bulk")
fig,ax=plt.subplots(figsize=(6.4,5.2))
yy=np.arange(len(d))
ax.hlines(yy,d.rho_tcga_bulk,d.rho_stromal_meta,color=GRID,lw=2.4,zorder=1)
ax.scatter(d.rho_tcga_bulk,yy,c=ORANGE,s=46,zorder=3,label="TCGA-BRCA bulk tumours (n=1,097)",
           edgecolor=SURF,linewidth=1.1)
ax.scatter(d.rho_stromal_meta,yy,c=BLUE,s=46,zorder=3,label="within stromal Visium spots (n=8,626)",
           edgecolor=SURF,linewidth=1.1)
ax.axvline(0,color=INK2,lw=1.0)
ax.set_yticks(yy); ax.set_yticklabels(d.tf,fontsize=8)
ax.set_xlabel("Spearman correlation with COL1A1")
ax.set_title("TF–COL1A1 correlation collapses once the compartment is held fixed",
             fontsize=10,loc="left")
ax.grid(axis="y",visible=False); ax.legend(loc="lower right",fontsize=7)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_spatial_correlation_collapse.png")
fig.savefig(f"{FIG}/fig_spatial_correlation_collapse.pdf"); plt.close(fig)

# ---------------------------------------------------------------- FIGURE 4
M=pd.read_csv(f"{RES}/pdx_matched_pairs_rank_change.csv").sort_values("median_paired_delta")
PRIOR={"COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12","MMP14","PLAU","MET","STAT5A",
       "CCND2","MYBL2","E2F1","E2F3","EZH2","DNMT1","BRCA1","GATA3","ESR1","JUN","EGR2","SREBF1"}
M["kind"]=np.where(M.gene.isin(PRIOR),"prioritised node","control marker")
fig,ax=plt.subplots(figsize=(7.0,8.4))
yy=np.arange(len(M))
for i,(_,r) in enumerate(M.iterrows()):
    c = ORANGE if r.median_paired_delta<0 else BLUE
    ax.annotate("",xy=(r.median_pctrank_pdx,i),xytext=(r.median_pctrank_originator,i),
                arrowprops=dict(arrowstyle="-|>",color=c,lw=1.4,alpha=0.85,
                                shrinkA=0,shrinkB=0,mutation_scale=8))
ax.scatter(M.median_pctrank_originator,yy,facecolor="none",edgecolor=INK2,s=28,lw=1.0,
           zorder=3,label="patient tumour (PDMR originator)")
ax.scatter(M.median_pctrank_pdx,yy,c=[ORANGE if v<0 else BLUE for v in M.median_paired_delta],
           s=34,zorder=4,edgecolor=SURF,linewidth=1.0,label="matched PDX (human reads)")
ax.set_yticks(yy)
ax.set_yticklabels([f"{g}{' *' if k=='control marker' else ''}" for g,k in zip(M.gene,M.kind)],
                   fontsize=8)
ax.set_xlabel("within-sample percentile rank among 2,090 reference genes")
ax.set_title("Human collagen collapses when the stroma becomes mouse\n14 matched patient/PDX pairs, NCI PDMR",
             fontsize=10,loc="left")
ax.grid(axis="y",visible=False); ax.set_xlim(-2,102)
ax.legend(loc="upper left",fontsize=7)
ax.text(1.0,-0.055,"* control marker, not a prioritised node",transform=ax.transAxes,
        ha="right",fontsize=6.5,color=INK2)
fig.tight_layout(); fig.savefig(f"{FIG}/fig_pdx_matched_rank_change.png")
fig.savefig(f"{FIG}/fig_pdx_matched_rank_change.pdf"); plt.close(fig)

# ---------------------------------------------------------------- FIGURE 5
J=pd.read_csv(f"{RES}/pdx_stromality_vs_rankdrop.csv")
x=J.log2_CAF_over_CancerEpi.values; y=J.median_log2FC_panelnormalised.values
rho,p=stats.spearmanr(x,y)
fig,ax=plt.subplots(figsize=(6.6,5.4))
ax.axhline(0,color=INK2,lw=0.9); ax.axvline(0,color=INK2,lw=0.9)
col=[ORANGE if xi>0 else BLUE for xi in x]
ax.scatter(x,y,c=col,s=46,edgecolor=SURF,linewidth=1.0,zorder=3)
b=np.polyfit(x,y,1); xs=np.linspace(x.min()-0.4,x.max()+0.4,50)
ax.plot(xs,np.polyval(b,xs),color=INK2,lw=1.2,ls=(0,(5,3)),zorder=2)
show={"COL1A1","COL3A1","FN1","POSTN","PDGFRB","CXCL12","DCN","LUM","CCND2","E2F1","EZH2",
      "MYBL2","GATA3","EPCAM","KRT8","ESR1","NFKB1","ETS1","MKI67","PTPRC"}
from adjustText import adjust_text
txts=[ax.text(xi,yi,g,fontsize=7.0,color=INK) for xi,yi,g in zip(x,y,J.gene) if g in show]
adjust_text(txts,ax=ax,expand=(1.15,1.3),
            arrowprops=dict(arrowstyle="-",color="#a9a8a4",lw=0.6))
ax.set_xlabel("Wu et al. 2021 breast atlas:  log2 (CAF CPM / carcinoma-cell CPM)")
ax.set_ylabel("matched PDX vs patient tumour, median log2 fold change")
ax.set_title(f"Single-cell stromality predicts the PDX collapse\nSpearman rho = {rho:.3f}, p = {p:.1e}, n = {len(x)} genes",
             fontsize=10,loc="left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig_pdx_stromality_vs_drop.png")
fig.savefig(f"{FIG}/fig_pdx_stromality_vs_drop.pdf"); plt.close(fig)

# ---------------------------------------------------------------- FIGURE 6
# Wu et al. 2021 sections with the authors' own pathology annotation
from scipy import io as spio
import anndata as ad
def load_wu(samp):
    mtx=f"{CACHE}/spatial_wu/filtered_count_matrices/{samp}_filtered_count_matrix"
    spd=f"{CACHE}/spatial_wu/spatial/{samp}_spatial"
    X=spio.mmread(f"{mtx}/matrix.mtx").T.tocsr().astype(np.float32)
    bc=pd.read_csv(f"{mtx}/barcodes.tsv",header=None)[0].astype(str).values
    gn=pd.read_csv(f"{mtx}/features.tsv",header=None,sep="\t")
    gn=(gn[1] if gn.shape[1]>1 else gn[0]).astype(str).values
    a=ad.AnnData(X); a.obs_names=bc; a.var_names=gn; a.var_names_make_unique()
    pos=pd.read_csv(f"{spd}/tissue_positions_list.csv",header=None,
        names=["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]).set_index("barcode")
    keep=[b for b in a.obs_names if b in pos.index]; a=a[keep].copy()
    a.obsm["spatial"]=pos.loc[a.obs_names,["pxl_col","pxl_row"]].values.astype(float)
    a.uns["spatial"]={samp:{"images":{"hires":plt.imread(f"{spd}/tissue_hires_image.png")},
                            "scalefactors":json.load(open(f"{spd}/scalefactors_json.json")),
                            "metadata":{}}}
    md=pd.read_csv(f"{CACHE}/spatial_wu/metadata/{samp}_metadata.csv",index_col=0)
    a.obs["Classification"]=md.reindex(a.obs_names)["Classification"].values
    sc.pp.calculate_qc_metrics(a,inplace=True,percent_top=None,log1p=False)
    a=a[a.obs.total_counts>=500,:].copy(); sc.pp.filter_genes(a,min_cells=3)
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    return a

AUTH_STROMA={"Stroma","Stroma + adipose tissue"}
AUTH_EPI={"Invasive cancer","Invasive cancer + lymphocytes","DCIS","Normal duct"}
wulibs=["CID4535","CID44971"]
fig,axes=plt.subplots(2,4,figsize=(13.4,6.8))
for r,samp in enumerate(wulibs):
    a=load_wu(samp)
    ax=axes[r,0]; x,y=tissue(ax,a,samp)
    cl=a.obs.Classification.astype(str).values
    grp=np.array(["stromal" if c in AUTH_STROMA else
                  "carcinoma" if c in AUTH_EPI else "other/mixed" for c in cl])
    cols={"stromal":ORANGE,"carcinoma":BLUE,"other/mixed":"#c9c8c4"}
    for k in ["other/mixed","carcinoma","stromal"]:
        m=grp==k
        ax.scatter(x[m],y[m],c=cols[k],s=7.0,linewidths=0,label=k,zorder=2)
    ax.set_title("pathologist annotation (Wu et al. 2021)",fontsize=9,pad=4)
    if r==0: ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.02),ncol=3,
                       fontsize=6.5,markerscale=2.0,handletextpad=0.2,columnspacing=1.0)
    ax.text(-0.06,0.5,f"Wu et al. 2021 section {samp}",transform=ax.transAxes,rotation=90,
            va="center",ha="center",fontsize=8,color=INK2)
    for c,g in enumerate(["COL1A1","NFKB1","ETS1"],start=1):
        if g in a.var_names: gene_panel(axes[r,c],a,samp,g,size=7.0)
fig.suptitle("Author-annotated stroma and carcinoma: COL1A1 follows the annotation, NFKB1 does not",
             fontsize=11,y=0.985)
fig.tight_layout(rect=[0.01,0,1,0.96])
fig.savefig(f"{FIG}/fig_spatial_wu_author_annotation.png")
fig.savefig(f"{FIG}/fig_spatial_wu_author_annotation.pdf"); plt.close(fig)

print("figures written:", sorted(os.listdir(FIG)))
