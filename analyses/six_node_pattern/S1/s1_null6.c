/* s1_null6.c -- S1 of analyses/six_node_pattern: are six-node patterns over-represented?
 *
 * Graph: the census graph without the legacy miRNA-miRNA arcs ('nolegacy'), in the null format of
 * scripts/v2/v2_null.c (uncontracted, 9,226 arcs, fine class = source type -> target type), plus one
 * layer label per arc (make_s1_inputs.py):
 *   0 TF_miRNA  1 miRNA_target  2 TF_target (analysed)  3 TF_target (TRRUST census layer)
 *   4 gene_gene (deposit, analysed)  5 gene_gene (STRING census layer, oriented low->high index =
 *   alphabetical)  6 miRNA_miRNA (10 kb census layer, both directions)
 *
 * RANDOMISATION
 *  NULL-A/B/C  verbatim copy of scripts/v2/v2_null.c: same xorshift64 stream, same groups, same order
 *              of draws, each replicate restarting from the observed graph.  Replicate r is therefore
 *              the graph behind row r of the stored runs, which the three-node counts check.
 *              Layer labels travel with the arc index (a Maslov-Sneppen swap changes only the target).
 *              The curveball groups of NULL-C rewrite their arcs in (row, column) order, which
 *              scrambles the index; there each source keeps its own multiset of layer labels, re-dealt
 *              uniformly over its new arcs from a SEPARATE stream (the structural stream is untouched).
 *  NULL-L      new, core-conditional.  All 6,829 analysed-network arcs (labels 0,1,2,4) are fixed, so
 *              every typed three-node core is preserved.  Only census-layer arcs are rewired, each
 *              within its own layer and node-type class, preserving every node's degree in it:
 *                TRRUST (3)  directed Maslov-Sneppen swaps within TF->TF, Gene->TF, Gene->Gene; a swap
 *                            is rejected if it changes the number of reciprocal TF<->TF pairs among
 *                            TF_target arcs (labels 2,3 between TF-typed nodes; 53 observed)
 *                STRING (5)  undirected double-edge swaps within Gene-Gene, Gene-TF, TF-TF; each new
 *                            link is oriented low->high index as in the census graph builder
 *                10 kb (6)   undirected double-edge swaps; each link stays in both directions
 *              No swap may create a self-loop or an arc (u,v) already present in the graph.
 *              Modes: L = all three layers; LT, LG, LM = TRRUST, STRING or 10 kb only.
 *              100 swap attempts per rewired item (the v2 setting), fresh start from the observed
 *              graph at every replicate, xorshift64 seeded from the command line.
 *
 * COUNTS per graph
 *  three-node  count3() of v2_null.c, verbatim (D1-D4 transitive triangles of the contracted graph)
 *  BHAT6       literal Fig. 1d six-node patterns of Bhat et al. 2024, exactly as enumerated by
 *              analyses/six_node_pattern_networks/build_bhat_networks.py (composite, miRNA-FFL, TF-FFL:
 *              role assignments and distinct six-node sets) and the three-node composite (M1<->T1,
 *              M1->G1, T1->G1)
 *  MODEL6      COMP_C2_toggle six-node architecture of scripts/v3/dyn_models.py, sign-agnostic,
 *              non-induced, TF1 TF2 TF-typed, miR1 miR2 miRNA, G1 G2 non-miRNA, all distinct:
 *              FULL (every interaction in the n6 equations):
 *                TF1->miR1, miR1->TF1, miR1->G1, TF1->G1                          (core)
 *                G1-G2, TF1->G2, miR1->G2                                         (gene-gene layer)
 *                miR1-miR2, TF1->miR2, miR2->G1, miR2->G2, miR2->TF1               (miRNA-miRNA layer)
 *                TF1->TF2, TF2->G1, TF2->G2, TF2->miR1, TF2->miR2                  (TF-TF layer)
 *              DEF (layer-defining edges of the dyn_models.py docstring only):
 *                core + G1-G2 + miR1-miR2 + TF1->TF2 + TF2->G1 + TF2->miR1
 *              arcs: TF->miRNA = label 0, miRNA->X = 1, TF->X = 2/3, G1-G2 = 4/5, miR1-miR2 = 6
 *  CENSUS6c    number of six-node D1-D4 modules of the contracted graph containing the composite core
 *              with TF2 regulating both the target and the miRNA ('c_full' of
 *              analyses/census_and_motif_nulls/first_rerun/ffl_census_composition.c: TF1->TF2, TF1->m, TF2->m,
 *              TF2->t, TF1->t, m->t in the contracted graph, m->TF1 in the uncontracted one, t
 *              non-miRNA).  Exact: every seed (TF1,TF2,m,t) is extended by every pair of further
 *              vertices that keeps the six-set weakly connected; each six-set is tested with the
 *              verbatim acyclic() + check_module() of v2_census2.c and counted once (hash set).
 *
 * usage: s1_null6 graph labels MODE R seed spe census_reps nproc proc [listprefix]
 *   MODE  A B C L LT LG LM  (randomised replicates), OBS or OBSNOSTRING (observed graph only)
 *   rows are written for replicates r with r % nproc == proc; CENSUS6c only for r < census_reps
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 640
#define W 10
#define MAXE 20000
#define MAXC 16
#define TMIR 0
#define TTF  1
#define TGEN 2

static int NV,NE,NCLS;
static int nodetype[MAXN];
static int eu[MAXE],ev[MAXE],ec[MAXE],lab[MAXE];
static uint64_t ADJ[MAXN][W];
static uint64_t rngstate=88172645463325252ULL;
static inline uint64_t xr(void){ rngstate^=rngstate<<13; rngstate^=rngstate>>7; rngstate^=rngstate<<17; return rngstate; }
static inline int ri(int n){ return (int)(xr()%(uint64_t)n); }
/* separate stream: re-dealing of layer labels inside the NULL-C curveball groups only */
static uint64_t rng2=0x2545F4914F6CDD1DULL;
static inline uint64_t xr2(void){ rng2^=rng2<<13; rng2^=rng2>>7; rng2^=rng2<<17; return rng2; }
static inline int ri2(int n){ return (int)(xr2()%(uint64_t)n); }

static inline int has(int u,int v){ return (ADJ[u][v>>6]>>(v&63))&1ULL; }
static inline void setE(int u,int v){ ADJ[u][v>>6]|=1ULL<<(v&63); }
static inline void clrE(int u,int v){ ADJ[u][v>>6]&=~(1ULL<<(v&63)); }

/* ===================== verbatim from scripts/v2/v2_null.c (randomisation) ===================== */
typedef struct { int *idx; int n; int kind; } Group;   /* kind 0=MS swap, 1=curveball */
static Group grp[MAXC*2]; static int ngrp=0;
static int frozen[MAXE];
static long long attempts_total=0, accepts_total=0;

static void ms_swap(Group *g,long long attempts){
    if(g->n<2) return;
    for(long long t=0;t<attempts;t++){
        int i=g->idx[ri(g->n)], j=g->idx[ri(g->n)];
        if(i==j) continue;
        int a=eu[i],b=ev[i],c=eu[j],d=ev[j];
        attempts_total++;
        if(a==d||c==b) continue;
        if(has(a,d)||has(c,b)) continue;
        clrE(a,b); clrE(c,d); setE(a,d); setE(c,b);
        ev[i]=d; ev[j]=b;
        accepts_total++;
    }
}

static int rows[MAXN],nrows;
static uint64_t rowset[MAXN][W];
static uint64_t colmask[W];
static int poolbuf[MAXN*4];
/* added for label bookkeeping: each row's labels before the trades */
static int rowlabs[MAXN][MAXN*2], nrowlabs[MAXN];

static void curveball(Group *g,long long trades){
    static int rowof[MAXN];
    memset(rowof,0xff,sizeof(rowof));
    nrows=0; memset(colmask,0,sizeof(colmask));
    for(int k=0;k<g->n;k++){
        int e=g->idx[k];
        if(rowof[eu[e]]<0){ rowof[eu[e]]=nrows; rows[nrows]=eu[e]; memset(rowset[nrows],0,sizeof(uint64_t)*W); nrowlabs[nrows]=0; nrows++; }
        rowset[rowof[eu[e]]][ev[e]>>6] |= 1ULL<<(ev[e]&63);
        colmask[ev[e]>>6] |= 1ULL<<(ev[e]&63);
        rowlabs[rowof[eu[e]]][nrowlabs[rowof[eu[e]]]++] = lab[e];          /* added */
    }
    if(nrows<2) return;
    for(long long t=0;t<trades;t++){
        int i=ri(nrows), j=ri(nrows);
        if(i==j) continue;
        uint64_t X[W],Y[W];
        int a=0,b=0;
        for(int k=0;k<W;k++){ X[k]=rowset[i][k]&~rowset[j][k]; Y[k]=rowset[j][k]&~rowset[i][k]; }
        int np=0;
        for(int k=0;k<W;k++){ uint64_t x=X[k]; while(x){ int p=__builtin_ctzll(x); x&=x-1; poolbuf[np++]=k*64+p; a++; } }
        for(int k=0;k<W;k++){ uint64_t y=Y[k]; while(y){ int p=__builtin_ctzll(y); y&=y-1; poolbuf[np++]=k*64+p; b++; } }
        if(np<2) continue;
        attempts_total++; accepts_total++;
        for(int k=np-1;k>0;k--){ int m=ri(k+1); int tmp=poolbuf[k]; poolbuf[k]=poolbuf[m]; poolbuf[m]=tmp; }
        for(int k=0;k<W;k++){ rowset[i][k]&=~X[k]; rowset[i][k]&=~Y[k]; rowset[j][k]&=~X[k]; rowset[j][k]&=~Y[k]; }
        for(int k=0;k<a;k++){ int c=poolbuf[k]; rowset[i][c>>6]|=1ULL<<(c&63); }
        for(int k=a;k<np;k++){ int c=poolbuf[k]; rowset[j][c>>6]|=1ULL<<(c&63); }
    }
    for(int k=0;k<g->n;k++){ int e=g->idx[k]; clrE(eu[e],ev[e]); }
    int pos=0;
    for(int r=0;r<nrows;r++){
        int start=pos;                                                      /* added */
        for(int k=0;k<W;k++){ uint64_t x=rowset[r][k]; while(x){ int p=__builtin_ctzll(x); x&=x-1;
            int e=g->idx[pos++]; eu[e]=rows[r]; ev[e]=k*64+p; setE(eu[e],ev[e]); } }
        /* added: re-deal this source's own labels over its new arcs (separate stream) */
        int nl=nrowlabs[r];
        if(pos-start!=nl){ fprintf(stderr,"row label count mismatch\n"); exit(1); }
        for(int k=nl-1;k>0;k--){ int m=ri2(k+1); int tmp=rowlabs[r][k]; rowlabs[r][k]=rowlabs[r][m]; rowlabs[r][m]=tmp; }
        for(int k=0;k<nl;k++) lab[g->idx[start+k]]=rowlabs[r][k];
    }
    if(pos!=g->n){ fprintf(stderr,"curveball size mismatch %d vs %d\n",pos,g->n); exit(1); }
}

static uint64_t OU[MAXN][W], IU[MAXN][W];
static int redu[MAXE],redv[MAXE],nred;
static int iscomp[MAXE];

static void count3(long long *out){   /* out: total, comp, tffl, mirfl, other, nrecip */
    memset(OU,0,sizeof(OU)); memset(IU,0,sizeof(IU));
    nred=0; long long nrecip=0;
    for(int e=0;e<NE;e++){
        int u=eu[e],v=ev[e];
        int rec = has(v,u) && ((nodetype[u]==TTF&&nodetype[v]==TMIR)||(nodetype[u]==TMIR&&nodetype[v]==TTF));
        if(rec && nodetype[u]==TMIR) continue;
        redu[nred]=u; redv[nred]=v; iscomp[nred]= (rec?1:0); nred++;
        if(rec) nrecip++;
        OU[u][v>>6]|=1ULL<<(v&63); IU[v][u>>6]|=1ULL<<(u&63);
    }
    long long tot=0,cp=0,tf=0,mr=0,ot=0;
    for(int k=0;k<nred;k++){
        int u=redu[k],v=redv[k];
        if((OU[v][u>>6]>>(u&63))&1ULL) continue;
        int c=0;
        for(int w=0;w<W;w++){
            uint64_t m = OU[u][w]&OU[v][w]&~IU[u][w]&~IU[v][w];
            c += __builtin_popcountll(m);
        }
        if(!c) continue;
        tot+=c;
        if(iscomp[k]) cp+=c;
        else if(nodetype[u]==TTF && nodetype[v]==TMIR) tf+=c;
        else if(nodetype[u]==TMIR && nodetype[v]==TTF) mr+=c;
        else ot+=c;
    }
    out[0]=tot; out[1]=cp; out[2]=tf; out[3]=mr; out[4]=ot; out[5]=nrecip;
}
/* ============================================================================================== */

/* ---------------- NULL-L ---------------- */
typedef struct { int e1[MAXE], e2[MAXE]; int n; int kind; int tA, tB; } LCls;  /* kind 0 directed, 1 STRING link, 2 10 kb link */
static LCls LC[8]; static int nLC=0;
static uint64_t TTADJ[MAXN][W];     /* TF_target arcs between TF-typed nodes (labels 2,3) */
static inline int hasTT(int u,int v){ return (TTADJ[u][v>>6]>>(v&63))&1ULL; }
static long long L_att=0, L_acc=0;

static void build_Lclasses(const char *mode){
    int doT = !strcmp(mode,"L")||!strcmp(mode,"LT"), doG = !strcmp(mode,"L")||!strcmp(mode,"LG"),
        doM = !strcmp(mode,"L")||!strcmp(mode,"LM");
    nLC=0;
    if(doT){ int tt[3][2]={{TTF,TTF},{TGEN,TTF},{TGEN,TGEN}};
        for(int c=0;c<3;c++){ LCls*L=&LC[nLC++]; L->n=0; L->kind=0; L->tA=tt[c][0]; L->tB=tt[c][1];
            for(int e=0;e<NE;e++) if(lab[e]==3 && nodetype[eu[e]]==L->tA && nodetype[ev[e]]==L->tB) L->e1[L->n++]=e; } }
    if(doG){ int tt[3][2]={{TGEN,TGEN},{TGEN,TTF},{TTF,TTF}};
        for(int c=0;c<3;c++){ LCls*L=&LC[nLC++]; L->n=0; L->kind=1; L->tA=tt[c][0]; L->tB=tt[c][1];
            for(int e=0;e<NE;e++) if(lab[e]==5){
                int x=nodetype[eu[e]], y=nodetype[ev[e]];
                if((x==L->tA&&y==L->tB)||(x==L->tB&&y==L->tA)) L->e1[L->n++]=e; } } }
    if(doM){ LCls*L=&LC[nLC++]; L->n=0; L->kind=2; L->tA=TMIR; L->tB=TMIR;
        for(int e=0;e<NE;e++) if(lab[e]==6 && eu[e]<ev[e]){
            int f=-1; for(int e2=0;e2<NE;e2++) if(lab[e2]==6 && eu[e2]==ev[e] && ev[e2]==eu[e]){ f=e2; break; }
            if(f<0){ fprintf(stderr,"10 kb arc without its reverse\n"); exit(1); }
            L->e1[L->n]=e; L->e2[L->n]=f; L->n++; } }
    for(int c=0;c<nLC;c++) fprintf(stderr,"NULL-L class %d kind %d types %d-%d items %d\n",c,LC[c].kind,LC[c].tA,LC[c].tB,LC[c].n);
}

static void rebuild_TT(void){
    memset(TTADJ,0,sizeof(TTADJ));
    for(int e=0;e<NE;e++) if((lab[e]==2||lab[e]==3) && nodetype[eu[e]]==TTF && nodetype[ev[e]]==TTF && eu[e]!=ev[e])
        TTADJ[eu[e]][ev[e]>>6]|=1ULL<<(ev[e]&63);
}
static int count_TT_recip(void){
    int r=0; for(int u=0;u<NV;u++) for(int v=u+1;v<NV;v++) if(hasTT(u,v)&&hasTT(v,u)) r++; return r;
}

static void nullL(long long spe){
    rebuild_TT();
    for(int c=0;c<nLC;c++){
        LCls*L=&LC[c]; if(L->n<2) continue;
        long long att=spe*(long long)L->n;
        for(long long t=0;t<att;t++){
            int i=ri(L->n), j=ri(L->n); if(i==j) continue;
            L_att++;
            if(L->kind==0){                                  /* directed TRRUST arcs */
                int ei=L->e1[i], ej=L->e1[j];
                int a=eu[ei],b=ev[ei],cc=eu[ej],d=ev[ej];
                if(a==d||cc==b) continue;
                if(has(a,d)||has(cc,b)) continue;
                if(L->tA==TTF && L->tB==TTF){
                    int before = hasTT(b,a) + hasTT(d,cc);
                    int after  = hasTT(d,a) + hasTT(b,cc);
                    if(before!=after) continue;
                    TTADJ[a][b>>6]&=~(1ULL<<(b&63)); TTADJ[cc][d>>6]&=~(1ULL<<(d&63));
                    TTADJ[a][d>>6]|=1ULL<<(d&63); TTADJ[cc][b>>6]|=1ULL<<(b&63);
                }
                clrE(a,b); clrE(cc,d); setE(a,d); setE(cc,b);
                ev[ei]=d; ev[ej]=b; L_acc++;
            } else {                                         /* undirected links */
                int ei=L->e1[i], ej=L->e1[j];
                int p=eu[ei], q=ev[ei], r=eu[ej], s=ev[ej];
                if(L->tA!=L->tB){                            /* bipartite Gene-TF: p,r = Gene ends */
                    if(nodetype[p]!=L->tA){ int tmp=p; p=q; q=tmp; }
                    if(nodetype[r]!=L->tA){ int tmp=r; r=s; s=tmp; }
                } else if(xr()&1ULL){ int tmp=r; r=s; s=tmp; }
                /* new links {p,s} and {r,q} */
                if(p==s||r==q) continue;
                if(L->kind==1){
                    int n1a = p<s?p:s, n1b = p<s?s:p, n2a = r<q?r:q, n2b = r<q?q:r;
                    if(n1a==n2a && n1b==n2b) continue;
                    clrE(eu[ei],ev[ei]); clrE(eu[ej],ev[ej]);
                    if(has(n1a,n1b)||has(n2a,n2b)){ setE(eu[ei],ev[ei]); setE(eu[ej],ev[ej]); continue; }
                    /* no second STRING link on the same pair in the other orientation (impossible: all low->high) */
                    setE(n1a,n1b); setE(n2a,n2b);
                    eu[ei]=n1a; ev[ei]=n1b; eu[ej]=n2a; ev[ej]=n2b; L_acc++;
                } else {
                    int fi=L->e2[i], fj=L->e2[j];
                    if((p==r&&s==q)||(p==q)) continue;
                    clrE(p,q); clrE(q,p); clrE(r,s); clrE(s,r);
                    if(has(p,s)||has(s,p)||has(r,q)||has(q,r)||(p==r&&s==q)||(p==q&&s==r)){
                        setE(p,q); setE(q,p); setE(r,s); setE(s,r); continue; }
                    setE(p,s); setE(s,p); setE(r,q); setE(q,r);
                    eu[ei]=p; ev[ei]=s; eu[fi]=s; ev[fi]=p;
                    eu[ej]=r; ev[ej]=q; eu[fj]=q; ev[fj]=r; L_acc++;
                }
            }
        }
    }
}

/* ---------------- structures for the six-node counts ---------------- */
static uint64_t CO[MAXN][W], CI[MAXN][W], NB[MAXN][W];      /* contracted graph */
static uint64_t TMb[MAXN][W], MXb[MAXN][W], TXb[MAXN][W], GGa[MAXN][W], MMb[MAXN][W];
static int is_type(int v,int t){ return nodetype[v]==t; }
static inline int bit(const uint64_t *b,int v){ return (b[v>>6]>>(v&63))&1ULL; }
static int tolist(const uint64_t *b,int *out){ int n=0; for(int k=0;k<W;k++){ uint64_t x=b[k]; while(x){ out[n++]=k*64+__builtin_ctzll(x); x&=x-1; } } return n; }

static void build_structs(int nostring){
    memset(CO,0,sizeof(CO)); memset(CI,0,sizeof(CI)); memset(NB,0,sizeof(NB));
    memset(TMb,0,sizeof(TMb)); memset(MXb,0,sizeof(MXb)); memset(TXb,0,sizeof(TXb));
    memset(GGa,0,sizeof(GGa)); memset(MMb,0,sizeof(MMb));
    for(int e=0;e<NE;e++){
        int u=eu[e], v=ev[e];
        if(nostring && lab[e]==5) continue;
        int rec = has(v,u) && ((nodetype[u]==TTF&&nodetype[v]==TMIR)||(nodetype[u]==TMIR&&nodetype[v]==TTF));
        if(!(rec && nodetype[u]==TMIR)){ CO[u][v>>6]|=1ULL<<(v&63); CI[v][u>>6]|=1ULL<<(u&63); }
        switch(lab[e]){
            case 0: TMb[u][v>>6]|=1ULL<<(v&63); break;
            case 1: MXb[u][v>>6]|=1ULL<<(v&63); break;
            case 2: case 3: if(u!=v) TXb[u][v>>6]|=1ULL<<(v&63); break;
            case 4: case 5: if(u!=v){ GGa[u][v>>6]|=1ULL<<(v&63); GGa[v][u>>6]|=1ULL<<(u&63); } break;
            case 6: if(u!=v){ MMb[u][v>>6]|=1ULL<<(v&63); MMb[v][u>>6]|=1ULL<<(u&63); } break;
        }
    }
    for(int u=0;u<NV;u++) for(int k=0;k<W;k++) NB[u][k]=CO[u][k]|CI[u][k];
}

/* ---------------- hash set of sorted six-sets ---------------- */
typedef struct { uint64_t *k; size_t cap, n; } HSet;
static void hs_init(HSet *h,size_t cap){ h->cap=cap; h->n=0; h->k=calloc(cap,sizeof(uint64_t)); }
static void hs_clear(HSet *h){ memset(h->k,0,h->cap*sizeof(uint64_t)); h->n=0; }
static int hs_add(HSet *h,uint64_t key);
static void hs_grow(HSet *h){ HSet g; hs_init(&g,h->cap*2); for(size_t i=0;i<h->cap;i++) if(h->k[i]) hs_add(&g,h->k[i]); free(h->k); *h=g; }
static int hs_add(HSet *h,uint64_t key){
    if(2*(h->n+1)>h->cap) hs_grow(h);
    uint64_t x=key*0x9E3779B97F4A7C15ULL; size_t i=(size_t)(x>>20)&(h->cap-1);
    while(h->k[i]){ if(h->k[i]==key) return 0; i=(i+1)&(h->cap-1); }
    h->k[i]=key; h->n++; return 1;
}
static uint64_t key6(const int *s,int n){
    int t[8]; memcpy(t,s,n*sizeof(int));
    for(int i=1;i<n;i++){ int v=t[i],j=i-1; while(j>=0&&t[j]>v){ t[j+1]=t[j]; j--; } t[j+1]=v; }
    uint64_t k=1ULL<<63; for(int i=0;i<n;i++) k|=((uint64_t)t[i])<<(10*i); return k;
}

/* ---------------- BHAT6 (build_bhat_networks.py, Fig. 1d literal) ---------------- */
static long long bh3_comp, bh_inst[3], bh_sets[3];      /* 0 miRNA_FFL, 1 TF_FFL, 2 Composite */
static HSet HB[3];
static FILE *LIST_BH=NULL;
static void bhat6(void){
    static int L1[MAXN],L2[MAXN],L3[MAXN],L4[MAXN],L5[MAXN];
    static uint64_t MG[MAXN][W], TG[MAXN][W], TT[MAXN][W], TTi[MAXN][W], GG[MAXN][W];
    memset(MG,0,sizeof(MG)); memset(TG,0,sizeof(TG)); memset(TT,0,sizeof(TT)); memset(TTi,0,sizeof(TTi)); memset(GG,0,sizeof(GG));
    for(int u=0;u<NV;u++){
        if(nodetype[u]==TMIR) for(int k=0;k<W;k++){ uint64_t x=MXb[u][k]; while(x){ int v=k*64+__builtin_ctzll(x); x&=x-1; if(nodetype[v]==TGEN) MG[u][v>>6]|=1ULL<<(v&63); } }
        if(nodetype[u]==TTF) for(int k=0;k<W;k++){ uint64_t x=TXb[u][k]; while(x){ int v=k*64+__builtin_ctzll(x); x&=x-1;
            if(nodetype[v]==TGEN) TG[u][v>>6]|=1ULL<<(v&63);
            if(nodetype[v]==TTF){ TT[u][v>>6]|=1ULL<<(v&63); TTi[v][u>>6]|=1ULL<<(u&63); } } }
        if(nodetype[u]==TGEN) for(int k=0;k<W;k++){ uint64_t x=GGa[u][k]; while(x){ int v=k*64+__builtin_ctzll(x); x&=x-1; if(nodetype[v]==TGEN) GG[u][v>>6]|=1ULL<<(v&63); } }
    }
    for(int c=0;c<3;c++){ hs_clear(&HB[c]); bh_inst[c]=0; }
    bh3_comp=0;
    for(int t1=0;t1<NV;t1++){ if(nodetype[t1]!=TTF) continue;
      for(int m1=0;m1<NV;m1++){ if(nodetype[m1]!=TMIR) continue;
        int mt = bit(MXb[m1],t1), tm = bit(TMb[t1],m1);
        if(!mt && !tm) continue;
        int cls_ok[3] = { mt, tm, mt&&tm };
        if(mt&&tm){ uint64_t z[W]; for(int k=0;k<W;k++) z[k]=MG[m1][k]&TG[t1][k]; for(int k=0;k<W;k++) bh3_comp+=__builtin_popcountll(z[k]); }
        int ng1=tolist(MG[m1],L1);
        for(int a=0;a<ng1;a++){ int g1=L1[a];
            int nm2=0; { int n=tolist(MMb[m1],L2); for(int b=0;b<n;b++) if(bit(MG[L2[b]],g1)) L3[nm2++]=L2[b]; }
            if(!nm2) continue;
            for(int c=0;c<3;c++){ if(!cls_ok[c]) continue;
                uint64_t t2s[W];
                for(int k=0;k<W;k++) t2s[k] = c==0? TT[t1][k] : c==1? TTi[t1][k] : (TT[t1][k]&TTi[t1][k]);
                int nt2=tolist(t2s,L4);
                for(int d=0;d<nt2;d++){ int t2=L4[d];
                    uint64_t g2s[W]; for(int k=0;k<W;k++) g2s[k]=GG[g1][k]&TG[t2][k];
                    int ng2=tolist(g2s,L5);
                    for(int e2=0;e2<nm2;e2++) for(int f=0;f<ng2;f++){
                        int s[6]={m1,L3[e2],t1,t2,g1,L5[f]};
                        bh_inst[c]++; hs_add(&HB[c],key6(s,6));
                        if(LIST_BH && c==2) fprintf(LIST_BH,"%d\t%d\t%d\t%d\t%d\t%d\n",m1,L3[e2],t1,t2,g1,L5[f]);
                    } } } } } }
    for(int c=0;c<3;c++) bh_sets[c]=HB[c].n;
}

/* ---------------- MODEL6 ---------------- */
static long long m6_inst[2], m6_sets[2];                  /* 0 FULL, 1 DEF */
static HSet HM[2];
static FILE *LIST_M6[2]={NULL,NULL};
static void model6(void){
    static int L1[MAXN],L2[MAXN],L3[MAXN],L4[MAXN];
    for(int c=0;c<2;c++){ hs_clear(&HM[c]); m6_inst[c]=0; }
    for(int t1=0;t1<NV;t1++){ if(nodetype[t1]!=TTF) continue;
      for(int m1=0;m1<NV;m1++){ if(nodetype[m1]!=TMIR) continue;
        if(!(bit(TMb[t1],m1) && bit(MXb[m1],t1))) continue;             /* TF1 <-> miR1 */
        uint64_t g1s[W]; for(int k=0;k<W;k++) g1s[k]=MXb[m1][k]&TXb[t1][k];
        int ng1=tolist(g1s,L1);
        for(int a=0;a<ng1;a++){ int g1=L1[a]; if(g1==t1||nodetype[g1]==TMIR) continue;
            int ng2=tolist(GGa[g1],L2);
            int nm2=tolist(MMb[m1],L3);
            uint64_t t2s[W]; for(int k=0;k<W;k++) t2s[k]=TXb[t1][k];
            int nt2=tolist(t2s,L4);
            for(int b=0;b<ng2;b++){ int g2=L2[b]; if(g2==t1||g2==g1||nodetype[g2]==TMIR) continue;
                int full_g2 = bit(TXb[t1],g2) && bit(MXb[m1],g2);
                for(int c2=0;c2<nm2;c2++){ int m2=L3[c2]; if(m2==m1) continue;
                    int full_m2 = full_g2 && bit(TMb[t1],m2) && bit(MXb[m2],g1) && bit(MXb[m2],g2) && bit(MXb[m2],t1);
                    for(int d=0;d<nt2;d++){ int t2=L4[d];
                        if(nodetype[t2]!=TTF||t2==t1||t2==g1||t2==g2) continue;
                        if(!(bit(TXb[t2],g1) && bit(TMb[t2],m1))) continue;       /* DEF: TF2->G1, TF2->miR1 */
                        int s[6]={t1,t2,m1,m2,g1,g2};
                        m6_inst[1]++; hs_add(&HM[1],key6(s,6));
                        if(LIST_M6[1]) fprintf(LIST_M6[1],"%d\t%d\t%d\t%d\t%d\t%d\n",t1,t2,m1,m2,g1,g2);
                        if(full_m2 && bit(TXb[t2],g2) && bit(TMb[t2],m2)){
                            m6_inst[0]++; hs_add(&HM[0],key6(s,6));
                            if(LIST_M6[0]) fprintf(LIST_M6[0],"%d\t%d\t%d\t%d\t%d\t%d\n",t1,t2,m1,m2,g1,g2);
                        }
                    } } } } } }
    for(int c=0;c<2;c++) m6_sets[c]=HM[c].n;
}

/* ---------------- CENSUS6c: verbatim D1-D4 of v2_census2.c ---------------- */
static int adj[8][8];
static int check_module(int n){
    int i, j, it;
    int indeg[8], outdeg[8];
    for(i=0;i<n;i++){indeg[i]=0;outdeg[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){outdeg[i]++;indeg[j]++;}
    int src=-1,snk=-1,ns=0,nk=0;
    for(i=0;i<n;i++){ if(indeg[i]==0){src=i;ns++;} if(outdeg[i]==0){snk=i;nk++;} }
    if(ns!=1||nk!=1) return 0;
    if(src==snk) return 0;
    int reach[8],r2[8],ch;
    for(i=0;i<n;i++){reach[i]=0;r2[i]=0;}
    reach[src]=1; do{ch=0; for(i=0;i<n;i++) if(reach[i]) for(j=0;j<n;j++) if(adj[i][j]&&!reach[j]){reach[j]=1;ch=1;} }while(ch);
    r2[snk]=1; do{ch=0; for(i=0;i<n;i++) if(r2[i]) for(j=0;j<n;j++) if(adj[j][i]&&!r2[j]){r2[j]=1;ch=1;} }while(ch);
    for(i=0;i<n;i++) if(!reach[i]||!r2[i]) return 0;
    {
        int NN=2*n, cap[16][16]; memset(cap,0,sizeof(cap));
        for(i=0;i<n;i++) cap[2*i][2*i+1] = (i==src||i==snk)? n : 1;
        for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]) cap[2*i+1][2*j]=1;
        int S=2*src+1, T=2*snk, flow=0;
        for(it=0; it<3; it++){
            int prev[16],qq[16],h=0,t2=0,u,v;
            for(u=0;u<NN;u++) prev[u]=-1;
            prev[S]=S; qq[t2++]=S;
            while(h<t2){ u=qq[h++]; for(v=0;v<NN;v++) if(prev[v]<0&&cap[u][v]>0){prev[v]=u;qq[t2++]=v;} }
            if(prev[T]<0) break;
            v=T; while(v!=S){ u=prev[v]; cap[u][v]--; cap[v][u]++; v=u; }
            flow++;
            if(flow>=2) break;
        }
        if(flow<2) return 0;
    }
    return 1;
}
static int acyclic(int n){
    int indeg[8],i,j,removed=0,ch;
    int gone[8];
    for(i=0;i<n;i++){indeg[i]=0;gone[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]) indeg[j]++;
    do{ ch=0;
        for(i=0;i<n;i++) if(!gone[i]&&indeg[i]==0){ gone[i]=1;removed++;ch=1;
            for(j=0;j<n;j++) if(adj[i][j]) indeg[j]--; }
    }while(ch);
    return removed==n;
}
static void load_adj(const int *S,int n){
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) adj[i][j] = (i!=j) && bit(CO[S[i]],S[j]);
}
static HSet HC;
static long long c6_seeds;
static FILE *LIST_C6=NULL;
static long long census6c(void){
    static int A[MAXN], Y[MAXN];
    hs_clear(&HC); c6_seeds=0;
    for(int a=0;a<NV;a++){ if(nodetype[a]!=TTF) continue;
      for(int m=0;m<NV;m++){ if(nodetype[m]!=TMIR || !bit(CO[a],m) || !has(m,a)) continue;
        for(int b=0;b<NV;b++){ if(b==a||nodetype[b]!=TTF||!bit(CO[a],b)||!bit(CO[b],m)) continue;
          for(int t=0;t<NV;t++){ if(t==a||t==b||t==m||nodetype[t]==TMIR) continue;
            if(!(bit(CO[a],t)&&bit(CO[b],t)&&bit(CO[m],t))) continue;
            c6_seeds++;
            int S[6]={a,b,m,t,0,0};
            uint64_t un[W], seedm[W]; memset(seedm,0,sizeof(seedm));
            for(int i=0;i<4;i++) seedm[S[i]>>6]|=1ULL<<(S[i]&63);
            for(int k=0;k<W;k++) un[k]=(NB[a][k]|NB[b][k]|NB[m][k]|NB[t][k])&~seedm[k];
            load_adj(S,4); if(!acyclic(4)) continue;
            int nA=tolist(un,A);
            for(int ix=0;ix<nA;ix++){ int x=A[ix];
                S[4]=x; load_adj(S,5); if(!acyclic(5)) continue;
                /* y: later members of A, or neighbours of x outside A and the seed */
                int ny=0;
                for(int iy=ix+1;iy<nA;iy++) Y[ny++]=A[iy];
                for(int k=0;k<W;k++){ uint64_t z=NB[x][k]&~un[k]&~seedm[k]; while(z){ int y=k*64+__builtin_ctzll(z); z&=z-1; if(y!=x) Y[ny++]=y; } }
                for(int iy=0;iy<ny;iy++){
                    S[5]=Y[iy]; load_adj(S,6);
                    if(!acyclic(6) || !check_module(6)) continue;
                    if(hs_add(&HC,key6(S,6)) && LIST_C6){
                        fprintf(LIST_C6,"%d\t%d\t%d\t%d\t%d\t%d\n",S[0],S[1],S[2],S[3],S[4],S[5]); }
                } } } } } }
    return (long long)HC.n;
}

/* ---------------- main ---------------- */
static int eu0[MAXE],ev0[MAXE],lab0[MAXE]; static uint64_t ADJ0[MAXN][W];

int main(int argc,char**argv){
    if(argc<10){ fprintf(stderr,"usage: %s graph labels MODE R seed spe census_reps nproc proc [listprefix]\n",argv[0]); return 1; }
    FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
    if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3) return 1;
    for(int i=0;i<NV;i++) if(fscanf(f,"%d",&nodetype[i])!=1) return 1;
    for(int i=0;i<NE;i++){ if(fscanf(f,"%d %d %d",&eu[i],&ev[i],&ec[i])!=3) return 1; setE(eu[i],ev[i]); }
    fclose(f);
    f=fopen(argv[2],"r"); if(!f){perror("labels");return 1;}
    for(int i=0;i<NE;i++) if(fscanf(f,"%d",&lab[i])!=1){ fprintf(stderr,"label file short\n"); return 1; }
    fclose(f);
    const char *mode=argv[3];
    int R=atoi(argv[4]);
    rngstate=strtoull(argv[5],NULL,10); if(!rngstate) rngstate=88172645463325252ULL;
    long long spe=atoll(argv[6]);
    int census_reps=atoi(argv[7]), nproc=atoi(argv[8]), proc=atoi(argv[9]);
    const char *lp = argc>10? argv[10] : NULL;
    for(int c=0;c<3;c++) hs_init(&HB[c],1<<16);
    for(int c=0;c<2;c++) hs_init(&HM[c],1<<14);
    hs_init(&HC,1<<20);

    printf("rep\ttotal3\tcomp3\ttf3\tmir3\tother3\trecip\tTTrecip\tbhat3_comp\tbhat6_mirFFL_inst\tbhat6_mirFFL_sets\t"
           "bhat6_TFFFL_inst\tbhat6_TFFFL_sets\tbhat6_comp_inst\tbhat6_comp_sets\tmodel6_full_inst\tmodel6_full_sets\t"
           "model6_def_inst\tmodel6_def_sets\tcensus6c_seeds\tcensus6c\tsec\n");

    if(!strcmp(mode,"OBS")||!strcmp(mode,"OBSNOSTRING")){
        int nostr = !strcmp(mode,"OBSNOSTRING");
        char fn[1024];
        if(lp){ sprintf(fn,"%s_bhat6_composite_instances.tsv",lp); LIST_BH=fopen(fn,"w");
                sprintf(fn,"%s_model6_full_instances.tsv",lp); LIST_M6[0]=fopen(fn,"w");
                sprintf(fn,"%s_model6_def_instances.tsv",lp); LIST_M6[1]=fopen(fn,"w");
                sprintf(fn,"%s_census6c_modules.tsv",lp); LIST_C6=fopen(fn,"w"); }
        if(nostr){                                  /* drop the STRING arcs from the arc list itself */
            int k=0;
            for(int e=0;e<NE;e++){ if(lab[e]==5){ clrE(eu[e],ev[e]); continue; }
                eu[k]=eu[e]; ev[k]=ev[e]; ec[k]=ec[e]; lab[k]=lab[e]; k++; }
            fprintf(stderr,"no-STRING graph: %d of %d arcs kept\n",k,NE); NE=k;
        }
        clock_t t0=clock();
        long long o[6]; count3(o);
        rebuild_TT();
        build_structs(0); bhat6(); model6();
        long long c6 = census_reps>0 ? census6c() : -1;
        printf("obs\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.1f\n",
               o[0],o[1],o[2],o[3],o[4],o[5],count_TT_recip(),bh3_comp,bh_inst[0],bh_sets[0],bh_inst[1],bh_sets[1],bh_inst[2],bh_sets[2],
               m6_inst[0],m6_sets[0],m6_inst[1],m6_sets[1],c6_seeds,c6,(double)(clock()-t0)/CLOCKS_PER_SEC);
        return 0;
    }

    memcpy(eu0,eu,sizeof(eu)); memcpy(ev0,ev,sizeof(ev)); memcpy(ADJ0,ADJ,sizeof(ADJ)); memcpy(lab0,lab,sizeof(lab));
    char model=mode[0];
    int isL = (model=='L');
    /* ---- v2_null.c set-up, verbatim ---- */
    static int gbuf[MAXC*2][MAXE]; static int gn[MAXC*2];
    if(!isL){
        memset(frozen,0,sizeof(frozen));
        if(model=='B'||model=='C'){
            for(int e=0;e<NE;e++){
                int u=eu[e],v=ev[e];
                if(has(v,u) && ((nodetype[u]==TTF&&nodetype[v]==TMIR)||(nodetype[u]==TMIR&&nodetype[v]==TTF)))
                    frozen[e]=1;
            }
        }
        memset(gn,0,sizeof(gn));
        ngrp=NCLS;
        if(model=='C') ngrp=NCLS+2;
        for(int e=0;e<NE;e++){
            if(frozen[e]) continue;
            int g=ec[e];
            if(model=='C'){
                int st=nodetype[eu[e]], tt=nodetype[ev[e]];
                if(st==TMIR && (tt==TGEN||tt==TTF)) g=NCLS;
                else if(st==TTF && (tt==TGEN||tt==TTF)) g=NCLS+1;
            }
            gbuf[g][gn[g]++]=e;
        }
        for(int g=0;g<ngrp;g++){ grp[g].idx=gbuf[g]; grp[g].n=gn[g];
            grp[g].kind = (model=='C' && g>=NCLS) ? 1 : 0; }
        fprintf(stderr,"groups:"); for(int g=0;g<ngrp;g++) fprintf(stderr," g%d:%d%s",g,gn[g],grp[g].kind?"(cb)":"");
        fprintf(stderr,"\n");
    } else build_Lclasses(mode);
    rebuild_TT(); int ttr0=count_TT_recip();
    fprintf(stderr,"observed reciprocal TF<->TF pairs among TF_target arcs: %d\n",ttr0);

    for(int r=0;r<R;r++){
        memcpy(eu,eu0,sizeof(eu)); memcpy(ev,ev0,sizeof(ev)); memcpy(ADJ,ADJ0,sizeof(ADJ)); memcpy(lab,lab0,sizeof(lab));
        if(!isL){
            for(int g=0;g<ngrp;g++){
                if(grp[g].n<2) continue;
                if(grp[g].kind==0) ms_swap(&grp[g], spe*(long long)grp[g].n);
                else curveball(&grp[g], spe*(long long)grp[g].n);
            }
        } else {
            nullL(spe);
            rebuild_TT();
            if(count_TT_recip()!=ttr0){ fprintf(stderr,"NULL-L changed the TF<->TF reciprocity\n"); return 1; }
        }
        if(r % nproc != proc) continue;
        clock_t t0=clock();
        long long o[6]; count3(o);
        rebuild_TT(); int ttr=count_TT_recip();
        build_structs(0); bhat6(); model6();
        long long c6=-1;
        if(r<census_reps) c6=census6c(); else c6_seeds=-1;
        printf("%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.1f\n",
               r,o[0],o[1],o[2],o[3],o[4],o[5],ttr,bh3_comp,bh_inst[0],bh_sets[0],bh_inst[1],bh_sets[1],bh_inst[2],bh_sets[2],
               m6_inst[0],m6_sets[0],m6_inst[1],m6_sets[1],c6_seeds,c6,(double)(clock()-t0)/CLOCKS_PER_SEC);
        fflush(stdout);
    }
    if(!isL) fprintf(stderr,"swap attempts=%lld accepted=%lld\n",attempts_total,accepts_total);
    else fprintf(stderr,"NULL-L attempts=%lld accepted=%lld\n",L_att,L_acc);
    return 0;
}
