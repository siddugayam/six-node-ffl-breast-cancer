/* motif_null_tier_sign.c -- three-node FFL motif test under the null models of motif_null_models.c
 * (analyses/census_and_motif_nulls/nulls/) and under tier- and sign-restricted versions of them.
 *
 * Randomisation, counting (count3) and output are copied from motif_null_models.c. Models A, B and C are unchanged, so
 * they reproduce its replicate files line for line. The other models use a per-arc attribute file (one line per arc,
 * "tier sign": tier 0 strong, 1 weak, 2 predicted_only, -1 not a miRNA->target arc; sign 0 "+", 1 "-",
 * 2 unannotated, -1 not a TF->{TF,Gene} arc):
 *   BT   NULL-B freezing; free miRNA->{Gene,TF} arcs randomised by curveball within each tier (3 bipartite groups);
 *        all other classes swapped within class as in NULL-B.
 *   BS   NULL-B with the TF->Gene and TF->TF classes split by annotated sign (Maslov-Sneppen within class x sign).
 *   BTS  BT and BS together.
 *   CT   NULL-C with its miRNA-side curveball group split by tier.
 *   CTS  CT with its TF-side curveball group split by sign.
 * Group order: the fine-class groups in class order, then the split groups (miRNA side by tier, then TF side by class
 * and sign).
 *
 * usage: motif_null_tier_sign graph.txt MODEL R seed swaps_per_edge [attributes.txt]
 *   MODEL one of A B C BT BS BTS CT CTS; the attribute file is required for the last five.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 640
#define W 10
#define MAXE 20000
#define MAXC 16
#define TMIR 0
#define TTF  1
#define TGEN 2

static int NV,NE,NCLS;
static int nodetype[MAXN];
static int eu[MAXE],ev[MAXE],ec[MAXE];
static uint64_t ADJ[MAXN][W];
static uint64_t rngstate=88172645463325252ULL;
static inline uint64_t xr(void){ rngstate^=rngstate<<13; rngstate^=rngstate>>7; rngstate^=rngstate<<17; return rngstate; }
static inline int ri(int n){ return (int)(xr()%(uint64_t)n); }

static inline int has(int u,int v){ return (ADJ[u][v>>6]>>(v&63))&1ULL; }
static inline void setE(int u,int v){ ADJ[u][v>>6]|=1ULL<<(v&63); }
static inline void clrE(int u,int v){ ADJ[u][v>>6]&=~(1ULL<<(v&63)); }

/* edge groups to randomise */
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

/* curveball on a group: rows = distinct sources, cols = the group's target universe */
static int rows[MAXN],nrows;
static uint64_t rowset[MAXN][W];
static uint64_t colmask[W];
static int poolbuf[MAXN*4];

static void curveball(Group *g,long long trades){
    /* build row sets */
    static int rowof[MAXN];
    memset(rowof,0xff,sizeof(rowof));
    nrows=0; memset(colmask,0,sizeof(colmask));
    for(int k=0;k<g->n;k++){
        int e=g->idx[k];
        if(rowof[eu[e]]<0){ rowof[eu[e]]=nrows; rows[nrows]=eu[e]; memset(rowset[nrows],0,sizeof(uint64_t)*W); nrows++; }
        rowset[rowof[eu[e]]][ev[e]>>6] |= 1ULL<<(ev[e]&63);
        colmask[ev[e]>>6] |= 1ULL<<(ev[e]&63);
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
    /* write back into ADJ and edge list */
    for(int k=0;k<g->n;k++){ int e=g->idx[k]; clrE(eu[e],ev[e]); }
    int pos=0;
    for(int r=0;r<nrows;r++){
        for(int k=0;k<W;k++){ uint64_t x=rowset[r][k]; while(x){ int p=__builtin_ctzll(x); x&=x-1;
            int e=g->idx[pos++]; eu[e]=rows[r]; ev[e]=k*64+p; setE(eu[e],ev[e]); } }
    }
    if(pos!=g->n){ fprintf(stderr,"curveball size mismatch %d vs %d\n",pos,g->n); exit(1); }
}

/* ---- counting ---- */
static uint64_t OU[MAXN][W], IU[MAXN][W];
static int redu[MAXE],redv[MAXE],nred;
static int iscomp[MAXE];

static void count3(long long *out){   /* out: total, comp, tffl, mirfl, other, nrecip */
    memset(OU,0,sizeof(OU)); memset(IU,0,sizeof(IU));
    nred=0; long long nrecip=0;
    for(int e=0;e<NE;e++){
        int u=eu[e],v=ev[e];
        int rec = has(v,u) && ((nodetype[u]==TTF&&nodetype[v]==TMIR)||(nodetype[u]==TMIR&&nodetype[v]==TTF));
        if(rec && nodetype[u]==TMIR) continue;             /* drop miRNA->TF arc of the pair */
        redu[nred]=u; redv[nred]=v; iscomp[nred]= (rec?1:0); nred++;
        if(rec) nrecip++;
        OU[u][v>>6]|=1ULL<<(v&63); IU[v][u>>6]|=1ULL<<(u&63);
    }
    long long tot=0,cp=0,tf=0,mr=0,ot=0;
    for(int k=0;k<nred;k++){
        int u=redu[k],v=redv[k];
        if((OU[v][u>>6]>>(u&63))&1ULL) continue;           /* residual 2-cycle -> not induced-acyclic */
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

int main(int argc,char**argv){
    if(argc<6){ fprintf(stderr,"usage: %s graph.txt MODEL R seed swaps_per_edge [attributes]\n",argv[0]); return 1; }
    FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
    if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3) return 1;
    for(int i=0;i<NV;i++) if(fscanf(f,"%d",&nodetype[i])!=1) return 1;
    for(int i=0;i<NE;i++){ if(fscanf(f,"%d %d %d",&eu[i],&ev[i],&ec[i])!=3) return 1; setE(eu[i],ev[i]); }
    fclose(f);
    const char *M=argv[2];
    char model=M[0];
    int tierB = !strcmp(M,"BT")||!strcmp(M,"BTS"), signB = !strcmp(M,"BS")||!strcmp(M,"BTS");
    int tierC = !strcmp(M,"CT")||!strcmp(M,"CTS"), signC = !strcmp(M,"CTS");
    int plain = !strcmp(M,"A")||!strcmp(M,"B")||!strcmp(M,"C");
    if(!plain && !tierB && !signB && !tierC){ fprintf(stderr,"unknown model %s\n",M); return 1; }
    static int tier[MAXE], sgn[MAXE];
    if(!plain){
        if(argc<7){ fprintf(stderr,"attribute file required for model %s\n",M); return 1; }
        FILE *fa=fopen(argv[6],"r"); if(!fa){ perror(argv[6]); return 1; }
        for(int i=0;i<NE;i++) if(fscanf(fa,"%d %d",&tier[i],&sgn[i])!=2){ fprintf(stderr,"attribute file short\n"); return 1; }
        fclose(fa);
    }
    int R=atoi(argv[3]);
    rngstate=strtoull(argv[4],NULL,10); if(!rngstate) rngstate=88172645463325252ULL;
    long long spe=atoll(argv[5]);

    /* observed */
    long long obs[6]; 
    static int eu0[MAXE],ev0[MAXE]; static uint64_t ADJ0[MAXN][W];
    memcpy(eu0,eu,sizeof(eu)); memcpy(ev0,ev,sizeof(ev)); memcpy(ADJ0,ADJ,sizeof(ADJ));
    count3(obs);
    fprintf(stderr,"OBSERVED total=%lld comp=%lld tf=%lld mir=%lld other=%lld recip=%lld\n",
            obs[0],obs[1],obs[2],obs[3],obs[4],obs[5]);

    /* freeze mask */
    memset(frozen,0,sizeof(frozen));
    if(model=='B'||model=='C'){
        for(int e=0;e<NE;e++){
            int u=eu[e],v=ev[e];
            if(has(v,u) && ((nodetype[u]==TTF&&nodetype[v]==TMIR)||(nodetype[u]==TMIR&&nodetype[v]==TTF)))
                frozen[e]=1;
        }
    }
    /* groups */
    static int gbuf[MAXC*4][MAXE]; static int gn[MAXC*4]; static int gkind[MAXC*4];
    memset(gn,0,sizeof(gn)); memset(gkind,0,sizeof(gkind));
    ngrp=NCLS;
    if(plain){
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
        for(int g=NCLS;g<ngrp;g++) gkind[g]=1;
    } else {
        /* split groups after the class groups:
           miRNA side: 3 tier groups (curveball) if tier-split, else 1 curveball group (C models) or class groups (B)
           TF side:    C models: 3 sign groups (curveball) if sign-split, else 1 curveball group
                       B models with sign: class x sign Maslov-Sneppen groups (classes TF->Gene, TF->TF) */
        int cMIR = model=='C' || tierB;            /* miRNA side leaves the class groups */
        int cTF  = model=='C' || signB;            /* TF side leaves the class groups */
        int nm = cMIR ? ((tierB||tierC)?3:1) : 0;
        int g0m = NCLS, g0t = NCLS+nm;
        int nt = 0;
        if(model=='C') nt = signC?3:1; else if(signB) nt = 2*3;
        ngrp = NCLS + nm + nt;
        for(int g=g0m; g<g0m+nm; g++) gkind[g]=1;
        if(model=='C') for(int g=g0t; g<g0t+nt; g++) gkind[g]=1;
        /* index of the TF->Gene and TF->TF classes for the B sign groups: order of class index */
        int tfcls[2]={-1,-1}, ntfc=0;
        for(int e=0;e<NE;e++){ int st=nodetype[eu[e]], tt=nodetype[ev[e]];
            if(st==TTF && (tt==TGEN||tt==TTF)){ int c=ec[e]; if(c!=tfcls[0] && c!=tfcls[1]){ if(ntfc>=2){ fprintf(stderr,"more than two TF->{TF,Gene} classes\n"); return 1; } tfcls[ntfc++]=c; } } }
        if(ntfc==2 && tfcls[0]>tfcls[1]){ int x=tfcls[0]; tfcls[0]=tfcls[1]; tfcls[1]=x; }
        for(int e=0;e<NE;e++){
            if(frozen[e]) continue;
            int g=ec[e];
            int st=nodetype[eu[e]], tt=nodetype[ev[e]];
            if(cMIR && st==TMIR && (tt==TGEN||tt==TTF)){
                if(tierB||tierC){ if(tier[e]<0||tier[e]>2){ fprintf(stderr,"arc %d without tier\n",e); return 1; } g=g0m+tier[e]; }
                else g=g0m;
            } else if(cTF && st==TTF && (tt==TGEN||tt==TTF)){
                if(model=='C') g = signC ? g0t+sgn[e] : g0t;
                else g = g0t + (ec[e]==tfcls[0]?0:3) + sgn[e];
                if((signC||signB) && (sgn[e]<0||sgn[e]>2)){ fprintf(stderr,"arc %d without sign\n",e); return 1; }
            }
            gbuf[g][gn[g]++]=e;
        }
    }
    for(int g=0;g<ngrp;g++){ grp[g].idx=gbuf[g]; grp[g].n=gn[g]; grp[g].kind=gkind[g]; }
    fprintf(stderr,"groups:"); for(int g=0;g<ngrp;g++) fprintf(stderr," g%d:%d%s",g,gn[g],grp[g].kind?"(cb)":"");
    fprintf(stderr,"  frozen=%d\n",(int){0}+ ({int s=0;for(int e=0;e<NE;e++)s+=frozen[e];s;}));

    printf("rep\ttotal\tcomp\ttf\tmir\tother\trecip\n");
    printf("obs\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",obs[0],obs[1],obs[2],obs[3],obs[4],obs[5]);
    for(int r=0;r<R;r++){
        memcpy(eu,eu0,sizeof(eu)); memcpy(ev,ev0,sizeof(ev)); memcpy(ADJ,ADJ0,sizeof(ADJ));
        for(int g=0;g<ngrp;g++){
            if(grp[g].n<2) continue;
            if(grp[g].kind==0) ms_swap(&grp[g], spe*(long long)grp[g].n);
            else {
                /* count rows for trade budget */
                curveball(&grp[g], spe*(long long)grp[g].n);
            }
        }
        long long o[6]; count3(o);
        printf("%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\n",r,o[0],o[1],o[2],o[3],o[4],o[5]);
    }
    fprintf(stderr,"swap attempts=%lld accepted=%lld (%.3f)\n",attempts_total,accepts_total,
            attempts_total? (double)accepts_total/attempts_total:0.0);
    return 0;
}
