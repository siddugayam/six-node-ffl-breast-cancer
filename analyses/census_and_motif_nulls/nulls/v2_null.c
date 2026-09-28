/* v2_null.c -- motif significance for 3-node FFLs under three null models.
 *
 * Input graph = the ORIGINAL directed graph G (reciprocal TF<->miRNA arcs still present),
 * edges labelled by fine class (source_type -> target_type).
 *
 * NULL-A : degree- and class-preserving Maslov-Sneppen double-edge swaps within every
 *          fine edge class independently.  Reciprocity is NOT protected.
 * NULL-B : identical, except the arcs belonging to a reciprocal TF<->miRNA pair are frozen;
 *          only the non-reciprocal TF->miRNA and miRNA->TF arcs are swapped.
 * NULL-C : NULL-B freezing, but the bipartite miRNA->target layer (miRNA->Gene plus free
 *          miRNA->TF, one joint column universe) and the bipartite TF->target layer are
 *          randomised with curveball trades (uniform, exact joint degree sequence).
 *
 * After randomisation the reciprocal TF<->miRNA pairs are contracted to their TF->miRNA arc
 * (D1) and the induced transitive triangles of the reduced graph are counted; at n=3 that is
 * exactly the D1-D4 definition.  Composite cores are counted ONCE per reciprocal pair.
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
    if(argc<6){ fprintf(stderr,"usage: %s graph.txt NULL{A,B,C} R seed swaps_per_edge\n",argv[0]); return 1; }
    FILE*f=fopen(argv[1],"r"); if(!f){perror("open");return 1;}
    if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3) return 1;
    for(int i=0;i<NV;i++) if(fscanf(f,"%d",&nodetype[i])!=1) return 1;
    for(int i=0;i<NE;i++){ if(fscanf(f,"%d %d %d",&eu[i],&ev[i],&ec[i])!=3) return 1; setE(eu[i],ev[i]); }
    fclose(f);
    char model=argv[2][0];
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
    static int gbuf[MAXC*2][MAXE]; static int gn[MAXC*2];
    memset(gn,0,sizeof(gn));
    int gid_of_class[MAXC];
    for(int c=0;c<NCLS;c++) gid_of_class[c]=c;
    ngrp=NCLS;
    if(model=='C'){
        /* joint bipartite layers: all miRNA->* free arcs into group NCLS, all TF->{TF,Gene} into NCLS+1 */
        ngrp=NCLS+2;
    }
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
