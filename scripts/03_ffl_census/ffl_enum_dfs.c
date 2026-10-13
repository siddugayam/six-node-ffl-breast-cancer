/* ffl_enum.c -- exact enumeration of n-node feed-forward-loop (FFL) modules
 *
 * Definition implemented (see 02_ffl_census.py for provenance):
 *   S is a vertex-induced connected subgraph on n nodes of the directed graph G
 *   such that there EXISTS an orientation of the reciprocal (mutual) arc pairs
 *   -- reciprocal TF<->miRNA pairs and undirected miRNA-miRNA co-transcription
 *   edges are stored as reciprocal arc pairs and are "collapsed", i.e. their
 *   effective direction is fixed by the circuit -- for which
 *     (ii)  S is a DAG,
 *     (iii) S has exactly one source and exactly one sink,
 *     (i)   S contains a 3-node FFL core R->M, R->T, M->T,
 *     (iv)  >=2 internally vertex-disjoint directed source->sink paths exist.
 *   Every directed path (core indirect arm included) may traverse at most ONE
 *   undirected edge.
 *
 * Enumerating orientations == enumerating topological orders with the source
 * first and the sink last, then orienting every edge forward.  This is complete
 * and sound: any valid orientation is recovered by its own topological order.
 *
 * build: gcc -O3 -march=native -fopenmp -o ffl_enum ffl_enum.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN   700
#define NW     11                 /* 11*64 = 704 >= MAXN */
#define MAXMOD 8

static int N;
static uint64_t OUTB[MAXN][NW], INB[MAXN][NW], UNDB[MAXN][NW];
static int ntype[MAXN];           /* 0=TF 1=Gene 2=miRNA */
static unsigned char *dF, *dB;    /* dF[s*N+v], dB[t*N+v] */

static inline int hasarc(int u,int v){ return (OUTB[u][v>>6]>>(v&63))&1ULL; }
static inline int isund (int u,int v){ return (UNDB[u][v>>6]>>(v&63))&1ULL; }

/* ---- interaction classes -------------------------------------------------
   0 TF-TF   1 TF-gene  2 TF-miRNA  3 miRNA-gene  4 miRNA-TF
   5 miRNA-miRNA        6 gene-TF   7 gene-miRNA  8 gene-gene              */
static const int ICL[3][3] = {
    /* u=TF   */ { 0, 1, 2 },
    /* u=Gene */ { 6, 8, 7 },
    /* u=miR  */ { 4, 3, 5 }
};

/* architecture key: bits0-8 iclass mask | nTF<<9 | nmiR<<13 | fflclass<<17 */
#define KEYBITS 19
#define NKEYS   (1<<KEYBITS)

static uint64_t *cnt_g;              /* [NKEYS] global counters   */
static uint64_t (*nodemask_g)[NW];   /* [NKEYS][NW] participating nodes */
static uint64_t allnodes_g[NW];
static uint64_t total_g = 0;

/* per-module scratch ------------------------------------------------------ */
typedef struct {
    int n;
    int mem[MAXMOD];
    unsigned char adj[MAXMOD];   /* out-mask, any direction present */
    unsigned char und[MAXMOD];   /* undirected out-mask */
    unsigned char inm[MAXMOD];   /* in-mask */
} Mod;

static void build_mod(Mod *m,const int *mem,int n){
    int i,j; m->n=n;
    for(i=0;i<n;i++) m->mem[i]=mem[i];
    for(i=0;i<n;i++){ m->adj[i]=0; m->und[i]=0; m->inm[i]=0; }
    for(i=0;i<n;i++) for(j=0;j<n;j++){
        if(i==j) continue;
        if(hasarc(mem[i],mem[j])){
            m->adj[i]|=(unsigned char)(1u<<j);
            m->inm[j]|=(unsigned char)(1u<<i);
            if(isund(mem[i],mem[j])) m->und[i]|=(unsigned char)(1u<<j);
        }
    }
}

/* enumerate simple s->t paths in the oriented DAG, <=1 undirected arc.
   store internal-node masks; return count (capped).                        */
static __thread int paths_int[128];
static int enum_paths(const Mod*m,const unsigned char*fw,const unsigned char*fwu,
                      int s,int t,int cur,int used_und,unsigned char intmask,
                      int *out,int *np,int cap){
    if(cur==t){ if(*np<cap) out[(*np)++]=intmask; return 0; }
    unsigned char nb=fw[cur];
    while(nb){
        int j=__builtin_ctz(nb); nb&=nb-1;
        int uu=used_und+((fwu[cur]>>j)&1);
        if(uu>1) continue;
        unsigned char im=intmask;
        if(j!=t) im|=(unsigned char)(1u<<j);
        enum_paths(m,fw,fwu,s,t,j,uu,im,out,np,cap);
        if(*np>=cap) return 0;
    }
    return 0;
}

/* ------------------------------------------------------------------------
 * Optimised admissibility test.
 *
 * Instead of enumerating all (n-2)! permutations, we enumerate TOPOLOGICAL
 * ORDERS directly by depth-first search, placing the source first and the sink
 * last.  Three conditions prune the search at every placement:
 *   - all STRICT predecessors of the candidate must already be placed, and it
 *     must have no strict arc back to a placed node (otherwise the order is
 *     cyclic);
 *   - it must have at least one in-neighbour among the placed nodes, otherwise
 *     it would be a second source, violating (iii).
 * A final sweep checks that every node except the sink has a later
 * out-neighbour (no second sink), then conditions (i) and (iv) are tested.
 * This is the same predicate as the permutation version in ffl_enum.c, and the
 * two are cross-validated against each other at n = 3, 4 and 5.
 * -------------------------------------------------------------------------- */
static int check_order(const Mod*m,const int*order,int *ndisj){
    int n=m->n,i,j,a,b,c;
    int pos[MAXMOD];
    unsigned char fw[MAXMOD],fwu[MAXMOD],bw[MAXMOD];
    for(i=0;i<n;i++) pos[order[i]]=i;
    for(i=0;i<n;i++){ fw[i]=0; fwu[i]=0; bw[i]=0; }
    for(i=0;i<n;i++){
        unsigned char nb=m->adj[i];
        while(nb){ j=__builtin_ctz(nb); nb&=nb-1;
            if(pos[i]<pos[j]){ fw[i]|=(unsigned char)(1u<<j); bw[j]|=(unsigned char)(1u<<i);
                               if((m->und[i]>>j)&1) fwu[i]|=(unsigned char)(1u<<j); } }
    }
    int S=order[0], T=order[n-1];
    for(i=0;i<n;i++){ if(i==T) continue; if(!fw[i]) return 0; }
    for(i=0;i<n;i++){ if(i==S) continue; if(!bw[i]) return 0; }
    /* (i) 3-node core, at most one undirected arc on the two-step arm */
    int core=0;
    for(a=0;a<n&&!core;a++){
        unsigned char oa=fw[a];
        unsigned char nb=oa;
        while(nb){ c=__builtin_ctz(nb); nb&=nb-1;
            unsigned char mid=oa & bw[c];
            while(mid){ b=__builtin_ctz(mid); mid&=mid-1;
                int nu=((fwu[a]>>b)&1)+((fwu[b]>>c)&1);
                if(nu<=1){ core=1; break; } }
            if(core) break; }
    }
    if(!core) return 0;
    /* (iv) >= 2 internally vertex-disjoint source->sink paths */
    int np=0;
    enum_paths(m,fw,fwu,S,T,S,0,0,paths_int,&np,64);
    if(np<2) return 0;
    int ok2=0;
    for(i=0;i<np&&!ok2;i++) for(j=i+1;j<np;j++) if(!(paths_int[i]&paths_int[j])){ ok2=1; break; }
    if(!ok2) return 0;
    if(ndisj){
        int bb=2;
        for(i=0;i<np;i++)for(j=i+1;j<np;j++){ if(paths_int[i]&paths_int[j])continue;
            for(a=j+1;a<np;a++){ if(paths_int[a]&(paths_int[i]|paths_int[j]))continue;
                if(bb<3)bb=3;
                for(b=a+1;b<np;b++){ if(paths_int[b]&(paths_int[i]|paths_int[j]|paths_int[a]))continue;
                    if(bb<4)bb=4; } } }
        *ndisj=bb;
    }
    return 1;
}

typedef struct {
    const Mod*m;
    unsigned char sIn[MAXMOD], sOut[MAXMOD];
    int order[MAXMOD];
    int *ndisj;
} TS;

static int topo_dfs(TS*t,int pos,unsigned char placed){
    const Mod*m=t->m; int n=m->n;
    unsigned char full=(unsigned char)((1u<<n)-1);
    unsigned char sinkbit=(unsigned char)(1u<<(n-1));
    if(pos==n-1){
        int T=n-1;
        if(t->sIn[T] & (unsigned char)(~placed & full & ~sinkbit)) return 0;
        if(t->sOut[T] & placed) return 0;
        if(!(m->inm[T] & placed)) return 0;
        t->order[n-1]=T;
        return check_order(m,t->order,t->ndisj);
    }
    unsigned char avail=(unsigned char)(full & ~placed & ~sinkbit);
    while(avail){
        int x=__builtin_ctz(avail); avail&=avail-1;
        if(t->sIn[x] & (unsigned char)(~placed & full)) continue;
        if(t->sOut[x] & placed) continue;
        if(!(m->inm[x] & placed)) continue;
        t->order[pos]=x;
        if(topo_dfs(t,pos+1,(unsigned char)(placed|(1u<<x)))) return 1;
    }
    return 0;
}

static int test_module_fixed(const Mod*m,int *ndisj){
    int n=m->n,i,j;
    unsigned char full=(unsigned char)((1u<<n)-1);
    /* cheap rejections first */
    for(i=0;i<n;i++){
        if(i!=n-1 && !m->adj[i]) return 0;      /* no out-neighbour at all */
        if(i!=0   && !m->inm[i]) return 0;      /* no in-neighbour at all  */
    }
    TS t; t.m=m; t.ndisj=ndisj;
    for(i=0;i<n;i++){ t.sIn[i]=0; t.sOut[i]=0; }
    for(i=0;i<n;i++){
        unsigned char nb=m->adj[i];
        while(nb){ j=__builtin_ctz(nb); nb&=nb-1;
            if(!((m->adj[j]>>i)&1)){ t.sOut[i]|=(unsigned char)(1u<<j);
                                     t.sIn[j]|=(unsigned char)(1u<<i); } }
    }
    if(t.sIn[0]) return 0;                       /* strict arc into the source */
    if(t.sOut[n-1]) return 0;                    /* strict arc out of the sink */
    (void)full;
    t.order[0]=0;
    return topo_dfs(&t,1,(unsigned char)1);
}

/* full test on an arbitrary node set: is it valid, and is (s,t)=(mem[0],mem[n-1])
   the lexicographically smallest valid (source,sink) pair?  (dedup)         */
static int test_and_canonical(const int *mem,int n,int *ndisj){
    Mod m; build_mod(&m,mem,n);
    if(!test_module_fixed(&m,ndisj)) return 0;
    /* Alternative (source,sink) choices can only exist for members that carry
       no STRICT incoming (resp. outgoing) arc inside S.  If index 0 is the only
       such candidate source and index n-1 the only candidate sink, this module
       cannot be reached from any other (s,t) pair and no dedup scan is needed. */
    int i,j,mem2[MAXMOD],q,p;
    unsigned char srcC=0,snkC=0;
    for(i=0;i<n;i++){
        int strict_in=0,strict_out=0;
        for(j=0;j<n;j++){ if(i==j) continue;
            if(((m.adj[j]>>i)&1) && !((m.adj[i]>>j)&1)) strict_in=1;
            if(((m.adj[i]>>j)&1) && !((m.adj[j]>>i)&1)) strict_out=1; }
        if(!strict_in) srcC|=(unsigned char)(1u<<i);
        if(!strict_out) snkC|=(unsigned char)(1u<<i);
    }
    if(srcC==(unsigned char)1 && snkC==(unsigned char)(1u<<(n-1))) return 1;
    int s0=mem[0],t0=mem[n-1];
    for(i=0;i<n;i++) for(j=0;j<n;j++){
        if(i==j) continue;
        if(!((srcC>>i)&1) || !((snkC>>j)&1)) continue;
        int s1=mem[i],t1=mem[j];
        if(s1>s0 || (s1==s0 && t1>=t0)) continue;
        p=0; mem2[p++]=s1;
        for(q=0;q<n;q++) if(q!=i&&q!=j) mem2[p++]=mem[q];
        mem2[p++]=t1;
        Mod m2; build_mod(&m2,mem2,n);
        if(test_module_fixed(&m2,NULL)) return 0;   /* counted elsewhere */
    }
    return 1;
}

static int arch_key(const int *mem,int n){
    int i,j,nTF=0,nmiR=0,icl=0,composite=0;
    for(i=0;i<n;i++){ if(ntype[mem[i]]==0)nTF++; else if(ntype[mem[i]]==2)nmiR++; }
    for(i=0;i<n;i++) for(j=0;j<n;j++){
        if(i==j) continue;
        if(hasarc(mem[i],mem[j])){
            icl|=1<<ICL[ntype[mem[i]]][ntype[mem[j]]];
            if(hasarc(mem[j],mem[i])){
                int a=ntype[mem[i]],b=ntype[mem[j]];
                if((a==0&&b==2)||(a==2&&b==0)) composite=1;
            }
        }
    }
    int fclass;   /* 0 composite, 1 TF-FFL, 2 miRNA-FFL, 3 gene-FFL */
    if(composite) fclass=0;
    else fclass = (ntype[mem[0]]==0)?1:((ntype[mem[0]]==2)?2:3);
    return (icl&0x1FF) | (nTF<<9) | (nmiR<<13) | (fclass<<17);
}

/* ------------------------------------------------------------------ main */
static int NMOD;
static FILE *finst=NULL;
static long MAXINST=0;
static char **name;

int main(int argc,char**argv){
    if(argc<4){ fprintf(stderr,"usage: %s graph.txt n out_prefix [maxinst]\n",argv[0]); return 1; }
    const char*gf=argv[1]; NMOD=atoi(argv[2]); const char*pref=argv[3];
    MAXINST = (argc>4)? atol(argv[4]) : 0;

    FILE*f=fopen(gf,"r"); if(!f){perror("graph");return 1;}
    int M,i,j,u,v,und;
    if(fscanf(f,"%d %d",&N,&M)!=2){fprintf(stderr,"bad header\n");return 1;}
    name=malloc(sizeof(char*)*N);
    char buf[256];
    for(i=0;i<N;i++){ if(fscanf(f,"%255s %d",buf,&ntype[i])!=2){fprintf(stderr,"bad node %d\n",i);return 1;}
                      name[i]=strdup(buf); }
    for(i=0;i<M;i++){ if(fscanf(f,"%d %d %d",&u,&v,&und)!=3){fprintf(stderr,"bad arc %d\n",i);return 1;}
        OUTB[u][v>>6]|=1ULL<<(v&63); INB[v][u>>6]|=1ULL<<(u&63);
        if(und) UNDB[u][v>>6]|=1ULL<<(v&63); }
    fclose(f);
    fprintf(stderr,"[C] N=%d M=%d n=%d\n",N,M,NMOD);

    /* BFS distances, depth cap n-1 */
    int L=NMOD-1;
    dF=malloc((size_t)N*N); dB=malloc((size_t)N*N);
    memset(dF,255,(size_t)N*N); memset(dB,255,(size_t)N*N);
    #pragma omp parallel for schedule(dynamic)
    for(int s=0;s<N;s++){
        unsigned char*d=dF+(size_t)s*N; d[s]=0;
        int *fr=malloc(sizeof(int)*N),*nx=malloc(sizeof(int)*N);
        int nf=0,nn; fr[nf++]=s;
        for(int st=1;st<=L&&nf;st++){ nn=0;
            for(int a=0;a<nf;a++){ int x=fr[a];
                for(int w=0;w<NW;w++){ uint64_t bb=OUTB[x][w];
                    while(bb){ int b=__builtin_ctzll(bb); bb&=bb-1; int y=w*64+b;
                        if(y<N&&d[y]==255){ d[y]=st; nx[nn++]=y; } } } }
            memcpy(fr,nx,sizeof(int)*nn); nf=nn; }
        free(fr);free(nx);
    }
    #pragma omp parallel for schedule(dynamic)
    for(int t=0;t<N;t++){
        unsigned char*d=dB+(size_t)t*N; d[t]=0;
        int *fr=malloc(sizeof(int)*N),*nx=malloc(sizeof(int)*N);
        int nf=0,nn; fr[nf++]=t;
        for(int st=1;st<=L&&nf;st++){ nn=0;
            for(int a=0;a<nf;a++){ int x=fr[a];
                for(int w=0;w<NW;w++){ uint64_t bb=INB[x][w];
                    while(bb){ int b=__builtin_ctzll(bb); bb&=bb-1; int y=w*64+b;
                        if(y<N&&d[y]==255){ d[y]=st; nx[nn++]=y; } } } }
            memcpy(fr,nx,sizeof(int)*nn); nf=nn; }
        free(fr);free(nx);
    }
    fprintf(stderr,"[C] distances done\n");

    cnt_g=calloc(NKEYS,sizeof(uint64_t));
    nodemask_g=calloc(NKEYS,sizeof(uint64_t)*NW);

    char instf[512]; snprintf(instf,sizeof(instf),"%s_instances.tsv",pref);
    if(MAXINST>0){ finst=fopen(instf,"w");
        fprintf(finst,"n\tsource\tsink\tmembers\tn_disjoint_paths\tarch_key\n"); }
    long inst_written=0;

    int K=NMOD-2;
    uint64_t grand=0;
    #pragma omp parallel reduction(+:grand)
    {
        int32_t *kidx=malloc(sizeof(int32_t)*NKEYS);
        for(int a=0;a<NKEYS;a++) kidx[a]=-1;
        int kcap=1024,kn=0;
        uint64_t *kcnt=calloc(kcap,sizeof(uint64_t));
        int *kkey=malloc(sizeof(int)*kcap);
        uint64_t (*knm)[NW]=calloc(kcap,sizeof(uint64_t)*NW);
        int *U=malloc(sizeof(int)*N);
        int mem[MAXMOD], sel[MAXMOD];
        #pragma omp for schedule(dynamic,1)
        for(int s=0;s<N;s++){
            const unsigned char*ds=dF+(size_t)s*N;
            for(int t=0;t<N;t++){
                if(t==s) continue;
                if(ds[t]==255) continue;
                if(hasarc(t,s)&&!hasarc(s,t)) continue;
                const unsigned char*dt=dB+(size_t)t*N;
                int nu=0;
                for(int w=0;w<N;w++){
                    if(w==s||w==t) continue;
                    if(ds[w]==255||dt[w]==255) continue;
                    if(ds[w]+dt[w]>L) continue;
                    if(hasarc(w,s)&&!hasarc(s,w)) continue;
                    if(hasarc(t,w)&&!hasarc(w,t)) continue;
                    U[nu++]=w;
                }
                if(nu<K) continue;
                int rmin = hasarc(s,t)?1:2;
                /* suffix availability for pruning */
                /* iterative K-nested loop */
                int idx[MAXMOD]; for(int a=0;a<K;a++) idx[a]=a;
                if(K==0) continue;
                int d=0; idx[0]=-1;
                int nsrc=0,nsnk=0;
                while(d>=0){
                    idx[d]++;
                    if(idx[d] > nu-(K-d)){ d--;
                        if(d>=0){ int w=U[idx[d]];
                                  if(hasarc(s,w)) nsrc--;
                                  if(hasarc(w,t)) nsnk--; }
                        continue; }
                    int w=U[idx[d]];
                    if(hasarc(s,w)) nsrc++;
                    if(hasarc(w,t)) nsnk++;
                    if(d==K-1){
                        if(nsrc>=rmin && nsnk>=rmin){
                            mem[0]=s;
                            for(int a=0;a<K;a++) mem[a+1]=U[idx[a]];
                            mem[K+1]=t;
                            int nd=0;
                            if(test_and_canonical(mem,NMOD,&nd)){
                                int key=arch_key(mem,NMOD);
                                int ki=kidx[key];
                                if(ki<0){
                                    if(kn==kcap){ int nc=kcap*2;
                                        kcnt=realloc(kcnt,sizeof(uint64_t)*nc);
                                        memset(kcnt+kcap,0,sizeof(uint64_t)*(nc-kcap));
                                        kkey=realloc(kkey,sizeof(int)*nc);
                                        knm=realloc(knm,sizeof(uint64_t)*NW*nc);
                                        memset(knm+kcap,0,sizeof(uint64_t)*NW*(nc-kcap));
                                        kcap=nc; }
                                    ki=kn++; kidx[key]=ki; kkey[ki]=key; }
                                kcnt[ki]++; grand++;
                                for(int a=0;a<NMOD;a++){ int x=mem[a];
                                    knm[ki][x>>6]|=1ULL<<(x&63); }
                                if(MAXINST>0 && inst_written<MAXINST){
                                    #pragma omp critical
                                    { if(inst_written<MAXINST){
                                        fprintf(finst,"%d\t%s\t%s\t",NMOD,name[s],name[t]);
                                        for(int a=1;a<=K;a++) fprintf(finst,"%s%s",name[mem[a]],a<K?";":"");
                                        fprintf(finst,"\t%d\t%d\n",nd,key);
                                        inst_written++; } }
                                }
                            }
                        }
                        if(hasarc(s,w)) nsrc--;
                        if(hasarc(w,t)) nsnk--;
                    } else {
                        /* prune: max reachable rmin counts */
                        int rem=K-1-d;
                        if(nsrc+rem<rmin || nsnk+rem<rmin){
                            if(hasarc(s,w)) nsrc--;
                            if(hasarc(w,t)) nsnk--;
                            continue;
                        }
                        d++; idx[d]=idx[d-1];
                    }
                }
            }
        }
        #pragma omp critical
        { for(int a=0;a<kn;a++){ int key=kkey[a]; cnt_g[key]+=kcnt[a];
            for(int w=0;w<NW;w++) nodemask_g[key][w]|=knm[a][w]; } }
        free(kidx); free(kcnt); free(kkey); free(knm); free(U);
    }
    if(finst) fclose(finst);

    char cf[512]; snprintf(cf,sizeof(cf),"%s_census.tsv",pref);
    FILE*fo=fopen(cf,"w");
    fprintf(fo,"n\tarch_key\ticlass_mask\tnTF\tnGene\tnmiRNA\tffl_class\tcount\tn_distinct_nodes\n");
    uint64_t tot=0;
    for(int a=0;a<NKEYS;a++) if(cnt_g[a]){
        int icl=a&0x1FF,nTF=(a>>9)&15,nmiR=(a>>13)&15,fc=(a>>17)&3;
        int nd=0; for(int w=0;w<NW;w++) nd+=__builtin_popcountll(nodemask_g[a][w]);
        fprintf(fo,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%llu\t%d\n",NMOD,a,icl,nTF,NMOD-nTF-nmiR,nmiR,fc,
                (unsigned long long)cnt_g[a],nd);
        tot+=cnt_g[a];
    }
    fclose(fo);
    printf("TOTAL\t%d\t%llu\n",NMOD,(unsigned long long)tot);
    fprintf(stderr,"[C] n=%d total modules = %llu\n",NMOD,(unsigned long long)tot);
    return 0;
}
