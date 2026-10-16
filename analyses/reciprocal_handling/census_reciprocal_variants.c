/* census_reciprocal_variants.c -- the topological definition census of the census graph under four treatments of the reciprocal
 * TF<->miRNA pairs (and, for the condensation variants, of every directed cycle).
 *
 * Enumeration (source-rooted ESU over out-neighbourhoods), the RAND-ESU stream (xorshift64, same order of draws) and the
 * module test (acyclic() + check_module(), i.e. D2 and D3-D4 in the paper's numbering; D1 holds by construction) are
 * copied from analyses/census_and_motif_nulls/census/ffl_census_composition.c, so mode 'base' reproduces its counts.
 *
 * Modes (structure graph used for the enumeration in brackets):
 *   base  [contracted census graph]    the published census. Each module found is also tested under the variants:
 *         U  valid on the uncontracted graph (both arcs of each reciprocal pair kept);
 *         S  retained after graph-level SCC condensation (its nodes lie in n distinct SCCs and their images form a
 *            topological-definition module of the condensation DAG);
 *         M  valid after module-level condensation (SCCs of its own induced subgraph in the uncontracted graph
 *            contracted, then D3-D4); the first failing condition is recorded.
 *   unc   [uncontracted census graph]  topological-definition modules.
 *   mod   [uncontracted census graph]  module-level condensation. No acyclicity pruning; a node set is counted from the
 *         root r = Vsub[0] only if r is the smallest-index vertex of the source component, so each set is counted once.
 *         Each module found is also tested on the contracted graph (shared with the baseline).
 *   scc   [condensation DAG]           topological-definition modules of the condensation DAG (vertices are SCCs).
 *
 * usage: census_reciprocal_variants MODE contracted.txt uncontracted.txt scc.txt K nparts part seed q1..q6 keyprefix
 *   graph files: "NV NE NCLS" / NV node types / NE lines "u v c" (the census formats of make_graphs.py)
 *   scc.txt: "NV NC NA" / NV component ids / NA lines "a b" (arcs of the condensation DAG)
 *   keyprefix: '-' for none; else for K <= 6 writes <prefix>.keys (uint64, sorted 10-bit vertex ids) and, in mode
 *              base, <prefix>.flags (uint8 per module: bit0 U valid, bit1 S retained, bit2 M valid, bits4-5 M failing
 *              condition 3 or 4) and <prefix>.images (uint64 image keys of the modules retained under S).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define W 16

static int MODE;            /* 0 base, 1 unc, 2 mod, 3 scc */
static int NV;
static uint64_t OUT[MAXN][W], IN[MAXN][W];          /* structure graph of the enumeration */
static uint64_t CADJ[MAXN][W], UADJ[MAXN][W];       /* contracted and uncontracted census graphs */
static uint64_t SADJ[MAXN][W];                      /* condensation DAG */
static int scc_of[MAXN], NC;

static int K;
static double q[8];
static uint64_t rngstate;
static inline double rnd(void){
    rngstate ^= rngstate<<13; rngstate ^= rngstate>>7; rngstate ^= rngstate<<17;
    return (double)(rngstate>>11) / 9007199254740992.0;
}
static inline int bit(uint64_t (*G)[W], int u, int v){ return (G[u][v>>6]>>(v&63))&1ULL; }

static long long found = 0, nodes_visited = 0;
static long long n_u = 0, n_s = 0, n_m = 0, n_mfail3 = 0, n_mfail4 = 0, n_shared_base = 0, n_rootskip = 0;
static FILE *FK = NULL, *FF = NULL, *FI = NULL;

static int adj[8][8];
static int Vsub[8];
static uint64_t forbid[W];

/* ---- D3-D4, verbatim logic of check_module (ffl_census_composition.c), on an explicit adjacency ---- */
static int check_module_A(int A[8][8], int n){
    int i, j, it;
    int indeg[8], outdeg[8];
    for(i=0;i<n;i++){indeg[i]=0;outdeg[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(A[i][j]){outdeg[i]++;indeg[j]++;}
    int src=-1,snk=-1,ns=0,nk=0;
    for(i=0;i<n;i++){ if(indeg[i]==0){src=i;ns++;} if(outdeg[i]==0){snk=i;nk++;} }
    if(ns!=1||nk!=1) return 3;
    if(src==snk) return 3;
    int reach[8],r2[8],ch;
    for(i=0;i<n;i++){reach[i]=0;r2[i]=0;}
    reach[src]=1; do{ch=0; for(i=0;i<n;i++) if(reach[i]) for(j=0;j<n;j++) if(A[i][j]&&!reach[j]){reach[j]=1;ch=1;} }while(ch);
    r2[snk]=1; do{ch=0; for(i=0;i<n;i++) if(r2[i]) for(j=0;j<n;j++) if(A[j][i]&&!r2[j]){r2[j]=1;ch=1;} }while(ch);
    for(i=0;i<n;i++) if(!reach[i]||!r2[i]) return 3;
    {
        int NN=2*n, cap[16][16]; memset(cap,0,sizeof(cap));
        for(i=0;i<n;i++) cap[2*i][2*i+1] = (i==src||i==snk)? n : 1;
        for(i=0;i<n;i++)for(j=0;j<n;j++) if(A[i][j]) cap[2*i+1][2*j]=1;
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
        if(flow<2) return 4;
    }
    return 0;
}
static int check_module(int n){ return check_module_A(adj, n) == 0; }

static int acyclic_A(int A[8][8], int n){
    int indeg[8],i,j,removed=0,ch;
    int gone[8];
    for(i=0;i<n;i++){indeg[i]=0;gone[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(A[i][j]) indeg[j]++;
    do{ ch=0;
        for(i=0;i<n;i++) if(!gone[i]&&indeg[i]==0){ gone[i]=1;removed++;ch=1;
            for(j=0;j<n;j++) if(A[i][j]) indeg[j]--; }
    }while(ch);
    return removed==n;
}
static int acyclic(int n){ return acyclic_A(adj, n); }

static void fill(int A[8][8], uint64_t (*G)[W], const int *S, int n){
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) A[i][j] = (i!=j) && bit(G,S[i],S[j]);
}

static uint64_t key_of(const int *S, int n){
    int s[8]; for(int i=0;i<n;i++) s[i]=S[i];
    for(int i=1;i<n;i++){ int x=s[i], j=i-1; while(j>=0 && s[j]>x){ s[j+1]=s[j]; j--; } s[j+1]=x; }
    uint64_t k=0; for(int i=0;i<n;i++) k=(k<<10)|(uint64_t)s[i];
    return k;
}

/* module-level condensation on the uncontracted graph.
 * returns 0 valid, 3 D3 fails, 4 D4 fails, 9 weakly disconnected (D1), 8 root rule fails (mode mod only) */
static int mtest(const int *S, int n, int rootrule){
    int A[8][8], R[8][8], i, j, k;
    fill(A, UADJ, S, n);
    /* D1: weak connectivity */
    { int seen[8]={0}, st[8], sp=0; seen[0]=1; st[sp++]=0;
      while(sp){ int u=st[--sp]; for(int v=0;v<n;v++) if(!seen[v] && (A[u][v]||A[v][u])){ seen[v]=1; st[sp++]=v; } }
      for(i=0;i<n;i++) if(!seen[i]) return 9; }
    for(i=0;i<n;i++) for(j=0;j<n;j++) R[i][j] = (i==j) || A[i][j];
    for(k=0;k<n;k++) for(i=0;i<n;i++) if(R[i][k]) for(j=0;j<n;j++) if(R[k][j]) R[i][j]=1;
    int comp[8], nc=0;
    for(i=0;i<n;i++) comp[i]=-1;
    for(i=0;i<n;i++) if(comp[i]<0){ for(j=i;j<n;j++) if(comp[j]<0 && R[i][j] && R[j][i]) comp[j]=nc; nc++; }
    int C[8][8]; memset(C,0,sizeof(C));
    for(i=0;i<n;i++) for(j=0;j<n;j++) if(A[i][j] && comp[i]!=comp[j]) C[comp[i]][comp[j]]=1;
    if(rootrule){
        /* the root's component must be the unique source component, and the root its smallest vertex */
        int cr=comp[0];
        for(i=0;i<nc;i++) if(C[i][cr]) return 8;
        for(i=1;i<n;i++) if(comp[i]==cr && S[i]<S[0]) return 8;
    }
    if(nc==1) return 3;
    return check_module_A(C, nc);
}

static void found_module(int n){
    if(MODE==0){
        if(!check_module(n)) return;
        found++;
        int A[8][8]; uint8_t fl=0;
        fill(A, UADJ, Vsub, n);
        if(acyclic_A(A,n) && check_module_A(A,n)==0){ fl|=1; n_u++; }
        int img[8], distinct=1;
        for(int i=0;i<n;i++){ img[i]=scc_of[Vsub[i]]; for(int j=0;j<i;j++) if(img[j]==img[i]) distinct=0; }
        if(distinct){ fill(A, SADJ, img, n); if(acyclic_A(A,n) && check_module_A(A,n)==0){ fl|=2; n_s++;
            if(FI){ uint64_t ik=key_of(img,n); fwrite(&ik,8,1,FI); } } }
        int m=mtest(Vsub,n,0);
        if(m==0){ fl|=4; n_m++; } else if(m==3){ fl|=(3<<4); n_mfail3++; } else if(m==4){ fl|=(1<<4); n_mfail4++; }
        else { fprintf(stderr,"unexpected M code %d for a baseline module\n",m); exit(3); }
        if(FK){ uint64_t k=key_of(Vsub,n); fwrite(&k,8,1,FK); fwrite(&fl,1,1,FF); }
    } else if(MODE==1 || MODE==3){
        if(!check_module(n)) return;
        found++;
        if(FK){ uint64_t k=key_of(Vsub,n); fwrite(&k,8,1,FK); }
    } else {
        int m=mtest(Vsub,n,1);
        if(m==8){ n_rootskip++; return; }
        if(m!=0) return;
        found++;
        int A[8][8]; fill(A, CADJ, Vsub, n);
        if(acyclic_A(A,n) && check_module_A(A,n)==0) n_shared_base++;
        if(FK){ uint64_t k=key_of(Vsub,n); fwrite(&k,8,1,FK); }
    }
}

static void extend(int n, uint64_t *Vext, uint64_t *seen){
    nodes_visited++;
    if(n==K){ found_module(n); return; }
    uint64_t ext[W]; memcpy(ext,Vext,sizeof(uint64_t)*W);
    for(;;){
        int w=-1,k;
        for(k=0;k<W;k++) if(ext[k]){ w = k*64 + __builtin_ctzll(ext[k]); break; }
        if(w<0) break;
        ext[w>>6] &= ~(1ULL<<(w&63));
        if(q[n]<1.0 && rnd()>q[n]) continue;
        int i, ok=1;
        for(i=0;i<n;i++){
            adj[i][n] = (OUT[Vsub[i]][w>>6]>>(w&63))&1ULL;
            adj[n][i] = (OUT[w][Vsub[i]>>6]>>(Vsub[i]&63))&1ULL;
        }
        adj[n][n]=0;
        Vsub[n]=w;
        if(MODE!=2 && !acyclic(n+1)) ok=0;
        if(ok){
            uint64_t nx[W], ns[W];
            for(k=0;k<W;k++){
                uint64_t excl = OUT[w][k] & ~seen[k] & ~forbid[k];
                nx[k] = ext[k] | excl;
                ns[k] = seen[k] | OUT[w][k];
            }
            ns[w>>6] |= 1ULL<<(w&63);
            nx[w>>6] &= ~(1ULL<<(w&63));
            extend(n+1,nx,ns);
        }
    }
}

static FILE *xopen(const char *p){ FILE *f=fopen(p,"r"); if(!f){perror(p); exit(1);} return f; }

static int read_graph(const char *p, uint64_t (*G)[W], int *types){
    int n, m, c, u, v, cc; FILE *f=xopen(p);
    if(fscanf(f,"%d %d %d",&n,&m,&c)!=3){ fprintf(stderr,"bad header %s\n",p); exit(1); }
    for(int i=0;i<n;i++){ int t; if(fscanf(f,"%d",&t)!=1) exit(1); if(types) types[i]=t; }
    for(int i=0;i<m;i++){ if(fscanf(f,"%d %d %d",&u,&v,&cc)!=3) exit(1); G[u][v>>6]|=1ULL<<(v&63); }
    fclose(f); return n;
}

int main(int argc,char**argv){
    if(argc<16){ fprintf(stderr,"usage: %s MODE contracted uncontracted scc K nparts part seed q1..q6 keyprefix\n",argv[0]); return 1; }
    const char *md=argv[1];
    MODE = !strcmp(md,"base")?0 : !strcmp(md,"unc")?1 : !strcmp(md,"mod")?2 : !strcmp(md,"scc")?3 : -1;
    if(MODE<0){ fprintf(stderr,"unknown mode\n"); return 1; }
    int i, t1[MAXN], t2[MAXN];
    int n1=read_graph(argv[2], CADJ, t1), n2=read_graph(argv[3], UADJ, t2);
    if(n1!=n2){ fprintf(stderr,"node count mismatch\n"); return 1; }
    for(i=0;i<n1;i++) if(t1[i]!=t2[i]){ fprintf(stderr,"node type mismatch\n"); return 1; }
    { FILE *f=xopen(argv[4]); int nv, na, a, b;
      if(fscanf(f,"%d %d %d",&nv,&NC,&na)!=3 || nv!=n1){ fprintf(stderr,"bad scc file\n"); return 1; }
      for(i=0;i<nv;i++) if(fscanf(f,"%d",&scc_of[i])!=1) return 1;
      for(i=0;i<na;i++){ if(fscanf(f,"%d %d",&a,&b)!=2) return 1; SADJ[a][b>>6]|=1ULL<<(b&63); }
      fclose(f); }
    if(MODE==0) { NV=n1; memcpy(OUT,CADJ,sizeof(OUT)); }
    else if(MODE==3) { NV=NC; memcpy(OUT,SADJ,sizeof(OUT)); }
    else { NV=n1; memcpy(OUT,UADJ,sizeof(OUT)); }
    for(int u=0;u<NV;u++) for(int v=0;v<NV;v++) if(bit(OUT,u,v)) IN[v][u>>6]|=1ULL<<(u&63);

    K=atoi(argv[5]);
    int nparts=atoi(argv[6]), part=atoi(argv[7]);
    rngstate=(uint64_t)strtoull(argv[8],NULL,10);
    if(rngstate==0) rngstate=88172645463325252ULL;
    for(i=0;i<8;i++) q[i]=1.0;
    for(i=1;i<=6;i++) q[i]=atof(argv[8+i]);
    double prod=1.0; for(i=1;i<K;i++) prod*=q[i];
    const char *kp=argv[15];
    if(strcmp(kp,"-") && K<=6){
        char fn[2048];
        snprintf(fn,sizeof fn,"%s.keys",kp); FK=fopen(fn,"wb");
        if(MODE==0){ snprintf(fn,sizeof fn,"%s.flags",kp); FF=fopen(fn,"wb");
                     snprintf(fn,sizeof fn,"%s.images",kp); FI=fopen(fn,"wb"); }
        if(!FK || (MODE==0 && (!FF||!FI))){ perror("key files"); return 1; }
    }

    for(int s=part; s<NV; s+=nparts){
        if(MODE==2){                       /* only in-neighbours of the root with a smaller index are excluded */
            for(int k=0;k<W;k++) forbid[k]=0;
            for(int u=0;u<s;u++) if(bit(IN,s,u)) forbid[u>>6]|=1ULL<<(u&63);
        } else memcpy(forbid,IN[s],sizeof(uint64_t)*W);
        uint64_t ext[W],seen[W];
        int k;
        for(k=0;k<W;k++){ ext[k]=OUT[s][k]&~forbid[k]; seen[k]=OUT[s][k]; }
        ext[s>>6]&=~(1ULL<<(s&63));
        seen[s>>6]|=1ULL<<(s&63);
        Vsub[0]=s; adj[0][0]=0;
        extend(1,ext,seen);
    }
    printf("mode %s\nK %d\nnparts %d\npart %d\nseed %s\nq", md, K, nparts, part, argv[8]);
    for(i=1;i<K;i++) printf(" %g",q[i]);
    printf("\nprod_q %.17g\nfound %lld\nvisited %lld\nrng_final %llu\n", prod, found, nodes_visited,
           (unsigned long long)rngstate);
    if(MODE==0) printf("valid_U %lld\nretained_S %lld\nvalid_M %lld\nM_fail_D3 %lld\nM_fail_D4 %lld\n",
                       n_u, n_s, n_m, n_mfail3, n_mfail4);
    if(MODE==2) printf("shared_with_baseline %lld\nroot_rule_skipped %lld\n", n_shared_base, n_rootskip);
    if(FK) fclose(FK); if(FF) fclose(FF); if(FI) fclose(FI);
    return 0;
}
