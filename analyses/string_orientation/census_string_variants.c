/* census_string_variants.c -- census of the four conditions of the topological definition of a census graph, with STRING-related counts per module.
 *
 * Enumeration (source-rooted ESU over out-neighbourhoods), the RAND-ESU stream (xorshift64, same order of draws) and the
 * module test (acyclic() + check_module(); acyclicity, single source and sink, two disjoint paths; connectedness by construction) are copied from
 * analyses/census_and_motif_nulls/census/ffl_census_composition.c, so the published counts are reproduced.
 * For every module found it also counts:
 *   no_string_arc   the induced subgraph contains none of the arcs listed in ARCS (the STRING arcs of the graph)
 *   no_string_pair  no two of its vertices form one of the STRING associations listed in PAIRS (either orientation)
 *   in_baseline     its vertex set is in BASEKEYS (sorted uint64 keys of the baseline modules, same packing)
 *
 * usage: census_string_variants graph.txt K nparts part seed q1..q6 PAIRS ARCS BASEKEYS KEYPREFIX
 *   graph.txt: "NV NE 1" / NV node types / NE lines "u v c" (census fine format)
 *   PAIRS, ARCS: first line the count, then "a b" per line; '-' for none
 *   BASEKEYS: binary file of sorted uint64 keys or '-'; KEYPREFIX: write <prefix>.keys (K <= 6) or '-'
 *   key = sorted vertex indices packed 10 bits each
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define W 16

static int NV;
static uint64_t OUT[MAXN][W], IN[MAXN][W], PAIR[MAXN][W], SARC[MAXN][W];
static int K;
static double q[8];
static uint64_t rngstate;
static inline double rnd(void){
    rngstate ^= rngstate<<13; rngstate ^= rngstate>>7; rngstate ^= rngstate<<17;
    return (double)(rngstate>>11) / 9007199254740992.0;
}
static inline int bit(uint64_t (*G)[W], int u, int v){ return (G[u][v>>6]>>(v&63))&1ULL; }

static long long found = 0, nodes_visited = 0, n_noarc = 0, n_nopair = 0, n_inbase = 0;
static uint64_t *BK = NULL; static long long NBK = 0;
static FILE *FK = NULL;
static int adj[8][8];
static int Vsub[8];
static uint64_t forbid[W];

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

static uint64_t key_of(const int *S, int n){
    int s[8]; for(int i=0;i<n;i++) s[i]=S[i];
    for(int i=1;i<n;i++){ int x=s[i], j=i-1; while(j>=0 && s[j]>x){ s[j+1]=s[j]; j--; } s[j+1]=x; }
    uint64_t k=0; for(int i=0;i<n;i++) k=(k<<10)|(uint64_t)s[i];
    return k;
}

static int in_base(uint64_t k){
    long long lo=0, hi=NBK-1;
    while(lo<=hi){ long long m=(lo+hi)/2; if(BK[m]==k) return 1; if(BK[m]<k) lo=m+1; else hi=m-1; }
    return 0;
}

static void extend(int n, uint64_t *Vext, uint64_t *seen){
    nodes_visited++;
    if(n==K){
        if(check_module(n)){
            found++;
            int arc=0, pair=0;
            for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(i!=j){
                if(bit(SARC,Vsub[i],Vsub[j])) arc=1;
                if(bit(PAIR,Vsub[i],Vsub[j])) pair=1; }
            if(!arc) n_noarc++;
            if(!pair) n_nopair++;
            if(BK || FK){ uint64_t k=key_of(Vsub,n); if(BK && in_base(k)) n_inbase++; if(FK) fwrite(&k,8,1,FK); }
        }
        return;
    }
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
        if(!acyclic(n+1)) ok=0;
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

static void read_pairs(const char *p, uint64_t (*G)[W], int both){
    if(!strcmp(p,"-")) return;
    FILE *f=xopen(p); int m, a, b;
    if(fscanf(f,"%d",&m)!=1) exit(1);
    for(int i=0;i<m;i++){ if(fscanf(f,"%d %d",&a,&b)!=2) exit(1);
        G[a][b>>6]|=1ULL<<(b&63); if(both) G[b][a>>6]|=1ULL<<(a&63); }
    fclose(f);
}

int main(int argc,char**argv){
    if(argc<16){ fprintf(stderr,"usage: %s graph K nparts part seed q1..q6 PAIRS ARCS BASEKEYS KEYPREFIX\n",argv[0]); return 1; }
    int i, m, c, u, v, t;
    FILE *f=xopen(argv[1]);
    if(fscanf(f,"%d %d %d",&NV,&m,&c)!=3) return 1;
    for(i=0;i<NV;i++) if(fscanf(f,"%d",&t)!=1) return 1;
    for(i=0;i<m;i++){ if(fscanf(f,"%d %d %d",&u,&v,&c)!=3) return 1;
        OUT[u][v>>6]|=1ULL<<(v&63); IN[v][u>>6]|=1ULL<<(u&63); }
    fclose(f);
    K=atoi(argv[2]);
    int nparts=atoi(argv[3]), part=atoi(argv[4]);
    rngstate=(uint64_t)strtoull(argv[5],NULL,10);
    if(rngstate==0) rngstate=88172645463325252ULL;
    for(i=0;i<8;i++) q[i]=1.0;
    for(i=1;i<=6;i++) q[i]=atof(argv[5+i]);
    double prod=1.0; for(i=1;i<K;i++) prod*=q[i];
    read_pairs(argv[12], PAIR, 1);
    read_pairs(argv[13], SARC, 0);
    if(strcmp(argv[14],"-")){
        FILE *g=fopen(argv[14],"rb"); if(!g){ perror(argv[14]); return 1; }
        fseek(g,0,SEEK_END); NBK=ftell(g)/8; fseek(g,0,SEEK_SET);
        BK=malloc(NBK*8); if(fread(BK,8,NBK,g)!=(size_t)NBK){ fprintf(stderr,"short key file\n"); return 1; }
        fclose(g);
        for(long long x=1;x<NBK;x++) if(BK[x]<=BK[x-1]){ fprintf(stderr,"key file not sorted and unique\n"); return 1; }
    }
    if(strcmp(argv[15],"-") && K<=6){ char fn[2048]; snprintf(fn,sizeof fn,"%s.keys",argv[15]); FK=fopen(fn,"wb");
        if(!FK){ perror(fn); return 1; } }
    for(int s=part; s<NV; s+=nparts){
        memcpy(forbid,IN[s],sizeof(uint64_t)*W);
        uint64_t ext[W],seen[W];
        int k;
        for(k=0;k<W;k++){ ext[k]=OUT[s][k]&~forbid[k]; seen[k]=OUT[s][k]; }
        ext[s>>6]&=~(1ULL<<(s&63));
        seen[s>>6]|=1ULL<<(s&63);
        Vsub[0]=s; adj[0][0]=0;
        extend(1,ext,seen);
    }
    printf("graph %s\nK %d\nnparts %d\npart %d\nseed %s\nq", argv[1], K, nparts, part, argv[5]);
    for(i=1;i<K;i++) printf(" %g",q[i]);
    printf("\nprod_q %.17g\nfound %lld\nvisited %lld\nrng_final %llu\nno_string_arc %lld\nno_string_pair %lld\n",
           prod, found, nodes_visited, (unsigned long long)rngstate, n_noarc, n_nopair);
    if(BK) printf("in_baseline %lld\nbaseline_keys %lld\n", n_inbase, NBK);
    if(FK) fclose(FK);
    return 0;
}
