/* Exhaustive / sampled ESU census of n-node FFLs under D1-D4, replicating
   03_ffl_census.py exactly (greedy 2-path D4), plus a rigorous max-flow D4 for
   comparison, plus the class-diversity histogram under the paper's own _cls map. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 600
#define W 10                 /* 600/64 */
static int N, M;
static uint64_t out_b[MAXN][W], in_b[MAXN][W], und_b[MAXN][W];
static unsigned char ecls[MAXN][MAXN];   /* class+1, 0 = no edge */
static int K;
static long long n_sub=0, n_ffl_greedy=0, n_ffl_flow=0;
static long long classdiv[9];
static int best_classes=0;
static int best_mod[8];
static long long divergent=0;

static inline int bget(uint64_t *b,int i){return (b[i>>6]>>(i&63))&1ULL;}
static inline void bset(uint64_t *b,int i){b[i>>6]|=1ULL<<(i&63);}

/* ---- FFL predicate on a set of K vertices ---- */
static int sub[8];
static int adjS[8][8];       /* induced adjacency, local indices */

static int reach_from(int s,int *vis){
    int st[8],sp=0; st[sp++]=s; vis[s]=1; int cnt=1;
    while(sp){int u=st[--sp];for(int v=0;v<K;v++) if(adjS[u][v]&&!vis[v]){vis[v]=1;cnt++;st[sp++]=v;}}
    return cnt;
}
static int reach_to(int t,int *vis){
    int st[8],sp=0; st[sp++]=t; vis[t]=1; int cnt=1;
    while(sp){int u=st[--sp];for(int v=0;v<K;v++) if(adjS[v][u]&&!vis[v]){vis[v]=1;cnt++;st[sp++]=v;}}
    return cnt;
}
static int colour[8];
static int dfs_acyc(int u){
    colour[u]=1;
    for(int v=0;v<K;v++) if(adjS[u][v]){ if(colour[v]==1) return 0; if(colour[v]==0&&!dfs_acyc(v)) return 0; }
    colour[u]=2; return 1;
}
/* greedy 2-path, mirroring the python find_path (DFS with a stack) */
static int prevv[8];
static int find_path(int s,int t,int *banned,int *path,int *plen){
    int st[8],sp=0; for(int i=0;i<K;i++) prevv[i]=-2;
    prevv[s]=-1; st[sp++]=s;
    while(sp){
        int u=st[--sp];
        if(u==t) break;
        for(int v=0;v<K;v++){
            if(!adjS[u][v]) continue;
            if(prevv[v]!=-2) continue;
            if(banned[v]&&v!=t) continue;
            prevv[v]=u; st[sp++]=v;
        }
    }
    if(prevv[t]==-2) return 0;
    int L=0,u=t; while(u!=-1){path[L++]=u;u=prevv[u];}
    *plen=L; return 1;
}
/* rigorous: >=2 internally vertex-disjoint s->t paths via unit-node-capacity max flow */
static int maxflow2(int s,int t){
    /* node split: in=i, out=i+K ; cap 1 for internal nodes, big for s,t */
    int nn=2*K; static int cap[16][16];
    memset(cap,0,sizeof(cap));
    for(int i=0;i<K;i++) cap[i][i+K]=(i==s||i==t)?9:1;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) if(adjS[i][j]) cap[i+K][j]=9;
    int flow=0;
    while(flow<2){
        int par[16]; for(int i=0;i<nn;i++)par[i]=-2; par[s]=-1;
        int q[16],qh=0,qt=0; q[qt++]=s;
        while(qh<qt){int u=q[qh++];for(int v=0;v<nn;v++) if(cap[u][v]>0&&par[v]==-2){par[v]=u;q[qt++]=v;}}
        if(par[t+K]==-2) break;
        int v=t+K; while(v!=s){int u=par[v];cap[u][v]--;cap[v][u]++;v=u;}
        flow++;
    }
    return flow>=2;
}
static void handle(int *S){
    n_sub++;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) adjS[i][j]= (i!=j) && (ecls[S[i]][S[j]]!=0);
    int nsrc=0,nsnk=0,s=-1,t=-1;
    for(int i=0;i<K;i++){
        int hi=0,ho=0;
        for(int j=0;j<K;j++){ if(adjS[j][i])hi=1; if(adjS[i][j])ho=1; }
        if(!hi){nsrc++;s=i;} if(!ho){nsnk++;t=i;}
    }
    if(nsrc!=1||nsnk!=1) return;
    for(int i=0;i<K;i++)colour[i]=0;
    for(int i=0;i<K;i++) if(colour[i]==0&&!dfs_acyc(i)) return;
    int v1[8]={0},v2[8]={0};
    reach_from(s,v1); reach_to(t,v2);
    for(int i=0;i<K;i++) if(!(v1[i]&&v2[i])) return;
    int banned[8]={0},path[8],pl;
    if(!find_path(s,t,banned,path,&pl)) return;
    int g=0;
    for(int i=1;i<pl-1;i++) banned[path[i]]=1;
    int p2[8],pl2;
    if(find_path(s,t,banned,p2,&pl2)) g=1;
    int f=maxflow2(s,t);
    if(g) n_ffl_greedy++;
    if(f) n_ffl_flow++;
    if(g!=f) divergent++;
    if(!g) return;
    int seen[8]={0},nc=0;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++) if(adjS[i][j]){int c=ecls[S[i]][S[j]]-1; if(!seen[c]){seen[c]=1;nc++;}}
    classdiv[nc]++;
    if(nc>best_classes){best_classes=nc;for(int i=0;i<K;i++)best_mod[i]=S[i];}
}
/* ---- ESU ---- */
static int Sbuf[8];
static void extend(int depth,int v,uint64_t *ext){
    if(depth==K){ handle(Sbuf); return; }
    uint64_t e[W]; memcpy(e,ext,sizeof(e));
    while(1){
        int w=-1;
        for(int i=0;i<W;i++) if(e[i]){ w=i*64+__builtin_ctzll(e[i]); break; }
        if(w<0) break;
        e[w>>6]&=~(1ULL<<(w&63));
        uint64_t excl[W]; memset(excl,0,sizeof(excl));
        for(int i=0;i<depth;i++){ bset(excl,Sbuf[i]); for(int j=0;j<W;j++) excl[j]|=und_b[Sbuf[i]][j]; }
        uint64_t ne[W];
        for(int i=0;i<W;i++) ne[i]= e[i] | (und_b[w][i] & ~excl[i]);
        for(int i=0;i<v/64;i++) ne[i]=0;
        ne[v/64] &= ~((v%64)==63?~0ULL:((1ULL<<((v%64)+1))-1));
        Sbuf[depth]=w;
        extend(depth+1,v,ne);
    }
}
int main(int argc,char**argv){
    const char*gf=argv[1]; K=atoi(argv[2]);
    FILE*f=fopen(gf,"r"); if(!f){perror("graph");return 1;}
    fscanf(f,"%d %d\n",&N,&M);
    static char nm[MAXN][64],ty[MAXN][16];
    for(int i=0;i<N;i++) fscanf(f,"%63[^\t]\t%15[^\n]\n",nm[i],ty[i]);
    for(int i=0;i<M;i++){int a,b,c;fscanf(f,"%d %d %d\n",&a,&b,&c);
        ecls[a][b]=c+1; bset(out_b[a],b); bset(in_b[b],a); bset(und_b[a],b); bset(und_b[b],a);}
    fclose(f);
    for(int v=0;v<N;v++){
        uint64_t ext[W]; memset(ext,0,sizeof(ext));
        for(int i=0;i<W;i++) ext[i]=und_b[v][i];
        for(int i=0;i<v/64;i++) ext[i]=0;
        ext[v/64] &= ~((v%64)==63?~0ULL:((1ULL<<((v%64)+1))-1));
        Sbuf[0]=v; extend(1,v,ext);
    }
    printf("graph=%s K=%d nodes=%d edges=%d\n",gf,K,N,M);
    printf("connected_subgraphs=%lld  FFL_greedy(paper)=%lld  FFL_maxflow=%lld  divergent=%lld\n",
           n_sub,n_ffl_greedy,n_ffl_flow,divergent);
    printf("class_diversity_histogram (paper's is_ffl):\n");
    for(int i=1;i<=7;i++) if(classdiv[i]) printf("   %d classes -> %lld\n",i,classdiv[i]);
    printf("max_classes=%d  example_module:",best_classes);
    for(int i=0;i<K;i++) printf(" %s(%s)",nm[best_mod[i]],ty[best_mod[i]]);
    printf("\n");
    return 0;
}
