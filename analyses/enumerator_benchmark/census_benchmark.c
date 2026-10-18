/* census_benchmark.c -- the census enumerator of analyses/census_and_motif_nulls/census/ffl_census_composition.c
 * (source-rooted ESU over out-neighbourhoods, acyclicity pruning, D2-D4 test of check_module(); exhaustive), writing the
 * vertex set of every module found, with two relaxed tests:
 *   strict   the topological definition (paper numbering) as in the census
 *   d3relax  exactly one source and one sink still required; "every vertex on a source-sink path" dropped
 *   d4relax  one source-sink path suffices instead of two internally vertex-disjoint paths
 * usage: census_benchmark graph.txt K MODE keys_out
 *   graph.txt: "NV NE 1" / NV node types / NE lines "u v c" (census fine format); keys_out: binary uint64 keys
 *   (sorted vertex indices packed 10 bits each) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define W 16
static int NV, K, MODE;   /* 0 strict, 1 d3relax, 2 d4relax */
static uint64_t OUT[MAXN][W], IN[MAXN][W];
static long long found = 0, visited = 0;
static FILE *FK;
static int adj[8][8], Vsub[8];
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
    if(MODE!=1){
        int reach[8],r2[8],ch;
        for(i=0;i<n;i++){reach[i]=0;r2[i]=0;}
        reach[src]=1; do{ch=0; for(i=0;i<n;i++) if(reach[i]) for(j=0;j<n;j++) if(adj[i][j]&&!reach[j]){reach[j]=1;ch=1;} }while(ch);
        r2[snk]=1; do{ch=0; for(i=0;i<n;i++) if(r2[i]) for(j=0;j<n;j++) if(adj[j][i]&&!r2[j]){r2[j]=1;ch=1;} }while(ch);
        for(i=0;i<n;i++) if(!reach[i]||!r2[i]) return 0;
    }
    if(MODE==2) return 1;
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
    int indeg[8],i,j,removed=0,ch,gone[8];
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
static void extend(int n, uint64_t *Vext, uint64_t *seen){
    visited++;
    if(n==K){ if(check_module(n)){ found++; uint64_t k=key_of(Vsub,n); fwrite(&k,8,1,FK); } return; }
    uint64_t ext[W]; memcpy(ext,Vext,sizeof(uint64_t)*W);
    for(;;){
        int w=-1,k;
        for(k=0;k<W;k++) if(ext[k]){ w = k*64 + __builtin_ctzll(ext[k]); break; }
        if(w<0) break;
        ext[w>>6] &= ~(1ULL<<(w&63));
        for(int i=0;i<n;i++){
            adj[i][n] = (OUT[Vsub[i]][w>>6]>>(w&63))&1ULL;
            adj[n][i] = (OUT[w][Vsub[i]>>6]>>(Vsub[i]&63))&1ULL;
        }
        adj[n][n]=0; Vsub[n]=w;
        if(!acyclic(n+1)) continue;
        uint64_t nx[W], ns[W];
        for(k=0;k<W;k++){
            uint64_t excl = OUT[w][k] & ~seen[k] & ~forbid[k];
            nx[k] = ext[k] | excl; ns[k] = seen[k] | OUT[w][k];
        }
        ns[w>>6] |= 1ULL<<(w&63); nx[w>>6] &= ~(1ULL<<(w&63));
        extend(n+1,nx,ns);
    }
}
int main(int argc, char **argv){
    if(argc<5){ fprintf(stderr,"usage: %s graph K strict|d3relax|d4relax keys_out\n",argv[0]); return 1; }
    int m, c, u, v, t;
    FILE *f=fopen(argv[1],"r"); if(!f){ perror(argv[1]); return 1; }
    if(fscanf(f,"%d %d %d",&NV,&m,&c)!=3) return 1;
    for(int i=0;i<NV;i++) if(fscanf(f,"%d",&t)!=1) return 1;
    for(int i=0;i<m;i++){ if(fscanf(f,"%d %d %d",&u,&v,&c)!=3) return 1; OUT[u][v>>6]|=1ULL<<(v&63); IN[v][u>>6]|=1ULL<<(u&63); }
    fclose(f);
    K=atoi(argv[2]);
    MODE = !strcmp(argv[3],"strict")?0 : !strcmp(argv[3],"d3relax")?1 : !strcmp(argv[3],"d4relax")?2 : -1;
    if(MODE<0){ fprintf(stderr,"unknown mode\n"); return 1; }
    FK=fopen(argv[4],"wb"); if(!FK){ perror(argv[4]); return 1; }
    for(int s=0; s<NV; s++){
        memcpy(forbid,IN[s],sizeof(uint64_t)*W);
        uint64_t ext[W],seen[W];
        for(int k=0;k<W;k++){ ext[k]=OUT[s][k]&~forbid[k]; seen[k]=OUT[s][k]; }
        ext[s>>6]&=~(1ULL<<(s&63)); seen[s>>6]|=1ULL<<(s&63);
        Vsub[0]=s; adj[0][0]=0;
        extend(1,ext,seen);
    }
    fclose(FK);
    printf("K %d\nmode %s\nfound %lld\nvisited %lld\n", K, argv[3], found, visited);
    return 0;
}
