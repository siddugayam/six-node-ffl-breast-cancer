/* census_09_esu_classmask.c -- independent n-node FFL census under the formal definition.
 *
 * D1 S acyclic after contracting each reciprocal TF<->miRNA pair to its TF->miRNA arc
 *    (the reduced graph G' is supplied on stdin already contracted)
 * D2 exactly one source (indeg 0 in S) and exactly one sink (outdeg 0 in S)
 * D3 every vertex lies on a directed source->sink path inside S
 * D4 >=2 internally vertex-disjoint source->sink paths
 *
 * Enumeration: source-rooted ESU.  For every candidate source s we grow vertex-induced
 * subgraphs using OUT-neighbourhoods only, with Wernicke's exclusive-neighbourhood rule
 * so that each out-reachable vertex set is generated exactly once.  A valid module has a
 * unique source and is out-reachable only from it, so each module is produced exactly once.
 * RAND-ESU: each extension candidate at depth d is kept with probability q[d];
 * unbiased estimate = found / prod(q).
 *
 * usage: census_09_esu_classmask graph.txt n q2 q3 q4 q5 q6 seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define W 16
static int NV, NE; static int lastmask=0;
static uint64_t OUT[MAXN][W], IN[MAXN][W];
static int ECLS[MAXN][MAXN];        /* too big; use hash instead */
static int NCLS;
static int nodetype[MAXN];

/* edge class lookup: sparse via per-node arrays */
static int *cls_idx[MAXN];          /* class of edge u->v, indexed by v, -1 if none */

static int K;                        /* subgraph size */
static double q[8];
static uint64_t rngstate;
static inline double rnd(void){
    rngstate ^= rngstate<<13; rngstate ^= rngstate>>7; rngstate ^= rngstate<<17;
    return (double)(rngstate>>11) / 9007199254740992.0;
}

static long long found=0;
static long long clshist[16]; static long long maskhist[512];
static int maxcls=0;
static long long nodes_visited=0;
static int examples_needed=0;
static int bestset[8]; static int bestn=0;

/* ---- small-graph tests ---- */
static int adj[8][8];

static int check_module(int *V,int n,int *pclasses){
    int i,j;
    int indeg[8],outdeg[8];
    for(i=0;i<n;i++){indeg[i]=0;outdeg[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){outdeg[i]++;indeg[j]++;}
    int src=-1,snk=-1,ns=0,nk=0;
    for(i=0;i<n;i++){ if(indeg[i]==0){src=i;ns++;} if(outdeg[i]==0){snk=i;nk++;} }
    if(ns!=1||nk!=1) return 0;
    if(src==snk) return 0;
    /* D3: reachable from src and reaches snk (acyclicity already enforced) */
    int reach[8],r2[8],ch,it;
    for(i=0;i<n;i++){reach[i]=0;r2[i]=0;}
    reach[src]=1; do{ch=0; for(i=0;i<n;i++) if(reach[i]) for(j=0;j<n;j++) if(adj[i][j]&&!reach[j]){reach[j]=1;ch=1;} }while(ch);
    r2[snk]=1; do{ch=0; for(i=0;i<n;i++) if(r2[i]) for(j=0;j<n;j++) if(adj[j][i]&&!r2[j]){r2[j]=1;ch=1;} }while(ch);
    for(i=0;i<n;i++) if(!reach[i]||!r2[i]) return 0;
    /* D4: >=2 internally vertex-disjoint src->snk paths, via unit-capacity vertex-split flow */
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
    /* distinct edge classes */
    {
        int seen[64]; memset(seen,0,sizeof(seen)); int nc=0; int mk=0;
        for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){
            int c=cls_idx[V[i]][V[j]];
            if(c>=0){ mk|=1<<c; if(!seen[c]){seen[c]=1;nc++;} }
        }
        *pclasses=nc; lastmask=mk;
    }
    return 1;
}

/* acyclicity of current induced subgraph on Vsub (size n) */
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

static int Vsub[8];
static uint64_t forbid[W];

static void extend(int n, uint64_t *Vext, uint64_t *seen){
    nodes_visited++;
    if(n==K){
        int nc=0;
        if(check_module(Vsub,n,&nc)){ found++; clshist[nc]++; maskhist[lastmask&511]++; if(nc>maxcls){maxcls=nc; bestn=n; memcpy(bestset,Vsub,sizeof(int)*n);} }
        return;
    }
    uint64_t ext[W]; memcpy(ext,Vext,sizeof(uint64_t)*W);
    int wi;
    for(;;){
        /* pop lowest set bit */
        int w=-1,k;
        for(k=0;k<W;k++) if(ext[k]){ w = k*64 + __builtin_ctzll(ext[k]); break; }
        if(w<0) break;
        ext[w>>6] &= ~(1ULL<<(w&63));
        if(q[n]<1.0 && rnd()>q[n]) continue;
        /* add w */
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
    (void)wi;
}

int main(int argc,char**argv){
    if(argc<4){fprintf(stderr,"usage: %s graph.txt K q2..q6 seed\n",argv[0]);return 1;}
    FILE*f=fopen(argv[1],"r");
    if(!f){perror("open");return 1;}
    if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3){fprintf(stderr,"bad header\n");return 1;}
    int i;
    for(i=0;i<NV;i++){ cls_idx[i]=malloc(sizeof(int)*NV); memset(cls_idx[i],0xff,sizeof(int)*NV); }
    for(i=0;i<NV;i++) fscanf(f,"%d",&nodetype[i]);
    for(i=0;i<NE;i++){
        int u,v,c; fscanf(f,"%d %d %d",&u,&v,&c);
        OUT[u][v>>6]|=1ULL<<(v&63);
        IN[v][u>>6]|=1ULL<<(u&63);
        cls_idx[u][v]=c;
    }
    fclose(f);
    K=atoi(argv[2]);
    for(i=0;i<8;i++) q[i]=1.0;
    int na=3;
    for(i=1;i<=6 && na<argc-1;i++){ q[i]=atof(argv[na++]); }
    rngstate = (argc>na)? (uint64_t)strtoull(argv[na],NULL,10) : 12345ULL;
    if(rngstate==0) rngstate=88172645463325252ULL;
    double prod=1.0; for(i=1;i<K;i++) prod*=q[i];

    for(int s=0;s<NV;s++){
        memcpy(forbid,IN[s],sizeof(uint64_t)*W);
        uint64_t ext[W],seen[W];
        int k;
        for(k=0;k<W;k++){ ext[k]=OUT[s][k]&~forbid[k]; seen[k]=OUT[s][k]; }
        ext[s>>6]&=~(1ULL<<(s&63));
        seen[s>>6]|=1ULL<<(s&63);
        Vsub[0]=s; adj[0][0]=0;
        extend(1,ext,seen);
    }
    printf("K=%d q=", K); for(i=1;i<K;i++) printf("%g,",q[i]);
    printf(" seed=%llu found=%lld est=%.1f visited=%lld maxclasses=%d\n",
           (unsigned long long)rngstate, found, (double)found/prod, nodes_visited, maxcls);
    printf("clshist:"); for(i=0;i<12;i++) if(clshist[i]) printf(" %d:%lld",i,clshist[i]);
    printf("\n");
    printf("maskhist:"); for(i=0;i<512;i++) if(maskhist[i]) printf(" %d:%lld",i,maskhist[i]); printf("\n");
    if(bestn){ printf("maxclass_example:"); for(i=0;i<bestn;i++) printf(" %d",bestset[i]); printf("\n"); }
    return 0;
}
