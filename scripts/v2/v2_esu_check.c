/* Independent cross-check: classical (undirected, min-label-root) ESU enumeration of ALL
   connected induced subgraphs of size K, each tested against D1-D4.  Completely different
   enumeration route from v2_census.c (which grows only along out-edges from the source). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 640
#define W 10
static int NV,NE,NCLS,K;
static uint64_t OUT[MAXN][W],UND[MAXN][W];
static int nodetype[MAXN];
static int *cls_idx[MAXN];
static long long found=0, visited=0; static int maxcls=0;
static long long clshist[16];
static int adj[8][8],Vsub[8];

static int test(int *V,int n,int*pnc){
    int i,j,indeg[8],outdeg[8];
    for(i=0;i<n;i++){indeg[i]=outdeg[i]=0;}
    for(i=0;i<n;i++)for(j=0;j<n;j++){ adj[i][j]= (i!=j) && ((OUT[V[i]][V[j]>>6]>>(V[j]&63))&1ULL); }
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){outdeg[i]++;indeg[j]++;}
    /* acyclic */
    { int gone[8],removed=0,ch,d[8];
      for(i=0;i<n;i++){gone[i]=0;d[i]=indeg[i];}
      do{ch=0; for(i=0;i<n;i++) if(!gone[i]&&d[i]==0){gone[i]=1;removed++;ch=1; for(j=0;j<n;j++) if(adj[i][j]) d[j]--;}}while(ch);
      if(removed!=n) return 0; }
    int src=-1,snk=-1,ns=0,nk=0;
    for(i=0;i<n;i++){ if(indeg[i]==0){src=i;ns++;} if(outdeg[i]==0){snk=i;nk++;} }
    if(ns!=1||nk!=1||src==snk) return 0;
    int reach[8],r2[8],ch;
    for(i=0;i<n;i++){reach[i]=0;r2[i]=0;}
    reach[src]=1; do{ch=0;for(i=0;i<n;i++) if(reach[i])for(j=0;j<n;j++) if(adj[i][j]&&!reach[j]){reach[j]=1;ch=1;}}while(ch);
    r2[snk]=1; do{ch=0;for(i=0;i<n;i++) if(r2[i])for(j=0;j<n;j++) if(adj[j][i]&&!r2[j]){r2[j]=1;ch=1;}}while(ch);
    for(i=0;i<n;i++) if(!reach[i]||!r2[i]) return 0;
    { int NN=2*n,cap[16][16],it; memset(cap,0,sizeof(cap));
      for(i=0;i<n;i++) cap[2*i][2*i+1]=(i==src||i==snk)?n:1;
      for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]) cap[2*i+1][2*j]=1;
      int S=2*src+1,T=2*snk,flow=0;
      for(it=0;it<3;it++){ int prev[16],q[16],h=0,t2=0,u,v;
        for(u=0;u<NN;u++)prev[u]=-1; prev[S]=S;q[t2++]=S;
        while(h<t2){u=q[h++];for(v=0;v<NN;v++) if(prev[v]<0&&cap[u][v]>0){prev[v]=u;q[t2++]=v;}}
        if(prev[T]<0)break; v=T; while(v!=S){u=prev[v];cap[u][v]--;cap[v][u]++;v=u;} flow++; if(flow>=2)break; }
      if(flow<2) return 0; }
    { int seen[64]; memset(seen,0,sizeof(seen)); int nc=0;
      for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){int c=cls_idx[V[i]][V[j]]; if(c>=0&&!seen[c]){seen[c]=1;nc++;}}
      *pnc=nc; }
    return 1;
}
static void ext(int n,uint64_t *Vext,uint64_t *seen,int root){
    visited++;
    if(n==K){ int nc=0; if(test(Vsub,n,&nc)){found++;clshist[nc]++; if(nc>maxcls)maxcls=nc;} return; }
    uint64_t e[W]; memcpy(e,Vext,sizeof(uint64_t)*W);
    for(;;){
        int w=-1,k;
        for(k=0;k<W;k++) if(e[k]){ w=k*64+__builtin_ctzll(e[k]); break; }
        if(w<0) break;
        e[w>>6]&=~(1ULL<<(w&63));
        uint64_t nx[W],ns[W];
        for(k=0;k<W;k++){ uint64_t excl=UND[w][k]&~seen[k];
            /* min-label root rule: only vertices > root */
            nx[k]=e[k]|excl; ns[k]=seen[k]|UND[w][k]; }
        for(k=0;k<=(root>>6);k++){ uint64_t m = (k<(root>>6))? ~0ULL : ((root&63)==63?~0ULL:((1ULL<<((root&63)+1))-1)); nx[k]&=~m; }
        ns[w>>6]|=1ULL<<(w&63); nx[w>>6]&=~(1ULL<<(w&63));
        Vsub[n]=w;
        ext(n+1,nx,ns,root);
    }
}
int main(int argc,char**argv){
    FILE*f=fopen(argv[1],"r"); if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3) return 1;
    int i; for(i=0;i<NV;i++){cls_idx[i]=malloc(sizeof(int)*NV);memset(cls_idx[i],0xff,sizeof(int)*NV);}
    for(i=0;i<NV;i++) if(fscanf(f,"%d",&nodetype[i])!=1) return 1;
    for(i=0;i<NE;i++){int u,v,c; if(fscanf(f,"%d %d %d",&u,&v,&c)!=3) return 1;
        OUT[u][v>>6]|=1ULL<<(v&63); UND[u][v>>6]|=1ULL<<(v&63); UND[v][u>>6]|=1ULL<<(u&63); cls_idx[u][v]=c;}
    fclose(f); K=atoi(argv[2]);
    for(int v=0;v<NV;v++){
        uint64_t e[W],s[W]; int k;
        for(k=0;k<W;k++){ e[k]=UND[v][k]; s[k]=UND[v][k]; }
        for(k=0;k<=(v>>6);k++){ uint64_t m=(k<(v>>6))?~0ULL:(((v&63)==63)?~0ULL:((1ULL<<((v&63)+1))-1)); e[k]&=~m; }
        s[v>>6]|=1ULL<<(v&63);
        Vsub[0]=v; ext(1,e,s,v);
    }
    printf("ESU-crosscheck K=%d found=%lld visited=%lld maxclasses=%d\nclshist:",K,found,visited,maxcls);
    for(i=0;i<12;i++) if(clshist[i]) printf(" %d:%lld",i,clshist[i]);
    printf("\n"); return 0;
}
