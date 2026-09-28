/* ffl_census_composition.c -- counting-only n-node FFL census for the first re-runs, Part B.
 *
 * Enumeration and module test are copied from scripts/v2/v2_census2.c (the program behind the
 * exhaustive counts in results/v2/verify_census.csv): source-rooted ESU over out-neighbourhoods,
 * acyclicity pruning, D1-D4 with a unit-capacity max-flow test for D4, and the same xorshift64
 * RAND-ESU with the same order of random draws.  With nparts=1 the 'found' and 'visited' counts
 * therefore reproduce logs/v2/census2_n*.log and logs/v2/census_n7_sampled*.log exactly.
 *
 * Nothing is listed.  For every module that passes D1-D4 the program only increments counters:
 *   composition   (#TF, #miRNA, #Gene)                         node types from graph_pub_fine.txt
 *   class mask    set of edge classes under the paper's 7-class _cls map (R4_graph_aug.txt codes:
 *                 0 TF->Gene, 1 TF->TF, 2 TF->miRNA, 3 miRNA->Gene, 4 miRNA->TF, 5 Gene->Gene,
 *                 6 miRNA-miRNA)
 *   flag bit 0 (a)      >=2 TF and >=2 miRNA and >=2 Gene nodes
 *   flag bit 1 (b_type) an arc TF-type->TF-type, an arc Gene-type->Gene-type, and a path
 *                       TF-type -> miRNA-type -> Gene-type, all inside the module
 *   flag bit 2 (b_cls)  the same, with TF->TF and Gene->Gene read as the paper's edge classes
 *                       (class 1 and class 5) and the path as class 2 followed by class 3
 *   flag bit 3 (c_lit)  TF1->TF2, TF2->t, TF2->m and TF1->m, with m->TF1 present in the
 *                       uncontracted network (orig_pub.txt), i.e. TF1<->m reciprocal; t is a
 *                       non-miRNA node distinct from TF1, TF2 and m
 *   flag bit 4 (c_full) c_lit plus TF1->t and m->t: the composite core {TF1, m, t} with TF2
 *                       regulating both the target and the miRNA (dyn_models.py tf_tf topology)
 * Up to three example modules per flag and per class mask with >= 6 classes are printed.
 *
 * usage: ffl_census_composition graph_pub_fine.txt R4_graph_aug.txt orig_pub.txt K nparts part seed
 *                               q1 q2 q3 q4 q5 q6
 *   q[d] is the probability of keeping a candidate when growing a subgraph of size d (v2 order).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define W 16
enum { T_MIR = 0, T_TF = 1, T_GENE = 2 };

static int NV, NE, NCLS;
static uint64_t OUT[MAXN][W], IN[MAXN][W];
static unsigned char PCLS[MAXN][MAXN];   /* paper class + 1, 0 = no arc */
static unsigned char ORIG[MAXN][MAXN];   /* arc in the uncontracted network */
static int nodetype[MAXN];
static char nodename[MAXN][64];

static int K;
static double q[8];
static uint64_t rngstate;
static inline double rnd(void){
    rngstate ^= rngstate<<13; rngstate ^= rngstate>>7; rngstate ^= rngstate<<17;
    return (double)(rngstate>>11) / 9007199254740992.0;
}

static long long found = 0, nodes_visited = 0;
static long long comp[512], maskhist[128], flaghist[32];
static int n_ex_flag[5], ex_flag[5][3][8], ex_role[5][3][4];
static int n_ex_mask[128], ex_mask[128][3][8];

static int adj[8][8];
static int Vsub[8];
static uint64_t forbid[W];

/* ---- D1-D4, verbatim logic of v2_census2.c check_module ---- */
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

/* composite TF2 architecture; full=1 also requires the composite core arcs TF1->t and m->t */
static int check_c(int n, int full, int *role){
    for(int a=0;a<n;a++){ if(nodetype[Vsub[a]]!=T_TF) continue;
      for(int b=0;b<n;b++){ if(b==a||nodetype[Vsub[b]]!=T_TF||!adj[a][b]) continue;
        for(int m=0;m<n;m++){ if(nodetype[Vsub[m]]!=T_MIR||!adj[a][m]||!adj[b][m]) continue;
          if(!ORIG[Vsub[m]][Vsub[a]]) continue;                 /* m -> TF1 in the original network */
          for(int t=0;t<n;t++){ if(t==a||t==b||t==m||nodetype[Vsub[t]]==T_MIR) continue;
            if(!adj[b][t]) continue;
            if(full && !(adj[a][t] && adj[m][t])) continue;
            role[0]=Vsub[a]; role[1]=Vsub[b]; role[2]=Vsub[m]; role[3]=Vsub[t];
            return 1; } } } }
    return 0;
}

static void tally(int n){
    int nt[3]={0,0,0}, i, j, k;
    for(i=0;i<n;i++) nt[nodetype[Vsub[i]]]++;
    comp[nt[T_TF]*64 + nt[T_MIR]*8 + nt[T_GENE]]++;
    int mask=0, tftf=0, gg=0, tftf_c=0, gg_c=0, path=0, path_c=0;
    for(i=0;i<n;i++)for(j=0;j<n;j++) if(adj[i][j]){
        int c = PCLS[Vsub[i]][Vsub[j]];
        if(!c){ fprintf(stderr,"arc without paper class %d->%d\n",Vsub[i],Vsub[j]); exit(2); }
        c--; mask |= 1<<c;
        int ti=nodetype[Vsub[i]], tj=nodetype[Vsub[j]];
        if(ti==T_TF && tj==T_TF) tftf=1;
        if(ti==T_GENE && tj==T_GENE) gg=1;
        if(c==1) tftf_c=1;
        if(c==5) gg_c=1;
        for(k=0;k<n;k++) if(adj[j][k]){
            if(ti==T_TF && tj==T_MIR && nodetype[Vsub[k]]==T_GENE) path=1;
            if(c==2 && PCLS[Vsub[j]][Vsub[k]]-1==3) path_c=1;
        }
    }
    maskhist[mask]++;
    int role_l[4], role_f[4];
    int fa = (nt[T_TF]>=2 && nt[T_MIR]>=2 && nt[T_GENE]>=2);
    int fb = tftf && gg && path;
    int fbc = tftf_c && gg_c && path_c;
    int fcl = check_c(n,0,role_l);
    int fcf = check_c(n,1,role_f);
    int fl = fa | (fb<<1) | (fbc<<2) | (fcl<<3) | (fcf<<4);
    flaghist[fl]++;
    for(int f=0; f<5; f++) if(((fl>>f)&1) && n_ex_flag[f]<3){
        int e=n_ex_flag[f]++;
        for(i=0;i<n;i++) ex_flag[f][e][i]=Vsub[i];
        for(i=0;i<4;i++) ex_role[f][e][i] = (f==3)? role_l[i] : (f==4)? role_f[i] : -1;
    }
    if(__builtin_popcount(mask)>=6 && n_ex_mask[mask]<3){
        int e=n_ex_mask[mask]++;
        for(i=0;i<n;i++) ex_mask[mask][e][i]=Vsub[i];
    }
}

static void extend(int n, uint64_t *Vext, uint64_t *seen){
    nodes_visited++;
    if(n==K){
        if(check_module(n)){ found++; tally(n); }
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

int main(int argc,char**argv){
    if(argc<14){fprintf(stderr,"usage: %s fine.txt R4.txt orig.txt K nparts part seed q1..q6\n",argv[0]);return 1;}
    int i, u, v, c;
    /* structure and node types: the graph used for the exhaustive counts */
    FILE *f=xopen(argv[1]);
    if(fscanf(f,"%d %d %d",&NV,&NE,&NCLS)!=3){fprintf(stderr,"bad header\n");return 1;}
    for(i=0;i<NV;i++) if(fscanf(f,"%d",&nodetype[i])!=1) return 1;
    for(i=0;i<NE;i++){ if(fscanf(f,"%d %d %d",&u,&v,&c)!=3) return 1;
        OUT[u][v>>6]|=1ULL<<(v&63); IN[v][u>>6]|=1ULL<<(u&63); }
    fclose(f);
    /* paper's 7-class codes and node names; must describe the same arc set */
    int n2, m2; f=xopen(argv[2]);
    if(fscanf(f,"%d %d\n",&n2,&m2)!=2 || n2!=NV || m2!=NE){fprintf(stderr,"R4 graph mismatch\n");return 1;}
    for(i=0;i<NV;i++){ char ty[16]; if(fscanf(f,"%63[^\t]\t%15[^\n]\n",nodename[i],ty)!=2) return 1;
        int t = !strcmp(ty,"miRNA")?T_MIR : !strcmp(ty,"TF")?T_TF : T_GENE;
        if(t!=nodetype[i]){fprintf(stderr,"node type mismatch at %d\n",i);return 1;} }
    for(i=0;i<NE;i++){ if(fscanf(f,"%d %d %d\n",&u,&v,&c)!=3) return 1;
        if(!((OUT[u][v>>6]>>(v&63))&1ULL)){fprintf(stderr,"R4 arc not in fine graph\n");return 1;}
        PCLS[u][v]=(unsigned char)(c+1); }
    fclose(f);
    /* uncontracted network, for reciprocity */
    int n3, m3, c3; f=xopen(argv[3]);
    if(fscanf(f,"%d %d %d",&n3,&m3,&c3)!=3 || n3!=NV){fprintf(stderr,"orig mismatch\n");return 1;}
    for(i=0;i<NV;i++){ int t; if(fscanf(f,"%d",&t)!=1) return 1; }
    for(i=0;i<m3;i++){ if(fscanf(f,"%d %d %d",&u,&v,&c)!=3) return 1; ORIG[u][v]=1; }
    fclose(f);

    K=atoi(argv[4]);
    int nparts=atoi(argv[5]), part=atoi(argv[6]);
    rngstate=(uint64_t)strtoull(argv[7],NULL,10);
    if(rngstate==0) rngstate=88172645463325252ULL;
    for(i=0;i<8;i++) q[i]=1.0;
    for(i=1;i<=6;i++) q[i]=atof(argv[7+i]);
    double prod=1.0; for(i=1;i<K;i++) prod*=q[i];

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
    printf("K %d\nnparts %d\npart %d\nseed %s\nq", K, nparts, part, argv[7]);
    for(i=1;i<K;i++) printf(" %g",q[i]);
    printf("\nprod_q %.17g\nfound %lld\nvisited %lld\nrng_final %llu\n", prod, found, nodes_visited,
           (unsigned long long)rngstate);
    for(i=0;i<512;i++) if(comp[i]) printf("comp %d %d %d %lld\n", i/64, (i/8)%8, i%8, comp[i]);
    for(i=0;i<128;i++) if(maskhist[i]) printf("mask %d %lld\n", i, maskhist[i]);
    for(i=0;i<32;i++) if(flaghist[i]) printf("flags %d %lld\n", i, flaghist[i]);
    for(int fl=0; fl<5; fl++) for(int e=0;e<n_ex_flag[fl];e++){
        printf("exflag %d", fl);
        for(i=0;i<K;i++) printf(" %s", nodename[ex_flag[fl][e][i]]);
        if(fl>=3){ printf(" | TF1=%s TF2=%s miRNA=%s target=%s", nodename[ex_role[fl][e][0]],
                   nodename[ex_role[fl][e][1]], nodename[ex_role[fl][e][2]], nodename[ex_role[fl][e][3]]); }
        printf("\n");
    }
    for(int mk=0; mk<128; mk++) for(int e=0;e<n_ex_mask[mk];e++){
        printf("exmask %d", mk);
        for(i=0;i<K;i++) printf(" %s", nodename[ex_mask[mk][e][i]]);
        printf("\n");
    }
    return 0;
}
