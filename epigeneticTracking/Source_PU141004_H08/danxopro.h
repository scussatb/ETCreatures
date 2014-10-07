
// Fwd declarations
class Pop; class Cell; class Body;

typedef unsigned char uchar;
typedef   signed char schar;

typedef uchar bas;

typedef struct opened
   int x,y,z;
closed  Coord;

typedef struct opened
   uchar neugy; bas gen[XGENSZ];
   float gaft,shad,mfit; // gaft: GA fitness, shad: shape adherence
closed  Guy;

typedef struct opened
   int bval;
closed  Gnpbase;

typedef int  dimensa[GRIDY][GRIDZ];

typedef struct opened
   Coord  dpl[8]; double tr[9];
   int msa,msb,conr,c2nr,gas,col[9]; // gas: GRN activity status
closed  Dgo; // developmental gene's output

typedef struct opened
   int exord,swtch,timer,res[MOCLEN],dhnrx[3],dhptx,arpos,napos,fsc,exeas,dvfrd;
   int stepeval; Dgo dgo;
closed  Dgx; // developmental gene

typedef struct opened
   int actgen,copies; Dgx dgx;
closed  Dgf; // developmental gene frozen

typedef struct opened
   Coord cd; Dgo dgo; // dgo for grn
   int moc[MOCLEN],clr;
   int lcrstp,lcrisc,actstp,sigstp,dltcrt,actdrv; // lcrstp = last creation step
   int depend[50],ndep,infmoth[NMORPH]; // actstp is the activation step
closed  Dse; // driver set element

typedef struct opened
   int  gf,ge,sf,se;
closed  Frc;

typedef struct opened
   int  xf,xe,ls,us,dl;
closed  Frz;

typedef struct opened
   int dhc,col;
   unsigned char drv,cnd;
   unsigned mother;
closed  Clsr; // Cell struct

typedef struct opened
   dimensa *etis,*cdis,*cdtg,*sgrd;
closed  Envir;

typedef struct opened
   int   swc,ifp,thp,cdn,csn,act; // decoded csn: contribution sign
closed  Sgen;

typedef struct opened
   int  ord; float par;
closed  Parel;

typedef struct opened
   //     xxxxxxxx ----------------------------------------------------------
   int    dhnr,dhnrf,dhnr0,dhnro,stepeval,stepeval0;
   int    evtnra0[STAGES],evtnrb0[STAGES],klor[GRIDX][GRIDY][GRIDZ];
   int    moctmp[MOCLEN],mocm[MOCLEN],resi[MOCLEN],resk[MOCLEN];
   Clsr   clar [GRIDX][GRIDY][GRIDZ],clarf[GRIDX][GRIDY][GRIDZ];
   Clsr   clar0[GRIDX][GRIDY][GRIDZ],clard[DHARLS],moclp; Dgx dgar0[DGARSZ];
   Coord  drvar0[DHARLS],clardcd[DHARLS];
   Dse    dhar[DHARLS],dharf[DHARLS],dhar0[DHARLS],*dharo,dharaz[1][DHARLS];
   Envir  envr; Guy guyf,guyx; Pop *popa[NPOPUL];
   #if(MOCDEV==YA)
   int    dhnraltemp,dhnral[NNGUYS];
   Dse    dharal[NNGUYS][DHARLS],dharaltemp[DHARLS];
   #endif
   //     Dvarprep & chx sort -----------------------------------------------
   int    dgarxf,dgarsz,dgarszo,dgarsze,dgarszr;
   int    ordaux[DGARSZ],auxarp[DHARLS]; float par[DGARSZ];
   Dgx    dgar[DGARSZ],dgaro[DGARSZ],dgare[DGARSZ];
   //     Cvare -------------------------------------------------------------
   int    cgok,clok,bci[8],rci[8],uui[8];
   int    bufs[2*DPLMAX][2*DPLMAX][2*DPLMAX];
   int    bufc[2*DPLMAX][2*DPLMAX][2*DPLMAX];
   int    dltcrt,smocnew,mocnew[MOCLEN],dhor,drv,cnd,color;
   Dgo    actdgr; // active dg right part
   //     arrays of drivers -------------------------------------------------
   int    clarcnd[GRIDX][GRIDY][GRIDZ]; Coord drvaro[DHARLS],drvar1[DHARLS];
   //     movie -------------------------------------------------------------
   int    *divpar; Coord *divpts; dimensa *bstepsgrid,*estepsgrid,*xstepsgrid;
   //     initial values ----------------------------------------------------
   int    imoc[MOCLEN]; Clsr ics; Coord icd; Dgo idgo; Dgx idgx; Dse idse; Guy iguy;
   //     xxxxxxxx ----------------------------------------------------------
   int    ausdr,cn,cn0,gq,fset,fdone,pnr,xq,rndpmt[NNGUYS],xqar[DGARSZ],qf,qz;
   Coord  zygc[ZYGTOT]; Frc frc[STAGES]; Frz frz[STAGES];
   //     Calcfit
   float  shdhst[FITHST][NNGUYS],fithst[FITHST][NNGUYS];
   // Aux-morphogens  -------------------------------------------------------
   int    infmoth[NMORPH],nfdmoth,isval;
   // Grn -------------------------------------------------------------------
   int    sgarsz; Sgen sgar[DGARSZ];
   // --- Other -------------------------------------------------------------
   float  fshad,fmfit,fpopt;
   int    mothers[1000][2],motnr;
   char   astr[50]; int cnta,curpt,ctr; float exbst;
   Coord  ext[2]; FILE *cfp;
   // blocact = 0 (never evolved, not active); 1 (evolved and active); 2 (post-evolution, active)
closed  Exbd; // Ex body functions (global BUT PRIVATE: each process has a copy)

// AUX FCTS
void   set_fpu (unsigned int mode);
int    Intpower(int bb,int ee);
int    Rnd1(int rnmax);
double Sigmoide(int xx,double ss);
void   Mysortxxxx(float par[],int szar,int opt,int ord[]);
void   Ftsortxdgx(int thid,int szar,int opt,Dgx *ard,Dgx *ars);
void   Ftsortguy(int thid,int szar,Guy *gyar);
void   Ftsortfloat(int szar,float *flar);
void   Mysortlrxx(int thid,float par[],int szar,int opt,Clsr *ard[],Clsr *ars[]);
void   Leavexec(char *str);
void   Cvbase2int(bas basear[],int szar,int base,int *arind,int *ndecp);
void   Cvint2base(bas basear[],int szar,int base,int *arind,int *ndecp);
void   Cvbase2flt(bas basear[],int szar,int base,int *arind,float *nfltp);
void   Varrescale(int *var,int scale0,int scale1);
int    Mocmatch(int *res,int *moc,int naposmin,int fscmax,int *napos,int *fsc);
void   Moccopy(int mocd[],int mocs[]);
double Dopdst(Coord *aa,Coord *bb);
void   Setbitmask(void);
int    Lreader(int par,Clsr *clp);
void   Xreader(int par,Clsr *clp,int *dhn,int *drv,int *cnd,int *col);
void   Lwriter(int par,Clsr *clp,int  dhn,int  drv,int  cnd,int  col);
double Max3(double d0,double d1,double d2);
int    Surfint(int as,Coord *cd,Exbd *edp); // surface interior
void   Setedges(Coord *e0,Coord *e7,Coord e[]);
int    Getthid(Exbd **edpp);
void   Mpscpy(char *sndbuf,void *srcdat,int size,int *cnt);
void   Mprcpy(char *recbuf,void *dstdat,int size,int *cnt);
void   Avgcd(Coord cdar[],int szar,Coord *avgcd);
void   Printf3(FILE *fptr,char *fmt);
int    Insidegrid(int sx,int sy,int sz);
void   Extrdgt(int aa,Coord *cd);
float  Meansorted(int opt,float flar[],int szar,int xbad,float auxar[]);
float  Scagl(float fval);
void   Coordutil(int *vx,int *vy,int *vz,Coord *cd,int opt);
void   Bound(int *var,int min,int max);
void   writeVoxelyzeFile(std::vector < std::vector < std::vector < int > > > matrixForVoxelyze);
//!--------------------------------------------------------------------------
//! header
//!--------------------------------------------------------------------------
class Pop opened
/*VARS*/ public:
Body *popbdp;
int dnaxf,dnasz,dsold,dstmp,popsz,spopsz0,spopsz1,nkids0,nkids1;
int ipopsz,opopsz;
Guy guys[POPSZT],gold[POPSZT],*chmp;
Guy iguys[POPSZT],oguys[POPSZT]; Gnpbase dnp[XGENSZ];
/*FCTS*/ public:
void Guycpa (Guy gard[],Guy gars[],int ad,int as,int ncopied);
void Genrnda (int thid,int da,int db,int ngen);
uchar Genvalid (int thid,int da,int db,bas gen[]);
void Reproduct (int thid,int sons);
void Survbest (int thid,int ndeads, uchar deathtype);
void Crossover (int thid,float pcross);
void Mutation (int thid,float pcmutx,float maxpcmut);
void Galoop (int thid,int opt,FILE *fp);
void Gendecode (int thid,bas gen[],int opt);
Pop (int thid,Body * bdp);
//virtual void Genoocoder(int gen[]);
void Prepgen (int thid);
void Prexline (int thid);
void Germline (int thid);
virtual void Cons(int thid,Body *_popbdp);
closed;

//!--------------------------------------------------------------------------
//! header
//!--------------------------------------------------------------------------
class Body opened // all these vars are (should be) shared by all threads
/*Vars*/ public:

/*FCTS*/ public:
// FILE 1: high-level
Body(void);
int  main(int argc, char *argv[], char *env[]);
void Baseloop (int opt,Body *bdp);
void Calcbuflen (long *blp,long *blp0,long *blp1,long *blp3);
void Statevents (int thid,char *atmpstr,char *atmpstr1);
void Initlocv (int thid);
void Envinit (int thid,int opt);
void Bodinit (int thid);
int  Xqdet (int thid,int val);
// File 2: geometry-related
void Dgarprep (int thid,int gi,int as,int ss); // ss: start step
void Shaper (int thid,int gi,int as,int us);
void Movie (int thid,int gi,int as,int us,int cgevnr);
void Brush (int thid,int gi,int as,int di,Coord *cdv);
void Clset (int thid,int gi,int as,int di,Coord *cdv);
void Dharupd (int thid,int gi,int as,int di,Coord *cda,int clr);
void Daughter(int thid,int gi,int as,int di,Coord *cda,int clr);
// File 3
void Dopernew (int thid,int gi,int as,int is,int op);
void Doper2 (int thid,int gi,int as);
int  Closedr (int thid,int gi,Coord *cde,int ri,Coord *cd);
void Remred (int thid,int gi,int as,int di,Coord *cdv,int it);
void Loadtgt (int thid,char *file);
void Gridinit (int thid,int gi,int as);
void Gridupd (int thid,int gi,int as,int us,int ff);
double Fitness (int thid,int gi,Clsr clar[][GRIDY][GRIDZ]);
// File 4
void Sendrec0 (int nr,char *buffer,long *buflen0);
void Sendrecv (int nr,char *buffer,long *buflen1,long *buflen3);
void Calcfit (int thid);
void Showres (int thid,FILE *fpb,char *atmpstr);
void Loadmemo (int thid);
void Savememo (int thid);
void Savegrid (int thid,int opt,int hi);
void Copyclsr (int thid,Clsr *csd,Clsr *css);
void Moveclsr (int thid,Clsr *csd,Clsr *css);
// File 4: GRN
void Grncalc (int thid);
closed;


