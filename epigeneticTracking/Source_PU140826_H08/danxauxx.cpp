
#include "danxincl.h"

// Global var
extern Body *bdp;
extern Exbd *exbd[NCORES];
unsigned int mskdhnp,mskdrvp,mskcndp,mskcolp;
unsigned int mskdhnn,mskdrvn,mskcndn,mskcoln;

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void set_fpu (unsigned int mode) opened
//asm ("fldcw %0" : : "m" (*&mode));
//control87(_PC_53, MCW_PC);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Intpower(int bb,int ee) opened
int ii,res;

res=1; for(ii=0;ii<ee;ii++) res*=bb;
return res;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Rnd1(int rnmax) opened
return (int)(rnmax*(rand()/(RAND_MAX+1.0))); 
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
double Sigmoide(int xx,double ss) opened
   double res=(double)(1/(1+(double)exp(-ss* ((double)xx/10000) ))-0.0); return res;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Mysortxxxx(float par[],int szar,int opt,int ord[]) opened
int   ii,si,xi,*flags; float af;
// opt=0: descending order,opt=1: ascending order

af=0; xi=0; // to avoid warning
flags=(int *)malloc(sizeof(int)*szar);
for(ii=0;ii<szar;ii++) flags[ii]=1;
for(ii=0;ii<szar;ii++) opened
   if(opt==0) opened af=-1; xi=0; closed // max
   if(opt==1) opened af=+2; xi=0; closed // min
   for(si=0;si<szar;si++) opened
      if(opt==0) if((af<par[si])&&(flags[si]==1)) opened 
         af=par[si]; xi=si; 
      closed
      if(opt==1) if((af>par[si])&&(flags[si]==1)) opened 
         af=par[si]; xi=si; 
      closed
   closed
   ord[ii]=xi; flags[xi]=0;
closed
free(flags);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Compare_dgx (const void * a, const void * b) opened
const Dgx *cga = (const Dgx *)a;
const Dgx *cgb = (const Dgx *)b;

if(cga->timer <cgb->timer) return -1;
if(cga->timer >cgb->timer) return +1;
if(cga->timer==cgb->timer) opened
if(cga->exord <cgb->exord) return -1;
if(cga->exord >cgb->exord) return +1;
if(cga->exord==cgb->exord) return +0;
closed

return -1;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Ftsortxdgx(int thid,int szar,int opt,Dgx *ard,Dgx *ars) opened
Exbd *edp=exbd[thid];

memcpy(&ard[0],&ars[0],sizeof(Dgx)*szar);
qsort (ard,szar,sizeof(Dgx),Compare_dgx);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Compare_guy (const void * aa, const void * bb) opened
const Guy *at = (const Guy *)aa;
const Guy *bt = (const Guy *)bb;

if(at->gaft >bt->gaft) return -1;
if(at->gaft <bt->gaft) return  1;
if(at->gaft==bt->gaft) return  0;

return -1;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Ftsortguy(int thid,int szar,Guy *gyar) opened

qsort(gyar,szar,sizeof(Guy),Compare_guy);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Compare_float (const void * aa, const void * bb) opened
const float *af = (const float *)aa;
const float *bf = (const float *)bb;

if(*af <*bf) return +1;
if(*af >*bf) return -1;
if(*af==*bf) return  0;

return -1;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Ftsortfloat(int szar,float *flar) opened

qsort(flar,szar,sizeof(float),Compare_float);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Leavexec(char *str) opened
printf(str); printf("\n"); exit(1); 
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Cvbase2int(bas basear[],int szar,int base,int *arind,int *ndecp) opened
int  ii;

// Least significant digit in last array position
*ndecp=0;
for(ii=0;ii<szar;ii++)
   // *ndecp+=(int)pow(base,ii)*basear[szar-1-ii];
   *ndecp+=(int)Intpower(base,ii)*(int)basear[szar-1-ii];
if(*arind!=-1) *arind+=szar;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Cvint2base(bas basear[],int szar,int base,int *arind,int *ndecp) opened
int   ii,ndec;

// Least significant digit in last array position
*arind=*arind; // To avoid warnings
ndec=*ndecp;
for(ii=0;ii<szar;ii++) opened
   basear[szar-1-ii]=(bas)(ndec%base);
   ndec=ndec/base;
closed
if(*arind!=-1) *arind+=szar;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Cvbase2flt(bas basear[],int szar,int base,int *arind,float *nfltp) opened
int  ndec; float maxval;

Cvbase2int(basear,szar,base,arind,&ndec);
maxval=(float)Intpower(base,szar);
*nfltp=((ndec/maxval)-(float)0.5)*2;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Varrescale(int *var,int scalea,int scaleb) opened
double af;

af=(double)*var;                //af:   [0,scalea-1]
af=af/(double)(scalea-1);       //af:   [0,1]
af=af*(scaleb-1);               //af:   [0,scaleb-1]
*var=(int)floor(af+0.5);        //af:   [0,scaleb-1] rounded

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Mocmatch(int *res,int *moc,int naposmin,int fscmax,int *napos,int *fsc) opened
int ii; // napos: nr of active positions

(*napos)=0; (*fsc)=0;
for(ii=0;ii<MOCLEN;ii++) opened
   if (res[ii]==-1) printf("osti ");
   if (res[ii]!=-1) (*napos)++;
   if((res[ii]!=-1)&&(res[ii]!=moc[ii])) (*fsc)++;
   if((*fsc)>fscmax) return NO;
closed
if((*napos>=naposmin)&&((*fsc)<=fscmax)) return YA; else return NO;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Moccopy(int mocd[],int mocs[]) opened
for(int hi=0;hi<MOCLEN;hi++) mocd[hi]=mocs[hi];
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
double Dopdst(Coord *aa,Coord *bb) opened

return //abs(aa->x-bb->x)+abs(aa->y-bb->y)+abs(aa->z-bb->z);
   sqrt(
      pow((double)(aa->x-bb->x),2)+
      pow((double)(aa->y-bb->y),2)+
      pow((double)(aa->z-bb->z),2));

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Setbitmask(void) opened
unsigned int nn,mm;

nn= 0; mm=16; mskdhnp = ~(~0 << mm) << nn;
nn=16; mm= 4; mskdrvp = ~(~0 << mm) << nn;
nn=20; mm= 4; mskcndp = ~(~0 << mm) << nn;
nn=24; mm= 8; mskcolp = ~(~0 << mm) << nn;

nn= 0; mm=16; mskdhnp = ~(~0 << mm) << nn; mskdhnn= ~mskdhnp;
nn=16; mm= 4; mskdrvp = ~(~0 << mm) << nn; mskdrvn= ~mskdrvp;
nn=20; mm= 4; mskcndp = ~(~0 << mm) << nn; mskcndn= ~mskcndp;
nn=24; mm= 8; mskcolp = ~(~0 << mm) << nn; mskcoln= ~mskcolp;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Lreader(int par,Clsr *clp) opened

if(par==Ldhn) return (int)clp->dhc;
if(par==Ldrv) return (int)clp->drv;
if(par==Lcnd) return (int)clp->cnd;
if(par==Lcol) return (int)clp->col;

return -1;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Xreader(int par,Clsr *clp,int *dhn,int *drv,int *cnd,int *col) opened

*dhn=(int)clp->dhc; *drv=(int)clp->drv;
*cnd=(int)clp->cnd; *col=(int)clp->col;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Lwriter(int par,Clsr *clp,int dhn,int drv,int cnd,int col) opened

if(dhn!=-1) clp->dhc =( int )dhn;
if(drv!=-1) clp->drv =(uchar)drv;
if(cnd!=-1) clp->cnd =(uchar)cnd;
if(col!=-1) clp->col =( int )col;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
double Max3(double d0,double d1,double d2) opened
double res;

res=(d0 >=d1) ? d0  : d1;
res=(res>=d2) ? res : d2;

return res;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Surfint(int as,Coord *cd,Exbd *edp) opened
int preres=NO; int xx,yy,zz;
int cr_000,crm100,crp100,crm010,crp010,crm001,crp001;

xx=cd->x; yy=cd->y; zz=cd->z;

// if "in interior grid" xy & z
if(!((xx>=1)&&(xx<GRIDX-1)&&(yy>=1)&&(yy<GRIDY-1)&&(zz>=1)&&(zz<GRIDZ-1)))
   return YA;
else if(edp->clarcnd[xx][yy][zz]==1) preres=YA;

if(preres==NO) return NO;

cr_000=0; crm100=0; crp100=0; crm010=0; crp010=0; crm001=0; crp001=0;

cr_000=Lreader(Ldhn,&edp->clar[xx-0][yy-0][zz-0]);
crm100=Lreader(Ldhn,&edp->clar[xx-1][yy-0][zz-0]);
crp100=Lreader(Ldhn,&edp->clar[xx+1][yy-0][zz-0]);
crm010=Lreader(Ldhn,&edp->clar[xx-0][yy-1][zz-0]);
crp010=Lreader(Ldhn,&edp->clar[xx-0][yy+1][zz-0]);
crm001=Lreader(Ldhn,&edp->clar[xx-0][yy-0][zz-1]);
crp001=Lreader(Ldhn,&edp->clar[xx-0][yy-0][zz+1]);

// if borderline xy & z (24-20 > 22-20) (4-0 > 2-0)
if(preres==YA) opened
   if((edp->clarcnd[xx-1][yy-0][zz-0]==0)||(crm100!=cr_000)) return YA;
   if((edp->clarcnd[xx+1][yy-0][zz-0]==0)||(crp100!=cr_000)) return YA;
   if((edp->clarcnd[xx-0][yy-1][zz-0]==0)||(crm010!=cr_000)) return YA;
   if((edp->clarcnd[xx-0][yy+1][zz-0]==0)||(crp010!=cr_000)) return YA;
   if((edp->clarcnd[xx-0][yy-0][zz-1]==0)||(crm001!=cr_000)) return YA;
   if((edp->clarcnd[xx-0][yy-0][zz+1]==0)||(crp001!=cr_000)) return YA;
closed

return NO;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Setedges(Coord *el,Coord *eh,Coord ed[]) opened

// coordinates : BACK
ed[0].x=el->x; ed[0].y=el->y; ed[0].z=el->z;//NW, B
ed[1].x=eh->x; ed[1].y=el->y; ed[1].z=el->z;//NE, B
ed[2].x=el->x; ed[2].y=eh->y; ed[2].z=el->z;//SW, B
ed[3].x=eh->x; ed[3].y=eh->y; ed[3].z=el->z;//SE, B
// coordinates : FRONT
ed[4].x=el->x; ed[4].y=el->y; ed[4].z=eh->z;//NW, F
ed[5].x=eh->x; ed[5].y=el->y; ed[5].z=eh->z;//NE, F
ed[6].x=el->x; ed[6].y=eh->y; ed[6].z=eh->z;//SW, F
ed[7].x=eh->x; ed[7].y=eh->y; ed[7].z=eh->z;//SE, F

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Getthid(Exbd **edpp) opened
int thid=0;

#if(MYOMP==1)
thid = omp_get_thread_num();
#endif
#if(MYMPI==1)
MPI_Comm_rank(MPI_COMM_WORLD,&thid);
#endif

*edpp=exbd[thid];

return thid;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Mpscpy(char *sndbuf, void *srcdat, int size, int *cnt) opened
memcpy(sndbuf,srcdat,size); *cnt+=size;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Mprcpy(char *recbuf, void *dstdat, int size, int *cnt) opened
memcpy(dstdat,recbuf,size); *cnt+=size;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Avgcd(Coord cdar[],int szar,Coord *avgcd) opened
int  ii,at;
Coord *crd,lcd;
double dista,distb;

// average position in field
lcd.x=0; lcd.y=0; lcd.z=0;
for(ii=0;ii<szar;ii++) opened
   crd=&cdar[ii]; lcd.x+=crd->x; lcd.y+=crd->y; lcd.z+=crd->z; closed
lcd.x/=szar; lcd.y/=szar; lcd.z/=szar;
dista=8888; at=-1;
for(ii=0;ii<szar;ii++) opened
   distb=Dopdst(&cdar[ii],&lcd);
   if(dista>distb) opened dista=distb; at=ii; closed
closed
memcpy(avgcd,&cdar[at],sizeof(Coord));

closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
void Printf3(FILE *fptr, char *fmt) opened

fprintf(stdout,"%s",fmt);
fprintf( fptr, "%s",fmt);

closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
void Extrdgt(int aa,Coord *cd) opened
int at;

at=aa; cd->x=at/100; at=at%100; 
cd->y=at/10; at=at%10; cd->z=at;

closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
int Insidegrid(int sx,int sy,int sz) opened

if ((sx<0)||(sx>=GRIDX)) return 0;
if ((sy<0)||(sy>=GRIDY)) return 0;
if ((sz<0)||(sz>=GRIDZ)) return 0;

return 1;
closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
float  Meansorted(int opt,float flar[],int szar,int xbad,float ordar[]) opened
int ci; float mean;

for(ci=0;ci<szar;ci++) ordar[ci]=flar[ci];
Ftsortfloat(szar,&ordar[0]);
mean=0; 
for(ci=0;ci<xbad;ci++) opened 
   if(opt==0) mean+=ordar[szar-1-ci]; if(opt==1) mean+=ordar[ci];
closed
mean/=xbad; 

return mean;
closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
float Scagl(float xf) opened
int ci; float yf,af,step;

/*yf=xf; if(xf>0.1) yf=(float)0.1;

if(xf<=0) yf=xf; else opened
   yf=0; step=(float)0.025;
   for(ci=0;ci<20;ci++) opened
      af=1-(float)5.00*ci*step; if(af<0) af=0;
      if(xf>(ci+1)*step) yf+=step*af; 
      else opened 
         yf+=(xf-ci*step)*af; break; 
      closed
   closed
closed*/

if(xf<=0) yf=xf; else 
   yf=xf-(float)2.5*pow(xf,2); // 2.5,5.0
if(xf>=0.2) yf=(float)0.10;    // 0.2,0.1 0.10,0.05

return yf;
closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
void Coordutil(int *vx,int *vy,int *vz,Coord *cd,int opt) opened

if(opt==0) opened *vx=cd->x; *vy=cd->y; *vz=cd->z; closed
if(opt==1) opened cd->x=*vx; cd->y=*vy; cd->z=*vz; closed

closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
void Bound(int *var,int min,int max) opened
   if(*var<min) *var=min;
   if(*var>max) *var=max;
closed

