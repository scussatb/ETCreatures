
#include "danxincl.h"

// Global var
extern Exbd *exbd[NCORES]; extern char cvar;

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Guycpa(Guy gard[],Guy gars[],int ad,int as,int ncopied) opened

for(int gy=0;gy<ncopied;gy++) opened
   memcpy(&gard[ad+gy].gen[0],&gars[as+gy].gen[0],sizeof(bas)*dnasz);
   memcpy(&gard[ad+gy].gaft  ,&gars[as+gy].gaft  ,sizeof(float));
closed
if(gard==iguys) ipopsz=ncopied;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
uchar Pop::Genvalid(int thid,int dl,int du,bas gen[]) opened
uchar res; int di,dx,dy,dz,xl,xu;
Exbd *edp=exbd[thid];

xl=(dl-BEGNSQ)/DGENSQ;
xu=(du-BEGNSQ)/DGENSQ;

res=YA;
Gendecode(0,&gen[0],1);
for(di=xl;di<xu;di++) opened
   if(edp->dgar[di].exord!=-1) opened
      dx=(edp->dgar[di].dgo.dpl[7].x-edp->dgar[di].dgo.dpl[0].x+1);
      dy=(edp->dgar[di].dgo.dpl[0].y-edp->dgar[di].dgo.dpl[7].y+1);
      dz=(edp->dgar[di].dgo.dpl[7].z-edp->dgar[di].dgo.dpl[0].z+1);
      if(dx<=0) res=NO; // useful for "external" proliferation
      if(dy<=0) res=NO; // useful for "external" proliferation
      if(dz<=0) res=NO; // useful for "external" proliferation
   closed
closed
return res;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Genrnda(int thid,int da,int db,int nn) opened
int gy,bi,genok=0;

for(gy=0;gy<popsz;gy++) memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz);
opopsz=nn;

gy=0;
while(gy<nn) opened
   genok=0;
   for(bi=0;bi< dnasz;bi++) oguys[gy].gen[bi]=(bas)Rnd1(dnp[bi].bval);
   // header always initialized to 0
   for(bi=0;bi<GHDRSZ;bi++) oguys[gy].gen[bi]=(bas)0; 
   if(Genvalid(thid,da,db,oguys[gy].gen)==YA) genok=1;
   else memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz);
   if (genok==1) gy++;
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Reproduct(int thid,int sons) opened
int gy,ii,rn0,son=0; float gaftsum,rn,*wheel;

// cmt: gaft must be >=0
wheel=(float *)malloc(sizeof(float)*popsz);
for(gy=0;gy<popsz;gy++) memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz);
opopsz=sons;

gaftsum=0; for(gy=0;gy<ipopsz;gy++) gaftsum+=iguys[gy].gaft;
if(gaftsum>0) opened
   wheel[0]=iguys[0].gaft / gaftsum;
   for(gy=0;gy<ipopsz-1;gy++) wheel[gy+1]=wheel[gy]+iguys[gy+1].gaft/gaftsum;
   for(ii=0;ii<sons;ii++) opened
      rn=Rnd1(1000)/(float)1000;
      if(rn < wheel[0]) son=0;
      for(gy=0;gy<ipopsz-1;gy++) if(rn>=wheel[gy] && rn<wheel[gy+1]) opened
         son=gy+1; break; 
      closed
      memcpy(&oguys[ii].gen[0],&iguys[son].gen[0],sizeof(bas)*dnasz);
   closed
closed
else
   for(ii=0;ii<sons;ii++) opened
      rn0=Rnd1(ipopsz);
      memcpy(&oguys[ii].gen[0],&iguys[rn0].gen[0],sizeof(bas)*dnasz);
   closed

free(wheel);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Survbest(int thid,int ndeads,unsigned char deathtype) opened
int gy,gb;

for(gy=0;gy<popsz;gy++) memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz); 
opopsz=ipopsz-ndeads;

Ftsortguy(thid,ipopsz,iguys);
for(gy=0;gy<ipopsz-ndeads;gy++) opened
   if(deathtype==0) gb=gy;
   if(deathtype==1) gb=ipopsz-1-gy;
   Guycpa(oguys,iguys,gy,gb,1);
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Crossover(int thid,float pcross) opened
uchar crossok=NO; bas *crossdna0,*crossdna1;
int *crosslist,gy,ii,nn,er,es,za,zb,at,bt,iter; float fu;
Exbd *edp=exbd[thid];

crossdna0=(bas *)malloc(sizeof(bas)*dnasz);
crossdna1=(bas *)malloc(sizeof(bas)*dnasz);
crosslist=(int *)malloc(sizeof(int)*popsz);

for(gy=0;gy<popsz;gy++) memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz);
opopsz=ipopsz;

for(gy=0;gy<ipopsz;gy++) crosslist[gy]=gy;

for(nn=ipopsz;nn>0;nn-=2) opened
   ii=ipopsz-nn;
   er=-1; es=-1; 
   while(er>=es) opened er=Rnd1(nn); es=Rnd1(nn); closed
   fu=Rnd1(1000)/(float)1000;
   // Crossover attempt
   crossok=NO; iter=0;
   if(fu<pcross) opened
      while((crossok==NO)&&(iter<=ATPMAX)) opened
         iter++;
         za=Rnd1(dnasz-0); // 1. crossover point
         zb=Rnd1(dnasz-0); // 2. crossover point
         if (za<=zb) opened at=za; bt=zb; closed
         if (za> zb) opened at=zb; bt=za; closed
         if((at>=dstmp)&&(bt>=dstmp)) opened
            memcpy(&crossdna0[ 0],&iguys[crosslist[er]].gen[ 0],sizeof(bas)*dnasz);
            memcpy(&crossdna1[ 0],&iguys[crosslist[es]].gen[ 0],sizeof(bas)*dnasz);
            memcpy(&crossdna0[at],&iguys[crosslist[es]].gen[at],sizeof(bas)*(bt-at+1));
            memcpy(&crossdna1[at],&iguys[crosslist[er]].gen[at],sizeof(bas)*(bt-at+1));
            if(Genvalid(thid,0,dnasz,&crossdna0[0])==YA)
            if(Genvalid(thid,0,dnasz,&crossdna1[0])==YA) opened
               memcpy(&oguys[ii+0].gen[0],&crossdna0[0],sizeof(bas)*dnasz);
               memcpy(&oguys[ii+1].gen[0],&crossdna1[0],sizeof(bas)*dnasz); 
               crossok=YA;
            closed
         closed
      closed
   closed
   // Crossover failure
   if((fu>=pcross)||(crossok==NO)) opened
      Guycpa(oguys,iguys,ii+0,crosslist[er],1);
      Guycpa(oguys,iguys,ii+1,crosslist[es],1);
   closed
   // crossar rearrangement
   for(gy=er-0;gy<=nn-1;gy++) crosslist[gy]=crosslist[gy+1];
   for(gy=es-1;gy<=nn-1;gy++) crosslist[gy]=crosslist[gy+1];
closed

free(crossdna0);
free(crossdna1);
free(crosslist);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Mutation(int thid,float pcmutx,float maxpcmut) opened
uchar dnaok; bas *neudna;
int gy,ii,nn,zp,iter,mxchgd; float uf;
Exbd *edp=exbd[thid];

neudna=(bas *)malloc(sizeof(bas)*dnasz);
mxchgd=(int)(maxpcmut*dnasz);

for(gy=0;gy<popsz;gy++) memset(&oguys[gy].gen[0],0,sizeof(bas)*dnasz);
opopsz=ipopsz;

for(gy=0;gy<ipopsz;gy++) opened
   uf=Rnd1(1000)/(float)1000;
   dnaok=NO; iter=0;
   if(uf<pcmutx) opened
      while((dnaok==NO)&&(iter<=ATPMAX)) opened
         if(0==0) opened
            memcpy(&neudna[0],&iguys[gy].gen[0],sizeof(bas)*dnasz);
            nn=Rnd1(mxchgd);
         closed
         for(ii=0;ii<nn;ii++) opened
            zp=Rnd1(dnasz);
            if(zp<GHDRSZ) zp=GHDRSZ; // no mutations in the header
            if(zp>=dstmp) 
            if(zp>=dnaxf) if(zp<dnasz) neudna[zp]=Rnd1(dnp[zp].bval);
         closed
         if(Genvalid(thid,0,dnasz,&neudna[0])==YA) opened
            memcpy(&oguys[gy].gen[0],&neudna[0],sizeof(bas)*dnasz);
            dnaok=YA;
         closed
         else iter++;
      closed
   closed
   if((uf>=pcmutx)||(dnaok==NO)) opened 
      memcpy(&oguys[gy].gen[0],&iguys[gy].gen[0],sizeof(bas)*dnasz); 
   closed
closed

free(neudna);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Galoop(int thid,int opt,FILE *fp) opened
char atmpstr[50]; int gy; float step;
Exbd *edp=exbd[thid];

#if(DEBUG0==YA) // non-monotony check
if(edp->cn>edp->cn0) if(edp->exbst>guys[0].shad) opened
   sprintf(atmpstr," exbst %d %f %f",edp->cn,edp->exbst,guys[0].shad); 
   Printf3(fp,atmpstr); 
   //Leavexec(" exbst error");
closed
edp->exbst=guys[0].shad;
#endif

// Inits iguys,oguys,guysold
for(gy=0;gy<popsz;gy++) memcpy(&iguys[gy],&edp->iguy,sizeof(Guy));
for(gy=0;gy<popsz;gy++) memcpy(&oguys[gy],&edp->iguy,sizeof(Guy));
for(gy=0;gy<popsz;gy++) memcpy(&gold [gy],&guys[gy], sizeof(Guy));
for(gy=0;gy<popsz;gy++) memcpy(&guys [gy],&edp->iguy,sizeof(Guy));

// nkids0 gen's are sons of spopsz0
Guycpa(iguys,gold,0,0,spopsz0); Reproduct(thid,nkids0);
Guycpa(iguys,oguys,0,0,nkids0); Crossover(thid,PCROSS);
Guycpa(iguys,oguys,0,0,nkids0); Mutation(thid,(float)PMUTAT,(float)MPCMUT);
Guycpa( guys,oguys,0,0,nkids0);

// nkids1 gen's are sons of spopsz1
Guycpa(iguys,gold,0,spopsz0,spopsz1); Reproduct(thid,nkids1);
Guycpa( guys,oguys,nkids0,0,nkids1); // once gaft not copied

// Replaces worst nkids0+nkids1 of spopsz0
Guycpa(iguys,gold,0,0,spopsz0); Survbest(thid,nkids0+nkids1,0);
Guycpa( guys,oguys,nkids0+nkids1,0,spopsz0-nkids0-nkids1); // once gaft not copied

// New spopsz1
Genrnda(thid,0,dnasz,spopsz1);
Guycpa(guys, oguys,spopsz0,0,spopsz1);

// Old chmp in pos 0
memcpy(&guys[0],&gold[0],sizeof(Guy));

// Init guys - Checks validity
for(gy=0;gy<popsz;gy++) opened
   if(Genvalid(thid,0,dnasz,guys[gy].gen)==NO) opened
      Genrnda (thid,0,dnasz,1);
      memcpy(&guys[gy].gen[0],&oguys[0].gen[0],sizeof(bas)*dnasz);
   closed
   guys[gy].neugy=YA;
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
Pop::Pop (int thid,Body * bdp) opened
Cons(thid,bdp);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Pop::Cons(int thid,Body * _popbdp) opened
int  gy,bi,at;

popbdp=_popbdp;  
spopsz0=POPSZ0; nkids0=NKIDS0; popsz=POPSZT; 
spopsz1=POPSZ1; nkids1=NKIDS1; dnasz=XGENSZ;

// Init dnp
for(bi=0;bi<XGENSZ;bi++) opened
   at=4; if(bi==0) at=DGARSZ; 
   if(bi==1) at=-1; if(bi==2) at=DPLMAX;
   dnp[bi].bval=at;
closed

if(thid==0) opened
   Genrnda(thid,0,dnasz,popsz);
   for(gy=0;gy<popsz;gy++) opened
      memcpy(&guys[gy].gen[0],&oguys[gy].gen[0],sizeof(bas)*dnasz);
      guys[gy].neugy=YA; guys[gy].gaft=0; guys[gy].shad=0;
   closed
   //Gendecoder(&guys[0].gen[0],&guys[0]);
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Pop::Gendecode(int thid,bas gen[],int opt) opened
int   dx,dy,dz,ii,hi,di,bi,at,bt,ct,dl,max,min,xf,*cds,maxdgd;
float coe; Dgx curdgx; Coord nwb,sef; Exbd *edp=exbd[thid];

bi=3; //first 3 numbers in the Genome are not decoded
// stepeval
Cvbase2int(&gen[bi],XTMRSQ,4,&bi,&edp->stepeval);
if(edp->stepeval<0) edp->stepeval=0;
if(edp->stepeval>STAGES-1) edp->stepeval=STAGES-1;

// Inits dgar
if (opt==0) xf=0; if (opt==1) xf=edp->dgarxf;
for(di=0;di<xf;di++)
   memcpy(&edp->dgar[di],&edp->dgar0[di],sizeof(Dgx));

for(di=xf;di<edp->dgarsz;di++) opened
   #if(GRNMODE==YA)
   Cvbase2int(&gen[bi],&edp->sgar[di].swc,SGSNSQ,4,&bi);
   for(ii=0;ii<SGIASZ;ii++) 
   Cvbase2int(&gen[bi],&edp->sgar[di].ifp,SGIFSQ,4,&bi);
   Cvbase2int(&gen[bi],&edp->sgar[di].thp,SGTHSQ,4,&bi);
   Cvbase2int(&gen[bi],&edp->sgar[di].cdn,SGCDSQ,4,&bi);
   Cvbase2int(&gen[bi],&edp->sgar[di].csn,SGSCSQ,4,&bi);
   if(edp->sgar[di].swc<=2) edp->sgar[di].swc =0; else edp->sgar[di].swc=1;
   if(edp->sgar[di].csn<=2) edp->sgar[di].csn=-1; else edp->sgar[di].csn=1;
   Varrescale(&edp->sgar[di].cdn,(int)Intpower(4,SGSCSQ),10);
   goto GTLEOC;
   #endif
   // inits curdgx
   memcpy(&curdgx,&edp->idgx,sizeof(Dgx));
   bi=BEGNSQ+di*DGENSQ;
   if((int)gen[bi++]<=1) curdgx.swtch=0; else curdgx.swtch=1;
   Cvbase2int(&gen[bi],XORDSQ,4,&bi,&curdgx.exord);
   #if(PGFREEZE==NO)
   Varrescale(&curdgx.exord,(int)pow(4,(float)XORDSQ),edp->dgarsz);
   #endif
   Cvbase2int(&gen[bi],XTMRSQ,4,&bi,&curdgx.timer);
   Varrescale(&curdgx.timer,(int)Intpower(4,XTMRSQ),STAGES);
   Cvbase2int(&gen[bi],RTMRSQ,4,&bi,&bt); bt=3; // FORCE timer relevance
   if(bt<=1) curdgx.timer=-1; // timer not relevant
   if(TMON==NO) curdgx.timer=-1;
   // res
   for(hi=0;hi<MOCLEN;hi++) opened
      Cvbase2int(&gen[bi],XRESSQ,4,&bi,&curdgx.res[hi]); 
   closed
   for(hi=0;hi<MOCLEN;hi++) opened
      Cvbase2int(&gen[bi],1,4,&bi,&bt); // ignored
   closed
   // dpl
   for(ii=0;ii<6;ii++) opened
      if(ii==0) opened cds=&nwb.x; at=-1; closed // - negative
      if(ii==1) opened cds=&nwb.y; at=+1; closed // + positive
      if(ii==2) opened cds=&nwb.z; at=-1; closed // - negative
      if(ii==3) opened cds=&sef.x; at=+1; closed // + positive
      if(ii==4) opened cds=&sef.y; at=-1; closed // - negative
      if(ii==5) opened cds=&sef.z; at=+1; closed // + positive
      Cvbase2int(&gen[bi],XDPLSQ,4,&bi,&bt);
      dl=edp->frz[edp->xqar[di]].dl;
      if(bt>dl) bt=dl; *cds=(at*bt)-at;
      // -at allows for "external proliferation"
      // maxval=0 -> maxval=+1, minval=0 -> minval=-1
   closed
   // msa,msb,shp,dgd
   Cvbase2int(&gen[bi],XMSASQ,4,&bi,&curdgx.dgo.msa);
   Cvbase2int(&gen[bi],XMSBSQ,4,&bi,&curdgx.dgo.msb);
   maxdgd=Intpower(4,XDGDSQ)-1;
   for(ii=0;ii<9;ii++) opened
      Cvbase2int(&gen[bi],XDGDSQ,4,&bi,&at);
      curdgx.dgo.tr[ii]=(((double)at/maxdgd)-0.5)*2; 
   closed
   if(NDIMS==2) for(ii=4;ii<8;ii++) curdgx.dgo.tr[ii]=0;
   // colours
   for(ii=0;ii<9;ii++) 
      Cvbase2int(&gen[bi],XCOLSQ,4,&bi,&curdgx.dgo.col[ii]);
   curdgx.dgo.conr=curdgx.dgo.col[0];
   curdgx.dgo.c2nr=curdgx.dgo.col[1];
   // pmov 
   Cvbase2int(&gen[bi],APOSSQ,4,&bi,&bt); curdgx.dgo.dpos.x=bt;
   Cvbase2int(&gen[bi],APOSSQ,4,&bi,&bt); curdgx.dgo.dpos.y=bt;
   Cvbase2int(&gen[bi],APOSSQ,4,&bi,&bt); curdgx.dgo.dpos.z=bt; 
   Varrescale( &curdgx.dgo.dpos.x,Intpower(4,APOSSQ),GRIDX);
   Varrescale( &curdgx.dgo.dpos.y,Intpower(4,APOSSQ),GRIDY);
   Varrescale( &curdgx.dgo.dpos.z,Intpower(4,APOSSQ),GRIDZ);
   if(NDIMS==2) curdgx.dgo.dpos.z=0;   
   // end of decoding -------------------------------------------------------
   // on-the-fly corrections and forced -------------------------------------
   // dx,dy,dz
   dx=(sef.x-nwb.x+1); // ACHTUNG
   dy=(nwb.y-sef.y+1); // ACHTUNG
   dz=(sef.z-nwb.z+1); // ACHTUNG
   if(dx<=0) opened nwb.x-=1; sef.x+=1; closed
   if(dy<=0) opened sef.y-=1; nwb.y+=1; closed
   if(dz<=0) opened nwb.z-=1; sef.z+=1; closed

   #if  (SPHERE==1)
   max=(dx >=dy) ? dx : dy; max=(max>=dz) ? max : dz;
   min=(dx <=dy) ? dx : dy; min=(min<=dz) ? min : dz;
   coe=(min/(float)max);
   if(coe<0.5) opened nwb.y=sef.x; sef.y=nwb.x; closed
   if(coe<0.5) opened nwb.z=nwb.x; sef.z=sef.x; closed
   #elif(SPHERE==2)
   nwb.x=-sef.x; sef.x=+sef.x;
   nwb.y=+sef.x; sef.y=-sef.x;
   nwb.z=-sef.x; sef.z=+sef.x;
   #endif
   Setedges(&nwb,&sef,&curdgx.dgo.dpl[0]);
   // end -------------------------------------------------------------------
   memcpy(&edp->dgar[di],&curdgx,sizeof(Dgx));
closed
// Last inits
for(di=0;di<edp->dgarsz;di++) edp->dgar[di].arpos=di;

// Sorts by exord ?
cvar=cvar; //anchor
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Pop::Prepgen(int thid) opened
int   pi,gy,at,bt,ct,dt,st,nt;
Guy  *gyp,*gzp,*gtp; Exbd *edp=exbd[thid];

if(edp->gq>0) ct=edp->frc[edp->gq-1].ge; else ct=0;
if(edp->gq>0) st=edp->frc[edp->gq-1].se; else st=0;
if((edp->cn-ct>=0)&&(edp->cn-ct<COLDGENS)) 
   dstmp=BEGNSQ+(DGENSQ*st);
else dstmp=-1;
nt=dstmp;

#if(PGFREEZE==YA)
for(pi=0;pi<NPOPUL;pi++) for(gy=0;gy<edp->popa[pi]->popsz;gy++) opened
   at=  BEGNSQ;               // lower boundary of genome
   bt=  edp->popa[ 0]->dnaxf; // upper boundary of sector frozen
   dt=  edp->popa[ 0]->dnasz; // upper boundary of sector evolved 
   ct=  edp->popa[ 0]->dsold; // upper boundary of sector evolved old
   gyp=&edp->popa[pi]->guys [gy];
   gzp=&edp->popa[ 0]->guys [ 0];
   gtp=&edp->popa[pi]->oguys[gy];
   // copies frozen sector into all genomes
   memcpy(&gyp->gen[at],&gzp->gen[at],sizeof(bas)*(bt-at));
 //memcpy(&gyp->gen[ 0],&gzp->gen[nt],sizeof(int)*(nt- 0));
   if((edp->fdone!=0)||(edp->cn==edp->cn0)) goto KERNOS;
   // generates new random section
   edp->popa[pi]->Genrnda(thid,ct,dt,edp->popa[pi]->popsz);
   memcpy(&gyp->gen[ct],&gtp->gen[ct],sizeof(bas)*(dt-ct));
   KERNOS:cvar=cvar;
closed
#endif
 
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Pop::Prexline(int thid) opened
int   gi,pi,gy,bi,bo,bx,ii,di,at,bt,ct,gt,ht,vx,vy,vz,pt,aa,bb,cc,sc;
int   val,tmr,force,add,neusect,done,dsold; Coord cdv;
bas  *gen; float af; Dgx curdgx; Exbd *edp=exbd[thid];

//     swc    novel  0+others
// ord tmr  evolved    others
// ord tmr    novel  0

for(pi=0;pi<NPOPUL;pi++) for(gy=0;gy<edp->popa[pi]->popsz;gy++) opened
   gi=pi*POPSZT+gy;
   gen=edp->popa[pi]->guys[gy].gen;

   #if(PGFREEZE==NO)
   // initialize stepeval
   if(edp->cn==0) opened
      edp->stepeval=0; //printf(" ekkor ");
      Cvint2base(&gen[3],&edp->stepeval,CLKMAX,4,&at);
   closed
   #endif

   //if(gi==0) goto EOCEOC;
   edp->popa[pi]->Gendecode(0,edp->popa[pi]->guys[gy].gen,1);

   #if(PGFREEZE==NO)
   // mutation to stepeval
   if(gi!=0) opened
      af=Rnd1(1001)/(float)1000;
      if (af<=PLSTEVAL) opened bb=Rnd1(2); add=1*Intpower(-1,bb); closed
      else add=0;
      edp->stepeval+=add;
      if(edp->stepeval<0) edp->stepeval=0;
      if(edp->stepeval>CLKMAX-1) edp->stepeval=CLKMAX-1;
      Cvint2base(&gen[3],&edp->stepeval,XTMRSQ,4,&at);
   closed
   #endif

   for(di=edp->dgarxf;di<edp->dgarsz;di++) opened
      dsold=(edp->popa[pi]->dsold-BEGNSQ)/DGENSQ;
      neusect=0;
      if((edp->fdone==0)&&(edp->cn!=edp->cn0))
      if((di>=dsold)&&(di<edp->dgarsz)) 
         neusect=1;
      // sets swc=0 if shift has just occurred
      if(neusect==1) opened
         bo=BEGNSQ+(di*DGENSQ); gen[bo+0]=0;
      closed
      if((neusect==1)||(gi!=0)) opened
         memcpy(&curdgx,&edp->dgar[di],sizeof(Dgx));
         bo=BEGNSQ+(curdgx.arpos*DGENSQ); //bi=bx+bj
         #if(PGFREEZE==YA)
         bx=1; force=0;
         // determines xq
         val=di;
         for(ii=-1;ii<STAGES-1;ii++) opened
            if(ii==-1) opened aa=0; bb=edp->frz[0].xe; closed
            else opened aa=edp->frz[ii].xe; bb=edp->frz[ii+1].xe; closed
            if(((aa<=val)&&(val<bb))) opened edp->xq=ii+1; break; closed
         closed
         // Sets timer
         force=0;
         aa=edp->frz[edp->xq].ls; bb=edp->frz[edp->xq].us; cc=bb-aa+1;
         if((curdgx.timer<aa)||(curdgx.timer>bb)) opened
            curdgx.timer=aa+Rnd1(cc); force=1; 
         closed
         #endif
         #if(EVFTEVAL==YA)
         force=1; if (PLGENETM==0) force=0;
         af=Rnd1(1001)/(float)1000;
         if (af<=PLGENETM) opened at=Rnd1(2); add=1*Intpower(-1,at); closed
         else add=0;
         curdgx.timer+=add;
         if(curdgx.timer<0) curdgx.timer=0;
         if(curdgx.timer>(CLKMAX-1)) curdgx.timer=(CLKMAX-1);
         #endif
         if(force==1) opened
            tmr=curdgx.timer;
            Varrescale(&tmr,STAGES,Intpower(4,XTMRSQ));
            bx=(1+XORDSQ); bt=1;
            Cvint2base(&gen[(bo+bx)],XTMRSQ,4,&bt,&tmr);
         closed
         #if(PAROPT!=0)
         // Sets apos
         force=0;
         Coordutil(&vx,&vy,&vz,&curdgx.dgo.dpos,0); done=0;
         do opened
            popbdp->Detplate(thid,&pt,vx,vy,vz);
            if(edp->platmr[pt]<=curdgx.timer) done=1; else
            opened vx=Rnd1(GRIDX); vy=Rnd1(GRIDY); vz=Rnd1(GRIDZ); closed
         closed while(done==0);
         Coordutil(&vx,&vy,&vz,&cdv,1);
         if(memcmp(&curdgx.dgo.dpos,&cdv,sizeof(Coord))) force=1;
         if(force==1) for(ii=0;ii<3;ii++) opened
            bx=DGLPSQ+DGRASQ; 
            sc=Intpower(4,APOSSQ); 
            bt=1; ht=(bo+bx+ii*APOSSQ);
            if(ii==0) opened ct=cdv.x; gt=GRIDX; closed
            if(ii==1) opened ct=cdv.y; gt=GRIDY; closed
            if(ii==2) opened ct=cdv.z; gt=GRIDZ; closed
            Varrescale(&ct,gt,sc);  
            Cvint2base(&gen[ht],APOSSQ,4,&bt,&ct);
         closed
         #endif
      closed
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Pop::Germline(int thid) opened
int   gi,pi,gy,ii,di,hi,n0,n1,cnt,nr,dt,rt,bo,bx,bk,pt,vx,vy,vz,qr;
int   mocy,sva,svb,oka,okb,rnok,lcrstp,pos;
bas  *gen; Dgx curdgx; Dse (*dharap)[DHARLS]; Sgen cursgn;
Exbd *edp=exbd[thid];

for(pi=0;pi<NPOPUL;pi++) for(gy=0;gy<edp->popa[pi]->popsz;gy++) opened
   gi=pi*POPSZT+gy;
   #if  (MOCDEV==NO)
   dt=edp->dhnr0;
   dharap=edp->dharaz; pos=0;
   #elif(MOCDEV==YA)
   dt=edp->dhnral[gi]; if(dt==1) dt=0; // ACHTUNG: ad-hoc!!!
   dharap=edp->dharal; pos=gi;
   #endif

   //if(gi==0) goto CEOCEO;
   edp->popa[pi]->Gendecode(0,edp->popa[pi]->guys[gy].gen,1);

   for(di=edp->dgarxf;di<edp->dgarsz;di++) opened
      nr=0;
      memcpy(&curdgx,&edp->dgar[di],sizeof(Dgx));
      memcpy(&cursgn,&edp->sgar[di],sizeof(Sgen));
      // checks whether sequence matches any MOC in MOC array
      mocy=0; rt=NO;
      for(ii=0;ii<=dt;ii++) opened
         #if  (GRNMODE==NO)
         rt=Mocmatch(&curdgx.res[0],&dharap [pos][ii].moc[0],
            NAPMIN,FSCMAX,&sva,&svb);
         #elif(GRNMODE==YA)
         for(hi=0;hi<CLKMAX;hi++)
            if(cursgn.ifp==dharap[pos][ii].moc[hi]) rt=YA;
         #endif
         if (rt==YA) mocy=1;
      closed
      if((mocy==0)&&(nr<=nr)) opened
         // selects randomly existing moc
         rnok=0; cnt=0;
         do opened
            n0=Rnd1(dt+1); if(n0>dt) n0=dt; cnt++;
            if(dharap[pos][n0].moc[0]!=-1) opened
               oka=1; okb=1; 
               lcrstp=dharap[pos][n0].lcrstp;
               qr=popbdp->Xqdet(thid,edp->dgar[di].arpos);
               if(lcrstp>=edp->frz[qr].us) oka=0;
               #if(PAROPT!=0)
               Coordutil(&vx,&vy,&vz,&edp->drvaro[n0],0);
               popbdp->Detplate(thid,&pt,vx,vy,vz);
               if(edp->platmr[pt]>curdgx.timer) okb=0;
               #endif
               if((oka==1)&&(okb==1)) rnok=1; else rnok=0; 
            closed
         closed while((rnok==0)&&(cnt<50));
         #if  (GRNMODE==NO)
         for(hi=0;hi<MOCLEN;hi++) curdgx.res[hi]=dharap[pos][n0].moc[hi];
         #elif(GRNMODE==YA)
         n1=Rnd1(CLKMAX); cursgn.ifp=dharap[pos][n0].moc[hi];
         #endif
      closed
      gen=&edp->popa[pi]->guys[gy].gen[0];
      if((mocy==0)&&(nr<=nr)&&(rnok==YA)) opened
         #if  (GRNMODE==NO)
         bo=BEGNSQ+(curdgx.arpos*DGENSQ); bx=0;
         // sets activation to off
         if(gy==0) gen[bo+bx]=(bas)0;
         // sets res
         bx=(1+XORDSQ+XTMRSQ+1); bk=0;
         for(hi=0;hi<MOCLEN;hi++)
            Cvint2base(&gen[(bo+bx+bk)],XRESSQ,4,&bk,&curdgx.res[hi]);
         // sets masking bits to active
         bx=(1+XORDSQ+XTMRSQ+1)+(XRESSQ*MOCLEN);
         for(bk=0;bk<MOCLEN;bk++) gen[bo+bx+bk]=(bas)3;
         nr++;
         #elif(GRNMODE==YA)
         bo=BEGNSQ+(curdgx.arpos*DGENSQ); bx=0;
         // sets activation to off
         if(gy==0) gen[bo+bx]=(bas)0;
         // sets ifp
         bx=SGSCSQ; bk=0;
         Cvint2base(&gen[(bo+bx+bk)],&cursgn.ifp,SGIFSQ,4,&bk);
         nr++;
         #endif
      closed
   closed
   //CEOCEO:cvar=cvar;
closed

#if(DEBUG7==YA)
Gpstatsx(this);
#endif

closed
