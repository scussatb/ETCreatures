
#include "danxincl.h"

// Global var
extern Exbd *exbd[NCORES]; extern char cvar;
char   atmpstr[100],btmpstr[100];

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Dopernew(int thid,int gi,int as,int is,int op) opened closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Doper2(int thid,int gi,int as) opened
int  vx,vy,vz,mx,my,mz,ex,ey,ez,rx,ry,rz,ri,ii;
int  cnd,dfound,nhx,nhy,nhz,ct,ot; Coord cd,vcd,mcd,scd,cde;
Exbd *edp=exbd[thid];

// compiles clarcnd (only cells active and sleep)
Forvxyz opened
   edp->clarcnd[vx][vy][vz]=0;
   cnd=Lreader(Lcnd,&edp->clar[vx][vy][vz]);
   if((cnd==CLALIVN)||(cnd==CLALIVE)) edp->clarcnd[vx][vy][vz]=1;
closed

// here considers only the first driver ([0])
ri=0; for(ii=0;ii<DHARLS;ii++) if(edp->drvaro[ii].x>=0) opened
   ex=edp->drvaro[ii].x;
   ey=edp->drvaro[ii].y;
   ez=edp->drvaro[ii].z;
   memcpy(&edp->clardcd[ri],&edp->drvaro[ii],sizeof(Coord));
   edp->clard[ri++]=edp->clar[ex][ey][ez]; 
 closed
if(ri==0) return;

// if clar[vx][vy][vz] is superficial and non-driver and in the vicinity
// no cell is superficial and driver ... then clar[vx][vy][vz] becomes driver
rx=-888; ry=-888; rz=-888;
edp->smocnew=SMOCMA;
nhx=DOP2NS; nhy=DOP2NS; nhz=DOP2NS; if(NDIMS==2) nhz=0;
for(vz=nhz;vz<GRIDZ-nhz;vz++)
for(vx=nhx;vx<GRIDX-nhx;vx++) 
for(vy=nhy;vy<GRIDY-nhy;vy++) opened
   vcd.x=vx; vcd.y=vy; vcd.z=vz;
   // if "in interior grid", non-void, borderline, NON-driver
   if(edp->clarcnd[vx][vy][vz]==1)
   if(Lreader(Ldrv,&edp->clar[vx][vy][vz])==0)
   if(Surfint(as,&vcd,edp)==YA) opened
      // search neighbourhood, start ----------------
      dfound=NO;
      if((abs(rx-vx)<=nhx)&&(abs(ry-vy)<=nhy)&&(abs(rz-vz)<=nhz)) 
         dfound=YA;
      if(dfound==NO) opened
         for(mz=vz-nhz;mz<=vz+nhz;mz++)
         for(mx=vx-nhx;mx<=vx+nhx;mx++) 
         for(my=vy-nhy;my<=vy+nhy;my++) opened
            mcd.x=vx; mcd.y=vy; mcd.z=vz;
            // if "in interior grid",non-void,borderline,driver
            if(edp->clarcnd[mx][my][mz]==1)
            if(Lreader(Ldrv,&edp->clar[mx][my][mz])==1)
            if(Surfint(as,&mcd,edp)==YA) opened
               dfound=YA; rx=mx; ry=my; rz=mz; goto GEOTOS;
            closed
         closed
      closed
      GEOTOS:cvar=cvar;

      if(dfound==NO) opened
         Coordutil(&vx,&vy,&vz,&cde,1);
         ct=Closedr(thid,gi,&cde,ri,&scd);
         Moccopy(&edp->mocnew[0],&edp->dhar[ct].moc[0]);
         edp->smocnew++;
         cd.x=vx; cd.y=vy; cd.z=vz;
         ot=Lreader(Ldhn,&edp->clar[vx][vy][vz]);
         edp->dltcrt=edp->dhar[ot].dltcrt;
         edp->moclp.dhc=0;
         edp->dhor=-1; edp->drv=YA; edp->cnd=CLALIVN; 
         edp->color=edp->clar[cd.x][cd.y][cd.z].col;
         Dharupd (thid,gi,as,-1,&cd,2);
         Daughter(thid,gi,as,-1,&cd,2);
      closed
      // search neighbourhood, end ------------------
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int Body::Closedr(int thid,int gi,Coord *cde,int ri,Coord *cd) opened
int  ex,ey,ez,si,rd,md,mda;
Exbd *edp=exbd[thid];

Coordutil(&ex,&ey,&ez,cde,0);
md=(GRIDX+GRIDY+GRIDZ); rd=-1;
for(si=0;si<ri;si++) opened
   mda=abs(edp->clardcd[si].x-ex)
      +abs(edp->clardcd[si].y-ey)
      +abs(edp->clardcd[si].z-ez);
   if(md>mda) opened md=mda; rd=si; closed
closed

return rd;
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Remred(int thid,int gi,int as,int di,Coord *cdv,int ss) opened
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Loadtgt(int thid,char *file) opened
int  vx,vy,vz,rx,ry,rz,val,ii,sr; FILE *fp; Exbd *edp=exbd[thid];

if((fp=fopen(file,"r"))==NULL) Leavexec("File not found"); ii=0;
for(rz=0;rz<GRIDZ;rz++) for(ry=0;ry<GRIDY;ry++) for(rx=0;rx<GRIDX;rx++) opened
   if((ii==0)||((ii+0)%GRIDX==0)) sr=fscanf(fp,"\n (vz:%d vy:%d vx:%d)",&vz,&vy,&vx);
   sr=fscanf(fp," %d",&val); ii++; vx=rx;
   edp->envr.cdtg[vx][vy][vz]=val;
closed
fclose(fp);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Gridinit(int thid,int gi,int as) opened
int vx,vy,vz,ii,zi,at,bt,dr,cnd,col; Coord cd;
Exbd *edp=exbd[thid];

if(as!=0) return;

zi=0; 
Forvxyz opened
   dr=0;
   if(ISGRID==0) opened
      at=0; bt=0; 
      cnd=CLCHOLE; col=0;
      for(ii=0;ii<ZYGTOT;ii++) 
         if((vx==edp->zygc[ii].x)&&(vy==edp->zygc[ii].y)&&(vz==edp->zygc[ii].z)) opened
         cd.x=vx; cd.y=vy; cd.z=vz; zi=ii; dr=1; break;
      closed
   closed
   if(ISGRID!=0) opened
      at=edp->envr.etis[vx][vy][vz];
      bt=edp->envr.cdis[vx][vy][vz]; 
      cnd=CLALIVE; col=bt;
      if(at==1) dr=1;
   closed
   Lwriter(Lxxx,&edp->clar[vx][vy][vz],0,0,cnd,col);
   if(dr==1) opened // driver
      for(ii=0;ii<MOCLEN;ii++) edp->mocnew[ii]=0;
      edp->smocnew=SMOCMA+zi;
      cd.x=vx; cd.y=vy; cd.z=vz;
      edp->moclp.dhc=0; // should be -1, but unsigned short
      edp->isval=0; edp->dltcrt=as;
      edp->dhor=-1; edp->drv=YA; edp->cnd=CLALIVN;
      edp->color=col;
      Dharupd (thid,gi,as,-1,&cd,0);
      Daughter(thid,gi,as,-1,&cd,0); zi++;
   closed
closed

if(PAROPT!=0) 
   Forvxyz edp->pmap[vx][vy][vz]=edp->envr.imap[vx][vy][vz];

closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
void  Body::Gridupd(int thid,int gi,int as,int us,int ff) opened
int   pi,gy,vx,vy,vz,drv,clcnd,hs,at,bt,ct;
Exbd *edp=exbd[thid];

pi=gi/POPSZT; gy=gi%POPSZT;

if(gi==0) if(ff==0) Forvxyz opened
   at=-1; bt=-1;
   if (edp->clar[vx][vy][vz].dhc!=-1) opened
      drv  =Lreader(Ldrv,&edp->clar[vx][vy][vz]);
      clcnd=Lreader(Lcnd,&edp->clar[vx][vy][vz]);
      if(clcnd==CLALIVE) at=Lreader(Ldhn,&edp->clar[vx][vy][vz]);
      if(clcnd==CLALIVE) bt=Lreader(Lcol,&edp->clar[vx][vy][vz]);
      if(clcnd==CLALIVN) at=Lreader(Ldhn,&edp->clar[vx][vy][vz])+30000;
      if(clcnd==CLALIVN) bt=Lreader(Lcol,&edp->clar[vx][vy][vz]);
      if(drv==1)         at=Lreader(Ldhn,&edp->clar[vx][vy][vz])+60000;
      if(clcnd==CLCHOLE) at=-1;
      if(clcnd==CLCHOLE) bt=-1;
      if(clcnd==CLCHOLN) at=-2;
      if(clcnd==CLCHOLN) bt=-1;
      // copies to sgrd
      #if(VTKOPT==0)
      if (at==-1)    ct=-1; // empty cells
      if (at!=-1)    ct=+1; // normal cells +1
      if (at>=30000) ct=+2; // new normal cells +2
      if (at>=60000) ct=+3; // driver cells
      if (at==-2)    ct=+4; // new holes
      #endif
      #if(VTKOPT==1)
      ct=bt;
      #endif
      #if(VTKOPT==2)
      ct=bt;
      if (bt>=15)    ct=14;
      if (at>=30000) ct=15; // new normal cells +2
      if (at>=60000) ct=14; // driver cells
      #endif
      #if(VTKOPT==3) // stem only
      if (at>=60000) ct=+3; else ct=-1;
      #endif
      #if(VTKOPT==4) // map
      at=edp->pmap[vx][vy][vz]; ct=at; if(at!=-1) ct=1;
      #endif
      edp->envr.sgrd[vx][vy][vz]=ct;  
   closed
closed
if(gi==0) if(ff==0) if(edp->ausdr==YA) opened
   hs=STAGES-1; if(edp->cn!=edp->cn0) hs=edp->frz[edp->qz].us;
   if(as<=hs) Savegrid(thid,0,as); 
closed

closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
float Body::Fitbox(int thid,int gi,Coord *cdl,Coord *cdu,float *shas) opened
int   pi,gy,lx,ly,lz,ux,uy,uz,vx,vy,vz,ci,ides,iins,at,bt;
float shad,bdes,bins,bous,cins,xins,cinsar[COLRNG],shco0,oc,ious,af,bf;
Clsr  *clp,*clf; Exbd *edp=exbd[thid];

Coordutil(&lx,&ly,&lz,cdl,0); Coordutil(&ux,&uy,&uz,cdu,0);
pi=gi/POPSZT; gy=gi%POPSZT;

ides=0; iins=0; ious=0; for(ci=0;ci<COLRNG;ci++) cinsar[ci]=0;
for(vx=lx;vx<=ux;vx++) for(vy=ly;vy<=uy;vy++) for(vz=lz;vz<=uz;vz++) opened
   clp=&edp->clar[vx][vy][vz]; clf=&edp->clarf[vx][vy][vz];
   at= edp->envr.cdtg[vx][vy][vz];
   if(clp->dhc==-1) bt=-1; // For speed
   else bt=((Lreader(Lcnd,clp)==CLALIVE) || (Lreader(Lcnd,clp)==CLALIVN)) ? Lreader(Lcol,clp):-1;
   #if(CFLAG0==1)
   oc=1.0; //if(vy< 65) oc=0.5;
   #endif
   if (at!=-1)            ides+= 1;
   if((at!=-1)&&(bt!=-1)) iins+= 1;
   if((at==-1)&&(bt!=-1)) ious+=oc;
   #if(FITOPA==0)
   for(ci=0;ci<COLRNG;ci++) if((at==ci)&&(bt==ci)) cinsar[ci]++;
   #endif
   #if(FITOPA==1)
   if((at!=-1)&&(bt!=-1)) opened
      bf=0; af=(float)(abs(at-bt));
      if(af<FITTOL) bf=1-(float)0.7*(af/FITTOL);
      if((at<=COLTOL)&&(bt<=COLTOL)) bf=1;
      cinsar[at]+=bf;
   closed
   #endif
closed

// float conv
bdes=(float)ides;
bins=(float)iins;
bous=(float)ious;

cvar=cvar;//anchor
cins=0;
for(ci=0;ci<COLRNG;ci++) opened cins+=cinsar[ci]; closed
if (cins<0) cins=0;
// calc
shco0=(float)SHCOE0;
xins=(float)((1-shco0)*bins+(shco0)*cins);

#if(FITOPB==0)
// calc shad
if (bdes==0) shad=-1;
else shad=(float)(xins-SHCOE2*bous)/bdes;
if(shad<0) shad=0; if(shad>1) shad=1;
// calc shas
af=(float)(xins-SHCOE2*bous); if (af<0) af=0;
*shas=(bdes-af)/SZGRID;
if(*shas<0) *shas=0; if(*shas>1) *shas=1;
#endif

#if(FITOPB==1)
// calc shad
float szgrid=(float)(ux-lx+1)*(uy-ly+1)*(uz-lz+1);
if((szgrid-bdes)==0) shad=-1;
else opened
   shad=(float)(szgrid-bous-bdes-SHCOE2*(bdes-xins))/(szgrid-bdes);
   if(shad<0) shad=0;
closed
// calc shas
af=(szgrid-bous-bdes-SHCOE2*(bdes-xins)); if (af<0) af=0;
*shas=((szgrid-bdes)-af)/SZGRID;
#endif

return shad;
closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
void Body::Smooth3d(int thid,int gi,Coord *exta,Coord *extb,int nhsize) opened
int ti,avgcol,cnd,col,znhsze,vx,vy,vz,ax,ay,az;
Exbd *edp=exbd[thid];

Forvxyz opened
   cnd=edp->clar[vx][vy][vz].cnd;
   edp->klor[vx][vy][vz]=edp->clar[vx][vy][vz].col;
   if((cnd!=CLALIVN)&&(cnd!=CLALIVE)) edp->klor[vx][vy][vz]=-1;
closed

znhsze=nhsize; if(NDIMS==2) znhsze=0;
for(vx=exta->x;vx<=extb->x;vx++)
for(vy=exta->y;vy<=extb->y;vy++)
for(vz=exta->z;vz<=extb->z;vz++) opened
   ti=0; avgcol=0;
   for(ax=vx-nhsize;ax<=vx+nhsize;ax++)
   for(ay=vy-nhsize;ay<=vy+nhsize;ay++)
   for(az=vz-znhsze;az<=vz+znhsze;az++) if(Insidegrid(ax,ay,az)==1) opened
      col=edp->klor[ax][ay][az];
      if(col!=-1) opened avgcol+=col; ti++; closed
   closed
   if(edp->klor[vx][vy][vz]!=-1)if(ti!=0) opened
      avgcol/=ti;
      edp->clar[vx][vy][vz].col=avgcol;
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Loadneur(int thid) opened
char *filename,str[5]; int pnr,ii,xi,yi,lht,at; float af; FILE *fp;
Exbd *edp=exbd[thid];

for(ii=0;ii<NCAT;ii++) for(pnr=0;pnr<NEXTOT;pnr++) opened
   lht = (int)strlen(SMBFN)+20;
   filename = (char *)malloc(lht);
   //printf("\npnr=%d",pnr);
   strcpy(filename,SMBFN); sprintf(str,"%d",ii);
   strcat(filename,str); strcat(filename,"_");
   if(pnr+1<1000) strcat(filename,"0");
   if(pnr+1< 100) strcat(filename,"0");
   if(pnr+1<  10) strcat(filename,"0"); sprintf(str,"%d",pnr+1);
   strcat(filename,str); strcat(filename,".txt");
   fp=fopen(filename,"r"); if(fp==NULL) return;
   for(yi=0;yi<28;yi++) for(xi=0;xi<28;xi++) opened 
      at=fscanf(fp,"%f",&af); //yr=27-yi;
      edp->envr.neurtot[ii][pnr][yi*28+xi]=af; 
   closed   
   fclose(fp);
closed

#if(DEBUG7==YA)
fp=fopen(XDBG3FN,"w");
for(ii=0;ii<4;ii++) for(pnr=0;pnr<100;pnr++) opened
   for(xi=0;xi<784;xi++) opened
      af=edp->envr.neurodt[ii][pnr][xi]; fprintf(fp," %3f",af);
   closed
   fprintf(fp," %3d\n",ii);
closed
fclose(fp);
#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Dsplneur(int thid) opened
int ii,pn,xi,at;
Exbd *edp=exbd[thid];

for(ii=0;ii<NCAT;ii++) for(pn=0;pn<NEXSPL;pn++) opened
   for(xi=0;xi<784;xi++)
      edp->envr.neurspl[ii][pn][xi]=
      edp->envr.neurtot[ii][pn][xi];
closed

#if(HRNDSPL==1)
for(ii=0;ii<NCAT;ii++) opened
   for(pn=NEXSPL/2;pn<NEXSPL;pn++) opened
      at=Rnd1(NEXTOT);
      for(xi=0;xi<784;xi++)
         edp->envr.neurspl[ii][pn][xi]=
         edp->envr.neurtot[ii][at][xi];
   closed
closed
#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Parlr(int thid,int oa,int ob,int da,int db,int *ii) opened
int   oi,di; float af;
Exbd *edp=exbd[thid];

for(oi=oa;oi<ob;oi++) opened
   if(0==0) opened // threshold
      af=edp->parset[(*ii)++];
      edp->whtsmat[oi][oi]=(int)(af*100);
   closed
   for(di=da;di<db;di++) opened // dendrites
      af=edp->parset[(*ii)++];
      edp->whtsmat[oi][di]=(int)(af*100);
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Parmapr(int thid,int gi,int aa,int oo) opened
int   oi,di,at,ii,vx,vy,vz,ci,oa,ob,os,da,db,lr,subnets; 
int   nquadr,n0,n1,n2,n4,li; 
float af; Clsr *clp; FILE *fp1; Exbd *edp=exbd[thid];

if(aa==1) goto SADLES;
for(ii=0;ii<PARSETSZ;ii++) opened edp->parset[ii]=0; edp->parcnt[ii]=0; closed
Forvxyz opened
   clp=&edp->clar[vx][vy][vz];
   at = edp->pmap[vx][vy][vz];
   if(at!=-1) opened
      af =(float)((clp->col/(float)(COLRNG-1)-0.5)*2);
      //if (clp->cnd!=CLALIVN) if (clp->cnd!=CLALIVE) af=0; 
      edp->parset[at]+=af; edp->parcnt[at]+=1; 
   closed
closed
//for(ii=0;ii<PARSETSZ;ii++) if(edp->parcnt[ii]!=0) 
// edp->parset[ii]/=edp->parcnt[ii]; // average
SADLES: cvar=cvar;
#if(PAROPT==1)
// regularisation of parset
for(ii=0;ii<PARSETSZ;ii++) opened
   if(edp->parset[ii]<-1.000) edp->parset[ii]=-1;
   if(edp->parset[ii]>+1.000) edp->parset[ii]=+1;
   if(edp->parset[ii]>-0.001) 
   if(edp->parset[ii]<+0.001) edp->parset[ii]=+0;
closed
// download into neurcnx 
for(oi=0;oi<NDSTOT;oi++) for(di=0;di<NDSTOT;di++) opened
 //edp->whtsmat[oi][di].cat=-1;
   edp->whtsmat[oi][di]=+0;
closed

if(NETMAP==1) opened
   ii=0; if(aa==1) printf("\nPARSUBS: ");
   subnets=4; // subnets=NCAT;
   for(ci=0;ci<subnets;ci++) opened
      if(aa==1) printf(",%d",ii);
      oa=0; ob=0; // to avoid warnings
      for(lr=1;lr<LRSTOT;lr++) opened
         n0= edp->ndsplr[lr-0]; n1=edp->ndsplr[lr-1]; 
         n2=0; if(lr>=2) n2=edp->ndsplr[lr-2];
         os=(edp->ndsplr[lr-0]-edp->ndsplr[lr-1])/subnets;
         if(lr==1) opened da=0; db=784; closed
         if(lr>=2) if(CNXTYPE==0) opened da=n2; db=n1; closed // tot cnx         
         if(lr>=2) if(CNXTYPE==1) opened da=n2; db=ob; closed // asymmetric tot cnx
         if(lr>=2) if(CNXTYPE==2) opened da=oa; db=ob; closed // compartmentalised cnx
         oa= edp->ndsplr[lr-1]+os*(ci+0); 
         ob= edp->ndsplr[lr-1]+os*(ci+1); 
         Parlr(thid,oa,ob,da,db,&ii);
         //if(aa==1) printf("\nsn %d lr %d %4d-%4d-%4d-%4d",ci,lr,oa,ob,da,db);
      closed
      if(aa==1) printf(",%d",ii);
      cvar=cvar;
   closed
   if(aa==1) printf("\n");
   cvar=cvar;
closed

if(NETMAP==2) opened
   ii=0; if(aa==1) printf("\nPARSUBS: ");
   for(li=1;li<LRSTOT;li++) opened
      lr=li-1; if(li==1) lr=LRSTOT-1; 
      nquadr=1;
      for(ci=0;ci<nquadr;ci++) opened
         if(aa==1) printf(",%d",ii);
         n0= edp->ndsplr[lr-0]; n1=edp->ndsplr[lr-1]; 
         n2=0; if(lr>=2) n2=edp->ndsplr[lr-2];
         n4=0; if(lr>=4) n4=edp->ndsplr[lr-4];
         os=(edp->ndsplr[lr-0]-edp->ndsplr[lr-1])/nquadr;
         if(lr==1) opened da=0; db=784; closed
         if(lr>=2) if(CNXTYPE==0) opened da=n2; db=n1; closed // tot cnx         
         if(lr>=2) if(CNXTYPE==1) opened da=n2; db=ob; closed // asymmetric tot cnx
         if(lr>=2) if(CNXTYPE==2) opened da=oa; db=ob; closed // compartmentalised cnx
         oa= edp->ndsplr[lr-1]+os*(ci+0);
         ob= edp->ndsplr[lr-1]+os*(ci+1);
         if(lr==LRSTOT-1) 
            opened da=n4; db=n1; closed 
         Parlr(thid,oa,ob,da,db,&ii);
         if(aa==1) printf(",%d",ii);
         cvar=cvar;
      closed
   closed
   if(aa==1) printf("\n");
   cvar=cvar;
closed

ii=0; if(gi==0) if(edp->ausdr==1) opened 
   fp1=fopen(XDBG1FN,"w");
   for(oi=0;oi<NDSTOT;oi++) for(di=0;di<NDSTOT;di++) opened
      fprintf(fp1," Nin: %4d Nout: %4d wht: %8.5f",di,oi,
         ((float)edp->whtsmat[oi][di])/100);
      ii++; if(ii%10==0) fprintf(fp1,"\n");
   closed
   fclose(fp1);
closed
#endif
#if(PAROPT==2)
float step,cstp,dlt;

// correction
//for(ii=0;ii<PARSETSZ;ii++) if(oo==1) edp->parset[ii]+=edp->pardlt[ii];
// regularisation (-1 < val < +1)
// decoding
for(ii=0;ii<PARSETSZ;ii++) opened 
   edp->citlist[ii].ord=ii;
   edp->citlist[ii].par=edp->parset[ii];
closed
cvar=cvar;
// calc delta
if(oo==0) opened
   Ftsortpar(thid,PARSETSZ,&edp->citlist[0]);
   step=1/(float)PARSETSZ;
   edp->pardlt[0]=-2;
   for(ii=0;ii<PARSETSZ;ii++) opened
      cstp=(ii*step);
      dlt=cstp-edp->citlist[ii].par;
      edp->pardlt[edp->citlist[ii].ord]=dlt;
   closed
closed
#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Parexcl(int thid,int as) opened
int vx,vy,vz,ai,bi,at; Clsr *clp;
Exbd *edp=exbd[thid];

Forvxyz opened
   clp=&edp->clar[vx][vy][vz];
   at = edp->pmap[vx][vy][vz];
   if(at!=-1) opened
      for(ai=0;ai<NPLATES;ai++) opened
         if(at>=edp->parsubt[ai][0])
         if(at< edp->parsubt[ai][1]) break;
      closed
      if(as<edp->platmr[ai]) clp->col=COLRNG/2; //clp->col=0; 
   closed
   if(at!=-1) opened
      Detplate(thid,&bi,vx,vy,vz);
      if(as<edp->platmr[bi]) clp->col=COLRNG/2; //clp->col=0;
      // plate bi is modified by genes whose tmr is  >= edp->platmr[bi] 
      // therefore, before the plate is inactive and its weights are all 0
   closed
closed

closed

#define CLCMEAN 2
//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
float Body::Clcneur(int thid,int gi,int ff,float *mfit,FILE *fp) opened
int   ii,ci,ai,qi,xi,ki,pnr,aa,bb,oc,lc,lr,la,st;
int   oneok,ncat,ncac,exout[NCAT],lastlr,texn,wexn;
float gmean[NCAT];
float pgdar[NCAT],ngdar[NCAT],mgdar[NCAT];
float dgdar[NCAT],sgdar[NCAT],cgdar[NCAT];
float pgood,ngood,mgood,dgood,sgood,cgood,ggood,tgood;
float shad,maxcurct,maxother,maxallct,af,hf,kf,gf,thres,adif,bdif;
FILE *fp0; Exbd *edp=exbd[thid];
int   gqtrga[3]= opened GQTRGA closed;
float mgdcoe,mgdcoa[3]= opened MGDCOA closed; 
float rshcoe,rshcoa[3]= opened RSHCOA closed; 
float rskcoe,rskcoa[3]= opened RSKCOA closed; 
float cgdcoe,cgdcoa[3]= opened CGDCOA closed;

//Drawex(thid,0);

// lastlr determination
int mstones[4],xlastar[4];
sscanf(LASTLR, "M %d %d %d %d A %d %d %d %d",
   &mstones[0],&mstones[1],&mstones[2],&mstones[3],
   &xlastar[0],&xlastar[1],&xlastar[2],&xlastar[3]);
// default value?
for(ii=0;ii<4;ii++) if(edp->gq>=mstones[ii]) lastlr=xlastar[ii];

// ncac determination
int xcacset[4]; 
sscanf(NCACSET,"M %d %d %d %d A %d %d %d %d",
   &mstones[0],&mstones[1],&mstones[2],&mstones[3],
   &xcacset[0],&xcacset[1],&xcacset[2],&xcacset[3]);
if(edp->gq>=0) ncac=xcacset[0]; // default value
for(ii=0;ii<4;ii++) if(edp->gq>=mstones[ii]) ncac=xcacset[ii];

// mbad,dbad,sbad determination
int mbad,xxmbada[3],dbad,xxdbada[3],sbad,xxsbada[3];
sscanf(XXBADA, "M %d %d %d A %d %d %d B %d %d %d C %d %d %d",
   &mstones[0],&mstones[1],&mstones[2],
   &xxmbada[0],&xxmbada[1],&xxmbada[2],
   &xxdbada[0],&xxdbada[1],&xxdbada[2],
   &xxsbada[0],&xxsbada[1],&xxsbada[2]);
// default value?
for(ii=0;ii<3;ii++) if(edp->gq>=mstones[ii]) opened
   mbad=xxmbada[ii]; dbad=xxdbada[ii]; sbad=xxsbada[ii]; 
closed

// mgdcoe et al determination
for(ii=0;ii<3;ii++) if(edp->gq>=gqtrga[ii]) qi=ii;
mgdcoe=mgdcoa[qi]; rshcoe=rshcoa[qi]; 
rskcoe=rskcoa[qi]; cgdcoe=cgdcoa[qi];

ncat=NCAT; thres=0.5;
for(ci=0;ci<ncat;ci++) opened mgdar[ci]=0; cgdar[ci]=0; closed
if(gi==0)if(edp->ausdr==1) fp0=fopen(XDBG0FN,"w");
for(ii=0;ii<ncat;ii++) for(pnr=0;pnr<NEXSPL;pnr++) opened
   for(xi=0;xi<784;xi++)
      edp->neurout[xi]=(int)(edp->envr.neurspl[ii][pnr][xi]*100);
   for(lr=1;lr<LRSTOT;lr++) opened
      la=edp->ndsplr[lr-2]; if(lr==1) la=0; 
      #if(NETMAP==2)
      if(lr==LRSTOT-1) la=edp->ndsplr[lr-4]; // achtung: customised!!!
      #endif
      Calclr(thid,la,edp->ndsplr[lr-1],edp->ndsplr[lr-1],edp->ndsplr[lr]);
   closed
   #if(CLCMEAN==1)
   for(ci=0;ci<ncat;ci++) opened
      avr[ci]=+0; st=(edp->ndsplr[LRSTOT-1]-edp->ndsplr[LRSTOT-2])/ncat;
      aa=edp->ndsplr[LRSTOT-2]+(ci+0)*st;
      bb=edp->ndsplr[LRSTOT-2]+(ci+1)*st;
      for(ki=aa;ki<bb;ki++) avr[ci]+=((float)edp->neurout[ki])/100;
      avr[ci]/=(bb-aa);
   closed
   #endif      
   #if(CLCMEAN==2)
   // calc mean
   for(ci=0;ci<ncat;ci++) gmean[ci]=0;
   lastlr=edp->ndsplr[LRSTOT-1]-edp->ndsplr[LRSTOT-2];
   aa=edp->ndsplr[LRSTOT-2]; bb=aa+lastlr;
   for(xi=aa;xi<bb;xi++) opened
      cvar=cvar;
      for(ci=0;ci<ncat;ci++) opened
         if((xi-aa)%ncat==ci) 
            gmean[ci]+=((float)edp->neurout[xi])/100;
      closed
   closed
   for(ci=0;ci<ncat;ci++) gmean[ci]/=(lastlr/ncat);
   #endif
   // calc outmat
   for(ci=0;ci<ncat;ci++) opened
      if(ci==ii) edp->outmat[ci][NEXSPL*ii+pnr]=+gmean[ci];
      if(ci!=ii) edp->outmat[ci][NEXSPL*ii+pnr]=-gmean[ci];
   closed
   // calc maxcurct,maxother
   maxother=-1; maxallct=-1; for(ci=0;ci<ncat;ci++) exout[ci]=0;
   maxcurct=gmean[ii];
   for(ci=0;ci<ncat;ci++) opened
      if(gmean[ci]>=thres) opened exout[ci]=1; closed
      if(ii!=ci) if(maxother<gmean[ci]) opened maxother=gmean[ci]; oc=ci; closed
      if(ii==ii) if(maxallct<gmean[ci]) opened maxallct=gmean[ci]; lc=ci; closed
   closed
   // calc nexok[ci]
   oneok=0;
   if(maxcurct>maxother) opened //if(maxcurct>=thres) if(maxother<thres)
      oneok=1;
   closed
   // calc cgdar et al.
   adif=Scagl(maxcurct-maxother);
   bdif=Scagl(maxcurct-thres)+Scagl(thres-maxother); 
   if(SCAGL2==YA) bdif=Scagl(bdif);
   edp->mscar[NEXSPL*ii+pnr]=bdif;
   if(lc==ii) mgdar[lc]+=Scagl(maxallct-maxother);
   if(lc!=ii) mgdar[lc]-=Scagl(maxallct-maxcurct);
   cgdar[ii]+=(float)oneok;
   cvar=cvar;
   #if(DEBUG0==YA)
   if((gi==0)&&(edp->ausdr==1)) opened
      fprintf(fp0,"\n%3d) ii=%3d %3d",gi,ii,pnr);
      fprintf(fp0," avr");
      for(ci=0;ci<ncat;ci++) fprintf(fp0," %6.3f",gmean[ci]);
      fprintf(fp0," out");
      for(ci=0;ci<ncat;ci++) fprintf(fp0," %3d"  ,exout[ci]);
      if(oneok==1) fprintf(fp0," ok");
   closed
   #endif
closed
for(ci=0;ci<ncat;ci++) mgdar[ci]/=NEXSPL;
for(ci=0;ci<ncat;ci++) cgdar[ci]/=NEXSPL;
if(gi==0)if(edp->ausdr==1) fclose(fp0);

#if(DEBUG0==YA)
if(edp->ausdr==1) if(gi==0) if(ff==0) opened
   fp0=fopen(XDBG5FN,"w"); ai=0; 
   for(ii=0;ii<ncat;ii++) opened
      fprintf(fp0,"\n\ncat=%2d\n",ii);
      for(pnr=0;pnr<NEXSPL;pnr++) opened
         af=edp->mscar[NEXSPL*ii+pnr]*100;
         if(ai%10==0) fprintf(fp0,"\n");
         fprintf(fp0," %6.1f",af); ai++;
      closed
   closed
   fclose(fp0);
closed
#endif

// calc p,n,d,s
for(ci=0;ci<ncat;ci++) opened pgdar[ci]=0; ngdar[ci]=0; closed
for(ci=0;ci<ncat;ci++) opened
   for(ai=0;ai<ncat*NEXSPL;ai++) opened
      af=edp->outmat[ci][ai];
      if(af>0) pgdar[ci]+=af; else ngdar[ci]+=af/(ncat-1);
   closed
   pgdar[ci]/=NEXSPL; 
   ngdar[ci]/=NEXSPL;
   cvar=cvar;
   dgdar[ci]=pgdar[ci]+ngdar[ci];
   sgdar[ci]=fabs(1-pgdar[ci]+ngdar[ci]);
closed

// calc totals
pgood=Meansorted(0,&pgdar [0],ncat,ncat,edp->ordar);
ngood=Meansorted(0,&ngdar [0],ncat,ncat,edp->ordar);
dgood=Meansorted(0,&dgdar [0],ncat,dbad,edp->ordar);
sgood=Meansorted(1,&sgdar [0],ncat,sbad,edp->ordar);
mgood=Meansorted(0,&mgdar [0],ncat,mbad,edp->ordar);
cgood=Meansorted(0,&cgdar [0],ncat,ncat,edp->ordar);
if(ff==0) opened 
   //if(cgood>=0.00) edp->gresid=0.75; if(cgood>=0.25) edp->gresid=0.50;
   //if(cgood>=0.50) edp->gresid=0.25; if(cgood>=0.75) edp->gresid=0.10;
   edp->gresid=(float)1.00;
closed
edp->gresid=1.00; // whole array
texn=NCAT*NEXSPL; wexn=(int)(edp->gresid*texn);
ggood=Meansorted(0,edp->mscar,texn,wexn,edp->mscar);
gf=(float)0.5+ggood; hf=(float)0.5+dgood; kf=(float)0.5+mgood; //kf=(1-sgood)/2;
tgood=mgdcoe*gf+rshcoe*hf+rskcoe*kf+cgdcoe*cgood;
if(fabs(dgood)<0.001) tgood=0; // to avoid local minimum

// print
if (gi==0)if(ff==0)if(edp->ausdr==1) opened 
   Printf3(fp,"\n");
   for(ii=0;ii<6;ii++) for(ci=0;ci<ncat+1;ci++) opened 
      if(ii==0) opened strcpy(atmpstr,"pgood"); af=(ci<ncat) ? pgdar[ci]:pgood; closed
      if(ii==1) opened strcpy(atmpstr,"ngood"); af=(ci<ncat) ? ngdar[ci]:ngood; closed
      if(ii==2) opened strcpy(atmpstr,"dgood"); af=(ci<ncat) ? dgdar[ci]:dgood; closed
      if(ii==3) opened strcpy(atmpstr,"sgood"); af=(ci<ncat) ? sgdar[ci]:sgood; closed
      if(ii==4) opened strcpy(atmpstr,"mgood"); af=(ci<ncat) ? mgdar[ci]:mgood; closed
      if(ii==5) opened strcpy(atmpstr,"cgood"); af=(ci<ncat) ? cgdar[ci]:cgood; closed
      strcpy(btmpstr,"\n     "); strcat(btmpstr,atmpstr); strcat(btmpstr," by cat:");
      if(ci==0) Printf3(fp,btmpstr);
      af*=100; sprintf(btmpstr," %6.2f",af); Printf3(fp,btmpstr);
      if(ci==ncat) Printf3(fp,"^"); 
   closed
   sprintf(atmpstr,"\ngf,hf,kf,c-,tgood:"); 
   Printf3(fp,atmpstr);
   sprintf(atmpstr,"  %5.2f*%6.2f  %5.2f*%6.2f  %5.2f*%6.2f  %5.2f*%6.2f %6.2f$ (%4.2f)",
      mgdcoe,gf*100,rshcoe,hf*100,rskcoe,kf*100,cgdcoe,cgood*100,tgood*100,edp->gresid); 
   Printf3(fp,atmpstr);
closed

// final assignment
*mfit=cgood; shad=tgood; 
//shad=(shad+1.0)/2.0; 
if (shad<0) shad=0;
if((shad<0)||(shad>1)) Leavexec("shad error");
return shad;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Calclr(int thid,int la,int lb,int ua,int ub) opened
int  oi,di,sum;
Exbd *edp=exbd[thid];

for(oi=ua;oi<ub;oi++) opened
   sum=(edp->whtsmat[oi][oi]*100); // threshold       
   for(di=la;di<lb;di++) opened 
      sum+=(edp->whtsmat[oi][di]*edp->neurout[di]);
   closed
   edp->neurout[oi]=(int)(Sigmoide(sum,4)*100);
   cvar=cvar;
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Drawex(int thid,int gi) opened
int ii,pnr,xi,yi,dx,dy,at; float af;
Exbd *edp=exbd[thid];

for(xi=0;xi<GRIDX;xi++) for(yi=0;yi<GRIDY;yi++)
   edp->envr.sgrd[xi][yi][0]=-1;

dx=0; dy=0;
for(ii=0;ii<NCAT;ii++) for(pnr=0;pnr<NEXSPL;pnr++) opened
   for(yi=0;yi<28;yi++) for(xi=0;xi<28;xi++) opened
      af=edp->envr.neurspl[ii][pnr][yi*28+xi];
      at=(int)((1-af)*COLRNG);
      edp->envr.sgrd[xi+dx][yi+dy][0]=at;
   closed
   dx+=33; if (dx>=825) opened dx=0; dy+=33; closed
closed
Savegrid(thid,0,0);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Tspinit (int thid) opened
/*int  ii;
Exbd *edp=exbd[thid];

for(ii=0;ii<PARSETSZ;ii++) opened 
   edp->envr.citpos[ii].x=Rnd1(10);
   edp->envr.citpos[ii].y=Rnd1(10);
   edp->envr.citpos[ii].z=0; 
closed*/ 

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
float Body::Calctsp(int thid,int gi) opened
/*long int ii, dst,expdst; float res;
Exbd *edp=exbd[thid];

expdst=(10/2*2)*(PARSETSZ-1);
Ftsortpar(thid,PARSETSZ,&edp->citlist[0]);
dst=0;
for(ii=1;ii<PARSETSZ;ii++) opened
   dst+=
      abs(edp->envr.citpos[edp->citlist[ii-0].ord].x-
          edp->envr.citpos[edp->citlist[ii-1].ord].x)+
      abs(edp->envr.citpos[edp->citlist[ii-0].ord].y-
          edp->envr.citpos[edp->citlist[ii-1].ord].y)+
      abs(edp->envr.citpos[edp->citlist[ii-0].ord].z-
          edp->envr.citpos[edp->citlist[ii-1].ord].z);
closed 

res=1-(float)dst/expdst; 
if(res<0) res=0;
return res;*/
return 0;
closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
void Body::Detplate(int thid,int *pi,int x,int y,int z) opened
int ti; Exbd *edp=exbd[thid];

//#if(PAROPT!=0)
ti=*pi;
for(ti=0;ti<NPLATES;ti++) opened
   if((x>=edp->platbnd[ti][0].x)&&(x<edp->platbnd[ti][1].x)) 
   if((y>=edp->platbnd[ti][0].y)&&(y<edp->platbnd[ti][1].y))
   if((z>=edp->platbnd[ti][0].z)&&(z<edp->platbnd[ti][1].z)) break; 
closed
*pi=ti;
//#endif

//if((vx<0)||(vx>=GRIDX)) return 0;
//if((vy<0)||(vy>=GRIDY)) return 0;
//if((vz<0)||(vz>=GRIDZ)) return 0;

closed
