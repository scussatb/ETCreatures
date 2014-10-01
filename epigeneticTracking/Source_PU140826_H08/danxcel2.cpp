
#include "danxincl.h"

// Global var
extern Exbd *exbd[NCORES]; extern char cvar;

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Dgarprep(int thid,int gi,int as,int ss) opened
int    vx,vy,vz,di,ii,ntypes,result,nap,fsc,fs,tmr,ot,cnd,dltmin,dltmax;
double ad; Clsr *clpv; Exbd *edp=exbd[thid];

for(ii=0;ii<DHARLS;ii++)
   memcpy(&edp->drvar1[ii],&edp->drvaro[ii],sizeof(Coord));

// Inits --------------------------------------------------------------------
for(di=0;di<edp->dgarsz;di++) memcpy(&edp->dgare[di],&edp->idgx,sizeof(Dgx));
for(ii=0;ii<DHARLS;ii++) edp->auxarp[ii]=-1;
if(as==ss) for(ii=0;ii<DHARLS;ii++) opened
   edp->drvaro[ii].x=-1; edp->drvaro[ii].y=-1; edp->drvaro[ii].z=-1;
closed

ot=1; //if((PGFREEZE==NO)&&(as>0)&&(edp->evtnr[as-1]==0)) ot=0;
if (ot==1)
Forvxyz opened
   clpv=&edp->clar[vx][vy][vz];
   if(clpv->dhc!=-1) opened
      cnd=Lreader(Lcnd,clpv);
      if(cnd==CLALIVN) Lwriter(Lxxx,clpv,-1,-1,CLALIVE,-1);
      if(cnd==CLCHOLN) Lwriter(Lxxx,clpv,-1,-1,CLCHOLE,-1);
      if(as==ss)
      if((cnd==CLALIVN)||(cnd==CLALIVE))
      if(Lreader(Ldrv,clpv)==1) opened
         ot=Lreader(Ldhn,clpv);
         edp->drvaro[ot].x=vx; edp->drvaro[ot].y=vy; edp->drvaro[ot].z=vz;
      closed
   closed
closed

// calc dltcrt
dltmin=8888; dltmax=0;
for(ii=0;ii<DHARLS;ii++) if(edp->drvaro[ii].x!=-1) opened
   vx=edp->drvaro[ii].x; vy=edp->drvaro[ii].y; vz=edp->drvaro[ii].z;
   ot=Lreader(Ldhn,&edp->clar[vx][vy][vz]);
   if(dltmax < edp->dhar[ot].dltcrt) dltmax = edp->dhar[ot].dltcrt;
   if(dltmin > edp->dhar[ot].dltcrt) dltmin = edp->dhar[ot].dltcrt;
closed

// Sorts by exord
ad=(float)Intpower(4,XORDSQ)-1;
for(di=0;di<edp->dgarsz;di++)
   edp->par[di]=(float)edp->dgar[di].exord/(float)ad;
//Mysortxdgx(thid,edp->par,edp->dgarsz,1,edp->dgare,edp->dgar);
Ftsortxdgx(thid,edp->dgarsz,-1,edp->dgare,edp->dgar);
edp->dgarsze=edp->dgarsz;

// eliminates inactive, undue & "colour" genes
memcpy(&edp->dgaro[0],&edp->dgare[0],sizeof(Dgx)*edp->dgarsz);
edp->dgarszo=edp->dgarsze; edp->dgarsze=0;
fs=ESTAGE;
for(di=0;di<edp->dgarszo;di++) opened
   tmr=as-edp->dgaro[di].timer;
   if((edp->dgaro[di].swtch!=-1)&&(edp->dgaro[di].swtch!=0))
   if((edp->dgaro[di].timer==-1)||((tmr>=dltmin)&&(tmr<=dltmax)))
   if((as-edp->frz[edp->xqar[edp->dgaro[di].arpos]].us<=dltmax)||(as>fs))
   if((edp->dgaro[di].dgo.msa>=0)&&(edp->dgaro[di].dgo.msa<=MSAMAX))
      memcpy(&edp->dgare[edp->dgarsze++],&edp->dgaro[di],sizeof(Dgx));
closed

ntypes=0;
for(ii=0;ii<DHARLS;ii++) if(edp->drvaro[ii].x>=0) edp->auxarp[ntypes++]=ii;
cvar=cvar;//anchor

#if(DEBUG0==YA)
for(ii=0;ii<DHARLS;ii++) opened
   if(as!=ss) if((edp->drvar1[ii].x!=-1)||(edp->drvaro[ii].x!=-1))
   if((edp->drvar1[ii].x!=edp->drvaro[ii].x)||
      (edp->drvar1[ii].y!=edp->drvaro[ii].y)||
      (edp->drvar1[ii].z!=edp->drvaro[ii].z)) printf("cxxxo");
closed
#endif

// associates matching mocs to genes
for(di=0;di<edp->dgarsze;di++) opened
   edp->dgare[di].dhptx=0;
   for(ii=0;ii<ntypes;ii++) opened
      result=Mocmatch(&edp->dgare[di].res[0],&edp->dhar[edp->auxarp[ii]].moc[0],
         NAPMIN,FSCMAX,&nap,&fsc);
      if(result==YA) opened
         // if it finds matching moc for gene, records the moc's dhnrx
         edp->dgare[di].dhnrx[edp->dgare[di].dhptx]=edp->auxarp[ii];
         edp->dgare[di].dhptx++;
      closed
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Shaper(int thid,int gi,int as,int us) opened
int   vx,vy,vz,di,ii,ri,zi,at,ot;
int   evtnra,evtnrb,goahd,matched;
char  str[20]; double coe,hex,hey,hez,max,min;
Clsr  curcl; Dgo *dgo; Coord cdv; Exbd *edp=exbd[thid];

strcpy(str,"ahi ahi restore");

#if(ZSTEPS==YA)
if(gi==0) if(thid==0) opened
   edp->bstepsgrid = new dimensa[GRIDX];
   edp->estepsgrid = new dimensa[GRIDX];
   edp->xstepsgrid = new dimensa[GRIDX];
closed

if(gi==0) if(thid==0) opened
   Gridupd(thid,gi,as,us,0);
   Forvxyz opened
      edp->bstepsgrid[vx][vy][vz]=edp->envr.sgrd[vx][vy][vz];
   closed
closed
#endif

edp->cgok=1; evtnra=0; evtnrb=0; // change event nr
// the limit evtnr<NEVENT is put here because some x could be discarded because ko
for(di=0;((di<edp->dgarsze)&&((evtnra+evtnrb)<NEVENT));di++)
//for(zi=0;(zi<edp->dgare[di].dhptx);zi++) opened
for(zi=0;zi<1;zi++) if(edp->dgare[di].dgo.msa<8) opened
   matched=NO;
   // scans the grid
   for(ri=0;ri<DHARLS;ri++) if(edp->drvaro[ri].x!=-1) opened
      memset(&curcl,0,sizeof(Clsr));
      Lwriter(Lxxx,&curcl,ri,+1,CLALIVE,-1);
      memcpy(&cdv,&edp->drvaro[ri],sizeof(Coord));
      Coordutil(&vx,&vy,&vz,&cdv,0);
      if(NDIMS==2) if(vz!=0) printf(" z error");
      
      goahd=1;
      #if(GRNMODE==NO)
      ot=Lreader(Ldhn,&edp->clar[vx][vy][vz]);
      at=edp->dhar[ot].dltcrt;
      edp->ctr=edp->dgare[di].timer; // ctr is current gene timer 
      Detplate(thid,&edp->curpt,vx,vy,vz);
      cvar=cvar;
      // ACHTUNG: with Matchet it can match more cells,with dhnrx JUST ONE!!!
      //if(((Matchet(&edp->dgare[di].res[0],&moc[0],NAPMIN,FSCMAX,&sva,&svb)==NO))) ....
      if((Lreader(Lcnd,&curcl)!=CLALIVE)||(Lreader(Ldrv,&curcl)!=1)) opened goahd=0; goto DECIDE; closed
      if((edp->cgok!=1)||(evtnra>evtnra))                            opened goahd=0; goto DECIDE; closed
      if (edp->dhar [ot].actdrv!=1)                                  opened goahd=0; goto DECIDE; closed
      if (edp->dgare[di].dhnrx[zi]!=Lreader(Ldhn,&curcl))            opened goahd=0; goto DECIDE; closed
      if((edp->dgare[di].timer!=-1)&&(edp->dgare[di].timer!=as-at))  opened goahd=0; goto DECIDE; closed
    //if((edp->dgare[di].dvfrd==0)&&(as>=REGFDS))                    opened goahd=0; goto DECIDE; closed
      if(PAROPT!=0) if(edp->platmr[edp->curpt]>edp->ctr)             opened goahd=0; goto DECIDE; closed
      // plate curpt can be worked on by genes with timer edp->platmr[curpt] onwards
      edp->clok=1; matched=YA;
      edp->platchgd[edp->curpt]=1;
      #endif
      #if(GRNMODE==YA)
      if(edp->dhar[ri].dgo.gas==0)
         opened goahead=0; goto DECIDE; closed
      #endif
      DECIDE:cvar=cvar;
      if(goahd==1) opened
         #if(GRNMODE==NO)
         dgo=&edp->dgare[di].dgo;
         #endif
         #if(GRNMODE==YA)
         dgo=&edp->dhar[ri].dgo;
         #endif
         // Calc geometry
         // coordinates are relative ((vx,vy) is in (0,0)), movement of pts makes a "z"
         memcpy(&edp->actdgr.dpl[0],&dgo->dpl[0],sizeof(Coord)*8);
         memcpy(&edp->actdgr.col[0],&dgo->col[0],sizeof( int )*9);
         memcpy(&edp->actdgr.  dpos,&dgo->  dpos,sizeof(Coord)*1);
         // compiles ellissoid parameters
         hex=(double)abs((edp->actdgr.dpl[7].x-edp->actdgr.dpl[0].x))/2;
         hey=(double)abs((edp->actdgr.dpl[7].y-edp->actdgr.dpl[0].y))/2;
         hez=(double)abs((edp->actdgr.dpl[7].z-edp->actdgr.dpl[0].z))/2;
         // max,min
         max=(hex>=hey) ? hex : hey; max=(max>=hez) ? max : hez;
         min=(hex<=hey) ? hex : hey; min=(min<=hez) ? min : hez;
         coe=(min/max);
         // compiles transformation matrix
         for(ii=0;ii<9;ii++) edp->actdgr.tr[ii]=dgo->tr[ii];
         #if(SPHERE==2)
         edp->tr[0]=1; edp->tr[1]=0; edp->tr[2]=0; // no rot
         edp->tr[3]=0; edp->tr[4]=1; edp->tr[5]=0; // no rot
         edp->tr[6]=0; edp->tr[7]=0; edp->tr[8]=1; // no rot
         #endif
         // calc for remred, checks coordinates?
         #if(NOAPOPT==YA)
         if(edp->dgare[di].dgo.msa/2==1) edp->dgare[di].dgo.msa=0;
         #endif
         if(edp->dgare[di].dgo.msa/2<=2) evtnra++;
         if(edp->dgare[di].dgo.msa/2==3) evtnrb++;
         if(edp->dgare[di].dgo.msa/2<=2) opened
            // Copies mother into tmp var
            memcpy(&edp->moclp,&edp->clar[vx][vy][vz],sizeof(Clsr));
            edp->mothers[edp->motnr][0]=ot; edp->motnr++;
            if(edp->dgare[di].dgo.msa/2==0) Remred(thid,gi,as,di,&cdv,0);
         closed
         #if(DEBUG7==YA)
         int hx,hy,hz;
         if(edp->dgare[di].dgo.msa/2==3)
         for(hx=0;hx<GRIDX;hx++) for(hy=0;hy<GRIDY;hy++) for(hz=0;hz<GRIDZ;hz++) 
            edp->envr.sgrd[hx][hy][hz]=edp->pmap[hx][hy][hz]; 
         Savegrid(thid,0,0);
         #endif
         Brush(thid,gi,as,di,&cdv);
         Clset(thid,gi,as,di,&cdv);
         #if(DEBUG7==YA)
         if(edp->dgare[di].dgo.msa/2==3)
         for(hx=0;hx<GRIDX;hx++) for(hy=0;hy<GRIDY;hy++) for(hz=0;hz<GRIDZ;hz++) 
            edp->envr.sgrd[hx][hy][hz]=edp->pmap[hx][hy][hz]; 
         Savegrid(thid,0,0);
         #endif
         if(edp->dgare[di].dgo.msa/2<=2) opened
            if(edp->dgare[di].dgo.msa/2==0) Remred(thid,gi,as,di,&cdv,1);
            // deactivates mother
            edp->dhar[ot].actdrv=0;

            // Puts back mother (mothers of apoptosis events are not put back)
            if((NOCANC!=0)||(edp->dgare[di].dgo.msa/2!=1)) opened 
               memcpy(&edp->clar[vx][vy][vz],&edp->moclp,sizeof(Clsr));
               // updates array of drivers (it could have been deleted during the event)
               ot=Lreader(Ldhn,&edp->clar[vx][vy][vz]);
               edp->drvaro[ot].x=vx;
            closed

            // updates dvfrd
            edp->dgar[edp->dgare[di].arpos].dvfrd=1;
            #if(TAGCHK==YA)
            // updates frgenesc
            if(as<=edp->stepeval) opened
               memcpy(&edp->frgenesc[gi][(NEVENT*as)+evtnra-1],&edp->dgare[di],sizeof(Dgx));
               edp->frgenesc[gi][(NEVENT*as)+evtnra-1].exeas=as;
               edp->frgenesc[gi][(NEVENT*as)+evtnra-1].stepeval=edp->stepeval;
            closed
            #endif
         closed
      closed
   closed
   if(GRNMODE==YA) goto KAPLOS; // does cycle on drivers just once
   cvar=cvar;//anchor
   // Movie steps
closed

KAPLOS:cvar=cvar;

Movie(thid,gi,as,us,evtnra);

if(gi==0) opened edp->evtnra0[as]=evtnra; edp->evtnrb0[as]=evtnrb; closed

#if(DEBUG0==YA)
if(edp->cgok!=1) opened
   printf(" gi: %d dhnr: %d",gi,edp->dhnr); Leavexec(str); 
closed
#endif

#if(ZSTEPS==YA)
if(gi==0) if(thid==0) opened
   delete edp->bstepsgrid; delete edp->estepsgrid; delete edp->xstepsgrid;
closed
#endif

cvar=cvar;//anchor for debugger
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Brush(int thid,int gi,int as,int di,Coord *cdpv) opened
int   vx,vy,vz,bx,by,bz,dx,dy,dz,kx,ky,kz,crtdrv,ot,pt,dltcrt,a,b,ii,xi,yi,zi;
double hex,hey,hez,fx,fy,fz,ad,dstcoe[9],sumcoe; 
Dgo *dgo; Coord pnt[9];
Exbd *edp=exbd[thid];

Coordutil(&vx,&vy,&vz,cdpv,0);
ot=Lreader(Ldhn,&edp->moclp);
dltcrt=edp->dhar[ot].dltcrt;

#if(GRNMODE==NO)
dgo=&edp->dgare[di].dgo;
#endif
#if(GRNMODE==YA)
dgo=&edp->dhar [ot].dgo;
#endif

#if  (DOPNEW==NO)
crtdrv=1; if(NOSCCH==2) crtdrv=0;
#elif(DOPNEW==YA)
if(as-dltcrt==1) at=1; else at=0;
#endif

// ellissoid: semi-axes, centre and displacement
hex=(double)abs((edp->actdgr.dpl[7].x-edp->actdgr.dpl[0].x))/2;
hey=(double)abs((edp->actdgr.dpl[7].y-edp->actdgr.dpl[0].y))/2;
hez=(double)abs((edp->actdgr.dpl[7].z-edp->actdgr.dpl[0].z))/2;
kx=2*DPLMAX/2; dx=(int)hex+edp->actdgr.dpl[0].x;
ky=2*DPLMAX/2; dy=(int)hey+edp->actdgr.dpl[7].y;
kz=2*DPLMAX/2; dz=(int)hez+edp->actdgr.dpl[0].z;
// xxxx
fx=1; //fx=(1.00+(0.33*(double)edp->dgare[di].dgo.dgd[7]/63));
fy=1; //fy=(1.00+(0.33*(double)edp->dgare[di].dgo.dgd[7]/63));
fz=1; //fz=(1.00+(0.33*(double)edp->dgare[di].dgo.dgd[7]/63));

// for speed
if(NDIMS==2) opened a=-1; b=-1; hez=1; kz=0; dz=0; closed
if(NDIMS==3) opened a=+1; b=+1; closed

// computes extremes
for(xi=0;xi<2;xi++) for(yi=0;yi<2;yi++) for(zi=0;zi<2;zi++) opened
   pnt[4*xi+2*yi+zi].x=kx-Intpower(-1,xi)*(int)hex; pnt[8].x=kx; 
   pnt[4*xi+2*yi+zi].y=ky-Intpower(-1,yi)*(int)hey; pnt[8].y=ky; 
   pnt[4*xi+2*yi+zi].z=kz-Intpower(-1,zi)*(int)hez; pnt[8].z=kz;
closed
// set color extremes
memset(&edp->bufc[0],0,sizeof(int)*(2*DPLMAX)*(2*DPLMAX)*(2*DPLMAX));
for(ii=0;ii<9;ii++) edp->bufc[pnt[ii].x][pnt[ii].y][pnt[ii].z]=edp->actdgr.col[ii];

// xxxx
for(bz=(kz-(int)hez-a);bz<=(kz+(int)hez+b);bz++)
for(by=(ky-(int)hey-1);by<=(ky+(int)hey+1);by++)
for(bx=(kx-(int)hex-1);bx<=(kx+(int)hex+1);bx++) opened
   edp->bufs[bx][by][bz]=6; // no signal
   // ellissoid
   if((pow((bx-kx)/(hex*fx),2)+pow((by-ky)/(hey*fy),2)+pow((bz-kz)/(hez*fz),2))<=1) opened
      if((dgo->msa/2==0)||(dgo->msa/2==2)) edp->bufs[bx][by][bz]=7; // normal signal
      if (dgo->msa/2==1) edp->bufs[bx][by][bz]=8; // delete signal
      if((dgo->msa/2==0)||(dgo->msa/2==2)) if(crtdrv==1)
      if((bx%N2DRAT==0)&&(by%N2DRAT==0)&&(bz%N2DRAT==0)) edp->bufs[bx][by][bz]=9; // driver signal
      if (dgo->msa/2==3) edp->bufs[bx][by][bz]=7; // for swap
      // colour
      sumcoe=0; pt=1;
      for(ii=0;ii<9;ii++) opened
         ad=double(abs(bx-pnt[ii].x)+abs(by-pnt[ii].y)+abs(bz-pnt[ii].z));
         if(ad!=0) dstcoe[ii]=1/ad; else opened pt=0; break; closed
         sumcoe+=dstcoe[ii];
      closed
      for(ii=0;ii<9;ii++) dstcoe[ii]/=sumcoe;
      if (pt==1) 
      for(ii=0;ii<9;ii++) opened
         edp->bufc[bx][by][bz]+=(int)(edp->actdgr.col[ii]*dstcoe[ii]);
      closed
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Clset(int thid,int gi,int as,int di,Coord *cdpv) opened
int   ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,kx,ky,kz,vx,vy,vz,ab,cb;
int   ix,iy,iz,a,b,done,bufval,first,ct,dt,ot,tmpvar,cola,colc;
double hex,hey,hez,tr[9],hem,rx,ry,rz,ix0,ix1,iy0,iy1,iz0,iz1;
Coord cd; Clsr *clpa,*clpc; Exbd *edp=exbd[thid];

Coordutil(&vx,&vy,&vz,cdpv,0);
Moccopy(&edp->mocm[0],&edp->dhar[Lreader(Ldhn,&edp->moclp)].moc[0]);
edp->smocnew=0; done=NO; edp->dhor=edp->dhnr+1;
memcpy(&tr[0],&edp->actdgr.tr[0],sizeof(double)*9);
ot=Lreader(Ldhn,&edp->moclp);
edp->dltcrt=edp->dhar[ot].dltcrt;

// ellissoid: semi-axes, centre and displacement
hex=(double)abs((edp->actdgr.dpl[7].x-edp->actdgr.dpl[0].x))/2;
hey=(double)abs((edp->actdgr.dpl[7].y-edp->actdgr.dpl[0].y))/2;
hez=(double)abs((edp->actdgr.dpl[7].z-edp->actdgr.dpl[0].z))/2;
kx=2*DPLMAX/2; dx=(int)hex+edp->actdgr.dpl[0].x;
ky=2*DPLMAX/2; dy=(int)hey+edp->actdgr.dpl[7].y;
kz=2*DPLMAX/2; dz=(int)hez+edp->actdgr.dpl[0].z;
hem=Max3(hex+1,hey+1,hez+1);

// for speed
if(NDIMS==2) opened a=-1; b=-1; hez=1; kz=0; dz=0; closed
if(NDIMS==3) opened a=+1; b=+1; closed

// xxxx
for(ix=0;ix<2;ix++) for(iy=0;iy<2;iy++) for(iz=0;iz<2;iz++)
for(bz=(kz-(int)hez-a);bz<=(kz+(int)hez+b);bz++)
for(by=(ky-(int)hey-1);by<=(ky+(int)hey+1);by++)
for(bx=(kx-(int)hex-1);bx<=(kx+(int)hex+1);bx++) opened
   first=NO;
   if((bz==(kz-(int)hez-a))&&(iz==0)) 
   if((by==(ky-(int)hey-1))&&(iy==0))
   if((bx==(kx-(int)hex-1))&&(ix==0)) first=YA;
   bufval=edp->bufs[bx][by][bz];
   if(bufval!=6) opened
      rx=(tr[0]*(bx-kx)+tr[1]*(by-ky)+tr[2]*(bz-kz));
      ry=(tr[3]*(bx-kx)+tr[4]*(by-ky)+tr[5]*(bz-kz));
      rz=(tr[6]*(bx-kx)+tr[7]*(by-ky)+tr[8]*(bz-kz));

      // to make nore robust: if difference is too big, takes the other boundary
      ix0=floor(rx); //if(fabs(ix0-rx)>0.98) ix0=ceil (rx);
      ix1=ceil (rx); //if(fabs(ix1-rx)>0.98) ix1=floor(rx);
      iy0=floor(ry); //if(fabs(iy0-ry)>0.98) iy0=ceil (ry);
      iy1=ceil (ry); //if(fabs(iy1-ry)>0.98) iy1=floor(ry);
      iz0=floor(rz); //if(fabs(iz0-rz)>0.98) iz0=ceil (rz);
      iz1=ceil (rz); //if(fabs(iz1-rz)>0.98) iz1=floor(rz);

      if(ix==0) ax=vx+dx+(int)ix0; if(ix==1) ax=vx+dx+(int)ix1;
      if(iy==0) ay=vy+dy+(int)iy0; if(iy==1) ay=vy+dy+(int)iy1;
      if(iz==0) az=vz+dz+(int)iz0; if(iz==1) az=vz+dz+(int)iz1;
      if(NDIMS==2) az=0;
      clpa=&edp->clar[ax][ay][az];
      Coordutil(&ax,&ay,&az,&cd,1);
      Detplate(thid,&ab,ax,ay,az);
      cvar=cvar;
      #if(PAROPT!=0)
      if(edp->ctr>=edp->platmr[ab])
      #endif
      if(Insidegrid(ax,ay,az)==1) opened
         // normal and delete
         if(edp->dgare[di].dgo.msa/2<=2)
         if((bufval==7)||(bufval==8)) opened  // normal (7) or delete (8) signal
            if(bufval==7) opened Moccopy(&edp->mocnew[0],&edp->mocm[0]); edp->cnd=CLALIVN; closed
            if(bufval==8) opened Moccopy(&edp->mocnew[0],&edp->imoc[0]); edp->cnd=CLCHOLN; closed
            if(done==YA) goto DENTOS;
            if(bufval==7) opened edp->smocnew=+1; closed  // normal
            if(bufval==8) opened edp->smocnew=-1; closed  // delete
            edp->isval=0; edp->drv=NO;
            Dharupd(thid,gi,as,di,&cd,1);
            done=YA;
            DENTOS:cvar=cvar;
            dt=1;
            #if(NOCANC==1)
            ct=Lreader(Lcnd,clpa);
            if(bufval==7)
            if((ct==CLALIVE)||(ct==CLALIVN)) dt=0;
            #endif
            #if(NOCANC==2)
            ct=Lreader(Lcnd,clpa);
            if((bufval==7)||(bufval==8))
            if((ct==CLALIVE)||(ct==CLALIVN)) dt=0;
            #endif
            if(dt==0) goto BELFOS;
            edp->drv=NO;
            edp->color=edp->bufc[bx][by][bz];
            Daughter(thid,gi,as,di,&cd,1);
            BELFOS:cvar=cvar;
         closed
         // drivers
         if(edp->dgare[di].dgo.msa/2<=2)
         if((ix==1)&&(iy==1)&&(iz==1)) // last cycle to avoid drv get cancelled
         if((bufval==9)||(first==YA)) opened         // driver signal
            edp->smocnew++; 
            Moccopy(&edp->mocnew[0],&edp->mocm[0]);
            //if(first==YA) Coordutil(&vx,&vy,&vz,&cd,1);
            edp->isval=0; edp->drv=YA;
            edp->dhor=-1; edp->cnd=CLALIVN; 
            edp->color=edp->bufc[bx][by][bz];
            Dharupd (thid,gi,as,di,&cd,1);
            Daughter(thid,gi,as,di,&cd,1);
         closed
         #if(PAROPT!=0)
         // swap operator
         if((ix==1)&&(iy==1)&&(iz==1)) // last cycle to avoid repetitions
         if((edp->dgare[di].dgo.msa==6)||(edp->dgare[di].dgo.msa==7)) opened
            cx=ax-vx+edp->actdgr.dpos.x;
            cy=ay-vy+edp->actdgr.dpos.y;
            cz=az-vz+edp->actdgr.dpos.z;
            clpc=&edp->clar[cx][cy][cz];
            Detplate(thid,&cb,cx,cy,cz);
            cvar=cvar;
            if(edp->ctr>=edp->platmr[cb])
            if(Insidegrid(cx,cy,cz)==1) opened
            //if(Cboundchk(thid,edp->curgb,cx,cy,cz)==1) opened
               tmpvar= edp->pmap[ax][ay][az];
               edp->pmap[ax][ay][az]=edp->pmap[cx][cy][cz];
               edp->pmap[cx][cy][cz]=tmpvar;
               //edp->pmap[ax][ay][az]=0; edp->pmap[cx][cy][cz]=0;
               if(edp->dgare[di].dgo.msa==7) opened
                  cola=Lreader(Lcol,clpa); 
                  colc=Lreader(Lcol,clpc);
                  Lwriter(Lcol,clpa,-1,-1,-1,colc); 
                  Lwriter(Lcol,clpc,-1,-1,-1,cola);
               closed         
            closed         
         closed
         #endif
      closed
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Daughter(int thid,int gi,int as,int di,Coord *cdpa,int clr) opened
int  ax,ay,az,xcol,pcol,pdrv,pcnd,pdhn,at,bt,ct,kt,conr,c2nr,msb;
Exbd *edp=exbd[thid]; Clsr *clpa;

Coordutil(&ax,&ay,&az,cdpa,0); clpa=&edp->clar[ax][ay][az];
if(NDIMS==2) opened if(az!=0) printf(" z error"); az=0; closed

msb =edp->dgare[di].dgo.msb;
conr=edp->dgare[di].dgo.conr;
c2nr=edp->dgare[di].dgo.c2nr;

pdhn=Lreader(Ldhn,clpa); pdrv=Lreader(Ldrv,clpa); // previous status
pcnd=Lreader(Lcnd,clpa); pcol=Lreader(Lcol,clpa);

if((pcnd==CLALIVE)||(pcnd==CLALIVN)) opened
   // mothers can't be deleted
   if(NOCANC==1) if(edp->dhar[pdhn].sigstp!=-1) return; 
   // all alive cells can't be deleted
   if(NOCANC==2) return; 
closed

if((NOSCCH!=2)||(clr==0)) opened
// updates drvaro
   if(pdrv==1) opened
      pdhn=Lreader(Ldhn,&edp->clar[ax][ay][az]);
      if(memcmp(&edp->drvaro[pdhn],cdpa,sizeof(Coord))==0) 
         edp->drvaro[pdhn].x=-1;
   closed
closed

// dhn,drv
if(edp->drv==0) opened at=edp->dhor; bt=0; closed
if(edp->drv==1) opened at=edp->dhnr; bt=1; closed
// cnd
if(edp->cnd==CLCHOLN) if(pcnd==CLCHOLE) edp->cnd=CLCHOLE;
ct=edp->cnd;
// col
#if(PAROPT==0)
kt=conr;
#endif
#if(PAROPT!=0)
xcol=edp->color; //xcol=conr;
if(di==-1) msb=0;
if((msb <0)||(msb >3)) msb=0;
if (msb>=2) if(c2nr<=COLRNG/2) xcol*=-1;
if((msb==0)||(msb==1)) kt=xcol; // absolute
if((msb==2)||(msb==3)) kt=xcol+pcol; // delta
if (kt<=0) kt=0; if(kt>=COLRNG) kt=COLRNG-1;
//if((msb<=1)||((msb>=2)&&(pcnd!=CLALIVN)))
#endif

if((NOSCCH==2)&&(clr!=0)) opened at=-1; bt=-1; closed
Lwriter(Lxxx,clpa,at,bt,ct,kt);
if(edp->moclp.dhc!=-1) clpa->mother=edp->moclp.dhc;

if((NOSCCH!=2)||(clr==0)) opened
// updates drvaro
   if(edp->drv==1)if(edp->cnd==CLALIVN) opened
      pdhn=edp->dhnr;
      memcpy(&edp->drvaro[pdhn],cdpa,sizeof(Coord));
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Dharupd(int thid,int gi,int as,int di,Coord *cda,int clr) opened
int  hi,at,ht,smocmax,min; Dgx *dgx; Exbd *edp=exbd[thid];

dgx=&edp->dgare[di]; if (di==-1) dgx=&edp->idgx;

min=SMOCMB+1;
for(hi=0;hi<MOCLEN;hi++) if(min>edp->mocnew[hi]) opened
   min=edp->mocnew[hi]; ht=hi;
closed
edp->mocnew[ht]=edp->smocnew;
if(ht>=MOCLEN) Leavexec("MOCLEN overflow");
//edp->mocnew[as-edp->dltcrt]=edp->smocnew;

at=1;
if(edp->dhnr>=0)
   at=memcmp(&edp->dhar[edp->dhnr].moc[0],&edp->mocnew[0],sizeof(int)*MOCLEN);

if(at) opened
   if((edp->dhnr+1)<=(DHARLS-1)) opened
      edp->dhnr++;
      for(hi=0;hi<MOCLEN;hi++) edp->dhar[edp->dhnr].moc[hi]=edp->mocnew[hi];
      edp->dhar[edp->dhnr].clr=clr;
      edp->dhar[edp->dhnr].lcrstp=as;
      edp->dhar[edp->dhnr].lcrisc=edp->isval;
      edp->dhar[edp->dhnr].dltcrt=edp->dltcrt;
      edp->dhar[edp->dhnr].actdrv= (edp->drv==1)? 1:0;
      memcpy(&edp->dhar[edp->dhnr].cd,cda,sizeof(Coord));
   closed
   else edp->cgok=0;
closed

if(clr==1) smocmax=SMOCMA; else smocmax=SMOCMB;
if(edp->smocnew>smocmax) opened
   printf("smocmax %d %d %d %d %d",thid,gi,as,clr,edp->smocnew); edp->cgok=0; 
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Movie(int thid,int gi,int as,int us,int cgevnr) opened
#if(ZSTEPS==YA)
int   act,nbs,vx,vy,vz,qx,qy,qz,xi,yi,zi,pi,qi,tmp,cc,at,bt;
long  zpts,zptstot,zptsres,cfound,cfoundmax;
double zprob,rprob,adv;
char  cvar=0;
Exbd  *edp=exbd[thid];

if((gi!=0)||(thid!=0)) return;

Gridupd(thid,gi,as,us,0);
Forvxyz opened
   edp->estepsgrid[vx][vy][vz]=edp->envr.sgrd[vx][vy][vz];
closed

//Calcolo different pnts
for(pi=0;pi<DVPTSMAX;pi++) edp->divpar[pi]=0;
zptstot=0;
Forvxyz opened
   act=NO;
   #if  (VTKOPT==0)
   at=edp->estepsgrid[vx][vy][vz];
   bt=edp->bstepsgrid[vx][vy][vz];
   if (at!=bt) opened 
      act=YA;
      if((at==-1)&&(bt!=-1)) edp->divpar[zptstot]=1; // first absent, last present
      if((at!=-1)&&(bt==-1)) edp->divpar[zptstot]=2; // first present,last absent 
   closed
   #elif(VTKOPT==1)
   at=edp->estepsgrid[vx][vy][vz];
   bt=edp->bstepsgrid[vx][vy][vz];
   if (at!=bt) opened 
      act=YA;
      if((at==-1)&&(bt!=-1)) edp->divpar[zptstot]=1; // first absent, last present
      if((at!=-1)&&(bt==-1)) edp->divpar[zptstot]=2; // first present, last absent
      if((at!=-1)&&(bt!=-1)) edp->divpar[zptstot]=3; // both present, diff colours 
   closed 
   #endif
   if(act==YA) opened
      edp->divpts[zptstot].x=vx;
      edp->divpts[zptstot].y=vy;
      edp->divpts[zptstot].z=vz;
      zptstot++;
   closed
closed
zpts=zptstot/NSTEPS;

Forvxyz opened
   edp->xstepsgrid[vx][vy][vz]=edp->bstepsgrid[vx][vy][vz];
closed

if(zptstot>0)
for(qi=0;qi<(NSTEPS-1);qi++) opened
   // calc probs
   cfoundmax=0;
   for(pi=0;pi<zptstot;pi++) if(edp->divpar[pi]!=-1) opened
      cfound=0; nbs=1;
      vx=edp->divpts[pi].x;
      vy=edp->divpts[pi].y;
      vz=edp->divpts[pi].z;
      for(qz=vz-nbs;qz<=vz+nbs;qz++)
      for(qx=vx-nbs;qx<=vx+nbs;qx++)
      for(qy=vy-nbs;qy<=vy+nbs;qy++)
         if(Insidegrid(qx,qy,qz)==1)
         if(edp->xstepsgrid[qx][qy][qz]!=-1) cfound++;
      if(cfound>cfoundmax) cfoundmax=cfound;
   closed

   zptsres=zpts; pi=0; cc=0;
   while((zptsres>=0)&&(cc<100)) opened
      if(edp->divpar[pi]!=-1) opened
         cfound=0; nbs=1;
         vx=edp->divpts[pi].x;
         vy=edp->divpts[pi].y;
         vz=edp->divpts[pi].z;
         for(qz=vz-nbs;qz<=vz+nbs;qz++)
         for(qx=vx-nbs;qx<=vx+nbs;qx++)
         for(qy=vy-nbs;qy<=vy+nbs;qy++)
            if(Insidegrid(sx,sy,sz)==1)
            if(edp->xstepsgrid[qx][qy][qz]!=-1) cfound++;
         adv=cfound/(double)cfoundmax;
         if(edp->divpar[pi]==1) zprob=0+adv; // if full,  dense volume is privileged
         if(edp->divpar[pi]==2) zprob=1-adv; // if hole, sparse volume is privileged
         if(edp->divpar[pi]==3) zprob=0+adv; // if none,  dense volume is privileged
         rprob=Rnd1(100)/(double)100;
         if(rprob<zprob) opened
            xi =edp->divpts[pi].x;
            yi =edp->divpts[pi].y;
            zi =edp->divpts[pi].z;
            edp->xstepsgrid[xi][yi][zi]=
            edp->estepsgrid[xi][yi][zi];
            edp->divpar[pi]=-1;
            zptsres--;
            //printf("%d ",zptsres);
         closed
      closed
      pi++; if(pi>=zptstot) opened pi=0; cc++; closed
   closed
   Forvxyz opened
      edp->envr.sgrd[vx][vy][vz]=edp->xstepsgrid[vx][vy][vz];
   closed
   Savegrid(thid,1,as);
closed

// End-of-step snapshot
Gridupd(thid,gi,as,us,0);
Savegrid(thid,1,as);
if(as==(CLKMAX-1)) exit(0);

#endif

closed
