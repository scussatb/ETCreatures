
#include "danxincl.h"

// Global var
extern Exbd *exbd[NCORES]; extern char cvar;
static Dgx dgarh[POPSZT][AMSDGN]; static float divar[POPSZT][AMSDGN]; 
static float mtar[POPSZT]; 
#if(DEBUG0==YA)
double voxgrid[GRIDX][GRIDY][GRIDZ];
int  ordar[DHARLS]; float dstar[DHARLS];
#endif

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Sendrec0(int nr,char *buffer,long *buflen0) opened
int  ei,vx,vy,vz,ii,pn,nd;
#if(MYMPI==1)
int  thid,si,rc; MPI_Status stat;
Exbd *edp;
#endif

#if(MYMPI==0)
for(ei=1;ei<NCORES;ei++) opened
   exbd[ei]->cn0=exbd[0]->cn0;
   Forvxyz opened
      exbd[ei]->envr.cdtg[vx][vy][vz]=exbd[0]->envr.cdtg[vx][vy][vz];
      exbd[ei]->envr.etis[vx][vy][vz]=exbd[0]->envr.etis[vx][vy][vz];
      exbd[ei]->envr.cdis[vx][vy][vz]=exbd[0]->envr.cdis[vx][vy][vz];
   closed
   if(PAROPT!=0) Forvxyz
      exbd[ei]->envr.imap[vx][vy][vz]=exbd[0]->envr.imap[vx][vy][vz];
   if(PAROPT==1) // neural
   for(ii=0;ii<NCAT;ii++) for(pn=0;pn<NEXTOT;pn++) for(nd=0;nd<784;nd++)
      exbd[ei]->envr.neurtot[ii][pn][nd]=exbd[0]->envr.neurtot[ii][pn][nd];
   if(PAROPT==2) opened closed // tsp
 //for(ii=0;ii<PARSETSZ;ii++) 
 //   memcpy(&exbd[ei]->envr.citpos[0],&exbd[0]->envr.citpos[0],sizeof(Coord));
closed
#endif
#if(MYMPI==1)
thid=Getthid(&edp);
if (thid==0) opened
   si=0; Mpscpy(&buffer[si],&edp->cn0,sizeof(int),&si);
   Forvxyz opened
      Mpscpy(&buffer[si],&edp->envr.cdtg[vx][vy][vz],sizeof(int),&si);
      Mpscpy(&buffer[si],&edp->envr.etis[vx][vy][vz],sizeof(int),&si);
      Mpscpy(&buffer[si],&edp->envr.cdis[vx][vy][vz],sizeof(int),&si);
   closed
   if(PAROPT!=0) Forvxyz
      Mpscpy(&buffer[si],&edp->envr.imap[vx][vy][vz],sizeof(int),&si);
   if(PAROPT==1) // neural
   for(ii=0;ii<NCAT;ii++) for(pn=0;pn<NEXTOT;pn++) for(nd=0;nd<784;nd++)
      Mpscpy(&buffer[si],&edp->envr.neurtot[ii][pn][nd],sizeof(float),&si);
   if(PAROPT==2) opened closed // tsp
 //for(ii=0;ii<PARSETSZ;ii++) 
 //   Mpscpy(&buffer[si],&edp->envr. citpos[ii],sizeof(Coord),&si);
   for(ei=1;ei<NCORES;ei++) rc=MPI_Send(buffer,
      *buflen0,MPI_BYTE,ei,0,MPI_COMM_WORLD);
closed
else opened
   rc=MPI_Recv(buffer,*buflen0,MPI_BYTE, 0,0,MPI_COMM_WORLD,&stat);
   si=0; Mprcpy(&buffer[si],&edp->cn0,sizeof(int),&si);
   Forvxyz opened
      Mprcpy(&buffer[si],&edp->envr.cdtg[vx][vy][vz],sizeof(int),&si);
      Mprcpy(&buffer[si],&edp->envr.etis[vx][vy][vz],sizeof(int),&si);
      Mprcpy(&buffer[si],&edp->envr.cdis[vx][vy][vz],sizeof(int),&si);
   closed
   if(PAROPT!=0) Forvxyz
      Mprcpy(&buffer[si],&edp->envr.imap[vx][vy][vz],sizeof(int),&si);
   if(PAROPT==1) // neural
   for(ii=0;ii<NCAT;ii++) for(pn=0;pn<NEXTOT;pn++) for(nd=0;nd<784;nd++)
      Mprcpy(&buffer[si],&edp->envr.neurtot[ii][pn][nd],sizeof(float),&si);
   if(PAROPT==2) opened closed // tsp
 //for(ii=0;ii<PARSETSZ;ii++) 
 //   Mprcpy(&buffer[si],&edp->envr. citpos[ii],sizeof(Coord),&si);
closed
#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Sendrecv(int nr,char *buffer,long *buflen1,long *buflen3) opened
int  thid,ei,pi,gy,gi,gl,gu,at,bt;
double bd; Exbd *edp; Guy *gyp,*gzp;
#if(MYMPI==1)
int si,rc; MPI_Status stat;
#endif

if(nr==1) opened
   thid=Getthid(&edp);
   at=edp->popa[0]->dnaxf; 
   bt=edp->popa[0]->dnasz;
   *buflen1 =sizeof(int)*6;
   *buflen1+=sizeof(int)*7;
   *buflen1+=sizeof(Guy);
   *buflen1+=sizeof(int)* NNGUYS;
   *buflen1+=sizeof(bas)*(NNGUYS/NCORES+1)*(bt-at);
   *buflen1+=sizeof(bas)*(NNGUYS/NCORES+1)*BEGNSQ;
   if (edp->ausdr==1) fprintf(edp->cfp," buflen1=%ld",*buflen1);
   #if(MYMPI==0)
   for(ei=1;ei<NCORES;ei++) opened
      exbd[ei]->cn=exbd[0]->cn; exbd[ei]->dgarxf=exbd[0]->dgarxf;
      exbd[ei]->gq=exbd[0]->gq; exbd[ei]->dgarsz=exbd[0]->dgarsz;
      memcpy(&exbd[ei]->guyf,&exbd[0]->guyf,sizeof(Guy));
      memcpy(&exbd[ei]->rndpmt[0],&exbd[0]->rndpmt[0],sizeof(int)*NNGUYS);
      for(pi=0;pi<NPOPUL;pi++) for(gy=0;gy<POPSZT;gy++) opened
         gyp =&exbd[ei]->popa[pi]->guys[gy];
         gzp =&exbd[ 0]->popa[pi]->guys[gy];
         memcpy(&gyp->gen[at],&gzp->gen[at],sizeof(bas)*(bt-at));
         memcpy(&gyp->gen[ 0],&gzp->gen[ 0],sizeof(bas)*BEGNSQ);
      closed
   closed
   #endif
   #if(MYMPI==1)
   if (thid==0) opened
      for(ei=1;ei<NCORES;ei++) opened
         bd=NNGUYS/(double)NCORES;
         gl=(int)(bd*ei); gu=(int)(bd*(ei+1));
         if (ei==NCORES-1) gu=NNGUYS;
         si=0;
         Mpscpy(&buffer[si],&edp->cn,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->gq,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->dgarxf,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->dgarsz,sizeof(int),&si);
         Mpscpy(&buffer[si],&gl,sizeof(int),&si); //boundaries
         Mpscpy(&buffer[si],&gu,sizeof(int),&si); //boundaries
         #if(CFLAG0==YA)
         Mpscpy(&buffer[si],&at,sizeof(int),&si);
         Mpscpy(&buffer[si],&bt,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->qf,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->qz,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->popa[0]->dsold,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->popa[0]->dnaxf,sizeof(int),&si);
         Mpscpy(&buffer[si],&edp->popa[0]->dnasz,sizeof(int),&si);
         #endif
         Mpscpy(&buffer[si],&edp->guyf,sizeof(Guy),&si);
         Mpscpy(&buffer[si],&edp->rndpmt[0],sizeof(int)*NNGUYS,&si);
         for(gi=gl;gi<gu;gi++) opened
            pi=gi/POPSZT; gy=gi%POPSZT;
            gyp=&edp->popa[pi]->guys[gy];
            Mpscpy(&buffer[si],&gyp->gen[at],sizeof(bas)*(bt-at),&si);
            Mpscpy(&buffer[si],&gyp->gen[ 0],sizeof(bas)*BEGNSQ,&si);
         closed
         rc=MPI_Send(buffer,*buflen1,MPI_BYTE,ei,1,MPI_COMM_WORLD);
      closed
   closed
   else opened
      if (0==0) opened // to keep indentation
         si=0;
         rc=MPI_Recv(buffer,*buflen1,MPI_BYTE,0,1,MPI_COMM_WORLD,&stat);
         Mprcpy(&buffer[si],&edp->cn,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->gq,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->dgarxf,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->dgarsz,sizeof(int),&si);
         Mprcpy(&buffer[si],&gl,sizeof(int),&si); //boundaries
         Mprcpy(&buffer[si],&gu,sizeof(int),&si); //boundaries
         #if(CFLAG0==YA)
         Mprcpy(&buffer[si],&at,sizeof(int),&si);
         Mprcpy(&buffer[si],&bt,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->qf,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->qz,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->popa[0]->dsold,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->popa[0]->dnaxf,sizeof(int),&si);
         Mprcpy(&buffer[si],&edp->popa[0]->dnasz,sizeof(int),&si);
         #endif
         Mprcpy(&buffer[si],&edp->guyf,sizeof(Guy),&si);
         Mprcpy(&buffer[si],&edp->rndpmt[0],sizeof(int)*NNGUYS,&si);
         for(gi=gl;gi<gu;gi++) opened
            pi=gi/POPSZT; gy=gi%POPSZT;
            gyp=&edp->popa[pi]->guys[gy];
            Mprcpy(&buffer[si],&gyp->gen[at],sizeof(bas)*(bt-at),&si);
            Mprcpy(&buffer[si],&gyp->gen[ 0],sizeof(bas)*BEGNSQ,&si);
         closed
      closed
   closed
   #endif
closed

if(nr==3) opened
   thid=Getthid(&edp);
   bd=NNGUYS/(double)NCORES;
   #if(MYMPI==0)
   // consolidates pop
   for(ei=1;ei<NCORES;ei++) opened
      gl=(int)(bd*ei); gu=(int)(bd*(ei+1));
      if (ei==NCORES-1) gu=NNGUYS;
      for(gi=gl;gi<gu;gi++) opened
         pi=gi/POPSZT; gy=gi%POPSZT;
         exbd[0]->popa[pi]->guys[gy].shad=exbd[ei]->popa[pi]->guys[gy].shad;
         exbd[0]->popa[pi]->guys[gy].mfit=exbd[ei]->popa[pi]->guys[gy].mfit;
         #if(MOCDEV==YA)
         memcpy(&exbd[0]->dharal[gi][0],&exbd[ei]->dharal[gi][0],sizeof(Dse)*DHARLS);
         exbd[0]->dhnral[gi]=exbd[ei]->dhnral[gi];
         #endif
         #if(TAGCHK==YA)
         memcpy(&exbd[0]->frgenesc[gi][0],
            &exbd[ei]->frgenesc[gi][0],sizeof(Dgx)*(NEVENT*STAGES));
         #endif
      closed
   closed
   #endif
   #if(MYMPI==1)
   gl=(int)(bd*thid); gu=(int)(bd*(thid+1));
   if (thid==NCORES-1) gu=NNGUYS;
   if (thid!=0) opened
      si=0;
      Mpscpy(&buffer[si],&gl,sizeof(int),&si); // boundaries
      Mpscpy(&buffer[si],&gu,sizeof(int),&si); // boundaries
      for(gi=gl;gi<gu;gi++) opened
         pi=gi/POPSZT; gy=gi%POPSZT;
         gyp=&edp->popa[pi]->guys[gy];
         Mpscpy(&buffer[si],&gyp->shad,sizeof(float),&si);
         Mpscpy(&buffer[si],&gyp->mfit,sizeof(float),&si);
         #if(MOCDEV==YA)
         Mpscpy(&buffer[si],&edp->dharal[gi][0],sizeof(Dse)*DHARLS,&si);
         Mpscpy(&buffer[si],&edp->dhnral[gi],sizeof(int),&si);
         #endif
         #if(TAGCHK==YA)
         Mpscpy(&buffer[si],&edp->frgenesc[gi][0],sizeof(Dgx)*(NEVENT*CLKMAX),&si);
         #endif
      closed
      rc=MPI_Send(buffer,*buflen3,MPI_BYTE,0,3,MPI_COMM_WORLD);
   closed
   else for(ei=1;ei<NCORES;ei++) opened
      si=0;
      rc=MPI_Recv(buffer,*buflen3,MPI_BYTE,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,&stat);
      Mprcpy(&buffer[si],&gl,sizeof(int),&si); // boundaries
      Mprcpy(&buffer[si],&gu,sizeof(int),&si); // boundaries
      for(gi=gl;gi<gu;gi++) opened
         pi=gi/POPSZT; gy=gi%POPSZT;
         gyp=&edp->popa[pi]->guys[gy];
         Mprcpy(&buffer[si],&gyp->shad,sizeof(float),&si);
         Mprcpy(&buffer[si],&gyp->mfit,sizeof(float),&si);
         #if(MOCDEV==YA)
         Mprcpy(&buffer[si],&edp->dharal[gi][0],sizeof(Dse)*DHARLS,&si);
         Mprcpy(&buffer[si],&edp->dhnral[gi],sizeof(int),&si);
         #endif
         #if(TAGCHK==YA)
         Mprcpy(&buffer[si],&edp->frgenesc[gi][0],sizeof(Dgx)*(NEVENT*CLKMAX),&si);
         #endif
      closed
   closed
   #endif
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Calcfit(int thid) opened
int   pi,gy,di,gz,dk,at,pres,popsz;
float eh,th,mt,shad,ftavg,ftmax,ftmin,af;
int   bstpos; float bstval; Guy *bstgyp;
Pop   *poptr; FILE *fp; Exbd *edp=exbd[thid];

popsz=edp->popa[0]->popsz;
// assigns gaft
for(pi=0;pi<NPOPUL;pi++) for(gy=0;gy<popsz;gy++) opened
   shad=edp->popa[pi]->guys[gy].shad;
   edp->popa[pi]->guys[gy].gaft=shad;
closed

// temporary crowning (each pop)
for(pi=0;pi<NPOPUL;pi++) opened
   bstval=-1; bstgyp=NULL;
   for(gy=0;gy<popsz;gy++) opened
      if(bstval< edp->popa[pi]->guys[gy].gaft) opened
         bstval= edp->popa[pi]->guys[gy].gaft;
         bstgyp=&edp->popa[pi]->guys[gy]; bstpos=gy;
      closed
   closed
   // Transplants best guys in pos 0
   memcpy(&edp->popa[pi]->guys[0],bstgyp,sizeof(Guy));
   #if(TAGCHK==YA)
   memcpy(&edp->frgenesc[pi*POPSZT+0][0],&edp->frgenesc[bstpos][0],
      sizeof(Dgx)*(NEVENT*STAGES));
   #endif
closed

// MOC-based niching
for(pi=0;pi<NPOPUL;pi++) opened
   for(gy=0;gy<popsz;gy++) opened
      edp->popa[pi]->Gendecode(0,edp->popa[pi]->guys[gy].gen,1);
      for(di=0;di<AMSDGN;di++) 
         memcpy(&dgarh[gy][di],&edp->dgar[edp->dgarxf+di],sizeof(Dgx));
      cvar=cvar;
   closed
   for(gy=0;gy<popsz;gy++) for(di=0;di<AMSDGN;di++) opened
      eh=0;
      if(dgarh[gy][di].swtch==1) opened         
         for(gz=0;gz<popsz;gz++) opened
            pres=NO;
            for(dk=0;dk<AMSDGN;dk++) opened
               if(dgarh[gz][dk].swtch==1) opened
                  pres=NO;
                  memcpy(&edp->resi[0],&dgarh[gy][di].res[0],sizeof(int)*MOCLEN);
                  memcpy(&edp->resk[0],&dgarh[gz][dk].res[0],sizeof(int)*MOCLEN);                  
                  if((memcmp(edp->resi,edp->imoc,sizeof(int)*MOCLEN))!=0) 
                  if((memcmp(edp->resi,edp->resk,sizeof(int)*MOCLEN))==0) if(gy!=gz) opened
                     pres=YA; break;
                  closed
               closed
            closed
            if(pres==YA) eh++;
         closed
      closed
      divar[gy][di]=eh/popsz;
   closed
   for(gy=0;gy<popsz;gy++) opened
      at=0; mt=0;
      for(di=0;di<AMSDGN;di++) if(dgarh[gy][di].swtch==1) opened
         at++; mt+=divar[gy][di];
      closed
      if(at!=0) mt/=at; else mt=1; mtar[gy]=(1-NICHECOE*mt);
      #if(NICHEX==YA)
      edp->popa[pi]->guys[gy].gaft=
      edp->popa[pi]->guys[gy].shad*mtar[gy];
      #endif
   closed
   #if(DEBUG0==YA)
   if(edp->ausdr==1) opened
      fp=fopen(XDBG6FN,"w");
      for(gy=0;gy<popsz;gy++) opened
         mt=(1-mtar[gy])/NICHECOE;
         fprintf(fp,"\n%3d) %6.3f %6.3f@",gy,mtar[gy],mt);
         for(di=0;di<AMSDGN;di++)
            fprintf(fp," %6.3f",divar[gy][di]);
      closed
      fclose(fp);
   closed
   #endif
closed

for(pi=0;pi<NPOPUL;pi++) opened
   poptr=edp->popa[pi];
   // equidistant fitness values
   //step=1/(float)(poptr->popsz-1);
   //for(gy=0;gy<poptr->popsz;gy++) poptr->guys[gy].gaft=(gy*step);
   // dynamic fitness values
   ftavg=0; ftmax=-1; ftmin=1; 
   for(gy=0;gy<poptr->popsz;gy++) opened
      af=poptr->guys[gy].gaft; ftavg+=af/poptr->popsz; 
      if(ftmax<af) ftmax=af; 
      if(ftmin>af) ftmin=af; 
   closed
   for(gy=0;gy<poptr->popsz;gy++) opened
      if(ftmax-ftmin!=0)
         af=(poptr->guys[gy].gaft-ftmin)/(ftmax-ftmin);
      else af=1;
      if(af<0) af=0; if(af>1) af=1;
      poptr->guys[gy].gaft=af;
   closed
   #if(DEBUG0==YA)
   if(edp->ausdr==1) opened
      fp=fopen(XDBG2FN,"w"); if(fp!=NULL) opened
         for(gy=0;gy<poptr->popsz;gy++) opened
            if(gy%5==0) fprintf(fp,"\n");
            at=fprintf(fp," %3d) %7.4f %7.4f",gy,
               poptr->guys[gy].shad,poptr->guys[gy].gaft);
         closed
      closed
      fclose(fp);
   closed
   #endif
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Showres(int thid,FILE *fptr00,char *atmpstr) opened
int   gi,pi,gy,ii,nshowed,*ord;
float *par,shad,gaft,mfit;
float gaftavg[NPOPUL],shadavg[NPOPUL],mfitavg[NPOPUL],shadbst;
Pop   *poptr; Exbd *edp=exbd[thid];

ord=( int *)malloc(sizeof( int )*edp->popa[0]->popsz);
par=(float*)malloc(sizeof(float)*edp->popa[0]->popsz);

nshowed=SHOWED;
for(pi=0;pi<NPOPUL;pi++) opened
   poptr=edp->popa[pi];
 //for(gy=0;gy<poptr->popsz;gy++) par[gy]=poptr->guys[gy].gaft;
   for(gy=0;gy<poptr->popsz;gy++) par[gy]=poptr->guys[gy].shad;
   Mysortxxxx(&par[0],poptr->popsz,0,&ord[0]);
   for(ii=0;ii<nshowed;ii++) opened
      gaft=poptr->guys[ord[ii]].gaft*100;
      shad=poptr->guys[ord[ii]].shad*100;
      mfit=poptr->guys[ord[ii]].mfit*100;
      if(pi==0)if(ii==0) shadbst=shad;
      strcpy(atmpstr,"\nBEST %3d) (%3d) gaft: %5.2f shad: %5.2f mfit: %5.2f ");
      fprintf(stdout,atmpstr,ii,ord[ii],gaft,shad,mfit);
      fprintf(fptr00,atmpstr,ii,ord[ii],gaft,shad,mfit);
   closed
   // avg
   gaftavg[pi]=0; shadavg[pi]=0; mfitavg[pi]=0;
   for(gy=0;gy<poptr->popsz;gy++) opened
      gaftavg[pi]+=poptr->guys[gy].gaft;
      shadavg[pi]+=poptr->guys[gy].shad;
      mfitavg[pi]+=poptr->guys[gy].mfit;
   closed
   gaftavg[pi]=gaftavg[pi]*100/POPSZT;
   shadavg[pi]=shadavg[pi]*100/POPSZT;
   mfitavg[pi]=mfitavg[pi]*100/POPSZT;
   strcpy(atmpstr,"   AVG gaft: %5.2f shad: %5.2f mfit: %5.2f ");
   fprintf(stdout,atmpstr,gaftavg[pi],shadavg[pi],mfitavg[pi]);
   fprintf(fptr00,atmpstr,gaftavg[pi],shadavg[pi],mfitavg[pi]);
closed

// History
if((edp->cn-1)%HSTSTP==0) opened
   for(gi=0;gi<NNGUYS;gi++) opened
      pi=gi/POPSZT; gy=gi%POPSZT;
      edp->fithst[(edp->cn-1)/HSTSTP][gi]=poptr->guys[gy].mfit;
      edp->shdhst[(edp->cn-1)/HSTSTP][gi]=poptr->guys[gy].shad;
   closed
closed

free(ord); free(par);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Loadmemo(int thid) opened
int  pi,gy,bi,ii,at,sr,gi; FILE *fp; float af;
Exbd *edp=exbd[thid];

// HSTFN
fp=fopen(HSTFN,"r"); if(fp!=NULL) opened
   for(ii=0;ii<FITHST;ii++) opened
      sr=fscanf(fp,"%d",&at);
      for(gi=0;gi<NNGUYS;gi++) sr=fscanf(fp,"%f",&edp->fithst[ii][gi]);
   closed
   for(ii=0;ii<FITHST;ii++) opened
      sr=fscanf(fp,"%d",&at);
      for(gi=0;gi<NNGUYS;gi++) sr=fscanf(fp,"%f",&edp->shdhst[ii][gi]);
   closed
closed
if(fp!=NULL) fclose(fp);

// Genome
fp=fopen(GENFN,"r"); if(fp!=NULL) opened
   sr=fscanf(fp,"\nCYCLE %d\n",&edp->cn0);
   for(pi=0;pi<NPOPUL;pi++) opened
      for(gy=0;gy<edp->popa[pi]->popsz;gy++) opened
         sr=fscanf(fp,"%d)",&at);
         sr=fscanf(fp,"%f", &af); edp->popa[pi]->guys[gy].shad=af;
         sr=fscanf(fp,"%f", &af); edp->popa[pi]->guys[gy].gaft=af;
         // developmental genes
         for(bi=BEGNSQ;bi<XGENSZ;bi++) opened
            if((bi-BEGNSQ)%DGENSQ==0) sr=fscanf(fp," \n[%d]",&at);
            if(GRNMODE==NO) opened
               sr=fscanf(fp,"%d",&at); edp->popa[pi]->guys[gy].gen[bi]=at; 
            closed
            if(GRNMODE==YA) if((bi-BEGNSQ)%DGENSQ<SGENSQ) opened
               sr=fscanf(fp,"%d",&at); edp->popa[pi]->guys[gy].gen[bi]=at; 
            closed
         closed
      closed
   closed
closed
if(fp!=NULL) fclose (fp);

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Savememo(int thid) opened
char  *filename,cn[10],buf[40],*ac;
int   gi,pi,gy,ii,bi,hi,ri,ci,ce,lht,ilastall[10],at,bt,vx,vy,vz;
float af,bf,cf,flastall[10];
FILE  *fp; Dgx *dgx; Dgo *dgr;
Exbd  *edp=exbd[thid];

/*char  buf[80],str[80],tall[10][10]; time_t tneu;
struct tm *cur;
tneu = time(0);
cur = localtime(&tneu);
sprintf(str,"d");
if(cur->tm_mday<10) strcat(str,"0"); sprintf(cnb,"%dh",cur->tm_mday); strcat(str,cnb);
if(cur->tm_hour<10) strcat(str,"0"); sprintf(cnb,"%d:",cur->tm_hour); strcat(str,cnb);
if(cur->tm_min <10) strcat(str,"0"); sprintf(cnb,"%d ",cur->tm_min ); strcat(str,cnb);*/

fp=fopen(LASTAFN,"r"); if(fp!=NULL) opened
   ac=fgets(buf,40,fp);
   for(ii=0;ii<10;ii++) opened
      at= fscanf(fp,"%d %f",
      &ilastall[ii],&flastall[ii]);
      if((ilastall[ii]<0)||(ilastall[ii]>999999)) ilastall[ii]=-1;
      if((flastall[ii]<0)||(flastall[ii]>999999)) flastall[ii]=-1;
   closed
closed
fp=fopen(LASTAFN,"w"); if(fp!=NULL) opened
   bt=RNSEED%100-1; ilastall[bt]=edp->cn; 
   flastall[bt]=edp->popa[0]->guys[0].mfit;
   fprintf(fp,"cycle nr  best gaft\n");
   for(ii=0;ii<10;ii++) opened
      at=fprintf(fp,"%6d %10.3f\n",ilastall[ii], flastall[ii]); 
   closed
   fclose(fp);
closed

// xxxx
fp=fopen(HSTFN,"w"); if(fp!=NULL) opened
   for(ii=0;ii<FITHST;ii++) opened
      fprintf(fp," %4d",ii);
      for(gi=0;gi<NNGUYS;gi++) fprintf(fp," %6.3f",edp->fithst[ii][gi]);
      fprintf(fp,"\n");
   closed
   for(ii=0;ii<FITHST;ii++) opened
      fprintf(fp," %4d",ii);
      for(gi=0;gi<NNGUYS;gi++) fprintf(fp," %6.3f",edp->shdhst[ii][gi]);
      fprintf(fp,"\n");
   closed
closed
if(fp!=NULL) fclose(fp);

#if(TAGCHK==YA)
// TAG2FN
fp=fopen(TAG2FN,"w"); if(fp!=NULL) opened
   fprintf(fp," %4d]",-1);
   for(ri=0;ri<STAGES;ri++) fprintf(fp," %5d",ri);
   fprintf(fp," avg_se count\n");
   for(ri=0;ri<320;ri++) opened
      fprintf(fp," %4d)",ri);
      for(ci=0;ci<STAGES;ci++) fprintf(fp," %5.2f",edp->taghst[ri][ci]);
      af=edp->taghst[ri][STAGES+0];
      bf=edp->taghst[ri][STAGES+1];
      if(bf!=0) cf=af/bf; else cf=-1;
      fprintf(fp," %5.2f %5.2f\n",cf,bf);
   closed
closed
if(fp!=NULL) fclose(fp);

// TAG3FN
fp=fopen(TAG3FN,"w"); if(fp!=NULL) opened
   fprintf(fp,"%6s%3s%5s%5s","   nr);","timer;","actgn;","cpies;");
   fprintf(fp,"\n");
   for(ii=0;ii<edp->hgnr;ii++) opened
      fprintf(fp," %4d)",ii);
      fprintf(fp," %5d",edp->hstgenes[ii].dgx.timer);
      fprintf(fp," %5d",edp->hstgenes[ii].actgen);
      fprintf(fp," %5d",edp->hstgenes[ii].copies);
      for(hi=0;hi<MOCLEN;hi++) fprintf(fp," %5d",edp->hstgenes[ii].dgx.res[hi]);
      fprintf(fp,"\n");
   closed
closed
if(fp!=NULL) fclose(fp);
#endif

// Genome
lht = (int)strlen(OUTDIRN)+25; filename = (char *)malloc(lht);
strcpy(filename,OUTDIRN); strcat(filename,"bxgen_");
if(edp->cn<100000) strcat(filename,"0");
sprintf(cn, "%d", edp->cn); strcat(filename,cn); strcat(filename,".txt");
at=0;
//for(ii=1;ii<=10;ii++) if(abs(edp->cn-1000*AMSGECOE*5*ii)<=SHOWPACE) aiv=1;
if(edp->cn!=1)if((edp->cn-1)%125000==0) at=1;
if(at==0) fp=fopen(GENFN,"w"); if(at==1) fp=fopen(filename,"w");
if(fp!=NULL) opened
   fprintf(fp,"\nCYCLE %6d\n",edp->cn);
   for(pi=0;pi<NPOPUL;pi++) opened
      for(gy=0;gy<edp->popa[pi]->popsz;gy++) opened
         fprintf(fp,"%4d) %9.6f %9.6f ",
            gy,edp->popa[pi]->guys[gy].shad,edp->popa[pi]->guys[gy].gaft);
         // dev. genes
         for(bi=BEGNSQ;bi<XGENSZ;bi++) opened
            if((bi-BEGNSQ)%DGENSQ==0) fprintf(fp," \n[%3d]",(bi-BEGNSQ)/DGENSQ);
            if(GRNMODE==NO)
               fprintf(fp," %1d",(int)edp->popa[pi]->guys[gy].gen[bi]);
            if(GRNMODE==YA) if((bi-BEGNSQ)%DGENSQ<SGENSQ)
               fprintf(fp," %1d",(int)edp->popa[pi]->guys[gy].gen[bi]);
         closed
         fprintf(fp," \n");
      closed
   closed
closed
if(fp!=NULL) fclose(fp);
free(filename);

// dgar and dhar
for(ii=0;ii<2;ii++) opened
   if(ii==0) fp=fopen(DGARFN,"w"); if(ii==1) fp=fopen(DHLSFN,"w");
   fprintf(fp,"%5s","nr)");
   if(ii==0) opened
      fprintf(fp,"%4s%5s%4s","swc","ord","ds");
      for(hi=0;hi<4;hi++) fprintf(fp,"  hi%2d",hi);
      fprintf(fp,"%4s%4s","msa","msb");
      fprintf(fp,"%4s%4s%4s","0.x","0.y","0.z");
      fprintf(fp,"%4s%4s%4s","7.x","7.y","7.z");
      fprintf(fp,"%6s%6s%6s","tr0","tr1","tr2");
      fprintf(fp,"%6s%6s%6s","tr3","tr4","tr5");
      fprintf(fp,"%6s%6s%6s","tr6","tr7","tr8");
      fprintf(fp,"%5s%5s","conr","dpsx"); 
      fprintf(fp,"%5s%5s","dpsy","dpsz"); 
   closed
   if(ii==1) opened
      fprintf(fp,"  xx xx xx xx ais act");
      fprintf(fp,"%5s%5s%5s%5s%5s","d0","d1","d2","d3","d4");
      for(hi=0;hi<MOCLEN;hi++) fprintf(fp,"  hi%2d",hi); 
   closed
 //for(pi=0;pi<1;pi++) for(gy=0;gy<POPSZT;gy++) opened
   for(pi=0;pi<1;pi++) for(gy=0;gy<4;gy++) opened
      gi=pi*POPSZT+gy;
      fprintf(fp,"\n %4d ; %4d ; %4.2f",pi,gy,edp->popa[pi]->guys[gy].gaft);
      edp->popa[pi]->Gendecode(0,edp->popa[pi]->guys[gy].gen,1);
      fprintf(fp,"\n stepeval: %d\n",edp->stepeval);
      for(ce=0;ce<edp->dgarsz;ce++)
         edp->par[ce]=(float)edp->dgar[ce].arpos/edp->dgarsz;
      //Ftsortxdgx(0,edp->dgarsz,1,edp->dgar,edp->dgar);
      edp->dhnro=edp->dhnr0; edp->dharo=edp->dhar0;
      #if(MOCDEV==YA)
      edp->dhnro=edp->dhnral[gi]; edp->dharo=edp->dharal[gi];
      #endif
      if(ii==0) at=edp->dgarsz-1; if(ii==1) at=edp->dhnro;
      for(ce=0;ce<=at;ce++) opened
         fprintf(fp,"\n%4d)",ce);
         if(ii==0) opened
            dgx=&edp->dgar[ce]; dgr=&edp->dgar[ce].dgo;
            fprintf(fp," %3d %4d %3d",dgx->swtch,dgx->exord,dgx->timer);
            for(hi=0;hi<4;hi++) fprintf(fp," %5d",dgx->res[hi]);
            fprintf(fp," %3d %3d",dgr->msa,dgr->msb);
            for(hi=0;hi<=7;hi+=7) fprintf(fp," %3d %3d %3d",
               dgx->dgo.dpl[hi].x,dgx->dgo.dpl[hi].y,dgx->dgo.dpl[hi].z);
            for(hi=0;hi<9;hi++) fprintf(fp," %5.2f",dgx->dgo.tr[hi]);
            fprintf(fp," %4d %4d", dgx->dgo.conr, dgx->dgo.dpos.x);
            fprintf(fp," %4d %4d",dgx->dgo.dpos.y,dgx->dgo.dpos.z);
         closed
         if(ii==1) opened
            fprintf(fp,"%4d",edp->dharo[ce].cd.x);
            fprintf(fp,"%4d",edp->dharo[ce].cd.y);
            fprintf(fp,"%4d",edp->dharo[ce].cd.z);
            fprintf(fp,"%3d",edp->dharo[ce].clr);
            fprintf(fp,"%4d",edp->dharo[ce].lcrstp);
            fprintf(fp,"%4d",edp->dharo[ce].actstp);
            for(hi=0;hi<5;hi++) fprintf(fp," %4d",edp->dharo[ce].depend[hi]);
            for(hi=0;hi<MOCLEN;hi++) fprintf(fp," %5d",edp->dharo[ce].moc[hi]); 
         closed
      closed
   closed
   fclose(fp);
closed

ii=0;
fp=fopen(XDBG3FN,"w"); if(fp!=NULL) opened
   Forvxyz opened
      if(ii%10==0) at=fprintf(fp,"\n");
      at=fprintf(fp," %7d",edp->pmap[vx][vy][vz]); ii++;
   closed
closed
fclose(fp);

if(PAROPT!=0) opened
   fp=fopen(XDBG4FN,"w"); if(fp!=NULL) opened
      for(ii=0;ii<PARSETSZ;ii++) opened
         if(ii%10==0) at=fprintf(fp,"\n");
         at=fprintf(fp," %6.3f",edp->parset[ii]);
      closed
   closed
   fclose(fp);
closed

closed

#define SPEC_SHAPE0 0
#define SPEC_SHAPE1 0
//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Savegrid(int thid,int opt,int hi) opened
char  *filename,str[5]; int vx,vy,vz,ii,lht,ct,bb; FILE  *fp;
Exbd  *edp=exbd[thid];

lht = (int)strlen(GRDFILE)+14;
filename = (char *)malloc(lht);
strcpy(filename,GRDFILE);

if(opt==0) opened
   sprintf(str,"%d",hi);
   if(hi <1000) strcat(filename,"0");
   if(hi  <100) strcat(filename,"0");
   if(hi   <10) strcat(filename,"0");
   strcat(filename, str);
   strcat(filename,".txt");
closed

if(opt==1) opened
   bb=edp->pnr;
   //bb=(edp->pnr/(NSTEPS)+1)*NSTEPS-edp->pnr%NSTEPS;
   //if(edp->pnr%NSTEPS==0) bb=edp->pnr;
   if(bb <1000) strcat(filename,"0");
   if(bb  <100) strcat(filename,"0");
   if(bb   <10) strcat(filename,"0");
   sprintf(str,"%d",bb);
   strcat(filename, str);
   strcat(filename,".txt");
   printf("\n%4d %4d",edp->pnr,bb);
   edp->pnr++;
closed

fp=fopen(filename,"w");
fprintf(fp,"%3d %3d %3d\n",GRIDX,GRIDY,GRIDZ);
ii=0;
for(vz=0;vz<GRIDZ;vz++) for(vy=0;vy<GRIDY;vy++) for(vx=0;vx<GRIDX;vx++) opened
// for(vy=GRIDY-1;vy>=0;vy--)
   ct=edp->envr.sgrd[vx][vy][vz];
   if((ii==0)||((ii+0)%(GRIDX)==0))
      fprintf(fp,"\n (vz:%3d vy:%3d vx:%3d)",vz,vy,vx);
   fprintf(fp," %4d",ct); ii++;
closed

if(fp!=NULL) fclose(fp);
free(filename);

closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
void Body::Moveclsr(int thid,Clsr *csd,Clsr *css) opened
Exbd *edp=exbd[thid];

memcpy(csd,css,sizeof(Clsr)); memcpy(css,&edp->ics,sizeof(Clsr));

closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
void Body::Copyclsr(int thid,Clsr *csd,Clsr *css) opened

memcpy(csd,css,sizeof(Clsr));

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Tagchck2(int thid) opened
#if(TAGCHK==YA)
int   pi,gy,presacta,presactb,ga,gb,gc,ca,cb,cc;
float *par=NULL,af;
Exbd *edp=exbd[thid];

if ((edp->cn-1)%TAGPACE!=0) return;

for(gc=0;gc<NNGUYS;gc++) for(cc=0;cc<(NEVENT*STAGES);cc++) opened // only activated genes
   if (edp->frgenesc[gc][cc].timer!=-1) opened
      if(0==0) opened
         presactb=0;
         for(gb=0;gb<NNGUYS;gb++) for(cb=0;cb<(NEVENT*STAGES);cb++) opened
            if (edp->frgenesc[gc][cc].timer==edp->frgenesb[gb][cb].timer)
            if (memcmp(
               &edp->frgenesc[gc][cc].res[0],
               &edp->frgenesb[gb][cb].res[0],sizeof(int)*STAGES)==0) opened
               presactb=1; break;
            closed
         closed
      closed
      if(presactb==1) opened
         presacta=0;
         for(ga=0;ga<NNGUYS;ga++) for(ca=0;ca<(NEVENT*STAGES);ca++) opened
            if (edp->frgenesc[gc][cc].timer==edp->frgenesa[ga][ca].timer)
            if (memcmp(
               &edp->frgenesc[gc][cc].res[0],
               &edp->frgenesa[ga][ca].res[0],sizeof(int)*STAGES)==0) opened
               presacta=1; break;
            closed
         closed
      closed
      if(presacta==0) if(presactb==1) opened
         // presacta=0 and presactb=1: accepted gene
         af=1;
         pi=gc/POPSZT; gy=gc%POPSZT;
         af=edp->popa[pi]->guys[gy].gaft;
         edp->taghst[(edp->cn-1)/SHOWPACE][edp->frgenesc[gc][cc].timer]+=af;
         edp->taghst[(edp->cn-1)/SHOWPACE][STAGES+0]+=(float)edp->frgenesc[gc][cc].stepeval;
         edp->taghst[(edp->cn-1)/SHOWPACE][STAGES+1]+=1;
      closed
   closed
closed

#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Tagchck3(int thid) opened
#if(TAGCHK==YA)
int   pres,gc,cc,ii,at; float *par=NULL; Dgx *dgxp;
Exbd *edp=exbd[thid];

if ((edp->cn-1)%TAGPACE!=0) return;

// assigns std values to var excluded from comparison
for(gc=0;gc<NNGUYS;gc++) for(cc=0;cc<(NEVENT*STAGES);cc++) opened
   dgxp=&edp->frgenesc[gc][cc];
   for(ii=0;ii<3;ii++) dgxp->dhnrx[ii]=-2;
   dgxp->arpos=-2; dgxp->exord=-2; dgxp->dhptx=-2;
   dgxp->napos=-2; dgxp-> fsc =-2; dgxp->exeas=-2;
   dgxp->dgo.msa =-2; dgxp->dgo.c2nr=-2;
closed

// are old genes still present?
for(ii=0;ii<edp->hgnr;ii++) opened
   edp->hstgenes[ii].copies=0;
   pres=0;
   for(gc=0;gc<NNGUYS;gc++) for(cc=0;cc<(NEVENT*STAGES);cc++) opened
      if (memcmp(&edp->frgenesc[gc][cc],&edp->hstgenes[ii].dgx,sizeof(Dgx))==0) opened
         edp->hstgenes[ii].copies++; if (pres==0) edp->hstgenes[ii].actgen++;
         pres=1;
      closed
   closed
   if (pres==0) opened // position is freed
      //memcpy(&edp->hstgenes[ii].dgx,&edp->dgx0,sizeof(Dgx));
      edp->hstgenes[ii].actgen=-1;
      edp->hstgenes[ii].copies=-1;
   closed
closed

// are there new genes?
for(gc=0;gc<NNGUYS;gc++) for(cc=0;cc<(NEVENT*STAGES);cc++) opened
   if (edp->frgenesc[gc][cc].res[0]!=-1) opened
      pres=0; at=-1;
      for(ii=0;ii<edp->hgnr;ii++) opened
         if (edp->hstgenes[ii].actgen!=-1) // position occupied
         if (memcmp(&edp->frgenesc[gc][cc],&edp->hstgenes[ii].dgx,sizeof(Dgx))==0) opened
            pres=1; break; // if pres is 1 aiv is not used anyway
         closed
         if (edp->hstgenes[ii].actgen==-1) // position free
         if (at==-1) at=ii;              // assigns only once
      closed
      if(pres==0) opened
         if(at==-1) opened at=edp->hgnr; edp->hgnr++; closed
         memcpy(&edp->hstgenes[at].dgx,&edp->frgenesc[gc][cc],sizeof(Dgx));
         edp->hstgenes[at].actgen=0;
         edp->hstgenes[at].copies=0;
      closed
   closed
closed

// sort
Ftsortxdgf(thid,par,edp->hgnr,1,edp->hstgenes);
#endif

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Grncalc(int thid) opened closed

