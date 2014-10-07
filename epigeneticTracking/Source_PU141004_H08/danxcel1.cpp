
#include "danxincl.h"

// Global var
Exbd *exbd[NCORES]; char cvar;

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
int  main(int argc, char *argv[], char *env[]) opened
int  ii,opt,mrar[20],mrtot; Body *bdp;

#if(DEBUG0==YA)
double a=3.0,b=7.0,c;

set_fpu (0x27F);  // use double-precision rounding

c=a/b;
if (c==a/b) printf ("comparison succeeds\n");
else        printf ("unexpected result");
#endif

#if(DEBUG7==YA)
int ei,gu; double bdv;
for(ii=1;ii<=20;ii++) opened
   printf("\n");
   for(ei=0;ei<ii;ei++) opened
      bdv=NNGUYS/(double)ii; gu=(int)(bdv*(ei+1));
      if (ei==ii-1) gu=NNGUYS; printf(" %3d",gu);
   closed
closed
#endif

#if(DEBUG0==YA)
mrar[ 0]=sizeof( Exbd)              /(1024*1024); // Exbd
mrar[ 1]=sizeof( Guy )*NNGUYS       /(1024*1024); // popa
mrar[ 2]=sizeof( int )*SZGRID*3     /(1024*1024); // envr
mrtot=0; for(ii=0;ii<=2;ii++) mrtot+=mrar[ii];

printf("\n\nMemory check, alternative (MB)\n");
printf("sizeof Exbd  x core: %5d\n",mrar[ 0]);
printf("sizeof popa  x core: %5d\n",mrar[ 1]);
printf("sizeof envr  x core: %5d\n",mrar[ 2]);
printf("tot          x core: %5d\n",mrtot);

printf("\n\nMemory check, other\n");
printf("sizeof uns char      %5d\n",(int)sizeof(unsigned char));
printf("sizeof uns int       %5d\n",(int)sizeof(unsigned int));
printf("sizeof     int       %5d\n",(int)sizeof(int));
printf("sizeof   float       %5d\n",(int)sizeof(float));
#endif

// HINT UNIX
srand(RNSEED);

bdp =new Body();

if(argv[1]!=NULL) opt=0; else opt=1;
#if(RELOAD==YA)
opt=0;
#endif

bdp->Baseloop(opt,bdp);

return 0;
closed

//!--------------------------------------------------------------------------
//! FCT
//!--------------------------------------------------------------------------
Body::Body(void) {}

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Baseloop(int opt,Body *bdp) opened
char  *atmpstr,*atmpstr1,*buffer;
int   thid,gi,pi,gy,as,ls,us,ri,ei,ce,vx,vy,vz,ii,hi,gl,gu,ff,eval,di;
int   rn,aa,bb,size,ddone;
long  t0,t1,buflen,buflen0,buflen1,buflen3;
float shas; double gph,bdv; Coord cdl,cdu;
Guy *gyp,*gzp,*gtp; Exbd *edp;
#if(MYMPI==1)
int rc; MPI_Status stat;
#endif

#if(MYOMP==1)
omp_set_num_threads(NCORES);
#endif

#if(MYMPI==1)
rc=MPI_Init(NULL,NULL);
if (rc!=MPI_SUCCESS) opened
   printf ("Error MPI: Abort.\n"); MPI_Abort(MPI_COMM_WORLD, rc);
closed
else printf ("\nMPI_Init ok\n");
MPI_Comm_size (MPI_COMM_WORLD,&size); printf ("size: %d\n",size);
#endif

// mallocs
gtp=(Guy*)malloc(sizeof(Guy));
atmpstr =(char *)malloc(200);
atmpstr1=(char *)malloc(200);
Calcbuflen (&buflen,&buflen0,&buflen1,&buflen3);
buffer=(char*)malloc(buflen);

// calc all -----------------------------------------------------------------
#if(MYOMP==1)
#pragma omp parallel private(ei,thid,pi)
#endif
for(ei=0;ei<NCORES;ei++) opened
   thid=Getthid(&edp);
   if (thid==ei) opened
      exbd[ei] = (Exbd *)malloc(sizeof(Exbd));
      for(pi=0;pi<NPOPUL;pi++) exbd[ei]->popa[pi]=new Pop(ei,this);
      bdp->Initlocv(ei); bdp->Bodinit(ei);
      bdp->Envinit(ei,0); exbd[ei]->fdone=0;
   closed
   printf("\nThread %d ", thid);
closed
// calc all_end -------------------------------------------------------------

// calc master --------------------------------------------------------------
thid=Getthid(&edp);
if (thid==0) opened
   thid=Getthid(&edp);
   edp->cnta=0;
   //Txt2vtk(thid,1,46); exit(0);
   if(opt==0) bdp->Loadmemo(0); 
   if(opt==1) edp->cn0=0;
closed
// calc master_end ----------------------------------------------------------

// calc master --------------------------------------------------------------
remove(CONSFN);
thid=Getthid(&edp);
if (thid==0) opened
   strcpy(atmpstr1,CDTGFN);
   Loadtgt(0,atmpstr1);
   for(ii=0;ii<STAGES;ii++) opened 
      edp->evtnra0[ii]=0; edp->evtnrb0[ii]=0; 
   closed
closed
// calc master_end ----------------------------------------------------------

Sendrec0(0,buffer,&buflen0);

t0=(long)time(0);
thid=Getthid(&edp);
for(edp->cn=edp->cn0;edp->cn<=CYCLES;edp->cn++) opened

   // calc master -----------------------------------------------------------
   thid=Getthid(&edp);
   if (thid==0) opened
      edp->cfp=fopen(CONSFN, "a");
      edp->ausdr=0; if(((edp->cn-1)%SHOWPACE)==0) edp->ausdr=1;
      sprintf(atmpstr,"\n%d ",edp->cn);
      if(edp->ausdr==1) Printf3(edp->cfp,atmpstr); else fprintf(stdout,"X");
      edp->pnr=1;
   closed
   // calc master_end -------------------------------------------------------

   // calc master -----------------------------------------------------------
   // rnd permutation
   thid=Getthid(&edp);
   if (thid==0) opened
      for(gi=0;gi<NNGUYS;gi++) edp->rndpmt[gi]=gi;
      for(gi=1;gi<NNGUYS;gi++) opened
         rn=Rnd1(NNGUYS-1)+1; // 0 is excluded
         aa=edp->rndpmt[gi]; bb=edp->rndpmt[rn];
         edp->rndpmt[gi]=bb; edp->rndpmt[rn]=aa;
      closed
   closed
   // calc master_end -------------------------------------------------------

   GSLAVE:cvar=cvar; // slave starts here
   #if(MYOMP==1)
   thid=Getthid(&edp);
   if (thid!=0) for(ei=0;ei<NCORES;ei++) exbd[ei]->cn=exbd[0]->cn;
   #endif
   #if(MYMPI==1)
   thid=Getthid(&edp);
   if (thid!=0) edp->cn++;
   #endif

   // calc all --------------------------------------------------------------
   #if(MYOMP==1)
   #pragma omp parallel private(ei,thid,pi)
   #endif
   for(ei=0;ei<NCORES;ei++) opened
      thid=Getthid(&edp);
      if (thid==ei) opened
         // gq
         if((0<=edp->cn)&&(edp->cn<edp->frc[0].ge)) edp->gq=0;
         for(ii=0;ii<STAGES-1;ii++)
            if((edp->frc[ii].ge<=edp->cn)&&(edp->cn<edp->frc[ii+1].ge)) opened
               edp->gq=ii+1; break; 
            closed
         // boundaries of evolved sector, as number of dev. genes
         edp->dgarxf=edp->frc[edp->gq].sf;
         edp->dgarsz=edp->frc[edp->gq].se;
         edp->qf=Xqdet(thid,edp->dgarxf-0); 
         edp->qz=Xqdet(thid,edp->dgarsz-1);
         // boundaries of evolved sector, as number of bases
         for(pi=0;pi<NPOPUL;pi++) opened
            edp->popa[pi]->dsold=edp->popa[pi]->dnasz;
            edp->popa[pi]->dnaxf=(BEGNSQ+(DGENSQ*edp->dgarxf));
            edp->popa[pi]->dnasz=(BEGNSQ+(DGENSQ*edp->dgarsz));
         closed
      closed
   closed
   // calc all_end ----------------------------------------------------------

   // calc master -----------------------------------------------------------
   thid=Getthid(&edp);
   if (thid==0) opened
      edp->popa[0]->Prepgen (0);
      edp->popa[0]->Prexline(0);
      memcpy(&edp->guyf,&edp->popa[0]->guys[0],sizeof(Guy));
   closed

   // calc master -----------------------------------------------------------
   thid=Getthid(&edp);
   if (thid==0) opened
      sprintf(atmpstr,"\n[gq: %d ge: %5d sf: %3d se: %3d ls: %d us: %d]",
         edp->gq,edp->frc[edp->gq].ge,
         edp->frc[edp->gq].sf,edp->frc[edp->gq].se,
         edp->frz[edp->qf].ls,edp->frz[edp->qz].us);         
      if(edp->ausdr==1) Printf3(edp->cfp,atmpstr);
   closed

   // calc master -----------------------------------------------------------
   // swaps guys
   thid=Getthid(&edp);
   if (thid==0) for(aa=1;aa<NNGUYS;aa++) opened
      bb=edp->rndpmt[aa];
      gyp=&edp->popa[aa/POPSZT]->guys[aa%POPSZT];
      gzp=&edp->popa[bb/POPSZT]->guys[bb%POPSZT];
      memcpy(gtp,gyp,sizeof(Guy));
      memcpy(gyp,gzp,sizeof(Guy));
      memcpy(gzp,gtp,sizeof(Guy));
   closed
   // calc master_end ---------------------------------------------------------

   Sendrecv(1,buffer,&buflen1,&buflen3);

   // calc all ----------------------------------------------------------------
   #if(MYOMP==1)
   #pragma omp parallel private(ei,thid,pi)
   #endif
   for(ei=0;ei<NCORES;ei++) opened
      thid=Getthid(&edp);
      if (thid==ei) opened
         //Envinit(ei,1); // Inits env (actual grids only)
         edp->fset=0;
      closed
   closed
   // calc all_end ------------------------------------------------------------

   // computes extremes
   gl=0; gu=NNGUYS;
   #if(MYMPI==1)
   thid=Getthid(&edp);
   bdv=NNGUYS/(double)NCORES;
   gl=(int)(bdv*thid); gu=(int)(bdv*(thid+1));
   if (thid==NCORES-1) gu=NNGUYS;
   #endif

   #if(MYOMP==1)
   #pragma omp parallel private(gl,gu,gi,thid,edp,pi,gy,as,ls,us,ce,vx,vy,vz,hi,ff,eval,ddone)
   #endif
   if (0==0) opened
   gl=0; gu=NNGUYS;
   thid=Getthid(&edp);
   bdv=NNGUYS/(double)NCORES;
   gl=(int)(bdv*thid); gu=(int)(bdv*(thid+1));
   if (thid==NCORES-1) gu=NNGUYS; ff=0;
   for(gi=gl;gi<gu;gi++) opened
      edp->fset=0; pi=gi/POPSZT; gy=gi%POPSZT; ii=0;
      thid=Getthid(&edp);
      if(thid==0) if(edp->ausdr==1) if(gy==0) 
         printf("\nThread %d pi %d",thid,pi);
      ii=ii+1; if(ii%100==0) if(ii==0) opened printf(" "); ii=1; closed
      // Inits clar & dhar
      Forvxyz memcpy(&edp->clar[vx][vy][vz],&edp->ics,sizeof(Clsr));
      for(ce=0;ce<DHARLS;ce++) memcpy(&edp->dhar[ce],&edp->idse,sizeof(Dse));
      edp->dhnr=-1;
      // Decodes gen
      memcpy(&edp->guyx,&edp->popa[pi]->guys[gy],sizeof(Guy));
      if(ff==0) edp->popa[pi]->Gendecode(thid,edp->guyf.gen,0);
      else      edp->popa[pi]->Gendecode(thid,edp->guyx.gen,1);
      // dgar0
      if(ff==0) for(ce=0;ce<edp->dgarsz;ce++)
          memcpy(&edp->dgar0[ce],&edp->dgar[ce],sizeof(Dgx));
      // inits dvfrd
      for(ce=0;ce<edp->dgarsz;ce++) edp->dgar[ce].dvfrd=0;
      // shad & mfit
      edp->popa[pi]->guys[gy].shad=0;
      edp->popa[pi]->guys[gy].mfit=0; 
      // inits mothers
      for(hi=0;hi<1000;hi++) opened 
         edp->mothers[hi][0]=-1; edp->mothers[hi][1]=-1; 
      closed
      edp->motnr=0;
      // development cycle
      ls=0; us=STAGES-1; /*us=CLKMAX-1;*/
      if(edp->fdone!=0) opened 
         ls=edp->frz[edp->qf].ls; //ls=edp->frz[edp->gq].ls;  
         us=edp->frz[edp->qz].us; //us=edp->frz[edp->gq].us;
      closed
      for(as=ls;as<=us;as++) opened
         // Inits zygotes and isgrid
         Gridinit(thid,gi,as);
         if ((edp->fdone==0)&&(as==edp->frz[edp->qf].ls)) opened
            Forvxyz memcpy(&edp->clarf[vx][vy][vz],&edp->clar [vx][vy][vz],sizeof(Clsr));
            memcpy(&edp->dharf[0],&edp->dhar[0],sizeof(Dse)*DHARLS);
            edp->dhnrf=edp->dhnr;
            edp->fset=1;
         closed
         if ((edp->fdone!=0)&&(as==ls)) opened
            Forvxyz memcpy(&edp->clar [vx][vy][vz],&edp->clarf[vx][vy][vz],sizeof(Clsr));
            memcpy(&edp->dhar[0],&edp->dharf[0],sizeof(Dse)*DHARLS);
            edp->dhnr=edp->dhnrf;
         closed
         //Disrupt (thid,gi,as);
         Dgarprep(thid,gi,as,ls);
         #if(GRNMODE==YA)
         Grncalc(thid);
         #endif
         Shaper  (thid,gi,as,us);
         //Doper2 (if gi==0 or it has other as ahead or is about to change ms)
         if((as!=us)||(DOPAON==YA)) opened
            if(NOSCCH==0) Doper2(thid,gi,as);
         closed
         if ((edp->cn<edp->frc[edp->gq].ge)&&((edp->cn+1)>=edp->frc[edp->gq].ge))
         if ((gi==0)&&(as==ls)) printf(" dop2");
         // Doper2_end
         // fitness evaluation
         if((PGFREEZE==YA)||(PGFTEVAL==YA)) opened
            eval=0;
            if((as==ESTAGE)||((as==us)&&(us<ESTAGE))) eval=1;
         closed
         if(EVFTEVAL==YA) opened
            eval=0; if(as==edp->stepeval) eval=1;
         closed
         #if(CFLAG7==YA)
         if(eval==1) opened
            edp->ext[0].x=0; edp->ext[1].x=GRIDX-1;
            edp->ext[0].y=0; edp->ext[1].y=GRIDY-1;
            edp->ext[0].z=0; edp->ext[1].z=GRIDZ-1;
            Smooth3d(thid,gi,&edp->ext[0],&edp->ext[1],1);
         closed
         #endif
         // fitness evaluation
         if(eval==1) opened
            cdl.x=0; cdl.y=0; cdl.z=0; cdu.x=GRIDX-1; cdu.y=GRIDY-1; cdu.z=GRIDZ-1;
            edp->fshad=(float)Fitness(thid,gi,edp->clar);
            edp->popa[pi]->guys[gy].shad=edp->fshad;
         closed
         Gridupd(thid,gi,as,us,ff);
         // for guy 0, does some operations
         if (ff==0) if (as==us) opened
            memcpy(&edp->dhar0,&edp->dhar,sizeof(Dse)*DHARLS);
            for(ce=0;ce<DHARLS;ce++) opened
               for(hi=0;hi<MOCLEN;hi++) edp->dharaz[0][ce].moc[hi]=edp->dhar[ce].moc[hi];
               edp->dharaz[0][ce].lcrstp=edp->dhar[ce].lcrstp; 
            closed
            edp->dhnr0=edp->dhnr;
            for(ri=0;ri<DHARLS;ri++)
               memcpy(&edp->drvar0[ri],&edp->drvaro[ri],sizeof(Coord));
            edp->stepeval0=edp->stepeval;
            Forvxyz memcpy(&edp->clar0[vx][vy][vz],&edp->clar[vx][vy][vz],sizeof(Clsr));
         closed
      closed
      #if(MOCDEV==YA)
      memcpy(&edp->dharal[gi][0],&edp->dhar[0],sizeof(Dse)*DHARLS);
      /*for(ce=0;ce<DHARLS;ce++)
         if((ce<=edp->dhnr)&&(edp->dhar[ce].asval<=edp->frz[edp->gq].fs))
            edp->dhnral[gi]=ce;*/
      edp->dhnral[gi]=edp->dhnr;
      #endif
      if (ff==0) gi--;
      if((ff==0)&&(edp->fset==1)) edp->fdone=1; ff++;
      if (gi==gu-1)
      if((edp->cn+0)< edp->frc[edp->gq].ge)
      if((edp->cn+1)>=edp->frc[edp->gq].ge) edp->fdone=0;
   closed
   closed

   thid=Getthid(&edp);
   Statevents (thid,atmpstr,atmpstr1);
   Sendrecv(3,buffer,&buflen1,&buflen3);

   // calc master -------------------------------------------------------------
   // restores order guys
   thid=Getthid(&edp);
   if (thid==0) for(aa=NNGUYS-1;aa>=1;aa--) opened
      bb=edp->rndpmt[aa];
      //printf(" aa=%d",aa); printf(" bb=%d",bb);
      gyp=&edp->popa[aa/POPSZT]->guys[aa%POPSZT];
      gzp=&edp->popa[bb/POPSZT]->guys[bb%POPSZT];
      memcpy(gtp,gyp,sizeof(Guy));
      memcpy(gyp,gzp,sizeof(Guy));
      memcpy(gzp,gtp,sizeof(Guy));
   closed
   // restores mocdev
   #if(MOCDEV==YA)
   thid=Getthid(&edp);
   if (thid==0) opened
      for(aa=NNGUYS-1;aa>=1;aa--) opened
         bb=edp->rndpmt[aa];
         memcpy(&edp->dharaltemp[0],&edp->dharal[aa][0],sizeof(Dse)*DHARLS);
         memcpy(&edp->dharal[aa][0],&edp->dharal[bb][0],sizeof(Dse)*DHARLS);
         memcpy(&edp->dharal[bb][0],&edp->dharaltemp[0],sizeof(Dse)*DHARLS);
         edp->dhnraltemp=edp->dhnral[aa];
         edp->dhnral[aa]=edp->dhnral[bb];
         edp->dhnral[bb]=edp->dhnraltemp;
      closed
   closed
   #endif
   // calc master_end -------------------------------------------------------

   #if(MYMPI==1)
   thid=Getthid(&edp);
   if((thid!=0)&&((edp->cn)<CYCLES)) goto GSLAVE;
   #endif

   // calc master -------------------------------------------------------------
   thid=Getthid(&edp);
   if (thid==0) opened
      Calcfit(0);
      if(edp->ausdr==YA) opened
         Showres(0,edp->cfp,atmpstr); Savememo(0); 
      closed
      if(edp->cn%GPSTEP==0) edp->popa[0]->Germline(0);
      for(pi=0;pi<NPOPUL;pi++) edp->popa[pi]->Galoop(0,0,edp->cfp);
      t1=(long)time(0); gph=(double)((edp->cn-edp->cn0+1)*3600)/(t1-t0);
      sprintf(atmpstr,"gph: %.0f \n",gph);
      if(edp->ausdr==YA) Printf3(edp->cfp,atmpstr);
      fclose(edp->cfp);
   closed
   if(REPLAY==YA) exit(0);
closed

#if(MYMPI==1)
MPI_Finalize();
#endif

free(gtp); free(atmpstr); free(atmpstr1);
closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Calcbuflen (long *blp,long *blp0,long *blp1,long *blp3) opened
long  buflen,buflen0,buflen1,buflen3;

// sendrec_0
buflen0 =sizeof(int);
buflen0+=sizeof(int)*SZGRID; 
buflen0+=2*sizeof(int)*SZGRID;
// sendrec_1
buflen1 =sizeof(int)*6;
buflen1+=sizeof(int)*7;
buflen1+=sizeof(Guy);
buflen1+=sizeof(int)*NNGUYS;
buflen1+=sizeof_gen_*(NNGUYS/NCORES+1);
buflen1+=sizeof(bas)*BEGNSQ*(NNGUYS/NCORES+1);
// sendrec_3
buflen3 =sizeof(int);
buflen3+=sizeof(int);
buflen3+=sizeof(float)*(NNGUYS/NCORES+1);
buflen3+=sizeof(float)*(NNGUYS/NCORES+1);
#if(MOCDEV==YA)
buflen3+=sizeof(Dse)*(NNGUYS/NCORES+1)*DHARLS;
buflen3+=sizeof(Dgx)*(NNGUYS/NCORES+1);
#endif
// calc
buflen=(buflen1>=buflen0) ? buflen1 :buflen0;
buflen=(buflen3>=buflen ) ? buflen3 :buflen ;
printf("%ld %ld %ld %ld",buflen0,buflen1,buflen3,buflen);
*blp=buflen; *blp0=buflen0; *blp1=buflen1; *blp3=buflen3; 

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Statevents (int thid,char *atmpstr,char *atmpstr1) opened
int ii,ri,si,at,bt,ct;
Exbd *edp=exbd[thid];
   
// calc master -------------------------------------------------------------
thid=Getthid(&edp);
if (thid==0) opened
   if(edp->ausdr==1) opened
      sprintf(atmpstr,"\n\nstepeval0: %2d\n ",edp->stepeval0); 
      Printf3(edp->cfp,atmpstr);
      for(ii=0;ii<STAGES;ii++) opened
         sprintf(atmpstr," %2d",ii); 
         if((ii+1)%25==0) atmpstr=strcat(atmpstr,"\n ");
         Printf3(edp->cfp,atmpstr); 
      closed
      at=0;
      sprintf(atmpstr,"\na"); Printf3(edp->cfp,atmpstr);
      for(ii=0;ii<STAGES;ii++) opened
         sprintf(atmpstr," %2d",edp->evtnra0[ii]); at+=edp->evtnra0[ii];
         if((ii+1)%25==0) atmpstr=strcat(atmpstr,"\na"); Printf3(edp->cfp,atmpstr);
      closed
      sprintf(atmpstr,"\nb"); Printf3(edp->cfp,atmpstr);
      for(ii=0;ii<STAGES;ii++) opened
         sprintf(atmpstr," %2d",edp->evtnrb0[ii]); at+=edp->evtnrb0[ii];
         if((ii+1)%25==0) atmpstr=strcat(atmpstr,"\nb"); Printf3(edp->cfp,atmpstr);
      closed
      sprintf(atmpstr," %2d\n",at);
      Printf3(edp->cfp,atmpstr);
      // statistics double mocs
      bt=0; ct=-1;
      for(ri=0;ri<edp->dhnr0;ri++) if (edp->dhar0[ri].moc[0]!=-1) opened
         memcpy(&edp->moctmp[0],&edp->dhar0[ri].moc[0],sizeof(int)*MOCLEN);
         at=0;
         for(si=0;si<edp->dhnr0;si++) if(memcmp(&edp->moctmp[0],
            &edp->dhar0[si].moc[0],sizeof(int)*MOCLEN)==0) at++;
         cvar=cvar; // anchor
         if (at>bt) opened bt=at; ct=ri; closed
      closed
      sprintf(atmpstr,"max nr of double moc: %2d %4d",bt,ct);
      Printf3(edp->cfp,atmpstr);
      // statistics double mocs_end
   closed
closed
// calc master_end ---------------------------------------------------------

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Initlocv (int thid) opened
int  ii,hi;
Exbd *edp=exbd[thid];

// Inits
for(hi=0;hi<MOCLEN;hi++) edp->imoc[hi]=-1;

// Inits
edp->icd.x=-1; edp->icd.y=-1; edp->icd.z=-1;

// Inits
edp->idgo.msa=-1; edp->idgo.msb=-1;
edp->idgo.conr=-1; edp->idgo.c2nr=-1;
for(ii=0;ii<2;ii++) opened
   edp->idgo.dpl[ii].x=-1;
   edp->idgo.dpl[ii].y=-1;
   edp->idgo.dpl[ii].z=-1;
closed

// Inits
edp->idgx.exord=-1; edp->idgx.swtch=-1; edp->idgx.timer=-1;
edp->idgx.arpos=-1; edp->idgx.dvfrd=-1;
for(ii=0;ii<3;ii++) edp->idgx.dhnrx[ii]=-1;
for(hi=0;hi<MOCLEN;hi++) edp->idgx.res[hi]=-1;
memcpy(&edp->idgx.dgo,&edp->idgo,sizeof(Dgo));

// Inits
edp->iguy.neugy=YA; edp->iguy.gaft=0; edp->iguy.shad=-1;
memset(&edp->iguy.gen,0,sizeof(bas)*XGENSZ);

// Inits
Lwriter(Lxxx,&edp->ics,-1,0,CLCHOLE,0);

// Inits
for(hi=0;hi<MOCLEN;hi++) edp->idse.moc[hi]=-1;
for(ii=0;ii<50;ii++) edp->idse.depend[ii]=-1; edp->idse.ndep=0;
for(ii=0;ii<NMORPH;ii++) edp->idse.infmoth[ii]=-1;
edp->idse.lcrstp=-1; edp->idse.actstp=-1; edp->idse.sigstp=-1;
edp->idse.dltcrt=0; edp->idse.actdrv=0;

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void Body::Envinit(int thid,int opt) opened
int  vx,vy,vz,ii,it,xi,yi,zi,dist,distmin,aa,bb,cnt; Coord cdclose[2],cdmin;
int  ndqblc[3]= opened NDQBLC closed;
Exbd *edp=exbd[thid];

if(opt==0) opened
   if(0==0) opened
      edp->envr.sgrd = new dimensa[GRIDX];
      edp->envr.cdtg = new dimensa[GRIDX];
      edp->envr.etis = new dimensa[GRIDX];
      edp->envr.cdis = new dimensa[GRIDX];
   closed
closed

// Inits
if(opt==0) Forvxyz opened
   if(0==0) opened
      edp->envr.sgrd[vx][vy][vz]=-1;
      edp->envr.cdtg[vx][vy][vz]=0;
      edp->envr.etis[vx][vy][vz]=0;
      edp->envr.cdis[vx][vy][vz]=0;
   closed
closed

// cdis
if(opt==0) opened
   Forvxyz edp->envr.cdis[vx][vy][vz]=-1;
   Forvxyz opened 
      if(vx%SMOOTHSZ==0)if(vy%SMOOTHSZ==0)if(vz%SMOOTHSZ==0)
         edp->envr.cdis[vx][vy][vz]=COLRNG/2-RNDRNG/2+Rnd1(RNDRNG);
       //edp->envr.cdis[vx][vy][vz]=Rnd1(RNDRNG);
   closed
   // determines dominating colour
   Forvxyz opened
      distmin=999999;
      cdmin.x=-1; cdmin.y=-1; cdmin.z=-1;
      cdclose[0].x=(vx/SMOOTHSZ+0)*SMOOTHSZ; cdclose[1].x=(vx/SMOOTHSZ+1)*SMOOTHSZ;
      cdclose[0].y=(vy/SMOOTHSZ+0)*SMOOTHSZ; cdclose[1].y=(vy/SMOOTHSZ+1)*SMOOTHSZ;
      cdclose[0].z=(vz/SMOOTHSZ+0)*SMOOTHSZ; cdclose[1].z=(vz/SMOOTHSZ+1)*SMOOTHSZ;
      if((cdclose[1].x<0)||(cdclose[1].x>=GRIDX)) cdclose[1].x=cdclose[0].x;
      if((cdclose[1].y<0)||(cdclose[1].y>=GRIDY)) cdclose[1].y=cdclose[0].y;
      if((cdclose[1].z<0)||(cdclose[1].z>=GRIDZ)) cdclose[1].z=cdclose[0].z;
      for(xi=0;xi<=1;xi++) for(yi=0;yi<=1;yi++) for(zi=0;zi<=1;zi++) opened
         dist=
         abs(vx-cdclose[xi].x)+abs(vy-cdclose[yi].y)+abs(vz-cdclose[zi].z);
         if(distmin>dist) opened
            distmin=dist;
            cdmin.x=cdclose[xi].x;
            cdmin.y=cdclose[yi].y;
            cdmin.z=cdclose[zi].z;
         closed
      closed
      cvar=cvar;
      edp->envr.cdis[vx][vy][vz]=edp->envr.cdis[cdmin.x][cdmin.y][cdmin.z];
   closed
closed

// etis
if(opt==0) opened
   // writes drv signals to one quadrant
   Forvxyz edp->envr.etis[vx][vy][vz]=0;
   Forvxyz
      if((vz<=GRIDZ/2)&&(vy<=GRIDY/2)&&(vz<=GRIDZ/2))
      if((vx%ndqblc[0]==0)&&(vy%ndqblc[1]==0)&&(vz%ndqblc[2]==0)) 
         edp->envr.etis[vx][vy][vz]=1;
   // copies symmetrically to other quadrants
   Forvxyz if(vx>GRIDX/2) opened
      aa=GRIDX-1-vx;
      edp->envr.etis[vx][vy][vz]=edp->envr.etis[aa][vy][vz]; 
   closed
   Forvxyz if(vy>GRIDY/2) opened
      aa=GRIDY-1-vy;
      edp->envr.etis[vx][vy][vz]=edp->envr.etis[vx][aa][vz]; 
   closed
   Forvxyz if(vz>GRIDZ/2) opened
      aa=GRIDZ-1-vz;
      edp->envr.etis[vx][vy][vz]=edp->envr.etis[vx][vy][aa]; 
   closed
closed

closed

//!--------------------------------------------------------------------------
//! fctheader
//!--------------------------------------------------------------------------
void  Body::Bodinit(int ei) opened
int   vx,vy,vz,cx,ii,gi,di,ri,ci,fi,ti,val,aa,bb,min,max; float cf;
Exbd  *edp=exbd[ei];

// grn
edp->sgarsz=DGARSZ;

for(ii=0;ii<STAGES;ii++) opened
   fi=ii-1; cf=AMSCOE;
   if(AMSLLT==0) ti=ii; if(AMSLLT==1) ti=1;
   edp->frc[ii].gf =(int)(fi*1000*cf); 
   edp->frc[ii].ge =(int)(ii*1000*cf);
   edp->frc[ii].sf =(int)(fi*AMSDGN);  
   edp->frc[ii].se =(int)(ii*AMSDGN+AMSDGO);
   // xxxx
   edp->frz[ii].xe =(int)(ii*AMSDGN);  
   edp->frz[ii].xf =(int)(fi*AMSDGN); 
   edp->frz[ii].ls =(int)(ti);         
   edp->frz[ii].us =(int)(ti);
   edp->frz[ii].dl =(int)(AMSDLZ);
   // corrections
   min=0; max=DGARSZ;
   edp->frc[0].gf=0; edp->frc[0].ge=0;
   edp->frc[0].sf=0; edp->frc[0].se=0;
   edp->frz[0].xf=0; edp->frz[0].xe=0;
   Bound(&edp->frc[ii].sf,0,DGARSZ);
   Bound(&edp->frc[ii].se,0,DGARSZ);
   Bound(&edp->frz[ii].xf,0,DGARSZ);
   Bound(&edp->frz[ii].xe,0,DGARSZ);
closed

#if(DEBUG0==YA)
if(ei==0) opened
   printf("\n");
   for(ii=0;ii<STAGES;ii++) opened
      printf("%3d",edp->frz[ii].dl); if((ii+1)%25==0) printf("\n");
   closed
   printf("\n");
closed
#endif

// xqar
for(di=0;di<DGARSZ;di++) opened
   val=di;
   for(ii=-1;ii<STAGES-1;ii++) opened
      if(ii==-1) opened aa=0; bb=edp->frz[0].xe; closed
      else opened aa=edp->frz[ii].xe; bb=edp->frz[ii+1].xe; closed
      if(((aa<=val)&&(val<bb))) opened edp->xqar[di]=ii+1; break; closed
   closed
closed

// dgarsz,cn,cn0,clarf
edp->dgarsz=DGARSZ; edp->cn=0; edp->cn0=0;
Forvxyz memcpy(&edp->clarf[vx][vy][vz],&edp->ics,sizeof(Clsr));

// fithst,zyg
int zygsss[12]=ZYGSSS;
for(ii=0;ii<FITHST;ii++) for(gi=0;gi<NNGUYS;gi++) opened
   edp->shdhst[ii][gi]=0; edp->fithst[ii][gi]=0;
closed

for(ii=0;ii<ZYGTOT;ii++) opened
   edp->zygc[ii].x=zygsss[3*ii+0];
   edp->zygc[ii].y=zygsss[3*ii+1];
   edp->zygc[ii].z=zygsss[3*ii+2];
closed

// Exdp various
Forvxyz memcpy(&edp->clar[vx][vy][vz],&edp->ics,sizeof(Clsr));
//Forvxyz Clsrinit(&edp->claro[vx][vy][vz],CLNOCEL);
for(cx=0;cx<DHARLS;cx++) memcpy(&edp->clard[cx],&edp->ics,sizeof(Clsr));
memcpy(&edp->moclp,&edp->ics,sizeof(Clsr));
for(di=0;di<edp->dgarsz;di++) opened
   memcpy(&edp->dgar  [di],&edp->idgx,sizeof(Dgx));
   memcpy(&edp->dgare [di],&edp->idgx,sizeof(Dgx));
   memcpy(&edp->dgaro [di],&edp->idgx,sizeof(Dgx));
closed

// rrorder was here

#if(CFLAG0==YA)
Setbitmask();
#endif

closed

//!----------------------------------------------------------------------------
//! fctheader
//!----------------------------------------------------------------------------
int Body::Xqdet(int thid,int val) opened
int ii,aa,bb,xq;
Exbd  *edp=exbd[thid];

for(ii=-1;ii<STAGES-1;ii++) opened
   if(ii==-1) opened aa=0; bb=edp->frz[0].xe; closed
   else opened aa=edp->frz[ii].xe; bb=edp->frz[ii+1].xe; closed
   if(((aa<=val)&&(val<bb))) opened xq=ii+1; break; closed
closed

return xq;
closed
