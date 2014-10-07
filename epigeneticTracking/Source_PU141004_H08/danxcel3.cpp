
#include "danxincl.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cstdlib>

using namespace std;

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
double Body::Fitness (int thid,int gi,Clsr clar[][GRIDY][GRIDZ]) opened
/* VOXCAD FITNESS */
/**/
bool debug=false;

if (debug) cout << "Entering evaluation " << GRIDX << " " << GRIDY << " " << GRIDZ << endl;


double fitness=0;
int num_x_voxels=10;
int num_y_voxels=10;
int num_z_voxels=10;
double voxelSize=0.01;

// number of ETcells per voxelyze voxel
int xratio=GRIDX/num_x_voxels;
int yratio=GRIDY/num_y_voxels;
int zratio=GRIDZ/num_z_voxels;

vector< vector< vector< int > > > matrixForVoxelyze (num_x_voxels, std::vector<std::vector<int> >(num_y_voxels, std::vector<int>(num_z_voxels, 0)));

// Number of cells per type (0 is no cell)
int ntype[5];
int nCellsPerType[]={0,0,0,0,0};

if (debug) cout << "Building voxelyze matrix...";
// Subdividing the matrix
for (int vx=1; vx<num_x_voxels-1; vx++) opened
   for (int vy=1; vy<num_y_voxels-1; vy++) opened
      for (int vz=1; vz<num_z_voxels-1; vz++) opened
	 // Getting the number of cells of each type
	 for (int i=0; i<5; i++) opened
	    ntype[i]=0;
	 closed
	 for (int vxi=0; vxi<xratio; vxi++) opened
	    for (int vyi=0; vyi<yratio; vyi++) opened
	       for (int vzi=0; vzi<zratio; vzi++) opened
 		  //cout << clar[vx*xratio+vxi][vy*yratio+vyi][vz*zratio+vzi].dhc << " "; 
                  switch (clar[vx*xratio+vxi][vy*yratio+vyi][vz*zratio+vzi].dhc) opened
		     case 0:  ntype[1]++; break;
		     case 1:  ntype[2]++; break;
                     case 2:  ntype[3]++; break;
                     case 3:  ntype[4]++; break;
		     default: ntype[0]++; break; // unknown type => no cell
                  closed
               closed
            closed
         closed
if (debug)         cout << "[" << ntype[0] << ",";
	 // looking for maximum
         int maxind=0;
	 for (int i=1; i<5; i++) opened
	    if (ntype[i]>ntype[maxind]) opened
	       maxind=i;
	    closed
if (debug)            cout << ntype[i] << ",";
         closed
if (debug)         cout << "] ";
	 // writing in the voxelyze matrix
if (debug) 	 cout << maxind << " ";
         matrixForVoxelyze[vx][vy][vz]=maxind;
	 nCellsPerType[maxind]++;
      closed
   closed
if (debug)    cout << endl;
closed

cout << "[";
for (int i=0; i<4; i++) opened
   cout << nCellsPerType[i] << ",";
closed
cout << nCellsPerType[4] << "]  ";

if (debug) cout << "done" << endl;
// writing voxelyze file
if (debug) cout << "Writing voxelyze file...";
writeVoxelyzeFile(matrixForVoxelyze);
if (debug) cout << "done" << endl;

                // Evaluate Individual =====================================================================================================================================
if (debug) cout << "Running voxelyze...";
        // just for timekeeping
                clock_t start;
                clock_t end;
                start = clock();
                // cout << "starting voxelyze evaluation now" << endl;

                // System call to voxelyze to run simulation -- make sure that executable is already built, and that it lives in the path specified here
                FILE* input = popen("./voxelyze -f voxelyzeInputFromCPPN.vxa","r");
                // FILE* input = popen(callVoxleyzeCmd.str().c_str(),"r");
if (debug) cout << "done" << endl;
if (debug) cout << "Waiting for voxelyze...";

                // Read Fitness File ===========================================================================================================================================

                // continually check for return of fitness file (non-optimized, for a cleaner approach pass with sockets)int exitCode0 = std::system("mkdir champVXAs");
                bool doneEval = false;
                while (not doneEval)
                {
                        end = clock();
                        std::ifstream infile("softbotsOutput.xml"); // this is the name of the fitness file created (as we specified it when writing the vxa)
                        // std::ifstream infile(outFileName.str().c_str());

                        // if file exists, note how long it took to complete, and exit this loop by switching flag
                        if (infile.is_open())
                        {
                                printf("voxelyze took %.6lf seconds\n", float(end-start)/CLOCKS_PER_SEC);
                                //usleep(1);
                                doneEval = true;
if (debug) 				cout << "done" << endl;
                        }
                        // if the file doesn't exist, wait for a small period of time and check for it again
                        else
                        {
                                //usleep(1);
                                // if the file hasn't been found in a long period of time, voxelyze may have become unstable and crashed, kill the simulations and assign minimum fitness so the whole program doesn't crash/hang.  
                                if ( double(end-start)/CLOCKS_PER_SEC > 120.0) // amount of time set arbitrarily.  For a more scalable value, make function of the number of voxels
                                {
                                        cout << "voxelyze hung after 120 seconds... assigning fitness of 0.000001"<<endl;

                                        // find process and kill it.  Please optimize this better.
                                        int exitCode3 = std::system("ps axu > /tmp/HnPsFile.txt"); // write system processes to file
                                        std::ifstream psfile("/tmp/HnPsFile.txt"); // read them from that file
                                        std::string thisLine;
                                        if (psfile.is_open())
                                        {
                                                while (std::getline(psfile, thisLine))
                                                {
                                                        std::size_t foundvox = thisLine.find("./voxelyze"); // if a thread is running voxelyze
                                                        if (foundvox!=std::string::npos)
                                                        {
                                                                if (atoi(thisLine.substr(foundvox-8,4).c_str()) >= 2) // if it has been running longer that 2 minutes
                                                                {
                                                                        std::size_t foundsp = thisLine.find(" ");
                                                                        cout << ("kill "+thisLine.substr(foundsp,11)).c_str() << endl;
                                                                        int exitCode4 = std::system(("kill "+thisLine.substr(foundsp,11)).c_str()); // kill it
                                                                        int exitCode5 = std::system("killall <defunct>"); // kill any other defunct processes while we are here (as this crude clean up method may cause them)
                                                                }
                                                        }
                                                }
                                        }

                                        cout << endl;
                                        return 0.000001; // if we kill the simulation, just return the minimum fitness
                                }
                        }
                        infile.close();
                }

                pclose(input); // clean up the process running the simulation

                // open the fitness file again (you could also just move your pointer back to the beginning of the file)
                std::ifstream infile("softbotsOutput.xml"); // this is the name of the fitness file created (as we specified it when writing the vxa)
                // std::ifstream infile(outFileName.str().c_str());

                std::string line;
                float FinalCOM_DistX;
                float FinalCOM_DistY;
                float FinalCOM_DistZ;

                // parse the file line by line to find the fitness values you are looking for (in this case diatance moved).  this could also be done by parsing the xml tree
                if (infile.is_open())
                {
                        while (std::getline(infile, line))
                        {
                            std::size_t foundx = line.find("<FinalCOM_DistX>");
                            if (foundx!=std::string::npos)
                            {
                                FinalCOM_DistX = atof(line.substr(foundx+strlen("<FinalCOM_DistX>"),line.find("</")-(foundx+strlen("<FinalCOM_DistX>"))).c_str());
                            }
                            std::size_t foundy = line.find("<FinalCOM_DistY>");
                            if (foundy!=std::string::npos)
                            {
                                FinalCOM_DistY = atof(line.substr(foundy+strlen("<FinalCOM_DistY>"),line.find("</")-(foundy+strlen("<FinalCOM_DistY>"))).c_str());
                            }
                            std::size_t foundz = line.find("<FinalCOM_DistZ>");
                            if (foundz!=std::string::npos)
                            {
                                FinalCOM_DistZ = atof(line.substr(foundz+strlen("<FinalCOM_DistZ>"),line.find("</")-(foundz+strlen("<FinalCOM_DistZ>"))).c_str());
                            }
                        }

                        fitness = pow(pow(FinalCOM_DistX,2)+pow(FinalCOM_DistY,2),0.5);  // in this example we are only looking for displacement in the X-Y plane

                        // normalize the fitness from an absolute distance measure to display in the unity of "body lengths" (the largest dimension of the array holding the robot)
                        fitness = fitness/(max(num_x_voxels,max(num_y_voxels,num_z_voxels))*voxelSize);

                        // cout << "Original Fitness from voxelyze: " << fitness << endl;

                        infile.close();
                }

                // if fitness is not greater than zero or is an absurdly large number, something probably went wrong.  Just assign minimum fitness and exit the evaluation.
                // SCB if (fitness < 0.000001) return 0.000001;
                // SCB if (fitness > 10000) return 0.000001;

                // adjust fitness by the penalty factor
                // SCB fitness = fitness * calculateFitnessAdjustment( matrixForVoxelyze );

                if (fitness < 0.000001) fitness = 0.000001;
                if (fitness > 10000) fitness = 0.000001;
		if (fitness != fitness) fitness = 0.000001;

		cout << "Fitness: " << fitness << endl;
                return fitness;


/**/
   


/*
// ALEX'S FITNESS 

int   pi,gy,lx,ly,lz,ux,uy,uz,vx,vy,vz,ci,ides,iins,at,bt;
float shad,bdes,bins,bous,cins,xins,cinsar[COLRNG],shco0,oc,ious,af,bf;
Clsr  *clp,*clf; Exbd *edp=exbd[thid];

pi=gi/POPSZT; gy=gi%POPSZT;

ides=0; iins=0; ious=0; for(ci=0;ci<COLRNG;ci++) cinsar[ci]=0;
Forvxyz opened
   clp=&clar[vx][vy][vz]; clf=&edp->clarf[vx][vy][vz];
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
std::cout << std::endl;

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
#endif

#if(FITOPB==1)
// calc shad
float szgrid=(float)(ux-lx+1)*(uy-ly+1)*(uz-lz+1);
if((szgrid-bdes)==0) shad=-1;
else opened
   shad=(float)(szgrid-bous-bdes-SHCOE2*(bdes-xins))/(szgrid-bdes);
   if(shad<0) shad=0;
closed
#endif

return (double)shad;
/**/
closed
