
#define _CRT_SECURE_NO_DEPRECATE

// LABELS & abbreviations ---------------------------------------------------
#define NO         0
#define YA         1
#define opened     {
#define closed     }

// Frequently modified ------------------------------------------------------
#define RNSEED     101
#define MYOMP      0
#define MYMPI      0
#define NCORES     1
#define CYCLES     40000
#define SHOWPACE   500
#define AMSCOE     1.0 // (Dvfreeze)
#define NEXSPL     10 // number of examples for neural sample
#define RELOAD     NO
#define REPLAY     NO
#define SCAGL2     NO
#define NICHEX     NO // for niching in fitness fct
#define NICHECOE   (float)0.01

// Key parameters, always ---------------------------------------------------
#define NEVENT     10 // (other)
#define COLDGENS   0 // span of gens for temp. freezing of past genes
#define N2DRAT     5 // normal to driver ratio
#define DOP2NS     5 // doper2 neighbourhood size
#define NDQBLC     7,7,5 // normal to driver quot for 'marble' block
// Key parameters, paropt ---------------------------------------------------
#define PAROPT     0 // parameter optimisation (1: pattern recogn., 2: TSP)
#define RNDRNG     (COLRNG/16) // initial rnd range for colour values
#define SMOOTHSZ   3 // smooth neighbourhood size
#define PARSETSZ   112860
#define PARSUBS    0,112860,0,112860,0,112860,0,112860
#define NPLATES    4 // 
#define CBSCALE    1.00
#define CBOUND0    -1, 0.00,0.10,0.00,0.10,0.00,0.04, -1, 0.00,0.10,0.10,0.20,0.00,0.04,
#define CBOUND1    -1, 0.10,0.20,0.00,0.10,0.00,0.04, -1, 0.10,0.20,0.10,0.20,0.00,0.04,
#define MAPCPN     1,1,1,1 // map copy number
#define NCAT       10 // number of categories for neural tot
#define NEXTOT     1000 // number of examples for neural tot
#define LRSTOT     4 // neural network layers
#define NDSPLR     opened 784,900,1000,1100, closed // nodes per layer (cumulative numbers)
#define NDSTOT     1100 // neural network tot number of nodes
#define LASTLR     "M 0 0 0 00 A 50 100 150 200"
#define NCACSET    "M 0 0 0 0 A 10 10 10 10" // number of categories for neural currently evolved
#define NETMAP     1
#define CNXTYPE    0
#define HRNDSPL    0 // half random sample
#define GQTRGA     0,3,6
#define MGDCOA     (float)0.00, (float)0.00, (float)0.00 // typical value: 0.50 
#define RSHCOA     (float)0.00, (float)0.00, (float)0.00 // typical value: 0.50
#define RSKCOA     (float)0.00, (float)0.00, (float)0.00 // typical value: 0.50
#define CGDCOA     (float)0.00, (float)0.00, (float)0.00 // typical value: 0.50
#define XXBADA     "M 0 0 0 A 10 10 10 B 10 10 10 C 4 4 4"
// Key parameters, anubs ----------------------------------------------------
#define TAGPACE    10
#define PLSTEVAL   0.015 // mutation prob. for stepeval
#define PLGENETM   0.015 // mutation prob. for gene timers

// Big options --------------------------------------------------------------
#define ISGRID     0 // option on initial grid status
#define NICHEX     NO // for niching in fitness fct
#define REMRED     NO // for remove-redeploy procedure
#define TMON       YA // to disable timers
#define PGFREEZE   YA // progressive freezing
#define PGFTEVAL   NO // progressive fitness evaluation
#define EVFTEVAL   NO // evolvable fitness evaluation
#define MOCDEV     NO // if NO it's used only moc set of the best (PF)
#define SPHERE     0 // uses spheres
#define TAGCHK     NO // tagchk
#define DOPAON     NO // Dop always ON with PGFREEZE=YA
#define GRNMODE    NO // 
#define NOAPOPT    YA // 
#define NOCANC     0  // No cancellation of cells, for val=1 (mothers), =2 (all)
#define NOSCCH     0  // No stem cells change, for val=1 (doper2), =2 (all)

// Space related ------------------------------------------------------------
#define NDIMS      3  // if() option. Number of dimensions
#define COLRNG     16
#define COLTOL    -1  // below colour differences are not counted
#define GRIDX      38
#define GRIDY      134
#define GRIDZ      94
#define SZGRID     (GRIDX*GRIDY*GRIDZ)
#define ZYGSSS     opened 19,67,21,-1,-1,-1,-1,-1,-1,-1,-1,-1 closed
#define ZYGTOT     1
#define MSAMAX     15
#define DPLMAX     20

// Fitness related ----------------------------------------------------------
#define FITOPA     1 // fitness option A
#define FITOPB     0 // fitness option B
#define FITEXP     2.5
#define FITTOL     2 // fitness tolerance
#define SHCOE0     1.00 // shad coef ([0-1] -> [b&w-color])
#define SHCOE2     1.00 // shad coef

// Big array size -----------------------------------------------------------
#define DGARSZ     2020  // dgar total elements
#define DHARLS     20000
#define DVPTSMAX   100000

// Other --------------------------------------------------------------------
#define SMOCMA     5000
#define SMOCMB     19999
#define NAPMIN     1
#define FSCMAX     0
#define FITHST     500
#define GPSTEP     2
#define SHOWED     1
#define ZSTEPS     NO
#define NSTEPS     20
#define VTKOPT     1
#define NMORPH     4
#define HSTSTP     500 // (other)
#define CFLAG0     YA
#define DEBUG0     YA

// Dvfreeze -----------------------------------------------------------------
#define AMSDGN     20  // gene span
#define AMSDGO     0  // gene overlap
#define AMSDLZ     13  // size of parallelepiped
#define AMSLLT     0   // timer option (0 = one by one, 1 = all together)
#define ESTAGE     100  // evaluation stage
#define STAGES     101  // number of clock values
#define MOCLEN     20  // MOC length //#define CLKMAX    80  // number of clock values

// Files
#define  CDTGFN   "Dinputx/aaanubi_038134094_voxgr.txt"
#define   SMBFN   "Dinputx/MNIST_E/mnist"
#define LASTAFN   "Dinputx/lastall.txt"
#define OUTDIRN   "Doutput/"
#define  CONSFN   "Doutput/bxconsole.txt"
#define   GENFN   "Doutput/bxgen.txt"
#define   HSTFN   "Doutput/bxhst.txt"
#define  TAG2FN   "Doutput/bxtag2.txt"
#define  TAG3FN   "Doutput/bxtag3.txt"
#define  DGARFN   "Doutput/bxdgar.txt"
#define  DHLSFN   "Doutput/bxdhls.txt"
#define GRDFILE   "Doutput/exgrdfn"
#define XDBG0FN   "Doutput/exdbg0.txt"
#define XDBG1FN   "Doutput/exdbg1.txt"
#define XDBG2FN   "Doutput/exdbg2.txt"
#define XDBG3FN   "Doutput/exdbg3.txt"
#define XDBG4FN   "Doutput/exdbg4.txt"
#define XDBG5FN   "Doutput/exdbg5.txt"
#define XDBG6FN   "Doutput/exdbg6.txt"

// Variables' size in the Genome
// Genes, left  part
#define XORDSQ    5 //order of precedence field
#define XTMRSQ    4 //timer               field
#define RTMRSQ    1 //timer relevance     field
#define XRESSQ    7 //regulatory set      field
// Genes, right part --------------------------------------------------------
#define XMSASQ    2 //master switch       field
#define XMSBSQ    1 //master switch       field
#define XDGDSQ    3 //diagonal distance   field
#define XDPLSQ    6 //displacement        field
#define XCOLSQ    2 //color number        field
#define APOSSQ    6 //parmap position     field
// Totals -------------------------------------------------------------------
#define DGLPSQ    (1+XORDSQ+XTMRSQ+RTMRSQ+(XRESSQ*MOCLEN)+MOCLEN)
#define DGRASQ    ((6*XDPLSQ)+XMSASQ+XMSBSQ+(9*XDGDSQ)+(9*XCOLSQ))
#define DGRBSQ    (3*APOSSQ)
#define DGENSQ    (DGLPSQ+DGRASQ+DGRBSQ)
#define DGARSQ    (DGENSQ*DGARSZ) // DGAR total bases
// grn structure ------------------------------------------------------------
#define SGIASZ    2 // if array size
#define SGSCSQ    1 // switch
#define SGIFSQ    7 // if part
#define SGTHSQ    7 // then part
#define SGSNSQ    1 // contribution sign
#define SGCDSQ    3 // contribution code
#define SGENSQ    (SGSCSQ+(SGIFSQ*SGIASZ)+SGTHSQ+SGSNSQ+SGCDSQ) // entire gene

// GA -----------------------------------------------------------------------
#define NPOPUL    1 //12
#define GHDRSZ    (3+XTMRSQ) // gen header size
#define BEGNSQ    (GHDRSZ)
#define XGENSZ    (GHDRSZ+DGARSQ)
#define PCROSS    0.5
#define PMUTAT    0.5
#define MPCMUT    0.005
#define ATPMAX    60
#define POPSZ0    256
#define POPSZ1    0
#define POPSZT    (POPSZ0+POPSZ1)
#define NKIDS0    224 // must be: POPSZ0 > NKIDS0+NKIDS1
#define NKIDS1    0
#define NICHESS   1
#define NICHEA    2
#define NNGUYS    (NPOPUL*POPSZT)

// Cell conditions
#define CLCHOLE  -1 // cellular hole as a result of cell died
#define CLCHOLN  -2 // cellular hole just created
#define CLALIVE   1 // cell alive
#define CLALIVN   2 // cell just created
#define CLDESTR   6 // cell destroyed
// xxxx
#define Ldhn      0
#define Ldrv      1
#define Lcnd      2
#define Lcol      3
#define Lxxx      4
#define Forvxyz   for(vx=0;vx<GRIDX;vx++) for(vy=0;vy<GRIDY;vy++) for(vz=0;vz<GRIDZ;vz++)
#define sizeof_gen_ (sizeof(bas)*XGENSZ)
