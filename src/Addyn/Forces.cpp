#include "Forces.h"
#include <stdarg.h>
void Forces::Shout(const char * s, ...){
#ifndef DIN_DEBUG
  va_list args;
  va_start(args, s);
  fprintf(stderr, ">");
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
#else  
  return;
#endif
}
Forces::Forces(int argc,char **argv,int NInEdge,char *ConfFileExt){
  Shout("Constructor/no starting configuration");
  InitConst();
  sprintf(ConfFile,"%s",ConfFileExt);
  char ConfF[40];
  sprintf(ConfF,ConfFile);
  if(ReadConfDinamica(ConfF)){
    Dx = 1./(double)(NEdge-1);
  }
  if(VAR_IF_TYPE(SysShape,SYS_2D)){
    nEdge[0] = NEdge;
    double Ratio = pEdge(1)*pInvEdge(0);
    nEdge[1] = (int)(nEdge[0]*Ratio + 0.0001);
    for(int i=0;;i++){
      if((int)(nEdge[1]/Ratio) == nEdge[0]) break;
      nEdge[0]++;
      nEdge[1] = (int)(nEdge[0]*Ratio + 0.0001);
      printf("using: nEdge[0] %d nEdge[1] %d\n",nEdge[0],nEdge[1]);
      if(i>=10){
	printf("Could not find the appropriate border partition %d %d\n",nEdge[0],nEdge[1]);
	exit(0);
      }
    }
    SetNLink(4);
    ReSetNPart(nEdge[0]*nEdge[1]);
    ReSetNChain(nEdge[1]);
    SetNPCh(nEdge[0]);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_3D)){
    SetNLink(6);
    ReSetNChain(NEdge*NEdge);
    ReSetNPart(NEdge*NEdge*NEdge);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    SetNLink(1);
    ReSetNChain(1);
    ReSetNPart(NEdge);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_STALK)){
    SetNLink(3);
    ReSetNPart(4*NEdge);
    ReSetNChain(4);
    //    Kf.El[0] *= Gen->NPart;
  }
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    SetNLink(3);
    ReSetNPart(NEdge);
    ReSetNChain(1);
    //    Kf.El[0] *= Gen->NPart;
  }
  else if(VAR_IF_TYPE(SysShape,SYS_PORE)){
    SetNLink(3);
    ReSetNPart(NEdge);
    ReSetNChain(1);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_1D)){
    SetNLink(2);
    ReSetNPart(NEdge);
    ReSetNChain(1);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_RIGID)){
    ReSetNPart(0);
    SetNLink(0);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_TRIAL)){
    ReSetNPart(3*3*3);
    SetNLink(0);
    ReSetNChain(1);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MD)){
    SetNLink(0);
    ReSetNPart(NEdge);
    ReSetNChain(1);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MC)){
    SetNLink(0);
    ReSetNPart(NEdge);
    ReSetNChain(1);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    SetNLink(4);
    ReSetNPart(NEdge+NSpline);
    ReSetNChain(1);
  }
  ReSetNPCh(pNPart()/pNChain());
  //SetNPCh(NEdge);
  SetDeltat(Deltat);
  SetStep(0);
  SetNType(2);
  CreateInitial();
  {
    MInt = new MatInt(3,3);
    MInt->SetCoeff(-24.33,0,0);
    MInt->SetCoeff(-7.22,0,1);
    MInt->SetCoeff(-24.33,0,2);
    MInt->SetCoeff(-0.1,1,1);
    MInt->SetCoeff(-7.22,1,2);
    MInt->SetCoeff(0.,2,2);
    MInt->SetCoeff(3.,0,0,0);
    MInt->SetCoeff(3.,0,0,1);
    MInt->SetCoeff(3.,0,0,2);
    MInt->SetCoeff(3.,0,1,1);
    MInt->SetCoeff(3.,0,1,2);
    MInt->SetCoeff(3.,0,2,2);
    MInt->SetCoeff(0.,1,1,1);
    MInt->SetCoeff(3.,1,1,2);
    MInt->SetCoeff(3.,1,2,2);
    MInt->SetCoeff(0.,2,2,2);
  }
  AllocMethod();
  PrepareSys();
  PrepareParallel(argc,argv);
  if(VAR_IF_TYPE(SysShape,SYS_LEAVES)) IfNano = 2;
  if(VAR_IF_TYPE(SysShape,SYS_PORE)) IfNano = 2;
  if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    Pm[NEdge].Typ = 2;
    Pm[NEdge+1].Typ = 2;
  }
  //Pc->PrintCells();
  //for(int p=0;p<pNPart();p++) CheckDomDec(p);
  //Interp();
}
Forces::~Forces(){
  Shout("Freeing");
  fclose(StatFile1);
  fclose(StatFile2);
  Shout("Freeing/md: forces domain decomposition");  
  if(VAR_IF_TYPE(SysAlloc,ALL_FORCES)){
    delete [] Fm;
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_MD)){
    delete Pc;
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_MC)){
    delete Pc;
    delete [] OldNrgBead;
    delete [] OldNrgCh;
    for(int p=0;p<pNPCh();p++){
      delete [] OldPos[p];
    }
    delete [] OldPos;
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_BIAS)){
    for(int t=0;t<NTrialBias;t++){
      delete [] BondPosBias[t];
    }
    delete [] BondPosBias;
    delete [] CumProbBias;
  }  
  if(VAR_IF_TYPE(SysAlloc,ALL_SPLINE)){
    free(Pl);
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_DENS)){
    free(LocDens2);
    free(LocDens3);
    free(Dens2);
    free(Dens3);
  }
  if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    delete IntMatrix;
  }
  if(VAR_IF_TYPE(SysShape,SYS_2D)){
    delete IntMatrix;
  }
  if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    delete IntMatrix;
  }
#ifdef __glut_h__
  free(Cylinder);
#endif
}
Forces::Forces(int argc,char **argv,char *ConfF,char *Snapshot){
  Shout("constructor with intial configuration");
  InitConst();
  if(ReadConfDinamica(ConfF)){
    Dx = 1./(double)(NEdge-1);
  }
  Open(Snapshot,BF_PART);
  if(VAR_IF_TYPE(SysShape,SYS_LEAVES)) IfNano = 2;
  if(VAR_IF_TYPE(SysShape,SYS_PORE)) IfNano = 2;
  AllocMethod();
  PrepareParallel(argc,argv);
  PrepareSys();
 //Interp();
}
void Forces::PrepareSys(){
  Shout("Preparing system");
  OldNrgSys = 0.;
  DefNanoForceParam();
  /*calculate the energies per chain or per particle for
    the calculation of the MC.
  */
  if(VAR_IF_TYPE(SysShape,SYS_MC)){
    if(VAR_IF_TYPE(SysAlloc,ALL_DENS)){
      Shout("Preparing system/adding densities\n");
      OldNrgSys = DensFuncNrgSys();
      ChemPotEx += NrgPBead;
      if(VAR_IF_TYPE(CalcMode,CALC_NcVT)){
	//CalcTotNrgCh();
      }
      if(VAR_IF_TYPE(CalcMode,CALC_NVT)){
	;//CalcTotNrgBead();
      }
      if(VAR_IF_TYPE(CalcMode,CALC_mcVT)){
	if(!VAR_IF_TYPE(CalcMode,CALC_CONF_BIAS))
	  ;//CalcTotNrgCh();
      }
    }
    //defining all the allocated chains
    for(int p=0;p<pNAllocP();p++){
      int c = (int)(p/(double)pNPCh());
      Pm[p].CId = c;
      Pm[p].Typ = 0;
      if( p%pNPCh() >= Block[0].Asym ) Pm[p].Typ = 1;
      if( p%pNPCh() == pNPCh() - 1) continue;
      Ln[p].NLink = 1;
      Ln[p].Link[0] = p+1;
    }
  }
  if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    for(int p=0;p<NEdge;p++){
      Pm[p].Typ = 2;
      Pm[p].Idx = p;
      Pm[p].CId = p;
      Ln[p].NLink = 0;
      SetBkf(p);
    }
    for(int p=NEdge;p<NEdge+NSpline;p++){
      Pm[p].Typ = 0;
      Pm[p].Idx = p;
      Pm[p].CId = NEdge;
      Ln[p].NLink = 1;
      Ln[p].Link[0] = p+1;
      if(p==NEdge+NSpline-1) Ln[p].NLink = 0;
      SetBkf(p);
    }
  }
  if(VAR_IF_TYPE(SysShape,SYS_MD)){
    Shout("Preparing system/calculating forces\n");
    OldNrgSys = SumForcesMD();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    IntMatrix = new Matrice(pNPart(),pNPart());
  }
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    IntMatrix = new Matrice(pNPCh(),pNPCh());
  }
  else if(VAR_IF_TYPE(SysShape,SYS_2D)){
    IntMatrix = new Matrice(pNPart(),pNPart());
  }
  /*for the chemical potential, to avoid the calculation of
    large exponential */
  NrgPBead = 0.;//2.*OldNrgSys/(double)pNPart();
  Shout("Prepared");
}
void Forces::FillMatrix(){
  if(!IfFillMatrix) return;
  IntMatrix->Clear();
  double Inter = pEdge(CLat1)/(double)NEdge;//fabs(pPos(1,0) - pPos(0,0));
  SPLINE Weight;
  Weight.a0 = Kf.El[2];
  Weight.a1 = 0./Inter;
  Weight.a2 = Kf.Lap/SQR(Inter);
  Weight.a3 = 0./(Inter*SQR(Inter));
  Weight.a4 = Kf.SLap/(SQR(Inter)*SQR(Inter));
  if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    int NDim = 1;
    Matrice *CoeffMatrix = new Matrice(Weight,NDim);
    CoeffMatrix->Print();
    for(int r=0;r<pNPart();r++){
      if(Pm[r].Typ != 0){ IntMatrix->Set(r,r,1.);continue;}
      if(r >= 2) IntMatrix->Set(r,r-2,CoeffMatrix->Val(2,0));
      if(r >= 1) IntMatrix->Set(r,r-1,CoeffMatrix->Val(2,1));
      if(r < pNPart()-1) IntMatrix->Set(r,r+1,CoeffMatrix->Val(2,3));
      if(r < pNPart()-2) IntMatrix->Set(r,r+2,CoeffMatrix->Val(2,4));
      IntMatrix->Set(r,r,CoeffMatrix->Val(2,2));
    }
    IntMatrix->Invert();
    //IntMatrix->Print();
    delete CoeffMatrix;
  }
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    int NDim = 1;
    Matrice *CoeffMatrix = new Matrice(Weight,NDim);
    CoeffMatrix->Print();
    for(int r=0;r<pNPCh();r++){
      if(Pm[r].Typ != 0){ IntMatrix->Set(r,r,1.);continue;}
      if(r >= 2) IntMatrix->Set(r,r-2,CoeffMatrix->Val(2,0));
      if(r >= 1) IntMatrix->Set(r,r-1,CoeffMatrix->Val(2,1));
      if(r < pNPart()-1) IntMatrix->Set(r,r+1,CoeffMatrix->Val(2,3));
      if(r < pNPart()-2) IntMatrix->Set(r,r+2,CoeffMatrix->Val(2,4));
      IntMatrix->Set(r,r,CoeffMatrix->Val(2,2));
    }
    //IntMatrix->Invert();
    IntMatrix->Print();
    delete CoeffMatrix;
  }
  else if(VAR_IF_TYPE(SysShape,SYS_2D)){
    int NDim = 2;
    Matrice *CoeffMatrix = new Matrice(Weight,NDim);
    CoeffMatrix->Print();
    for(int p=0;p<pNPart();p++){
      if(Pm[p].Typ != 0){
	IntMatrix->Set(p,p,1.);
	continue;
      }
      int pym1 = Ln[p].Link[0];
      int pyp1 = Ln[p].Link[1];
      int pym2 = Ln[pym1].Link[0];
      int pyp2 = Ln[pyp1].Link[1];
      int pxm1 = Ln[p].Link[2];
      int pxp1 = Ln[p].Link[3];
      int pxm2 = Ln[pxm1].Link[2];
      int pxp2 = Ln[pxp1].Link[3];
      // printf("%d)\n",p);
      //printf("%d %d %d %d\n",pym2,pym1,pyp1,pyp2);
      // printf("%d %d %d %d\n",pxm2,pxm1,pxp1,pxp2);
      if(PeriodicImage[0]){
	IntMatrix->Set(p,pxm2,CoeffMatrix->Val(2,0));
	IntMatrix->Set(p,pxm1,CoeffMatrix->Val(2,1));
	IntMatrix->Set(p,pxp1,CoeffMatrix->Val(2,3));
	IntMatrix->Set(p,pxp2,CoeffMatrix->Val(2,4));
      }
      else{
	if(pxm2 == p-2*nEdge[1])
	  IntMatrix->Set(p,pxm2,CoeffMatrix->Val(2,0));
	if(pxm1 == p-nEdge[1])
	  IntMatrix->Set(p,pxm1,CoeffMatrix->Val(2,1));
	if(pxp1 == p+nEdge[1])
	  IntMatrix->Set(p,pxp1,CoeffMatrix->Val(2,3));
	if(pxp2 == p+2*nEdge[1])
	  IntMatrix->Set(p,pxp2,CoeffMatrix->Val(2,4));
      }
      IntMatrix->Add(p,p,CoeffMatrix->Val(2,2));
      if(PeriodicImage[1]){
	IntMatrix->Set(p,pym2,CoeffMatrix->Val(2,0));
	IntMatrix->Set(p,pym1,CoeffMatrix->Val(2,1));
	IntMatrix->Set(p,pyp1,CoeffMatrix->Val(2,3));
	IntMatrix->Set(p,pyp2,CoeffMatrix->Val(2,4));
      }
      else{
	if(pym2 == p-2)
	  IntMatrix->Set(p,pym2,CoeffMatrix->Val(2,0));
	if(pym1 == p-1)
	  IntMatrix->Set(p,pym1,CoeffMatrix->Val(2,1));
	if(pyp1 == p+1)
	  IntMatrix->Set(p,pyp1,CoeffMatrix->Val(2,3));
	if(pyp2 == p+2)
	  IntMatrix->Set(p,pyp2,CoeffMatrix->Val(2,4));
      }
    }
    IntMatrix->Invert();
    //IntMatrix->Print();
    delete CoeffMatrix;
  }
  IfFillMatrix = 0;
}
void Forces::PrepareParallel(int argc,char **argv){
#ifdef OMPI_MPI_H 
  Shout("Preparing parallelisation");
  MPI_Init(&argc,&argv);
  int Rank=0,Size=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
  MPI_Comm_size(MPI_COMM_WORLD, &Size);
  int Partition = (int)(argc/(double)Size);
  NFile[0] = Partition*Rank;
  NFile[1] = Partition*(Rank+1);
  if(Rank==Size-1) NFile[1] += argc%Size;
  Proc = new SingProc(Size,Rank);
#endif
}
void Forces::AllocMethod(){
  Shout("Allocation");
  StatFile1 = fopen("StatDyn1.dat","w");
  StatFile2 = fopen("StatDyn2.dat","w");
  ChooseCalcMode(CalcMode);
  ChoosePot(CalcMode);
  ChooseThermostat(ThermMode);
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  double Pos[3];
  //Pc->PrintCells();
  //CheckPairList();
  SetDeltat(Deltat);
  if(VAR_IF_TYPE(SysShape,SYS_MD)){
    Shout("Allocating/md: forces, domain decomposition");
    Fm = (FORCES *)calloc(pNAllocP(),sizeof(FORCES));
    VAR_ADD_TYPE(SysAlloc,ALL_MD);
    VAR_ADD_TYPE(SysAlloc,ALL_FORCES);
    Pc = new DomDec(Edge,pNPart(),sqrt(Kf.CutOff2));
    for(int p=0;p<pNPart();p++){pPos(p,Pos);Pc->AddPart(p,Pos);}
  }
  if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    Shout("Allocating/rod: forces");
    Fm = (FORCES *)calloc(pNAllocP(),sizeof(FORCES));
    VAR_ADD_TYPE(SysAlloc,ALL_FORCES);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MC)){
    Shout("Allocating/mc: domain decomposition old chain positions, old energies for particles and chains, first bead distribution, bias (cumulative probabilities and bead positions)");  
    Pc = new DomDec(Edge,pNPart(),sqrt(Kf.CutOff2));
    double Pos[3];
    for(int p=0;p<pNPart();p++){pPos(p,Pos);Pc->AddPart(p,Pos);}
    //Pc->PrintCells();
    if(VAR_IF_TYPE(CalcMode,CALC_NVT)){
      OldNrgBead = new double[3*pNAllocP()];
      if(!OldNrgBead){
	printf("Could not allocate OldNrgPm \n");
	exit(1);
      }
    }
    OldNrgCh = new double[3*pNAllocC()];
    FirstBeadDistr = new double[NBin];
    if(!OldNrgCh){
      printf("Could not allocate OldNrgCh\n");
      exit(1);
    }
    OldPos = (double **)calloc(pNPCh(),sizeof(double));
    if( !OldPos){
      printf("Could not allocate OldPos\n");
      exit(1);
    }
    for(int p=0;p<pNPCh();p++){
      OldPos[p] = (double *)calloc(3,sizeof(double));
    }
    CumProbBias = new double[NTrialBias];
    BondPosBias = new double*[NTrialBias];
    for(int t=0;t<NTrialBias;t++){
      BondPosBias[t] = new double[3];
    }
    VAR_ADD_TYPE(SysAlloc,ALL_MC);
    VAR_ADD_TYPE(SysAlloc,ALL_BIAS);
    GaussVar = sqrt(pReOverCutOff()/(3.*(pNPCh()-1)))/2.;
    GaussVar = sqrt(1./pkSpr());
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    Shout("Allocating/mc: domain decomposition old chain positions, old energies for particles and chains, first bead distribution, bias (cumulative probabilities and bead positions)");  
    Pc = new DomDec(Edge,pNPart(),sqrt(Kf.CutOff2));
    double Pos[3];
    for(int p=0;p<pNPart();p++){pPos(p,Pos);Pc->AddPart(p,Pos);}
    //Pc->PrintCells();
    if(VAR_IF_TYPE(CalcMode,CALC_NVT)){
      OldNrgBead = new double[3*pNAllocP()];
      if(!OldNrgBead){
	printf("Could not allocate OldNrgPm \n");
	exit(1);
      }
    }
    OldNrgCh = new double[3*pNAllocC()];
    FirstBeadDistr = new double[NBin];
    if(!OldNrgCh){
      printf("Could not allocate OldNrgCh\n");
      exit(1);
    }
    OldPos = (double **)calloc(pNPCh(),sizeof(double));
    if( !OldPos){
      printf("Could not allocate OldPos\n");
      exit(1);
    }
    for(int p=0;p<pNPCh();p++){
      OldPos[p] = (double *)calloc(3,sizeof(double));
    }
    CumProbBias = new double[NTrialBias];
    BondPosBias = new double*[NTrialBias];
    for(int t=0;t<NTrialBias;t++){
      BondPosBias[t] = new double[3];
    }
    VAR_ADD_TYPE(SysAlloc,ALL_MC);
    VAR_ADD_TYPE(SysAlloc,ALL_BIAS);
  }
  if(VAR_IF_TYPE(SysShape,SYS_1D) || VAR_IF_TYPE(SysShape,SYS_TRIAL)){
    Shout("Allocating/splines");  
    Pl = (PART *)calloc(NSpline,sizeof(PART));
    VAR_ADD_TYPE(SysAlloc,ALL_SPLINE);
  }
  if(VAR_IF_TYPE(CalcMode,CALC_DENS)){
    Shout("Allocating/particle densities, sum of local densities");  
    LocDens2 = (double *)calloc(pNPCh()*pNType(),sizeof(double));
    LocDens3 = (double *)calloc(pNPCh()*pNType(),sizeof(double));
    Dens2 = (double *)calloc(pNAllocP()*pNType(),sizeof(double));
    Dens3 = (double *)calloc(pNAllocP()*pNType(),sizeof(double));
    if( !Dens2 || !Dens3){
      printf("Could not allocate Dens2 Dens3\n");
      exit(1);
    }
    VAR_ADD_TYPE(SysAlloc,ALL_DENS);
  }
#ifdef __glut_h__
  Cylinder = (GLuint *)calloc(pNNano(),sizeof(GLuint));
#endif  
  ChemPotId = log(NChemPotId/pVol());
}
void Forces::InitConst(){
  Shout("Init constants");  
  NEdge = 80+1;//NInEdge;
  Dx = 1./(double)(NEdge-1);
  Kf.El[0] = 11.;//2.*NEdge/10.;//Elastic coupling
  Kf.El[1] = 11.;//Elastic coupling
  Kf.El[2] = 11.;//Elastic coupling
  Kf.Lap = 20.*pow(Dx,2.);//bending rigidity
  Kf.SLap = 0.;//.015*pow(Dx,2.);//Surface tension
  Kf.Ext = 1.;
  Kf.LJ = 1.;
  Kf.LJMin = .5;
  Kf.CutOff2 = 1.;
  Kf.Cont = 0.;//.001;
  Kf.Elong[0] = 0.1;
  Kf.Elong[1] = 0.1;
  Kf.Elong[2] = 0.1;
  Kf.ForThr = 100.;
  IntMax = 100;
  IncrDist = 0.01;
  IfInterp=FIT_FORTH;
  IfFillMatrix=1;
  Nano->Rad = .05;
  Nano->Height = .3;
  ChemPotId = 0.;
  ChemPotEx = 0.;
  Nano->Pos[0] = 0.;
  Nano->Pos[1] = 0.;
  Nano->Pos[2] = .5;
  BoundCond[0] = 1;BoundCond[1] = 1;
  BoundCond[2] = 1;BoundCond[3] = 1;
  BoundCond[4] = 1;BoundCond[5] = 1;
  PeriodicImage[0] = 1;PeriodicImage[1] = 1;PeriodicImage[2] = 1;
  Time = 0.;
  NUpdate = 100;
  NWrite = 1000;
  NBin = 100;
  NTrialBias = 10;
  //  GaussVar = sqrt(1./pkSpr());
  ThermMode = THERM_LANG;
  CalcMode = CALC_LJ39;
  //  Part2Move = NEdge/2+NEdge;
#ifdef __glut_h__
  NSpline = 100;
  IfSphere = 0;
  IfLine = 1;
  IfSpline = 0;
  NShow = 1;
  IfMovie=0;
  Frame = 0;
  IfExt = 0;
  IfRot = 0;
  BeadType = 0;
#endif
  IfExit = 0;
  IfNano = 0;
  DynFlag = 0;
  SetNBlock(1);
}
void Forces::Info(){
  printf("-----INFO-------\n");
  if(VAR_IF_TYPE(SysShape,SYS_MC)){
    // Pc->Erase();
    // ClearDens();
    // for(int p=0;p<pNPart();p++)
    //   Pc->AddPart(p,Pm[p].Pos);
    // AddDens(0,pNPart());
    printf("Energy: Sys %lf Ch  Part \n",OldNrgSys);
    // CalcTotNrgCh();
    // for(int c=0;c<pNChain();c++)
    //   printf("%d %lf %lf %lf\n",c,OldNrgCh[c*3],OldNrgCh[c*3+1],OldNrgCh[c*3+2]);
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MD)){
    for(int p=0;p<pNPart();p++){
      if(p == Bead2Move) printf("p      ");
      printf("%d) (%lf,%lf,%lf) (%lf %lf %lf) %lf-%lf-%lf\n",p,pPos(p,0),pPos(p,1),pPos(p,2),pVel(p,0),pVel(p,1),pVel(p,2),Fm[p].Dir[2],Fm[p].Dir[2],Fm[p].Dir[2]);
    }
  }
  printf("------------\n");
}
int Forces::ReadConfDinamica(char *InFile){
  SysType = 0;
  SysShape = 0;
  CalcMode = 0;
  FILE *FileToRead = fopen(InFile,"r");
  if(FileToRead == NULL){
    printf("The conf file %s is missing\n",InFile);
    return 1;
  }
  double buff[12];
  char SysInit[20];
  char *Line = (char *)malloc(256*sizeof(char));
  //fgets(Line,256,FileToRead);
  int NNano = 0;
  for(int k=0;!(fgets(Line,256,FileToRead)==NULL);k++){
    if(strstr(Line, "Rigid") == Line) NNano++;
  }
  SetNNano(NNano);
  NNano = 0;
  rewind(FileToRead);
  for(int k=0;!(fgets(Line,256,FileToRead)==NULL);k++){
    //printf("%s",Line);
    if(1 == sscanf(Line,"NEdge %lf",buff) )
      NEdge = (int)*buff;
    else if(1 == sscanf(Line,"SysShape %s",SysInit) ){
      if(!strcmp(SysInit,"leaves") )
	VAR_ADD_TYPE(SysShape,SYS_LEAVES);
      else if(!strcmp(SysInit,"pore") )
	VAR_ADD_TYPE(SysShape,SYS_PORE);
      else if(!strcmp(SysInit,"1d") )
	VAR_ADD_TYPE(SysShape,SYS_1D);
      else if(!strcmp(SysInit,"2d") )
	VAR_ADD_TYPE(SysShape,SYS_2D);
      else if(!strcmp(SysInit,"3d") )
	VAR_ADD_TYPE(SysShape,SYS_3D);
      else if(!strcmp(SysInit,"rod") )
	VAR_ADD_TYPE(SysShape,SYS_ROD);
      else if(!strcmp(SysInit,"trial") )
	VAR_ADD_TYPE(SysShape,SYS_TRIAL);
      else if(!strcmp(SysInit,"rigid") )
	VAR_ADD_TYPE(SysShape,SYS_RIGID);
      else if(!strcmp(SysInit,"stalk") )
	VAR_ADD_TYPE(SysShape,SYS_STALK);
      else if(!strcmp(SysInit,"md") )
	VAR_ADD_TYPE(SysShape,SYS_MD);
      else if(!strcmp(SysInit,"mc") )
	VAR_ADD_TYPE(SysShape,SYS_MC); 
      else if(!strcmp(SysInit,"electro") ){
	VAR_ADD_TYPE(SysShape,SYS_ELECTRO);
	//VAR_ADD_TYPE(SysShape,SYS_MC); 
      }
      else{
	printf("system type not recognized\n");
	exit(1);
      }
    }
    else if(1 == sscanf(Line,"CalcMode %s",SysInit) ){
      if(!strcmp(SysInit,"NVT") )
	VAR_ADD_TYPE(CalcMode,CALC_NVT);
      else if(!strcmp(SysInit,"NcVT") )
	VAR_ADD_TYPE(CalcMode,CALC_NcVT);
      else if(!strcmp(SysInit,"mcVT") )
	VAR_ADD_TYPE(CalcMode,CALC_mcVT);
      else if(!strcmp(SysInit,"mVT") )
	VAR_ADD_TYPE(CalcMode,CALC_mVT);
      else{
	printf("calculation type not recognized\n");
	exit(1);
      }
    }
    else if(1 == sscanf(Line,"Thermostat %s",SysInit) ){
      if(!strcmp(SysInit,"Langevin") )
	ThermMode = THERM_LANG;
      else if(!strcmp(SysInit,"Andersen") )
	ThermMode = THERM_AND;
      else if(!strcmp(SysInit,"Berendsen") )
	ThermMode = THERM_BERE;
      else if(!strcmp(SysInit,"no") )
	ThermMode = THERM_NO;
      else{
	printf("thermostat not recognized\n");
	exit(1);
      }
    }
    else if(1 == sscanf(Line,"PotentialMode %s",SysInit) ){
      if(!strcmp(SysInit,"Pair") )
	VAR_ADD_TYPE(CalcMode,CALC_PAIR);
      else if(!strcmp(SysInit,"DensFunc") )
	VAR_ADD_TYPE(CalcMode,CALC_DENS);
      else if(!strcmp(SysInit,"DensFuncCh") )
	VAR_ADD_TYPE(CalcMode,CALC_DENS);
    }
    else if(1 == sscanf(Line,"Potential %s",SysInit) ){
      if(!strcmp(SysInit,"LJ") )
	VAR_ADD_TYPE(CalcMode,CALC_LJ);
      else if(!strcmp(SysInit,"LJ39") )
	VAR_ADD_TYPE(CalcMode,CALC_LJ39);
      else if(!strcmp(SysInit,"Harmonic") )
	VAR_ADD_TYPE(CalcMode,CALC_HARM);
      else if(!strcmp(SysInit,"Step") )
	VAR_ADD_TYPE(CalcMode,CALC_STEP);
      else if(!strcmp(SysInit,"Electro") )
	VAR_ADD_TYPE(CalcMode,CALC_ELECTRO);
      else{
	printf("interaction potential not recognized\n");
	exit(1);
      }
    }
    else if(1 == sscanf(Line,"IfInterp %lf",buff) )
      IfInterp = (int)*buff;
    else if(1 == sscanf(Line,"Lap %lf",buff) )
      Kf.Lap = *buff;
    else if(1 == sscanf(Line,"SLap %lf",buff) )
      Kf.SLap = *buff;
    else if(1 == sscanf(Line,"Ext %lf",buff) )
      Kf.Ext = *buff;
    else if(1 == sscanf(Line,"LJ %lf",buff) )
      Kf.LJ = *buff;
    else if(1 == sscanf(Line,"LJMin %lf",buff) )
      Kf.LJMin = *buff;
    else if(1 == sscanf(Line,"kSpr %lf",buff) )
      SetkSpr(*buff);
    else if(1 == sscanf(Line,"SprRest %lf",buff) )
      SetSprRest(*buff);
    else if(1 == sscanf(Line,"kBen %lf",buff) )
      SetkBen(*buff);
    else if(1 == sscanf(Line,"SimLimit %lf",buff) )
      SimLimit = (int)*buff;
    else if(1 == sscanf(Line,"CutOff %lf",buff) )
      Kf.CutOff2 = SQR(*buff);
    else if(1 == sscanf(Line,"Cont %lf",buff) )
      Kf.Cont = *buff;
    else if(1 == sscanf(Line,"NWrite %lf",buff) )
      NWrite = (int)*buff;
    else if(1 == sscanf(Line,"NBin %lf",buff) )
      NBin = (int)*buff;
    else if(1 == sscanf(Line,"NGrid %lf",buff) )
      NGrid = (int)*buff;
    else if(1 == sscanf(Line,"NUpdate %lf",buff) )
      NUpdate = (int)*buff;
    else if(strstr(Line, "Rigid") == Line){
      NanoString(Line,NNano++);
    }
#ifdef __glut_h__
    else if(1 == sscanf(Line,"IfMovie %lf",buff) )
      IfMovie = (int)*buff;
    else if(1 == sscanf(Line,"IfLine %lf",buff) )
      IfLine = (int)*buff;
    else if(1 == sscanf(Line,"NSpline %lf",buff) )
      NSpline = *buff;
    else if(1 == sscanf(Line,"IfSphere %lf",buff) )
      IfSphere = (int)*buff;
#endif
    else if(3 == sscanf(Line,"Edge %lf %lf %lf",buff,buff+1,buff+2) ){
      SetEdge(buff[0],0);
      SetEdge(buff[1],1);
      SetEdge(buff[2],2);
    }
    else if(3 == sscanf(Line,"El %lf %lf %lf",buff,buff+1,buff+2) ){
      Kf.El[0] = *buff;
      Kf.El[1] = *(buff+1);
      Kf.El[2] = *(buff+2);
    }
    else if(3 == sscanf(Line,"Elong %lf %lf %lf",buff,buff+1,buff+2) ){
      Kf.Elong[0] = *buff;
      Kf.Elong[1] = *(buff+1);
      Kf.Elong[2] = *(buff+2);
    }
    else if(6 == sscanf(Line,"Boundary %lf %lf %lf %lf %lf %lf",buff,buff+1,buff+2,buff+3,buff+4,buff+5)){
      BoundCond[0] = (int)buff[0];
      BoundCond[1] = (int)buff[1];
      BoundCond[2] = (int)buff[2];
      BoundCond[3] = (int)buff[3];
      BoundCond[4] = (int)buff[4];
      BoundCond[5] = (int)buff[5];
    }
    else if(3 == sscanf(Line,"Periodic  %lf %lf %lf",buff,buff+1,buff+2)){
      PeriodicImage[0] = (int)*buff;
      PeriodicImage[1] = (int)*(buff+1);
      PeriodicImage[2] = (int)*(buff+2);
    }
    else if(1 == sscanf(Line,"Deltat %lf",buff) )
      Deltat = *buff;
    else if(1 == sscanf(Line,"Temp %lf",buff) )
      SetTemp(buff[0]);
    else if(1 == sscanf(Line,"NChemPotId %lf",buff) )
      NChemPotId = *buff;
    else if(1 == sscanf(Line,"NTrialBias %lf",buff) ){
      NTrialBias =(int) *buff;
    }
    else if(1 == sscanf(Line,"ChemPotEx %lf",buff) )
      ChemPotEx = *buff;
    else if(1 == sscanf(Line,"IfConfBias %lf",buff) ){
      if((int)*buff == 1)
    	VAR_ADD_TYPE(CalcMode,CALC_CONF_BIAS);
    }
    else if(1 == sscanf(Line,"IfBilBias %lf",buff) ){
      if((int)*buff == 1){
    	VAR_ADD_TYPE(CalcMode,CALC_BIL_BIAS);
	StudySys();
      }
    }
    else if(1 == sscanf(Line,"IfSphBias %lf",buff) ){
      if((int)*buff == 1)
    	VAR_ADD_TYPE(CalcMode,CALC_SPH_BIAS);
    }
    else if(1 == sscanf(Line,"Viscosity %lf",buff) )
      Viscosity = *buff;
    else if(1 == sscanf(Line,"TNSlab %lf",buff) )
      Tens.NSlab = (int)*buff;
    else if(1 == sscanf(Line,"TNComp %lf",buff) )
      Tens.NComp = (int)*buff;
    else if(1 == sscanf(Line,"TCalcMode %s",SysInit) ){
      if(!strcmp(SysInit,"2d") ){
	VAR_ADD_TYPE(Tens.CalcMode,CALC_2d);
	Tens.NDim = 2;
      }
      else if(!strcmp(SysInit,"3d") ){
	VAR_ADD_TYPE(Tens.CalcMode,CALC_3d);
	Tens.NDim = 3;
      }
      else{
	printf("Pressure summation not recognized\n");
	exit(1);
      }
    }
  }
  //for(int n=0;n<pNNano();n++) for(int d=0;d<3;d++) Nano[n].Pos[d] *= pEdge(d);
  ChemPotId = log(NChemPotId/pVol());
  printf("Sys] NEdge %d SysShape %d %s Interp %d  \n",NEdge,SysShape,SysInit,IfInterp);
  printf("Forces] Lap %lf SLap %lf Ext %lf LJ %lf Cont %lf El %lf %lf %lf Elong %lf %lf %lf\n",Kf.Lap,Kf.SLap,Kf.Ext,Kf.LJ,Kf.Cont,Kf.El[0],Kf.El[1],Kf.El[2],Kf.Elong[0],Kf.Elong[1],Kf.Elong[2]);
  printf("External] Center %lf %lf %lf Rad %lf Hei %lf\n",Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],Nano->Rad,Nano->Height);
  fclose(FileToRead);
  return 0;
}
int Forces::ReSetNPart(int NewNPart){
  int OldNPart = pNPart();
  //if the number of allocated particle is less then the new particles returns
  Block[0].NPart = NewNPart;
  if(SetNPart(NewNPart)) return 1;
  if(VAR_IF_TYPE(SysAlloc,ALL_MC)){
    double *Dens2T = (double *)calloc(pNPart()*pNType(),sizeof(double));
    double *Dens3T = (double *)calloc(pNPart()*pNType(),sizeof(double));
    for(int p=0;p<OldNPart*pNType();p++){
      Dens2T[p] = Dens2[p];
      Dens3T[p] = Dens3[p];
    }
    //memcpy(Dens2,Dens2T,OldNPart*pNType());
    //memcpy(Dens3,Dens3T,OldNPart*pNType());
    free(Dens2);
    free(Dens3);
    Dens2 = (double *)calloc(pNAllocP()*pNType(),sizeof(double));
    Dens3 = (double *)calloc(pNAllocP()*pNType(),sizeof(double));
    // memcpy(Dens2T,Dens2,OldNPart*pNType());
    // memcpy(Dens3T,Dens3,OldNPart*pNType());
    for(int p=0;p<OldNPart*pNType();p++){
      Dens2[p] = Dens2T[p];
      Dens3[p] = Dens3T[p];
    }
    free(Dens2T);
    free(Dens3T);
    double *TempNrgPart = new double[OldNPart];
    for(int p=0;p<OldNPart;p++){
      TempNrgPart[p] = OldNrgBead[p];
    }
    delete[] OldNrgBead;
    OldNrgBead = new double[pNAllocP()];
    for(int p=0;p<OldNPart;p++){
      OldNrgBead[p] = TempNrgPart[p];
    }
    delete []TempNrgPart;
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_MD)){
    FORCES *TmpFm = new FORCES [OldNPart];
    for(int p=0;p<OldNPart;p++){
      for(int d=0;d<3;d++){
	TmpFm[p].Dir[d] = Fm[p].Dir[d];
	TmpFm[p].Ext[d] = Fm[p].Ext[d];
      }
    }
    free(Fm);
    Fm = (FORCES *)calloc(pNAllocP(),sizeof(FORCES));
    for(int p=0;p<OldNPart;p++){
      for(int d=0;d<3;d++){
	Fm[p].Dir[d] = TmpFm[p].Dir[d];
	Fm[p].Ext[d] = TmpFm[p].Ext[d];
      }
    }
    free(TmpFm);
  }
  return 0;
}
int Forces::ReSetNChain(int NewNChain){
  int OldNChain = pNChain();
  Block[0].NChain = NewNChain;
  if(SetNChain(NewNChain))return 1;
  if(VAR_IF_TYPE(SysAlloc,ALL_MC)){
    double *TempNrgCh = new double[OldNChain];
    for(int c=0;c<3*OldNChain;c++){
      TempNrgCh[c] = OldNrgCh[c];
    }
    delete[] OldNrgCh;
    OldNrgCh = new double[pNAllocC()];
    for(int c=0;c<3*OldNChain;c++){
      OldNrgCh[c] = TempNrgCh[c];
    }
    delete []TempNrgCh;
  }
  if(VAR_IF_TYPE(SysAlloc,ALL_MD)){
  }
  return 0;
}
void Forces::ReSetNPCh(int NewNPCh){
  Block[0].NPCh = NewNPCh;
  SetNPCh(NewNPCh);
}
void Forces::ReOpen(char *FName,int Bf){
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  Pc->Erase();
  delete Pc;
  Open(FName,Bf);
  Pc = new DomDec(Edge,pNPart(),sqrt(Kf.CutOff2));
  double Pos[3];
  for(int p=0;p<pNPart();p++){pPos(p,Pos);Pc->AddPart(p,Pos);}
  PrepareSys();
}
