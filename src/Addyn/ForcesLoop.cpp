#include "Forces.h"
/** main loop */
void Forces::ChooseSimMode(){
  if(VAR_IF_TYPE(SysShape,SYS_2D)){
    //CalcUpdate = &Forces::Sim2d;
  }
}
/** Solve a differential equation */
void Forces::Solve(){
  if(VAR_IF_TYPE(SysShape,SYS_2D))
    SolveLinksIterative();
    //SolveLinks();
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES))
    SolveLeaves();
  else if(VAR_IF_TYPE(SysShape,SYS_PORE))
    SolveLeaves();
  else if(VAR_IF_TYPE(SysShape,SYS_ROD))
    SolveRod();
}
void Forces::Dynamics(){
  if(VAR_IF_TYPE(SysShape,SYS_MD)){
    VelVerlet1();
  }
  if(VAR_IF_TYPE(SysShape,SYS_2D)){
    Wave();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_LEAVES)){
    ForceFieldLeaves();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_1D)){
    ForceFieldLine();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_3D)){
    ForceFieldBulk();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ROD)){
    ForceFieldRod();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_RIGID)){
    ForceFieldRigid();
    VelVerletRigid();
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MC)){
    if(VAR_IF_TYPE(CalcMode,CALC_NVT)){
      //for(int p=0;p<2*pNPCh();p++) 
      NInsertion += TryMove();
    }
    else if(VAR_IF_TYPE(CalcMode,CALC_NcVT)){
      NInsertion += TryMoveCh();
    }
    else if(VAR_IF_TYPE(CalcMode,CALC_mVT)){
      if(Mat->Casuale() < .5) 
	NRemoval += TryRemove();
      else 
	NInsertion += TryInsert();
    }
    else if(VAR_IF_TYPE(CalcMode,CALC_mcVT)){
      if(!VAR_IF_TYPE(CalcMode,CALC_CONF_BIAS)){
	if(Mat->Casuale() < .5) NRemoval += TryRemoveCh();
	else NInsertion += TryInsertCh();
      }
      else {
	if(Mat->Casuale() < .5) 
	  NRemoval += TryRemoveChBias();
	else 
	  NInsertion += TryInsertChBias();
      }
    }
  }
  else if(VAR_IF_TYPE(SysShape,SYS_ELECTRO)){
    if(VAR_IF_TYPE(CalcMode,CALC_NVT)){
      NInsertion += TryMove();
    }
    else if(VAR_IF_TYPE(CalcMode,CALC_mVT)){
      if(Mat->Casuale() < .5) 
	NRemoval += TryRemove();
      else 
	NInsertion += TryInsert();
    }
  }
  else if(VAR_IF_TYPE(SysShape,SYS_MD)){
    OldNrgSys = SumForcesMD();
    Pc->Erase();
    for(int p=0;p<pNPart();p++){
      Pc->AddPart(p,Pm[p].Pos);
    }
  }
  if(VAR_IF_TYPE(SysShape,SYS_MD)){
    ApplyTherm();
    VelVerlet2();
  }
  Time += pDeltat();
  IncrStep();
  SetTime(Time);
  Task();
}
void Forces::ExplorePepSize(){
  Shout("Explore pep size\n");
  if(SysShape != SYS_LEAVES) return ;
  char FName[120];
  double HeiMin = .15;
  double HeiMax = .4;
  double HeiStep = (HeiMax-HeiMin)/(double)NGrid;
  double AngMin = 5.;
  double AngMax = 35.;
  double AngStep = (AngMax-AngMin)/(double)NGrid;
  double SLapMin = 0.001;
  double SLapMax = 10.0;
  double SLapStep = 10.;
  for(double SLap=SLapMin;SLap<SLapMax;SLap*=SLapStep){
    sprintf(FName,"PepHeiElThickSLap%lf.dat",SLap);
    FILE *FWrite = fopen(FName,"w");
    for(double Hei=HeiMin;Hei<HeiMax;Hei+=HeiStep){
      for(double Angle=AngMin;Angle<AngMax;Angle+=AngStep){
	fprintf(stderr,"%lf %lf %lf\r",SLap,Hei,Angle);
	Nano->Height = Hei;
	Kf.SLap = SLap;
	Nano->Hamaker = Angle;
	Solve();
	int c = 0;
	double Min = pEdge(2);
	double Max = 0.;
	double xMin = 0.;
	double xMax = 0.;
	for(int p = c*pNPCh();p<(c+1)*pNPCh();p++){
	  if(Pm[p].Pos[2] > Max && Pm[p].Typ == 0){
	    Max = Pm[p].Pos[2];
	    xMax = Pm[p].Pos[0];
	  }
	}
	c = 1;
	for(int p = c*pNPCh();p<(c+1)*pNPCh();p++){
	  if(Pm[p].Pos[2] < Min && Pm[p].Typ == 0){
	    Min = Pm[p].Pos[2];
	    xMin = Pm[p].Pos[0];
	  }
	}
	fprintf(FWrite,"%lf %lf %lf %lf %lf\n",Hei/.2,Angle,(Min-Max)/.2,SLap,.5*(xMin+xMax)/.2);
      }
    }
    fclose(FWrite);
  }
  printf("\n");
}
void Forces::ExplorePepSize2d(){
  Shout("Explore pep size 2d\n");
  if(SysShape != SYS_2D) return ;
  SetEdge(.5*MIN(pEdge(0),pEdge(1)),3);
  char FName[120];
  double HeiMin = 0.2;
  double HeiMax = 4.2;
  double HeiStep = (HeiMax-HeiMin)/(double)NGrid;
  double AngMin = 0.;
  double AngMax = 80.;
  double AngStep = (AngMax-AngMin)/(double)NGrid;
  double InvNBin = 1./(double)NBin;
  SysFormat = VAR_SYS_TXVL;
  for(double Angle=AngMin;Angle<AngMax;Angle+=AngStep){
    for(double Hei=HeiMin;Hei<HeiMax;Hei+=HeiStep){
      fprintf(stderr,"%lf %lf\r",Hei,Angle);
      Nano->Height = Hei;
      Nano->Hamaker = Angle;
      Solve();
      sprintf(FName,"Sol2dHei%.1fAng%02d.dat",Hei,(int)Angle);
      Write(FName);
    }
  }
  printf("\n");
}
void Forces::ExploreDoubleMin(){
  Shout("Explore pep size\n");
  char String[120];
  if(SysShape != SYS_2D) return;
  if(pNNano() < 2) return;
  double DistMin = pNanoPos(1,0);
  double DistDelta = (pEdge(0)-pNanoPos(1,0))/(double)NGrid;
  double DeltaGrid = pEdge(0)/(double)nEdge[0];
  // DistDelta = floor(2.*DistDelta/DeltaGrid)*DeltaGrid;
  // SigErr(DistDelta <= 0.,"The grid step is too short for %d distances ",NGrid);
  char FName[60];
  FILE *FProf = fopen("DistThinMin.dat","w");
  for(int i=0;i<NGrid;i++){
    fprintf(stderr,"done %.2f %%\r",100.*i/(double)NGrid);
    IfFillMatrix = 1;
    Solve();
    sprintf(FName,"MinProf%03d.dat",i);
    FILE *FWrite = fopen(FName,"w");
    StringNano(String,0);
    fprintf(FWrite,String);
    StringNano(String,1);
    fprintf(FWrite,String);
    double Min = 100000.;
    double Dist = 0.;
    double zOld = 100000.;
    int IfContinue = 0;
    int nBin = 0;
    for(int p=nEdge[1]/2+1;p<pNPart();p+=nEdge[1]){
      double x = pPos(p,0) - pNanoPos(0,0);
      double z = pPos(p,2);
      if(pType(p) == 1) z = pNanoPos(0,2) + .5*Nano[0].Height;
      if( x > (pNanoPos(1,0) - pNanoPos(0,0) )) continue;
      fprintf(FWrite,"%lf %lf\n",x,z);
      if( x > .6*(pNanoPos(1,0) - pNanoPos(0,0) )) continue;
      if(Min > z){// && IfContinue){
	Min = z;
	Dist = x;
	nBin = p;
      }
      // if(z > zOld){IfContinue = 1;}
      // zOld = z;
    }
    fprintf(FProf,"%lf %lf %lf\n",Nano[1].Pos[0]-Nano[0].Pos[0],Dist,Min);
    fflush(FProf);
    fclose(FWrite);
    sprintf(FName,"Sol2dDist%.3f.dat",Nano[1].Pos[0]-Nano[0].Pos[0]);
    WriteTxvl(FName);
    Nano[1].Pos[0] += DistDelta;
    if(Nano[1].Pos[0] > pEdge(0)) return;
  }
  printf("\n");
  fclose(FProf);
}
void Forces::RunWidom(char *File2Open,int f){
  Shout("Widom on particles\n");
  int NInt = 10000;
  double NrgDiff = 0.;
  char File2Write[60];
  sprintf(File2Write,"WidomOut%05d.dat",f);
  ReOpen(File2Open,BF_PART);
  FILE *WidomOut = fopen(File2Write,"w");
  SigErr(WidomOut==NULL,"Cannot allocate %s\n",File2Write);
  for(int p=0;p<pNPart();p++){
    WidomRemove(&NrgDiff,p);
    fprintf(WidomOut,"%lf\n",NrgDiff);
  }
  fclose(WidomOut);
  sprintf(File2Write,"WidomIn%05d.dat",f);
  FILE *WidomIn = fopen(File2Write,"w");
  for(int p=0;p<NInt;p++){
    WidomInsert(&NrgDiff);
    fprintf(WidomIn,"%lf\n",NrgDiff);
  }
  fclose(WidomIn);
}
void Forces::RosenIn(FILE *WidomIn){
  Shout("Widom insertion biased/Rosenbluth\n");
  int NInt = 10000;
  double Weight = 0.;
  for(int p=0;p<NInt;p++){
    WidomBiasChIn(&Weight);
    fprintf(WidomIn,"%g\n",Weight);
  }
}
void Forces::RosenOut(FILE *WidomIn){
  Shout("Widom deletion biased/Rosenbluth\n");
  int NInt = 10000;
  double Weight = 0.;
  for(int c=0;c<pNChain();c++){
    WidomBiasChOut(&Weight,c);
    fprintf(WidomIn,"%g\n",Weight);
  }
}
void Forces::RunWidomChIn(char *File2Open,int f){
  Shout("Widom chain in\n");
  int NInt = 5*pNChain();
  double NrgDiff[3];
  char File2Write[60];
  ReOpen(File2Open,BF_PART);
  sprintf(File2Write,"WidomIn%05d.dat",f);
  FILE *WidomIn = fopen(File2Write,"w");
  //StudySys();
  for(int p=0;p<NInt;p++){
    IncrStep();
    WidomInsertCh(NrgDiff);
    fprintf(WidomIn ,"%lf\n",NrgDiff[2]);
  }
  fclose(WidomIn);
}
void Forces::RunWidomChOut(char *File2Open,int f){
  Shout("Widom chain out\n");
  double NrgDiff[3];
  char File2Write[60];
  sprintf(File2Write,"WidomOut%05d.dat",f);
  ReOpen(File2Open,BF_PART);
  FILE *WidomOut = fopen(File2Write,"w");
  CalcTotNrgCh();
  for(int c=0;c<pNChain();c++){
    WidomRemoveCh(NrgDiff,c);
    fprintf(WidomOut,"%lf\n",NrgDiff[2]);
  }
  fclose(WidomOut);
}
void Forces::CalcTotNrg(char *File2Open,int f){
  Shout("Calculating total energy\n");
  ReOpen(File2Open,BF_PART);
  double NrgNb = 0.;
  double NrgSpr = 0.;
  double NrgSpr2 = 0.;
  double NrgBend = 0.;
  for(int c=0;c<pNChain();c++){
    for(int p=c*pNPCh();p<(c+1)*pNPCh()-1;p++){
      double DistBA[3];
      double DistCB[3];
      double DistBA2 = 0.;
      double DistCB2 = 0.;
      double CosAngle = 0.;
      for(int d=0;d<3;d++){
	DistBA[d] = Pm[p].Pos[d] - Pm[p-1].Pos[d];
	DistBA[d] -= floor(DistBA[d]/(pEdge(d)) + .5)*pEdge(d);
	DistCB[d] = Pm[p+1].Pos[d] - Pm[p].Pos[d];
	DistCB[d] -= floor(DistCB[d]/(pEdge(d)) + .5)*pEdge(d);
	DistBA2 += SQR(DistBA[d]);
	DistCB2 += SQR(DistCB[d]);
	CosAngle += DistBA[d]*DistCB[d];
      }
      DistCB2 = sqrt(DistCB2);
      NrgSpr += .5*pkSpr()*SQR(DistCB2 - pSprRest());
      if(p == c*pNPCh()) continue;
      DistBA2 = sqrt(DistBA2);
      CosAngle /= (DistBA2*DistCB2);
      NrgBend += pkBen()*(1.-CosAngle);
   }
  }
  double NrgNano = 0.;
  for(int p=0;p<pNPart();p++){
    NrgNano += NanoNrg(p);
  }
  CalcTotNrgCh();
  printf("Nb %lf + Nano %lf = %lf Spr %lf + Ben %lf = %lf\n",OldNrgSys,NrgNano,OldNrgSys+NrgNano,NrgSpr,NrgBend,NrgSpr+NrgBend);
  fprintf(StatFile1,"%d %lf %lf %lf %lf\n",f,OldNrgSys,NrgSpr,NrgBend,NrgNano);
  fflush(StatFile1);
}
void Forces::CalcNrgPep(char *File2Open,int f){
  Shout("Calculating total energy\n");
  ReOpen(File2Open,BF_PART);
  AddDens(0,pNPart());
  OldNrgSys = SumDens(0,pNPart());  
  ClearDens();
  for(int b=0,NPep=0,cOff=0,pOff=0;b<pNBlock();cOff+=pNChain(b++)){
    if(!strncmp(Block[b].Name,"PEP",3)){
      continue;
    }
    for(int c=cOff;c<cOff+pNChain(b);c++,pOff+=pNPCh(b)){
      for(int p=pOff,link=0;p<MIN(pOff+pNPCh(b),pNPart());p++){
	AddDens(p,p+1);
      }
    }
  }
  double NrgPep = OldNrgSys - SumDens(0,pNPart());  
  printf("Nrg pep %lf\n",NrgPep);
  fprintf(StatFile1,"%d %lf %lf\n",f,NrgPep,OldNrgSys);
  fflush(StatFile1);
}
void Forces::CalcTens(char **argv,int *FilePos,int NFile){
  Shout("Calculating tension\n");
  AllocTens();
  char FName[60];
  for(int f=0;f<NFile;f++){
    fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
    ReOpen(argv[FilePos[f]],BF_PART);
    CalcTens();
    CalcDens();
    if( !(f%(NWrite)) ){
      for(int c=0;c<Tens.NComp;c++){
	if(VAR_IF_TYPE(Tens.CalcMode,CALC_2d))
	  sprintf(FName,"Tension2d%05dL%d.dat",f,c);
	if(VAR_IF_TYPE(Tens.CalcMode,CALC_3d))
	  sprintf(FName,"Tension3d%05dL%d.dat",f,c);
	WriteTens(FName,c,1./(double)NFile);
      }
    }
  }
  printf("\n");
}
void Forces::AvForces(char **argv,int *FilePos,int NFile){
  Shout("Calculating force field/average force\n");
  AllocTens();
  char FName[60];
  const int NBin = 120;
  double *Profile = (double *)calloc(3*NBin,sizeof(double));
  for(int f=0;f<NFile;f++){
    fprintf(stderr,"Elaborating file %s %.3f %%\r",argv[FilePos[f]],f/(double)NFile*100.);
    ReOpen(argv[FilePos[f]],BF_PART);
    CalcForcesDensFunc();
    FILE *FForce = fopen("ForceProfile.dat","w");
    double Force[3] = {0.,0.,0.};
    int cBin[3];
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	cBin[d] = (int)(Pm[p].Pos[d]*pInvEdge(d)*NBin);
	if(cBin[d] < 0 || cBin[d] >= NBin) continue;
	Profile[cBin[d]*3+d] += Fm[p].Dir[d];
	Force[d] += Fm[p].Dir[d];
	fprintf(FForce,"%d %g %g %g\n",p,Force[0],Force[1],Force[2]);
      }
    }
    printf("%g %g %g\n",Force[0],Force[1],Force[2]);
    fclose(FForce);
  }
  FILE *FProf = fopen("ForceProfile.dat","w");
  for(int b=0;b<NBin;b++){
    fprintf(FProf,"%d %lf %lf %lf\n",b,Profile[b*3],Profile[b*3+1],Profile[b*3+2]);
  }
  FILE *FField = fopen("ForceField.dat","w");
  double Delta = sqrt(Kf.CutOff2)/(double)NTab;
  for(int b=0;b<NTab;b++){
    double x = b*Delta;
    fprintf(FField,"%lf ",x);
    for(int t1=0;t1<pNType();t1++){
      for(int t2=0;t2<pNType();t2++){
	fprintf(FField,"%lf ",FTab[(b*pNType()+t1)*pNType()+t2]);
      }
    }
    fprintf(FField,"\n");
  }
  free(Profile);
  fclose(FProf);
  fclose(FField);
  printf("\n");
}
void Forces::Task(){
  if(VAR_IF_TYPE(SysAlloc,ALL_MC)){
    return;
    if((pStep()%NUpdate)!=0) return;
    OldNrgSys = 0.;
    double Pot[3];
    for(int c=0;c<pNChain();c++) OldNrgSys += CalcNrgCh(c,Pot);
    //fprintf(TempFile,"%d %d %lf \n",pStep(),pNChain(),NrgSys/(double)pNChain());
    //fprintf(StatFile1,"%d %d %lf \n",pStep(),pNChain(),OldNrgSys/(double)pNChain());
    //fflush(StatFile1);
    // int c = 0;//p%pNPCh();
    // double Pot[3];
    // CalcNrgCh(c,Pot);
    // for(int e=0;e<3;e++) OldNrgCh[c*3+e] = Pot[e];
    // double Dist[4];
    // TwoPartDist(0,pNPCh()-1,Dist);
    // fprintf(TempFile,"%d %lf %lf %lf %lf\n",pStep(),OldNrgCh[3*0],OldNrgCh[3*0+1],OldNrgCh[3*0+2],Dist[3]);
  }
  else if(VAR_IF_TYPE(SysAlloc,ALL_MD)){
    if( !(pStep()%NUpdate) ){
      double v2 = 0.;
      for(int p=0;p<pNPart();p++){
	for(int d=0;d<3;d++){
	  v2 += SQR(Pm[p].Vel[d]);
	}
      }
      //fprintf(StatFile1,"%d %lf %lf\n",pStep(),v2/(double)(3*pNPart()),OldNrgSys);
      //fflush(StatFile1);
    }
  }
}
void Forces::Trial(){
  CheckPairList();
}
void Forces::RunDynamics(){
  Shout("run dynamics\n");
  double NLoopSec = 0.;
  char FileName[60];
  DefForceParam();
  //MinimalNrg();
  for(int s=0;s<SimLimit;s++){
    Dynamics();
    //MinimalMD();
    if( !(s%NUpdate) ){
      CurrTime = time(NULL);
      NLoopSec += (s)/(double)(CurrTime-InitTime);
      double Temp = 0.;
      for(int p=0;p<pNPart();p++){
	for(int d=0;d<3;d++){
	  Temp += SQR(Pm[p].Vel[d]);
	}
      }
      Temp = Temp/(double)(3*pNPart());
      fprintf(stderr,"NPart %d loop/ms %.3g acc/step %.3f in/out %.4f T %.4f accomplished %.3f %% Nrg %lf\n",pNPart(),NLoopSec,(NRemoval+NInsertion)/(double)s,NInsertion/(double)NRemoval,Temp,s/(double)(SimLimit)*100.,OldNrgSys);
      fprintf(StatFile1,"%d %lf %lf\n",s,OldNrgSys,Temp);
    }
    if( !(s%(NWrite)) ){
      sprintf(FileName,"Trajectory%09d.dat",pStep());
      Write(FileName);
    }
  }
  printf("\n");
  CurrTime = time(NULL);
  NLoopSec += (SimLimit)/(double)(CurrTime-InitTime);
  double v2 = 0.;
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      v2 += SQR(Pm[p].Vel[d]);
    }
  }
  fprintf(StatFile1,"##NPart %d loop/ms %.3g acc/step %.3f in/out %.4f T %.4f accomplished Nrg %lf\r",pNPart(),NLoopSec,(NRemoval+NInsertion)/(double)pStep(),NInsertion/(double)NRemoval,v2/(double)(3*pNPart()),OldNrgSys);
}
double Forces::MinimalNrg(){
  double Pos[3];
  double Pot[3];
  double DistRel[4];
  OldNrgSys = 0.;
  Pc->Erase();
  for(int p=0;p<pNPart();p++){pPos(p,Pos);Pc->AddPart(p,Pos);}
  for(int p1=0;p1<pNPart();p1++){
    if(VAR_IF_TYPE(SysAlloc,ALL_FORCES)){
      for(int d=0;d<3;d++) Fm[p1].Dir[d] = 0.;
    }
    for(Pc->SetCurr(p1);Pc->IfCurr();Pc->NextCurr()){
      int p2 = Pc->p2Curr;
      if(p2 <= p1) continue;
      Pc->Dist2Curr(DistRel);
      if(DistRel[3] > Kf.CutOff2) continue;
      double Cons = Potential(DistRel[3],pType(p1),pType(p2),Pot);
      double InvDist = DistRel[3];
      if(VAR_IF_TYPE(SysAlloc,ALL_FORCES)){
	for(int d=0;d<3;d++){
	  Fm[p1].Dir[d] += Cons*DistRel[d]*InvDist;
	  Fm[p2].Dir[d] -= Cons*DistRel[d]*InvDist;
	}
      }
      OldNrgSys += Pot[0];
    }
  }
}
void Forces::MinimalMD(){
  double Sigma = sqrt(pTemp());
  double Pos[3];
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] += .5*Fm[p].Dir[d]*pDeltat();
      Pm[p].Pos[d] += Pm[p].Vel[d]*pDeltat();
      Pm[p].Pos[d] -= floor(Pm[p].Pos[d]*pInvEdge(d))*pEdge(d);
    }
  }
  MinimalNrg();
  //vv2
  double Temp = 0.;
  for(int p=0;p<pNPart();p++){
    for(int d=0;d<3;d++){
      Pm[p].Vel[d] += .5*Fm[p].Dir[d]*pDeltat();
      Temp += SQR(Pm[p].Vel[d]);
    }
  }
  //Andersen
  for(int p=0;p<pNPart();p++){
    if(Mat->Casuale() < pDeltat()){
      for(int d=0;d<3;d++){
	Pm[p].Vel[d] = Mat->Gaussiano(0.,Sigma);//Mat->Gaussiano(0.,Sigma);
      }
    }
  }
  Temp = Temp/(3*pNPart());
  double Norm = 1./sqrt(Temp);
  // for(int p=0;p<pNPart();p++){
  //   for(int d=0;d<3;d++){
  //     Pm[p].Vel[d] *= Norm;
  //   }
  // }
}
void Forces::MinimizeSol(){
  double *Sol2d = (double *)calloc(SQR(NEdge),sizeof(double));
  double *Count = (double *)calloc(SQR(NEdge),sizeof(double));
  double NInvEdge = 1./(double)NEdge;
  for(int p=0;p<pNPart();p++){
    if(pType(p) != 0) continue;
    int vx = (Pm[p].Pos[0]*NEdge*pInvEdge(0));
    if(vx < 0 || vx >= NEdge) continue;
    int vy = (Pm[p].Pos[1]*NEdge*pInvEdge(1));
    if(vy < 0 || vy >= NEdge) continue;
    Sol2d[vx*NEdge+vx] = Pm[p].Pos[2];
    Count[vx*NEdge+vx] = 1.;
  }
  for(int vx=0;vx<SQR(NEdge);vx++){
    double Norm = Count[vx] > 0. ?  1./Count[vx] : 1.;
    Sol2d[vx] *= Norm;
  }
  Create2d();
  int NInt = 1000;
  for(int n=0;n<NInt;n++){
    
  }
}
