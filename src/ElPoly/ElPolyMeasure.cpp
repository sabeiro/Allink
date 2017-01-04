#include "ElPoly.h"

int ElPoly::Diffusivity(){
  int NDim = 2;
  NFileTot = NFile[1] - NFile[0];
  int NChain = 0;
  for(int b=0;b<pNBlock();b++){
    if(!strcmp(Block[b].Name,"ADDED")) continue;
    if(!strcmp(Block[b].Name,"CHOL")) continue;
    NChain += Block[b].NChain;
  }
  double *FileAverage = (double *)calloc(NDim*NChain*NFileTot,sizeof(double));
  double *ChainDiff = (double *)calloc(NFileTot,sizeof(*ChainDiff));
  double *ChainOrigin = (double *)calloc(NFileTot,sizeof(*ChainDiff));
  double *Count = (double *)calloc(NFileTot,sizeof(*Count));
  double *Time = (double *)calloc(NFileTot,sizeof(*Time));
  double CmInit[3] = {pCm(0),pCm(1),pCm(2)};
  double InitTime = pTime();
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO))return 1;
    BfDefChain();
    Time[f] = pTime() - InitTime;
    if(Time[f]<=0.) Time[f] = (double)f;
    for(int c=0;c<NChain;c++){
      int ff = f - NFile[0];
      for(int d=0;d<NDim;d++){
	FileAverage[(c*NFileTot+ff)*NDim+d] = Ch[c].Pos[d];
      }
    }
  }
  printf("\n");
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,FileAverage,NDim*NChain*NFileTot, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    for(int f=0;f<NFileTot;f++){
      for(int c=0;c<NChain;c++){
	for(int d=0;d<NDim;d++){
	  ChainOrigin[f] += 
	    SQR((FileAverage[(c*NFileTot)*NDim+d]-FileAverage[(c*NFileTot+f)*NDim+d]));//-(pCm(d)-CmInit[d]) );
	}
      }
    }
    for(int f=0;f<NFileTot;f++){
      for(int ff = f; ff<NFileTot;ff++){
	int fff = ff - f;
	for(int c=0;c<NChain;c++){
	  for(int d=0;d<NDim;d++){
	    ChainDiff[fff] += 
	      SQR((FileAverage[(c*NFileTot+ff)*NDim+d]-FileAverage[(c*NFileTot+f)*NDim+d]));//-(pCm(d)-CmInit[d]) );
	  }
	  Count[fff] += 1.;
	}
	//printf("%d %d %d %lf\n",f,ff,fff,ChainDiff[fff]);
      }
    }
    char FileName[60];
    sprintf(FileName,"Diffusivity.dat");
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    fprintf(FileToWrite,"%s",cSystem);
    fprintf(FileToWrite,"#Time ChainOrigin OriginSum DiffSum ChainDiff Count\n");
    for(int f=0;f<NFileTot;f++){
      double Inv = Count[f] > 0. ? Count[f] : 1.;
      fprintf(FileToWrite,"%lf %.5g %.5g %.5g\n",Time[f],ChainOrigin[f]/(double)(NChain),ChainDiff[f]/Inv,ChainDiff[f]/(Inv*Time[f]) );
    }
#ifdef OMPI_MPI_H
  }
#endif
  free(FileAverage);
  free(ChainDiff);
  free(ChainOrigin);
  free(Count);
  free(Time);
  return 0;
}
void ElPoly::DiffSlab(int NSlab){
  int NBin = NFile[1] - NFile[0];
  double Volume1 = 0.;
  double HeightCyl = pEdge(CNorm);
  double LenMin = MIN(pEdge(CLat1),pEdge(CLat2));
  double LenMax = MAX(pEdge(CLat1),pEdge(CLat2));
  double BoxRad = .5*MIN(pEdge(CLat1),pEdge(CLat2));
  double InvBoxRad = 1./BoxRad;
  double *OldPos = (double *)calloc(3*pNPart(),sizeof(double));
  double *Diff = (double *)calloc(NBin*pNPart(),sizeof(double));
  double *Count = (double *)calloc(NSlab,sizeof(double));
  int *PartOr = (int *)calloc(pNPart(),sizeof(int));
  double dr[3];
  for(int p = 0;p<pNPart();p++){
    double dr2 = 0.;
    for(int d=0;d<3;d++){
      OldPos[3*p+d] = Pm[p].Pos[d];
      dr[d] = Nano->Pos[d] - Pm[p].Pos[d];
      dr[d] -= floor(dr[d]*pInvEdge(d))*pEdge(d);
    }
    dr2 = sqrt(SQR(dr[CLat1]) + SQR(dr[CLat2]));
    int vr = (int)(dr2/BoxRad*NSlab);
    PartOr[p] = vr;
    if(vr < 0 || vr >= NSlab) continue;
    Count[vr] += 1.;
  }
  char FName[60];
  if(1==1){//particles position and slab circles
    FILE *FRepr = fopen("SystemSlab.dat","w");
    for(int p=0;p<pNPart();p++){
      fprintf(FRepr,"%lf %lf %lf 0\n",Pm[p].Pos[CLat1],Pm[p].Pos[CLat2],.5);
    }
    for(int v=0;v<NSlab;v++){
      double Dist = v*BoxRad/(double)NSlab;
      for(int i=0;i<100;i++){
	double ang = i*.01*2.*M_PI;
	double x = Dist*cos(ang) + Nano->Pos[CLat1];
	double y = Dist*sin(ang) + Nano->Pos[CLat2];
	fprintf(FRepr,"%lf %lf %lf 2\n",x,y,.51);
      }
    }
    fclose(FRepr);
    double Norm = 1./(double)pNPart();
    for(int v=0;v<NSlab;v++){
      double Dist = v*BoxRad/(double)NSlab;
      sprintf(FName,"MeanSqDispl%.3f.dat",Dist);
      FILE *FWrite = fopen(FName,"w");
      fclose(FWrite);
    }
  }
  int NOut = 0;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART));
    double OutRatio = 100.*NOut/(double)(f*pNPart());
    fprintf(stderr,"step %d)/%d fraction of number of beads out of the slab %lf %%\r",f,NBin,OutRatio);
    double Rad2Mean = 0.;
    for(int p = 0;p<pNPart();p++){
      int vr1 = PartOr[p];
      if(vr1 < 0 || vr1 >= NSlab) continue;
      double dr2 = 0.;
      for(int d=0;d<3;d++){
	dr[d] = Nano->Pos[d] - Pm[p].Pos[d];
	dr[d] -= floor(dr[d]*pInvEdge(d))*pEdge(d);
      }
      dr2 = sqrt(SQR(dr[CLat1]) + SQR(dr[CLat2]));
      int vr = (int)(dr2*InvBoxRad*NSlab);
      if(vr != vr1){
	NOut++;
      }
      double Rad2 = 0.;
      for(int d=0;d<2;d++){
	double Dist = Pm[p].Pos[d] - OldPos[3*p+d];
	Dist -= floor(Dist*pInvEdge(d))*pEdge(d);
	Rad2 += SQR(Dist);
      }
      // vr or vr1?
      Diff[f*NSlab+vr1] += Rad2;
      Rad2Mean += Rad2;
    }
    if(1==1){//particle position and mean displacement
      Rad2Mean = pNPart()/Rad2Mean;
      sprintf(FName,"DisplPos%05d.dat",f);
      FILE *FWrite = fopen(FName,"w");
      fprintf(FWrite,"# l(%lf %lf %lf) d[part]\n",pEdge(CLat1),pEdge(CLat2),.5);
      for(int p=0;p<pNPart();p++){
	double Rad2 = 0.;
	for(int d=0;d<2;d++){
	  double Dist = Pm[p].Pos[d] - OldPos[3*p+d];
	  Dist -= floor(Dist*pInvEdge(d))*pEdge(d);
	  Rad2 += SQR(Dist);
	}
	int vr1 = PartOr[p];
	double Vel = Rad2*Rad2Mean;
	fprintf(FWrite,"{x(%lf %lf %lf) v(%lf %lf %lf)\n",Pm[p].Pos[CLat1],Pm[p].Pos[CLat2],Vel,Vel,Vel,Vel);
      }
      fclose(FWrite);
    }
  }
  printf("\n");
  double DisplMean = 0.;
  //normalize
  for(int v=0;v<NSlab;v++){
    Count[v] = Count[v] > 0. ? 1./Count[v] : 1.;
    for(int b=0;b<NBin;b++){
      Diff[b*NSlab+v] *= Count[v];
    }
  }
  for(int v=0;v<NSlab;v++){
    DisplMean += Diff[(NBin-1)*NSlab+v] ;
    double Dist = v*BoxRad/(double)NSlab;
    sprintf(FName,"MeanSqDispl%.3f.dat",Dist);
    FILE *FWrite = fopen(FName,"w");
    for(int b=1;b<NBin;b++){
      fprintf(FWrite,"%d %lf\n",b,Diff[b*NSlab+v]);
    }
    fclose(FWrite);
  }
  DisplMean /= (double) NSlab;
  if(1==1){//slab mean displ
    sprintf(FName,"SlabMeadDispl.dat");
    FILE *FWrite = fopen(FName,"w");
    fprintf(FWrite,"# l(%lf %lf %lf) d[part]\n",pEdge(CLat1),pEdge(CLat2),.5);
    for(int sx=0;sx<100;sx++){
      for(int sy=0;sy<100;sy++){
	double x = sx*.01*pEdge(CLat1);
	double y = sy*.01*pEdge(CLat2);
	double Rad2 = sqrt(SQR(x-Nano->Pos[CLat1])+SQR(y-Nano->Pos[CLat2]));
	int vr = (int)(Rad2*InvBoxRad*NSlab);
	if(vr < 0 || vr >= NSlab) continue;
	double Vel = Diff[(NBin-1)*NSlab+vr];
	fprintf(FWrite,"{x(%lf %lf %lf) v(%lf %lf %lf)\n",x,y,Vel,Vel,Vel,Vel);
      }
    }
    fclose(FWrite);
  }
  free(OldPos);
  free(Diff);
  free(Count);
  free(PartOr);
}
#include <Cubo.h>
typedef DdLinkedList DomDec;
void ElPoly::PairCorr(int NBin,int NDim){
  int NAlloc = NBin;
  double CutOff = 4.;
  if(NDim == 2) NAlloc = NBin*NBin;
  double *Prob = (double *)calloc(NAlloc,sizeof(double));
  double Edge[3] = {pEdge(0),pEdge(1),pEdge(2)};
  double Norm = 0.;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    // DomDec *Pc = new DomDec(Edge,pNPart(),1.);
    // for(int p=0;p<pNPart();p++) Pc->AddPart(p,Pm[p].Pos);
    // int NeiList[27];
    double DistRel[4];
    for(int p1=0;p1<pNPart();p1++){
      for(int p2=p1+1;p2<pNPart();p2++){
	// int NNei = Pc->GetNei(Pm[p1].Pos,NeiList);
	// for(int i=0;i<NNei;i++){
	//   int c1 = NeiList[i];
	//   for(Pc->SetCounters(c1);Pc->IfItCell(c1);Pc->IncrCurr(c1)){
	// 	int p2 = Pc->ItCell(c1);
	// 	if(p2 < p1) continue;
	if(!TwoPartDist(p1,p2,DistRel,CutOff)) continue;
	int r = (int)(DistRel[3]/CutOff*NBin);
	if(r < 0 || r >= NBin) continue;
	Prob[r] += 1.;
	Norm += 1.;
      }
    }
    // }
    // delete Pc;
  }
  printf("\n");
  Norm *= NFile[1] - NFile[0];
  double rInv = 1.*CutOff/(double)NBin;
  FILE *FPair = fopen("PairCorrelation.dat","w");
  fprintf(FPair,"%lf %lf \n",0.,Prob[0]/Norm);
  for(int r=1;r<NBin;r++){
    double rr = r*CutOff/(double)NBin;
    Prob[r] /= Norm*(CUBE(rr)-CUBE(rr-rInv));
    fprintf(FPair,"%lf %lf \n",rr,Prob[r]);
  }
  fclose(FPair);
  free(Prob);
}
int ElPoly::PairCorrelationF(int NBin,int How){
  double **Plot;
  Plot = (double **)calloc(NBin,sizeof(double));
  for(int i=0;i<NBin;i++){
    *(Plot+i) = (double *)calloc(NBin,sizeof(double));
  }
  double *dDensity;
  dDensity = (double *)calloc(NBin,sizeof(double));
  double *dDensity1;
  dDensity1 = (double *)calloc(NBin,sizeof(double));
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //printf("Opening: %s\n",cFile[f]);
    if(OpenRisk(cFile[f],BF_CHAIN))return 1;
    if(How == 0 || How == 1)
      PairCorrelation(dDensity,NBin,How,NChType);
    else if (How == 2)
      PairCorrelationRound(Plot,NBin,NChType);
    else if (How == 3)
      PairCorrelationSquare(Plot,NBin,NChType);
    else if (How == 4)
      PairCorrelationPep(Plot,NBin,NChType);
    for(int v=0;v<NBin;v++){
      //printf("%d %lf\n",v,dDensity[v]);
      dDensity1[v] += dDensity[v];
    }
  }
  printf("\n");
  if(How == 0 || How == 1){
    FILE *FileToWrite = fopen("PairCorrelation.dat","w");
    fprintf(FileToWrite,"# RadDist Density \n");
    for(int v=0;v<NBin;v++){
      dDensity1[v] /= NFile[1] - NFile[0];
      fprintf(FileToWrite,"%lf %lf\n",(double)v/(double)NBin*pEdge(3),dDensity1[v]);
    }
    fclose(FileToWrite);
  }
  else {
    double Max=0.;
    double InvNBin=1./(double)NBin;
    for(int v=0;v<NBin;v++){
      for(int vv=0;vv<NBin;vv++){
	//	Plot[v][vv] /= DUE_PI*( QUAD((Gen->Edge[3]/InvNBin))*(QUAD((v+1)) - QUAD(v)) );
	if(Max < Plot[v][vv])
	  Max = Plot[v][vv];
      }
    }
    FILE *FileToWrite = fopen("PairCorrelation.xvl","w");
    fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",2.*pEdge(3),2.*pEdge(3),pEdge(3),NBin,ChooseDraw(EL_QUAD1));
    int vRef = (int)(NBin/2.);
    for(int v=0;v<NBin;v++){
      for(int vv=0;vv<NBin;vv++){
	if(Plot[v][vv] > 0.){
	  if(How == 2){
	    double PosX = pEdge(3) + InvNBin*(double)v*pEdge(3)*cos((double)vv*InvNBin*DUE_PI);
	    double PosY = pEdge(3) + InvNBin*(double)v*pEdge(3)*sin((double)vv*InvNBin*DUE_PI);
	    fprintf(FileToWrite,"{x(%lf %lf %lf)}\n",
		    PosX,PosY,Plot[v][vv]/Max*pEdge(3));
	  }
	  else{
	    double PosX = pEdge(CLat1)*.5 + (v - vRef)*InvNBin*pEdge(CLat1);
	    double PosY = pEdge(CLat2)*.5 + (vv - vRef)*InvNBin*pEdge(CLat2);
	    fprintf(FileToWrite,"{x(%lf %lf %lf)}\n",
		    PosX,PosY,Plot[v][vv]/Max*pEdge(3));
	  }
	  //printf("%d %d %lf %lf %lf\n",v,vv,PosX,PosY,Plot[v][vv]);
	}
      }
    }
    fclose(FileToWrite);
  }
  free(Plot);
  free(dDensity);
  free(dDensity1);
  return 0;
}
int ElPoly::ScatteringF(int NBin,int How){
  double **Plot;
  Plot = (double **)calloc(NBin,sizeof(double));
  for(int i=0;i<NBin;i++){
    *(Plot+i) = (double *)calloc(NBin,sizeof(double));
  }
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    //printf("Opening: %s\n",cFile[f]);
    if(OpenRisk(cFile[f],BF_CHAIN))return 1;
    Scattering2d(Plot,NBin,NChType);
  }
  printf("\n");
  double Max=0.;
  double InvNBin=1./(double)NBin;
  for(int v=0;v<NBin;v++){
    for(int vv=0;vv<NBin;vv++){
      Plot[v][vv] = POS(Plot[v][vv]);
      if(Max < Plot[v][vv])
  	Max = Plot[v][vv];
    }
  }
  Max = 1./Max;
  FILE *FileToWrite = NULL;
  if(How == 0)
    FileToWrite = fopen("Scattering.xvl","w");
  else 
    FileToWrite = fopen("Scattering2.xvl","w");    
  fprintf(FileToWrite,"#l(%lf %lf %lf) v(%d) d[%s]\n",2.*pEdge(3),2.*pEdge(3),pEdge(3),NBin,ChooseDraw(EL_QUAD));
  int vRef = NBin/2;
  for(int q=0;q<4;q++){
    double Signx = !(q%2) ? 1. : -1.;
    double Signy =  q > 1 ? 1. : -1.;
    for(int vx=0;vx<NBin;vx++){
      for(int vy=0;vy<NBin;vy++){
	//if(Plot[v][vv] > 0.)
	{
	  double PosX = Signx*vx*InvNBin*pEdge(CLat1);
	  double PosY = Signy*vy*InvNBin*pEdge(CLat2);
	  //printf("%d %d %lf %lf %lf\n",v,vv,PosX,PosY,Plot[v][vv]);
	  fprintf(FileToWrite,"{x(%lf %lf %lf) v(%lf %lf %lf)}\n",
		  PosX,PosY,Plot[vx][vy]*.001,Plot[vx][vy],0.,0.);
	}
      } 
    }
  }
  fclose(FileToWrite);
  free(Plot);
  return 0;
}
//obsolete?
int ElPoly::NChainPSquareF(){
  int NTimes = 20;
  double *ChainPArea = (double *)calloc(3*NTimes,sizeof(double));
  double MeanSigma = 0.;
  int SubDiv[3] = {NTimes,NTimes,1};
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return 1;
    //    OrderPos();
    SubDiv[CLat1] = 1;
    SubDiv[CLat2] = 1;
    for(int t=0;t<NTimes;t++){
      if( (t%2)==0) 
	SubDiv[CLat1] += 1;
      if( (t%2-1)==0)
	SubDiv[CLat2] += 1;
      double NEdge = (double)(SubDiv[CLat1]*SubDiv[CLat2]);
      double *Plot = (double *)calloc((int)NEdge,sizeof(double));
      double Area = pEdge(CLat1)*pEdge(CLat2)/NEdge;
      double SumX = 0.;
      double SumX2 = 0.;
      for(int p=0;p<pNPart();p++){
	int vx = (int)(pPos(p,CLat1)*pInvEdge(CLat1)*SubDiv[CLat1]);
	if(vx > SubDiv[CLat1]) continue;
	int vy = (int)(pPos(p,CLat2)*pInvEdge(CLat2)*SubDiv[CLat2]);
	if(vy > SubDiv[CLat2]) continue;
	Plot[vx*SubDiv[CLat2] + vy] += 1.;
      }
      for(int n=0;n<NEdge;n++){
	SumX += Plot[n]/Area;
	SumX2 += QUAD((Plot[n]/Area));
	Plot[n] /= Area;
      }
      int Valori = 20;
      double *Intervalli = (double *)calloc(Valori,sizeof(double));
      int IfNorm = 0;
      MOMENTI m1 = Mat->Distribuzione(Plot,(int)NEdge,Intervalli,Valori,IfNorm);
      ChainPArea[3*t+0] = NEdge;
      MeanSigma += m1.Uno;
      double Media = SumX/NEdge;
      ChainPArea[3*t+1] = sqrt((SumX2 - Media*Media*NEdge)/(NEdge-1));
      printf("%lf  %lf-%lf\n",m1.Uno,ChainPArea[3*t+1],m1.Due);
      ChainPArea[3*t+2] = Area*QUAD(m1.Due)/QUAD(m1.Uno);
      //QUAD(Area)*ChainPArea[4*t+2]/QUAD(ChainPArea[4*t+1]);
      free(Plot);
      free(Intervalli);
    }
  }
  char *FileName = (char *)calloc(60,sizeof(char));
  sprintf(FileName,"ChainPArea%.0fKappa%.0fRho%.0f.dat",pchiN(),pkappaN(),prho());
  FILE *FileToWrite = fopen(FileName,"w");
  char cSystem[STRSIZE];
  SysDef(cSystem);
  fprintf(FileToWrite,"# %s",cSystem);
  fprintf(FileToWrite,"#%lf\n",MeanSigma/(double)(NFile[1]-NFile[0]));
  fprintf(FileToWrite,"# NEdge SDeviation Area*Dev/NEdge\n");
  for(int t=0;t<NTimes;t++)
    fprintf(FileToWrite,"%lf %lf %lf\n",ChainPArea[3*t+0],ChainPArea[3*t+1],ChainPArea[3*t+2]);
  fclose(FileToWrite);
  free(FileName);
  free(ChainPArea);
  return 0;
}
int ElPoly::SpectrumF(int NSample){
  int NHalf = NSample/2;
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  double InvNSample = 1./(double)NSample;
  double InvNSample2 = 1./(double)SQR(NSample);
  double *Plot = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotUp = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotDown = (double *) calloc(SQR(NSample),sizeof(double));
  double *CountUp = (double *) calloc(SQR(NSample),sizeof(double));
  double *CountDown = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotBS = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotA = (double *)calloc(SQR(NSample), sizeof(double));
  double *PlotAUp = (double *)calloc(SQR(NSample), sizeof(double));
  double *PlotADown = (double *)calloc(SQR(NSample), sizeof(double));
  double *PlotThin = (double *)calloc(SQR(NSample), sizeof(double));
  double *Count = (double *) calloc(NSample*NSample,sizeof(double));
#ifdef  USE_FFTW
  fftw_complex *out = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_complex *in = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_plan plan = fftw_plan_dft_2d(NSample,NSample,
				     in, out,FFTW_FORWARD,FFTW_PATIENT);
  fftw_complex *outUp = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_complex *inUp = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_plan planUp = fftw_plan_dft_2d(NSample,NSample,
				     inUp, outUp,FFTW_FORWARD,FFTW_PATIENT);
  fftw_complex *outDown = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_complex *inDown = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_plan planDown = fftw_plan_dft_2d(NSample,NSample,
				     inDown, outDown,FFTW_FORWARD,FFTW_PATIENT);
					//FFTW_MEASURE);
  fftw_complex *outThin = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_complex *inThin = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_plan planThin = fftw_plan_dft_2d(NSample,NSample,
				     inThin, outThin,FFTW_FORWARD,FFTW_PATIENT);
					//FFTW_MEASURE);
#endif
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    // memset(Plot,0.,NSample*NSample*sizeof(double));
    // memset(Count,0.,NSample*NSample*sizeof(double));
    for(int s=0;s<SQR(NSample);s++){
      Plot[s] = 0.;
      Count[s] = 0.;
      PlotUp[s] = 0.;
      CountUp[s] = 0.;
      PlotDown[s] = 0.;
      CountDown[s] = 0.;
    }
    if(OpenRisk(cFile[f],BF_CHAIN))return 0;
    for(int p=0;p<pNPart();p++){
      int sx = (int)(pPos(p,CLat1)*pInvEdge(CLat1)*NSample);
      int sy = (int)(pPos(p,CLat2)*pInvEdge(CLat2)*NSample);
      if(sx < 0 || sx >= NSample) continue;
      if(sy < 0 || sy >= NSample) continue;
      Plot[sx*NSample+sy] += pPos(p,CNorm);// - pCm(CNorm);
      Count[sx*NSample+sy] += 1.;
      //if(Pm[p].Typ != 1) continue;
      if(VAR_IF_TYPE(Ch[pChain(p)].Type,CHAIN_UP)){
	PlotUp[sx*NSample+sy] += pPos(p,CNorm);
	CountUp[sx*NSample+sy] += 1.;
      }
      else {
	PlotDown[sx*NSample+sy] += pPos(p,CNorm);
	CountDown[sx*NSample+sy] += 1.;
      }
    }
    double MeanThick = 0.;
    for(int s=0;s<SQR(NSample);s++){
      Plot[s] /= Count[s] > 0. ? Count[s] : 1.;
      PlotUp[s] /= CountUp[s] > 0. ? CountUp[s] : 1.;
      PlotDown[s] /= CountDown[s] > 0. ? CountDown[s] : 1.;
      MeanThick += CountUp[s] - CountDown[s];
    }
    MeanThick /= (double)SQR(NSample);
    //InterBSpline2D(Plot,PlotBS,NSample,NSample);
    //InterBSpline2D(PlotBS,Plot,NSample,NSample);
#ifdef USE_FFTW
    if(1==1){
      for(int s=0;s<SQR(NSample);s++){
	in[s][0] = Plot[s];
	in[s][1] = 0.;
	inThin[s][0] = PlotUp[s]-PlotDown[s];
	inThin[s][1] = 0.;
	inUp[s][0] = PlotUp[s];
	inUp[s][1] = 0.;
	inDown[s][0] = PlotDown[s];
	inDown[s][1] = 0.;
      }
      fftw_execute(plan);
      fftw_execute(planThin);
      fftw_execute(planUp);
      fftw_execute(planDown);
      for(int sx=0;sx<NSample;sx++){
	int ssx = sx < NHalf ? sx + NHalf : sx - NHalf;
	for(int sy=0;sy<NSample;sy++){
	  int ssy = sy < NHalf ? sy + NHalf : sy - NHalf;
	  PlotA[ssx*NSample+ssy] += (SQR(out[sx*NSample+sy][0])+SQR(out[sx*NSample+sy][1]))*SQR(InvNSample2);
	  PlotThin[ssx*NSample+ssy] += (SQR(outThin[sx*NSample+sy][0])+SQR(outThin[sx*NSample+sy][1]))*SQR(InvNSample2);
	  PlotAUp[ssx*NSample+ssy] += (SQR(outUp[sx*NSample+sy][0])+SQR(outUp[sx*NSample+sy][1]))*SQR(InvNSample2);
	  PlotADown[ssx*NSample+ssy] += (SQR(outDown[sx*NSample+sy][0])+SQR(outDown[sx*NSample+sy][1]))*SQR(InvNSample2);
	}
      }
    }
    else {
      int NMax = NSample;
      double *st = Plot;
      double *sw = PlotA;
      double dNMax = 1./(double)NMax;
      int NHalf = (int)(NMax/2.);
      for(int kx=-NHalf;kx<NHalf;kx++){
	double qx = kx*dNMax;
	for(int ky=-NHalf;ky<NHalf;ky++){
	  double qy = ky*dNMax;
	  double Re2=0.,Im2=0.;
	  double Re1=0.,Im1=0.;
	  for(int lx=0;lx<NMax;lx++){
	    for(int ly=0;ly<NMax;ly++){
	      double Arg = 2.*M_PI*(lx*qx + ly*qy);
	      double cy = cos(Arg);
	      double sy = sin(Arg);
	      Re1 += st[lx*NMax+ly]*cy;
	      Im1 += st[lx*NMax+ly]*sy;
	    }
	  }
	  int kkx = kx + NHalf;
	  int kky = ky + NHalf;
	  sw[kkx*NMax+kky] += SQR(Re1*dNMax) + SQR(Im1*dNMax);
	}
      }
    }
    if(1==0){
      FILE *FMidplane = fopen("Midplane.dat","w");
      FILE *FThickness = fopen("Thickness.dat","w");
      for(int sx=0;sx<NSample;sx++){
	double x = sx*pEdge(CLat1)*InvNSample;
	for(int sy=0;sy<NSample;sy++){
	  double y = sy*pEdge(CLat2)*InvNSample;
	  fprintf(FMidplane,"%lf %lf %lf \n",x,y, Plot[sx*NSample+sy]);
	  fprintf(FThickness,"%lf %lf %lf \n",x,y,PlotUp[sx*NSample+sy]-PlotDown[sx*NSample+sy]);
	}
      }
      fclose(FMidplane);
      fclose(FThickness);
    }
#else
    Mat->Spettro2d(Plot,PlotA,NSample);
#endif //FFTW3_H
  }
  printf("\n");
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,PlotA,NSample*NSample, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    char *NomeFile = (char*)calloc(60,sizeof(char));
    sprintf(NomeFile,"Spectrum.xvl",pchiN(),pkappaN(),prho());
    FILE *FileToWrite = fopen(NomeFile,"w");
    fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[squaremesh]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NSample);
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
	double x = pEdge(CLat1)*sx*InvNSample;
	double y = pEdge(CLat2)*sy*InvNSample;
	fprintf(FileToWrite,"{x(%lf %lf %lf) v(%lf %lf %lf)}\n",x,y,PlotA[sx*NSample+sy]*InvNFile,fabs(log10(PlotA[sx*NSample+sy]*InvNFile)),0.,0.);
      }
    }
    printf("Write\n");
    fclose(FileToWrite);
    sprintf(NomeFile,"Spectrum.dat",pchiN(),pkappaN(),prho());
    FileToWrite = fopen(NomeFile,"w");
    fprintf(FileToWrite,"# 1) q^2 2) h(q)^2*L_xL_y | q_x=s2pi/(L_x) NSample %d\n",NSample);
    for(int kx=0;kx<NSample;kx++){
      int kx1 = kx - NSample/2;
      double qx  = 2.*M_PI*kx*pInvEdge(CLat1);
      double qxt = 2.*M_PI*kx1*pInvEdge(CLat1);
      for(int ky=0;ky<NSample;ky++){
	int ky1 = ky - NSample/2;
	double qy  = 2.*M_PI*ky*pInvEdge(CLat2);
	double qyt = 2.*M_PI*ky1*pInvEdge(CLat2);
	double qq = SQR(qxt) + SQR(qyt);
	double Spe = PlotA[kx*NSample+ky]*InvNFile*pEdge(CLat1)*pEdge(CLat2);
	double SpeUp = PlotAUp[kx*NSample+ky]*InvNFile*pEdge(CLat1)*pEdge(CLat2);
	double SpeDown = PlotADown[kx*NSample+ky]*InvNFile*pEdge(CLat1)*pEdge(CLat2);
	double SpeThin = PlotThin[kx*NSample+ky]*InvNFile*pEdge(CLat1)*pEdge(CLat2);
	fprintf(FileToWrite,"%lf %g %g %g %g\n",qq,Spe,SpeUp,SpeDown,SpeThin);
      }
    }
    fclose(FileToWrite);
    free(NomeFile);
#ifdef OMPI_MPI_H
  }
#endif
  free(Plot);
  free(PlotA);
  free(PlotBS);
  free(Count);
#ifdef USE_FFTW
  fftw_destroy_plan(plan);
#endif
  return 0;
}
void ElPoly::Midplane(int NSample){
  int NHalf = NSample/2;
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  double InvNSample = 1./(double)NSample;
  double *PlotBS  = (double *) calloc(SQR(NSample),sizeof(double));
  double *Count = (double *) calloc(NSample*NSample,sizeof(double));
  double *Plot  = (double *) calloc(SQR(NSample),sizeof(double));
  double *PlotA  = (double *) calloc(SQR(NSample),sizeof(double));
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    // memset(Plot,0.,NSample*NSample*sizeof(double));
    // memset(Count,0.,NSample*NSample*sizeof(double));
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
    	Plot[sx*NSample+sy] = 0.;
      }
    }
    if(OpenRisk(cFile[f],NBackFold))return;
    for(int p=0;p<pNPart();p++){
      int sx = (int)(pPos(p,CLat1)*pInvEdge(CLat1)*NSample);
      int sy = (int)(pPos(p,CLat2)*pInvEdge(CLat2)*NSample);
      if(sx < 0 || sx >= NSample) continue;
      if(sy < 0 || sy >= NSample) continue;
      Plot[sx*NSample+sy] += pPos(p,CNorm);// - pCm(CNorm);
      Count[sx*NSample+sy] += 1.;
    }
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
	Plot[sx*NSample+sy] /= Count[sx*NSample+sy] > 0. ? Count[sx*NSample+sy] : 1.;
	PlotA[sx*NSample+sy] += Plot[sx*NSample+sy];
      }
    }
    // InterBSpline2D(Plot,PlotBS,NSample,NSample);
    // InterBSpline2D(PlotBS,Plot,NSample,NSample);
    char FName[60];
    sprintf(FName,"Midplane%05d.dat",f);
    FILE *FMid = fopen(FName,"w");
    for(int sy=0;sy<NSample;sy++){
      for(int sx=0;sx<NSample;sx++){
	double x = pEdge(CLat1)*sx*InvNSample;
	double y = pEdge(CLat2)*sy*InvNSample;
	fprintf(FMid,"%lf %lf %lf\n",x,y,Plot[sx*NSample+sy]);
      }
    }
    fclose(FMid);
  }
  printf("\n");
  FILE *FMid = fopen("average_middplane.dat","w");
  for(int sy=0;sy<NSample;sy++){
    for(int sx=0;sx<NSample;sx++){
      double x = pEdge(CLat1)*sx*InvNSample;
      double y = pEdge(CLat2)*sy*InvNSample;
      fprintf(FMid,"%lf %lf %lf\n",x,y,PlotA[sx*NSample+sy]*InvNFile);
    }
  }
  fclose(FMid);
  free(Plot);
  free(Count);
  free(PlotBS);
}
void ElPoly::SpectrumMidplane(int NSample){
  // FILE *Ciccia = fopen("Ciccia.dat","w");
  // for(int p=0;p<pNPart();p++){
  //   fprintf(Ciccia,"%lf %lf\n",SQR(pPos(p,0)-.5*pEdge(0))+SQR(pPos(p,1)-.5*pEdge(1)),pPos(p,2));
  // }
  // fclose(Ciccia);
  // return;
  int NHalf = NSample/2;
  double InvNSample = 1./(double)NSample;
  double *Plot = (double *) calloc(NSample*NSample,sizeof(double));
  double *Count = (double *) calloc(NSample*NSample,sizeof(double));
  double *PlotA = (double *) calloc(NSample*NSample,sizeof(double));
  int NPoint = (int)(NSample*sqrt(2.));
  double *Density = (double *) calloc(NPoint,sizeof(double));
#ifdef  USE_FFTW
  fftw_complex *out = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_complex *in = (fftw_complex *)fftw_malloc(SQR(NSample)*sizeof(fftw_complex));
  fftw_plan plan1 = fftw_plan_dft_2d(NSample,NSample,
  				     in, out,FFTW_FORWARD,FFTW_PATIENT);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
    	Plot[sx*NSample+sy] = 0.;
    	Count[sx*NSample+sy] = 0.;
      }
    }
    FILE *Ciccia = fopen(cFile[f],"r");
    char cLine[256];
    for(int p=0;!(fgets(cLine,STRSIZE,Ciccia)==NULL);p++){
      double x = 0., y = 0., z = 0.;
      sscanf(cLine,"%lf %lf %lf\n",&x,&y,&z);
      int v  = (int)(x*NSample*pInvEdge(CLat1));//(pEdge(0)+.4));
      int vv = (int)(y*NSample*pInvEdge(CLat2));//(pEdge(1)+.4));
      if(v  < 0 || v  >= NSample) continue;
      if(vv < 0 || vv >= NSample) continue;
      Plot[v*NSample+vv] += z;
      Count[v*NSample+vv] += 1.;
    }
    for(int s=0;s<NSample;s++){
      for(int ss=0;ss<NSample;ss++){
	Plot[s*NSample+ss] /= Count[s*NSample+ss] > 0. ? Count[s*NSample+ss] : 1.;
      }
    }    
    fclose(Ciccia);
    if(1==0){
      for(int sx=0;sx<NSample;sx++){
	for(int sy=0;sy<NSample;sy++){
	  in[sx*NSample+sy][0] = Plot[sx*NSample+sy];
	  in[sx*NSample+sy][1] = 0.;
	}
      }
      fftw_execute(plan1);
      for(int sx=0;sx<NSample;sx++){
	int sx1 = sx + NHalf;
	if(sx1 >= NSample) sx1 -= NSample;
	for(int sy=0;sy<NSample;sy++){
	  int sy1 = sy + NHalf;
	  if(sy1 >= NSample) sy1 -= NSample;
	  PlotA[sx*NSample+sy] += (SQR(out[sx*NSample+sy][0])+SQR(out[sx*NSample+sy][1]))*SQR(InvNSample);
	  //PlotA[sx*NSample+sy] += Plot[sx][sy];
	}
      }
    }    
    else {
      int NMax = NSample;
      double *st = Plot;
      double *sw = PlotA;
      double dNMax = 1./(double)NMax;
      int NHalf = (int)(NMax/2.);
      for(int kx=-NHalf;kx<NHalf;kx++){
	double qx = kx*dNMax;
	for(int ky=-NHalf;ky<NHalf;ky++){
	  double qy = ky*dNMax;
	  double Re2=0.,Im2=0.;
	  double Re1=0.,Im1=0.;
	  for(int lx=0;lx<NMax;lx++){
	    for(int ly=0;ly<NMax;ly++){
	      double Arg = 2.*M_PI*(lx*qx + ly*qy);
	      double cy = cos(Arg);
	      double sy = sin(Arg);
	      Re1 += st[lx*NMax+ly]*cy;
	      Im1 += st[lx*NMax+ly]*sy;
	    }
	  }
	  int kkx = kx + NMax/2;
	  int kky = ky + NMax/2;
	  sw[kkx*NMax+kky] += SQR(Re1*dNMax) + SQR(Im1*dNMax);
	}
      }
    }
  }
#endif // USE_FFTW
  printf("\n");
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,PlotA,NSample*NSample, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    char *NomeFile = (char*)calloc(60,sizeof(char));
    sprintf(NomeFile,"SpectrumChi%.0fKappa%.0fRho%.0f.xvl",pchiN(),pkappaN(),prho());
    FILE *FileToWrite = fopen(NomeFile,"w");
    // fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[squaremesh]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NSample);
    double DivInv = 1./(double)(NFile[1]-NFile[0]);
    for(int s=0;s<NSample;s++){
      for(int ss=0;ss<NSample;ss++){
	fprintf(FileToWrite,"%lf %lf %lf\n",pEdge(CLat1)*s*InvNSample,pEdge(CLat2)*ss*InvNSample,PlotA[s*NSample+ss]*DivInv);
	// fprintf(FileToWrite,"{x(%lf %lf %lf) v(%lf %lf %lf)}\n",pEdge(CLat1)*s*InvNSample,pEdge(CLat2)*ss*InvNSample,PlotA[s*NSample+ss]*DivInv,0.,0.,log10(PlotA[s*NSample+ss]*DivInv));
      }
    }
    fclose(FileToWrite);
    sprintf(NomeFile,"SpectrumChi%.0fKappa%.0fRho%.0f.dat",pchiN(),pkappaN(),prho());
    FileToWrite = fopen(NomeFile,"w");
    double qHalf = QUAD(pEdge(CLat1)*(NHalf)*InvNSample)+QUAD(pEdge(CLat2)*(NHalf)*InvNSample);
    double dx = pEdge(CLat1)*InvNSample;
    double dy = pEdge(CLat2)*InvNSample;
    for(int kx=0;kx<NSample;kx++){
      int kx1 = kx - NSample/2;
      double qx = kx*InvNSample;
      double qxt = kx1*InvNSample;
      for(int ky=0;ky<NSample;ky++){
	int ky1 = ky - NSample/2;
	double qy = ky*InvNSample;
	double qyt = ky1*InvNSample;
	double qq = SQR(qx/dx - .5/dx) + SQR(qy/dy - .5/dy);
	//double qq = SQR(qx) + SQR(qy);
	fprintf(FileToWrite,"%lf %g\n",qq,PlotA[kx*NSample+ky]*DivInv);
	//fprintf(FileToWrite,"%lf %g\n",QUAD(pEdge(CLat1)*(s-NHalf)*InvSample)+QUAD(pEdge(CLat2)*(ss-NHalf)*InvSample),PlotA[s*NSample+ss]*DivInv);
      }
    }
    // for(int p=1;p<NPoint;p++)
    //   fprintf(FileToWrite,"%d %lf\n",p,Density[p]*DivInv);
    fclose(FileToWrite);
    free(NomeFile);
#ifdef OMPI_MPI_H
  }
#endif
  free(Density);
  free(PlotA);
  free(Plot);
}
void ElPoly::HeaderAverage(int nNano){
  if(nNano >= pNNano()){
    printf("The specified nanoparticle doesn't exist\n");
    nNano = 0;
    // return;
  }
  int NDim = 2;
  int NumFile = NFile[1]-NFile[0];
  double *Area = (double *)calloc(NumFile,sizeof(double));
  double *NanoPos = (double *)calloc(6*NumFile,sizeof(double));
  double *NanoDiff = (double *)calloc(2*NumFile,sizeof(double));
  double *NanoCount = (double *)calloc(NumFile,sizeof(double));
  double *NanoDist = (double *)calloc(2*NumFile,sizeof(double));
  double *Time = (double *)calloc(NumFile,sizeof(double));
  double NanoInit[3] = {Nano->Pos[0],Nano->Pos[1],Nano->Pos[2]};
  double AreaAv = 0.;
  double AreaErr = 0.;
  double InitTime = pTime();
  //indentation
  double InitPos = pPos(0,CNorm) - pPos(1,CNorm) - MAX(Nano[0].Rad,Nano[0].Height);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    FILE *File2Read;
    if((File2Read = fopen(cFile[f],"r"))==0){
      printf("The file is missing\n");
      return ;
    }
    ReadHeader(File2Read);
    fclose(File2Read);
    AreaAv  += pEdge(CLat1)*pEdge(CLat2);
    AreaErr += SQR(pEdge(CLat1)*pEdge(CLat2));
    Area[f] = pEdge(CLat1)*pEdge(CLat2);
    Time[f] = pTime() - InitTime;
    for(int d=0;d<3;d++) NanoPos[f*6+d] = Nano[nNano].Pos[d];
    NanoPos[f*6+3] = Nano[nNano].Rad;
    NanoPos[f*6+4] = Nano[nNano].Height;
    NanoPos[f*6+5] = Nano[nNano].Area;
    for(int ff=f;ff>=0;ff--){
      int fff = f-ff;
      NanoDiff[fff*2] += SQR(NanoPos[f*6+0]-NanoPos[ff*6+0]) + SQR(NanoPos[f*6+1]-NanoPos[ff*6+1]);
      NanoCount[fff] += 1.;
    }
    NanoDiff[f*2+1] = SQR(NanoPos[f*6+0]-NanoInit[0]) + SQR(NanoPos[f*6+1]-NanoInit[1]);
    if(pNNano() > 1){
      //indentation
      NanoDist[f*2] = InitPos - pNanoPos(0,CNorm) + pNanoPos(1,CNorm) + MAX(Nano[nNano].Rad,Nano[nNano].Height);
      double Dist = 0.;
      for(int d=0;d<3;d++){
	double Pos1 = pNanoPos(0,d) - floor(pNanoPos(0,d)*pInvEdge(d))*pEdge(d);
	Dist += SQR(pNanoPos(0,d)-pNanoPos(1,d));
      }
      Dist = sqrt(Dist);
      //two np dist
      NanoDist[f*2+1] = Dist;
    }
  }
  AreaAv /= (double)(NFile[1] - NFile[0]);
  AreaErr = sqrt( AreaErr - (NFile[1] - NFile[0])*SQR(AreaAv))/(double)(NFile[1] - NFile[0]);
  printf("\n");
  printf("Area %lf \\pm %lf\n",AreaAv,AreaErr);
  FILE *File2Write = fopen("AreaTime.dat","w");
  fprintf(File2Write,"# Time Area %lf pm %lf\n",AreaAv,AreaErr);
  char FName[120];
  sprintf(FName,"NanoDiffR%.1fS%.1fH%.1f.dat",Nano[nNano].Rad,Nano[nNano].Hamaker,Nano[nNano].Height);
  FILE *FNano = fopen(FName,"w");
  fprintf(FNano,"# t <msd> msd Dz NpDist rad hei\n");
  for(int f=0;f<NumFile;f++){
    NanoDiff[f*2] /= NanoCount[f] > 0. ? NanoCount[f] : 1.;
    fprintf(File2Write,"%lf %lf\n",Time[f],Area[f]);
    fprintf(FNano,"%lf %lf %lf %lf %lf %lf %lf\n",Time[f],NanoDiff[f*2  ],NanoDiff[f*2+1],NanoDist[f*2],NanoDist[f*2+1],NanoPos[f*6+3],NanoPos[f*6+4]);
  }
  fclose(File2Write);
  fclose(FNano);
  free(NanoDiff);
  free(NanoDist);
  free(Area);
  free(Time);
}
//only one block
void ElPoly::WidomOut(char *NrgFile,int NBin){
  int NChains = 1;
  Block[0].NChain -= NChains;
  WriteXvt("SnapOut.dat");
  FILE *FileLast = fopen("LastChain.dat","w");
  for(int p=pNPart()-(NChains)*pNPCh();p<pNPart();p++){
    fprintf(FileLast,"%lf %lf %lf %lf %lf %lf %d\n",pPos(p,0),pPos(p,1),pPos(p,2),pVel(p,0),pVel(p,1),pVel(p,2),pType(p));
  }
  fclose(FileLast);
}
//only one block
void ElPoly::WidomOut(){
  int NChains = 1;
  Block[0].NChain -= NChains;
  char FileName[60];
  sprintf(FileName,"SnapOut%05d.dat",pNChain()-1);
  WriteXvt(FileName);
  for(int c=0;c<pNChain()-1;c++){
    SwapChain(c,pNChain()-1);
    sprintf(FileName,"SnapOut%05d.dat",c);
    WriteXvt(FileName);
  }
}
//the system has an additional fake chain
void ElPoly::WidomIn(){
  int NFile = 1000;
  char FileName[60];
  SetNPart(pNPart()+pNPCh());
  SetNChain(pNChain()+1);
  Block[0].NChain += 1;
  int p1 = pNPCh()*(pNChain()-1);
  double Var = sqrt(SQR(pReOverCutOff())/(double)(pNPCh()-1)/3.)/2.;
  for(int f=0;f<NFile;f++){
    for(int d=0;d<3;d++){
      SetPos(p1,d,Mat->Casuale()*pEdge(d));
    }
    for(int p=p1+1;p<pNPart();p++){
      for(int d=0;d<3;d++){
	SetPos(p1,d,pPos(p-1,d) + Mat->Gaussiano(0.,Var));
      }
      if( p%pNPCh() >= Block[0].Asym )
	SetType(p,1);
    }
    sprintf(FileName,"SnapIn%05d.dat",f);
    Write(FileName);
  }
}
//the system has an additional fake chain
void ElPoly::WidomIn(char *NrgFile,int NBin){
  int NChains = 1;
  Block[0].NChain += NChains;
  WriteXvt("SnapIn.dat");
  int NAddPart = NChains*pNPCh();
  FILE *FileLast = fopen("SnapIn.dat","a");
  PART *Pn = (PART *)calloc(NAddPart,sizeof(PART));
  double GaussVar = sqrt(SQR(pReOverCutOff())/(double)(pNPCh()-1)/3.)*2./3.;
  for(int c=0;c<NChains;c++){
    int p1 = c*pNPCh();
    for(int d=0;d<3;d++){
      SetPos(p1,d,Mat->Casuale()*pEdge(d));
    }
    for(int p=1;p<pNPCh();p++){
      for(int d=0;d<3;d++){
	Pn[p+p1].Pos[d] = Pn[p-1+p1].Pos[d] + Mat->Gaussiano(0,GaussVar);
      }
      Pn[p+p1].Typ = 0;
      if(p >= Block[0].Asym){
	Pn[p+p1].Typ = 1;
      }
    }
  }
  for(int p=0;p<NAddPart;p++){
    fprintf(FileLast,"%lf %lf %lf %lf %lf %lf %d\n",Pn[p].Pos[0],Pn[p].Pos[1],Pn[p].Pos[2],Pn[p].Vel[0],Pn[p].Vel[1],Pn[p].Vel[2],Pn[p].Typ);
  }
  free(Pn);
  fclose(FileLast);
}
void ElPoly::End2EndDistr(char *OutFile){
  char FileName[256];
  for(int f=NFile[0];f<NFile[1];f++){
    sprintf(FileName,"%s%04d.dat",OutFile,f);
    FILE *File2Write = fopen(FileName,"w");
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO))return;
    for(int c=0;c<pNChain();c++){
      int p1 = c*pNPCh();
      int p2 = c*pNPCh() + pNPCh() - 1;
      double End2End = 0.;
      for(int d=0;d<3;d++){
	End2End += SQR(pPos(p1,d) - pPos(p2,d));
      }
      fprintf(File2Write,"%lf ",End2End);
      fprintf(File2Write,"\n");
    }
    fclose(File2Write);
  }
}
void ElPoly::Decoupling(int What){
  NFileTot = NFile[1] - NFile[0];
  double *Angles = (double *)calloc((NFile[1]-NFile[0])*pNChain(),sizeof(double));
  double *Count = (double *)calloc((NFile[1]-NFile[0]),sizeof(double));
  double *AnglesNano = (double *)calloc((NFile[1]-NFile[0])*pNNano(),sizeof(double));
  double *Time = (double *)calloc((NFile[1]-NFile[0]),sizeof(double));
  double *AngleDiff = (double *)calloc(pNChain(),sizeof(double));
  double *AngleDiff2 = (double *)calloc(pNChain(),sizeof(double));
  Vettore Ax0(1.,0.,0.);
  Vettore Ax2(0.,0.,1.);
  double InvPhob = Block[0].Asym > 0 ? 1./(double)Block[0].Asym : 1.;
  double InvPhil = Block[0].Asym > 0 ? 1./(double)(pNPCh() - Block[0].Asym) : 1.;
  double DistRel[4];
  Vettore ChDir(0.,0.,0.);
  Vettore NanoAx(0.,0.,0.);
  double Average = 0.;
  double Variance = 0.;
  double InitTime = pTime();
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    int ff = f - NFile[0];
    if(OpenRisk(cFile[f],BF_CHAIN))return ;
    Time[f] = pTime() - InitTime;
    for(int n=0;n<pNNano();n++){
      for(int d=0;d<3;d++){
	NanoAx.Set(Nano[n].Axis[d],d);
      }
      NanoAx.Set(0.,CNorm);
      AnglesNano[ff*pNNano()+n] = Ax0.Angle(&NanoAx);
    }
    for(int c=0,b=0;c<pNChain();c++){
      if(What == 0){
	ChDir.Set(Ch[c].Dir[0],0);
	ChDir.Set(Ch[c].Dir[1],1);
	double Angle = Ax0.Angle(&ChDir);
	if(isnan(Angle)) Angles[ff*pNChain()+c] = Angles[(ff-1)*pNChain()+c];
	else Angles[ff*pNChain()+c] = Ax0.Angle(&ChDir);
      }
      else if(What == 1){
	Angles[ff*pNChain()+c] = TwoPartDist(c*pNPCh(),(c+1)*pNPCh()-1,DistRel);
      }
      Average += Angles[ff*pNChain()+c];
      Variance += SQR(Angles[ff*pNChain()+c]);
    }
  }
  printf("\n");
  double NTot = (double)(pNChain()*(NFile[1]-NFile[0]));
  Average /= NTot > 0. ? NTot : 1.;
  double PreVar = Variance;
  Variance = (Variance - SQR(Average)*NTot)/(NTot-1);
  if(isnan(Variance)){
    printf("Variance is not a number %lf %lf %lf\n",PreVar,Average,NTot);
    Variance = 1.;
  }
  for(int f=0;f<NFileTot;f++){
    for(int c=0;c<pNChain();c++){
      AngleDiff2[f] += (Angles[f*pNChain()+c]-Average)*(Angles[0*pNChain()+c]-Average);
    }
    for(int ff = f+1; ff<NFileTot;ff++){
      int fff = ff - f;
      for(int c=0;c<pNChain();c++){
	AngleDiff[fff]  += (Angles[ff*pNChain()+c]-Average)*(Angles[f*pNChain()+c]-Average);
	Count[fff] += 1.;
      }
    }
  }
  char FileName[256];
  if(What == 0)
    sprintf(FileName,"DirDecoupling.dat");
  else if(What == 1)
    sprintf(FileName,"E2EDecoupling.dat");    
  FILE *FileToWrite = fopen(FileName,"w");
  char cSystem[STRSIZE];
  SysDef(cSystem);
  //  fprintf(FileToWrite,"%s",cSystem);
  fprintf(FileToWrite,"#Time AngleProd\n");
  for(int f=1;f<NFileTot;f++){
    double Inv = Count[f] > 0. ? Count[f] : 1.;
    fprintf(FileToWrite,"%lf %lf %lf\n",Time[f],AngleDiff[f]/(Inv*Variance),AngleDiff2[f]/(Variance*pNChain()) );
  }
  fclose(FileToWrite);
  for(int n=0;n<pNNano();n++){
    sprintf(FileName,"NanoAngles%d.dat",n);
    FILE *FNanoAng = fopen(FileName,"w");
    for(int f=0;f<NFileTot;f++)
      fprintf(FNanoAng,"%lf %lf \n",Time[f],AnglesNano[f*pNNano()+n]*360./DUE_PI);
    fclose(FNanoAng);
  }
  // FileToWrite = fopen("RawAngles.dat","w");
  // for(int f=0;f<NFileTot;f++){
  //   for(int c=0;c<pNChain();c++){
  //     fprintf(FileToWrite,"%lf\n",Angles[f*pNChain()+c]);
  //   }
  // }
  // fclose(FileToWrite);
  free(Count);
  free(Time);
  free(Angles);
  free(AngleDiff);
  free(AngleDiff2);
}
/** Subdivide the simulation box in SubDiv[CLat1]xSubDiv[CLat2]
    squares and calculate the average number of chains and its
    stadard deviation. The compressibility is in the limit
    lim_l\to\infty = l^2<N_c^2>/<N_c>^2
    where l is the dimension of the patch and NTimes different 
    patches size are taken in account. It works in NVT.
*/
void ElPoly::AreaCompr(int NSample){
  int SubDiv[3] = {NSample,NSample,1};
  int NEdge = (SubDiv[CLat1]*SubDiv[CLat2]);
  //mean StdDev Area*StdDev/Mean^2
  double *ChainPArea = (double *)calloc(3*NSample,sizeof(double));
  double *Plot = (double *)calloc(NEdge,sizeof(double));
  double *Area = (double *)calloc(NSample,sizeof(double));
  double *NDiv = (double *)calloc(NSample,sizeof(double));
  char FileName[256];
  //for the simulation in NPtT a smaller box is used
  double Edge[3] = {pEdge(0)-1.,pEdge(1)-1.,pEdge(2)};
  double InvEdge[3] = {1./Edge[0],1./Edge[1],1./Edge[2]};
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_CHAIN);
    for(int t=0;t<NSample;t++){
      SubDiv[CLat1] = NSample-t;
      SubDiv[CLat2] = NSample-t;
      if(t == NSample -1) SubDiv[CLat2] = 2;
      NDiv[t] = SubDiv[CLat1]*SubDiv[CLat2];
      Area[t] = Edge[CLat1]*Edge[CLat2]/(SubDiv[CLat1]*SubDiv[CLat2]);
      for(int n=0;n<NEdge;n++) Plot[n] = 0.;
      for(int c=0;c<pNChain();c++){
	double Posx = pChPos(c,CLat1) - floor(pChPos(c,CLat1)*pInvEdge(CLat1))*pEdge(CLat1);
	double Posy = pChPos(c,CLat2) - floor(pChPos(c,CLat2)*pInvEdge(CLat2))*pEdge(CLat2);
	int vx = (int)(Posx*InvEdge[CLat1]*SubDiv[CLat1]);
	int vy = (int)(Posy*InvEdge[CLat2]*SubDiv[CLat2]);
	if(vx < 0 || vx >= SubDiv[CLat1]) continue;
	if(vy < 0 || vy >= SubDiv[CLat2]) continue;
	Plot[vx*NSample+vy] += 1.;
      }
      double SumX = 0.;
      double SumX2 = 0.;
      double Count = 0.;
      for(int sx = 0;sx<SubDiv[CLat1];sx++){
	for(int sy = 0;sy<SubDiv[CLat2];sy++){
	  int n = sx*NSample+sy;
	  Plot[n] /= Area[t];
	  SumX += Plot[n];
	  Count += 1.;
	}
      }
      double NInv = Count > 0. ? 1./Count : 1.;
      double MeanVal = SumX*NInv;
      for(int sx = 0;sx < SubDiv[CLat1];sx++){
	for(int sy = 0;sy < SubDiv[CLat2];sy++){
	  int n = sx*NSample+sy;
	  SumX2 += SQR(Plot[n] - MeanVal);
	}
      }
      ChainPArea[3*t  ] += MeanVal;
      ChainPArea[3*t+1] += SumX2/(Count - 1.);
      ChainPArea[3*t+2] += 1.;
      if(1==0){//writes the single distribution files
	sprintf(FileName,"Distr%.0fKappa%02d.dat",pkappaN(),t);
	FILE *FileToWrite = fopen(FileName,"a");
	for(int l1=0;l1<SubDiv[CLat1];l1++){
	  for(int l2=0;l2<SubDiv[CLat2];l2++){
	    fprintf(FileToWrite,"%lf\n",Plot[l1*NSample+l2]);
	  }
	}
	fclose(FileToWrite);
      }
    }
  }
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,ChainPArea,3*NSample, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    free(Plot);
    //--------normalizing--------------------
    double Norm = 1./(double)(NFile[1]-NFile[0]);
    for(int t=0;t<NSample;t++){
      ChainPArea[3*t  ] *= Norm;
      ChainPArea[3*t+1] *= Norm;
      ChainPArea[3*t+2] *= Norm;
    }
    sprintf(FileName,"ChainPArea%.0fRho%.0fChiN%.0fKappa.dat",prho(),pchiN(),pkappaN());
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    fprintf(FileToWrite,"# %s",cSystem);
    fprintf(FileToWrite,"# Area MeanChPArea SDeviation Area*Dev/MeanChPArea\n");
    for(int t=0;t<NSample;t++){
      double Compr = Area[t]*ChainPArea[3*t+1]/SQR(ChainPArea[3*t  ]);
      fprintf(FileToWrite,"%lf %lf %lf %lf\n",Area[t],(ChainPArea[3*t  ]),ChainPArea[3*t+1],Compr);
      //fprintf(FileToWrite,"%lf %lf %lf %lf\n",ChainPArea[3*t+0],Area[t]*ChainPArea[3*t+1],ChainPArea[3*t+1],ChainPArea[3*t+2]);
    }
    fclose(FileToWrite);
#ifdef OMPI_MPI_H
  }
#endif
  free(ChainPArea);
}
/**
   Divides the box size in small patches and calculate the thinning for every patch.
   For every patch size are stored the mean values and variance of the thickness.
   As output it is printed the area per patch, the correspondent thickness, the standard deviation and the elastic coupling modulus.
*/
void ElPoly::ElasticCoupling(int NSample){
  double Edge[3] = {pEdge(0)-1.,pEdge(1)-1.,pEdge(2)};
  int SubDiv[3] = {NSample,NSample,1};
  //mean StdDev Area*StdDev/Mean^2
  double *ElCoup = (double *)calloc(3*(NSample+1),sizeof(double));
  double *Area = (double *)calloc((NSample+1),sizeof(double));
  double *NDiv = (double *)calloc((NSample+1),sizeof(double));
  double MeanSigma = 0.;
  char FileName[256];
  //for the simulation in NPtT a smaller box is used
  double InvEdge[3] = {1./Edge[0],1./Edge[1],1./Edge[2]};
  double *PosUp   = (double *)calloc(SQR(NSample),sizeof(double));
  double *PosDown = (double *)calloc(SQR(NSample),sizeof(double));
  double *CountUp  = (double *)calloc(SQR(NSample),sizeof(double));
  double *CountDown= (double *)calloc(SQR(NSample),sizeof(double));
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_CHAIN);
    for(int t=0;t<NSample+1;t++){
      //defining the partitions
      // SubDiv[CLat1] = (int)pow(2.,t);
      // SubDiv[CLat2] = (int)pow(2.,t);
      SubDiv[CLat1] = NSample - t;
      SubDiv[CLat2] = NSample - t;
      if(t == NSample-1) SubDiv[CLat2] = 2;
      if(t == NSample){
	SubDiv[CLat1] = 1;
	SubDiv[CLat2] = 1;
      }	
      NDiv[t] = SubDiv[CLat1]*SubDiv[CLat2];
      Area[t] = Edge[CLat1]*Edge[CLat2]/(double)(SubDiv[CLat1]*SubDiv[CLat2]);
      for(int n=0;n<SQR(NSample);n++){
	PosUp[n] = 0.;
	PosDown[n] = 0.;
	CountUp[n] = 0.;
	CountDown[n] = 0.;
      }
      for(int p=0;p<pNPart();p++){
	if(pType(p) != 1) continue;
	double x = pPos(p,CLat1) - floor(pPos(p,CLat1)*pInvEdge(CLat1))*pEdge(CLat1);
	double y = pPos(p,CLat2) - floor(pPos(p,CLat2)*pInvEdge(CLat2))*pEdge(CLat2);
	int sx = (int)(x*InvEdge[CLat1]*SubDiv[CLat1]);
	if(sx < 0 || sx >= SubDiv[CLat1]) continue;
	int sy = (int)(y*InvEdge[CLat2]*SubDiv[CLat2]);
	if(sy < 0 || sy >= SubDiv[CLat2]) continue;
	if(VAR_IF_TYPE(Ch[pChain(p)].Type,CHAIN_UP)){
	  PosUp[sx*NSample+sy] += pPos(p,CNorm);
	  CountUp[sx*NSample+sy] += 1.;
	}
	else {
	  PosDown[sx*NSample+sy] += pPos(p,CNorm);
	  CountDown[sx*NSample+sy] += 1.;
	}
      }
      double SumX = 0.;
      double SumX2 = 0.;
      double Count = 0.;
      for(int sx = 0;sx<SubDiv[CLat1];sx++){
	for(int sy = 0;sy<SubDiv[CLat2];sy++){
	  int n = sx*NSample+sy;
	  PosUp[n] /= CountUp[n] > 1. ? CountUp[n] : 1.;
	  PosDown[n] /= CountDown[n] > 1. ? CountDown[n] : 1.;
	  double Val = (PosUp[n]-PosDown[n]);
	  SumX += Val;
	  Count += 1.;
	  ElCoup[3*t+2] += 1.;
	}
      }
      double NInv = Count > 0. ? 1./Count : 1.;
      double MeanVal = SumX*NInv;
      for(int sx = 0;sx < SubDiv[CLat1];sx++){
	for(int sy = 0;sy < SubDiv[CLat2];sy++){
	  int n = sx*NSample+sy;
	  double Val = (PosUp[n]-PosDown[n]);
	  SumX2 += SQR(Val - MeanVal);
	}
      }
      ElCoup[3*t  ] += MeanVal;
      ElCoup[3*t+1] += SumX2/(Count - 1.);
      if(1==0){//writes the single distribution files
	sprintf(FileName,"Distr%.0fKappa%02d.dat",pkappaN(),t);
	FILE *FileToWrite = fopen(FileName,"a");
	for(int l1=0;l1<SubDiv[CLat1];l1++){
	  for(int l2=0;l2<SubDiv[CLat2];l2++){
	    int n = l1*NSample+l2;
	    fprintf(FileToWrite,"%lf\n",(PosUp[n]-PosDown[n]));
	  }
	}
	fclose(FileToWrite);
      }
    }
  }
 #ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,ElCoup,3*NSample, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    //--------normalizing--------------------
    double Norm = 1./(double)(NFile[1]-NFile[0]);
    for(int t=0;t<NSample+1;t++){
      ElCoup[3*t  ] *= Norm;
      ElCoup[3*t+1] *= Norm;
      ElCoup[3*t+2] *= Norm;
    }
    ElCoup[3*NSample+1] = 0.;
    sprintf(FileName,"ElasticCoupling%.0fRho%.0fChiN%.0fKappa.dat",prho(),pchiN(),pkappaN());
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    fprintf(FileToWrite,"# %s",cSystem);
    fprintf(FileToWrite,"# Area Thickness ThickStdDev Area*ThickStdDev\n");
    for(int t=0;t<NSample+1;t++){
      fprintf(FileToWrite,"%lf %lf %lf %lf \n",Area[t],ElCoup[3*t  ],ElCoup[3*t+1],Area[t]*ElCoup[3*t+1]/SQR(ElCoup[3*t  ]));
    }
    fclose(FileToWrite);
#ifdef OMPI_MPI_H
  }
#endif
  free(ElCoup);
  free(PosUp);
  free(PosDown);
  free(CountUp);
  free(CountDown);
}
/** 
    Elastic coupling modulus calculation for a fixed box size simulation.
 */
void ElPoly::ElasticCouplingNVT(){
  //mean StdDev Area*StdDev/Mean^2
  double ElCoup[3] = {0.,0.,0.};
  char FileName[256];
  char cSystem[STRSIZE];
  //for the simulation in NPtT a smaller box is used
  sprintf(FileName,"ElCoupNVT%.0fRho%.0fChiN%.0fKappa.dat",prho(),pchiN(),pkappaN());
  FILE *FileToWrite = fopen(FileName,"w");
  SysDef(cSystem);
  fprintf(FileToWrite,"# %s",cSystem);
  fprintf(FileToWrite,"# step Thickness ThickStdDev\n");
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_CHAIN);
    double PosUp = 0.;
    double PosDown = 0.;
    double CountUp = 0.;
    double CountDown = 0.;
    for(int p=0;p<pNPart();p++){
      if(pType(p) != 1) continue;
      if(VAR_IF_TYPE(Ch[pChain(p)].Type,CHAIN_UP)){
	PosUp += pPos(p,CNorm);
	CountUp += 1.;
      }
      else {
	PosDown += pPos(p,CNorm);
	CountDown += 1.;
      }
    }
    PosUp /= CountUp > 1. ? CountUp : 1.;
    PosDown /= CountDown > 1. ? CountDown : 1.;
    ElCoup[0] += (PosUp-PosDown);
    ElCoup[1] += SQR(PosUp-PosDown);
    fprintf(FileToWrite,"%d %lf\n",f,PosUp-PosDown);
  }
  printf("\n");
  ElCoup[1] = sqrt((ElCoup[1] - SQR(ElCoup[0])*(NFile[1]-NFile[1]))/((double)(NFile[0]-NFile[1])-1.));
  fprintf(FileToWrite,"# %lf %lf %lf\n",pEdge(CLat1)*pEdge(CLat2),ElCoup[0],ElCoup[1]);
  fclose(FileToWrite);
}
void ElPoly::BilayerDistance(char *OutFile,int NSample){
  FILE *File2Write = fopen(OutFile,"w");
  int NLayer = 2;
  double PosBf[3];
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_CHAIN))return;
    double *PosLay   = (double *)calloc(4*NSample*NSample,sizeof(double));
    double *CountLay = (double *)calloc(4*NSample*NSample,sizeof(double));
    for(int b=0,cOff=0,pOff=0;b<pNBlock();cOff+=pNChain(b++)){
      int Level = 0;
      if(!strcmp(Block[b].Name,"LIPID1")) Level = 2;
      for(int c=cOff;c<pNChain(b)+cOff;c++,pOff+=pNPCh(b)){
	for(int d=0;d<3;d++){
	  PosBf[d] = Ch[c].Pos[d] - floor(Ch[c].Pos[d]*pInvEdge(d))*pEdge(d);
	}
	int cLevel = Level;
	if(CHAIN_IF_TYPE(Ch[c].Type,CHAIN_UP)) cLevel += 1;
	int sx = (int)(PosBf[CLat1]*pInvEdge(CLat1)*NSample);
	int sy = (int)(PosBf[CLat2]*pInvEdge(CLat2)*NSample);
	if(sx < 0 || sx >= NSample) continue;
	if(sy < 0 || sy >= NSample) continue;
	PosLay[(sx*NSample+sy)*4+cLevel] += Ch[c].Pos[CNorm];//PosBf[CNorm];
	CountLay[(sx*NSample+sy)*4+cLevel] += 1.;
      }
    }
    for(int sx = 0;sx < NSample;sx++){
      for(int sy = 0;sy < NSample;sy++){
	for(int l=0;l<4;l++){
	  PosLay[(sx*NSample+sy)*4+l] /= CountLay[(sx*NSample+sy)*4+l] > 0. ? CountLay[(sx*NSample+sy)*4+l] : 1.;
	}
	double h1 = PosLay[(sx*NSample+sy)*4+2]-PosLay[(sx*NSample+sy)*4+1];
	double h2 = PosLay[(sx*NSample+sy)*4+3]-PosLay[(sx*NSample+sy)*4+0];
	h2 -= floor(h2*pInvEdge(CNorm))*pEdge(CNorm);
	//if(PosLay[(sx*NSample+sy)*4+0]*PosLay[(sx*NSample+sy)*4+1]*PosLay[(sx*NSample+sy)*4+2]*PosLay[(sx*NSample+sy)*4+3] <= 0.) continue;
	fprintf(File2Write,"%lf %lf %lf %lf %lf\n",MIN(h1,h2),PosLay[(sx*NSample+sy)*4+0],PosLay[(sx*NSample+sy)*4+1],PosLay[(sx*NSample+sy)*4+2],PosLay[(sx*NSample+sy)*4+3]);
      }
    }
    free(PosLay);
    free(CountLay);
  }
  printf("\n");
  fclose(File2Write);
}
void ElPoly::BondDistr(char *OutFile,int NBin){
  FILE *File2Write = fopen(OutFile,"w");
  // double **Histo = (double **)calloc(3,sizeof(double));
  // double **Raw = (double **)calloc(3,sizeof(double));
  // for(int d=0;d<3;d++){
  //   Histo[d] = (double *)calloc(NBin,sizeof(double));
  //   Raw[d] = (double *)calloc(pNPart,sizeof(double));
  // }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO)) return;
    double Dist[3];
    for(int c=0;c<pNChain();c++){
      for(int p=0;p<pNPCh()-1;p++){
	int p1 = p + c*pNPCh();
	for(int d=0;d<3;d++){
	  Dist[d] = remainder(pPos(p1,d) - pPos(p1+1,d),pEdge(d));
	}
	fprintf(File2Write,"%lf %lf %lf %lf\n",Dist[0],Dist[1],Dist[2],sqrt(SQR(Dist[0])+SQR(Dist[1])+SQR(Dist[2])));
      }
    }
  }
  // MOMENTI Mom[3];
  // int xMin = Raw[0][0];
  // int xMax = Raw[0][0];
  // for(int d=0;d<3;d++){
  //   Mom[d] = Mat->Distribuzione(Raw[d],pNPart-pNChain,Histo[d],NBin);
  //   Mat->Normalize(Histo[d],NBin);
  //   if(xMin > Mom[d].Minimo) xMin = Mom[d].Minimo;
  //   if(xMax > Mom[d].Massimo) xMax = Mom[d].Massimo;
  // }
  // for(int v=0;v<NBin;v++){
  //   double x = v/(double)(NBin)*(xMax-xMin)+xMin;
  //   fprintf(File2Write,"%lf %lf %lf %lf\n",x,Histo[0][v],Histo[1][v],Histo[2][v]);
  // }
  // for(int d=0;d<3;d++){
  //   free(Histo[d]);
  //   free(Raw[d]);
  // }
  // free(Histo);
  // free(Raw);
  fclose(File2Write);
}
void ElPoly::EndToEndDist(){
  double *EndToEnd = (double *) calloc(pNPart()*(NFile[1]-NFile[0]),sizeof(double));
  double *Distr = (double *) calloc((NFile[1]-NFile[0]),sizeof(double));
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO))return;
    

  }

}
