#include "ElPoly.h"

int ElPoly::Temperature(int NBin,int Coord){
  double *dDensity;
  dDensity = (double *)malloc(NBin*sizeof(double));
  for(int v=0;v<NBin;v++){
    dDensity[v] = 0.;
  }
  FILE *FileToWrite = fopen("Temperature.dat","w");
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NANO)==0)return 0;
    //TemperatureProfile(NBin,Coord);
    for(int v=0;v<NBin;v++){
      //dDensity[v] += dPoint[v];
    }
  }
  printf("\n"); 
  for(int v=0;v<NBin;v++){
    fprintf(FileToWrite,"%lf \n",dDensity[v]/(double)(NFile[1] - NFile[0]));
  }
  //FileToWrite = fopen(OutFile,"w");
  return 0;
}
int ElPoly::NanoParticle(int NBin){
  //SetEdge(.5*sqrt(SQR(pEdge(CLat1))+SQR(pEdge(CLat2))),3);
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  int NType = 2;
  double *dDensity = (double *)calloc(NType*NBin,sizeof(double));
  double *dRad = (double *)calloc(NType*NBin,sizeof(double));
  double *ThickProf = (double *)calloc(NType*NBin,sizeof(*ThickProf));
  double *CountCh = (double *)calloc(NType*NBin,sizeof(double));
  double *CountPart = (double *)calloc(NBin,sizeof(double));
  double *AngleProf = (double *)calloc(2*NBin,sizeof(double));
  double *VolEl = (double *)calloc(NBin,sizeof(double));
  double *VelDistr = (double *)calloc(NBin*2,sizeof(double));
  double *DisplRad = (double *)calloc(NBin*2,sizeof(double));
  double *DisplTheta = (double *)calloc(NBin*2,sizeof(double));
  double *E2E = (double *)calloc(NBin*2,sizeof(double));
  double *HFluct = (double *)calloc(NBin*2,sizeof(double));
  CHAIN *Ch2 = (CHAIN *)calloc(pNChain(),sizeof(CHAIN));
  Vettore Directive(3);
  Vettore Axis(3);
  Vettore ZedDir(3);
  for(int d=0;d<3;d++) ZedDir.Set(0.,d);
  ZedDir.Set(1.,CNorm);
  VolumeCircSlab(VolEl,NBin);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_CHAIN)) return 0;
    //ChainDiffusion(DIFF_LATERAL,ChDiff);
    //Sum radial profile------------------------------------
    for(int b=0;b<pNBlock();b++){
      int IfAdded = 0;
      if(!strncmp(Block[b].Name,"PEP",3) ) continue;
      for(int p=Block[b].InitIdx,link=0;p<Block[b].EndIdx;p++){
	if(pPos(p,CNorm) > pNanoPos(0,CNorm) + Nano->Rad*.5) continue;
	if(pPos(p,CNorm) < pNanoPos(0,CNorm) - Nano->Rad*.5) continue;
	double RadDist = sqrt( SQR(pPos(p,CLat1)-pNanoPos(0,CLat1)) + SQR(pPos(p,CLat2)-pNanoPos(0,CLat2)) );
	double Zed = pPos(p,CNorm);
	int vr = (int)(RadDist/pEdge(3)*NBin);
	if(vr < 0 || vr >= NBin) continue;
	int t = 0;
	if(CHAIN_IF_TYPE(Ch[pType(p)].Type,CHAIN_ADDED) ) t = 1;
	dRad[vr*NType+t] += 1.;
	CountPart[vr] += 1.;
      }
    }
    // DensityProfile(3,NBin,NType,dDensity);
    // for(int v=0;v<NBin;v++){
    //   for(int t=0;t<NType;t++){
    // 	dRad[v*NType+t] += dDensity[v*NType+t];
    // 	dDensity[v*NType+t] = 0.;
    //   }
    // }
    //Sum thickness profile--------------------------------
    for(int b=0,cOff=0;b<pNBlock();cOff+=pNChain(b++)){
      if(!strncmp(Block[b].Name,"PEP",3)){
	continue;
      }
      for(int c=cOff;c<cOff+pNChain(b);c++){
	double RelPos[3];
	for(int d=0;d<3;d++){
	  RelPos[d] = Ch[c].Pos[d] - pNanoPos(0,d);
	  RelPos[d] -= floor(RelPos[d]*pInvEdge(d) + .5)*pEdge(d);
	  Directive.Set(RelPos[d],d);
	  Axis.Set(Ch[c].Dir[d],d);
	}
	double RadDist = sqrt( SQR(RelPos[CLat1]) + SQR(RelPos[CLat2]) );
	double Zed = (Ch[c].Pos[CNorm] - pCm(CNorm));
	int vr = (int)(RadDist*pInvEdge(3)*NBin);
	if(vr < 0 || vr >= NBin) continue;
	int t = 0;
	if(VAR_IF_TYPE(Ch[c].Type,CHAIN_DOWN)) t = 1;
	CountCh[(vr)*2+t] += 1;
	ThickProf[(vr)*2+t] += Zed;
	HFluct[vr*2  ] += Zed;
	HFluct[vr*2+1] += SQR(Zed);
	double Angle = ZedDir.Angle(&Axis,&Directive);
	Angle = Angle < .5*M_PI ? Angle : Angle - .5*M_PI;
	if(!isnan(Angle)){
	  AngleProf[vr*2  ] += Angle;
	  AngleProf[vr*2+1] += SQR(Angle);
	}
	double VelRad = sqrt( SQR(Ch[c].Vel[CLat1]) + SQR(Ch[c].Vel[CLat2]) );
	double VelTheta = fabs(atan2(Ch[c].Vel[CLat2],Ch[c].Vel[CLat1]));
	double DRad = sqrt(SQR(Ch[c].Pos[CLat1]-Ch2[c].Pos[CLat1]) + SQR(Ch[c].Pos[CLat2]-Ch2[c].Pos[CLat2]));
	double DTheta = atan2(Ch[c].Pos[CLat1]-Ch2[c].Pos[CLat1],Ch[c].Pos[CLat2]-Ch2[c].Pos[CLat2]);
	VelDistr[vr*2  ] += VelRad;
	VelDistr[vr*2+1] += VelTheta < 100. ? VelTheta : 0.;
	DisplRad[vr*2  ] += DRad;
	DisplRad[vr*2+1] += SQR(DRad);
	DisplTheta[vr*2  ] += DTheta;
	DisplTheta[vr*2+1] += SQR(DTheta);
	double DistE2E = 0.;
	int pE1 = c*pNPCh();
	int pE2 = (c+1)*pNPCh()-1;
	for(int d=0;d<3;d++){
	  Ch2[c].Pos[d] = Ch[c].Pos[d];
	  DistE2E += SQR(pPos(pE1,d)-pPos(pE2,d));
	}
	E2E[vr*2  ] += sqrt(DistE2E);
	E2E[vr*2+1] += DistE2E;
      }
    }
  }
#ifdef OMPI_MPI_H
  MPI_Allreduce(MPI_IN_PLACE,dDensity,NType*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,ThickProf,NType*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,CountCh,NType*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,CountPart,NType*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,VelDistr,2*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,DisplRad,2*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,DisplTheta,2*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  MPI_Allreduce(MPI_IN_PLACE,E2E,2*NBin, MPI_DOUBLE, MPI_SUM, Proc->CommGrid);
  int Rank=0;
  MPI_Comm_rank(Proc->CommGrid, &Rank);
  if(Rank==0){
#endif
    printf("\n");
    //Normalize-----------------------------------
    double InvNFile = 1./(double)(NFile[1]-NFile[0]);
    double VolumeZed = pEdge(0)*pEdge(1)*pEdge(2)/(double)SQR(NBin);
    double VolumeRad = pEdge(3)*pEdge(3)*Nano->Rad/(double)SQR(NBin);
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      double NormType = 0.;
      double NormVol = (VolumeRad*VolEl[v])/InvNFile;
      NormVol = NormVol > 0. ? 1./NormVol : 1.;
      for(int t=0;t<NType;t++){
	dRad[v*NType+t] *= NormVol;
	double Norm = CountCh[v*NType+t] > 0. ? 1./CountCh[v*NType+t] : 1.;
	ThickProf[v*NType+t] *= Norm;
	NormType += CountCh[v*NType+t];
      }
      NormType = NormType > 0. ? 1./NormType : 1.;
      if(r < Nano->Rad) NormType = 0.;
      if(r < Nano->Rad) NormVol = 0.;
      VelDistr[v*2  ]   *= NormType;
      VelDistr[v*2+1]   *= NormType;
      DisplRad[v*2  ]   *= NormType;
      DisplRad[v*2+1]    = sqrt(DisplRad[v*2+1]*NormType - SQR(DisplRad[v*2  ]));
      DisplTheta[v*2  ] *= NormType;
      DisplTheta[v*2+1]  = sqrt(DisplTheta[v*2+1]*NormType - SQR(DisplTheta[v*2  ]));
      HFluct[v*2  ]     *= NormType;
      HFluct[v*2+1]      = sqrt(HFluct[v*2+1]*NormType - SQR(HFluct[v*2  ]));
      AngleProf[v*2  ]  *= NormType;
      AngleProf[v*2+1]   = sqrt(AngleProf[v*2+1]*NormType - SQR(AngleProf[v*2  ]));
      E2E[v*2  ]        *= NormType;
      E2E[v*2+1]         = sqrt( (E2E[v*2+1]*NormType - SQR(E2E[v*2  ])) );
      //E2E[v*2+1]        *= NormVol;
    }
    //--------------radial density-------------------------------
    char FileName[256];
    sprintf(FileName,"NanoDensR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FILE *FileToWrite = fopen(FileName,"w");
    char cSystem[STRSIZE];
    SysDef(cSystem);
    fprintf(FileToWrite,"# %s",cSystem);
    fprintf(FileToWrite,"#r DensPhob DensPhil\n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %lf %lf\n",r,dRad[v*NType],dRad[v*NType+1]);
    }
    fclose(FileToWrite);
    // //---------------thickness profile-----------------------------
    // sprintf(FileName,"NanoThickR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    // FileToWrite = fopen(FileName,"w");
    // fprintf(FileToWrite,"#r ThickUp ThickDown\n");
    // for(int vr=0;vr<NBin;vr++){
    //   double r = vr*pEdge(3)/(double)NBin;
    //   if(r < Nano->Rad) continue;
    //   if(ThickProf[vr*NType+0]*ThickProf[vr*NType+1] > 0.)
    //     fprintf(FileToWrite,"%lf %lf %lf\n",r,ThickProf[vr*NType+0],ThickProf[vr*NType+1]);
    // }
    // fclose(FileToWrite);
    //----------------velocity profile-----------------------------
    sprintf(FileName,"NanoRadVelR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FileToWrite = fopen(FileName,"w");
    fprintf(FileToWrite,"#r VelRad VelTheta\n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %lf %lf\n",r,VelDistr[v*2+0],VelDistr[v*2+1]);
    }
    fclose(FileToWrite);
    //----------------z fluctuation profile-----------------------------
    sprintf(FileName,"NanoZedFluctR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FileToWrite = fopen(FileName,"w");
    fprintf(FileToWrite,"#r <z> <(z - <z>)^2> \n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %.5g %.5g\n",r,HFluct[v*2+0],HFluct[v*2+1]);
    }
    fclose(FileToWrite);
    //----------------angle profile-----------------------------
    sprintf(FileName,"NanoAngleR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FileToWrite = fopen(FileName,"w");
    fprintf(FileToWrite,"#r <a> <(a - <a>)^2> \n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %.5g %.5g\n",r,AngleProf[v*2+0],AngleProf[v*2+1]);
    }
    fclose(FileToWrite);
    //----------------end to end-----------------------------
    sprintf(FileName,"NanoE2ER%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FileToWrite = fopen(FileName,"w");
    fprintf(FileToWrite,"#r <e2e> <(e2e - <e2e>)^2> \n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %.5g %.5g\n",r,E2E[v*2+0],E2E[v*2+1]);
    }
    fclose(FileToWrite);
    //----------------radial displacement-----------------------------
    sprintf(FileName,"NanoRadDisplR%.1fS%.1fH%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
    FileToWrite = fopen(FileName,"w");
    fprintf(FileToWrite,"#r JumpRad JumpTheta \n");
    for(int v=0;v<NBin;v++){
      double r = v*pEdge(3)/(double)NBin;
      fprintf(FileToWrite,"%lf %lf %lf %lf %lf\n",r,DisplRad[v*2  ],DisplRad[v*2+1],DisplTheta[v*2  ],DisplTheta[v*2+1]);
    }
    fclose(FileToWrite);
#ifdef OMPI_MPI_H
  }
#endif
  //freeing-------------------------------------------------
  free(CountCh);
  free(CountPart);
  free(DisplRad);
  free(DisplTheta);
  free(dDensity);
  free(dRad);
  free(ThickProf);
  free(VelDistr);
  free(AngleProf);
  return 0;
}
int ElPoly::WormF(int Partition,int NBin){
  double *dDensity = (double *)calloc(2*NBin,sizeof(double));
  double Norma;
  double Border[2] = {-10.,10.};
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return 0;
    Worm(Partition,NBin,Border,dDensity);
  }
  printf("\n");
  FILE *FileToWrite = fopen("WormDensity.dat","w");
  for(int v=0;v<NBin;v++){
    fprintf(FileToWrite,"%lf %lf %lf\n",(double)v/(double)NBin*(Border[1]-Border[0])+Border[0],dDensity[v*2]/(NFile[1]-NFile[0]),dDensity[v*2+1]/(NFile[1]-NFile[0]));
  }
  fclose(FileToWrite);
  free(dDensity);
  return 1;
}
//The skin is not closing, better use the marching cubes
void ElPoly::StalkF(int NSample){
  double Threshold = 0.;
  int NLevel = 4;
  double **Plot = (double **) calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++){
    Plot[l] = (double *)calloc(NSample*NSample,sizeof(double));
  }
  double **PlotA = (double **) calloc(NLevel,sizeof(double));
  for(int l=0;l<NLevel;l++){
    PlotA[l] = (double *)calloc(NSample*NSample,sizeof(double));
  }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return ;
    Stalk(NSample,NLevel,Plot,Threshold);
    for(int s=0;s<NSample;s++)
      for(int ss=0;ss<NSample;ss++)
	for(int l=0;l<NLevel;l++)
	  PlotA[l][s*NSample+ss] += Plot[l][s*NSample+ss];
  }
  printf("\n");
  Vettore v00(3);
  Vettore v01(3);
  Vettore v11(3);
  Vettore v10(3);
  Vettore vN(3);
  char FileName[256];
  sprintf(FileName,"Stalk.xvl");
  FILE *FileToWrite = fopen(FileName,"w");
  fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NSample,ChooseDraw(EL_SKIN));
  double FactC1 = pEdge(CLat1)/(double)NSample;
  double FactC2 = pEdge(CLat2)/(double)NSample;
  double FactCN = 1./(double)(NFile[1]-NFile[0]);
  for(int l=0,p=0,c=0;l<NLevel;l++){
    for(int ss=0;ss<NSample-1;ss++){
      for(int s=0;s<NSample-1;s++){
	v00.x[CLat1]=(s+.5)*FactC1;
	v00.x[CLat2]=(ss+.5)*FactC2;
	v00.x[CNorm]=PlotA[l][s*NSample+ss]*FactCN;
	v01.x[CLat1]=(s+1.5)*FactC1;
	v01.x[CLat2]=(ss+.5)*FactC2;
	v01.x[CNorm]=PlotA[l][(s+1)*NSample+ss]*FactCN;
	v11.x[CLat1]=(s+1.5)*FactC1;
	v11.x[CLat2]=(ss+1.5)*FactC2;
	v11.x[CNorm]=PlotA[l][(s+1)*NSample+ss+1]*FactCN;
	v10.x[CLat1]=(s+.5)*FactC1;
	v10.x[CLat2]=(ss+1.5)*FactC2;
	v10.x[CNorm]=PlotA[l][s*NSample+ss+1]*FactCN;
	//------------defines-the-squares---------------------
	if(PlotA[l][(s)*NSample+ss] <= Threshold){
	  if(l==0)
	    v00.x[CNorm]=PlotA[3][(s)*NSample+ss+1]*FactCN;
	  if(l==1)
	    v00.x[CNorm]=PlotA[2][(s)*NSample+ss+1]*FactCN;
	  if(l==2) 
	    v00.x[CNorm]=PlotA[1][(s)*NSample+ss+1]*FactCN;
	  if(l==3) 
	    v00.x[CNorm]=PlotA[0][(s)*NSample+ss+1]*FactCN;
	  if(v00.x[CNorm] <= 0.) continue;
	}
	if(PlotA[l][(s+1)*NSample+ss] <= Threshold){
	  if(l==0)
	    v01.x[CNorm]=PlotA[3][(s)*NSample+ss]*FactCN;
	  if(l==1)
	    v01.x[CNorm]=PlotA[2][(s)*NSample+ss]*FactCN;
	  if(l==2) 
	    v01.x[CNorm]=PlotA[1][(s)*NSample+ss]*FactCN;
	  if(l==3) 
	    v01.x[CNorm]=PlotA[0][(s)*NSample+ss]*FactCN;
	  if(v01.x[CNorm] <= 0.) continue;
	}
	if(PlotA[l][s*NSample+ss+1] <= Threshold){
	  if(l==0)
	    v10.x[CNorm]=PlotA[3][(s)*NSample+ss]*FactCN;
	  if(l==1)
	    v10.x[CNorm]=PlotA[2][(s)*NSample+ss]*FactCN;
	  if(l==2) 
	    v10.x[CNorm]=PlotA[1][(s)*NSample+ss]*FactCN;
	  if(l==3) 
	    v10.x[CNorm]=PlotA[0][(s)*NSample+ss]*FactCN;
	  if(v10.x[CNorm] <= 0.) continue;
	}
	if(PlotA[l][(s+1)*NSample+(ss+1)] <= Threshold){
	  if(l==0)
	    v11.x[CNorm]=PlotA[3][(s)*NSample+ss+1]*FactCN;
	  if(l==1)
	    v11.x[CNorm]=PlotA[2][(s)*NSample+ss+1]*FactCN;
	  if(l==2) 
	    v11.x[CNorm]=PlotA[1][(s)*NSample+ss+1]*FactCN;
	  if(l==3) 
	    v11.x[CNorm]=PlotA[0][(s)*NSample+ss+1]*FactCN;
	  if(v11.x[CNorm] <= 0.) continue;
	}
	fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p+0,c,l,v00.x[CLat1],v00.x[CLat2],v00.x[CNorm],(double)l*.25,.5,.8,p+1,p+2,p+3);
	fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p+1,c,l,v01.x[CLat1],v01.x[CLat2],v01.x[CNorm],(double)l*.25,.5,.8,p+2,p+3,p+0);
	fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p+2,c,l,v11.x[CLat1],v11.x[CLat2],v11.x[CNorm],(double)l*.25,.5,.8,p+3,p+0,p+1);
	fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf) l[%d] l[%d] l[%d] }\n",p+3,c,l,v10.x[CLat1],v10.x[CLat2],v10.x[CNorm],(double)l*.25,.5,.8,p+0,p+1,p+2);
	p += 4;
	c++;
      }
    }
  }
  for(int l=0;l<NLevel;l++){
    free(Plot[l]);
    free(PlotA[l]);
  }
  free(Plot);
  free(PlotA);
}
/**
   Call the function StalkLineProf to reconstruct the linear shape of the stalk and calculate the power spectrum.
 */
void ElPoly::StalkLineProfF(int NBin){
#ifdef USE_FFTW
  double *Line = (double *)calloc(NBin,sizeof(double));
  double *Spe = (double *)calloc(NBin,sizeof(double));
  char FName[120];
  int NHalf = (int)(NBin/2);
  // Excluding the two ends that are at most di
  double InvNBin = 1./(double)(NBin-2);
  double InvNFile = 1./(double)(NFile[1]-NFile[0]);
  fftw_complex *out = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_complex *in = (fftw_complex *)fftw_malloc(NBin*sizeof(fftw_complex));
  fftw_plan plan = fftw_plan_dft_1d(NBin-2,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan plan2 = fftw_plan_dft_1d(NBin-2,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_CHAIN)) return;
    StalkLineProf(Line,NBin);
    for(int vx=0;vx<NBin-2;vx++) in[vx][0] = Line[vx+1];
    fftw_execute(plan);
    for(int vx=0;vx<NBin-2;vx++){
      //int vvx = vx < NHalf ? vx + NHalf : vx - NHalf;
      Spe[vx] += (SQR(out[vx][0]) + SQR(out[vx][1]))*SQR(InvNBin);
    }
    //-------------------------write-temp-files
    if(1==1){
      sprintf(FName,"StalkLine%05d.dat",f);
      FILE *FWrite = fopen(FName,"w");
      for(int vx=0;vx<NBin;vx++){
	fprintf(FWrite,"%lf %lf \n",vx*pEdge(CLat2)/(double)NBin,Line[vx]);
      }
      fclose(FWrite);
    }
  }
  printf("\n");
  double dx = pEdge(CLat1)*InvNBin;
  FILE *FWrite = fopen("StalkLineSpectrum.dat","w");
  for(int vx=0;vx<NBin/2-2;vx++){
    double qx  = 2.*M_PI*vx*pInvEdge(CLat1);
    double qq = SQR(qx/dx - .5/dx);
    //fprintf(FWrite,"%lf %lf \n",qx,Spe[vx]*InvNFile*pEdge(CLat2)*pEdge(CLat2));
    fprintf(FWrite,"%lf %lf \n",qx,Spe[vx]*InvNFile);
  }
  fclose(FWrite);
  free(Line);
  free(Spe);
  fftw_free(out);
  fftw_free(in);
#else
  printf("fftw not implemented\n");
#endif // USE_FFTW
}
void ElPoly::Prova(){
  char FName[60];
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_NO))return;
    Nano[0].Pos[0] = 12.0;
    Nano[0].Pos[1] = 12.0;
    Nano[0].Pos[2] = 9.0;
    sprintf(FName,"new%07d00.dat",f*2);
    Write(cFile[f]);
  }
  printf("\n");
}
void ElPoly::RadNormPos(int NBin,int NGrid){
  //SigErr((NFile[1]-NFile[0]) != SQR(NGrid),"RadNormPos: the number of files is not compatible with the parameter span");
  char FName[120];
  int NWeight = 5;
  double *RadProf = (double *)calloc(NBin,sizeof(double));
  double *RadProfS = (double *)calloc(NBin,sizeof(double));
  double *Count = (double *)calloc(NBin,sizeof(double));
  double *PlotMin = (double *)calloc(NGrid*NGrid,sizeof(double));
  double *PlotThin = (double *)calloc(NGrid*NGrid,sizeof(double));
  double *PlotTmp = (double *)calloc(NGrid*NGrid,sizeof(double));
  double *AngList = (double *)calloc((NFile[1]-NFile[0]),sizeof(double));
  double *HeiList = (double *)calloc((NFile[1]-NFile[0]),sizeof(double));
  double *Weight = (double *)calloc(NWeight,sizeof(double));
  int *WIndex = (int *)calloc(NWeight,sizeof(int));
  double InvNBin = 1./(NBin);
  Mat->FillWeightGauss(Weight,WIndex,NWeight,0.0,2.0);
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_PART)) return;
    double Angle = Nano->Hamaker;
    double Hei = Nano->Height;
    for(int r=0;r<NBin;r++){
      RadProf[r] = 0.;
      Count[r] = 0.;
    }
    for(int p=0;p<pNPart();p++){
      double Rad = sqrt(SQR(Pm[p].Pos[CLat1]-Nano->Pos[CLat1])+SQR(Pm[p].Pos[CLat2]-Nano->Pos[CLat2]));
      //if(Rad < 1.5*Nano->Rad) continue;
      int r = (int)(Rad*pInvEdge(3)*NBin);
      if(r < 0 || r >= NBin) continue;
      RadProf[r] += Pm[p].Pos[CNorm];
      Count[r] += 1.;
    }
    for(int r=0;r<NBin;r++){
      if(Count[r] <= 0.) continue;
      RadProf[r] /= Count[r];
    }
    //smooth
    Matrice Mask1(5);
    Mask1.FillGaussian(.5,3.);
    int NDim = 1;
    int IfMinImConv = 0;
    Mask1.ConvoluteMatrix(RadProf,NBin,NDim,IfMinImConv);
    // for(int v=0;v<NBin;v++) RadProfS[v] = RadProf[v];
    // InterBSpline1D(RadProfS,RadProf,NBin,NBin);
    // for(int v=0;v<NBin;v++) RadProfS[v] = RadProf[v];
    // InterBSpline1D(RadProfS,RadProf,NBin,NBin);
    // for(int v=0;v<NBin;v++) RadProfS[v] = RadProf[v];
    // InterBSpline1D(RadProfS,RadProf,NBin,NBin);
    // Mat->ConvWeight(RadProf,NBin,Weight,WIndex,NWeight);
    double Thinning = 100000.;
    double MinDist = 0.;
    int rLower = 0;
    for(int r=3;r<NBin;r++){
      if(Count[r] <= 0.) continue;
      if(Thinning > RadProf[r]){
	Thinning = RadProf[r];
	MinDist = r*pEdge(3)/(double)NBin;
	rLower = r;
      }
    }
    //printf("%lf %lf %lf %lf\n",Angle,Hei,MinDist,Thinning);
    AngList[f] = Angle;
    HeiList[f] = Hei;
    PlotMin[f] = MinDist;
    PlotThin[f] = Thinning;
    sprintf(FName,"ThinProfHei%.1fAng%02d.dat",Hei,(int)Angle);
    FILE *FProf = fopen(FName,"w");
    for(int r=0;r<NBin;r++){
      if(Count[r] <= 0.) continue;	
      double Rad = r*pEdge(3)*InvNBin;
      fprintf(FProf,"%lf %lf\n",Rad,RadProf[r]);
    }
    fclose(FProf);
  }
  //smooth and write dens
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 2;
  int IfMinImConv = 1;
  for(int t=0;t<0;t++){
    Mask.ConvoluteMatrix(PlotMin,NGrid,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(PlotThin,NGrid,NDim,IfMinImConv);
  }
  for(int s=0;s<NGrid*NGrid;s++) PlotTmp[s] = PlotMin[s];
  InterBSpline2D(PlotTmp,PlotMin,NGrid,NGrid);
  for(int s=0;s<NGrid*NGrid;s++) PlotTmp[s] = PlotThin[s];
  InterBSpline2D(PlotTmp,PlotThin,NGrid,NGrid);
  FILE *FMinDist = fopen("HeiAngMinDist.dat","w");
  FILE *FThin = fopen("HeiAngThin.dat","w");
  for(int f=NFile[0];f<NFile[1];f++){
    double Angle = AngList[f];
    double Hei = HeiList[f];
    fprintf(FMinDist,"%lf %lf %lf\n",Hei,Angle,PlotMin[f]);
    fprintf(FThin,"%lf %lf %lf\n",Hei,Angle,PlotThin[f]);
  }
  free(Count);
  free(RadProf);
  free(RadProfS);
  free(PlotMin);
  free(PlotThin);
  fclose(FThin);
  fclose(FMinDist);
  printf("\n");
}
#ifdef USE_CGAL
void ElPoly::AreaDistrF(int NBin){
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  double *Distr = (double *)calloc(NBin,sizeof(double));
  double *DistrRadAv = (double *)calloc(NBin,sizeof(double));
  double *DistrRad = (double *)calloc(NBin,sizeof(double));
  double *VolContr = (double *)calloc(NBin,sizeof(double));
  VolumeCircSlab(VolContr,NBin);
  double FNormaInv = 1./(double)(NFile[1]-NFile[0]);
  double BoxRadInv = .5*sqrt( SQR(pEdge(CLat1)) + SQR(pEdge(CLat2)) )/(double)NBin;
  double AreaMean = 3.*.5*pNChain()/(pEdge(CLat1)*pEdge(CLat2));
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_CHAIN))return ;
    BfDefChain();
    AreaDistr(Distr,DistrRad,NBin);
    for(int v=0;v<NBin;v++)
      DistrRadAv[v] += DistrRad[v];
  }
  printf("\n");
  //---------normalize------------------
  double NormF = NFile[1] - NFile[0];
  for(int v=0;v<NBin;v++){
    DistrRadAv[v] /= VolContr[v]*NormF;
    Distr[v] /= VolContr[v]*NormF;
  }
  //---------write----------------------
  char FileName[120];
  sprintf(FileName,"NanoAreaRadDistrR%.1fH%.1fh%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
  FILE *File2Write = fopen(FileName,"w");
  for(int v=0;v<NBin;v++){
    fprintf(File2Write,"%lf %lf\n",v*BoxRadInv,DistrRadAv[v]*FNormaInv);
  }
  fclose(File2Write);
  sprintf(FileName,"NanoAreaDistrR%.1fH%.1fh%.1f.dat",Nano->Rad,Nano->Hamaker,Nano->Height);
  File2Write = fopen(FileName,"w");
  for(int v=0;v<NBin;v++){
    fprintf(File2Write,"%lf %lf\n",v*AreaMean/(double)NBin,Distr[v]*FNormaInv);
  }
  fclose(File2Write);
  free(Distr);
  free(DistrRadAv);
  free(DistrRad);
}
#else
void ElPoly::AreaDistrF(int NBin){
  printf("CGAL libraries not implemented\n");
}
#endif
