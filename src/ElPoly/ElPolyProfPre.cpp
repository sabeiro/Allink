#include "ElPoly.h"

void ElPoly::SumTens(){
  int NType = 3;
  double NumDiff = 0.001;
  double **Plot = (double **)calloc(NType,sizeof(double));
  int v[3];
  double FNorma = 1./(double)(NFile[1]-NFile[0]);
  char FileName[256];
  for(int t=0;t<3;t++) Plot[t] = (double *)calloc(CUB(NEdge),sizeof(double));
  double LatLim[3] = {pEdge(0),pEdge(1),pEdge(2)};
  double InvLatLim[3] = {1./pEdge(0),1./pEdge(1),1./pEdge(2)};
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_NO))return;
    //ShiftSys(SHIFT_CM);
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
    	v[d] = (int)((pPos(p,d)+NumDiff)*InvLatLim[d]*NEdge);
	if(v[d] < 0 || v[d] >=NEdge) continue;
      }
      int vTot = v[0] + NEdge*(v[1] + NEdge*v[2]);
      for(int t=0;t<NType;t++)
    	Plot[t][vTot] += pVel(p,t);
    }
  }
  sprintf(FileName,"Av%s",cFile[NFile[0]]);
  FILE *FileToWrite = fopen(FileName,"w");
  if(FileToWrite == NULL){printf("Can't open %s\n",FileName);return;}
  fprintf(FileToWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",LatLim[0],LatLim[1],LatLim[2],Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],NEdge);
  double NEdgeInv = 1./(double)NEdge;
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<NType;t++){
    Mask.ConvoluteMatrix(Plot[t],NEdge,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot[t],NEdge,NDim,IfMinImConv);
  }
  for(int vx=0;vx<NEdge;vx++)
    for(int vy=0;vy<NEdge;vy++)
      for(int vz=0;vz<NEdge;vz++){
	int vTot = vx + NEdge*(vy + NEdge*vz);
	double x = vx*LatLim[0]*NEdgeInv;
	double y = vy*LatLim[1]*NEdgeInv;
	double z = vz*LatLim[2]*NEdgeInv;
	if(fabs(Plot[0][vTot] + Plot[1][vTot] + Plot[2][vTot]) > 0.){
	  fprintf(FileToWrite,"{x(%.3f %.3f %.3f) v(%lf %.2f %.2f)}\n",
		  x,y,z,
		  Plot[0][vTot]*FNorma,Plot[1][vTot]*FNorma,Plot[2][vTot]*FNorma);
	}
      }
  for(int t=0;t<NType;t++)
    free(Plot[t]);
  free(Plot);
  fclose(FileToWrite);
}
int ElPoly::PressRadial(){
  int NZed = NEdge;
  int NRad = NEdge;
  int NType = 3;
  double *Plot = (double *)calloc(NZed*NRad*NType,sizeof(double));
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  //Vettore NanoAxis(0,0,1);
  //Vettore PosRel(3);
  //Vettore Dist(3);
  for(int p=0;p<pNPart();p++){
    double Rad = 0.;
    for(int d=0;d<2;d++){
      Rad += SQR(remainder(Nano->Pos[d] - pPos(p,d) - 0.001,pEdge(d)));
      //PosRel.Set(remainder(Nano->PosBf[d] - pPos(p,d),pEdge(d)),d);
    }
    //double RadDist = ABS(Dist.PerpTo3(&PosRel,&NanoAxis));
    double RadDist = sqrt(Rad);
    int vr = (int)(RadDist/pEdge(3)*NRad);
    int vz = (int)((pPos(p,CNorm)+.01)/pEdge(CNorm)*NZed);
    if(vr < 0 || vr >= NRad) continue;
    if(vz < 0 || vz >= NZed) continue;
    for(int t=0;t<3;t++){
      Plot[(vr*NZed+vz)*NType+t] += pVel(p,t);
    }
  }
  //------------Normalize----------------------
  double *VolContr = (double *)calloc(NRad,sizeof(double));
  VolumeCircSlab(VolContr,NRad);
  double Bound[3];
  double InvNZed = 1./(double)NZed;
  for(int t=0,n=0;t<NType;t++){
    Bound[t] = 0.;
    for(int vz=0;vz<NZed;vz++){
      for(int vr=0;vr<NRad;vr++){
	Plot[(vr*NZed+vz)*NType+t] /= VolContr[vr]*3.;
	if(Bound[t] < Plot[(vr*NZed+vz)*NType+t])
	  Bound[t] = Plot[(vr*NZed+vz)*NType+t];
      }
    }
    if(Bound[t] < 0.) Bound[t] = 1.;
  }
  free(VolContr);  
  //-------------------Write---------------------------
  FILE *PRadial = fopen("ContourPress.xvl","w");
  fprintf(PRadial,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(3),pEdge(CNorm),1.,NZed,ChooseDraw(EL_QUAD1));
  FILE *PNormal = fopen("PressNormal.dat","w");
  for(int vr=0;vr<NRad;vr++){
    double PRad[3] = {0.,0.,0.};
    double r = vr*pEdge(3)/(double)NRad;
    for(int vz=0;vz<NZed;vz++){
      //if(Plot[(v*NBin+vv)*NType+0] + Plot[(v*NBin+vv)*NType+1] + Plot[(v*NBin+vv)*NType+2] + Plot[(v*NBin+vv)*NType+3] < .1) continue;
      double z = vz*pEdge(CNorm)/(double)(NZed);
      double Press = Plot[(vr*NZed+vz)*NType+0];
      double Phob  = Plot[(vr*NZed+vz)*NType+1];
      double Phil  = Plot[(vr*NZed+vz)*NType+2];
      if(ABS( Press + Phob + Phil) > 0.)
      fprintf(PRadial,"{x(%lf %lf 0.) v(%lf %lf %lf)}\n",r,z,Press,Phob,Phil);
      PRad[0] += Press;PRad[1] += Phob;PRad[2] += Phil;
    }
    fprintf(PNormal,"%lf %lf %lf %lf\n",r,PRad[0]*InvNZed,PRad[1]*InvNZed,PRad[2]*InvNZed);
  }
  fclose(PRadial);
  fclose(PNormal);
  free(Plot);
  free(VolContr);
}
int ElPoly::Tens2dCartRad(){
  int NType = 3;
  // Press, phob, phil
  double *TensRad = (double *)calloc(NEdge*NEdge*NType,sizeof(double));
  double *TensNorm = (double *)calloc(NEdge*NEdge*NType,sizeof(double));
  double *TensAng = (double *)calloc(NEdge*NEdge*NType,sizeof(double));
  //x y z 
  double *TensRadTemp = (double *)calloc(NEdge*NEdge*3,sizeof(double));
  if(NFile[1]-NFile[0] != 3) return 1;
  for(int f=NFile[0];f<NFile[1];f++){
    if(Open(cFile[f],BF_PART))return 1;
    for(int p=0;p<pNPart();p++){
      int vr = (int)((pPos(p,0)+.0001)*NEdge*pInvEdge(0));
      int vz = (int)((pPos(p,1)+.0001)*NEdge*pInvEdge(1));
      int vTot = (vr*NEdge+vz)*NType;
      //pressure
      TensRadTemp[vTot+f] += pVel(p,0);
      //phob density
      TensRad[vTot+1] += pVel(p,1)/3.;
      //phil density
      TensRad[vTot+2] += pVel(p,2)/3.;
    }
  }
  for(int vr=0;vr<NEdge;vr++){
    for(int vz=0;vz<NEdge;vz++){
      int vTot = (vr*NEdge+vz)*NType;
      double Rad = sqrt(SQR(TensRadTemp[vTot+0])+SQR(TensRadTemp[vTot+1]));
      double Trace = TensRadTemp[vTot+0] + TensRadTemp[vTot+1] + TensRadTemp[vTot+2];
      double Ang = .5*(Trace - Rad);
      //double Ang = atan2(TensRadTemp[vTot+1],TensRadTemp[vTot+0]);
      TensRad[vTot+0] = Rad - .5*(Ang+TensRadTemp[vTot+2]);
      //TensNorm[vTot] = - TensRadTemp[vTot+2] + .5*(Ang+Rad);
      TensNorm[vTot] = TensRadTemp[vTot+2] - .5*(TensRadTemp[vTot+0]+TensRadTemp[vTot+1]);
      TensAng[vTot+0] = Ang - .5*(Rad+TensRadTemp[vTot+2]);
    }
  }
  FILE *RadSurfTens = fopen("SurfTensRad.xvl","w");
  FILE *NormSurfTens = fopen("SurfTensNorm.xvl","w");
  FILE *AngSurfTens = fopen("SurfTensAng.xvl","w");
  fprintf(RadSurfTens,"# l(%.2f %.2f %.2f) r(%.2f %.2f %.2f) v[%d] d[color]\n",pEdge(0),pEdge(1),pEdge(2),Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],NEdge);
  fprintf(NormSurfTens,"# l(%.2f %.2f %.2f) r(%.2f %.2f %.2f) v[%d] d[color]\n",pEdge(0),pEdge(1),pEdge(2),Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],NEdge);
  fprintf(AngSurfTens,"# l(%.2f %.2f %.2f) r(%.2f %.2f %.2f) v[%d] d[color]\n",pEdge(0),pEdge(1),pEdge(2),Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],NEdge);
  for(int vr=0;vr<NEdge;vr++){
    for(int vz=0;vz<NEdge;vz++){
      int vTot = (vr*NEdge+vz)*NType;
      double r = vr*pEdge(0)/NEdge;
      double z = vz*pEdge(1)/NEdge;
      if(fabs(TensRad[vTot+0]+TensRad[vTot+1]+TensRad[vTot+2])>0.){
	fprintf(RadSurfTens,"{x(%.4f %.4f %.4f) v(%lf %lf %lf)}\n",
	      r,z,0.,TensRad[vTot+0],TensRad[vTot+1],TensRad[vTot+2]);
	fprintf(NormSurfTens,"{x(%.4f %.4f %.4f) v(%lf %lf %lf)}\n",
	      r,z,0.,TensNorm[vTot],TensRad[vTot+1],TensRad[vTot+2]);
	fprintf(AngSurfTens,"{x(%.4f %.4f %.4f) v(%lf %lf %lf)}\n",
	      r,z,0.,TensAng[vTot],TensRad[vTot+1],TensRad[vTot+2]);
      }
    }
  }
  free(TensRad);
  free(TensNorm);
  free(TensAng);
  free(TensRadTemp);
  fclose(RadSurfTens);
  fclose(NormSurfTens);
  fclose(AngSurfTens);
}
int ElPoly::SurfTens(int NBin){
  int NType = 3;
  double **TensRad = (double **)calloc(NBin,sizeof(double));
  double **NormRad = (double **)calloc(NBin,sizeof(double));
  double **TensAng = (double **)calloc(NBin,sizeof(double));
  double **NormAng = (double **)calloc(NBin,sizeof(double));
  double **TensNorm= (double **)calloc(NEdge,sizeof(double));
  double **NormNorm= (double **)calloc(NEdge,sizeof(double));
  double **TensCart = (double **)calloc(NBin,sizeof(double));
  double **NormCart = (double **)calloc(NBin,sizeof(double));
  for(int v=0;v<NBin;v++){
    TensRad[v] = (double *)calloc(NType,sizeof(double));
    NormRad[v] = (double *)calloc(NType,sizeof(double));
    TensNorm[v] = (double *)calloc(NType,sizeof(double));
    NormNorm[v] = (double *)calloc(NType,sizeof(double));
    TensAng[v] = (double *)calloc(NType,sizeof(double));
    NormAng[v] = (double *)calloc(NType,sizeof(double));
    TensCart[v] = (double *)calloc(NType,sizeof(double));
    NormCart[v] = (double *)calloc(NType,sizeof(double));
  }
  //double *Plot = (double *)calloc(NBin*NBin,sizeof(double));
  SetEdge(MIN((.5*pEdge(CLat1)),(.5*pEdge(CLat2))),3);
  double *VolContr = (double *)calloc(NBin,sizeof(double));
  VolumeCircSlab(VolContr,NBin);
  double Prex = 0.;
  double Prey = 0.;
  int IfTens = 0;
  if(NFile[1]-NFile[0] == 3) IfTens = 1;
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_PART))return 0;
    //BackFold(BF_PART);
    //Vettore NanoAxis(Nano->Axis[0],Nano->Axis[1],Nano->Axis[2]);
    Vettore NanoAxis(0,0,1);
    Vettore PosRel(3);
    Vettore Dist(3);
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	PosRel.Set(remainder(pNanoPos(0,d) - pPos(p,d),pEdge(d)),d);
      }
      double RadDist = ABS(Dist.PerpTo3(&PosRel,&NanoAxis));
      int vr = (int)(RadDist/pEdge(3)*NBin);
      int vz = (int)((pPos(p,CNorm)+.01)/pEdge(CNorm)*NEdge);
      if(vr < 0 || vr >= NBin){continue;}
      if(vz < 0 || vz >= NEdge){printf("%d out of %d\n",vz,NBin);continue;}
      for(int t=1;t<3;t++){
	TensRad[vr][t] += pVel(p,t);
	NormRad[vr][t] += 1.;
	TensNorm[vz][t] += pVel(p,t);
	NormNorm[vz][t] += 1.;
      }
      if(IfTens == 0){
	TensRad[vr][0] += pVel(p,0);
	NormRad[vr][0] += 1.;
	TensNorm[vz][0] += pVel(p,0);
	NormNorm[vz][0] += 1.;
      }
      else if(IfTens){
	TensCart[vr][f] += pVel(p,0);
	NormCart[vr][f] += 1.;
	double vPre = (pVel(p,1)+pVel(p,2)-pVel(p,0));
      	if(f==CLat1 || f == CLat2){
      	  TensNorm[vz][0] -= .5*vPre;
	  if(NType>3){
	    TensNorm[vz][3] -= .5*vPre*pPos(p,CNorm);
	    NormNorm[vz][3] += 1.;
	  }
	  if(NType>4){
	    TensNorm[vz][4] -= .5*vPre*pPos(p,CNorm)*pPos(p,CNorm);
	    NormNorm[vz][4] += 1.;
	  }
	  if(f==CLat1)
	    Prex += vPre;
	  else 
	    Prey += vPre;
      	}
     	if(f==CNorm){
      	  TensNorm[vz][0] += vPre;
	  if(NType>3){
	    TensNorm[vz][3] += vPre*pPos(p,CNorm);
	    NormNorm[vz][3] += 1.;
	  }
	  if(NType>4){
	    TensNorm[vz][4] += vPre*pPos(p,CNorm)*pPos(p,CNorm);
	    NormNorm[vz][4] += 1.;
	  }
      	}
      	NormNorm[vz][0] += 1.;
      }
    } 
  }
  printf("\n");
  double VolElm = 1.;//pEdge(CNorm)/(double)NBin;
  FILE *File2Write = fopen("TensRadial.dat","w");
  for(int v=0;v<NBin;v++){
    for(int d=0;d<3;d++){
      TensCart[v][d] /= NormCart[v][d] > 0. ? VolContr[v]*NormCart[v][d] : 1.;
    }
    TensRad[v][0] = sqrt( SQR(TensCart[v][CLat1])+SQR(TensCart[v][CLat2]) );
    TensRad[v][0] -= 5.*.5*(TensCart[v][0]+TensCart[v][1]+TensCart[v][2]-sqrt( SQR(TensCart[v][CLat1])+SQR(TensCart[v][CLat2])) );
    //TensRad[v][0] += atan2(TensCart[v][CLat2],TensCart[v][CLat1]);
    TensRad[v][0] -= .5*TensCart[v][CNorm];
  }
  for(int v=0;v<NBin;v++){
    fprintf(File2Write,"%lf ",v/(double)NBin*pEdge(3));
    for(int t=0;t<NType;t++){
      TensRad[v][t] /= NormRad[v][t] > 0. ? NormRad[v][t]: 1.;
      fprintf(File2Write,"%lf ",TensRad[v][t]);
    }
    fprintf(File2Write,"\n");
  }
  fclose(File2Write);
  File2Write = fopen("TensNormal.dat","w");
  fprintf(File2Write,"#Press DensPhob DensPhil SponCurv SaddleSplay \n");
  for(int v=0;v<NEdge;v++){
    fprintf(File2Write,"%lf ",v/(double)NEdge*pEdge(CNorm));
    for(int t=0;t<NType;t++){
      TensNorm[v][t] /= NormNorm[v][t] > 0. ? VolElm*NormNorm[v][t] : 1.;
      fprintf(File2Write,"%lf ",TensNorm[v][t]);
    }
    fprintf(File2Write,"\n");
  }
  fclose(File2Write);
  File2Write = fopen("TensAngle.dat","w");
  fprintf(File2Write,"#Press DensPhob DensPhil SponCurv SaddleSplay \n");
  for(int v=0;v<NEdge;v++){
    fprintf(File2Write,"%lf ",v/(double)NEdge*pEdge(CNorm));
    for(int t=0;t<NType;t++){
      TensAng[v][t] /= NormNorm[v][t] > 0. ? VolElm*NormNorm[v][t] : 1.;
      fprintf(File2Write,"%lf ",TensNorm[v][t]);
    }
    fprintf(File2Write,"\n");
  }
  fclose(File2Write);
  for(int v=0;v<NBin;v++){
    free(TensRad[v]);
    free(NormRad[v]);
  }
  for(int v=0;v<NEdge;v++){
    free(TensNorm[v]);
    free(NormNorm[v]);
  }
  free(TensRad);
  free(NormRad);
  free(TensNorm);
  free(NormNorm);
  free(VolContr);
  //free(Plot);
}
int ElPoly::PressTrace(){
  int NType = 3;
  double *Sum = (double *)calloc(NType*CUBE(NEdge),sizeof(double));
  if(NFile[1]-NFile[0] != 3){
    printf("Just three files\n");
    return 1;
  }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(Open(cFile[f],BF_NO))return 1;
    for(int p=0;p<pNPart();p++){
      int vx = (int)((pPos(p,CLat1)+.1)/pEdge(CLat1)*NEdge);
      int vy = (int)((pPos(p,CLat2)+.1)/pEdge(CLat2)*NEdge);
      int vz = (int)((pPos(p,CNorm)+.1)/pEdge(CNorm)*NEdge);
      if(vx < 0 || vx >= NEdge) continue;
      if(vy < 0 || vy >= NEdge) continue;
      if(vz < 0 || vz >= NEdge) continue;
      int vTot = (vx*NEdge+vy)*NEdge+vz;
      Sum[vTot*NType+0] += pVel(p,0)/3.;
      Sum[vTot*NType+1] += pVel(p,1)/3.;
      Sum[vTot*NType+2] += pVel(p,2)/3.;
    }
  }
  printf("\n");
  FILE *FileToWrite = fopen("PressTrace.xvl","w");
  fprintf(FileToWrite,"# l(%.1f %.1f %.1f) r(%.2f %.2f %.2f) v[%d] d[color]\n",pEdge(0),pEdge(1),pEdge(2),Nano->Pos[0],Nano->Pos[1],Nano->Pos[2],NEdge);
  for(int vx=0;vx<NEdge;vx++){
    double x = vx*pEdge(CLat1)/(double)NEdge;
    for(int vy=0;vy<NEdge;vy++){
      double y = vy*pEdge(CLat2)/(double)NEdge;
      for(int vz=0;vz<NEdge;vz++){
	double z = vz*pEdge(CNorm)/(double)NEdge;
	int v = (vx*NEdge+vy)*NEdge+vz;
	if(Sum[v*NType+0] + Sum[v*NType+1] + Sum[v*NType+2] < 0.1) continue;
	fprintf(FileToWrite,"{x(%.2f %.2f %.2f)",x,y,z);
	fprintf(FileToWrite," v( %lf %.2f %.2f)}\n",Sum[v*NType+0],Sum[v*NType+1],Sum[v*NType+2]);
      }
    }
  }
  fclose(FileToWrite);
  free(Sum);
  //free(Plot);
}
void ElPoly::RestPress(int NBin){
  int NType = 3;
  double Round = 0.00001;
  double InvNBin = 1./(double)NBin;
  double **Plot1  = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot1[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  double *Count1 = (double *)calloc(NType*NBin*NBin,sizeof(double));
  double *Count2 = (double *)calloc(NType*NBin*NBin,sizeof(double));
  double LatDim[3] = {pEdge(CLat1),pEdge(CLat2),pEdge(CNorm)};
  double InvLatDim[2] = {1./LatDim[0],1./LatDim[1]};
  double FirstCm[3] = {pCm(CLat1),pCm(CLat2),pCm(CNorm)};
  for(int p=0;p<pNPart();p++){
    int vr = (int)((pPos(p,CLat1)+Round)*InvLatDim[0]*NBin);
    if (vr < 0 || vr >= NBin) continue;
    int vz = (int)((pPos(p,CLat2)+Round)*InvLatDim[1]*NBin);
    if (vz < 0 || vz >= NBin) continue;
    int t = Pm[p].Typ;
    for(int t=0;t<3;t++){
      Plot1[t][(vr*NBin+vz)] += pVel(p,t);//pPos(p,2);
      Count1[(vr*NBin+vz)*NType+t] += 1.;
    }
  }
  FILE *FTecPlot = fopen("TecPlotPressDiff.dat","w");
  fprintf(FTecPlot,"VARIABLES = \"R\", \"Z\", \"diff\", \"d1\",\"d2\"\n");
  fprintf(FTecPlot,"ZONE J=%d, K=%d, F=POINT\n",NBin,NBin);
  for(int vx=0;vx<NBin;vx++){
    for(int vy=0;vy<NBin;vy++){
      double r = vx*pEdge(0)/(double)NBin;
      double z = vy*pEdge(1)/(double)NBin - pEdge(1)*.5;
      //fprintf(FTecPlot,"%lf %lf %lf %lf %lf\n",r,z,diff0,diff1,diff2);
      double Diff = Plot1[0][vx*NBin+vy] - Plot1[1][vx*NBin+vy] - Plot1[2][vx*NBin+vy];
      fprintf(FTecPlot,"%lf %lf %lf %lf %lf\n",r,z,Diff,Plot1[1][vx*NBin+vy],Plot1[2][vx*NBin+vy]);
    }
  }
  fclose(FTecPlot);
  //difference
  FILE *Difference = fopen("PressDifference.xvl","w");
  fprintf(Difference,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NBin,ChooseDraw(EL_COLOR));
  int link[4] = {0,0,0,0};
  for(int t=0,p=0,c=0;t<NType;t++){
    for(int vx=0;vx<NBin-1;vx++){
      for(int vy=0;vy<NBin-1;vy++){
	double Diff = Plot1[0][vx*NBin+vy] - Plot1[1][vx*NBin+vy] - Plot1[2][vx*NBin+vy];
	fprintf(Difference,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf)}\n",p,c,t,vx*InvNBin*pEdge(0),vy*InvNBin*pEdge(1),0.,Diff,Plot1[1][vx*NBin+vy],Plot1[2][vx*NBin+vy]);
	p++;
      }
    }
  }
  fclose(Difference);
  for(int t=0;t<NType;t++){
    free(Plot1[t]);
  }
  free(Plot1);
  free(Count1);
}
int ElPoly::Diff2Files(int NBin,int How){
  if(NFile[1] - NFile[0] != 2){
    printf("The number of files must be two\n");
    return 1;
  }
  int NType = 3;
  double Round = 0.00001;
  NBin = MIN(NBin,NEdge);
  double InvNBin = 1./(double)NBin;
  double **Plot1  = (double **)calloc(NType,sizeof(double));
  double **Plot2  = (double **)calloc(NType,sizeof(double));
  double **PlotDiff  = (double **)calloc(NType,sizeof(double));
  for(int t=0;t<NType;t++){
    Plot1[t] = (double *)calloc(NBin*NBin,sizeof(double));
    Plot2[t] = (double *)calloc(NBin*NBin,sizeof(double));
    PlotDiff[t] = (double *)calloc(NBin*NBin,sizeof(double));
  }
  double *Count1 = (double *)calloc(NType*NBin*NBin,sizeof(double));
  double *Count2 = (double *)calloc(NType*NBin*NBin,sizeof(double));
  double LatDim[3] = {pEdge(CLat1),pEdge(CLat2),pEdge(CNorm)};
  double InvLatDim[2] = {1./LatDim[0],1./LatDim[1]};
  double FirstCm[3] = {pCm(CLat1),pCm(CLat2),pCm(CNorm)};
  for(int p=0;p<pNPart();p++){
    int vr = (int)((pPos(p,CLat1)+Round)*InvLatDim[0]*NBin);
    if(vr < 0 || vr >= NBin) continue;
    int vz = (int)((pPos(p,CLat2)+Round)*InvLatDim[1]*NBin);
    if(vz < 0 || vz >= NBin) continue;
    if(How==0){//Press
      for(int t=0;t<3;t++){
	Plot1[t][vr*NBin+vz] += pVel(p,t);
	Count1[(vr*NBin+vz)*NType+t] += 1.;
      }
    }
    else{//Dens
      int t = Pm[p].Typ;
      Plot1[t][vr*NBin+vz] += pPos(p,2);
      Count1[(vr*NBin+vz)*NType+t] += 1.;
    }
  }
  if(Open(cFile[NFile[1]-1],BF_NO))return 1;
  for(int p=0;p<pNPart();p++){
    int vr = (int)((pPos(p,CLat1)+Round+1.)*InvLatDim[0]*NBin);
    if (vr < 0 || vr >= NBin) continue;
    int vz = (int)((pPos(p,CLat2)+Round)*InvLatDim[1]*NBin);
    if (vz < 0 || vz >= NBin) continue;
    if(How==0){//Press
      for(int t=0;t<3;t++){
	Plot2[t][vr*NBin+vz] += pVel(p,t);
	Count2[(vr*NBin+vz)*NType+t] += 1.;
      }
    }
    else{//Dens
      int t = Pm[p].Typ;
      Plot2[t][vr*NBin+vz] += pPos(p,2);
      Count2[(vr*NBin+vz)*NType+t] += 1.;
    }
  }
  Matrice Mask(5,5);
  Mask.FillGaussian(.5,3.);
  Mask.Print();
  for(int t=0;t<NType;t++){
    //ConvoluteMatrix(Plot1[t],NBin,&Mask,2);
    // ConvoluteMatrix(Plot1[t],NBin,&Mask,2);
    //ConvoluteMatrix(Plot2[t],NBin,&Mask,2);
    // ConvoluteMatrix(Plot2[t],NBin,&Mask,2);
  }
  for(int vr=0;vr<NBin;vr++){
    for(int vz=0;vz<NBin;vz++){
      for(int t=0;t<NType;t++){
	Plot1[t][(vr*NBin+vz)] /= Count1[(vr*NBin+vz)*NType+t] > 0. ? Count1[(vr*NBin+vz)*NType+t] : 1.;
	Plot2[t][(vr*NBin+vz)] /= Count2[(vr*NBin+vz)*NType+t] > 0. ? Count2[(vr*NBin+vz)*NType+t] : 1.;
	PlotDiff[t][(vr*NBin+vz)] = Plot1[t][(vr*NBin+vz)] - Plot2[t][(vr*NBin+vz)];
      }
      int t = 0;
    }
  }
  int NDim = 2;
  int IfMinImConv = 0;
  for(int t=0;t<NType;t++){
    Mask.ConvoluteMatrix(PlotDiff[t],NBin,NDim,IfMinImConv);
    // Mask.ConvoluteMatrix(PlotDiff[t],NBin,2);
  }
  //tecplot
  FILE *FTecPlot = fopen("TecPlotDiff.dat","w");
  fprintf(FTecPlot,"VARIABLES = \"R\", \"Z\", \"diff\", \"d1\",\"d2\"\n");
  fprintf(FTecPlot,"ZONE J=%d, K=%d, F=POINT\n",NBin,NBin);
  //difference
  FILE *Difference;
  if(How) Difference = fopen("DensDifference.rzd","w");
  else Difference = fopen("PressDifference.xvl","w"); 
  if(How==1){//Dens
    PrintDens(Difference,PlotDiff,LatDim,NBin);
    for(int vx=0;vx<NBin;vx++){
      for(int vy=0;vy<NBin;vy++){
	double r = vx*LatDim[0]/(double)NBin;
	double z = vy*LatDim[1]/(double)NBin - LatDim[1]*.5;
	fprintf(FTecPlot,"%lf %lf %lf %lf %lf\n",r,z,PlotDiff[0][vx*NBin+vy],Plot1[0][vx*NBin+vy],Plot2[0][vx*NBin+vy]);
      }
    }
  }
  else{//Pre
    fprintf(Difference,"# l(%.1f %.1f %.1f) v[%d] d[color]\n",LatDim[0],LatDim[1],LatDim[2],NBin);
    for(int vr=0;vr<NBin;vr++){
      for(int vz=0;vz<NBin;vz++){
	double NanoAdded = Plot1[2][(vr*NBin+vz)];
	double Phob = PlotDiff[0][(vr*NBin+vz)];
	double Phil = PlotDiff[1][(vr*NBin+vz)];
	double r = (vr)*InvNBin*LatDim[0];
	double z = (vz)*InvNBin*LatDim[1];
	double dens = (PlotDiff[2][(vr*NBin+vz)]);
	fprintf(Difference,"{x(%lf %lf %lf) v(%lf %lf %lf)}\n",r,z,dens,NanoAdded,Phob,Phil);
	double z1 = vz*LatDim[1]/(double)NBin - LatDim[1]*.5;
	fprintf(FTecPlot,"%lf %lf %lf %lf %lf\n",r,z1,PlotDiff[0][vr*NBin+vz],Plot1[1][vr*NBin+vz],Plot2[1][vr*NBin+vz]);
      }
    }
  }
  fclose(FTecPlot);
  fclose(Difference);
  for(int t=0;t<NType;t++){
    free(Plot1[t]);
    free(Plot2[t]);
  }
  free(Plot1);
  free(Plot2);
  free(Count1);
  free(Count2);
  return 0;
}
