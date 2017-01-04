#include "ElPoly.h"
int ElPoly::ProjectionF(int NBin,int Coord){
  if(Coord > 4 || Coord <0) return 1;
  int NType = 5;
  double *Plot = (double *)calloc(NBin*NBin*NType,sizeof(double));
  double InvNBin = 1./(double)NBin;
  double RefPos[3] = {0.,0.,0.};
  for(int d=0;d<3;d++){
    RefPos[d] = Nano->Pos[d]-.5*pEdge(d);
  }
  if(Coord == 3){
    RefPos[0]=pCm(0);RefPos[1]=pCm(1);RefPos[2]=pCm(2);
  }
  SetEdge(.5*MIN(pEdge(CLat1),pEdge(CLat2)),3);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    OpenRisk(cFile[f],BF_PART);
    ShiftSys(SHIFT_CM);
    int NPlot = 0;
    //---Projects-against-one-coordinate--
    if(Coord < 3){
      int coord1 = (Coord+1)%3;
      int coord2 = (Coord+2)%3;
      for(int p=0;p<pNPart();p++){
	double x = pPos(p,coord1) - RefPos[coord1];
	x -= floor(x*pInvEdge(coord1))*pEdge(coord1);
	double y = pPos(p,coord2) - RefPos[coord2];
	y -= floor(y*pInvEdge(coord2))*pEdge(coord2);
	int v = (int)(NBin*x*pInvEdge(coord1));
	if( v < 0 || v >= NBin) continue;
	int vv = (int)(NBin*y*pInvEdge(coord2));
	if( vv < 0 || vv >= NBin) continue;
	int t = pType(p);
	if( t < 0 || t > 3) continue;
	if( CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED) )
	  Plot[(v*NBin+vv)*NType+3] += 1.;
	Plot[(v*NBin+vv)*NType+t] += 1.;
	if(p<pNPart()-1)
	  if(pType(p+1) == 1 && pType(p) == 0)
	    Plot[(v*NBin+vv)*NType+4] += 1.;	  
      }
    }
    //---Projects-against-the-radial-coordinate--
    else if(Coord == 3){
      SetEdge(.5*MAX((pEdge(CLat1)),(pEdge(CLat2))),3);
      for(int p=0;p<pNPart();p++){
	double x = pPos(p,CLat1) - RefPos[CLat1];
	x -= floor(x*pInvEdge(CLat1))*pEdge(CLat1);
	double y = pPos(p,CLat2) - RefPos[CLat2];
	y -= floor(y*pInvEdge(CLat2))*pEdge(CLat2);
	double z = pPos(p,CNorm) - RefPos[CNorm];
	z -= floor(z*pInvEdge(CNorm))*pEdge(CNorm);
	double r = sqrt(SQR(x)+SQR(y));
	int v = (int)(NBin*r*pInvEdge(3));
	if( v < 0 || v >= NBin) continue;
	int vv = (int)(NBin*pPos(p,CNorm)/pEdge(CNorm));
	if( vv < 0 || vv >= NBin) continue;
	int t = pType(p);
	if( t < 0 || t > 3) continue;
	if( CHAIN_IF_TYPE(Ch[pChain(p)].Type,CHAIN_ADDED) )
	  Plot[(v*NBin+vv)*NType+3] += 1.;
	Plot[(v*NBin+vv)*NType+t] += 1.;
	if(p<pNPart()-1)
	  if(pType(p+1) == 1 && pType(p) == 0)
	    Plot[(v*NBin+vv)*NType+4] += 1.;	  
      }
    }
  }
  printf("\n");
  //-----writes-the-output-file-------------------
  FILE *FileToWrite = NULL;
  FileToWrite = fopen("Projection.xyd","w");
  fprintf(FileToWrite,"#l(%lf %lf %lf) v[%d] d[%s]\n",pEdge(CLat1),pEdge(CLat2),pEdge(CNorm),NBin,ChooseDraw(EL_QUAD1));
  int coord1 = (Coord+1)%3;
  int coord2 = (Coord+2)%3;
  if(Coord == 3){
    coord1 = 3;
    coord2 = CNorm;
  }
  double Max[NType];
  for(int t=0;t<NType;t++){
    Max[t] = Plot[t];
    for(int v=0;v<NBin;v++)
      for(int vv=0;vv<NBin;vv++)
	if(Max[t] < Plot[(v*NBin+vv)*NType+t]) Max[t] = Plot[(v*NBin+vv)*NType+t];
    Max[t] = Max[t] <= 0. ? 1. : Max[t];
  }
  //for(int t=0;t<NType-1;t++){
  for(int t=0;t<1;t++){
    for(int v=0;v<NBin;v++){
      for(int vv=0;vv<NBin;vv++){
	int p = (v*NBin+vv)*NType+t;
	int c = 0;
	if(Plot[p] < .1) continue;
	double x = (v)*InvNBin*pEdge(CLat1);
	double y = (vv)*InvNBin*pEdge(CLat2);
	double dens = Plot[p]/Max[t]*5.+.5*pEdge(CNorm);
	double NanoAdded = 0.;//Plot[p]/Max[t]+Plot[((v*NBin+vv)*NType+3]/Max[3];
	double Phob = t == 0 ? Plot[(v*NBin+vv)*NType+0]/Max[0] : 0.;
	double Phil = t == 1 ? Plot[(v*NBin+vv)*NType+1]/Max[1] : 0.;
	fprintf(FileToWrite,"{t[%d %d %d] x(%lf %lf %lf) v(%lf %lf %lf)}\n",p,c,t,x,y,dens,NanoAdded,Phob,Phil);
      }
    } 
  }
  free(Plot);
  fclose(FileToWrite);
  return 0;
}
// int ElPoly::LateralDistr(int NBin){
// }
int ElPoly::Surface(int NBin,int Threshold){
  double *PlotAv = (double *)calloc(NBin*NBin,sizeof(double));
  int NPoint=0;
  double AreaAv = 0.,SurfAv = 0.;
  double Max=0.;
  printf("Edge %lf Threshold %d\n",pEdge(3),Threshold);
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return 1;
    double Surf=0.,Area =0.;
    double *Plot = (double *)calloc(NBin*NBin,sizeof(double));
    for(int p=0;p<pNPart();p++){
      int v = (int)(pPos(p,CLat1)/pEdge(CLat1)*NBin);
      if( v < 0 || v > NBin){ printf("Oi 0 < %d < %d\n",v,NBin);return 1;} 
      int vv = (int)(pPos(p,CLat2)/pEdge(CLat2)*NBin);
      if( vv < 0 || vv > NBin){ printf("Oi 0 < %d < %d\n",vv,NBin);return 1;} 
      Plot[v*NBin+vv] += 1.;
    }
    for(int v=0;v<NBin;v++)
      for(int vv=0,n=0;vv<NBin;vv++)
	PlotAv[v*NBin+vv] += Plot[v*NBin+vv];
    for(int v=1;v<NBin-1;v++){
      for(int vv=1,n=0;vv<NBin-1;vv++){
	if(Plot[v*NBin+vv] > Threshold){
	  NPoint++;
	  //printf("(%d %d) %d %d %d %d\n",v,vv,Plot[v-1][vv],Plot[v+1][vv],Plot[v][vv-1],Plot[v][vv+1]);
	  if(Plot[(v-1)*NBin+vv] > Threshold && v > 0)
	    n++;
	  if(Plot[(v+1)*NBin+vv] > Threshold && v < NBin -1)
	    n++;
	  if(Plot[v*NBin+vv-1] > Threshold && vv > 0)
	    n++;
	  if(Plot[v*NBin+vv+1] > Threshold && vv < NBin -1)
	    n++;
	  if(n == 4){
	    Area += 1.;
	  }
	  else if(n < 4  && n != 0){
	    if(n == 3){
	      Surf += 1.;
	      Area += 1.;
	    }
	    else if(n == 2){
	      Surf += 2;
	    }
	    else if(n == 1){
	      Surf += sqrt(2);
	    }	   
	  }
	  n = 0;
	}
      }
    }
    AreaAv += Area;
    SurfAv += Surf;
    free(Plot);
  }
  printf("\n");
  //-----------normalize-------------------
  for(int v=0;v<NBin;v++)
    for(int vv=0,n=0;vv<NBin;vv++)
      if(PlotAv[v*NBin+vv] > Max)
	Max = PlotAv[v*NBin+vv];
  //------------write-output-----------------
  FILE *FileToWrite = NULL;
  FileToWrite = fopen("Surface.xyz","w");
  double FNorma = 1./(double)(NFile[1] - NFile[0]);
  double LengthConv = pEdge(0)*pEdge(1)/(double)(NBin*NBin);
  fprintf(FileToWrite,"#AreaTot %lf Area %lf Surf %lf Threshold %d NChain %d NBin %d\n",pEdge(CLat1)*pEdge(CLat2),AreaAv*FNorma*LengthConv,SurfAv*FNorma*LengthConv,Threshold,pNChain(),NBin);
  for(int v=0;v<NBin;v++)
    for(int vv=0,n=0;vv<NBin;vv++)
      if(PlotAv[v*NBin+vv] > 0.)
	fprintf(FileToWrite,"%lf %lf %lf\n",(double)v/NBin*pEdge(CLat1),(double)vv/NBin*pEdge(CLat2),PlotAv[v*NBin+vv]/Max);
  fclose(FileToWrite);
  free(PlotAv);
  return 0;
}
int ElPoly::From3To2d(int NSample,double Param){
  char FName[120];
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NANO));
    sprintf(FName,"ProjOnSurf%05d.dat",f);
    FILE *FLine = fopen(FName,"w");
    fprintf(FLine,"# l(%lf %lf %lf) d[part]\n",1.,pEdge(CLat2)*pInvEdge(CLat1),.5);
    for(int p=0;p<pNPart();p++){
      if(pType(p) != 0) continue;
      fprintf(FLine,"{x(%lf %lf %lf) v(%lf %lf %lf)}\n",pPos(p,CLat1)*pInvEdge(CLat1),pPos(p,CLat2)*pInvEdge(CLat1),.5,pVel(p,0),pVel(p,1),pVel(p,2));
    }
    fclose(FLine);
  }
  printf("\n");
}
int ElPoly::From2To1d(int Coord){
  if(Coord != 0 && Coord != 1){
    printf("Coordinate accepted are 0 or 1 not %d\n",Coord);
    return 1;
  }
  int SecCoord = (Coord+1)%2;
  double *Line = (double *) calloc(3*NEdge,sizeof(double));
  double *Count = (double *) calloc(3*NEdge,sizeof(double));
  for(int p=0;p<pNPart();p++){
    int vz = (int)((pPos(p,Coord)+0.001)*pInvEdge(Coord)*NEdge);
    if(vz < 0 || vz >= NEdge) continue;
    for(int d=0;d<3;d++){
      Line[vz*3+d] += pVel(p,d);
      Count[vz*3+d] += 1.;
    }
  }
  for(int vz=0;vz<NEdge;vz++){
    for(int d=0;d<3;d++){
      double Norm = Count[vz*3+d] > 0. ? Count[vz*3+d] : 1.;
      //Line[vz*3+d] /= Norm;
      Line[vz*3+d] *= pEdge(SecCoord)/(double)NEdge;
    }
  }
  FILE *FLine = fopen("ProjOnLine.dat","w");
  for(int v=0;v<NEdge;v++){
    fprintf(FLine,"%lf %lf %lf %lf\n",v*pEdge(Coord)/(double)NEdge,Line[v*3+0],Line[v*3+1],Line[v*3+2]);
  }
  fclose(FLine);
  free(Line);
  free(Count);
}
int ElPoly::From3To1d(int Coord){


}
int ElPoly::RadialShell(int NBin){
  double Volume=0;//Global constant
  double Area=0.;
  double Norm=0.;
  double **Plot;
  double Ypsilon = 0.;
  double InvNFile  = 1./(double)(NFile[1]-NFile[0]);
  Plot = (double **)calloc(NBin,sizeof(double));
  for(int i=0;i<NBin;i++){
    *(Plot+i) = (double *)calloc(NBin,sizeof(double));
  }
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART))return 0;
    SetEdge(MIN((pEdge(CLat1)*.5),(pEdge(CLat2)*.5)),3);
    for(int p=0;p<pNPart();p++){
      double Rad = sqrt(SQR((pPos(p,CLat1)-pCm(CLat1))) 
			+ SQR((pPos(p,CLat2)-pCm(CLat2))) );
      int vr = (int)(NBin*Rad*pInvEdge(3));
      int vz = (int)(NBin*pPos(p,CNorm)*pInvEdge(CNorm));
      //printf("%lf %lf %d %d \n",Rad,pPos(p,CNorm),v,vv);
      if( vr < 0 || vr >= NBin) continue;
      if( vz < 0 || vz >= NBin) continue;
      Plot[vr][vz] += InvNFile;
    }
    Ypsilon += pEdge(2)-pCm(2)-1.;
  }
  printf("\n");
  FILE *FileToWrite = fopen("RadialShell.xye","w");
  //fprintf(FileToWrite,"# l(%.2f %.2f %.2f) v[60] d[chain]\n",Gen->Edge[0],Gen->Edge[1],Gen->Edge[2]);
  fprintf(FileToWrite,"%lf %lf %lf\n",0.,0.,0.);
  fprintf(FileToWrite,"%lf %lf %lf\n",pEdge(0),pEdge(1),1.);
  double Max=0.;
  for(int i=0;i<NBin;i++){
    for(int j=0;j<NBin;j++){
      if(Plot[i][j] > Max) Max = Plot[i][j];
    }
  }
  int th=0;
  for(int vz=0;vz<NBin;vz++){
    //for(int vr=NBin-1;vr>=0;vr--){
    for(int vr=0;vr<NBin;vr++){
      if(Plot[vr][vz] > 0.){
	fprintf(FileToWrite,"%lf %lf %lf\n",(double)vr/NBin*pEdge(3),(double)vz/NBin*pEdge(2),Plot[vr][vz]/Max);
	Norm += Plot[vr][vz];
	th++;
	if(th > 4){
	  th = 0;
	  break;
	}
      }
    }
  }
  fclose(FileToWrite);
  return 0;
  Mat->Ypsilon = Ypsilon*InvNFile;
  double Vol = 1.;//pNPart()/(Gen->rho/(double)pNPCh()*CUB(5.));
  Mat->PreFact = 3./8.*pow((3.*Vol)/DUE_PI,1./3.);
  Mat->Func = &Matematica::ContactAngle;
  int NRadici = 4;
  printf("Volume  %lf Cm %lf Area %lf # to Invert %lf\n",(double)pNPart()/10.,pCm(2),Area,pow(DUE_PI/(2.*pNPart()/10.),.25));
  double *Radici;
  Radici = (double *)malloc(NRadici*sizeof(double));
  int NZeri = Mat->Zeri(0.,DUE_PI/2.,Radici,NRadici);
  for(int i=0;i<NZeri;i++){
    Radici[i] *= 360./DUE_PI;
    printf("Angle %lf\n",Radici[i]);
  }
  if(NZeri == 0){
    printf("The function has no real roots\n");
  }
  free(Plot);
  return 0;
}
void ElPoly::IsoSurf(int NSample,double *IsoLevel,int NIso){
  int NPairF = NFile[1]-NFile[0];
  double OldPos[3] = {pNanoPos(0,0),pNanoPos(0,1),pNanoPos(0,2)};
  double DensEl = CUB(NSample)/(pVol()*NPairF);
  double Dens  = 13.3;
  FILE *FNano = fopen("NanoPos.txt","w");
  //IsoLevel *= NPairF;
  for(int ff=NFile[0];ff<NFile[1];ff+=NPairF){
    double *Plot = (double *)calloc(CUBE(NSample),sizeof(double));
    double Min = 0.;
    double Max = 0.;
    VAR_TRIANGLE *Triang = NULL;
    double Pos[3];
    for(int f=ff;f<ff+NPairF&&f<NFile[1];f++){
      Processing(f);
      if(OpenRisk(cFile[f],BF_PART))return;
      for(int b=0;b<pNBlock();b++){
	for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	  if(pType(p) != 0)continue;
	  for(int d=0;d<3;d++){
	    Pos[d] = pPos(p,d) - (pNanoPos(0,d)-.5*pEdge(d));
	    Pos[d] -= floor(Pos[d]*pInvEdge(d))*pEdge(d);
	  }
	  int sx = (int)(Pos[0]*pInvEdge(0)*NSample);
	  int sy = (int)(Pos[1]*pInvEdge(1)*NSample);
	  int sz = (int)(Pos[2]*pInvEdge(2)*NSample);
	  int sTot = (sx*NSample+sy)*NSample+sz;
	  Plot[sTot] += DensEl;
	  if(Max < Plot[sTot]) Max = Plot[sTot];
	  if(Min > Plot[sTot]) Min = Plot[sTot];
	}
      }
    }
    Matrice Mask(3,3,3);
    Mask.FillGaussian(.5,3.);
    Mask.Print();
    int NDim = 3;
    int IfMinImConv = 1;
    Mask.ConvoluteMatrix(Plot,NSample,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot,NSample,NDim,IfMinImConv);
    // ConvoluteMatrix(Plot,NSample,&Mask,3);
    // ConvoluteMatrix(Plot,NSample,&Mask,3);
    char FString[256];
    sprintf(FString,"IsoSurf%05d.dat",ff);
    FILE *F2Write = fopen(FString,"w");
    fprintf(F2Write,"#l(%lf %lf %lf) v[%d] d[polygon]\n",pEdge(0),pEdge(1),pEdge(2),NSample);
    HeaderNano(F2Write);
    int NTri = 0;
    for(int i=0;i<NIso;i++){
      printf("Min %lf Max %lf IsoLevel %lf DensEl %lf\n",Min,Max,IsoLevel[i],DensEl);
      Triang = MarchingCubes(Plot,NSample,IsoLevel[i],&NTri);
      for(int t=0;t<NTri;t++){
	for(int i=0;i<3;i++){
	  int l1 = t*3 + (i+1)%3;
	  int l2 = t*3 + (i+2)%3;
	  for(int d=0;d<3;d++) Pos[d] = Triang[t].p[i].x[d];
	  int sx = (int)(Pos[0]*pInvEdge(0)*NSample);
	  int sy = (int)(Pos[1]*pInvEdge(1)*NSample);
	  int sz = (int)(Pos[2]*pInvEdge(2)*NSample);
	  int sTot = (sx*NSample+sy)*NSample+sz;
	  fprintf(F2Write,"{t[%d %d %d] ",sTot,t,0);
	  fprintf(F2Write,"x(%lf %lf %lf) ",Pos[0],Pos[1],Pos[2]);
	  fprintf(F2Write,"v(%lf %lf %lf) ",Triang[t].n[i].x[0],Triang[t].n[i].x[1],Triang[t].n[i].x[2]);
	  fprintf(F2Write,"l[%d] l[%d]}\n",l1,l2);
	}
      }
    }
    fclose(F2Write);
    free(Plot);
    free(Triang);continue;
    int NVertex = CUBE(2*NSample-1);
    double *WeightL = (double *) calloc(NVertex,sizeof(double));
    NormalWeight(Triang,WeightL,NSample,NTri);
    double CmStalk[3] = {0.,0.,0.};//OldPos[0],OldPos[1],OldPos[2]};
    double CountStalk = 0.;
    for(int t=0;t<NTri;t++){
      for(int v=0;v<3;v++){
	int vRef = Triang[t].v[v];
	for(int d=0;d<3;d++){
	  CmStalk[d] = Triang[t].p[v].x[d]*WeightL[vRef];
	}
	CountStalk += WeightL[vRef];
      }
    }
    free(WeightL);
    free(Triang);
    if(CountStalk <= 0.) CountStalk = 1.;
    for(int d=0;d<3;d++){
      CmStalk[d] /= CountStalk;
    }
    pPos(CmStalk);
    SetNNano(1);
    for(int d=0;d<3;d++){
      Nano->Pos[d]   = CmStalk[d];
      OldPos[d] = Nano->Pos[d];
    }
    SetNanoBkf(0);
    Nano->Axis[0] = 0.;
    Nano->Axis[1] = 0.;
    Nano->Axis[2] = 1.;
    Nano->Rad     = .5;
    Nano->Height  = 4.;
    Nano->Hamaker = 1.;
    Nano->OffSet  = 1.;
    Nano->Shape = SHAPE_STALK;    
    for(int f=ff;f<ff+NPairF&&f<NFile[1];f++){
      char String[120];
      StringNano(String,0);
      fprintf(FNano,"sed '/Rigid/c\\%s' %s > Ciccia.dat\n",String,cFile[f]);
      fprintf(FNano,"mv Ciccia.dat %s\n",cFile[f]);
      //command("chmod u+x %s\n","NanoPos.txt");
      //SubNanoHeader(cFile[f]);
    }
    printf("\n");
  }
  fclose(FNano);
}
void ElPoly::IsoLine(int NSample,double *IsoLevel,int NIso,int How){
  int NPairF = 1;//NFile[1]-NFile[0];
  double OldPos[3] = {pNanoPos(0,0),pNanoPos(0,1),pNanoPos(0,2)};
  double DensEl = SQR(NSample)/(pVol()*NPairF);
  double Dens  = 13.3;
  //IsoLevel *= NPairF;
  for(int ff=NFile[0];ff<NFile[1];ff+=NPairF){
    double *Plot = (double *)calloc(SQR(NSample),sizeof(double));
    double Min = 0.;
    double Max = 0.;
    double Pos[3];
    for(int f=ff;f<ff+NPairF&&f<NFile[1];f++){
      Processing(f);
      if(OpenRisk(cFile[f],BF_PART))return;
      for(int b=0;b<pNBlock();b++){
	for(int p=Block[b].InitIdx;p<Block[b].EndIdx;p++){
	  if(pType(p) != 0)continue;
	  for(int d=0;d<3;d++){
	    Pos[d] = pPos(p,d) - (pNanoPos(0,d)-.5*pEdge(d));
	    Pos[d] -= floor(Pos[d]*pInvEdge(d))*pEdge(d);
	  }
	  int sx = (int)(Pos[CLat1]*pInvEdge(CLat1)*NSample);
	  int sy = (int)(Pos[CLat2]*pInvEdge(CLat2)*NSample);
	  int sTot = sx*NSample+sy;
	  if(How==0)//particle file
	    Plot[sTot] += DensEl;
	  else
	    Plot[sTot] += pPos(p,CNorm);
	  if(Max < Plot[sTot]) Max = Plot[sTot];
	  if(Min > Plot[sTot]) Min = Plot[sTot];
	}
      }
      printf("Min %lf Max %lf DensEl %lf\n",Min,Max,DensEl);
    }
    Matrice Mask(5,5);
    Mask.FillGaussian(.5,3.);
    Mask.Print();
    int NDim = 2;
    int IfMinImConv = 1;
    Mask.ConvoluteMatrix(Plot,NSample,NDim,IfMinImConv);
    Mask.ConvoluteMatrix(Plot,NSample,NDim,IfMinImConv);
    char FString[60];
    sprintf(FString,"IsoSurf%05d.dat",ff);
    FILE *F2Write = fopen(FString,"w");
    fprintf(F2Write,"#l(%lf %lf %lf) v[%d] d[part]\n",pEdge(0),pEdge(1),pEdge(2),NSample);
    HeaderNano(F2Write);
    IsoLine(F2Write,Plot,NSample,IsoLevel,NIso);
    free(Plot);
    fclose(F2Write);
  }
}
void ElPoly::IsoLine(FILE *F2Write,double *Plot,int NSample,double *IsoLevel,int NIso){
  VAR_LINE *Triang = NULL;
  int NTri = 0;
  int p=0;
  char FName[60];
  double Pos[3];
  for(int i=0;i<NIso;i++){
    Triang = MarchingSquares(Plot,NSample,IsoLevel[i],&NTri);
    ConnectLineChain(Triang,NSample,NTri);
    sprintf(FName,"IsoLineChain%d.dat",i);
    sprintf(cWhat2Draw,"part");
    Write(FName);
    // for(int c=0;c<pNChain();c++){
    // 	sprintf(FString,"Line%05dChain%02d.dat",ff,c);
    // 	FILE *FLine = fopen(FString,"w");
    // while(
    // for(int p=0;p<pNChain();p++){
    //   if
    //   fprintf(FString,"%lf %lf\n"
    // 	fclose(FLine);
    // }
    if(1==1){//Isoline without ordering
      for(int t=0;t<NTri;t++){
	for(int v=0;v<2;v++){
	  for(int d=0;d<3;d++) Pos[d] = Triang[t].p[v].x[d];
	  int sTot = Triang[t].v[v];
	  fprintf(F2Write,"{t[%d %d %d] ",sTot,t,0);
	  fprintf(F2Write,"x(%lf %lf %lf) ",Pos[0],Pos[1],Pos[2]);
	  if(v==0)
	    fprintf(F2Write,"l[%d]}\n",p+1);
	  else
	    fprintf(F2Write,"\n");
	  p++;
	}
      }
    }
    if(1==0){//Isoline separated between inside and outside an elips
      double Radx = 20.;
      double Rady = 12.;
      double Center[3] = {0.,.5*pEdge(CLat2),0.};
      int NIt = 100;
      double DistElips = SQR(Radx)+SQR(Rady);
      for(int i=0;i<NIt;i++){
	double Ang = i/(double)NIt*M_PI - M_PI*.5;
	double x = Center[0] + Radx*cos(Ang);
	double y = Center[1] + Rady*sin(Ang);
	fprintf(F2Write,"{t[%d %d %d] ",p++,2,2);
	fprintf(F2Write,"x(%lf %lf %lf) }\n",x,y,0.);	  
      }
      for(int t=0;t<NTri;t++){
	for(int v=0;v<2;v++){
	  for(int d=0;d<3;d++) Pos[d] = Triang[t].p[v].x[d];
	  int sTot = Triang[t].v[v];
	  int type = 0;
	  if(fabs(Pos[0] - Center[0]) < Radx && fabs(Pos[1] - Center[1]) < Rady) type = 1;
	  fprintf(F2Write,"{t[%d %d %d] ",sTot,t,type);
	  fprintf(F2Write,"x(%lf %lf %lf) ",Pos[0],Pos[1],Pos[2]);
	  if(v==0)
	    fprintf(F2Write,"l[%d]}\n",p+1);
	  else
	    fprintf(F2Write,"\n");
	  p++;
	}
      }
    }
  }
  if(1==1){//to write the continuum density values
    for(int sx=0;sx<NSample;sx++){
      for(int sy=0;sy<NSample;sy++){
	double x = pEdge(CLat1)*sx/(double)NSample;
	double y = pEdge(CLat2)*sy/(double)NSample;
	int sTot = sx*NSample+sy;
	fprintf(F2Write,"{t[%d %d %d] ",sTot,NTri+sTot/2,1);
	fprintf(F2Write,"x(%lf %lf %lf) }\n",x,y,Plot[sx*NSample+sy]);
      }
    }
  }
  if(1==0){//Write the chains separately
    for(int c=0;c<pNChain();c++){
      sprintf(FName,"Chain%d.dat",c);
      FILE *FChain = fopen(FName,"w");
      for(int p=0;p<pNPart();p++){
	if(Pm[p].CId != c) continue;
	fprintf(FChain,"%lf %lf\n",Pm[p].Pos[0],Pm[p].Pos[1]);
      }
    }
  }
  free(Triang);
}
void ElPoly::FetchPore(){
  FILE *FNano = fopen("NanoPos.sh","w");
  double OldPos[5] = {pNanoPos(0,0),pNanoPos(0,1),pNanoPos(0,2),Nano->Rad,Nano->Height};
  FILE *StalkArea = fopen("PoreArea.dat","w");
  fprintf(StalkArea,"#time rad asy\n");
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_PART)) return;
    double Asy = PorePos();
    Nano[0].Shape = ShapeId("pore");
    Nano->Height = .2;
    char String[120];
    StringNano(String,0);
    fprintf(StalkArea,"%lf %lf %lf\n",pTime(),Nano->Rad,Asy);
    fprintf(FNano,"sed '/pore/c\\%s' %s > Ciccia.dat\n",String,cFile[f]);
    fprintf(FNano,"mv Ciccia.dat %s\n",cFile[f]);
    //command("chmod u+x %s\n","NanoPos.txt");
    //SubNanoHeader(cFile[f]);
  }
  fclose(FNano);
  fclose(StalkArea);
  printf("\n");
}
/**
   Simple Monte Carlo to find the best position and radius of the osculating torus.
 */
void ElPoly::FetchStalk(){
  FILE *FNano = fopen("NanoPos.sh","w");
  double OldPos[5] = {pNanoPos(0,0),pNanoPos(0,1),pNanoPos(0,2),Nano->Rad,Nano->Height};
  FILE *StalkArea = fopen("StalkArea.dat","w");
  char FName[60];
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO)) return;
    SetNNano(1);
    if(StalkPos(OldPos)) continue;
    char String[120];
    StringNano(String,0);
    fprintf(StalkArea,"%lf %lf\n",pTime(),Nano->Area);
    fprintf(FNano,"sed '/Rigid/c\\%s' %s > Ciccia.dat\n",String,cFile[f]);
    fprintf(FNano,"mv Ciccia.dat %s\n",cFile[f]);
    sprintf(FName,"Centered%05d.dat",f);
    //Write(FName);
    //command("chmod u+x %s\n","NanoPos.txt");
    //SubNanoHeader(cFile[f]);
  }
  fclose(FNano);
  fclose(StalkArea);
  printf("\n");
}
/**
   Simple Monte Carlo to find the best position and radius of the osculating torus.
   The area and position of the torus are hence redifined counting how many hydrophilic beads are inside the torus.
 */
void ElPoly::StalkArea(){
  FILE *FNano = fopen("NanoPosA.sh","w");
  FILE *StalkArea = fopen("StalkArea.dat","w");
  FILE *AreaS = fopen("AreaS.dat","w");
  double OldPos[5] = {pNanoPos(0,0),pNanoPos(0,1),pNanoPos(0,2),Nano->Rad,Nano->Height};
  int NBin = 36;
  double *Count = (double *)calloc(NBin*NBin,sizeof(double));
  double RadTorus = 1.;//Nano->Rad;
  //fprintf(AreaS,"#l(%lf %lf %lf) \n",2.*Nano->Height,2.*Nano->Height,10.);
  fprintf(AreaS,"#l(%lf %lf %lf) \n",pEdge(0),pEdge(1),pEdge(2));
  HeaderNano(AreaS);
  char FName[60];
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO)) return;
    SetNNano(1);
    if(StalkPos(OldPos)) continue;
    //Nano->Pos[CNorm] = pCm(CNorm);
    //Nano->Pos[CNorm] -= floor(Nano->Pos[CNorm]*pInvEdge(CNorm))*pEdge(CNorm);
    double Cm[3] = {0.,0.,0.};
    double CountCm = 0;
    int nPart = 0;
    double Pos[3];
    //counts the particles inside the torus
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	Pos[d] = pPos(p,d) - Nano->Pos[d];
	Pos[d] -= floor(pPos(p,d)*pInvEdge(d))*pEdge(d);
      }
      double Rad = hypot(Pos[CLat1],Pos[CLat2]);
      if(Rad > Nano->Height) continue;
      if(fabs(Pos[CNorm]) > RadTorus) continue;
      // fprintf(AreaS,"{t[%d 0 %d] x(%lf %lf %lf)",nPart++,pType(p),pPos(p,0),pPos(p,1),pPos(p,2));
      // if(Ln[p].NLink > 0) fprintf(AreaS,"l[%d]",p-Ln[p].Link[0]+nPart+1);
      // fprintf(AreaS,"}\n");
      if(pType(p) == 1) continue;
      for(int d=0;d<3;d++){
	Cm[d] += pPos(p,d);
      }
      CountCm += 1.;
      int vx = (int)(Pos[CLat1]/(2.*Nano->Height)*NBin);
      vx += NBin/2;
      if(vx < 0 || vx >= NBin) continue;
      int vy = (int)(Pos[CLat2]/(2.*Nano->Height)*NBin);
      vy += NBin/2;
      if(vy < 0 || vy >= NBin) continue;
      Count[vx*NBin+vy] += 1.;
    }
    double Area = 0.;
    double Norm = 0.;
    for(int vx=0;vx<NBin;vx++){
      for(int vy=0;vy<NBin;vy++){
    	if(Count[vx*NBin+vy] < 1.) continue;
	// double x = vx*Nano->Height*2./(double)NBin;
	// double y = vy*Nano->Height*2./(double)NBin;
	// fprintf(AreaS,"{t[%d 0 2] x(%lf %lf %lf)}\n",nPart++,x+pNanoPos(0,0)-Nano->Height,y+pNanoPos(0,1)-Nano->Height,pNanoPos(0,2));
    	Area += 1.;
      }
    }
    if(CountCm <= 0.){
      printf("No particles in the torus %s\n",cFile[f]);
      return;
    }
    Nano->Area = SQR(2.*Nano->Height)*Area/(double)(SQR(NBin));
    for(int d=0;d<3;d++){
      Cm[d] /= CountCm;
    }
    Cm[CNorm] = pCm(CNorm);
    //fprintf(AreaS,"%lf %lf %lf\n",Cm[0]-Nano->Pos[0],Cm[1]-Nano->Pos[1],Cm[2]-Nano->Pos[2]);
    for(int d=0;d<3;d++){
      Nano->Pos[d] = Cm[d] - floor(Cm[d]*pInvEdge(d))*pEdge(d);
      fprintf(StalkArea,"%lf %lf\n",pTime(),Nano->Area);
    }
    SetNanoBkf(0);
    Nano->Height = Nano->Rad + sqrt(Nano->Area/DUE_PI);
    char String[120];
    StringNano(String,0);
    fprintf(StalkArea,"%lf %lf\n",pTime(),Nano->Area);
    fprintf(FNano,"sed '/Rigid/c\\%s' %s > Ciccia.dat\n",String,cFile[f]);
    fprintf(FNano,"mv Ciccia.dat %s\n",cFile[f]);
    sprintf(FName,"Centered%05d.dat",f);
    //Write(FName);
    //HeaderNano(AreaS);
  }
  fclose(AreaS);
  fclose(StalkArea);
  fclose(FNano);
  free(Count);
  printf("\n");
}
void ElPoly::AvSnap(){
  int NPairF = 10;
  double InvAv = 1./(double)NPairF;
  double Norm = 1./(double)(NFile[1]-NFile[0]);
  char FName[256];
  // for(int ff=NFile[0];ff<NFile[1];ff+=NPairF){
  PART *Pn = (PART *)calloc(pNPart(),sizeof(PART));
  //for(int f=ff;f<ff+NPairF&&f<NFile[1];f++){
  for(int f=NFile[0];f<NFile[1];f++){
    Processing(f);
    if(OpenRisk(cFile[f],BF_NO))return;
    for(int p=0;p<pNPart();p++){
      for(int d=0;d<3;d++){
	Pn[p].Pos[d] += pPos(p,d);
      }
    }
  }
  for(int p=0;p<pNPart();p++){
    double Pos[3] = {Pn[p].Pos[0]*Norm,Pn[p].Pos[1]*Norm,Pn[p].Pos[2]*Norm};
    SetPos(p,Pos);
  }
  printf("\n");
  free(Pn);
  sprintf(FName,"Av%s",cFile[NFile[0]]);
  Write(FName);
}
void ElPoly::Sample(int NSample){
  MOMENTI m1 = SampleSurfaceMem(NSample);
  FILE *FWrite = fopen("SolSim.dat","w");
  double InvNSample = 1./(double)NSample;
  fprintf(FWrite,"# l(%.2f %.2f %.2f) d[part]\n",pEdge(0),pEdge(1),pEdge(2));
  for(int sx=0;sx<NSample;sx++){
    int sy = 0;
    PlotMem[sx*NSample + sy] = m1.Uno;
    PlotMem[sy*NSample + sx] = m1.Uno;
    sy = NSample - 1;
    PlotMem[sx*NSample + sy] = m1.Uno;
    PlotMem[sy*NSample + sx] = m1.Uno;
  }
  printf("sample average: %lf\n",m1.Uno);
  double ConvFact = 1.45/m1.Uno;
  for(int sx=0;sx<NSample;sx++){
    double x = sx*InvNSample*pEdge(CLat1);
    int sx1 = sx + 1 == NSample ? 0 : sx + 1;
    for(int sy=0;sy<NSample;sy++){
      double y = sy*InvNSample*pEdge(CLat2);
      int sy1 = sy + 1 == NSample ? 0 : sy + 1;
      int p = sx*NSample+sy;
      int l1 = sx1*NSample+sy + SQR(NSample);
      int l2 = sx*NSample+sy1 + SQR(NSample);
      //fprintf(FWrite,"{t[%d 0 2] x(%lf %lf %lf) l[%d] l[%d]}\n",p,x,y,PlotMem[p],sx1*NSample+sy,sx*NSample+sy1);
      fprintf(FWrite,"{t[%d 0 3] x(%lf %lf %lf) l[%d] l[%d]}\n",p,x,y,(PlotMem[p]-m1.Uno)*ConvFact,l1,l2);
      //fprintf(FWrite,"{t[%d 0 2] x(%lf %lf %lf)}\n",p,x,y,PlotMem[p]-m1.Uno);
    }
  }
}
